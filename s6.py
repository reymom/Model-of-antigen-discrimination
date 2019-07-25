import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""
def two_phos(y, t, tao, w):
    tm, c0, c1, c2, c3, c4, c5, c6 = y

    dtmdt = -w*tm + tao*(c0+c1+c2+c3+c4+c5+c6)
    dc0dt = w*tm-(tao+w)*c0
    dc1dt = w*c0-(tao+w)*c1
    dc2dt = w*c1-(tao+w)*c2
    dc3dt = w*c2-(tao+w)*c3
    dc4dt = w*c3-(tao+w)*c4
    dc5dt = w*c4-(tao+w)*c5
    dc6dt = w*c5-tao*c6

    return [dtmdt, dc0dt, dc1dt, dc2dt, dc3dt, dc4dt, dc5dt, dc6dt]
"""


def two_phos(y, t, tao, w, act_shp, dis_shp, act_mapk, dis_mapk):
    tm, c0, c1, c2, c3, c4, c5, c6, shp, mapk = y

    dtmdt = -w*tm + tao*(c0+c1+c2+c3+c4+c5+c6)
    dc0dt = w*tm-(tao+w)*c0 + shp*c1 - mapk*c0
    dc1dt = w*c0-(tao+w)*c1 + shp*(c2-c1) + mapk * \
        (c0-c1) - act_shp*c1 + dis_shp*shp
    dc2dt = w*c1-(tao+w)*c2 + shp*(c3-c2) + mapk * \
        (c1-c2) - act_shp*c2 + dis_shp*shp
    dc3dt = w*c2-(tao+w)*c3 + shp*(c4-c3) + mapk * \
        (c2-c3) - act_shp*c3 + dis_shp*shp
    dc4dt = w*c3-(tao+w)*c4 + shp*(c5-c4) + mapk * \
        (c3-c4) - act_shp*c4 + dis_shp*shp
    dc5dt = w*c4-(tao+w)*c5 + shp*(c6-c5) + mapk * \
        (c4-c5) - act_shp*c5 + dis_shp*shp
    dc6dt = w*c5-tao*c6 - shp*c6 + mapk*c5 - act_mapk*c6 + dis_mapk*mapk

    dshpdt = act_shp*(c1+c2+c3+c4+c5)-5*dis_shp*shp
    dmapkdt = act_mapk*c6-dis_mapk*mapk

    return [dtmdt, dc0dt, dc1dt, dc2dt, dc3dt, dc4dt, dc5dt, dc6dt, dshpdt, dmapkdt]


w = 1./12.
t = np.linspace(0, 180, 100)

tao_ago = 1./18.
tao_noago = 1./3.

act_shp = 0.00001*w
dis_shp = 0.00001*w

act_mapk = 0.5*w
dis_mapk = 0.5*w

fullyphos_ago = []
fullyphos_noago = []
num_pmhc = []
ini = 1.
for i in range(27):
    ini = 10**(0.1*i)
    y0 = [ini, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    #sols1 = odeint(two_phos, y0, t, args=(tao_ago, w))
    sols1 = odeint(two_phos, y0, t, args=(
        tao_ago, w, act_shp, dis_shp, act_mapk, dis_mapk))
    #sols2 = odeint(two_phos, y0, t, args=(tao_noago, w))
    sols2 = odeint(two_phos, y0, t, args=(
        tao_noago, w, act_shp, dis_shp, act_mapk, dis_mapk))
    fullyphos_ago.append(sols1[:, 7][-1])
    fullyphos_noago.append(sols2[:, 7][-1])
    num_pmhc.append(ini)

fullyphos_ago = np.asarray(fullyphos_ago)
fullyphos_noago = np.asarray(fullyphos_noago)
num_pmhc = np.asarray(num_pmhc)

# interseccion
for i in range(len(fullyphos_ago)-2):
    if fullyphos_ago[i+1] < 1 and fullyphos_ago[i+2] > 1:
        inter_ago = float(num_pmhc[i+1]+num_pmhc[i+2])/2.
    if fullyphos_noago[i+1] < 1 and fullyphos_noago[i+2] > 1:
        inter_noago = float(num_pmhc[i+1]+num_pmhc[i+2])/2.


plt.figure(1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("# pMHC")
plt.ylabel(r"$6P\;at\;180s$")
plt.plot(num_pmhc, fullyphos_ago, color='red', label=r"$agonist$")
plt.plot(num_pmhc, fullyphos_noago, color='green', label=r"$non-agonist$")
plt.axhline(y=1, linestyle='--', color='black', label="threshold")
# plt.vlines(
#    x=inter_ago, ymin=fullyphos_noago[0], ymax=1, colors='red', linestyles='--')
# plt.vlines(x=inter_noago,
#           ymin=fullyphos_noago[0], ymax=1, colors='green', linestyles='--')
plt.legend()
plt.show()
