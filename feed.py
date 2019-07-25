import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


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


t = np.linspace(0, 180, 100)


w = 1./12.
act_shp = 0.1*w
dis_shp = 0.5*w

act_mapk = 0.5*w
dis_mapk = 0.6*w


fullyphos1 = []
fullyphos2 = []
fullyphos3 = []
timescales = []
rango = 41
for i in range(rango):
    print("j=%i de %i" % (i+1, rango))
    timescale = 10**(-0.5+0.1*i)
    timescales.append(timescale)
    tao = 1./float(timescale)
    ini1 = 1.
    # ini2 = 100
    # ini3 = 100000
    y0_1 = [ini1, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    # y0_2 = [ini2, 0., 0., 0., 0., 0.]
    # y0_3 = [ini3, 0., 0., 0., 0., 0.]
    sols1 = odeint(two_phos, y0_1, t, args=(
        tao, w, act_shp, dis_shp, act_mapk, dis_mapk))
    # sols2 = odeint(two_phos, y0_2, t, args=(
    #    tao, w, act_shp, dis_shp, act_mapk, dis_mapk))
    # sols3 = odeint(two_phos, y0_3, t, args=(
    #    tao, w, act_shp, dis_shp, act_mapk, dis_mapk))
    fin = sols1[:, 7][-1]/0.56
    print(fin)
    fullyphos1.append(fin)
    # fullyphos2.append(sols2[:, 3][-1]/ini2)
    # fullyphos3.append(sols3[:, 3][-1]/ini3)

act_shp = 0*w
dis_shp = 0*w

act_mapk = 0*w
dis_mapk = 0*w

fullyphos2 = []
rango = 41
for i in range(rango):
    print("j=%i de %i" % (i+1, rango))
    timescale = 10**(-0.5+0.1*i)
    tao = 1./float(timescale)
    ini1 = 1.
    y0_1 = [ini1, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
    sols1 = odeint(two_phos, y0_1, t, args=(
        tao, w, act_shp, dis_shp, act_mapk, dis_mapk))
    fin = sols1[:, 7][-1]
    print(fin)
    fullyphos2.append(fin)

"""
w = 1./12.
act_shp = 0.01*w
act_mapk = 0.01*w

tao = 1./18.
# tao_noago = 1./3.
y0 = [1., 0., 0., 0., 0., 0.]
sols1 = odeint(two_phos, y0, t, args=(
    tao, w, act_shp, act_mapk))
plt.figure(0)
plt.xlabel("t")
plt.ylabel(r"$C_{i}$")
plt.plot(t, sols1[:, 0], label='tm')
plt.plot(t, sols1[:, 1], label='c0')
plt.plot(t, sols1[:, 2], label='c1')
plt.plot(t, sols1[:, 3], label='c2')
plt.plot(t, sols1[:, 4], label='shp')
plt.plot(t, sols1[:, 5], label='mapk')
plt.legend(loc='best')
plt.show()
"""
# what i want to graphicate

plt.figure(1)
plt.xscale('log')
plt.xlabel("timescale")
plt.ylabel(r"$2P\;at\;180s$")
plt.plot(timescales, fullyphos2, label=r"$Simple\;KP$")
plt.plot(timescales, fullyphos1, label=r"$KP\;with\;feedbacks$")
# plt.plot(timescales, fullyphos3, label=r"$10^5\;pMHC$")
plt.legend()
plt.show()
