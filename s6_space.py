import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors

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


tao_ago = 1./18.
tao_noago = 1./3.

act_shp = 0.00001/12.
dis_shp = 0.00001/12.

act_mapk = 0.5/12.
dis_mapk = 0.5/12.

t = np.linspace(0, 180, 100)

tao_ago = 1./18.
tao_noago = 1./3.

dimx = 27
dimy = 298
space = np.zeros((dimy, dimx))
timescales = []
for j in range(dimy):
    print("j=%i de %i" % (j,dimy))
    timescale = 0.3 + 0.1*j
    timescales.append(timescale)
    w = 1./float(timescale)
    for i in range(dimx):
        ini = 10**(0.1*i)
        y0 = [ini, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        sols1 = odeint(two_phos, y0, t, args=(
            tao_ago, w, act_shp, dis_shp, act_mapk, dis_mapk))
        sols2 = odeint(two_phos, y0, t, args=(
            tao_noago, w, act_shp, dis_shp, act_mapk, dis_mapk))
        fullyphos_ago = sols1[:, 7][-1]
        fullyphos_noago = sols2[:, 7][-1]
        if fullyphos_ago < 1 and fullyphos_noago <= 1:
            space[j, i] = 2  # black
        if fullyphos_ago >= 1 and fullyphos_noago < 1:
            space[j, i] = 1  # grey
        if fullyphos_ago >= 1 and fullyphos_noago >= 1:
            space[j, i] = 0  # white
timescales = np.asarray(timescales)
num_pmhc = []
for i in range(dimx):
    ini = 10**(0.1*i)
    num_pmhc.append(ini)
num_pmhc = np.asarray(num_pmhc)

cmap = colors.ListedColormap(['white', 'gray', 'black'])

"""
fig, ax = plt.subplots()
plt.xscale('log')
plt.pcolor(space, cmap=cmap)
plt.yticks([1,3,10,30])
plt.ylabel('Timescale',rotation='vertical')
plt.show()
"""

print(len(num_pmhc))
print(num_pmhc)
print(len(timescales))
print(timescales)
print(len(space))
print(space)

fig, ax = plt.subplots()
cax = ax.pcolor(num_pmhc, timescales, space, cmap=cmap)
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Timescale(s)')
plt.xlabel('# pMHC')
#plt.axis([1, 0.0001, 2, 21])
# plt.yticks([2,4,6,8,10,12,14,16,18,20])
# plt.ylabel('<k>',rotation='vertical')
ax.minorticks_on()
# ax.tick_params(axis='y',which='minor',direction='out')
#cbar = fig.colorbar(cax,ticks=[0,25,50,75,100,125,150,175])
# cbar.ax.set_yticklabels([0,25,50,75,100,125,150,'N.C'])
plt.show()