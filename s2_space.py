import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors


def two_phos(y, t, tao, w):
    tm, c0, c1, c2 = y

    dtmdt = -w*tm + tao*(c0+c1+c2)
    dc0dt = w*tm-(tao+w)*c0
    dc1dt = w*c0-(tao+w)*c1
    dc2dt = w*c1-tao*c2

    return [dtmdt, dc0dt, dc1dt, dc2dt]


t = np.linspace(0, 180, 100)

tao_ago = 1./18.
tao_noago = 1./3.

dimx = 40
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
        y0 = [ini, 0., 0., 0.]
        sols1 = odeint(two_phos, y0, t, args=(tao_ago, w))
        sols2 = odeint(two_phos, y0, t, args=(tao_noago, w))
        fullyphos_ago = sols1[:, 3][-1]
        fullyphos_noago = sols2[:, 3][-1]
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
