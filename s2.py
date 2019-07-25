import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def two_phos(y, t, tao, w):
    tm, c0, c1, c2 = y

    dtmdt = -w*tm + tao*(c0+c1+c2)
    dc0dt = w*tm-(tao+w)*c0
    dc1dt = w*c0-(tao+w)*c1
    dc2dt = w*c1-tao*c2

    return [dtmdt, dc0dt, dc1dt, dc2dt]


w = 1./12.
t = np.linspace(0, 180, 100)

tao_ago = 1./18.
tao_noago = 1./3.

fullyphos_ago = []
fullyphos_noago = []
num_pmhc = []
ini = 1.
for i in range(40):
    ini = 10**(0.1*i)
    y0 = [ini, 0., 0., 0.]
    sols1 = odeint(two_phos, y0, t, args=(tao_ago, w))
    sols2 = odeint(two_phos, y0, t, args=(tao_noago, w))
    fullyphos_ago.append(sols1[:, 3][-1])
    fullyphos_noago.append(sols2[:, 3][-1])
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
plt.ylabel(r"$2P\;at\;180s$")
plt.plot(num_pmhc, fullyphos_ago, color='red', label=r"$agonist$")
plt.plot(num_pmhc, fullyphos_noago, color='green', label=r"$non-agonist$")
plt.axhline(y=1, linestyle='--', color='black', label="threshold")
plt.vlines(
    x=inter_ago, ymin=fullyphos_noago[0], ymax=1, colors='red', linestyles='--')
plt.vlines(x=inter_noago,
           ymin=fullyphos_noago[0], ymax=1, colors='green', linestyles='--')
plt.legend()
plt.show()
