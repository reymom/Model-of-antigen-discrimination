import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def complexes(y, t, k1, kp, kminus1):
    tm, c0, c1, c2, c3, c4, c5 = y

    dtmdt = -k1*tm + kminus1*(c0+c1+c2+c3+c4+c5)
    dc0dt = k1*tm-(kminus1+kp)*c0
    dc1dt = kp*c0-(kminus1+kp)*c1
    dc2dt = kp*c1-(kminus1+kp)*c2
    dc3dt = kp*c2-(kminus1+kp)*c3
    dc4dt = kp*c3-(kminus1+kp)*c4
    dc5dt = kp*c4-kminus1*c5

    return [dtmdt, dc0dt, dc1dt, dc2dt, dc3dt, dc4dt, dc5dt]


k1 = 0.1
kp = 0.1
y0 = [1., 0., 0., 0., 0., 0., 0.]
t = np.linspace(0, 125, 1000)

kminus1 = 0.5*kp
sols1 = odeint(complexes, y0, t, args=(k1, kp, kminus1))

kminus1 = kp
sols2 = odeint(complexes, y0, t, args=(k1, kp, kminus1))

tm1 = np.empty([len(sols1[:, 0])])
tm2 = np.empty([len(sols2[:, 0])])
c01 = np.empty([len(sols1[:, 1])])
c02 = np.empty([len(sols2[:, 1])])
c11 = np.empty([len(sols1[:, 2])])
c12 = np.empty([len(sols2[:, 2])])
c21 = np.empty([len(sols1[:, 3])])
c22 = np.empty([len(sols2[:, 3])])
c31 = np.empty([len(sols1[:, 4])])
c32 = np.empty([len(sols2[:, 4])])
c41 = np.empty([len(sols1[:, 5])])
c42 = np.empty([len(sols2[:, 5])])
active1 = np.empty([len(sols1[:, 6])])
active2 = np.empty([len(sols2[:, 6])])

print(sols1[:, 6][-1])
print(sols2[:, 6][-1])

for i in range(len(active1)):
    tm1[i] = 0
    tm2[i] = 0
    c01[i] = 0
    c02[i] = 0
    c11[i] = 0
    c12[i] = 0
    c21[i] = 0
    c22[i] = 0
    c31[i] = 0
    c32[i] = 0
    c41[i] = 0
    c42[i] = 0
    active1[i] = 0
    active2[i] = 0
    if i > 0:
        if sols1[:, 0][i-1] < sols1[:, 0][i]:
            tm1[i] = sols1[:, 0][i]-sols1[:, 0][i-1]
        else:
            tm1[i] = 0
        if sols2[:, 0][i-1] < sols2[:, 0][i]:
            tm2[i] = sols2[:, 0][i]-sols2[:, 0][i-1]
        else:
            tm2[i] = 0
        if sols1[:, 1][i-1] < sols1[:, 1][i]:
            c01[i] = sols1[:, 1][i]-sols1[:, 1][i-1]
        else:
            c01[i] = 0
        if sols2[:, 1][i-1] < sols2[:, 1][i]:
            c02[i] = sols2[:, 1][i]-sols2[:, 1][i-1]
        else:
            c02[i] = 0
        if sols1[:, 2][i-1] < sols1[:, 2][i]:
            c11[i] = sols1[:, 2][i]-sols1[:, 2][i-1]
        else:
            c11[i] = 0
        if sols2[:, 2][i-1] < sols2[:, 2][i]:
            c12[i] = sols2[:, 2][i]-sols2[:, 2][i-1]
        else:
            c12[i] = 0
        if sols1[:, 3][i-1] < sols1[:, 3][i]:
            c21[i] = sols1[:, 3][i]-sols1[:, 3][i-1]
        else:
            c21[i] = 0
        if sols2[:, 3][i-1] < sols2[:, 3][i]:
            c22[i] = sols2[:, 3][i]-sols2[:, 3][i-1]
        else:
            c22[i] = 0
        if sols1[:, 4][i-1] < sols1[:, 4][i]:
            c31[i] = sols1[:, 4][i]-sols1[:, 4][i-1]
        else:
            c31[i] = 0
        if sols2[:, 4][i-1] < sols2[:, 4][i]:
            c32[i] = sols2[:, 4][i]-sols2[:, 4][i-1]
        else:
            c32[i] = 0
        if sols1[:, 5][i-1] < sols1[:, 5][i]:
            c41[i] = sols1[:, 5][i]-sols1[:, 5][i-1]
        else:
            c41[i] = 0
        if sols2[:, 5][i-1] < sols2[:, 5][i]:
            c42[i] = sols2[:, 5][i]-sols2[:, 5][i-1]
        else:
            c42[i] = 0
        if sols1[:, 6][i-1] < sols1[:, 6][i]:
            active1[i] = sols1[:, 6][i]-sols1[:, 6][i-1]
        else:
            active1[i] = 0
        if sols2[:, 6][i-1] < sols2[:, 6][i]:
            active2[i] = sols2[:, 6][i]-sols2[:, 6][i-1]
        else:
            active2[i] = 0

# veo como evolucionan todos los complejos (poner sols1, sols2)
"""
plt.figure(0)
plt.xlabel("t")
plt.ylabel(r"$C_{i}$")
plt.plot(t, sols2[:, 0], label='tm')
plt.plot(t, sols2[:, 1], label='c0')
plt.plot(t, sols2[:, 2], label='c1')
plt.plot(t, sols2[:, 3], label='c2')
plt.plot(t, sols2[:, 4], label='c3')
plt.plot(t, sols2[:, 5], label='c4')
plt.plot(t, sols2[:, 6], label='c5')
plt.legend(loc='best')
plt.show()
"""
# comparacion activos para diferente k-1
"""
plt.figure(1)
plt.xlabel("t")
plt.ylabel(r"$C_{N}$")
plt.plot(t, sols1[:, 6], label=r"$k_{-1}=0.1·k_{p}$")
plt.plot(t, sols2[:, 6], label=r"$k_{-1}=k_{p}$")
plt.legend(loc='best')
plt.show()
"""
plt.figure(1)
plt.xlabel("t")
plt.ylabel(r"$density\;variation$")
plt.xlim(1, 100)
plt.ylim(0, 0.0008)
plt.plot(t, active1, label=r"$C_{N}\;for\;k_{-1}=0.5·k_{p}$")
plt.plot(t, tm1+c01+c11+c21+c31+c41+active1,
         label=r"$All\;for\;k_{-1}=0.5·k_{p}$")
plt.plot(t, active2, label=r"$C_{N}\;for\;k_{-1}=k_{p}$")
plt.plot(t, tm2+c02+c12+c22+c32+c42+active2, label=r"$All\;for\;k_{-1}=k_{p}$")
plt.legend(loc='best')
plt.show()