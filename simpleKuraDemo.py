'''
file simpleKuraDemo.py

Created on Nov 12, 2015

@author: bertalan@princeton.edu
'''
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
n = 64  # Make a random "network".
A = np.random.rand(n, n)
p = .5
A[A>p] = 1
A[A<=p] = 0  # Make the network symmetric (undirected).
for i in range(n):
    for j in range(n):
        A[i,j] = A[j,i]

# Make a uniform natural frequency distribution and set the coupling strength.
om = np.random.rand(n, 1) * 1.5
om = om - om.mean()
K = 2.4

def dThetaDt(thetaVec, unused_t):
    out = np.zeros((n,))
    for i in range(n):
        sumsini = 0
        for j in range(n):
            sumsini += A[i,j] * np.sin(thetaVec[j] - thetaVec[i])
        out[i] = om[i] + K/n * sumsini
    return out

# Integrate the IVP.
theta0 = np.random.rand(n) * 3.1415
tmax = 4
times = np.arange(0, tmax, .01)
history = odeint(
                 dThetaDt,
                 theta0,
                 times,
                 )

# Plot the trajectory.
fig, ax = plt.subplots()
ax.plot(times, history, 'k')
ax.set_xlabel('time')
ax.set_ylabel(r'$\theta_i(t)$');

# Plot a couple low-dimensional states vs omegas.
tis = [
        round(.1*len(times)),
        round(.4*len(times)),
        round(.7*len(times)),
        round(.9*times.size),  # times.size gives same value as len(times)
       ]

fig = plt.figure(figsize=(16,9))
for i in range(len(tis)):
    ti = tis[i]
    ax = fig.add_subplot(1, len(tis), i+1)
    ax.scatter(om, history[ti, :])
    ax.set_xlabel(r'$\omega_i$')
    ax.set_ylabel(r'$\theta_i(t=%f)$' % times[ti])
    thetaMin = history.min()
    thetaMax = history.max()
    ax.set_ylim([thetaMin, thetaMax])

# Make figures visible.
plt.show()
