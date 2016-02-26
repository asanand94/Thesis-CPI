'''
file kuraRegressionLeg

Created on Feb 19, 2016

@author: asanand@princeton.edu
'''
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib.cm as cm

n = 64  # Make a random "network".
A = np.random.rand(n, n)
p = .5
A[A>p] = 1
A[A<=p] = 0  # Make the network symmetric (undirected).
for i in range(n):
    for j in range(n):
        A[i,j] = A[j,i]

# Make a gaussian natural frequency distribution and set the coupling strength.
om = np.random.normal(0,.3,n)
om -= om.mean()
om.sort()  # This makes plotting easier, without actually affecting dynamics.
K = 5.0

#kuramoto oscillators
def dThetaDt(thetaVec, unused_t):
    out = np.zeros((n,))
    for i in range(n):
        sumsini = 0
        for j in range(n):
            sumsini += A[i,j] * np.sin(thetaVec[j] - thetaVec[i])
        out[i] = om[i] + K/n * sumsini
    return out


#fine variable calculation    
def fineVariables(initial, tmin, tmax, nsteps):
    dt = float(tmax-tmin)/nsteps
    times = np.arange(tmin,tmax,dt)
    return times, odeint(
                  dThetaDt,
                  initial,
                  times
                  )


# theta0 = np.random.rand(n) * 3.1415 #initial distribution of theta
theta0 = np.zeros((n,))
times, states = fineVariables(theta0, 0, 4, 10000)

## Restrict the fine states.
# Create empty array size nx3 to fill with the omegas multiplied into 1st and 3rd order
# Hermite polynomials, [H1(omega), H3(omega), H5(omega)]
nHermite = 3
hermitePoly = np.zeros((n,nHermite))

# Fill hermitePoly matrix with H0(omega), H2(omega), and H3(omega) values
hermitePoly[:, 0] = (om).ravel()
hermitePoly[:, 1] = (om ** 3 - 3 * om).ravel()
hermitePoly[:, 2] = (om**5 - 10*(om**3)+15*om).ravel()

alphaHistory = np.zeros((times.size, nHermite))
for i in range(len(times)):
    # Define residual function
    rhs = states[i, :]
    def residual(alpha):
        return rhs - np.dot(hermitePoly, alpha)
            
    #use Least Squares to calculate the best-fit alphas (coarse variables).
    initialIterate = np.zeros((nHermite, 1))# Empty initial vector for alpha values of given time step
    answer, return_code = leastsq(residual, initialIterate)
    alphaHistory[i,:] = answer
## Plot the theta and alpha histories.
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('times')
ax.set_ylabel(r'$\alpha$')  # The r'... (for "raw") is to prevent Python from interpreting \t as an escape code for a tab character.

ax.plot(times, alphaHistory)


fig, ax = plt.subplots()
for i in np.round(np.linspace(0, times.size-1, 10)):
    ax.plot(om, states[i, :])

ax.set_xlabel('$\omega$')
ax.set_ylabel('$\theta_i(t)$')
# Show the plots.
plt.show()
