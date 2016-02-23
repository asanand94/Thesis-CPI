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
om.sort()  # This makes plotting easier, without actually affecting dynamics.
K = 2.4

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
    dt = (tmax-tmin)/nsteps
    times = np.arange(tmin,tmax,dt)
    return odeint(
                  dThetaDt,
                  initial,
                  times
                  )

#coarse function    
def coarseIntegrate(initial, dadt, tStep):
    states = np.empty((initial.size,2))
    print states.shape
    states[0,:] = initial[0,:]
    tSteps = np.empty(range(initial.size)); tSteps.fill(tStep)
    states[1,:] = initial[0,:] + dadt[0,:]*tSteps[0,:]
    return states

# theta0 = np.random.rand(n) * 3.1415 #initial distribution of theta
theta0 = np.zeros((n,))

#define small and large time steps as well as number of fine steps taken
smallTstep = .1 #total time we fine integrate for
bigTstep = .1 #jump for the coarse integration
initialT = 0
nsteps = 5 #number of fine steps taken (dtSmall = smallTstep/nsteps)
ncsteps = 12 #number of coarse steps taken
nHermite = 3 #number of hermite Polynomials used



# Create empty array size nx3 to fill with the omegas multiplied into 1st and 3rd order
# Hermite polynomials, [H1(omega), H3(omega), H5(omega)]
hermitePoly = np.zeros((n,nHermite))

# Fill hermitePoly matrix with H0(omega), H2(omega), and H3(omega) values
hermitePoly[:, 0] = (om).ravel()
hermitePoly[:, 1] = (om ** 3 - 3 * om).ravel()
hermitePoly[:, 2] = (om**5 - 10*(om**3)+15*om).ravel()


alphas = np.zeros((nHermite, nsteps)) #Empty vector to store the fine alphas for each time step

    
#Loop for CPI
totalThetaHistory = np.empty((n, ncsteps))
totalAlphaHistory = np.empty((nHermite, ncsteps))
times = []
t = 0
for i in range(ncsteps):
    
    # Creat a trajectory of the fine variable (theta) for a short burst of time.
    history = fineVariables(theta0, initialT, initialT+smallTstep, nsteps)
    t += smallTstep  # Track time just for later plotting purposes.

    # Calculate best-fit alpha values at each fine timestep
    # (though we really only need the last two).
    for j in range(nsteps):
    
        rhs = history[j, :]
            
        # Define residual function
        def residual(alpha):
            return rhs - np.dot(hermitePoly, alpha)
        
        #use Least Squares to calculate the best-fit alphas (coarse variables).
        initialIterate = np.zeros((nHermite, 1))# Empty initial vector for alpha values of given time step
        answer, return_code = leastsq(residual, initialIterate)
    
        # enter alpha values for given time step into matrix and then iterate again
        alphas[0,j] = answer[0]
        alphas[1,j] = answer[1]
        alphas[2,j] = answer[2]
        
    # Extract the last two alpha vectors.
    a1 = alphas[:, -2]
    a2 = alphas[:, -1]
    
    # Save the "current" state (at the end of the fine burst).
    totalThetaHistory[:, i] = history[-1, :]
    totalAlphaHistory[:, i] = a2
    times.append(t)
        
    # Create array to calculate d(alpha)/dt to do coarse Euler integration.
    alphaDerivatives = np.zeros((3,))
    # Fill array with d(alpha)/dt values for each alphas.
    alphaDerivatives[:] = (a2-a1)/smallTstep
        
    # Find next alpha by an Euler step.
    nextAlpha = a2 + alphaDerivatives * bigTstep
    t += bigTstep
    
    # Lift back to the fine fine (theta) space to resume the fine integration.
    theta0 = np.dot(hermitePoly, nextAlpha)
    

## Plot the theta and alpha histories.
fig = plt.figure()
colors = cm.hot(np.linspace(0, 1, ncsteps))

ax = fig.add_subplot(2, 1, 1)
ax.set_xlabel('$\omega_i(t)$')
ax.set_ylabel(r'$\theta_i(t)$')  # The r'... (for "raw") is to prevent Python from interpreting \t as an escape code for a tab character.

bx = fig.add_subplot(2, 1, 2)
bx.set_xlabel('$t$')
bx.set_ylabel(r'$\alpha_k(t)$')

for i, c in zip(range(ncsteps), colors):
    ax.scatter(om, totalThetaHistory[:, i], color=c)
    ax.plot(om, totalThetaHistory[:, i], color='black')
    bx.scatter([times[i]]*nHermite, totalAlphaHistory[:, i], color=c)
    
# Draw lines through coarse points.
bx.plot(times, totalAlphaHistory.T, color='black')

fig.suptitle('$t\in[%.2f, %.2f]$' % (times[0], times[-1]))

# Show the plots.
plt.show()

'''    
Euler Integration
                          
def eulerIntegrate(initial, dxdt, tmin, tmax, nstep):
    nstep = int(nstep)
    states = np.empty((nstep, initial.size))
    states[0, :] = initial

    times = np.linspace(tmin, tmax, nstep)

    for step in bar(xrange(nstep-1)):
        state = states[step, :]
        time = times[step]
        dt = times[step+1] - time
        rhs = dxdt(state, time)
        states[step+1, :] = state + rhs * dt

    return states, times
'''