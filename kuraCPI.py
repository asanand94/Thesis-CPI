'''
file kuraRegressionLeg

Created on Feb 19, 2016

@author: asanand@princeton.edu
'''
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
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
    
# Define residual function
def residual(alpha):
	return rhs - np.dot(hermitePoly, alpha)

#fine variable calculation	
def fineVariables(initial, tmin, tmax, nsteps):
	dt = (tmax-tmin)/nsteps
	times = np.arange(tmin,tmax,dt)
	return	  odeint(
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

theta0 = np.random.rand(n) * 3.1415 #initial distribution of theta

#define small and large time steps as well as number of fine steps taken
smallTstep = .1 #total time we fine integrate for
bigTstep = .1 #jump for the coarse integration
initialT = 0
nsteps = 5 #number of fine steps taken (dtSmall = smallTstep/nsteps)
nHermite = 3 #number of hermite Polynomials used



# Create empty array size nx3 to fill with the omegas multiplied into 1st and 3rd order
# Hermite polynomials, [H1(omega), H3(omega), H5(omega)]
hermitePoly = np.zeros((n,nHermite))

# Fill hermitePoly matrix with H1(omega) and H3(omega) values
hermitePoly[:, 0] = (om).ravel()
hermitePoly[:, 1] = (om ** 3 - 3 * om).ravel()
hermitePoly[:, 2] = (om**5 - 10*(om**3)+15*om).ravel()


alphas = np.zeros((nHermite,nsteps)) #Empty vector to store the fine alphas for each time step

initialIterate = np.zeros((nHermite,1))# Empty initial vector for alpha values of given time step


	
#Loop for CPI		 
for i in range(1):
	
	#begin by calculating fine variable (theta) for a short burst of time
	history = fineVariables(theta0, initialT, initialT+smallTstep, nsteps)
	allTheta = np.transpose(history)
	theta0 = allTheta[:,nsteps-1]

	#create loop to calculate the fine alpha values
	for j in range(nsteps):
	
		rhs = theta0[:]
		thetaMatrix = history.transpose()
		
		#use Least Squares to calculate the fine alphas
		answer, return_code = leastsq(residual, initialIterate)
	
		# enter alpha values for given time step into matrix and then iterate again
		alphas[0,j] = answer[0]
		alphas[1,j] = answer[1]
		alphas[2,j] = answer[2]
		
	#enter last two fine alphas into separate matrix	
	fineAlphas = np.zeros((3,2))
	fineAlphas[:,0] = alphas[:,nsteps*i-1]
	fineAlphas[:,1] = alphas[:,nsteps*i-2]
	
	#create initial alpha vector for coarse integration
	initialAlphas = np.zeros((3,1))
	initialAlphas[:,0] = fineAlphas[:,0]
	
	#create array to calculate d(alpha)/dt to coarse grain
	alphaDerivatives = np.zeros((3,1))
	
	#fill array with d(alpha)/dt values for each alphas
	for k in range(3):
		alphaDerivatives[k,0] = (fineAlphas[k,0]-fineAlphas[k,1])/smallTstep
	
	#send it to coarse function to approximate alphas at t+bigTstep
	coarseIntegrate(initialAlphas, alphaDerivatives, bigTstep)
	

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