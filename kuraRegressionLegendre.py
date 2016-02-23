'''
file kuraRegressionLeg

Created on Feb 17, 2016

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
tmax = 30
times = np.arange(0, tmax, .01)
history = odeint(
                 dThetaDt,
                 theta0,
                 times,
                 )
               
# Create empty array size nx2 to fill with the omegas multiplied into 1st and 3rd order
# Hermite polynomials, [H1(omega), H3(omega)]
hermitePoly = np.zeros((n,5))

# Fill hermitePoly matrix with H1(omega) and H3(omega) values
hermitePoly[:, 0] = (om * 0 + 1).ravel()
hermitePoly[:, 1] = (om*2.3094).ravel()
hermitePoly[:, 2] = (5.96285*(om**2) - 1.11803).ravel()
hermitePoly[:, 3] = (15.6785*(om**3) - 5.2915*om).ravel()
hermitePoly[:, 4] = (41.4815*(om**4) - 20*(om**2) + 1.125).ravel()
# Transpose history matrix to get an n x 400 matrix where each column is a time step
thetaMatrix = history.transpose()

t = len(history) # Total number of time steps

alphas = np.zeros((5,t)) # Empty vector for alpha values for ALL time steps

initialIterate = np.zeros((5,1))# Empty initial vector for alpha values of given time step

# Define residual function
def residual(alpha):
	return rhs - np.dot(hermitePoly, alpha)

# Calculate the alpha value for all time steps
for i in range(t):

# Create matrix of theta values for given time step
	rhs = thetaMatrix[:,i]

# Use leastsq to get alpha values for give time step
	answer, return_code = leastsq(residual, initialIterate)
	
	# enter alpha values for given time step into matrix and then iterate again
	alphas[0,i] = answer[0]
	alphas[1,i] = answer[1]
	alphas[2,i] = answer[2]
	alphas[3,i] = answer[3]
	alphas[4,i] = answer[4]
	
# Transpose the alpha vector so it can be plotted against time.
alphaTrans = alphas.transpose()

# Plot the alpha trajectory and show plot.
fig, ax = plt.subplots()
ax.plot(times, alphaTrans, 'k')
ax.set_xlabel('time')
ax.set_ylabel('alpha');
plt.show()