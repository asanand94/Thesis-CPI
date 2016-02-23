'''
file kuraCoeffSum

Created on Feb 14, 2016

@author: asanand@princeton.edu
'''
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
n = 256  # Make a random "network".
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
tmax = 30
times = np.arange(0, tmax, .01)
history = odeint(
                 dThetaDt,
                 theta0,
                 times,
                 )
            
t = len(history) # Total number of time steps
alphas = np.zeros((2,t)) # Empty vector for alpha values for ALL time steps

# Transpose history matrix to get an n x 400 matrix where each column is a time step
thetaMatrix = history.transpose()

# Create empty array size nx2 to fill with the omegas multiplied into 1st and 3rd order
# Hermite polynomials, [H1(omega), H3(omega)]
hermitePoly = np.zeros((n,2))

# Fill hermitePoly matrix with H1(omega) and H3(omega) values
hermitePoly[:, 0] = (om).ravel()
hermitePoly[:, 1] = (om ** 3 - 3 * om).ravel()

#for loop to calculate alphas
for i in range(t):

	sum1 = 0
	sum2 = 0
	thetaAtT = thetaMatrix[:,i] #thetas for given t
	
	#loop to sum H1/H3*theta for all thetas at that time step to get alphas
	for j in range(n):
		
		sum1 = sum1 + thetaAtT[j]*hermitePoly[j,0] #calculate alpha for H1
		sum2 = sum2 + thetaAtT[j]*hermitePoly[j,1] #calculate alpha for H3
		
	# enter alpha values for given time step into matrix and then iterate again
	alphas[0,i] = sum1
	alphas[1,i] = sum2

# Transpose the alpha vector so it can be plotted against time.
alphaTrans = alphas.transpose()

# Plot the alpha trajectory and show plot.
fig, ax = plt.subplots()
ax.plot(times, alphaTrans, 'k')
ax.set_xlabel('time')
ax.set_ylabel('alpha');
plt.show()
