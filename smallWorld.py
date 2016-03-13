'''Create a Small World Network using the Watts-Strogatz method

Created on Mar 08, 2016

@author: asanand@princeton.edu
'''
import numpy as np

N = 10
# Small World Network (Watts Strogatz)
A = np.zeros((N,N))
b = .2 #probability of rewiring a node
    
# Make regular ring matrix. Undirected
# K = 4
for i in range(N):
    for j in range(N):
        if i == j: #No self-connection
            A[i,j] = 0
        elif i < j: #Lower triangle matrix (easier to make symmetric)
            A[i,j] = 0        
        elif abs(i-j) == 1: #immediate neighbor
            A[i,j] = 1
        elif abs(i-j) == 2: #neighbor once removed
            A[i,j] = 1
        elif i == N-1: #loops to bottom of matrix (1st column)
            A[i,0] = 1
            A[i-1,0] = 1
        elif i == N-2: #loops to bottom of matrix (2nd column)
            A[i+1,1] = 1
        else:
            A[i,j] = 0
        
A = A+A.transpose() #Make symmetric

#probability of rewiring node i,j
for i in range(N):
    for j in range(N):
        if np.random.random() < b:
            A[i,j] = 0
            A[j,i] = 0
            k = np.random.randint(N)
            while (k == i) or (k == j): #prevent self-connection
                k = np.random.randint(N)
            A[i,k] = 1
            A[k,i] = 1
print A
print A+A.transpose()