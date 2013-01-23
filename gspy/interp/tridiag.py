'''
Author: Chase Coleman

This program takes a near diagonal matrix and uses reduction to quickly solve for the
solution of Ax = b.  For the methodology used see Judd 1998 pg 66
'''
from __future__ import division
import scipy as sp
import numpy as np
import scipy.linalg as la

def tridiagsolve(A, b):
	'''
	This definition takes the matrix A and solution vector b and solves for the
	solution of Ax = b.
	
	A is a square matrix
	b is the solutions vector to Ax = b
	l is the dimension of the near diagonal.
	i.e. how many non zero terms are in the middle rows
	^First I will program the tridiagonal example.  May generalize later
	
	ex: l = 3 would give tridiagonal of form
	[a11 a12 0   0
	 a21 a22 a23 0
	 0   a32 a33 a34
	 0   0   a43 a44]
	 
	 returns x
	 
	 TODO: Figure out why after about 16 values it starts to fall apart.
	'''
	
	if A.shape[0] != A.shape[1]:
		raise ValueError('A must be a square matrix.')
	
	if A.shape[0] != b.size:
		raise ValueError('A and b must have same number of rows')
	
	#Define the size of the matrix and then initialize our c vector that will be filled
	#with values calculated in the way that Judd demonstrated
	n = b.size
	c = np.zeros(n)
	d = np.zeros(n)
	x = np.zeros(n)
	
	#Need to solve for rows 2-3 before we can solve the last n-2
	#We solve them s.t. x_i = c_i + d_i * x1
	c[1] = (b[0] / A[0,1])
	d[1] = (-A[0,0] / A[0,1])
	
	c[2] = (b[1] - A[1,1] * c[1]) / A[1,2]
	d[2] = (-A[1,1]*d[1] - A[1,0]) / A[1,2]
	
	c = np.squeeze(c)
	d = np.squeeze(d)
	x = np.squeeze(d)
	
	for i in range(3,n):
		c[i] = (b[i-1] - A[i-1,i-2]*c[i-2] - A[i-1,i-1]*c[i-1])/A[i-1,i]
		d[i] = (-A[i-1,i-2]*d[i-2] - A[i-1,i-1]*d[i-1])/A[i-1,i]
	
	#Solving like this leaves us the last equation in terms of x1
	#Thus we can solve for x1 and then use it to solve for the other
	#terms
	x[0] = (b[n-1] - A[n-1,n-2]*c[n-2] - A[n-1,n-1]*c[n-1]) \
			/ (A[n-1,n-2]*d[n-2] + A[n-1,n-1]*d[n-1])
	
	for i in range(1,n):
		x[i] = c[i] + (d[i]*x[0])
	
	return x