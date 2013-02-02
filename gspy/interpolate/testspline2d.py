'''
Author: Chase Coleman

Use this to test the 2d cubicspline
'''
from __future__ import division
import numpy as np
import scipy as sp
import scipy.linalg as la
import math
import matplotlib.pyplot as plt
from tridiag import tridiagsolve
from cubicspline2d import *
import mpl_toolkits.mplot3d.axes3d as p3

#This section sets up the left endpoints and right endpoints for
#the x and y values and defines how many intervals we want.
ax = 0.
bx = 4.

ay = 0.
by = 4.

nx = 25
ny = 25


#Sets the values of alpha and beta which are the values of 2nd deriv at
#endpoints.  Alpha=beta=0 is the natural spline
alpha = 0.
beta = 0.

#Define a set of x and y points we will use to interpolate a function z = f(x,y)
x = np.linspace(ax, bx, nx + 1)
y = np.linspace(ay, by, ny + 1)

#Find stepsize
hx = x[1] - x[0]
hy = y[1] - y[0]

#This puts x and y in a form that we can calculate each of the Z_ij points
X, Y = np.meshgrid(x,y)
#We will attempt to match the function Z = sin(x) - cos(y^2).  Thus
#we give ourselves the points evaluated at the given x and y.
Z = np.sin(X) - np.cos(Y**2)

#Calls the interpolation function and finds the coefficients
c_mat = calccubicspline2d(x, y, Z, alpha, beta, hx, hy)

#We want to test the function so we evaluate over many more points than we
#Initially had.
xtest = np.r_[1.0:4.0:100j]
ytest = np.r_[1.0:4.0:100j]

#Initialize a matrix that we will fill with z values
ztest = np.zeros((xtest.size, ytest.size))

#Fill our z matrix with values using feval2 that uses the coefficients from
#earlier to calculate the value at some x_i and y_j to find Z_ij
for i in range(xtest.size):
    for j in range(ytest.size):
        ztest[i, j] = feval2(xtest[i], ax, hx, nx, ytest[j], ay, hy, ny, c_mat)
        # ztest[i, j] = cs2d.eval([xtest[i], ytest[j]])

#Set up x and y to calculate the true values of z for our test
xx, yy = np.meshgrid(xtest,ytest)
zz = np.sin(xx) - np.cos(yy**2)

#Calculate difference between our interpolation values and the true values
err = zz - ztest

#Gives us the maximum error
maxerr = np.max(err)
print maxerr

#Plot each of the graphs
truefig = plt.figure(1)
trueax = p3.Axes3D(truefig)
trueax.plot_surface(xx, yy, zz)
plt.show()

testfig = plt.figure(2)
testax = p3.Axes3D(testfig)
testax.plot_surface(xx, yy, ztest)
plt.show()

#^It doesn't seem to want to plot both at the same time.
#Need to save first one and then exit and save second one
#in order to see both at same time just open them up.