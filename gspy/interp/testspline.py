'''
Author: Chase Coleman

Use this to test 1dcubicspline
'''
from __future__ import division
import numpy as np
import scipy as sp
import scipy.linalg as la
import math
import matplotlib.pyplot as plt
from tridiag import tridiagsolve
from cubicspline1d import *

#Here we initialize the left and right endpoint and the number of intervals we want
a = 0
b = 4
n = 25

#This calls gridbuild.  Gridbuild builds an evenly spaced interval
x,h = gridbuild(a,b,n)

#We will try and estimate the sin function for each of our x's so
#we calculate the real values for the x's that we are given.
y = np.sin(x)

#2nd derivatives at the endpoints.  alpha,beta = 0 results in a natural spline
alpha = 0
beta = 0

#Call the calccubicspline1d function that calculates the cubic spline for 1d
coeffs = calccubicspline1d(x,y,alpha,beta,h)

#We want to interpolate to 500 points with the 25 from earlier.
#This initializes the x values and creates a y mat to fill.
xtest = np.r_[1.0:4.0:500j]
ytest = np.zeros(500)

#Evaluates each of the x values from xtest using the coeffs found by the interpolation
#function
for i in range(500):
	ytest[i] = feval(xtest[i], a, h, n, coeffs)

#Calculates the error between the real sin(x) values and the values our interpolation gave us
err = np.sin(xtest) - ytest

#plots both functions on one graph.  It will be hard to notice the difference, but
#if you zoom in enough you will see both a blue and green graph.
#Green is true, blue is the feval values.
plt.plot(xtest, ytest, xtest, np.sin(xtest))
plt.show()