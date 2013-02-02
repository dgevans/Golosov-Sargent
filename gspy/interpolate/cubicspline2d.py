'''
Authors: Chase Coleman

This file computes the 2 dimensional cubic spline
'''
from __future__ import division
import numpy as np
import scipy as sp
import scipy.linalg as la
from tridiag import tridiagsolve
from cubicspline1d import *


def gridbuild2(ax, ay, bx, by, nx, ny):
    '''
    This function takes a(xy),b(xy),n(xy) as inputs and will return
    grid and has output
    a: the left endpoint of the interval you want to interpolate over.
    b: the right endpoint of the interval you want to interpolate over
    n: the number of equally spaced points you want to use


    grid: an evenly spaced interval in R1 with n+1 nodes
    h: the size of your steps

    x and y denote which variable it refers to for all of the above.
    '''

    xgrid = np.linspace(ax, bx, nx + 1)
    ygrid = np.linspace(ay, by, ny + 1)
    hx = xgrid[1] - xgrid[0]
    hy = ygrid[1] - ygrid[0]
    return xgrid, ygrid, hx, hy


def calcuv_i(x, a, h, i=None):
    '''
    This function calculates the value of u_i and v_i by applying the
    function phi to the t   The u_i and v_iare the basis functions
    used to evaluate the spline

    x: Value of x (or y) we are using
    a: left endpoint
    h: step size
    i: which step you are on
    '''
    #Determine t
    t = ((x - a) / h) + 2. - i

    #Calculate the value of the function phi(t)
    if np.abs(t) <= 1:
        out = 4. - 6.*np.abs(t)**2 + 3.*np.abs(t)**3
    elif 1 < np.abs(t) <= 2:
        out = (2 - np.abs(t))**3
    else:
        out = 0

    return out


def calccubicspline2d(x, y, z, alpha, beta, hx, hy):
    '''
    This function calculates a cubic spline for a 2 dimensional case
    We proceed on the method explained in Habermann and Kindermann 2007
    It should run faster than typical cubic B-Spline interpolation and
    maintain or even improve the accuracy.

    Inputs:
    x: Your x values.  Should be evenly spaced.  Use gridbuild func to choose x
    y: Your y values.  Should be evenly spaced.  Use gridbuild func to choose y
    z: Value of f(x,y) for the f that you are trying to approximate
    alpha: value of 2nd derivative of f at alpha
    beta: value of 2nd derivative of f at beta

    output:
    c: coefficients s.t. s(x) = sum(c_ij * u_i(x) * v_j(x))
    '''
    #Pull out necessary information from the arrays.
    ax = x[0]
    bx = x[-1]
    ay = y[0]
    by = y[-1]

    x = np.squeeze(x)
    y = np.squeeze(y)

    #Python only indexes up to n-1, but there are n terms in x
    #ie the nth term is x[n-1] b/c python include x[0]
    nx = x.size - 1
    ny = y.size - 1

    #Check simple errors.  x and y should be same size and x should be evenly spaced
    if x.size*y.size != z.size:
        raise ValueError('z should have a size of  nx * ny. Cannot proceed')
    if x[2] - 2*x[1] + x[0] > 1e-6:
        raise ValueError('x should be evenly spaced.  Cannot proceed')
    if y[2] - 2*y[1] + y[0] > 1e-6:
        raise ValueError('y should be evenly spaced.  Cannot proceed')
    #Initialize the u matrices

    #Create a cstar matrix.  Then fill it with coefficients calculated by running the 1d
    #cubic spline over every possible y and letting x change
    cstar = np.zeros((nx+3, ny + 1))

    for i in xrange(y.size):
        ctemp = calccubicspline1d(x, z[:,i], alpha, beta, hx)
        cstar[:, i] = ctemp.squeeze()

    #Create a true coeff mat.  Then fill it by running the 1d cubic spline fixed at all of the points
    #in y.

    c_mat = np.zeros((nx + 3, ny + 3))

    for i in xrange(nx + 3):
        ctemp = calccubicspline1d(y, cstar[i,:], alpha, beta, hy)
        c_mat[i, :] = ctemp

    return c_mat



def feval2(x, ax, hx, nx, y, ay, hy, ny, coeffs):
    '''
    This function takes a value x and the coefficients obtained by calccubicspline and
    evaluates the spline at x.

    Inputs:
    x,y: The point at which you are evaluating the spline
    a(xy): The left endpoint
    h(xy): The step size obtained by gridbuild2
    coeffs: The coefficients obtained in calccubicspline2

    Outputs:
    z_spline: The estimated value of z.  (Value of the interpolation s(x).)
    '''
    #We want to build a vector of possible u_ks
    uk_temp = 0.
    vk_temp = 0.
    tempsum = 0.

    lx = math.floor((x-ax)/hx) + 1
    ly = math.floor((y-ay)/hy) + 1
    mx = np.min([lx + 3, nx + 3])
    my = np.min([ly + 3, ny + 3])

    #We want to calculate the values of the functions from lx to mx and ly to my
    #We set it up like this to make it easier to see range
    for row in xrange(int(mx - lx +1)):
        for column in xrange(int(my - ly + 1)):
            #we take coeffs_row,column for each of the possible values in U_x X V_y
            #And multiply it by the basis function value.  Same as in 1d, but in 2d
            tempsum += coeffs[row + lx - 1, column + ly - 1] \
            * calcuv_i(x, ax, hx, row + lx) * calcuv_i(y, ay, hy, column + ly)

    z_spline = tempsum

    return z_spline