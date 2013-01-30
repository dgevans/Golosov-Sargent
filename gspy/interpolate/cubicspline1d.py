'''
Authors: Chase Coleman and Spencer Lyon

This file computes the 1 dimensional cubic spline
'''
from __future__ import division
import numpy as np
import scipy.linalg as la
import math
from numba import jit, autojit


def gridbuild(a, b, n):
    '''
    This function takes a,b,n as inputs and will return
    grid and has output
    a: the left endpoint of the interval you want to interpolate over
    b: the right endpoint of the interval you want to interpolate over
    n: the number of equally spaced points you want to use

    grid: an evenly spaced interval in R1 with n+1 nodes
    h: the size of your steps
    '''

    grid = np.linspace(a, b, n + 1)
    h = grid[1] - grid[0]
    return grid, h


# @jit('f8(f8,f8,f8,i8)')
def calcu_k(x, a, h, k=None):
    '''
    This function calculates the value of u_k by applying the
    function phi to the t that is passed in from get_t.  The
    u_k are the basis functions used to evaluate the spline

    x: Value of x we are using
    a: left endpoint
    h: step size
    k: which k you are on
    '''
    t = ((x - a) / h) - (k - 2)

    if np.abs(t) <= 1:
        out = 4. - 6. * np.abs(t) ** 2 + 3. * np.abs(t) ** 3
    elif 1 < np.abs(t) and np.abs(t) <= 2:
        out = (2 - np.abs(t)) ** 3
    else:
        out = 0

    return out


# @autojit
def calccubicspline1d(x, y, alpha, beta, h):
    '''
    This function calculates a cubic spline for a 1 dimensional case
    We proceed on the method explained in Habermann and Kindermann 2007
    It should run faster than typical cubic B-Spline interpolation and
    maintain the accuracy.

    Inputs:
    x: Your x values.  Should be evenly spaced.  Use gridbuild func to choose x
    y: Value of f(x) for the f that you are trying to approximate
    alpha: value of 2nd derivative of f at alpha
    beta: value of 2nd derivative of f at beta

    output:
    c: coefficients s.t. s(x) = sum(c_k * u_k(x))
    '''
    a = x[0]
    b = x[-1]
    x = np.squeeze(x)
    y = np.squeeze(y)
    n = x.size - 1   # Python only indexes up to n-1, but there are n terms in x
    # ie the nth term is x[n-1] b/c python include x[0]

    # Check simple errors.  x and y should be same size and x should be evenly spaced
    if x.size != y.size:
        raise ValueError('x and y should be the same size. Cannot proceed')
    if x[2] - 2 * x[1] + x[0] > 1e-6:
        raise ValueError('x should be evenly spaced.  Cannot proceed')

    # Initialize the u matrices
    umat = np.zeros((n + 3, n + 3))
    for row in xrange(1, x.size + 1):
        for column in xrange(umat.shape[1]):
            umat[row, column] = calcu_k(x[row - 1], a, h, column + 1)

    umat[0, :3] = [1, -2, 1]
    umat[-1, -3:] = [1, -2, 1]

    # We reduce our umat to a smaller tridiagonal matrix by simple mat algebra
    # We will fulfill conditions for c1,c2,cn+2,cn+3
    lumat = umat[2:-2, 2:-2]

    c1 = 1 / 6 * (y[0] - (alpha * h ** 2) / 6)  # put in c[1]. This is c_{2}
    cm2 = 1 / 6 * (y[-1] - (beta * h ** 2) / 6)  # put in c[-2]. this is c_{n + 2}

    # We need to satisfy n+3 conditions so y points and derivs. b does this
    # b = np.zeros(n + 3)
    # b[0] = alpha * h ** 2 / 6.
    # b[-1] = beta * h ** 2 / 6.
    # b[1:-1] = y

    # We have umat * c_n = b
    # To solve this we will first apply simple linear algebra and obtain c2,cn+2
    c_n = np.zeros(n + 3)  # c ranges from 1 to n+3 so we need n+3 spots
    c_n[1] = c1
    c_n[-2] = cm2

    lb = np.zeros(n - 1)
    lb[0] = y[1] - c1
    lb[-1] = y[-2] - cm2
    lb[1:-1] = y[2:-2]

    #We then will solve the middle piece which is lumat * c[middle part] = byvec[middlepart]
    c_n[2:-2] = la.solve(lumat, lb)

    #Then we obtain c_1 and c_n+3 by lin alg
    c_n[0] = ((alpha * h ** 2) / 6.) + 2. * c_n[1] - c_n[2]
    c_n[-1] = ((beta * h ** 2) / 6.) + 2. * c_n[-2] - c_n[-3]

    return c_n


def feval(x, a, h, n, coeffs):
    '''
    This function takes a value x and the coefficients obtained by calccubicspline and
    evaluates the spline at x.

    Inputs:
    x: The point at which you are evaluating the spline
    a: The left endpoint
    h: The step size obtained by gridbuild
    coeffs: The coefficients obtained in calccubicspline

    Outputs:
    y_spline: The estimated value of y.  (Value of the interpolation s(x).)
    '''

    l = math.floor((x - a) / h) + 1
    m = np.min([l + 3, n + 3])

    #We want to build a vector of possible u_ks
    uk_vec = np.zeros(m - l + 1)

    for column in xrange(int(m - l + 1)):
        uk_vec[column] = calcu_k(x, a, h, column + l)

    y_spline = np.dot(uk_vec, coeffs[l - 1:m])

    return y_spline
