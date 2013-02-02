'''
Author: Chase Coleman

This file contains the functions necessary to calculate first and second
derivatives for 1 dimensional cubic splines.
'''

def fderiv(val1, val2, h):
    '''
    This function uses the centered difference formula for calculating a
    derivative.  f'(x) = (f(x+h) - f(x-h))/2h.
    
    Inputs:
    val1: this is the value of f(x+h)
    val2: this is the value of f(x-h)
    h: Step size on grid
    
    output:
    fprime: value of derivative at point x
    '''
    
    fprime = (val1 - val2) / (2. * h)
    
    return fprime

def fderiv2(val1, val2, valx, h):
    '''
    This function uses the centered difference formula for calculating a
    derivative.  f''(x) = (f(x+h) - f(x) + f(x-h)) / (h^2).
    
    Inputs:
    val1: this is the value of f(x+h)
    val2: this is the value of f(x-h)
    valx: this is the value of f(x)
    h: Step size on grid
    
    output:
    f2prime: value of 2nd derivative at point x
    '''
    
    f2prime = (val1 - 2*valx + val2) / (h * h)
    
    return fprime

