'''
Chase Coleman and Spencer Lyon

This file mimics the matlab folder RunSimulations.m
'''

from __future__ import division
from itertools import product
import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
import scipy.interpolate as interp
import time
from compeconpy import fundefn, funfitxy
from steady_state import steady_state_res, find_steady_state
from inner_opt import uAlt
from set_params import DotDict

def runsimulations(coeffs, btild0, c10guess, c20guess, numsim, params, rhist0 = None):
    
    if rhist0 == None:
        print 'Not using existing shocks'
    else:
        print 'Using existing shocks'
        
    olddatapath = params.datapath
    oldtexpath = params.texpath
    oldplotpath =params.plotpath
    plotpath = oldplotpath
    datapath = olddatapath
    
    print 'Govt Exp'
    g = params.g
    n1 = params.n1
    n2 = params.n2
    alpha_1 = params.alpha_1
    alpha_2 = params.alpha_2
    theta_1 = params.theta_1
    theta_2 = params.theta_2
    psi = params.psi
    beta = params.beta
    sigma = params.sigma
    # SOLVE THE T-0 PROBLEM given btild(-1)
    btild_1 = btild0
    _s = 1
    
    print 'Computed V... Now solving V0(btild_1) where btild_1 is', btild_1
    
    #c1 and c2 solve
    

def getvalue(x, btild, params, c, V):
    '''This function mimics the file getValue.m from the Auxiliary file'''
    c1 = x[0]
    c2 = x[1]
    n1 = params.n1
    n2 = params.n2
    alpha_1 = params.alpha_1
    alpha_2 = params.alpha_2
    g = params.g[0]
    theta_1 = params.theta_1
    theta_2 = params.theta_2
    sigma = params.sigma
    gamma = params.gamma
    
    TotalResoucrces = (c1 * n1 + c2 * n2 + g[0])
    DenL1 = theta_1 * n1 + (c2 / c1)**(-sigma / gamma) * \
    (theta_2 / theta_1)**(1. / gamma) * n2 * theta_2
    l1 = TotalResoucrces / DenL1
    l2 = (c2 / c1)**(-sigma / gamma) * (theta_2 / theta_1)**(1. / gamma) * l1
    
    btildprime = c2**(sigma) * btild / beta2 - (c2 - c1 - l2**(1. + gamma)\
                * c2**(sigma) + l1**(1. + gamma) * c1**(sigma))
    u2btildprime = c2**(-sigma) * btildprime
    Rprime = c2**(-sigma) / c1**(-sigma)
    #v = alpha_1 * u(c1,l1,sigma,gamma)+alpha_2*u(c2,l2,sigma,gamma)+beta*beta*funeval(c(1,:),V(1),[u2btildprime Rprime]);
    #^^Need to find the function u(c1, l1, sigma, gamma)
    v=-v