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
    try:
        x = opt.fmin(getValue, np.array([c10 guess, c20guess]), (btild_1, _s, params, c, V))
    except ValueError:
        print 'Optimization failed for V0 once.  Trying with a constrained minimizer'
        lb = np.array([.001, .001])
        ub = np.array([10, 10])
        x = opt.fmin_l_bfgs_b(getValue, x, args = (btild_1, _s, params, c, V), \
            bounds = (lb, ub)
    
    #Establish values
    c10 = x[0]
    c20 = x[1]
    R0 = (c10 / c20)**(sigma)
    totalresources = (c10 * n1 + c20 * n2 + g[_s - 1]
    denl2 = theta_2 * R0 * n1 + theta_2 * n2
    l20 = (totalresources - theta_1 * n1 + theta_2 * n1 *R0) / (DenL2)
    l10 = 1. - (1. - l20) * theta_2 / theta_1 * R0
    xprime0 = -(c20 - c10) * (psi * c20**(-sigma)) - ((l10 / (1. - l10)) * R0\
        - l20 / (1. - l20)) * (1. - psi) + btild_1 * psi * c20**(-sigma)
    Rprime0 = c20**(-sigma) / c10**(-sigma)

    #Initialize matrices
    gHist = np.zeros((NumSim,1))
    xHist = np.zeros((NumSim,1))
    btildHist = np.zeros((NumSim,1))
    RHist = np.zeros((NumSim,1))
    btildHist[0] = btild_1
    TauHist = np.zeros((NumSim,1))
    YHist = np.zeros((NumSim,1))
    TransHist = np.zeros((NumSim,1))
    GMul = np.zeros((NumSim,1))
    c1Hist = np.zeros((NumSim,1))
    c2Hist = np.zeros((NumSim,1))
    l1Hist = np.zeros((NumSim,1))
    l2Hist = np.zeros((NumSim,1))
    sHist = np.zeros((NumSim,1))
    GiniCoeffHist = np.zeros((NumSim,1))
    IntHist = np.zeros((NumSim-1,1))
    IncomeFromAssets_Agent1Hist = np.zeros((NumSim-1,1)
    AfterTaxWageIncome_Agent1Hist = np.zeros((NumSim,1))
    AfterTaxWageIncome_Agent2Hist = np.zeros((NumSim,1))
    GShockDiffHist = np.zeros((NumSim-1,1)
    TransDiffHist = np.zeros((NumSim-1,1)
    LaborTaxAgent1DiffHist = np.zeros((NumSim-1,1))
    LaborTaxAgent2DiffHist = np.zeros((NumSim-1,1)
    DebtDiffHist = np.zeros((NumSim-1,1))

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
    v = 0
    v=-v
    return v
