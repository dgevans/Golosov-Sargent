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

    #Initializes values in matrices for t=0
    xHist[0] = xprime0
    ul20 = (1. - psi) / (1. - l20)
    ul10 = (1. - psi) / (1. - l10)
    uc20 = psi / (c20**(sigma))
    uc10 = psi / (c10**(sigma))
    c1Hist(1) = c10
    c2Hist(1) = c20
    l1Hist(1) = l10
    l2Hist(1) = l20
    btildHist(1) = xprime0 / uc20
    TauHist(1) = 1. - (ul10 / (theta_1 * uc10))
    TransHist(1) = c20 - l20 * ul20 / uc20
    RHist(1) = Rprime0
    YHist(1) = n1 * c10 + n2 * c20 + g[_s - 1]
    sHist(1) = _s
    gHist(1) = g[sHist[0]]
    AfterTaxWageIncome_Agent1Hist[0] = l10 * ul10 / uc10
    AfterTaxWageIncome_Agent2Hist[0] = l20 * ul10 / uc20
    
    for i in xrange(numsim - 1):
        if i % 500 == 0:
            print 'On simulation, t=', i
        x = xHist[i]
        R = RHist[i]
        _s = sHist[i]

        policyrulesinit = getinitialapproxpolicy(np.array([x, R, _s]), domain, policyrulesstore)
        policyrules = CheckGradNAG(inputs)

        #Check to see whether it is already 1d.  If so this is unnecessary
        policyrules = policyrules.squeeze() 

        c1 = policyrules[0:2]
        c2 = policyrules[2:4]
        l1 = policyrules[4:6]
        l2 = policyrules[6:8]

        ul2 = (1. - psi) / (1. - l2)
        uc2 = psi / (c2.**(sigma))
        ul1 = (1. - psi) / (1. - l1)
        uc1 = psi / (c1**(sigma))
        Rprime = policyrules[end - 4:end-2]
        xprime = policyrules[end - 2:]
        btildprime = policyrules[8:10]

        #Interest Rates
        # marginal utility of consumption (Agent 2) in s_
        intnum = psi / (c2Hist[i]**(sigma)) 
        # expected marginal utility of consumption (Agent 2)
        dennum = beta * np.sum(params.P(sHist[i],:) * uc2)
        IntHist[i] = intnum / dennum

        #Tau - From the wage optimality of Agent 2
        tau = 1. - (ul1 / (theta_1 * uc1))

        #Output
        y[0] = c1[0] * n1 + c2[0] * n2 + g[0]
        y[1] = c1[1] * n1 + c2[1] * n2 + g[1]

        # Transfers
        # These are transfers computed on the assumption that Agent 2 cannot
        # borrow and lend.  The transfers are the difference between his
        # consumption and after tax earning (l. U_l / U_c)
        trans = c2 - l2 * ul2 / uc2

        #Income
        aftertaxwageincome_agent2 = l2 * ul2 / uc2 + trans
        aftertaxwageincome_agent1 = l1 * ul1 / uc1 + trans

        #gini coeff
        ginicoeff = (aftertaxwageincome_agent2 + 2. * aftertaxwageincome_agent1) / \
                    (aftertaxwageincome_agent2 + afterwageincome_agent1) - 3./2

        #TODO: Continue programming the runs simulation.  Continue on line 185


def getvalue(x, btild, params, c, V):
    '''This function mimics the file getValue.m from the bellman_equation file'''
    c1 = x[0]
    c2 = x[1]
    n1 = params.n1
    n2 = params.n2
    alpha_1 = params.alpha_1
    alpha_2 = params.alpha_2
    g = params.g[0]
    theta_1 = params.theta_1
    theta_2 = params.theta_2
    psi = params.psi
    beta = params.beta
    sigma = params.sigma
    
    if np.min(x) > 0:
        R = (c1 / c2)**(sigma)
        totalresources = (c1 * n1 + c2 * n2 + g[0])
        denl2 = theta_2 * R * n1 + theta_2 * n2
        l2 = (totalresources - theta_1 * n1 + theta_2 * n1 * R) / (DenL2)
        l1 = 1. - (1. - l2) * theta_2 / theta_1 * R
        xprime = -(c2 - c1) * (psi * c2**(-sigma)) - \
        ((l1 / (1. - l1)) * R - l2 / (1. - l2)) * \
        (1. - psi) + btild * psi * c2**(-sigma)
        Rprime = c2**(-sigma) / c1**(-sigma)
        #Compute the value at time=0 using the policies and continuation values at 
        #T=1 thru C,V
        v=alpha_1 * uAlt(c1,l1,psi,sigma) + alpha_2 * uAlt(c2,l2,psi,sigma)\
        + beta * funeval(c(s_,:).T, V(s_), [xprime Rprime])
        v=-v;
        #funeval evaluates the current points by calling the interpolation function
        #If we change interpolation routines then change funeval to appropriate func

    else:
        v = np.max(abs(x)) * 100 + 10

    return v
