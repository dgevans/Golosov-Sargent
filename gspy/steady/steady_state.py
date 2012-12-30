"""
Created Dec 22, 2012

Python code that matches matlab code from directory ./SteadyStateCode
"""
from __future__ import division
import numpy as np
from inneropt.inner_opt import computeC2_2, computeR, computeL, compute_X_prime


def SteadyState_residuals(x, u2bdiff, rr, params, s):
    """
    Mimics the file ./SteadyStateCode/SteadyStateResiduals.m
    """
    p = params
    p.theta = np.array([p.theta_1, p.theta_2])
    p.alpha = np.array([p.alpha_1, p.alpha_2])

    u2btild = u2bdiff
    r = rr
    _s = s
    n1 = p.n1
    n2 = p.n2
    ctol = p.ctol

    psi = p.psi
    beta = p.beta
    P = p.P
    th_1 = p.theta_1
    th_2 = p.theta_2
    g = p.g
    alpha = p.alpha
    sigma = p.sigma

    frac = (r * P[_s, 0] * x[0] ** (-sigma) + r * P[_s, 1] * x[1] ** (-sigma)\
            - P[_s, 0] * x[2] ** (-sigma)) / (P[_s, 1])

    if x.min() > 0 and frac > 0:
        c1_1 = x[0]
        c1_2 = x[1]
        c2_1 = x[2]

        c1, c2, gradc1, gradc2 = computeC2_2(c1_1, c1_2, c2_1, r, _s, P, sigma)
        Rprime, gradRprime = computeR(c1, c2, gradc1, gradc2, sigma)
        l1, gradl1, l2, gradl2 = computeL(c1, gradc1, c2, gradc2, Rprime,
                                     gradRprime, th_1, th_2, g, n1, n2)
        xprime, gradxprime = compute_X_prime(c1, gradc1, c2, gradc2, Rprime,
                                            gradRprime, l1, gradl1, l2, gradl2,
                                            P, sigma, psi, beta, _s, u2btild)

        xprime = xprime[0, :]
        Rprime = Rprime[0, :]

        res = np.empty(3)
        res[0] = xprime[0] - xprime[1]
        res[1] = Rprime[0] - Rprime[1]
        res[2] = xprime[0] - u2btild

        c1 = c1[0, :]
        c2 = c2[0, :]
        l1 = l1[0, :]
        c2 = l2[0, :]

        if np.concatenate([l1, l2]).flatten().max > 1:
            res = np.abs(x) + 100

        if not np.isreal(res).all():
            res = np.abs(res) + 100

    else:
        res = np.abs(x) + 100

    return res, c1, c2, l1, l2

