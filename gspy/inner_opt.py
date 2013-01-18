"""
Created Dec 22, 2012

Python code that matches MatLab code from  ./InnerOptimizationCode
"""
from __future__ import division
import numpy as np


def computeC2_2(c1_1, c1_2, c2_1, R, s, P, sigma):
    """
    Mimics the file ./InnerOptimizationCode/computeC2_2.m

    Computes c_2(2) and return c1 and c2 as 3x2 matrix.  Here c1 and c2
    will be of the form:
        c_1(1)   c_1(2)
        c_1(1)   c_1(2)
        c_1(1)   c_1(2)

    The 3x2 matrix format is useful for computing future
    gradients.

    gradc1 will a matrix containing the derivative of c1 with respect to
    the various z's.  For instance the ith row and jth column will be
    the derivative of c_1(j) with respect to x(i).  Thus gradc1 will be
           1   0
           0   1
           0   0

    Similarly for gradc2
    """
    # is frac the same as ./steady/steady_state.py line 37 - 38?
    frac = (R * P[s, 0] * c1_1 ** (-sigma) + R * P[s, 1] * c1_2 ** (-sigma) -\
            P[s, 0] * c2_1 ** (-sigma)) / (P[s, 1])
    c2_2 = frac ** (-1. / sigma)

    c1 = np.kron(np.ones((3, 1)), np.array([c1_1, c1_2]))
    c2 = np.kron(np.ones((3, 1)), np.array([c2_1, c2_2]))

    gradc2_2 = np.zeros(3)
    gradc2_2[0] = c1_1 ** (-sigma - 1) * frac ** (-1. / sigma - 1) * R *\
                     P[s, 0] / P[s, 1]
    gradc2_2[1] = c1_2 ** (-sigma - 1) * frac ** (-1. / sigma - 1) * R
    gradc2_2[2] = -c2_1 ** (-sigma - 1) * frac ** (-1. / sigma - 1) * \
                    P[s, 0] / P[s, 1]

    grad_c1 = np.array([[1, 0], [0, 1], [0, 0]])
    grad_c2 = np.empty((3, 2))
    grad_c2[:, 0] = [0, 0, 1]
    grad_c2[:, 1] = gradc2_2

    return c1, c2, grad_c1, grad_c2


def computeR(c1, c2, gradc1, gradc2, sigma):
    """
    Mimics the file ./InnerOptimizationCode/computeR.m

    Computes Rprime c_2(s)^(-sigma)/c_1(s)^(-sigma) in the 3x2 matrix
    format as well as its gradient with respect to z.

    Uses the variables c1, c2 computed in computeC2 as well as their
        gradients.

    Also needed is the primitive sigma.
    """
    r_prime = c2 ** (-sigma) / (c1 ** (-sigma))

    grad_r_prime = sigma * c2 ** (-sigma) * c1 ** (sigma - 1) * gradc1 - \
                    sigma * c2 ** (-sigma - 1) * c1 ** (sigma) * gradc2

    return r_prime, grad_r_prime


def computeL(c1, gradc1, c2, gradc2, Rprime, gradRprime, theta_1,
             theta_2, g, n1, n2):
    """
    Mimics the file ./InnerOptimizationCode/computeL.m

    Computes l_1 and l_2, the labor supply  of agent 1 and 2 in the
    standard 3x2 format, along with their gradients with respect to z.

    Uses c1, c2, Rprime computed using computeC2 and computeRprime as
    well as their gradients.

    Also passed are the primitives theta_1,theta_2, n_1, n_2 and the
    vector of government expenditures g.
    """
    if g.shape[0] > 1:
        g = g.T

    g = np.kron(np.ones((3, 1)), g)

    # Compute l2 first
    l2 = (n1 * c1 + n2 * c2 + g + n1 * theta_2 * Rprime - n1 * theta_1) / \
            (theta_2 * (n2 + Rprime * n1))

    # Now gradl2
    gradl2 = n1 * gradRprime / (n2 + n1 * Rprime) - \
             n1 * gradRprime * l2 / (n2 + n1 * Rprime) + \
             n1 * gradc1 / (theta_2 * (n2 + n1 * Rprime)) + \
             n2 * gradc2 / (theta_2 * (n2 + n1 * Rprime))

    # now l1 and gradl1
    l1 = 1 - (1 - l2) * Rprime * theta_2 / theta_1
    gradl1 = gradl2 * Rprime * theta_2 / \
             theta_1 - (1 - l2) * gradRprime * theta_2 / theta_1

    return l1, gradl1, l2, gradl2


def compute_X_prime(c1, gradc1, c2, gradc2, Rprime, gradRprime, l1,
                    gradl1, l2, gradl2, P, sigma, psi, beta, s, x):
    """
    Mimics the file ./InnerOptimizationCode/computeXprime.m

    Computes the choice of the state variable xprime tomorrow in the
    standard 3x2 format as well as gradient with respect to z (note this
    is unfortunated notation, xprime is refering to xprime
    ($u_{c,2}\tilde b'$), while x is the vector [c_1(1), c_1(2), c_2(1)])


    First create c2 alt.  Here c2alt is a matrix of the form
       c_2(2)  c_2(1)
       c_2(2)  c_2(1)
       c_2(2)  c_2(1)
    This is so when we have multiplications like c2.*c2alt we get
       c_2(2)c_2(1) c_2(1)c_2(2)
       c_2(2)c_2(1) c_2(1)c_2(2)
       c_2(2)c_2(1) c_2(1)c_2(2)
    """
    c2alt = np.fliplr(c2)
    gradc2alt = np.fliplr(gradc2)

    # Euc2 = np.kron(np.ones((1, 2)), (psi * c2 ** (-sigma)).dot(P[s, :].T))
    Euc2 = np.outer(np.ones((1, 2)), (psi * c2 ** (-sigma)).dot(P[s, :].T)).T

    P = np.kron(np.ones((3, 1)), P[s, :])
    Palt = P[:, ::-1]

    # TODO: Check what lines 26-27 do in MatLab version

    xprime = x * psi * c2 ** (-sigma) / (beta * Euc2) + \
            (1 - psi) * l2 / (1 - l2) - \
            (1 - psi) * Rprime * l1 / (1 - l1) + \
            psi * c1 * c2 ** (-sigma) - psi * c2 ** (1 - sigma)

    # TODO: check this and figure out how to follow PEP8
    # SL: Checked on 1/7/13
    gradxprime = (-sigma * x * psi * c2 ** (-sigma - 1) / (beta * Euc2) + \
     (sigma * x * psi ** 2 * c2 ** (-2 * sigma - 1) * P * beta) / ((beta * Euc2) ** 2) - \
     sigma * psi * c2 ** (-sigma - 1) * c1 - (1 - sigma) * psi * c2 ** (-sigma)) * gradc2 + \
     (sigma * x * psi ** 2 * c2 ** (-sigma) * c2alt ** (-sigma - 1) * beta * Palt) / \
     ((beta * Euc2) ** 2) * gradc2alt + psi * c2 ** (-sigma) * gradc1 + \
     (1 - psi) * gradl2 / ((1 - l2) ** 2) - (1 - psi) * Rprime * gradl1 / ((1 - l1) ** 2) -\
     (1 - psi) * l1 * gradRprime / (1 - l1)

    return xprime, gradxprime


def uAlt(c, l, psi, sigma):
    """
    This is a utility function that allows for different values of
    risk aversion
    """
    if sigma == 1:
        u = psi * np.log(c) + (1. - psi) * np.log(1. - l)
    else:
        u = psi * c ** (1. - sigma) / (1. - sigma) + \
            (1. - psi) * np.log(1. - l)

    return u
