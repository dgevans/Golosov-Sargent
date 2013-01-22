"""
Created Dec 22, 2012

Python code that matches MatLab code from  ./InnerOptimizationCode
"""
from __future__ import division
import numpy as np
import numexpr as ne
from set_params import DotDict


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
    l2 = ne.evaluate("(n1 * c1 + n2 * c2 + g + n1 * theta_2 * Rprime - n1 * theta_1) / \
            (theta_2 * (n2 + Rprime * n1))")

    # Now gradl2
    gradl2 = ne.evaluate("n1 * gradRprime / (n2 + n1 * Rprime) - \
             n1 * gradRprime * l2 / (n2 + n1 * Rprime) + \
             n1 * gradc1 / (theta_2 * (n2 + n1 * Rprime)) + \
             n2 * gradc2 / (theta_2 * (n2 + n1 * Rprime))")

    # now l1 and gradl1
    l1 = 1 - (1 - l2) * Rprime * theta_2 / theta_1
    gradl1 = ne.evaluate("gradl2 * Rprime * theta_2 / \
             theta_1 - (1 - l2) * gradRprime * theta_2 / theta_1")

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

    xprime = ne.evaluate("x * psi * c2 ** (-sigma) / (beta * Euc2) + \
            (1 - psi) * l2 / (1 - l2) - \
            (1 - psi) * Rprime * l1 / (1 - l1) + \
            psi * c1 * c2 ** (-sigma) - psi * c2 ** (1 - sigma)")

    # TODO: check this and figure out how to follow PEP8
    # SL: Checked on 1/7/13
    gradxprime = ne.evaluate("(-sigma * x * psi * c2 ** (-sigma - 1) / (beta * Euc2) + \
     (sigma * x * psi ** 2 * c2 ** (-2 * sigma - 1) * P * beta) / ((beta * Euc2) ** 2) - \
     sigma * psi * c2 ** (-sigma - 1) * c1 - (1 - sigma) * psi * c2 ** (-sigma)) * gradc2 + \
     (sigma * x * psi ** 2 * c2 ** (-sigma) * c2alt ** (-sigma - 1) * beta * Palt) / \
     ((beta * Euc2) ** 2) * gradc2alt + psi * c2 ** (-sigma) * gradc1 + \
     (1 - psi) * gradl2 / ((1 - l2) ** 2) - (1 - psi) * Rprime * gradl1 / ((1 - l1) ** 2) -\
     (1 - psi) * l1 * gradRprime / (1 - l1)")

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


def check_grad(xx, rr, ss, c, vv, z_init, params):
    """
    Mimics the file ./InnerOptimizationCode/CheckGradNag.m

    THIS FUNCTION PERFORMS THE INNER OPTIMIZATION USING THE NAG LIBRARY
    The arguments are explained as follows

    Parameters
    ----------
    xx ,RR, ss: scalars, dtype=float
        A point in the domain

    c: array-like, dtype=float
        Current guess at coefficients

    vv:
        The functional space created in the CompEcon toolbox.

    zInit: array-like, dtype=float
         Initial guess for optimal policies

    params:
        Dictionary containing parameters for the mode

    """
    z_init = z_init[:3]
    params.theta = [params.theta_1, params.theta_2]
    params.alpha = [params.alpha_1, params.alpha_2]
    par = params
    x = xx
    r = rr
    v_coef = [c[0, :].T, c[1, :].T]
    v = vv
    _s = ss
    xLL = par.xLL
    xUL = par.xUL
    n1 = par.n1
    n2 = par.n2
    ctol = par.ctol

    # NOTE: I will pass the globs object through scipy.optimize.root to
    #       bel_obj_uncond_grad

    globs = DotDict
    globs.V = v
    globs.Vcoef = v_coef
    globs.r = r
    globs.x = x
    globs.params = par
    globs._s = _s

    # NOTE: As of scipy version 0.11.0 scipy.optimize.root can be used to
    #       call the same routine used by the function c05qb within NAG.
    #       This function is HYBRD1.f and is a modification of the
    #       hyrbid powell method.


def bel_obj_uncond_grad(n, z, user, iflag, globs):
    """
    Mimics ./InnerOptimizationCode/BelObjectiveUncondGradNAGBGP.m

    Computes the gradient of the bellman equation objective with
    respect to c_1(1), c_1(2) and c_2(1).  Substitues out for the rest
    of the variables using their respective gradients.

    z is the vector containing these three variables

    TODO: Fix this docstring to be in numpydoc format

    Notes
    -----
    We pass global variables as well (because NAG did not allow us
    to pass user defined information) V is a struct containing the
    interpolation information.  Vcoef has coefficents for the value
    function. R, x and s_ are the previous period states.  Par is a
    struct containig all the relevent parameters.  We return the
    gradient grad.  user and iflag our variables that nag requires
    but we don't use.

    the parameter n is just the number of equations. This can simply
    be inferred from z.size.

    The user and iflag parameters are simply NAG things that aren't
    used
    """
    # TODO: Put these in as args, but for now I will have them be like
    #       this.
    V = globs.V
    Vcoef = globs.Vcoef
    r = globs.r
    x = globs.x
    params = globs.params
    _s = globs._s

    psi = params.psi
    sigma = params.sigma
    beta = params.beta
    P = params.P
    th_1 = params.theta[0]
    th_2 = params.theta[1]
    g = params.g
    alpha = params.alpha
    n1 = params.n1
    n2 = params.n2

    frac = (r * P[_s, 0] * z[0] ** (-sigma) +
            r * P[_s, 1] * z[1] ** (-sigma) -
            P[_s, 0] * z[2] ** (-sigma)) / P[_s, 1]

    # frac must be positive
    if z.min() > 0 and frac > 0:
        c1_1 = z[0]
        c1_2 = z[1]
        c2_1 = z[2]

    # NOTE: If there is an error in passing arguments to these other
    #       functions it is because I compared the calls in
    #       BelObjectiveUncondGradNAGBGP.m to those in SteadyStateResiduals.m
    #       and the only difference was passing the x as the last arg
    #       of compute_X_prime instead of u2btild.

    c1, c2, gradc1, gradc2 = computeC2_2(c1_1, c1_2, c2_1, r, _s, P, sigma)
    Rprime, gradRprime = computeR(c1, c2, gradc1, gradc2, sigma)
    l1, gradl1, l2, gradl2 = computeL(c1, gradc1, c2, gradc2, Rprime,
                                 gradRprime, th_1, th_2, g, n1, n2)
    xprime, gradxprime = compute_X_prime(c1, gradc1, c2, gradc2, Rprime,
                                        gradRprime, l1, gradl1, l2, gradl2,
                                        P, sigma, psi, beta, _s, x)

    # TODO: Stopping on line 80 of the MatLab. Need to write funeval from
    #       compecon
