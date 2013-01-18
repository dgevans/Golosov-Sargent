"""
Created Dec 22, 2012

Python code that matches MatLab code from directory ./SteadyStateCode
"""
from __future__ import division
import numpy as np
import scipy.optimize as opt
import scipy.linalg as la
from inner_opt import computeC2_2, computeR, computeL, compute_X_prime


def steady_state_res(x, u2bdiff, rr, params, s, from_fsolve=0):
    """
    Mimics the file ./SteadyStateCode/SteadyStateResiduals.m

    Notes
    -----
    from_fsolve is a bool that tells whether or not this function is
    being called from scipy.optimize.fsolve. Default value is 0.
    """
    #Initialize the Parameters
    p = params
    p.theta = np.array([p.theta_1, p.theta_2])  # Why do this?
    p.alpha = np.array([p.alpha_1, p.alpha_2])

    u2btild = u2bdiff
    r = rr
    _s = s
    n1 = p.n1
    n2 = p.n2
    ctol = p.ctol

    #Get the Policy Rules
    psi = p.psi
    beta = p.beta
    P = p.P
    th_1 = p.theta_1  # Getting it directly from params.
    th_2 = p.theta_2
    g = p.g
    alpha = p.alpha
    sigma = p.sigma

    frac = (r * P[_s, 0] * x[0] ** (-sigma) + r * P[_s, 1] * x[1] ** (-sigma)\
            - P[_s, 0] * x[2] ** (-sigma)) / P[_s, 1]

    if x.min() > 0 and frac > 0:
        c1_1 = x[0]
        c1_2 = x[1]
        c2_1 = x[2]

        #Compute components from unconstrained guess
        # TODO: check if passing _s = 0 at first screws things up. I think
            # it is ok because as far as I can tell s_ is only used as index
            # in MatLab code. Ask David/Anmol to verify.
        c1, c2, gradc1, gradc2 = computeC2_2(c1_1, c1_2, c2_1, r, _s, P, sigma)
        Rprime, gradRprime = computeR(c1, c2, gradc1, gradc2, sigma)
        l1, gradl1, l2, gradl2 = computeL(c1, gradc1, c2, gradc2, Rprime,
                                     gradRprime, th_1, th_2, g, n1, n2)
        xprime, gradxprime = compute_X_prime(c1, gradc1, c2, gradc2, Rprime,
                                            gradRprime, l1, gradl1, l2, gradl2,
                                            P, sigma, psi, beta, _s, u2btild)

        #State Next Period
        xprime = xprime[0, :]
        Rprime = Rprime[0, :]

        res = np.empty(3)
        res[0] = xprime[0] - xprime[1]
        res[1] = Rprime[0] - Rprime[1]
        res[2] = xprime[0] - u2btild

        c1 = c1[0, :]
        c2 = c2[0, :]
        l1 = l1[0, :]
        l2 = l2[0, :]

        if np.concatenate([l1, l2]).flatten().max() > 1:
            res = np.abs(x) + 100

        if not np.isreal(res).all():
            res = np.abs(res) + 100

    else:
        res = np.abs(x) + 100

    # TODO: Python expects only res to be returned. Currently it gives some
        # error about a mismatch in the input and output shape of this func
        # when it is called by opt.fsolve in find_steady_state. I need to
        # figure out how to get c1, c2, l1, l2 out of this.

    if from_fsolve == 1:
        return res
    else:

        return res, c1, c2, l1, l2


def ss_residuals(X, Params):
    '''
    Mimics the  file ./SteadyStateCode/SSResiduals.m
    '''
    psi = Params.psi
    beta = Params.beta
    alpha_1 = Params.alpha_1
    alpha_2 = Params.alpha_2
    theta_1 = Params.theta_1
    theta_2 = Params.theta_2
    n1 = Params.n1
    n2 = Params.n2
    g = Params.g
    sigma = Params.sigma
    P = Params.P
    P = P[0, :]
    # Palt = np.fliplr(P)
    Palt = P[::-1]

    c1 = X[0:2]
    c2 = X[2:4]
    l1 = X[4:6]
    l2 = X[6:8]
    R = X[8]
    x = X[9]
    mu = X[10]
    lamb = X[11]
    phi = X[12:14]
    xi = X[14:16]
    rho = X[16:18]

    uc1 = psi * (c1 ** (-sigma))
    uc2 = psi * (c2 ** (-sigma))
    ul1 = -(1.0 - psi) / (1.0 - l1)
    ul2 = -(1.0 - psi) / (1.0 - l2)
    uc2Alt = uc2[::-1]

    Euc2 = np.dot(P, uc2)
    Euc1 = np.dot(P, uc1)

    #Initialize res matrix to fill
    res = np.empty(18)

    res[0:2] = (x * uc2) / (beta * Euc2) \
    - (uc2 * (c2 - c1)) \
    - x + R * l1 * ul1 \
    - l2 * ul2

    res[2:4] = theta_2 * R * (1.0 - l2) \
    - theta_1 * (1.0 - l1)

    res[4:6] = n1 * theta_1 * l1 \
    + n2 * theta_2 * l2 \
    - n1 * c1 \
    - n2 * c2 \
    - g

    res[6:8] = uc1 * R \
    - uc2

    res[8:10] = alpha_1 * P * uc1 \
    + mu * P * uc2 \
    - n1 * xi \
    - sigma * R * rho * uc1 / c1

    res[10:12] = alpha_2 * P * uc2\
    + mu * P * uc2 *( (sigma * x * (P *uc2 - Euc2)) / (beta * c2 * Euc2**2) \
    + sigma - 1.0 - sigma * c1 / c2 ) \
    + (sigma * x * mu * Palt * uc2Alt * P * uc2) / (beta * c2 * Euc2**2) \
    - n2 * xi + sigma * rho * uc2 / c2

    res[12:14] = alpha_1 * P * ul1\
    - mu * (1.0 - psi) * R * P / ( (1.0 - l1)**2 ) \
    + (phi + n1 * xi) * theta_1

    res[14:16] = alpha_2 * P * ul2 \
    + mu * (1.0 - psi) * P / ( (1.0 - l2)**2 ) \
    + (n2 * xi - R * phi) * theta_2

    res[16:18] = -beta * lamb * P * Euc1 \
    + mu * P * ul1 * l1 \
    + lamb * P *uc1 \
    + theta_2 * phi * (1.0 - l2) \
    + rho * uc1

    #Appeared to be an alternate method to build res, didn't translate.  Put here in case
    #we need it in the future

    ##% res(1:2) =x*uc2/(beta*Euc2)-( psi*uc2.*(c2-c1)+x+(1-psi)*( R*l1./(1-l1) - l2./(1-l2) ) )
    ##%
    ##% res(3:4) = theta_1*(1-l1)-theta_2*(1-l2)*R
    ##%
    ##% res(5:6) = theta_1*l1+theta_2*l2-c1-c2-g
    ##%
    ##% res(7:8) = psi*alpha_1*P.*uc1-psi*mu*P.*uc2-xi
    ##%
    ##% res(9:10) = alpha_2*psi*P.*uc2-mu*P.*( psi*(sigma*c1-(sigma-1)*c2)./(c2.^(sigma+1)) ...
    ##%     -sigma*x*(Euc2-P.*uc2)./(c2.^(sigma+1)*Euc2^2) )+sigma*x*mu*Palt.*P.*uc2Alt./(beta*c2.^(sigma+1).*Euc2^2)
    ##%
    ##% res(11:12) = theta_1*(xi+phi)-alpha_1*(1-psi)*P./(1-l1) -(1-psi)*R*mu*P./( (1-l1).^2 )
    ##%
    ##% res(13:14) = theta_2*(xi-R*phi) - alpha_2*(1-psi)*P./(1-l2) + mu.*P./( (1-l2).^2 )
    ##%
    ##% res(15:16) = beta*P*lambda*Euc1 - mu*(1-psi).*P.*l1./(1-l1)...
    ##%     -lambda*P.*uc1+phi.*(1-l2)

    return res


def find_steady_state(x0, R0, Params):
    '''
    Mimics file ./SteadyStateCode/findSteadyState.m
    '''
    p = Params
    cRat = R0 ** (-1.0 / p.sigma)

    c1_1 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[0]) / \
           (p.n1 + cRat * p.n2)

    c1_2 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[1]) / \
           (p.n1 + cRat * p.n2)

    c2_1 = cRat * c1_1

    guess = np.array([c1_1, c1_2, c2_1])
    xSS = opt.fsolve(steady_state_res, guess, args=(x0, R0, p, 1, True))
    # This syntax error is OK.

    [res, c1_, c2_, l1_, l2_] = steady_state_res(xSS, x0, R0, p, 1.0)

    X_temp = np.row_stack([c1_, c2_, l1_, l2_]).flatten()
    X = np.append(X_temp, [R0, x0])

    f = lambda Mult: findMultipliers(X, Mult, p)
    I = np.eye(8)

    b = - f(np.zeros((1, 8))).T

    A = np.zeros((8, 8))
    for i in xrange(8):
        A[:, i] = f(I[i, :]).T + b

    Mult = la.inv(A).dot(b)

    X = np.append(X, Mult)

    PolicyRule = opt.fsolve(ss_residuals, X, p)

    x = PolicyRule[9]
    R = PolicyRule[8]

    return x, R, PolicyRule


def findMultipliers(X, Mult, Params):
    X = np.append(X, Mult)
    resTemp = ss_residuals(X, Params)
    res = resTemp[10:18]
    return res
