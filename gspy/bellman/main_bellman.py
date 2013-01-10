"""
This is the equivalent to MainBellman.m from the original code
This is the main file that computes the value function via
time iteration for the parameters passed in the structure parameters

Notation:
    x = u_2 btild
    R = u_2/u_1
"""

from __future__ import division
import numpy as np
import scipy.optimize as opt
# from CompEcon import compeconpy
from steady.steady_state import steady_state_res, find_steady_state
from inneropt.inner_opt import uAlt
from itertools import product
import scipy.interpolate as interp

# Not sure what to do with the section DEFAULT PARAMETERS


def funfitxy(info_dict, dom, vals):
    """
    Do all of what funfitxy does

    Parameters
    ----------
    dom: array-like, dtype=float
        The domain over which the interpolating function is to be
        created: the matrix of (x, y) pairs if z = f(x, y)

    vals: array-like, dtype=float
        The range or points on the interpolating surface: z points if
        z = f(x, y)

    Returns
    -------
    i_dont_know_yet:
    """
    m = vals.size   # Number of data points


#Build Grid
def build_grid(params):
    '''
    This is the function that executes the equivalent of BuildGrid.m.
    This function defines the grid and defines the value function.  There
    are two alternatives.  First, the user could input either the x or Rgrid.
    This should supercede any other option.  Otherwise we use the SS
    computation and use efault DeltaX,DeltaR parameters to set the deviation
    from the SS.

    params.flagSetxGrid sets the flag for either using the default grid (Value = 0)
     or using the user defined grid (Value = 1)

    params.flagsetRdgrid sets the flag for either using the default grid (Value = 0)
    or using the user defined grid (Value = 1)
    '''
    #???Do we want to name all of the elements in params like params[0] = ___
    #How are we importing it?

    #Find the SS
    [xSS, RSS, NA] = find_steady_state(0, 3, params)

    #Check the flags and define Grid Endpoints
    if 'flagSetxGrid' in params:
        flagSetxGrid = params.flagSetxGrid
    else:
        flagSetxGrid = 0

    if flagSetxGrid == 1:
        xMin = params.xMin
        xMax = params.xMax
        print('Msg: Using user defined grid on x')
        #^Will proceed using user defined gridpoints
    else:
        xMin = xSS - params.DeltaX
        xMax = xSS + params.DeltaX
        print('Msg: Using default grid around SS')
        #^Will proceed using the default gridpoints
    #Uniformly space the gridpoints
    xGrid = np.linspace(xMin, xMax, params.xGridSize, endpoint=True)

    #Update the params struct
    params.xGrid = xGrid
    params.xLL = xMin
    params.xUL = xMax

    #Check the flags and define Grid Endpoints
    if 'flagSetRGrid' in params:
        flagSetRGrid = params.flagSetRGrid
    else:
        flagSetRGrid = 0

    if flagSetRGrid == 1:
        RMin = params.RMin
        RMax = params.RMax
        print('Proceeding with RGrid based on user inputs')
        #^Will proceed using user defined grid
    else:
        RMin = RSS - params.DeltaR
        RMax = RSS + params.DeltaR
        print('Proceeding with default RGrid')
        #^Will proceed using default grid

    #Uniformly space the gridpoints
    RGrid = np.linspace(RMin, RMax, params.RGridSize, endpoint=True)
    params.RGrid = RGrid
    #Gives the total gridsize
    GridSize = params.xGridSize * params.RGridSize * params.sSize

    #Update params
    params.Gridsize = GridSize
    params.xMin = xMin
    params.xMax = xMax
    params.RMax = RMax
    params.RMin = RMin

    ##Define Functional Space
    #Priorly used the CompEcon Library routine 'fundefn' to create a functional
    #space which stores key settings for the basis polynomials, domand and nodes.
    #V(s) is he functional space for the value function given the discrete shock
    #Need to look up documentation on Compecon toolbox function fundefn

    V = np.zeros(2)
    # V[0] = fundefn(params.ApproxMethod,[params.orderofappx_x, params.orderofappx_R],
    #         [xMin, RMin], [xMax, RMax])

    # V[1] = V[0]
    #We return the updated params and the functional space

    #It seems that we can just use our own grid that was built using linspace
    #One difference in the grid that is created by fundef is that it would appear
    #they account for degrees of freedom of some sort b/c the number of gridpoints
    #is n + 1 - k.  Don't think we need to worry about it much, but it's worth noting

    return params


def init_coef(params):
    #This function is the equivalent of InitializeCoeff.m
    '''
    INITIALIZE THE COEFF
    This section uses the stationary policies to initialze the value
    functions. The block prdoduces three outcomes
    1. domain which is the vectorized domain
    2. PolicyRulesStore : This serves as a matrix of guess for optimal
    policies at the interpolation nodes
    3. c0 : initial coeffecients
    '''
    xGrid = params.xGrid
    RGrid = params.RGrid
    gTrue = params.g
    params.g = gTrue.mean() * np.ones((2, 1))
    #Need to initialize arrays before we fill them.
    #TODO: Determine size of c0 and V0 after funfitxy is written
    n_size = params.xGridSize * params.RGridSize
    c0 = np.zeros((params.sSize, n_size))
    V0 = np.zeros((params.sSize, n_size))
    xInit_0 = np.zeros((2, params.xGridSize * params.RGridSize, 7))


    #Testing scipy.interpolate.RectBivariate.... 
    #Reshaping the matrices so that they pass in to function
    xx = np.linspace(params.xMin, params.xMax, 19, endpoint=True)
    xGrid = xx

    rr = np.linspace(params.RMin, params.RMax, 19, endpoint=True)
    RGrid = rr
    
    V0 = np.zeros((xGrid.size,RGrid.size))

    p = params

    _domain = np.array(list(product(params.xGrid, params.RGrid)))
    for _s in range(params.sSize):
        n = 0
        if _s == 0:
            for xctr in xrange(xGrid.size):
                for Rctr in xrange(RGrid.size):
                    _x = xGrid[xctr]
                    _R = RGrid[Rctr]

                    # NOTE: Just becomes product(xGrid, RGrid). See above loop
                    # _domain[_s, n, :] = [_x, _R]
                    #Initialize the guess for Stationary Policies
                    cRat = _R ** (-1. / p.sigma)

                    c1_1 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[0])\
                                / (p.n1 + cRat * p.n2)

                    c1_2 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[1])\
                                / (p.n1 + cRat * p.n2)

                    c2_1 = cRat * c1_1

                    guess = np.array([c1_1, c1_2, c2_1]).flatten()
                    [xSS, info, exitFlag, msg] = opt.fsolve(steady_state_res,
                                                       guess, full_output=1,
                                                        args=(_x, _R, p, _s, True))

                    [res, c1_, c2_, l1_, l2_] = steady_state_res(xSS, _x, _R,
                                                                 p, _s)

                    #Present Discounted value for Stationary policies
                    #change back to [_s,n]                   
                    V0[xctr, Rctr] = (p.alpha_1 * uAlt(c1_, l1_, p.psi, p.sigma) +
                                p.alpha_2 * uAlt(c2_, l2_, p.psi, p.sigma)).dot( \
                                p.P[_s, :].T) / (1 - p.beta)  # TODO: Check this

                    xInit_0[_s, n, :] = [c1_[_s], c2_[_s], l1_[_s],
                                         l2_[_s], _x, _R, _x]

                    n += 1

                    #Then need to initialize the coeffs by a routine that is equivalent to funfitxy
            coefs = interp.RectBivariateSpline(xGrid, RGrid, V0).get_coeffs()
            c0[_s, :] = 0  # TODO: Fix this after we figure out what funfitxy does
        else:
            c0[_s, :] = c0[_s - 1, :]
            V0[_s, :] = V0[_s - 1, :]
            xInit_0[_s, :, :] = xInit_0[_s - 1, :, :]

    domain_1 = np.column_stack((_domain, np.ones((n_size, 1))))
    domain_2 = np.column_stack((_domain, np.ones((n_size, 1)) * 2))
    domain = np.concatenate((domain_1, domain_2), axis=0)

    pass


def main(params):
    """
    This is the main file computes the value function via time iteration for
    the parameters passsed in the structure Para.

    NOTATION:
    --------
    x = u_2 btild
    R = u_2/u_1
    """
    #BUILD GRID
    params = build_grid(params)
    print('Msg: Completed definition of functional space')

    #INITIALIZE THE COEFF
    print('Msg: Initializing the Value Function...')
    [domain, c, PolicyRulesStore] = init_coef(params)
    print('Msg: ... Completed')
