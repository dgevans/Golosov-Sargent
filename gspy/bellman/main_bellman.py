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

# Not sure what to do with the section DEFAULT PARAMETERS


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

    if params.flagSetxGrid == 1:
        xMin = params.xMin
        xMax = params.xMax
        disp('Msg: Using user defined grid on x')
        #^Will proceed using user defined gridpoints
    else:
        xMin = xSS - params.DeltaX
        xMax = xSS + params.DeltaX
        disp('Msg: Using default grid around SS')
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

    if params.flagSetRGrid == 1:
        RMin = params.RMin
        Rmax = params.Rmax
        disp('Proceeding with RGrid based on user inputs')
        #^Will proceed using user defined grid
    else:
        RMin = RSS - params.DeltaR
        RMax = RSS + params.DeltaR
        disp('Proceeding with default RGrid')
        #^Will proceed using default grid

    #Uniformly space the gridpoints
    RGrid = np.linspace(RMin, RMax, params.RGridSize, endpoint=True)
    params.RGrid = RGrid
    #Gives the total gridsize
    GridSize = params.xGridSize * params.RGridSize * params.sSize

    #Update params struct
    param.Gridsize = GridSize
    param.xMin = xMin
    params.xMax = xMax
    params.RMax = RMax
    params.RMin = RMin

    ##Define Functional Space
    #Priorly used the CompEcon Library routine 'fundefn' to create a functional
    #space which stores key settings for the basis polynomials, domand and nodes.
    #V(s) is he functional space for the value function given the discrete shock
    #Need to look up documentation on Compecon toolbox function fundefn

    V = np.zeros(2)
    V[0] = fundefn(params.ApproxMethod,[params.orderofappx_x, params.orderofappx_R],
            [xMin, RMin], [xMax, RMax])

    V[1] = V[0]
    #We return the updated params and the funcitonal space
    return params, V


def init_coef(params, V):
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
    n_size = params.xGridSize + params.RGridSize
    c0 = np.empty(params.sSize, n_size)
    V0 = np.empty(params.sSize, n_size)
    xInit_0 = np.empty(params.sSize, n_size, 7)

    p = params

    domain_ = np.empty((1, p.xGridSize * p.RGridSize, 2))
    for _s in range(params.sSize):
        n = 0
        if _s == 0:
            for xctr in xrange(params.xGridSize):
                for Rctr in xrange(p.RGridSize):
                    _x = xGrid[xctr]
                    _R = RGrid[Rctr]
                    domain_[_s, n, :] = [_x, _R]
                    #Initialize the guess for Stationary Policies
                    cRat = _R ** (-1. / p.sigma)

                    c1_1 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[0])\
                                / (p.n1 + cRat * p.n2)

                    c1_2 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[1])\
                                / (p.n1 + cRat * p.n2)

                    c2_1 = cRat * c1_1
                    #SteadyStateResiduals Routine
                    guess = np.array([c1_1, c1_2, c2_1])
                    [xSS, info, exitFlag, msg] = opt.fsolve(steady_state_res,
                                                       guess, full_output=1,
                                                        args=(_x, _R, p, _s))

                    [res, c1_, c2_, l1_, l2_] = steady_state_res(xSS, _x, _R,
                                                                 p, _s)

                    #Present Discounted value for Stationary policies
                    V0[_s, n] = (p.alpha_1 * uAlt(c1_, l2_, p.psi, p.sigma) +
                                p.alpha_2 * uAlt(c1_, l2_, p.psi, p.sigma)) * \
                                p.P[_s, :].T / (1 - p.beta)  # TODO: Check this

                    xInit_0[_s, n, :] = [c1_[_s, 0], c2_[_s, 0], l1_[_s, 0],
                                         l2_[_s, 0], _x, _R, _x]

                    n += 1

                    #Then need to initialize the coeffs by a routine that is equivalent to funfitxy

            c0[_s, :] = 0  # TODO: Fix this after we figure out what funfitxy does
        else:
            c0[_s, :] = c0[_s - 1, :]
            V0[_s, :] = V0[_s - 1, :]
            xInit_0[_s, :] = xInit_0[_s - 1, :]
    #ends for _s in...

    pass


def main(params):
    #BUILD GRID
    [params, V] = build_grid(params)
    disp('Msg: Completed definition of functional space')

    #INITIALIZE THE COEFF
    disp('Msg: Initializing the Value Function...')
    [domain, c, PolicyRulesStore] = init_coef(params, V)
    disp('Msg: ... Completed')