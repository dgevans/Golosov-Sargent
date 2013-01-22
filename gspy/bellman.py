"""
This is the equivalent to MainBellman.m from the original code
This is the main file that computes the value function via
time iteration for the parameters passed in the structure parameters

Notation:
    x = u_2 btild
    R = u_2/u_1
"""

from __future__ import division
from itertools import product
import time
import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
from scipy.io import savemat
from compeconpy import fundefn, funfitxy
from steady_state import steady_state_res, find_steady_state
from inner_opt import uAlt, check_grad
from set_params import DotDict

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
    params.GridSize = GridSize
    params.xMin = xMin
    params.xMax = xMax
    params.RMax = RMax
    params.RMin = RMin

    ##Define Functional Space
    #Priorly used the CompEcon Library routine 'fundefn' to create a functional
    #space which stores key settings for the basis polynomials, domand and nodes.
    #V(s) is he functional space for the value function given the discrete shock
    #Need to look up documentation on Compecon toolbox function fundefn

    V = np.zeros(2, dtype=object)
    V[0] = fundefn([params.orderofappx_x, params.orderofappx_R],
                   [xMin, RMin], [xMax, RMax])
    V[1] = fundefn([params.orderofappx_x, params.orderofappx_R],
                   [xMin, RMin], [xMax, RMax])

    #It seems that we can just use our own grid that was built using linspace
    #One difference in the grid that is created by fundef is that it would appear
    #they account for degrees of freedom of some sort b/c the number of gridpoints
    #is n + 1 - k.  Don't think we need to worry about it much, but it's worth noting

    return params, V


def init_coef(params, info_dict):
    """
    Mimics ./BellmanEquationCode/InitializeCoeff.m

    INITIALIZE THE COEFF
    This section uses the stationary policies to initialze the value
    functions. The block prdoduces three outcomes
    1. domain which is the vectorized domain
    2. PolicyRulesStore : This serves as a matrix of guess for optimal
    policies at the interpolation nodes
    3. c0 : initial coeffecients
    """
    p = DotDict(params)
    xGrid = p.xGrid
    RGrid = p.RGrid
    gTrue = p.g
    p.g = gTrue.mean() * np.ones((2, 1))

    #Need to initialize arrays before we fill them.
    n_size = p.xGridSize * p.RGridSize
    n_coefs = (p.xGridSize - 1) * (p.RGridSize - 1)
    c0 = np.zeros((p.sSize, n_coefs))
    V0 = np.zeros((p.sSize, n_size))
    xInit_0 = np.zeros((2, p.xGridSize * p.RGridSize, 7))



    _domain = np.array(list(product(p.xGrid, p.RGrid)))
    for _s in range(p.sSize):
        n = 0
        if _s == 0:
            for xctr in xrange(p.xGridSize):
                for Rctr in xrange(p.RGridSize):
                    _x = xGrid[xctr]
                    _R = RGrid[Rctr]

                    # NOTE: _domain is product(xGrid, RGrid). See above loop
                    # _domain[_s, n, :] = [_x, _R]

                    # Initialize the guess for Stationary Policies
                    cRat = _R ** (-1. / p.sigma)

                    c1_1 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[0])\
                                / (p.n1 + cRat * p.n2)

                    c1_2 = (0.8 * (p.n1 * p.theta_1 + p.n2 * p.theta_2) - p.g[1])\
                                / (p.n1 + cRat * p.n2)

                    c2_1 = cRat * c1_1

                    guess = np.array([c1_1, c1_2, c2_1]).flatten()
                    [xSS, info, exitFlag, msg] = opt.fsolve(steady_state_res,
                                                       guess, full_output=1,
                                                        args=(_x, _R, p, _s, True),
                                                        xtol=1.0e-13)

                    [res, c1_, c2_, l1_, l2_] = steady_state_res(xSS, _x, _R,
                                                                 p, _s)

                    # Present Discounted value for Stationary policies
                    V0[_s, n] = (p.alpha_1 * uAlt(c1_, l1_, p.psi, p.sigma) +
                                p.alpha_2 * uAlt(c2_, l2_, p.psi, p.sigma)).dot( \
                                p.P[_s, :].T) / (1 - p.beta)

                    xInit_0[_s, n, :] = [c1_[_s], c2_[_s], l1_[_s],
                                         l2_[_s], _x, _R, _x]

                    n += 1

            c0[_s, :], B = funfitxy(info_dict[_s], _domain, V0[_s, :])
        else:
            c0[_s, :] = c0[_s - 1, :]
            V0[_s, :] = V0[_s - 1, :]
            xInit_0[_s, :, :] = xInit_0[_s - 1, :, :]

    domain_1 = np.column_stack((_domain, np.ones((n_size, 1))))
    domain_2 = np.column_stack((_domain, np.ones((n_size, 1)) * 2))
    domain = np.concatenate((domain_1, domain_2), axis=0)

    c = c0
    data = {'c': c}
    savemat(p.datapath + 'c1.mat', data)

    x_shape = xInit_0.shape
    policy_rules_store = np.empty((x_shape[0] * x_shape[1], x_shape[2] * 2)).T

    for i in range(x_shape[2]):
        cols = range(2 * i, 2 * i + 2)
        policy_rules_store[cols, :400] = xInit_0[0, :, i].squeeze()
        policy_rules_store[cols, 400:] = xInit_0[1, :, i].squeeze()

    policy_rules_store = policy_rules_store.T

    # NOTE: this is just here for comparison with the MatLab objects
    # data = {'my_c0': c0, 'my_v0': V0, 'my_xinit': xInit_0,
    #         'my_domain': domain, 'my_policyrules': policy_rules_store}
    # savemat('/users/spencerlyon2/documents/research/golosov-sargent/' + \
    #           'gspy/data/debugging/init_coefs.mat', data)

    # NOTE: All objects in data are within 1e-13 of corresponding MatLab objs

    return [domain, c, policy_rules_store]


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
    params, info_dict = build_grid(params)
    print('Msg: Completed definition of functional space')

    #INITIALIZE THE COEFF
    print('Msg: Initializing the Value Function...')
    [domain, c, policy_rules_store] = init_coef(params.copy(), info_dict)
    print('Msg: ... Completed')

    ctest = c  # NOTE: THis is just so the object c is available in pdb

    #Iterate on the Value Function
    #This block iterates on the bellman equation
    x_slice = domain[:, 0]
    R_slice = domain[:, 1]
    s_slice = domain[:, 2]
    GridSize = int(params.GridSize)

    #Initialize The Sup Norm Error matrix and variables needed in loop
    errorinsupnorm = np.ones(params.Niter)
    policy_rules_old = np.zeros(policy_rules_store.shape)
    xInit = np.zeros((GridSize / 2, policy_rules_store.shape[1]))  # Double Check Size
    vnew = np.zeros((GridSize / 2, 1))

    #Begin the for loops
    for iter in xrange(1, params.Niter):
        #Record Start Time.  Total time will be starttime-endtime
        starttime = time.time()

        #Clear Records of arrays that store index of failure
        #^We can add if we need to

        #Initialize the initial guess for the policy rules that the inneropt
        #will solve
        policy_rules_old = policy_rules_store

        for ctr in xrange(int(GridSize / 2)):
            #Here they use a parfor loop invoking parallel type for loops
            #Will make this parallel when we speed up program
            x = x_slice[ctr]
            R = R_slice[ctr]
            s = int(s_slice[ctr] - 1)  # TODO: These will be indices so I need to -1

            #Initalize Guess for the inneropt
            xInit = policy_rules_store[ctr, :]

            #Inner Optimization
            policyrules, v_new = check_grad(x, R, s, c, info_dict, xInit, params)
            vnew[ctr, 0] = v_new

            #Update Policy rules
            policy_rules_store[ctr, :] = policyrules
