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
import numpy as np
import scipy.optimize as opt
import scipy.interpolate as interp
# from CompEcon import compeconpy
from steady.steady_state import steady_state_res, find_steady_state
from inneropt.inner_opt import uAlt
from set_params import DotDict

# Not sure what to do with the section DEFAULT PARAMETERS


def lookup(tabvals, x, endadj):
    """
    Just a function that matches lookup.c or lookup.m from compecon
    """
    import pandas as pd
    n = np.prod(x.shape)
    m = tabvals.size
    if endadj >= 2:
        m = m - np.where(tabvals == tabvals[-1])[0].size

    temp_series = pd.Series(np.append(tabvals[:m], x))
    temp_series.sort()
    ind = temp_series.index.values
    temp = np.where(ind > m - 1)[0]
    j = ind[temp] - m
    ind = (temp - range(1, n + 1)).reshape(x.shape)
    ind = ind[j]
    if endadj == 1 or endadj == 2:
        ind[ind == 0] = np.where(tabvals == tabvals[0])[0].size
    return ind


def fundefn(n, lb, ub, interp_type='spli', order=3):
    """
    Mimics the file ./CompEcon/cetools/fundefn.m

    Parameters
    ----------
    n: array-like, dtype=int
        The number of points along each dimension of the approximation

    lb: array-like, dtype=float
        The lower bounds for each variable

    ub: array-like, dtype=float
        The upper bounds for each variable

    interp_type: string
        The type of interpolation to use

    order: int, optional (default=3)
        The order of the spline along each dimension.

    Returns
    -------
    info_dict: DotDict
        A dictionary containing n, lb, ub, interp_type, order, and
        spine knots for each dimension.
    """
    parms = ['', '']
    parms[0] = [interp_type, [lb[0], ub[0]], n[0] - order + 1, order]
    parms[1] = [interp_type, [lb[1], ub[1]], n[1] - order + 1, order]
    d = len(lb)
    info_dict = DotDict()
    info_dict.d = d
    info_dict.n = n
    info_dict.a = lb
    info_dict.b = ub
    info_dict.interptype = interp_type

    params = np.array(['', ''], dtype=object)
    for i in xrange(d):
        this_space = np.linspace(parms[i][1][0], parms[i][1][1], parms[i][2])
        params[i] = np.array([this_space, parms[i][2], parms[i][3]])

    info_dict.params = params

    return info_dict


def funfitxy(info_dict, dom, vals):
    """
    Do all of what ./CompEcon/funfitxy.m does

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

    def splibas(breaks, evennum, k, x, order):
        """
        Mimics the file ./CompEcon/splibas.m
        """
        # Skipping 28-52 (just error checking)

        p = breaks.size
        m = x.size
        minorder = min(order)
        n = p + k - 1
        a = breaks[0]
        b = breaks[-1]
        augbreaks = np.append(a * np.ones(k - minorder), breaks)
        augbreaks = np.append(augbreaks, b * np.ones(k - minorder))

        ind = lookup(augbreaks, x, 3)

        bas = np.zeros((m, k - minorder + 1))
        bas[:, 0] = 1

        # Skipping 76-77
        for i in range(1, k - minorder + 1):
            for ii in range(i, 0, -1):
                b0 = augbreaks[ind + ii - i]
                b1 = augbreaks[ind + ii]
                temp = bas[:, ii - 1] / (b1 - b0)
                bas[:, ii] = (x - b0) * temp + bas[:, ii]
                bas[:, ii - 1] = (b1 - x) * temp

            iii = np.where((k - i) == order)[0]
            if not iii.size == 0:
                ii = iii[0]
                r = np.arange(m).reshape(m, 1)
                r = np.tile(r, (1, 4))
                c1 = np.arange(order[ii] - k, 1) - (order[ii] - minorder)
                c1 = np.tile(c1, (m, 1)).T + ind
                c = c1.T
                B = np.zeros((m, n))
                for row in range(m):
                    B[r[row, :], c[row, :]] = bas[row, :]

        return B

    def funbasx(info_dict, x, order, bformat):
        """
        Mimics ./CompEcon/funbasx.m
        """
        order = np.ascontiguousarray(order)
        d = max(info_dict['n'].shape)  # Number of dimensions

        # Expand order if needed
        if order.size != d:
            order = np.tile(order, (1, d))

        m = order.shape[0]
        # Skipping 62-64

        minorder = order + np.zeros((1, d))
        numbases = np.ones((1, d))
        B = DotDict()
        B.vals = ['', '']
        B.order = minorder
        B.format = bformat

        # Skipping 73, 75-79 checks if x and bformat are empty: they aren't

        # Skipping 81-91 (not doing tensor or direct method)

        B.format = 'direct'  # Line 93 changes B.format away from bformat?

        # Skipping 98-110 (not doing tensor)
        if B.format == 'direct':  # of course it is!
            for j in xrange(d):
                orderj = np.unique(order[:, j]) if m > 1 else order[0, j]
                orderj = np.ascontiguousarray(orderj)
                if orderj.size == 1:
                    breaks = info_dict.params[j][0]
                    evennum = info_dict.params[j][1]
                    k = info_dict.params[j][2]
                    xj = x[:, j]
                    B.vals[j] = splibas(breaks, evennum, k, xj, orderj)

    def funbconv(b, order, format):
        """
        Mimics the file ./Compecon/funbconv.m
        """
        d = b.order.shape[1]
        if format == 'expanded':
            order = np.zeros((1, d))
        else:
            order = b.order

        numbas, d1 = order.shape

        if format == 'expanded':
            B = DotDict()
            B.vals = ['' for i in range(numbas)]
            B.order = order
            B.format = format

            # Skipping 55-63. Not doing 'tensor'

            if b.format == 'direct':
                n = 1
                for j in range(d):
                    n *= b.vals[j].shape[1]
                bsize = np.array([1, n])  # NOTE: This is different than 66!
                for i in range(numbas):
                    B.vals[i] = b.vals[order[i, d - 1] - b.order[d - 1], d - 1]

                # Stopping for the night at 68. Some things might be fishy with
                # order and b.order. Check on that.


    m = vals.size   # Number of data

    if np.prod(info_dict['n']) > m:
        raise TypeError(' Number of data points must be bigger than prod(n)')

    if dom.shape[0] != m:
        raise ValueError('dom and vals must have the same number of data points')

    B = funbasx(info_dict, dom, 0, 'expanded')


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

    V = np.zeros(2, dtype=object)
    V[0] = fundefn([params.orderofappx_x, params.orderofappx_R],
                   [xMin, RMin], [xMax, RMax])
    V[1] = fundefn([params.orderofappx_x, params.orderofappx_R],
                   [xMin, RMin], [xMax, RMax])

    # V[0] = fundefn(params.ApproxMethod,[params.orderofappx_x, params.orderofappx_R],
    #         [xMin, RMin], [xMax, RMax])

    # V[1] = V[0]
    #We return the updated params and the functional space

    #It seems that we can just use our own grid that was built using linspace
    #One difference in the grid that is created by fundef is that it would appear
    #they account for degrees of freedom of some sort b/c the number of gridpoints
    #is n + 1 - k.  Don't think we need to worry about it much, but it's worth noting

    return params, V


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
