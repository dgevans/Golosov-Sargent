"""
Created Jan 15, 2013

Author: Spencer Lyon and Chase Coleman

Translation of the compecon code used to solve the Golosov-Sargent economy
"""
import scipy.linalg as la
import numpy as np
import pandas as pd


class DotDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, in_dct=None):
        """
        A dictionary that supports dot notation as well as standard
        dictionary notation.

        Examples
        ---------
        >>> d = DotDict() or d = DotDict({'val1':'first'})
        >>> d.val2 = 'second'  # use dot notation to set attributes
        >>> d['val2'] = 'second'  # use dict notation to set attributes
        >>> d.val2  # use dot notation get attributes
        >>> d['val2']  # use dict notation get attributes
        """
        dct = in_dct if in_dct else dict()
        for key, value in dct.items():
            if hasattr(value, 'keys'):
                value = DotDict(value)
            self[key] = value


def lookup(tabvals, x, endadj):
    """
    Just a function that matches lookup.c or lookup.m from compecon

    This is correct, at least for init_coefs case
    """
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


def dprod(a, b):
    """
    Mimics the file ./CompEcon/dprod.m
    """
    ra, ca = a.shape
    rb, cb = b.shape
    if ra != rb:
        raise ValueError('a and b must have the same number of rows')
    c = np.zeros((ra, ca * cb))
    for i in range(ra):
        c[i, :] = np.kron(a[i, :], b[i, :])
    return c


def fundefn(n, lb, ub, interp_type='spli', order=3):
    """
    Mimics the file ./CompEcon/cetools/fundefn.m and is correct

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


def spdiags(B, d, a3, a4=None):
    """
    Mimics the file ./CompEcon/spdiags.m
    """
    moda = 1 if a4 == None else 0
    if moda:
        A = a3
    else:
        A = np.zeros((a3, a4))

    p = len(d)

    i, j = np.nonzero(A)
    a = A[i, j]
    a = np.column_stack((i, j, a))
    m, n = A.shape

    if moda:  # Skipping 96-105
        pass

    else:
        leng = np.zeros(p + 1)
        for k in xrange(p):
            leng[k + 1] = leng[k] + len(range(max(1, 1 - d[k]),
                                              min(m, n - d[k]) + 1))

        leng = np.array(leng, dtype=int)  # cast as int for indexing later

        a = np.zeros((leng[p], 3))
        for k in range(p):
            i = np.c_[range(max(1, 1 - d[k]), min(m, n - d[k]) + 1)] - 1
            a[np.arange(leng[k], leng[k + 1]), :] = np.column_stack((
                                                i,
                                                i + d[k],
                                                B[i + (m >= n) * d[k], k]))

        res1 = np.zeros((m, n))

        r = np.array(a[:, 0], dtype=int)
        c = np.array(a[:, 1], dtype=int)
        val = a[:, 2]

        if len(r.shape) == len(c.shape) == len(val.shape) == 1:
            res1[r, c] = val
        else:
            raise ValueError("Haven't implemented this yet. Come fix it")

    return res1


def splidop(breaks, evennum, k, order):
    """
    Mimics the file ./CompEcon/splidop.m
    """
    n = breaks.size + k - 1
    kk = max(k - 1, k - order - 1)
    a = breaks[0]
    b = breaks[-1]
    augbreaks = np.append(a + np.zeros(kk), breaks)
    augbreaks = np.append(augbreaks, b + np.zeros(kk))

    D = ['' for i in range(abs(order))]

    if order > 0:  # Doing derivative
        temp = k / (augbreaks[k: n + k - 1] - augbreaks[:n - 1])
        D[0] = spdiags(np.c_[-temp, temp], np.arange(2), n - 1, n)

        for i in range(1, order):  # TODO: This is probably wrong, but it isn't called
            temp = (k + 1 - i) / (augbreaks[k: + k - i] - augbreaks[i:n - 1])
            D[i] = spdiags(np.c_[-temp, temp], np.arange(2), n - i, n - i + 1)

    return D


def splibas(breaks, evennum, k, x, order):
    """
    Mimics the file ./CompEcon/splibas.m

    This is correct (at least for the init_coef case)
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

    if order.max() > 0:
        D = splidop(breaks, evennum, k, order.max())
    elif minorder < 0:
        I = splidop(breaks, evennum, k, minorder)

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
            r = np.tile(r, (1, k - order[ii] + 1))
            c1 = np.arange(order[iii] - k, 1) - (order[iii] - minorder)
            c1 = np.tile(c1, (m, 1)).T + ind
            c = c1.T
            B = np.zeros((m, n - order[iii]))
            for row in range(m):
                B[r[row, :], c[row, :]] = bas[row, :]

            if order[iii] > 0:
                B = B.dot(D[iii])

    return B


def funbasx(info_dict, x, order, bformat=None):
    """
    Mimics ./CompEcon/funbasx.m
    """
    order = np.ascontiguousarray(order)
    d = len(info_dict['n'])  # Number of dimensions
    x = np.atleast_2d(x)  # make sure x is correct shape for line xj = below

    # Expand order if needed
    if order.size != d:
        order = np.tile(order, (1, d))

    if len(order.shape) == 1:
        order = np.atleast_2d(order)

    if bformat == None:
        bformat = []

    m = order.shape[0]
    # Skipping 62-64

    minorder = order + np.zeros((1, d))
    numbases = np.ones((1, d))
    B = DotDict()
    B.vals = ['' for i in max(numbases)]
    B.order = minorder
    B.format = bformat

    # Skipping 73, checks if x is empty: it isn't

    if len(bformat) == 0:  # Skipped line 76
        bformat = 'direct'

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
    if bformat == 'expanded':
        B = funbconv(B, order, 'expanded')

    return B


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
            the_vals = np.asarray(b.vals)
            n = 1
            for j in range(d):
                n *= b.vals[j].shape[1]

            bsize = np.array([1, n])  # NOTE: This is different than 66!
            for i in range(numbas):
                if the_vals.shape[0] == d:
                    B.vals[i] = the_vals[d - 1]
                    for j in xrange(d - 1, 0, -1):
                        B.vals[i] = dprod(B.vals[i],
                                          the_vals[j - 1])
                else:
                    raise ValueError('I have not implemented this yet. \
                                     See line 68 for more info')

    return B

            # Stopping for the night at 68. Some things might be fishy with
            # order and b.order. Check on that.


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

    m = vals.size   # Number of data

    if np.prod(info_dict['n']) > m:
        raise TypeError(' Number of data points must be bigger than prod(n)')

    if dom.shape[0] != m:
        raise ValueError('dom and vals must have the same number of data points')

    B = funbasx(info_dict, dom, 0, 'expanded')
    c = la.lstsq(B.vals[0], vals)[0]
    return c, B


def funeval2(c, B, order):
    """
    Mimics the file ./CompEcon/funeval2.m

    # TODO: Fill in docs
    """
    order = np.atleast_2d(order)
    kk, d = order.shape

    # NOTE: This isn't generally true, but it works here.
    order2 = order + 1
    f = np.zeros((np.atleast_2d(B.vals[0]).shape[0],
                 np.atleast_2d(c).shape[0],
                 kk))

    for i in range(kk):
        # Putting code for cdprodx.m here
        b = B.vals
        ind = order2[i, :]
        d = ind.size
        a = b[ind[d - 1] - 1]
        for i in range(d - 2, -1, -1):
            a2 = dprod(b[ind[i] - 1], a)
        f[:, :, i] = a2.dot(c)

    return f


def funeval(c, info_dict, B, order):
    """
    Mimics the file ./CompEcon/funeval.m

    Parameters
    ----------
    c: array-like, dtype=float
        The matrix of coefficients for the interpoland

    info_dict: dictionary
        The python dictionary describing the functional space and the
        interpoland. Contains important information such as the number
        of dimensions, the number of coefficients in each dimension,
        the type of interpoland (spline, chebyshev, ect.) and the break
        points in each dimension.

    B:array-like, dtype=float
        A basis structure or an mxd matrix or a 1xd array of vectors.

    order: array-like, dtype=int
        An array describing the the order of the differential operator
        along each dimension of the interpoland. For example order =
        [0, 1, 1] corresponds to a mixed partial derivative with
        repect to the vairables in the 2nd and 3rd dimension.

    Notes
    -----
    When called from gspy, info_dict is is going to be V[0] or V[1].
    """
    # SKIPPING 107-115

    d = info_dict.d

    B2 = funbasx(info_dict, B, order)

    # SKIPPING 119-132

    # SKIPPING 139-141
    if B2.format == 'direct':
        y = funeval2(c, B2, order)

    return y.squeeze()  # NOTE: not returning B2 like they do.
