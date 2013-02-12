"""
Created Feb 11, 2013

Author: Spencer Lyon and Chase Coleman

Cython version of the compecon routines
"""
from __future__ import division
import numpy as np
import cython
import pandas as pd
import scipy.linalg as la

cimport cython
cimport numpy as np

DTYPE = np.float

ctypedef np.float_t DTYPE_t


@cython.boundscheck(False)
def lookup(np.ndarray[DTYPE_t, ndim=1] tabvals,
           np.ndarray[DTYPE_t, ndim=1] x,
           int endadj=0):
    """
    Just a function that matches lookup.c or lookup.m from compecon

    Performs a table lookup.

    Parameters
    ----------
    tabvals: array-like, dtype=float
        A sorted array of n values.

    x: array-like, dtype=float
        An array of values to find a position for.

    endadj: int, optional(default=0)
        Optional endpoint adjustment. Must be equal to 0, 1, 2, or 3.

    Returns
    -------
    ind: array-like, dtype=int
        An array of x.shape values where :math: `a_{i, j} = max k:
        x_{i, j} >= tabvals k`
    """
    cdef int n = np.prod(x.shape)
    cdef int m = tabvals.size
    if endadj >= 2:
        m = m - (tabvals == tabvals[-1]).sum()

    temp_series = pd.Series(np.append(tabvals[:m], x))
    temp_series.sort()

    cdef np.ndarray ind = temp_series.index.values
    cdef np.ndarray temp = np.where(ind > m - 1)[0]  # ind is 1d so return only rows.
    cdef np.ndarray j = ind[temp] - m
    ind = (temp - range(1, n + 1)).reshape(x.shape)
    ind[j] = ind.flatten()
    if endadj == 1 or endadj == 3:
        ind[ind == -1] = (tabvals == tabvals[0]).sum() - 1
    return ind


@cython.boundscheck(False)
def splidop(np.ndarray[DTYPE_t, ndim=1] breaks,
            int evennum,
            int k,
            int order):
    """
    Mimics the file ./CompEcon/splidop.m

    Computes the matrix operator that maps a vector of spline
    coefficients into the vector of coefficients of the derivative.
    Let g(x) be the spline evaluated at x, B be the coefficients,
    k be the order of each component in the spline. Then this function
    does the following:

        g(x) = np.dot(B ** k(x), c)
        g'(x) = np.dot(B ** {k-1}(x) * splidop(n, a, b, 1), c)
        g"(x) = np.dot(B ** {k-1}(x) * splidop(n, a, b, 2), c)

    Integrals are computed with order < 1

    Parameters
    ----------
    breaks: array-like, dtype=float
        The break points in the basis structure for the spline.

    evennum: int
        The number of unique items in breaks

    k: int
        The order of the spline. Generally, this is 3

    order: int
        The order of the derivative or anti-derivative to be evaluated.

    Notes
    -----
    You shouldn't call this function directly, rather you should use
    funeval which calls this function when needed.

    TODO
    ----
    spdiags: Note, I might be able to use the scipy.sparse version.
    """
    cdef int n = breaks.size + k - 1
    cdef int kk = max(k - 1, k - order - 1)
    cdef double a = breaks[0]
    cdef double b = breaks[-1]
    cdef np.ndarray augbreaks = np.append(a + np.zeros(kk), breaks)
    cdef np.ndarray augbreaks = np.append(augbreaks, b + np.zeros(kk))

    # No type declaration here b/c looping over lists is already optimized.
    D = ['' for i in range(abs(order))]

    cdef np.ndarray temp
    cdef int i

    if order > 0:  # Doing derivative
        temp = k / (augbreaks[k: n + k - 1] - augbreaks[:n - 1])
        D[0] = spdiags(np.c_[-temp, temp], np.arange(2), n - 1, n)

        for i in range(1, order):  # TODO: This is probably wrong, but it isn't called
            temp = (k + 1 - i) / (augbreaks[k: + k - i] - augbreaks[i:n - 1])
            D[i] = spdiags(np.c_[-temp, temp], np.arange(2), n - i, n - i + 1)

    return D


@cython.boundscheck(False)
def splibas(np.ndarray[DTYPE_t, ndim=1] breaks,
            int evennum,
            int k,
            np.ndarray[DTYPE_t, ndim=1] x,
            int order=0):

    cdef int p = breaks.size
    cdef int m = x.size
    cdef int minorder = min(order)
    cdef int n = p + k - 1
    cdef double a = breaks[0]
    cdef double b = breaks[-1]

    cdef np.ndarray augbreaks = np.append(a * np.ones(k - minorder), breaks)
    augbreaks = np.append(augbreaks, b * np.ones(k - minorder))

    cdef int ind = lookup(augbreaks, x, 3)

    cdef np.ndarray bas = np.zeros((m, k - minorder + 1))
    bas[:, 0] = 1

    if order.max() > 0:
        cdef np.ndarray D = splidop(breaks, evennum, k, order.max())
    elif minorder < 0:
        cdef np.ndarray I = splidop(breaks, evennum, k, minorder)

    cdef int i, ii, iii
    cdef np.ndarray temp, r, c1, c, B
    cdef dobule b0, b1

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
