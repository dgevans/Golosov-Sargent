"""
Created Jan 23, 2013

Author: Chase Coleman and Spencer Lyon

Cubic spline interpolation routines

TODO: Add re-scaling of nodes for more accuracy

TODO: Use algopy or numdifftools to add methods for calculating gradient

TODO: Cythonize or numaize this
"""
from __future__ import division
import numpy as np
from interpolate.cubicspline1d import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

__all__ = ['CubicSpline2d']


class CubicSpline2d(object):
    """
    TODO: Have Chase fill in the docstring that describes what the class
    does / references.
    """

    def __init__(self, lx, ly, ux, uy, nx, ny, alpha, beta):
        self.lx = lx
        self.ly = ly
        self.ux = ux
        self.uy = uy
        self.nx = nx
        self.ny = ny
        self.alpha = alpha
        self.beta = beta

        self.make_grid()

    def make_grid(self):
        """
        This function takes a(xy),b(xy),n(xy) as inputs and will return
        grid and has output
        a: the left endpoint of the interval you want to interpolate over.
        b: the right endpoint of the interval you want to interpolate over
        n: the number of equally spaced points you want to use


        grid: an evenly spaced interval in R1 with n+1 nodes
        h: the size of your steps

        x and y denote which variable it refers to for all of the above.
        """
        xgrid = np.linspace(self.lx, self.ux, self.nx + 1)
        ygrid = np.linspace(self.ly, self.uy, self.ny + 1)
        self.hx = xgrid[1] - xgrid[0]
        self.hy = ygrid[1] - ygrid[0]
        self.xgrid = xgrid
        self.ygrid = ygrid

        return xgrid, ygrid

    def _calcuv_i(self, x, x_or_y, i=None):
        """
        This function calculates the value of u_i and v_i by applying
        the function phi to the t   The u_i and v_iare the basis
        functions used to evaluate the spline

        x: Value of x (or y) we are using
        a: left endpoint
        h: step size
        i: which step you are on
        """
        # Get values of a or h depending on whether we are using x or y
        if x_or_y == 'x':
            a = self.lx
            h = self.hx
        elif x_or_y == 'y':
            a = self.ly
            h = self.hy
        else:
            raise ValueError("Variable x_or_y must have value equal to \
                             either 'x' or 'y'")

        #Determine t
        t = ((x - a) / h) + 2. - i

        #Calculate the value of the function phi(t)
        if np.abs(t) <= 1:
            out = 4. - 6. * np.abs(t) ** 2 + 3. * np.abs(t) ** 3
        elif 1 < np.abs(t) <= 2:
            out = (2 - np.abs(t)) ** 3
        else:
            out = 0

        return out

    def coefs(self, z):
        """
        This function calculates a cubic spline for a 2 dimensional case
        We proceed on the method explained in Habermann and Kindermann 2007
        It should run faster than typical cubic B-Spline interpolation and
        maintain or even improve the accuracy.

        Parameters
        ----------
        z: array-like, dtype=float
            Value of f(x, y) for the f that you are trying to
            approximate

        Returns
        -------
        c: coefficients s.t. s(x) = sum(c_ij * u_i(x) * v_j(x))
        """
        x = np.squeeze(self.xgrid)
        y = np.squeeze(self.ygrid)

        nx = self.nx
        ny = self.ny

        # Check simple errors.  x and y should be same size and x should be
        # evenly spaced
        if x.size * y.size != z.size:
            raise ValueError('z should have a size of  nx * ny. Cannot proceed')
        if x[2] - 2 * x[1] + x[0] > 1e-6:
            raise ValueError('x should be evenly spaced.  Cannot proceed')
        if y[2] - 2 * y[1] + y[0] > 1e-6:
            raise ValueError('y should be evenly spaced.  Cannot proceed')
        #Initialize the u matrices

        # Create a cstar matrix.  Then fill it with coefficients calculated by
        # running the 1d cubic spline over every possible y and letting x
        #  change
        cstar = np.zeros((nx + 3, ny + 1))

        for i in xrange(y.size):
            ctemp = calccubicspline1d(x, z[:, i], self.alpha, self.beta,
                                      self.hx)
            cstar[:, i] = ctemp.squeeze()

        # Create a true coeff mat.  Then fill it by running the 1d cubic
        # spline fixed at all of the points in y.

        c_mat = np.zeros((nx + 3, ny + 3))

        for i in xrange(nx + 3):
            ctemp = calccubicspline1d(y, cstar[i, :], self.alpha,
                                      self.beta, self.hy)
            c_mat[i, :] = ctemp

        self.c_mat = c_mat

        return c_mat

    def eval(self, point):
        """
        This function takes a value x and the coefficients obtained by calccubicspline and
        evaluates the spline at x.

        Inputs:
        x,y: The point at which you are evaluating the spline
        a(xy): The left endpoint
        h(xy): The step size obtained by gridbuild2
        coeffs: The coefficients obtained in calccubicspline2

        Outputs:
        z_spline: The estimated value of z.  (Value of the interpolation s(x).)
        """

        def find_z_point(x, y):
            tempsum = 0

            lx = math.floor((x - self.lx) / self.hx) + 1
            ly = math.floor((y - self.ly) / self.hy) + 1
            mx = np.min([lx + 3, self.nx + 3])
            my = np.min([ly + 3, self.ny + 3])

            #We want to calculate the values of the functions from lx to mx and ly to my
            #We set it up like this to make it easier to see range
            for row in xrange(int(mx - lx + 1)):
                for column in xrange(int(my - ly + 1)):

                    #we take coeffs_row,column for each of the possible values in U_x X V_y
                    #And multiply it by the basis function value.  Same as in 1d, but in 2d

                    tempsum += self.c_mat[row + lx - 1, column + ly - 1] \
                    * self._calcuv_i(x, 'x', row + lx) * self._calcuv_i(y, 'y', column + ly)

            z_spline = tempsum

            return z_spline

        if type(point) == np.ndarray:
            point = point if point.shape[0] == 2 else point.T
            z = np.zeros((point.shape[1], point.shape[1]))
            for ix in range(point.shape[1]):
                for iy in range(point.shape[1]):
                    x = point[0, ix]
                    y = point[1, iy]
                    z[ix, iy] = find_z_point(x, y)

        else:
            x = point[0]
            y = point[1]
            z = find_z_point(x, y)

        return z


if __name__ == '__main__':

    def test_2d():
        cs2d = CubicSpline2d(0, 0, 4, 4, 50, 50, 0, 0)
        Y, X = np.meshgrid(cs2d.ygrid, cs2d.xgrid)
        Z = np.sin(X) - np.cos(Y ** 2)
        cs2d.coefs(Z)
        xtest = np.r_[0.0:4.0:100j]
        ytest = np.r_[0.0:4.0:100j]
        ztest = cs2d.eval(np.row_stack([xtest, ytest]))

        # Build grid and get exact solution
        yy, xx = np.meshgrid(ytest, xtest)
        zz = np.sin(xx) - np.cos(yy ** 2)

        # Compute errors
        max_abs_err = np.abs(ztest - zz).max()

        print 'Max absolute error is ', max_abs_err
        return max_abs_err



    test_2d()
