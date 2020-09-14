#!/usr/bin/env pythonu
# -*- coding: utf-8 -*-
"""
Module containing Bicubic class and EfitData class for handling all interpolation
related computations.
"""
from __future__ import division, print_function, absolute_import
import numpy as np
import matplotlib.pyplot as plt


def Bicubic(f, fx: 'array-like', fy: 'array-like', fxy: 'array-like', x0: float, y0: float, derivs: str = '') -> float:
    """ Bicubic interpolation on unit square [0,1]x[0,1].
    Returns the value of the interpolated polynomial at the given interior
    point (x0,y0).
    f, fx, fy, fxy can be arrays, lists, tulpes etc.
    They must each have four values and are ordered:

    .. math::

        f_{i,j}, f_{i,j+1}, f_{i+1,j}, f_{i+1,j+1}.

    This is the natural order produced by the locate cell method,
    in the EfitData.EFit_Data class, which how
    this function is usually called.

    Parameters
    ----------
    f : array-like
        Contains the value of the function of each cell coordinate.
    fx : array-like
        Contains the x or R derivative
    fy : array-like
        Contains the y or Z derivative
    fxy : array-like
        Contains the cross derivative, xy or RZ
    x0 : float
        Normalized x or R coordinate of the point of interest
    y0 : float
        Normalized y  or Z coordinate of the point
    derivs : str, optional
        Contains the tag for the type of derivative we want.
        Accepts 'v', 'vr', 'vz', 'vrz'

    Returns
    -------
        The value of the function of its derivative at a given location
    """
    # All known quantities on 4 vertices.
    xall = np.array([[f[0], f[2], fy[0], fy[2]],
                     [f[1], f[3], fy[1], fy[3]],
                     [fx[0], fx[2], fxy[0], fxy[2]],
                     [fx[1], fx[3], fxy[1], fxy[3]]])

    # use two coefficient matrices
    A1 = np.array([[1, 0, 0, 0],
                   [0, 0, 1, 0],
                   [-3, 3, -2, -1],
                   [2, -2, 1, 1]])

    A2 = np.array([[1, 0, -3, 2],
                   [0, 0, 3, -2],
                   [0, 1, -2, 1],
                   [0, 0, -1, 1]])

    alp = np.matmul(A1, np.matmul(xall, A2))

    # build the polynomial using the relative coordinates
    if derivs == 'v' or derivs == '':
        res = np.matmul([1, x0, x0**2, x0**3], np.matmul(alp, [1, y0, y0**2, y0**3]), dtype=np.ndarray)

    elif derivs == 'vr':
        res = np.matmul([0, 1, 2 * x0, 3 * x0**2], np.matmul(alp, [1, y0, y0**2, y0**3]))

    elif derivs == 'vz':
        res = np.matmul([1, x0, x0**2, x0**3], np.matmul(alp, [0, 1, 2 * y0, 3 * y0**2]))

    elif derivs == 'vrz':
        res = np.matmul([0, 1, 2 * x0, 3 * x0**2], np.matmul(alp, [0, 1, 2 * y0, 3 * y0**2]))

    elif derivs == 'vrr':
        res = np.matmul([0, 0, 2, 6 * x0], np.matmul(alp, [1, y0, y0**2, y0**3]))

    elif derivs == 'vzz':
        res = np.matmul([1, x0, x0**2, x0**3], np.matmul(alp, [0, 0, 2, 6 * y0]))

    return res


class EfitData:
    """
    Structure to store the rectangular grid of psi data. It uses
    cylindrical coordinates, where R and Z are similar to the cartesian
    x and y. The phi components goes away due to the symmetry of a
    tokamak.


    Parameters
    ----------
    rmin : float, optional
        left boundary of the grid

    rmax : float, optional
        right boundary

    nr : int, optional
        number of grid points in the R direction

    zmin : float, optional
        bottom boundary for the grid

    zmax : float, optional
        top boundary

    nz : int, optional
        number of grid points in the Z direction

    name : str, optional
        Specify the title of the figure the data will be plotted on.
    """

    def __init__(self, rmin=0.0, rmax=1.0, nr=10, zmin=0.0, zmax=2.0, nz=20,
                 rcenter=1.6955000, bcenter=-2.1094041, rlimiter=None, zlimiter=None,
                 rmagx=0.0, zmagx=0.0, name='unnamed', parent=None):
        r, dr = np.linspace(rmin, rmax, nr, retstep=True)
        z, dz = np.linspace(zmin, zmax, nz, retstep=True)
        rgrid, zgrid = np.meshgrid(r, z, indexing='ij')
        value = np.zeros([nr, nz])

        self.nr = nr
        self.nz = nz
        self.rmin = rmin
        self.rmax = rmax
        self.zmin = zmin
        self.zmax = zmax
        self.r = rgrid
        self.z = zgrid
        self.v = value
        self.vr = value
        self.vz = value
        self.vrz = value
        self.dr = dr
        self.dz = dz
        self.rcenter = rcenter
        self.bcenter = bcenter
        self.rmagx = rmagx
        self.zmagx = zmagx
        self.rlimiter = rlimiter
        self.zlimiter = zlimiter
        self.name = name
        self.parent = parent
        self.psi_levels = {}

    def Gradient(self, xy: tuple) -> 'np.ndarray':
        """ Combines the first partial derivatives to solve the system for
        maximum, minimum, and saddle locations.

        Parameters
        ----------
        xy : array-like
            Contains x and y. Ex: xy = (x0, y0).

        Returns
        -------
        F : array
            Vector function to be used in find root.
        """
        # combine the deriv functions to solve the system
        x, y = xy
        F = np.zeros(2)
        F[0] = self.get_psi(x, y, tag='vr')
        F[1] = self.get_psi(x, y, tag='vz')
        return F

    def Hessian(self, xy: tuple) -> 'np.ndarray':
        """ Compute the Hessian at a point.

        Parameters
        ----------
        xy : array-like
            Contains x and y. Ex: xy = (x0, y0).

        Returns
        -------
        H : array
            Numpy array of shape (2, 2) representing the Hessian at xy.
        """
        x, y = xy
        H = np.zeros((2, 2))
        H[0, 0] = self.get_psi(x, y, 'vrr')
        H[1, 1] = self.get_psi(x, y, 'vzz')
        H[0, 1] = self.get_psi(x, y, 'vrz')
        H[1, 0] = self.get_psi(x, y, 'vrz')
        return H

    def PsiFunction(self, xy):
        x, y = xy
        return self.get_psi(x, y)

    def get_v(self, tag: str = 'v') -> 'np.ndarray':
        """ Returns the entire array of v, vr, vz, or vrz.
        If you want a single value use self.get_psi

        Parameters
        ----------
        tag : str, optional
            Specify the type of derivative. 'v', 'vr', 'vz', 'vrz'

        Returns
        -------
        ndarray
            value of the function or its derivative over the entire
            grid.
        """
        if tag == 'v':
            return self.v
        elif tag == 'vr':
            return self.vr
        elif tag == 'vz':
            return self.vz
        elif tag == 'vrz':
            return self.vrz

    def set_v(self, value, coords=None, tag='v'):
        """ sets a value for v, vr, vz, or vrz.

        Parameters
        ----------
        value : ndarray, float
            new set of values for the function. Must be the same shape
            as the grid if you are setting every value. Also accepts a
            single float for setting spevific values.

        coords : array-like, optional
            The coordinates of a single value, if you are setting one value.
            if set to none, it will set the entire value
        """
        if coords is not None:
            if tag == 'v':
                self.v[coords[0], coords[1]] = value
            elif tag == 'vr':
                self.vr[coords[0], coords[1]] = value
            elif tag == 'vz':
                self.vz[coords[0], coords[1]] = value
            elif tag == 'vrz':
                self.vrz[coords[0], coords[1]] = value
        else:
            if tag == 'v':
                self.v = value
            elif tag == 'vr':
                self.vr = value
            elif tag == 'vz':
                self.vz = value
            elif tag == 'vrz':
                self.vrz = value

    def Calculate_PDeriv(self, unit_spacing=True):
        """ Calculate partial derivatives at grid nodes.

        Use finite differences to increase accuracy at the
        boundaries.

        These formulas are derived from Taylor Series representations.
        Values for vr, vz, and vrz are produced and saved within the
        grid structure.

        Parameters
        ----------
        unit_spacing : bool, optional
            Allow for the derivatives to be calculated on a unit cell,
            or on a grid with any other spacing. Default is true.
        """
        # planning to make this more efficient using array slicing
        # need a place to store the new values of the grid
        from time import time
        print("Beginning Derivative calculation")
        start = time()

        f = self.v

        vr = np.zeros_like(self.vr)
        vz = np.zeros_like(self.vz)
        vrz = np.zeros_like(self.vrz)
        nr = self.nr
        nz = self.nz

        # step size becomes one because of the unit grid
        if unit_spacing:
            dr = 1
            dz = 1
        else:
            dr = self.dr
            dz = self.dz

        # inner square - can use centered difference
        for i in range(1, self.nr - 1):
            for j in range(1, self.nz - 1):
                vr[i, j] = (self.v[i + 1, j] - self.v[i - 1, j]) / 2 / dr

                vz[i, j] = (self.v[i, j + 1] - self.v[i, j - 1]) / 2 / dz

                vrz[i, j] = (self.v[i + 1, j + 1] + self.v[i - 1, j - 1]
                             - self.v[i - 1, j + 1] - self.v[i + 1, j - 1]) / 4 / dr / dz

        # missed a row in x
        for i in range(1, self.nr - 1):
            for j in [0, self.nz - 1]:
                vr[i, j] = (self.v[i + 1, j] - self.v[i - 1, j]) / 2 / dr

        # and in y
        for i in [0, self.nr - 1]:
            for j in range(1, self.nz - 1):
                vz[i, j] = (self.v[i, j + 1] - self.v[i, j - 1]) / 2 / dz

        # forward difference accuracy h^2
        for j in range(self.nz):
            vr[0, j] = (4 * self.v[1, j] - self.v[2, j]
                        - 3 * self.v[0, j]) / 2 / dr
            vr[-1, j] = -(4 * self.v[-2, j] - self.v[-3, j]
                          - 3 * self.v[-1, j]) / 2 / dr

        for i in range(self.nr):
            vz[i, 0] = (4 * self.v[i, 1] - self.v[i, 2]
                        - 3 * self.v[i, 0]) / 2 / dz
            vz[i, -1] = -(4 * self.v[i, -2] - self.v[i, -3]
                          - 3 * self.v[i, -1]) / 2 / dz

        # cross derivative on the edges
        for j in range(2, nz - 2):
            i = 0  # left Edge
            vrz[i, j] = (f[i + 2, j + 2] - 2 * f[i + 1, j + 1] + 2 * f[i + 1, j - 1]
                         - f[i + 2, j - 2]) / 4 / dr / dz
            i = nr - 1  # right edge
            vrz[i, j] = (f[i - 2, j + 2] - 2 * f[i - 2, j + 1] + 2 * f[i - 1, j - 1]
                         - f[i - 2, j - 2]) / 4 / dr / dz

        for i in range(2, nr - 2):
            j = 0  # bottom edge
            vrz[i, j] = (f[i - 1, j + 2] - 2 * f[i - 1, j + 1] + 2 * f[i + 1, j + 1]
                         - f[i + 2, j + 2]) / 4 / dr / dz
            j = nz - 1  # top edge
            vrz[i, j] = (f[i - 1, j - 2] - 2 * f[i - 1, j - 1] + 2 * f[i + 1, j - 1]
                         - f[i + 2, j - 2]) / 4 / dr / dz

        # cross derivatives at the corners
        for i, j in [[0, 0], [0, 1], [1, 0]]:
            # bottom left
            vrz[i, j] = (f[i, j] - f[i + 1, j] + f[i + 1, j + 1] - f[i, j + 1]) / dr / dz

        for i, j in [[nr - 2, 0], [nr - 1, 0], [nr - 1, 1]]:
            # bottom right
            vrz[i, j] = - (f[i, j] - f[i, j + 1] + f[i - 1, j + 1] - f[i - 1, j]) / dr / dz

        for i, j in [[0, nz - 2], [0, nz - 1], [1, nz - 1]]:
            # top left
            vrz[i, j] = - (f[i, j] - f[i, j - 1] + f[i + 1, j - 1] - f[i + 1, j]) / dr / dz

        for i, j in [[nr - 2, nz - 1], [nr - 1, nr - 1], [nr - 1, nz - 2]]:
            # top right
            vrz[i, j] = (f[i, j] - f[i - 1, j] + f[i - 1, j - 1] - f[i, j - 1]) / dr / dz

        # reload the new derivative values
        self.vr = vr
        self.vz = vz
        self.vrz = vrz
        end = time()
        print("Time calculating derivatives", end - start)

    def locate_cell(self, r0: float, z0: float) -> dict:
        """
        Locate the cell on the rectangular grid that surrounds
        a point of interest.

        Parameters
        ----------
        r0 : float
            R or x coordinate of the point
        z0 : float
            Z or y coordinate of the tested point

        Returns
        -------
            A dict with node indices of grid cell encompassing tested point
        """
        # indices of lower left vertex of the cell
        # (and make sure they are within [0,nr-1])
        ir = int(max([min([np.floor((r0 - self.rmin) / self.dr),
                           self.nr - 2]), 0]))
        iz = int(max([min([np.floor((z0 - self.zmin) / self.dz),
                           self.nz - 2]), 0]))

        ircell = [ir, ir + 1, ir, ir + 1]
        izcell = [iz, iz, iz + 1, iz + 1]

        return {'ir': ircell, 'iz': izcell}

    def get_psi(self, r0, z0, tag='v'):
        """ find grid cell encompassing (r0,z0)
        note: grid is the crude grid. Uses Bicubic Interpolation
        to calculate the exact value at the point. Useful for
        finding information inbetween grid points.

        Parameters
        ----------
        r0 : float
            R coordinate of the point of interest
        z0 : float
            Z coordinate of same point.
        tag : str, optional
            tag is the type of derivative we want: v, vr, vz, vrz
            if nothing is provided, it assumes no derivative (v).

        Returns
        -------
        float
            Value of psi or its derviative at the coordinate specified.
        """
        cell = self.locate_cell(r0, z0)
        rcell = self.r[cell['ir'], cell['iz']]
        zcell = self.z[cell['ir'], cell['iz']]
        fcell = self.v[cell['ir'], cell['iz']]
        frcell = self.vr[cell['ir'], cell['iz']]
        fzcell = self.vz[cell['ir'], cell['iz']]
        frzcell = self.vrz[cell['ir'], cell['iz']]

        # ==BICUBIC INTERPOLATION==
        # NOTE: assumes unit cell
        r0norm = (r0 - rcell[0]) / self.dr
        z0norm = (z0 - zcell[0]) / self.dz
        res = Bicubic(fcell, frcell, fzcell, frzcell,
                      x0=r0norm, y0=z0norm, derivs=tag)

        if tag == 'vr':
            res /= self.dr
        elif tag == 'vz':
            res /= self.dz
        elif tag == 'vrz':
            res /= self.dr * self.dz
        elif tag == 'vrr':
            res /= self.dr * self.dr
        elif tag == 'vzz':
            res /= self.dz * self.dz

        return res

    def plot_levels(self, level=1.0, color='red'):
        """
        This function is useful if you need to quickly see
        where a particular line of constant psi is. It in't able to store
        points of intersection, and cannot be generalized. If you
        need just a segment of psi, use the draw_lines method in the
        line tracing class.

        Parameters
        ----------
        level : float, optional
            Value of psi you wish to see
        color : str, optional
            color of the line.

        """
        # draw contour line on top of existing figure
        level = float(level)
        self.ax.contour(self.r, self.z, self.v, level, colors=color)

    def PlotLevel(self: object, level: float = 1.0, color: str = 'red', label: str = '', linestyles: str = 'solid') -> None:
        """
        Plot a psi level and provide it a label.

        This function is useful for management of psi boundaries
        such as 'psi_pf', 'psi_core', etc and ensuring the contour will
        be properly replotted (no duplicate of same label).

        Parameters
        ----------
        level : float, optional
            Psi level to plot. Default to 1.0 (separatrix of normalized psi)
        color : str, optional
            Color to pass to matplotlib contour function
        label : str, optional
            Label to associate with the psi level
        linestyles : str, optional
            Line style to pass to matplotlib contour function
        """
        try:
            self.psi_levels[label].collections[0].remove()
            self.psi_levels[label] = plt.contour(self.r, self.z, self.v, [float(level)], colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)
        except:
            self.psi_levels[label] = plt.contour(self.r, self.z, self.v, [float(level)], colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)

    def plot_data(self, nlevs=30, interactive=True, fig=None, ax=None, view_mode='filled'):
        """ generates the plot that we will be able to manipulate
        using the root finder

        Parameters
        ----------
        nlev : int, optional
            number of levels we want to be plotted
        """
        lev = (self.v.min() + (self.v.max()
                - self.v.min()) * np.arange(nlevs) / (nlevs - 1))
        self.fig = fig if fig is not None else plt.figure('INGRID: ' + self.name, figsize=(6, 10))
        self.fig.subplots_adjust(bottom=0.2)
        self.ax = ax if ax is not None else self.fig.add_subplot(111)
        if view_mode == 'lines':
            self.ax.contour(self.r, self.z, self.v, lev, cmap='gist_gray')
        elif view_mode == 'filled':
            self.ax.contourf(self.r, self.z, self.v, lev, cmap='gist_gray')
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_title(f'{self.name}')
        self.ax.set_xlabel('R')
        self.ax.set_ylabel('Z')
        self.ax.set_xlim(self.rmin, self.rmax)
        self.ax.set_ylim(self.zmin, self.zmax)
        if interactive:
            plt.ion()
        fig.show()

    def clear_plot(self):
        if plt.get_fignums():
            plt.clf()
        else:
            pass
