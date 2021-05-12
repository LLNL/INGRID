#!/usr/bin/env pythonu
# -*- coding: utf-8 -*-
"""
Module containing EfitData class for handling all interpolation
related computations.
"""
from __future__ import division, print_function, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from scipy.interpolate import RectBivariateSpline as rbs

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

    def init_bivariate_spline(self, r: 'np.ndarray', 
                                    z: 'np.ndarray', 
                                    v: 'np.ndarray') -> None:
        """ Initialize scipy.interpolate.RectBivariateSpline
        object for Bicubic interpolation.

        Sets class member v to crude EFIT grid.

        Parameters
        ----------
        r : array-like
            1-D array of r coordinates in strictly ascending order.
        z : array-like
            1-D array of z coordinates in strictly ascending order.
        v : array-like
            2-D array of EFIT data with shape (r.shape, z.shape)

        """
        self.v = v  # Crude EFIT grid.
        self.rbs = rbs(r, z, v)  # RectBivariateSpline object.

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

        lookup = {'v': (0, 0), 'vr': (1, 0), 'vrr': (2, 0),
                  'vz': (0, 1), 'vzz': (0, 2), 'vrz': (1, 1)
        }

        dx, dy = lookup[tag]
        return self.rbs(r0, z0, dx, dy)[0]

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

    def PlotLevel(self: object, level: float = 1.0, color: str = 'red', label: str = '', linestyles: str = 'solid',
                  refined: bool = True, refine_factor: int = 10) -> None:
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
        refined : bool, optional
            Plot level with hi-resolution cubic spline representation
        refine_factor: int, optional
            Refinement factor for to be passed to SciPy zoom method
        """

        data = self.v
        rgrid = self.r
        zgrid = self.z

        if refined is True:
            data = zoom(input=self.v, zoom=refine_factor)
            rgrid, zgrid = np.meshgrid(np.linspace(self.rmin, self.rmax, data.shape[0]),
                                       np.linspace(self.zmin, self.zmax, data.shape[1]),
                                       indexing='ij')
        try:
            self.psi_levels[label].collections[0].remove()
            self.psi_levels[label] = plt.contour(rgrid, zgrid, data, [float(level)], colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)
        except:
            self.psi_levels[label] = plt.contour(rgrid, zgrid, data, [float(level)], colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)

    def plot_data(self: object, nlevs: int = 30, interactive: bool = True, fig: object = None,
                  ax: object = None, view_mode: str = 'filled', refined: bool = True, refine_factor: int = 10):
        """
        Plot the EFIT data.

        Visualizes eqdsk file with either contour lines or filled contours.

        Parameters
        ----------
        nlev : int, optional
            number of levels we want to be plotted
        interactive : bool, optional
            Set matplotlib interactive mode on or off
        fig : object, optional
            Matplotlib figure handle
        ax : object, optional
            Matplotlib axes handle
        view_mode : str, optional
            Represent EFIT data with standard contour lines or filled contour lines.
            String value of 'filled' enables filled contours, whereas 'lines'
            omits filling of contours.
        refined : bool, optional
            Plot level with hi-resolution cubic spline representation
        refine_factor: int, optional
            Refinement factor for to be passed to SciPy zoom method
        """

        lev = self.v.min() + (self.v.max() - self.v.min()) * np.arange(nlevs) / (nlevs - 1)
        self.fig = fig if fig is not None else plt.figure('INGRID: ' + self.name, figsize=(8, 10))
        self.fig.subplots_adjust(bottom=0.075)
        self.ax = ax if ax is not None else self.fig.add_subplot(111)

        data = self.v
        rgrid = self.r
        zgrid = self.z

        if refined is True:
            data = zoom(input=self.v, zoom=refine_factor)
            rgrid, zgrid = np.meshgrid(np.linspace(self.rmin, self.rmax, data.shape[0]),
                                       np.linspace(self.zmin, self.zmax, data.shape[1]),
                                       indexing='ij')
        if view_mode == 'lines':
            self.ax.contour(rgrid, zgrid, data, lev, cmap='gist_gray')
        elif view_mode == 'filled':
            self.ax.contourf(rgrid, zgrid, data, lev, cmap='gist_gray')
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlabel('R')
        self.ax.set_ylabel('Z')
        self.ax.set_xlim(self.rmin, self.rmax)
        self.ax.set_ylim(self.zmin, self.zmax)
        if interactive:
            plt.ion()
        self.fig.show()

    def clear_plot(self):
        if plt.get_fignums():
            plt.clf()
        else:
            pass
