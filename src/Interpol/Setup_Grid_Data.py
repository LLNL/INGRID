#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:31:15 2019

@author: watkins35
"""
from __future__ import division, print_function
import numpy as np
from bicubicf import bicubic
import matplotlib.pyplot as plt


class Efit_Data:

    def __init__(self, rmin=0.0, rmax=1.0, nr=10, zmin=0.0, zmax=2.0, nz=20,
                 name='unnamed'):
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

        self.name = name

    def get_v(self, tag='v'):
        """ returns the entire array of v, vr, vz, or vrz.
        If you want a single value use self.get_psi
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
        """ sets a value for v, vr, vz, or vrz. Coords must be a tuple or list.
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
        Attempts to use finite differences to increase accuracy at the
        boundaries.
        These formulas are derived from Taylor Series representations"""
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
        for i in range(1, self.nr-1):
            for j in range(1, self.nz-1):
                vr[i, j] = (self.v[i+1, j] - self.v[i-1, j])/2/dr

                vz[i, j] = (self.v[i, j+1] - self.v[i, j-1])/2/dz

                vrz[i, j] = (self.v[i+1, j+1] + self.v[i-1, j-1]
                             - self.v[i-1, j+1] - self.v[i+1, j-1])/4/dr/dz
        # missed a row in x
        for i in range(1, self.nz-1):
            for j in [0, self.nr-1]:
                vr[i, j] = (self.v[i+1, j] - self.v[i-1, j])/2/dr
        # and in y
        for i in [0, self.nz-1]:
            for j in range(1, self.nr-1):
                vz[i, j] = (self.v[i, j+1] - self.v[i, j-1])/2/dz

        # forward difference accuracy h^2
        for j in range(self.nz):
            vr[0, j] = (4*self.v[1, j] - self.v[2, j]
                        - 3*self.v[0, j])/2/dr
            vr[-1, j] = -(4*self.v[-2, j] - self.v[-3, j]
                          - 3*self.v[-1, j])/2/dr

        for i in range(self.nr):
            vz[i, 0] = (4*self.v[i, 1] - self.v[i, 2]
                        - 3*self.v[i, 0])/2/dz
            vz[i, -1] = -(4*self.v[i, -2] - self.v[i, -3]
                          - 3*self.v[i, -1])/2/dz

        # cross derivative on the edges
        for j in range(2, nz-2):
            i = 0  # left Edge
            vrz[i, j] = (f[i+2, j+2] - 2 * f[i+1, j+1] + 2 * f[i+1, j-1]
                         - f[i+2, j-2]) / 4/dr/dz
            i = nr-1  # right edge
            vrz[i, j] = (f[i-2, j+2] - 2 * f[i-2, j+1] + 2 * f[i-1, j-1]
                         - f[i-2, j-2]) / 4/dr/dz

        for i in range(2, nr-2):
            j = 0  # bottom edge
            vrz[i, j] = (f[i-1, j+2] - 2*f[i-1, j+1] + 2*f[i+1, j+1]
                         - f[i+2, j+2]) / 4/dr/dz
            j = nz-1  # top edge
            vrz[i, j] = (f[i-1, j-2] - 2*f[i-1, j-1] + 2*f[i+1, j-1]
                         - f[i+2, j-2]) / 4/dr/dz

        # cross derivatives at the corners
        for i, j in [[0, 0], [0, 1], [1, 0]]:
            # bottom left
            vrz[i, j] = (f[i, j] - f[i+1, j] + f[i+1, j+1] - f[i, j+1])/dr/dz

        for i, j in [[nr-2, 0], [nr-1, 0], [nr-1, 1]]:
            # bottom right
            vrz[i, j] = - (f[i, j] - f[i, j+1] + f[i-1, j+1] - f[i-1, j])/dr/dz

        for i, j in [[0, nz-2], [0, nz-1], [1, nz-1]]:
            # top left
            vrz[i, j] = - (f[i, j] - f[i, j-1] + f[i+1, j-1] - f[i+1, j])/dr/dz

        for i, j in [[nr-2, nz-1], [nr-1, nr-1], [nr-1, nz-2]]:
            # top right
            vrz[i, j] = (f[i, j] - f[i-1, j] + f[i-1, j-1] - f[i, j-1])/dr/dz

        # reload the new derivative values
        self.vr = vr
        self.vz = vz
        self.vrz = vrz
        end = time()
        print("Time calculating derivatives", end-start)

    def locate_cell(self, r0, z0):
        """ Inputs:
               r0,z0 - tested point

               Outputs:
                  node indices of grid cell encompassing tested point
        """
        # indices of lower left vertex of the cell
        # (and make sure they are within [0,nr-1])
        ir = int(max([min([np.floor((r0-self.rmin)/self.dr),
                           self.nr-2]), 0]))
        iz = int(max([min([np.floor((z0-self.zmin)/self.dz),
                           self.nz-2]), 0]))

        ircell = [ir, ir+1, ir, ir+1]
        izcell = [iz, iz, iz+1, iz+1]

        return {'ir': ircell, 'iz': izcell}

    def get_psi(self, r0, z0, tag='v'):
        """ find grid cell encompassing (r0,z0)
        note: grid is the crude grid
        tag is the type of derivative we want: v, vr, vz, vrz
        if nothing is provided, it assumes no derivative (v).
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
        r0norm = (r0-rcell[0])/self.dr
        z0norm = (z0-zcell[0])/self.dz
        res = bicubic(fcell, frcell, fzcell, frzcell,
                      x0=r0norm, y0=z0norm, derivs=tag)

        if tag == 'vr':
            res /= self.dr
        elif tag == 'vz':
            res /= self.dz
        elif tag == 'vrz':
            res /= self.dr*self.dz

        return res

    def plot_levels(self, level=1.0, color='red'):
        """ level must be a float """
        # draw contour line on top of existing figure
        level = float(level)
        self.ax.contour(self.r, self.z, self.v, level, colors=color)

    def plot_data(self, nlev=30):
        """ generates the plot that we will be able to manipulate
        using the root finder """
        fig = plt.figure(self.name, figsize=(6, 10))
        self.ax = fig.add_subplot(111)
        lev = (self.v.min() + (self.v.max()
               - self.v.min()) * np.arange(nlev) / (nlev-1))
        self.ax.contour(self.r, self.z, self.v, lev)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title('psi')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.xlim(self.rmin, self.rmax)
        plt.ylim(self.zmin, self.zmax)
        plt.ion()
        plt.show()


if __name__ == "__main__":
    grid = Efit_Data()
