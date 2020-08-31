#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:29:28 2019

@author: watkins35
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from time import time
from Test_Functions import get_f
from EfitData import EfitData


def test_interpol(option=2, nfine=100, ncrude=10, tag='v'):
    """ Drives a simple test of a bicubic interpolation function.
    
    Parameters
    ----------
    option : int, optional
        Specifies which of several functions to demonstrate over.
        Accepts 1, 11, 12, 13, 2, 3, 4, 5, 51, 52
        See Test_Functions.get_f for more detail
    nfine : int, optional
        Density of the fine grid we will interpolate onto.
    ncrude : int, optional
        Density of the crude grid the sample data will be generated
        onto.
    tag : str, optional
        Specify if it is wanted to test the derivative interpolation
        methods. Accepts 'v', 'vr', 'vz', 'vrz'.

    """
    rmin = 0.0
    rmax = 2.0
    zmin = 0.0
    zmax = 1.0

    # accurate "exact" results for reference on fine grid
    grid0 = EfitData(nr=nfine, nz=nfine,
                     rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax)

    # crude mesh data for testing interpolation
    grid1 = EfitData(nr=ncrude, nz=ncrude,
                     rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax)

    # interpolated results on relatively fine grid
    grid2 = EfitData(nr=nfine, nz=nfine,
                     rmin=rmin, rmax=rmax, zmin=zmin, zmax=zmax)

    # get the function value, and analytic derivs
    grid0.v = get_f(grid0.r, grid0.z, option)

    # retrieve the data for our function on the crude grid
    grid1.v = get_f(grid1.r, grid1.z, option)

    # on the fine grid calculate derivatives by finite-difference
    # replace anayltic derivs by finite-difference
    grid0.Calculate_PDeriv(unit_spacing=False)

    # on the crude grid calculate derivatives by finite-difference
    grid1.Calculate_PDeriv()
    # don't adjust for dx/dy because of properties of the unit cell

    # determines where the contour lines are drawn
    nlev = 30
    lev = (grid0.get_v(tag).min() + (grid0.get_v(tag).max()
           - grid0.get_v(tag).min()) * np.arange(nlev) / (nlev - 1))

    # plot of what the function is supposed to look like
    plt.figure("interpol Demo {}-{}.".format(option, tag))
    plt.contour(grid0.r, grid0.z, grid0.get_v(tag), lev, colors='black')
    plt.gca().set_aspect('equal', adjustable='box')  # make it a square
    plt.title(tag)
    plt.ion()  # interactive mode. Let's us see the change in the plot later.
    plt.draw()
    plt.pause(.4)  # help us see the difference

    # interpolate to grid2 - the fine grid
    start = time()
    print("Beginning interpolation")
    # use same levels as for the exact contour plot
    for i in range(grid2.nr):
        for j in range(grid2.nz):
            x0 = grid2.r[i, j]
            y0 = grid2.z[i, j]
            # this is where bicubic interpolation comes into play
            grid2.set_v(grid1.get_psi(x0, y0, tag=tag), (i, j), tag)

    plt.contour(grid2.r, grid2.z, grid2.get_v(tag), lev, colors='red')

    sfactor = np.max(grid0.get_v(tag)) / np.max(grid2.get_v(tag))
    print("Needed scaling factor: {}".format(sfactor))
    end = time()
    print("Interpolation took {} seconds".format(end - start))
    plt.show()


if __name__ == "__main__":
    plt.close('all')
    test_interpol(option=1, tag='vr')
