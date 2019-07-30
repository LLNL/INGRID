#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:53:16 2019

@author: watkins35
"""
import numpy as np


def get_f(x, y, option=11, x0=0, y0=0):
    """ some functions for testing.
    Parameters
    ----------
    x, y: 1d or 2d arrays that define the grid.
    option: determines which function to return
    x0, y0: options to shift the centers of a certain functions. Default is 0
    """
    if option == 1:
        f = np.sin(np.pi*x)*np.cos(np.pi*y)

    elif option == 11:
        # diagonal lines
        f = x+y

    elif option == 12:
        # horizontal lines
        f = x

    elif option == 13:
        # vertical lines
        f = y

    elif option == 2:
        f = np.sin(np.pi*x)+np.cos(3*np.pi*y)

    elif option == 3:
        f = (np.cos(np.pi*(x+2*y))
             + np.cos(3*np.pi*y))

    elif option == 4:
        f = (np.exp(-(x-np.mean(x))**2
             - (y-np.mean(y))**2))

    elif option == 5:
        # gaussian funciton
        f = np.exp(-(x-x0)**2 - (y-y0)**2)

    elif option == 51:
        # x derivative of the gaussian
        f = np.exp(-2*(x-x0)**2 - (y-y0)**2)

    elif option == 52:
        # two gaussians?
        f = (np.exp(-2*(x-x0)**2 - (y-y0)**2) +
             np.exp(-2*(x+x0)**2 - (y+y0)**2))

    return f
