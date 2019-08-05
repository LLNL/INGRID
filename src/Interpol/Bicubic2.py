#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 13:24:06 2019

@author: Joey
"""

#from numpy import array, reshape, matmul, transpose
import numpy as np


def bicubic(f, fx, fy, fxy, x0, y0, derivs=None):
    """ Bicubic interpolation on unit square [0,1]x[0,1].
    Returns the value of the interpolated polynomial at the given interior
    point (x0,y0).
    """
    # Make sure that f, fx, fy, fxy are arrays.
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
    if derivs == 'v' or derivs is None:
        res = np.matmul([1, x0, x0**2, x0**3], np.matmul(alp, [1, y0, y0**2, y0**3]))

    elif derivs == 'vr':
        res = np.matmul([0, 1, 2*x0, 3*x0**2], np.matmul(alp, [1, y0, y0**2, y0**3]))

    elif derivs == 'vz':
        res = np.matmul([1, x0, x0**2, x0**3], np.matmul(alp, [0, 1, 2*y0, 3*y0**2]))

    elif derivs == 'vrz':
        res = np.matmul([0, 1, 2*x0, 3*x0**2], np.matmul(alp, [0, 1, 2*y0, 3*y0**2]))

    return res