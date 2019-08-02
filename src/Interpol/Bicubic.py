#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:24:53 2019

@author: watkins35
"""
from numpy import array, reshape, matmul, transpose


def bicubic(f, fx, fy, fxy, x0, y0, derivs=None):
    """ Bicubic interpolation on unit square [0,1]x[0,1].
    Returns the value of the interpolated polynomial at the given interior
    point (x0,y0).
    """
    # Make sure that f, fx, fy, fxy are arrays.
    # All known quantities on 4 vertices.
    xall = array([f, fx, fy, fxy]).flatten()

    # inverse of the matrix for determining coef
    A = array([[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [ 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
               [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
               [ 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0],
               [ 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0],
               [-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0],
               [ 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0],
               [ 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1],
               [-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1],
               [ 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
               [ 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0],
               [-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1],
               [ 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1]])

    alp = matmul(A, xall)
    alp = transpose(reshape(alp, (4, 4)))

    res = 0.0
    # build the polynomial using the relative coordinates
    if derivs == 'v' or derivs is None:
        for i in range(4):
            for j in range(4):
                res += alp[i, j] * (x0**i) * (y0**j)

    elif derivs == 'vr':
        for i in range(1, 4):
            for j in range(4):
                res += alp[i, j]*(i*x0**(i-1))*(y0**j)

    elif derivs == 'vz':
        for i in range(4):
            for j in range(1, 4):
                res += alp[i, j]*(x0**i)*(j*y0**(j-1))

    elif derivs == 'vrz':
        for i in range(1, 4):
            for j in range(1, 4):
                res += alp[i, j]*i*x0**(i-1)*j*y0**(j-1)

    return res


if __name__ == "__main__":
    import numpy as np
    f = np.array([1, 1, 1, 1])
    fx = np.zeros(4)
    fy = np.zeros(4)
    fxy = np.zeros(4)

    x0, y0 = .5, .5

    res = bicubic(f, fx, fy, fxy, x0, y0)
    print("({},{}): {}".format(x0, y0, res))
