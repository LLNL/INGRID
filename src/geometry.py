#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:00:04 2019

@author: watkins35
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


class Vector:
    def __init__(self, xy, origin):
        """ xy is of the form xy = (x, y) """
        self.x, self.y = xy
        self.xorigin = origin[0]
        self.yorigin = origin[1]
        self.xnorm = self.x - self.xorigin
        self.ynorm = self.y - self.yorigin
        self.quadrant = (int(np.sign(self.xnorm)), int(np.sign(self.ynorm)))

    def arr(self):
        """ returns the vector as an array """
        return np.array([self.xnorm, self.ynorm])

    def mag(self):
        """ computes the magnitude of the vector """
        return np.linalg.norm(self.arr())


class Point:
    def __init__(self, *pts):
        if np.shape(pts) == (2,):
            self.x, self.y = pts
        elif np.shape(pts) == (1, 2):
            self.x, self.y = pts[0]
        else:
            print('incompatible form')
            print(np.shape(pts), pts)

    def psi(self, grid, tag='v'):
        """ must pass in the grid upon which the value of psi is to be
        calculated on. Must be the Efit grid object.
        also accepts 'tag' as an optional parameter. This is to specify the
        type of psi derivative, if desired. accepts 'v', 'vr', 'vz', 'vrz'.
        default is 'v', or no derivative.
        """
        return grid.get_psi(self.x, self.y, tag)

    def plot(self):
        plt.plot(self.x, self.y, 'x')


#class Line:
#    """ line object """
#
#    def __init__(self, *points):
#        """ points are of the form p = (x, y)
#        # an ordered set of points defines a line
#        line = (p1, p2, ...)
#        """
#        if isinstance(points[0], Point):
#            self.p = points  # list of objects
#            self.xval = [p.x for p in points]
#            self.yval = [p.y for p in points]
#        else:
#            self.p = np.array(points)
#            self.xval = self.p[:, 0]
#            self.yval = self.p[:, 1]
#
#    def plot(self, color='#1f77b4'):
#        """ defaults to a light blue """
#        plt.plot(self.xval, self.yval, color=color)
#
#    def print_points(self):
#        print([(p.x, p.y) for p in self.p])
#
#    def divisions(self, num):
#        self.xs = np.linspace(self.p[0].x, self.p[-1].x, num)
#        self.ys = np.linspace(self.p[0].y, self.p[-1].y, num)

class Line:
    """ line object """

    def __init__(self, points, psi=None):
        """ points are of the form p = (x, y)
        must pass in a list, tuple, array...
        # an ordered set of points defines a line
        line = (p1, p2, ...)
        """
        self.p = points
        self.xval = [p.x for p in points]
        self.yval = [p.y for p in points]
        
        if psi is not None:
            self.psi = psi
            
    def reverse(self):    
        # points in the other direction
        self.p = self.p[::-1]
        return self
        
#    def straighten(self):
#        # removes interior points
#        self.p = [self.p[0], self.p[-1]]
#        return self

    def plot(self, color='#1f77b4'):
        """ defaults to a light blue """
        plt.plot(self.xval, self.yval, color=color, linewidth='2')

    def print_points(self):
        print([(p.x, p.y) for p in self.p])

    def divisions(self, num):
        self.xs = np.linspace(self.p[0].x, self.p[-1].x, num)
        self.ys = np.linspace(self.p[0].y, self.p[-1].y, num)
        
    def points(self):
        return ((self.p[0].x, self.p[0].y), (self.p[-1].x, self.p[-1].y))


class Patch:
    """ each patch contains a grid, and has it's own shape
    accepts Point objects and lists of Points"""
#    def __init__(self, *points):
#        self.p = []
#        for pt in points:
#            if type(pt) == list:
#                [self.p.append(p) for p in pt]
#            else:
#                self.p.append(pt)
    
    def __init__(self, lines):
        """ Each patch defined by four lines in order - 
        Nline, Eline, Sline, Wline
        order points to go clockwise
        """
        self.lines = lines
        self.N = lines[0]
        self.E = lines[1]
        self.S = lines[2]
        self.W = lines[3]
        
        self.p = list(self.N.p) + list(self.S.p)
        

    def plot_border(self, color='red'):
        """ draw solid borders around the patch """
        for line in self.lines:
            line.plot(color)

#    def plot_borders(self, color='#1f77b4'):
#        """ draw solid borders around the patch """
#        pts = self.p + [self.p[0]]
#        for i in range(len(self.p)):
#            line = Line(pts[i], pts[i+1])
#            line.plot(color)

    def fill(self, color='lightsalmon'):
        """ shades in the patch with a given color """
        x = [p.x for p in self.p]
        y = [p.y for p in self.p]
        plt.fill(x, y, facecolor=color)


def calc_mid_point(v1, v2):
    """
    v1 must be furthest right in a counter clockwise direction
    v1 and v2 are vectors of equal length
    """
    # bisection
    theta = np.arccos(np.dot(v1.arr(), v2.arr())/(v1.mag() * v2.mag())) / 2.

    # check with quadrant the vectors are in and
    # compute the angles appropriately
    if v1.quadrant == (1, 1):
        # NE quadrant
        angle = np.arccos(v1.xnorm/v1.mag())
    elif v1.quadrant == (-1, 1):
        # NW quadrant
        angle = np.pi - np.arcsin(v1.ynorm/v1.mag())
    elif v1.quadrant == (-1, -1):
        # SW quadrant
        angle = np.pi + np.arctan(v1.ynorm/v1.xnorm)
    elif v1.quadrant == (1, -1):
        # SE quadrant
        angle = - np.arccos(v1.xnorm/v1.mag())
    else:
        print("Something went wrong")

    x = v1.xorigin + v1.mag() * np.cos(theta+angle)
    y = v1.yorigin + v1.mag() * np.sin(theta+angle)
    return x, y


def test2points(p1, p2, line):
    """ check if two points are on opposite sides of a line.
    the signs are different if on opposite sides of the line
    """
#    x1, y1 = line.p[0].x, line.p[0].y
#    x2, y2 = line.p[-1].x, line.p[-1].y
    
    (x1, y1), (x2, y2) = line
    x, y = p1
    a, b = p2
    d1 = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)
    d2 = (a-x1)*(y2-y1)-(b-y1)*(x2-x1)
    return np.sign(d1),  np.sign(d2)


def intersect(line1, line2):
    """ lines of the form A = ((x, y), (x, y)) """
    def line(x, line):
        """ point slope form """
        (x1, y1), (x2, y2) = line
        return (y2-y1)/(x2-x1) * (x - x1) + y1

    def f(xy):
        # parse the line, access any point
        # normalized to help solve for zero
        x, y = xy
        return np.array([y - line(x, line1),
                         y - line(x, line2)])

    # use the mean for initial guess
    (a, b), (c, d) = line1
    (i, j), (p, q) = line2
    guess = (np.mean([a, c, i, p]), np.mean([b, d, j, q]))
    r, z = fsolve(f, guess)
    return r, z
    

#def intersect(line1, line2):
#    """ accepts two Line objects """
#    def line(x, line):
#        """ point slope form """
##        (x1, y1), (x2, y2) = line.p[0].x, , line.p[-1]
#        x1 = line.p[0].x
#        y1 = line.p[0].y
#        x2 = line.p[-1].x
#        y2 = line.p[-1].y
#        return (y2-y1)/(x2-x1) * (x - x1) + y1
#
#    def f(xy):
#        # parse the line, access any point
#        # normalized to help solve for zero
#        x, y = xy
#        return np.array([y - line(x, line1),
#                         y - line(x, line2)])
#
#    # use the mean for initial guess
##    (a, b), (c, d) = line1
##    (i, j), (p, q) = line2
#    a, b = line1.p[0].x, line1.p[0].y
#    c, d = line1.p[-1].x, line1.p[-1].y
#    i, j = line2.p[0].x, line2.p[0].y
#    p, q = line2.p[-1].x, line2.p[-1].y
#    
#    guess = (np.mean([a, c, i, p]), np.mean([b, d, j, q]))
#    r, z = fsolve(f, guess)
#    return r, z



def truncate_list(old_x, old_y, x, y):
    """ accepts lists of x and y coordinates, and shortens them.
    ends at the given x and y values. assumes they are on the same line
    as the rest of the list.
    """

    # find the nearest index
    # this index may be an upper or lower bound
    ix = (abs(old_x - x)).argmin()
    iy = (abs(old_y - y)).argmin()

    # case ascending x
    if old_x[0] - old_x[1] < 0:
        if old_x[ix] > x:
            # upper bound
            new_x = np.copy(old_x[:ix+1])
        else:
            # lower bound
            new_x = np.copy(old_x[:ix+2])

    # case descending x
    else:  # if old_x[0] - old_x[1] > 0:
        if old_x[ix] > x:
            new_x = np.copy(old_x[:ix+2])
        else:
            new_x = np.copy(old_x[:ix+1])
    # now we can set the last value
    new_x[-1] = x

    # case ascending y
    if old_y[0] - old_y[1] < 0:
        if old_y[iy] > y:
            new_y = np.copy(old_y[:iy+1])
        else:
            new_y = np.copy(old_y[:iy+2])

    # case descending
    else:  # if old_y[0] - old_y[1] > 0:
        if old_y[iy] > y:
            new_y = np.copy(old_y[:iy+2])
        else:
            new_y = np.copy(old_y[:iy+1])
    # set the last value
    new_y[-1] = y
    return new_x, new_y


if __name__ == "__main__":
    p1 = Point(1, 2)
    p2 = Point(3, 2)
    p3 = Point(4, 3)
    p4 = Point(2, 5)

    # functionality of a patch
    fig = plt.figure('patches', figsize=(6, 10))
    plt.clf()
    # pass in list of points - typically four
    patch = Patch([p1, p2, p3, p4])
    # plot the boundary

#    patch.plot_bounds()
    patch.fill('lightsalmon')

    # borders
    patch.plot_borders('red')

    plt.xlim(0, 6)
    plt.ylim(0, 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Example Patch')
    plt.show()
