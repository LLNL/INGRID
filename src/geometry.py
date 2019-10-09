#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:00:04 2019

@author: watkins35
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit, root_scalar


class Vector:
    """
    Defines a vector from a nontrivial origin.
    
    Parameters
    ----------
    xy : array-like
        Location of the vector. It if of the form (x, y).
    origin : array-like
        Location of the origin. This is to adjust for not being at the
        origin of the axes. Of the form (x, y).
    """
    def __init__(self, xy, origin):
        self.x, self.y = xy
        self.xorigin = origin[0]
        self.yorigin = origin[1]
        self.xnorm = self.x - self.xorigin
        self.ynorm = self.y - self.yorigin
        self.quadrant = (int(np.sign(self.xnorm)), int(np.sign(self.ynorm)))

    def arr(self):
        """ 
        Returns
        -------
        ndarray
            Returns the vector as an array.
        """
        return np.array([self.xnorm, self.ynorm])

    def mag(self):
        """  
        Returns
        -------
        float
            Computes the magnitude, or length of the vector.
        """
        return np.linalg.norm(self.arr())


class Point:
    """
    Point object 
    
    Parameters
    ----------
    pts : array-like
        Accepts either two values x, y as floats, or 
        a single tuple/list value (x, y).
    """
    
    def __init__(self, *pts):
        if np.shape(pts) == (2,):
            self.x, self.y = pts
        elif np.shape(pts) == (1, 2):
            self.x, self.y = pts[0]
        else:
            print('incompatible form')
            print(np.shape(pts), pts)

    def psi(self, grid, tag='v'):
        """ 
        Parameters
        ----------
        grid : Setup_Grid_Data.Efit_Data
            Must pass in the grid upon which the value of psi is to be
            calculated on. Must be the Efit grid object.
        tag : str, optional
            This is to specify the type of psi derivative, if desired. 
            Accepts 'v', 'vr', 'vz', 'vrz'.
            The default is 'v', or no derivative.
         
        Returns
        -------
        float
            Calculate the value of psi at the point.
        """
        return grid.psi_norm.get_psi(self.x, self.y, tag)

    def plot(self):
        """ Places an x on the location of the point. """
        plt.plot(self.x, self.y, 'x')


class Line:
    """ 
    Line object, in which an ordered set of points defines a line.
    
    Parameters
    ----------
    points : list
        Points are of the form p = (x, y), and the list should be
        made up of multiple points. [p, p, p]...  
    """

    def __init__(self, points):
        self.p = points
        self.xval = [p.x for p in points]
        self.yval = [p.y for p in points]

            
    def reverse(self):    
        """ Points the line in the other direction.
        It is intended to be used right after generating a line  
        using the draw_line funciton from the line tracer. 
        For example; LineTracing.draw_line(args...).reverse().
        
        Returns
        -------
        self
            geometry.Line
        """
        self.p = self.p[::-1]
        return self

    def reverse_copy(self):
        """
        Returns a copy of the line in the reversed direction.
        Does not overwrite the current line.
        """

        return Line(self.p[::-1])

    def plot(self, color='#1f77b4'):
        """ Plots the line of the current figure.
        
        Parameters
        ----------
        color : str, optional
            Defaults to a light blue.
        """
        plt.plot(self.xval, self.yval, color=color, linewidth='2')

    def print_points(self):
        """ Prints each point in the line to the terminal. """
        print([(p.x, p.y) for p in self.p])

    def divisions(self, num):
        """
        Splits the line into discrete segments.
        
        Parameters
        ----------
        num : int
            Number of points in the segmented line.
        """
        self.xs = np.linspace(self.p[0].x, self.p[-1].x, num)
        self.ys = np.linspace(self.p[0].y, self.p[-1].y, num)
        
    def points(self):
        """ Returns the points in the line as a tuple. """
        return ((self.p[0].x, self.p[0].y), (self.p[-1].x, self.p[-1].y))


class Patch:
    """
    Each patch contains a grid, and has it's own shape.
    
    Parameters
    ----------
    lines : geometry.Line
        Each patch defined by four lines in order - 
        Nline, Eline, Sline, Wline - order points to go clockwise.

    """
    
    def __init__(self, lines):
        self.lines = lines
        self.N = lines[0]
        self.E = lines[1]
        self.S = lines[2]
        self.W = lines[3]
        
        # This is the border for the fill function
        # It need to only include N and S lines
        self.p = list(self.N.p) + list(self.S.p)
        

    def plot_border(self, color='red'):
        """
        Draw solid borders around the patch.
        
        Parameters
        ----------
        color : str, optional
            Defaults to red.
        """
        for line in self.lines:
            line.plot(color)

    def fill(self, color='lightsalmon'):
        """
        Shades in the patch with a given color
        
        Parameters
        ----------
        color : str, optional
            Defaults to a light salmon.
        """
        x = [p.x for p in self.p]
        y = [p.y for p in self.p]
        plt.fill(x, y, facecolor=color)
        
    def refine(self, grid):
        """ Divides a patch into smaller cells based on N and S lines,
        and the psi levels of E and W lines.
        
        Parameters
        ----------
        grid : Setup_Grid_Data.Efit_Data
            Requires the grid the patches were calculated on.
        """
        
        # TODO: need a more universal fit for the function
        def f(x, a, b, c, d):
            # fit to a cubic polynomial
            return a + b*x + c*x**2 + d*x**3
        
        # curve fit return optimal Parameters and 
        # the covariance of those parameters
        poptN, pcovN = curve_fit(f, self.N.xval, self.N.yval)
        poptS, pcovS = curve_fit(f, self.S.xval, self.S.yval)
        
        x1 = np.linspace(self.N.xval[0], self.N.xval[-1])
        x2 = np.linspace(self.S.xval[0], self.S.xval[-1])
        
        plt.plot(x1, f(x1, *poptN), color='green')
        plt.plot(x2, f(x2, *poptS), color='magenta')
        plt.draw()
        
        
        # split horizontally 
        psiN = self.N.p[0].psi(grid)
        psiS = self.S.p[-1].psi(grid)
        
        alp = .5 # test split at half psi
        
        psiAlp = psiS + alp*(psiN - psiS)
                
        # TODO in the Ingrid.construct_SNL_patches method there must be some
        # inconsistency in the definition of lines, on the order isn't being 
        # maintained, because the below definition for the west endpoints
        # of the north and south lines is correct for most of the patches
        # but a few has one point on the wrong end.
        # these are: IDL, IST, OCB, and OPF
        x1 = self.S.p[-1].x
        x2 = self.N.p[0].x
        y1 = self.S.p[-1].y
        y2 = self.N.p[0].y
        
        plt.plot(x2, y2, 'X', color='blue')  # north
        plt.plot(x1, y1, 'X', color='red')  # south
        plt.draw()
                
        def fpsi(x):
            # line in 3d space
            # must manually calculate y each time we stick it into
            # the line of interest
            y = (y2-y1)/(x2-x1)*(x-x1)+y1
            return grid.psi_norm.get_psi(x, y) - psiAlp

        sol = root_scalar(fpsi, bracket=[x1, x2])
        r_psi = sol.root
        z_psi = (y2-y1)/(x2-x1)*(r_psi-x1)+y1
        
        plt.plot(r_psi, z_psi, 'x')
        plt.draw()
        
        mid_line = grid.eq.draw_line((r_psi, z_psi), {'line': self.E},option='theta', direction='cw', show_plot=True)
        
        self.lines.append(mid_line)
 

def calc_mid_point(v1, v2):
    """
    Calculates the bisection of two vectors of equal length, 
    and returns the point on the circle at that angle.
    
    
    Parameters
    ----------
    v1 : geometry.Vector
        v1 must be furthest right in a counter clockwise direction.
    v2 : geometry.Vector
        Vector on the left.
        
    Returns
    -------
    tuple
        The point at the bisection of two vectors.
    
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
    """
    Check if two points are on opposite sides of a given line.
        
    Parameters
    ----------
    p1 : tuple
        First point, (x, y)
    p2 : tuple
        Second point, (x, y)
    line : array-like
        The line is comprised of two points ((x, y), (x, y)).

    Returns
    -------
    tuple
        Returns two numbers, if the signs are different 
        the points are on opposite sides of the line.
    
    """    
    (x1, y1), (x2, y2) = line
    x, y = p1
    a, b = p2
    d1 = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)
    d2 = (a-x1)*(y2-y1)-(b-y1)*(x2-x1)
    
    return np.sign(d1),  np.sign(d2)


def intersect(line1, line2):
    """ Finds the intersection of two line segments
    
    
    Parameters
    ----------
    line1 : array-like
    line2 : array-like
        Both lines of the form A = ((x, y), (x, y)).
        
    Returns
    -------
    tuple
        Coordinates of the intersection.
    
    """ 
    # TODO: Imporove this method so it can handle verticle lines.
    # there is a division by zero error in some DIII-D data caused by this.
    def line(x, line):
        """ Point slope form. """
        (x1, y1), (x2, y2) = line
        if x2-x1 == 0:
            x1 += 1e-4
        
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

def segment_intersect(line1, line2):
    """ Finds the intersection of two FINITE line segments.

    Parameters
    ----------
    line1 : array-like
    line2 : array-like
        Both lines of the form line1 = (P1, P2), line2 = (P3, P4)

    Returns
    -------
    bool, tuple
        True/False of whether the segments intersect
        Coordinates of the intersection
    """

    (xa, ya), (xb, yb) = (line1.xval[0], line1.yval[0]), (line2.xval[0], line2.yval[0])
    (xc, yc), (xd, yd) = (line1.xval[1], line1.yval[1]), (line2.xval[1], line2.yval[1])


    print(xb - xa)
    print(yb - ya)
    M = np.array([[xb - xa, -xd + xc],\
                 [yb - ya, -yd + yc]])

    print(M)
    r = np.array([xc-xa, yc - ya])

    print(r)

    sol = np.linalg.solve(M, r)

    if sol[0] <= 1 and sol[1] <= 1\
       and sol[0] >= 0 and sol[1] >= 0:
           return True, sol
    else:
        False, (None,None)


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
    patch.plot_border('red')

    plt.xlim(0, 6)
    plt.ylim(0, 6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title('Example Patch')
    plt.show()
