#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 16:00:04 2019

@author: watkins35, garcia299
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import linecache
from matplotlib.patches import Polygon
from scipy.optimize import fsolve, curve_fit, root_scalar, brentq
from scipy.interpolate import splprep, splev
def DrawLine(data):
    (i,W_vertices,E,dynamic_step,eq,itot)=data   
    #if Verbose: 
    #return DrawLineBase(eq,W_vertices, {'line' : E}, option = {}, \
    #direction = 'cw', show_plot = 0, dynamic_step = dynamic_step)
    print('Tracing Radial Line in the poloidal direction {}/{}'.format(i,itot))
    return eq.draw_line(eq,W_vertices, {'line' : E},
                  option={}, direction='cw', show_plot=0, dynamic_step=dynamic_step) 
def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def which_non_increasing(L):
    return [i for i,(x, y) in enumerate(zip(L, L[1:])) if x>y]

def which_increasing(L):
    return [i for i,(x, y) in enumerate(zip(L, L[1:])) if x<y]

def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))

def find_split_index(split_point, line):
    same_line_split = False
    for i in range(len(line.p) - 1):
        # Split point is exactly on a Line object's point. Occurs often
        # when splitting a Line object with itself.

        if split_point.y == line.p[i].y and split_point.x == line.p[i].x:
            same_line_split = True
            return i, same_line_split
        # Create two vectors.
        end_u = np.array([line.p[i+1].x - line.p[i].x, line.p[i+1].y - line.p[i].y])
        split_v = np.array([split_point.x - line.p[i].x, split_point.y - line.p[i].y])

        if is_between(end_u, split_v):
            # store index corresponding to the start of the segment containing the split_point.
            return i, same_line_split
        else:
            continue
    return None, same_line_split

def is_between(end_u, split_v):
    eps = 1e-9
    # check cross product vector norm against eps.
    if np.linalg.norm(np.cross(end_u, split_v)) < eps:
        # Colinear up to eps.
        # Ensure dot product is positive and vector v lies in distance of u norm.
        if (np.dot(end_u, split_v) > 0) and (np.linalg.norm(end_u) > np.linalg.norm(split_v)):
            return True
    else:
        return False

def rotmatrix(theta):
    rot = np.zeros((2, 2))
    rot[0, 0] = np.cos(theta)
    rot[1, 1] = np.cos(theta)
    rot[0, 1] = -np.sin(theta)
    rot[1, 0] = np.sin(theta)
    return rot

def rotate(vec, theta, origin):
    return np.matmul(rotmatrix(theta), vec - origin) + origin

def unit_vector(v):
    """
    Returns unit vector
    """
    return v / np.linalg.norm(v)

def angle_between(u, v, origin, relative = False):
    """
    Compute angle in radians between vectors u and v
    """
    u_norm = unit_vector(u-origin)
    v_norm = unit_vector(v-origin)
    angle = np.arccos(np.clip(np.dot(u_norm, v_norm), -1, 1))
    return orientation_between(u, v, origin) * angle if relative else angle

def orientation_between(u, v, origin):
    """
    Compute angle in radians between vectors u and v
    """
    u_norm = unit_vector(u-origin)
    v_norm = unit_vector(v-origin)
    return np.sign(np.arctan2(u_norm[0] * v_norm[1] - u_norm[1] * v_norm[0], \
        u_norm[0] * v_norm[0] + u_norm[0] * v_norm[1]))


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
            self.coor = (self.x, self.y)
        elif np.shape(pts) == (1, 2):
            self.x, self.y = pts[0]
            self.coor = (self.x, self.y)
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
        if len(points)>1:
            self.X=self.xval[1]-self.xval[0]
            self.Y=self.yval[1]-self.yval[0]
        else:
            self.X=0
            self.Y=0
    def Norm(self):
        """
        Return norm of the lines

        Returns:
            None.

        """
        
        return np.sqrt(self.X**2+self.Y**2)
       
            
    def GetAngle(self,Line):
        """
        Return the angle between two lines in degree (between 0 and 180 degrees) 

        Args:
            Line (TYPE): DESCRIPTION.

        Returns:
            None.

        """
        if Line.Norm()!=0 and self.Norm()!=0:
            return np.rad2deg(np.arccos(((self.X)*(Line.X)+(self.Y)*(Line.Y))/(Line.Norm()*self.Norm())))
        else:
            return None
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

    def copy(self):
        """
        Returns a copy of the line.
        """
        return Line(self.p[::])

    def reverse_copy(self):
        """
        Returns a copy of the line in the reversed direction.
        Does not overwrite the current line.
        """

        return Line(self.p[::-1])

    def straighten(self):
        """
        Returns a Line instance consisting of the caller's
        first point and final point. To be used with 'curves'
        in order to generate chords.
        """

        return Line([self.p[0], self.p[-1]])

    def plot(self, color='#1f77b4', label=None):
        """ Plots the line of the current figure.

        Parameters
        ----------
        color : str, optional
            Defaults to a light blue.
        """

        plt.plot(self.xval, self.yval, color=color, zorder = 5, label=label)


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

    def fluff(self, num = 1000,verbose=False):
        if verbose: print('# Fluffing with n=',1000)
        x_fluff = np.empty(0)
        y_fluff = np.empty(0)
        if verbose: print('# fluff: len(self.xval)=',len(self.xval))
        for i in range(len(self.xval) - 1):
            x_fluff = np.append(x_fluff, np.linspace(self.xval[i], self.xval[i+1], num, endpoint = False))
            y_fluff = np.append(y_fluff, np.linspace(self.yval[i], self.yval[i+1], num, endpoint = False))
        x_fluff = np.append(x_fluff, self.xval[-1])
        y_fluff = np.append(y_fluff, self.yval[-1])

        return x_fluff, y_fluff

    def fluff_copy(self, num = 5):
        pts = []
        xf, yf = self.fluff(num)
        for i in zip(xf, yf):
            pts.append(Point(i))
        return Line(pts)


    def split(self, split_point, add_split_point = False):
        """
        Split a line object into two line objects at a particular point.
        Returns two Line objects Segment A and Segment B (corresponding to both subsets of Points)
        split_point: Point object
            - Point that determines splitting location.
            - split_point is always included in Segment B
        add_split_point: Boolean
            - Append the split point to Segment A.
        """
        d_arr = []

        # ind was not defined if the conditions are not met, e.g. when magnetic field lines do not intersect the targets
        # Safety check added with an exception
        ind, same_line_split = find_split_index(split_point, self)
        if ind==None:
            for i in range(len(self.p) - 1):
                plt.plot(self.p[i].x,self.p[i].y,'.',color='black',ms=8)
                plt.plot(split_point.x,split_point.y,'s',color='red',ms=8)
            plt.show()
            raise ValueError("Value of ind has not been set. This probably means that one of the targets does not intersect one of the field lines.\
                             Check that targets are wide enough to emcompass the SOL and PFR region")

        start__split = self.p[:ind + (1 if not same_line_split else 0)]
        start__split += [split_point] if add_split_point else []
        split__end = [split_point] + self.p[ind + (1 if not same_line_split else 2):]

        return Line(start__split), Line(split__end)


        """
        for i in range(len(self.p) - 1):
            v1 = np.array([self.p[i+1].x - self.p[i].x, self.p[i+1].y - self.p[i].y])
            v2 = np.array([split_point.x - self.p[i].x, split_point.y - self.p[i].y])

            # Split point is exactly on a Line object's point. Occurs often
            # when splitting a Line object with itself.
            if split_point.y == self.p[i].y and split_point.x == self.p[i].x:
                same_line_split = True
                d_arr.append((0.0, i))
                break
            # Vertical line segment.
            elif self.p[i+1].x == self.p[i].x:
                    # Append x-direction offset and the corresponding i index.
                    d_arr.append((np.abs(split_point.x - self.p[i].x), i))
                else:
                    continue
            # Horizontal line segment.
            elif self.p[i+1].y == self.p[i].y:
                    # Append y-direction offset and the corresponding i index.
                    d_arr.append((np.abs(split_point.y - self.p[i].y), i))
                else:
                    continue
            # Sloped line segment.
            else:
                # Compute slope and check if point lies on segment. Append difference and index.
                m = (self.p[i+1].y - self.p[i].y)/(self.p[i+1].x - self.p[i].x)
                y = self.p[i+1].y - split_point.y
                x = m * (self.p[i].x - split_point.x)
                d_arr.append((np.abs(x-y), i))
        """
        # Get index of minimal value.
        ind = np.asarray([dist[0] for dist in d_arr]).argmin()
        # Recover the index location we stored earlier corresponding to the target segment.
        """
        end_dist = np.sqrt((self.p[ind + 1].x - split_point.x) ** 2 + (self.p[ind + 1].y - split_point.y) ** 2)
        start_dist = np.sqrt((self.p[ind].x - split_point.x) ** 2 + (self.p[ind].y - split_point.y) ** 2)
        if start_dist == end_dist:
            same_line_split = True
        elif end_dist < start_dist:
            ind += -1

        start__split = self.p[:ind + (1 if not same_line_split else 0)]
        split__end = [split_point] + new_line.p[ind + (1 if not same_line_split else 2):]

        start__split += [split_point] if add_split_point else []
        return Line(start__split), Line(split__end)
        """
    def points(self):
        """ Returns the points in the line as a tuple. """
        return [(p.x, p.y) for p in self.p]


class Cell:
    """
    Each Cell is contained within a patch

    Parameters:
    -----------
    vertices : geometry.Point
        Each cell is defined by four points in a clockwise order.
        NW -> NE -> SE -> SW
    """
    def __init__(self, lines):
        self.lines = lines
        N = lines[0]
        S = lines[1]
        E = lines[2]
        W = lines[3]

        self.vertices = {'NW' : N.p[0], 'NE' : N.p[-1], \
                         'SW' : S.p[0], 'SE' : S.p[-1]}

        self.center = Point((np.mean([p.x for p in [self.vertices[coor] for coor in ['NW', 'NE', 'SE', 'SW']]]),\
                np.mean([p.y for p in [self.vertices[coor] for coor in ['NW', 'NE', 'SE', 'SW']]])))

        self.vertices.update({'CENTER' : self.center})

        self.p = list(N.p) + list(S.p)

    def plot_border(self, color = 'red'):
        plt.plot([self.vertices['NW'].x, self.vertices['NE'].x], \
                [self.vertices['NW'].y, self.vertices['NE'].y], \
                linewidth = 1, color = color,label='cell')

        plt.plot([self.vertices['NE'].x, self.vertices['SE'].x], \
                [self.vertices['NE'].y, self.vertices['SE'].y], \
                linewidth = 1, color = color,label='cell')

        plt.plot([self.vertices['SE'].x, self.vertices['SW'].x], \
                [self.vertices['SE'].y, self.vertices['SW'].y], \
                linewidth = 1, color = color,label='cell')

        plt.plot([self.vertices['SW'].x, self.vertices['NW'].x], \
                [self.vertices['SW'].y, self.vertices['NW'].y], \
                linewidth = 1, color = color,label='cell')

    def CollectBorder(self, color = 'red'):
        segs[0,0,0]=self.vertices['NW'].x
        segs[0,1,0]=self.vertices['NW'].y
        segs[1,0,0]=self.vertices['NE'].x
        segs[1,1,0]=self.vertices['NE'].y



        plt.plot([self.vertices['NE'].x, self.vertices['SE'].x], \
                [self.vertices['NE'].y, self.vertices['SE'].y], \
                linewidth = 1, color = color)

        plt.plot([self.vertices['SE'].x, self.vertices['SW'].x], \
                [self.vertices['SE'].y, self.vertices['SW'].y], \
                linewidth = 1, color = color)

        plt.plot([self.vertices['SW'].x, self.vertices['NW'].x], \
                [self.vertices['SW'].y, self.vertices['NW'].y], \
                linewidth = 1, color = color)

    def plot_center(self, color = 'black'):
        plt.plot(self.center.x, self.center.y, '.', markersize = 1, color = color,label='cell')
        """
        Line([self.vertices['NW'], self.vertices['NE']]).plot(color)
        Line([self.vertices['NE'], self.vertices['SE']]).plot(color)
        Line([self.vertices['SE'], self.vertices['SW']]).plot(color)
        Line([self.vertices['SW'], self.vertices['NW']]).plot(color)
        """

    def fill(self, color = 'salmon'):
        x = [p.x for p in self.p]
        y = [p.y for p in self.p]
        plt.fill(y, x, facecolor = color)


class Patch:
    """
    Each patch contains a grid, and has it's own shape.

    Parameters
    ----------
    lines : geometry.Line
        Each patch defined by four lines in order -
        Nline, Eline, Sline, Wline - order points to go clockwise.

    """

    def __init__(self, lines, patchName = '', PatchTagMap=None, platePatch = False, plateLocation = None, color='blue'):
        self.lines = lines
        self.N = lines[0]
        self.E = lines[1]
        self.S = lines[2]
        self.W = lines[3]
        self.BoundaryPoints={}

        # This is the border for the fill function
        # It need to only include N and S lines
        self.p = list(self.N.p) + list(self.E.p) + list(self.S.p) + list(self.W.p)
        self.PatchTagMap = PatchTagMap
        self.platePatch = platePatch
        self.plateLocation = plateLocation
        self.patchName = patchName
        self.color=color
        self.Verbose=True

        # TODO: Populate these attributes when generating subgrid.
        #       Currently this is being done in the concat_grid method.
        self.npol = 0
        self.nrad = 0

        self.PatchLabelDoc={
            'I': 'Inner',
            'O':' Outer',
            'D':'Divertor',
            'L':'Leg',
           'P': 'Private',
           'F':'Flux',
            'T': 'Top',
            'B': 'Bottom',
            'S': 'SOL',
            'C': 'Core'
            }

    def plot_border(self, ax=None,color='red'):
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
        x = np.array([p.x for p in self.p])
        y = np.array([p.y for p in self.p])
        arr = np.column_stack((x,y))
        PatchLabel=self.patchName+' (' +UnfoldLabel(self.PatchLabelDoc,self.patchName)+')'
        patch = Polygon(arr, fill = True, closed = True, color = color,label=PatchLabel)
        ax = plt.gca()
        ax.add_patch(patch)
        ax.plot()

        plt.show()

    def plot_subgrid(self, ax=None,color = 'blue'):
        for row in self.cell_grid:
            for cell in row:
                cell.plot_border(color)
        plt.pause(0.1)


    def adjust_corner(self, point, corner):
        if corner == 'NW':
            self.cell_grid[0][0].vertices[corner] = point
            self.cell_grid[0][0].vertices[corner] = point
        elif corner == 'NE':
            self.cell_grid[0][-1].vertices[corner] = point
            self.cell_grid[0][-1].vertices[corner] = point
        elif corner == 'SE':
            self.cell_grid[-1][-1].vertices[corner] = point
            self.cell_grid[-1][-1].vertices[corner] = point
        elif corner == 'SW':
            self.cell_grid[-1][0].vertices[corner] = point
            self.cell_grid[-1][0].vertices[corner] = point

    def get_tag(self):
        name = self.patchName
        return self.PatchTagMap[name] if len(name) == 3 else self.PatchTagMap[name[1:]]

    def make_subgrid(self, grid, np_cells = 2, nr_cells = 2, _poloidal_f=lambda x:x, _radial_f=lambda x:x,verbose = False, visual = False,ShowVertices=False,OptionTrace='theta'):
        """
        Generate a refined grid within a patch.
        This 'refined-grid' within a Patch is a collection
        of num x num Cell objects

        Parameters:
        ----------
        grid : Ingrid object
                To be used for obtaining Efit data and all
                other class information.
        num  : int, optional
                Number to be used to generate num x num
                cells within our Patch.
        """

        def psi_parameterize(grid, r, z):
            """
            Helper function to be used to generate a
            list of values to parameterize a spline
            in Psi. Returns a list to be used for splprep only
            """
            vmax = grid.psi_norm.get_psi(r[-1], z[-1])
            vmin = grid.psi_norm.get_psi(r[0], z[0])

            vparameterization = np.empty(0)
            for i in range(len(r)):
                vcurr = grid.psi_norm.get_psi(r[i], z[i])
                vparameterization = np.append(vparameterization, abs((vcurr - vmin) / (vmax - vmin)))

            return vparameterization

        def psi_test(leg, grid):
            """
            Determine whether a leg's orientation is increasing in psi
            by checking the endpoints of a given leg.

            Returns True if final Point entry is greater in psi than the
            first Point entry. Returns False otherwise.
            """
            p1 = leg.p[0]
            p2 = leg.p[-1]

            return True if grid.psi_norm.get_psi(p2.x, p2.y) > grid.psi_norm.get_psi(p1.x, p1.y) else False

        def set_dynamic_step(ratio = 0.01):
            """
            Compute a dt value for line_tracing that occurs during subgrid generation of
            this Patch object. The value is taken to be a ratio of the average length of
            a Patch in the radial direction.
            """
            d = 0
            for i in range(nr_lines):
                d += np.sqrt((E_vertices[i].x - W_vertices[i].x)**2 + (E_vertices[i].y - W_vertices[i].y)**2)
            dynamic_step = d / nr_lines
            print('Dynamic-Step value came out to be: {}\n'.format(dynamic_step * ratio))
            return dynamic_step * ratio


        # Allocate space for collection of cell objects.
        # Arbitrary 2D container for now.
        cell_grid = []


        if verbose:
            print('Constructing grid for patch "{}" with dimensions (np, nr) = ({}, {})'.format(self.patchName, np_cells, nr_cells))
            print(np_cells)
            print(nr_cells)
        np_lines = np_cells + 1
        nr_lines = nr_cells + 1
        if verbose: print(' # Create B-Splines along the North and South boundaries.')
        # Create B-Splines along the North and South boundaries.
        N_vals = self.N.fluff()

        self.N_spl, uN = splprep([N_vals[0], N_vals[1]], s = 0)
        # Reverse the orientation of the South line to line up with the North.

        S_vals = self.S.reverse_copy().fluff()
        self.S_spl, uS = splprep([S_vals[0], S_vals[1]], s = 0)
        if verbose: print(' # Create B-Splines along West boundaries.')
        # Create B-Splines along the East and West boundaries.
        # Parameterize EW splines in Psi
        try:
            #Cannot fluff with too many points
            n = 500 if len(self.W.p) < 50 else 100
           # W_vals = self.W.reverse_copy().fluff(num = n)
            W_vals = self.W.reverse_copy().fluff(n,verbose=verbose)
            Psi=psi_parameterize(grid, W_vals[0], W_vals[1])
            self.W_spl, uW = splprep([W_vals[0], W_vals[1]], u = Psi, s = 0)
        except Exception as e:
            exc_type, exc_obj, tb = sys.exc_info()
            f = tb.tb_frame
            lineno = tb.tb_lineno
            filename = f.f_code.co_filename
            linecache.checkcache(filename)
            line = linecache.getline(filename, lineno, f.f_globals)
            print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

        if verbose: print(' # Create B-Splines along the East boundaries.')
        try :
            n = 500 if len(self.E.p) < 50 else 100
            E_vals = self.E.fluff(num = n)
            self.E_spl, uE = splprep([E_vals[0], E_vals[1]], u = psi_parameterize(grid, E_vals[0], E_vals[1]), s = 10)
        except Exception as e:
            print(' Number of points on the boundary:', len(self.E.p))
            plt.plot(E_vals[0],E_vals[1],'.',color='black')
            print(repr(e))
        if verbose: print(' #check platePatch')
        # ACCURACY ON PSI_MIN SIDE OF W IDL LINE ISSUE.

        if self.platePatch:

            def f(u, *args):
                _y = splev(u, args[0])
                return grid.psi_norm.get_psi(args[1], args[2]) - grid.psi_norm.get_psi(_y[0], _y[1])

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return array[idx]

            if self.plateLocation == 'W':
                _u = uW
                U_vals = W_vals
                U_spl = self.W_spl
                plate_north = [N_vals[0][0], N_vals[1][0]]
                plate_south = [S_vals[0][0], S_vals[1][0]]
            elif self.plateLocation == 'E':
                _u = uE
                U_vals = E_vals
                U_spl = self.E_spl
                plate_north = [N_vals[0][-1], N_vals[1][-1]]
                plate_south = [S_vals[0][-1], S_vals[1][-1]]

            if verbose:
                print('=' * 80 + '\n')
                print('Parameterization of Target Plate in PSI:')
                print(_u)
                print('=' * 80 + '\n')

            lookup = {}
            for i in range(len(_u)):
                lookup[_u[i]] = i
            try:
                plate_north_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args = (U_spl, plate_north[0], plate_north[1])))]
            except ValueError:
                print('NorthIndex: brentq failed. attempting fsolve...')
                try:
                    plate_north_index = lookup[find_nearest(_u, fsolve(f, 0, args = (U_spl, plate_north[0], plate_north[1])))]
                    print('NorthIndex: ...fsolve success!')
                except:
                    print('NorthIndex: ERROR IN PARAMETERIZATION IN PSI')
            try:
                plate_south_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args = (U_spl, plate_south[0], plate_south[1])))]
            except ValueError:
                print('SouthIndex: brentq failed. attempting fsolve...')
                try:
                    plate_south_index = lookup[find_nearest(_u, fsolve(f, 0, args = (U_spl, plate_south[0], plate_south[1])))]
                    print('SouthIndex: ...fsolve success!')
                except:
                    print('SouthIndex: ERROR IN PARAMETERIZATION IN PSI')
            if plate_south_index > plate_north_index:
                plate_north_index, plate_south_index = plate_south_index, plate_north_index

            U_vals = [U_vals[0][plate_south_index:plate_north_index+1], U_vals[1][plate_south_index:plate_north_index+1]]
            U_spl, _u = splprep([U_vals[0], U_vals[1]], u = psi_parameterize(grid, U_vals[0], U_vals[1]), s = 0)

            if self.plateLocation == 'W':
                W_vals = U_vals
                self.W_spl = U_spl
                uW = _u
            elif self.plateLocation == 'E':
                E_vals = U_vals
                self.E_spl = U_spl
                uE = _u

        # Generate our sub-grid anchor points along the North
        # and South boundaries of our patch.
        self.N_vertices = []
        self.S_vertices = []
        E_vertices = []
        W_vertices = []
        
        if verbose: print('# Generate our sub-grid anchor points along the North and South boundaries of our patch.')
        # and South boundaries of our patch')
        
        if self.BoundaryPoints.get('N') is None:
            for i in range(np_lines):
                _n = splev(_poloidal_f(i / (np_lines-1)), self.N_spl)
                self.N_vertices.append(Point((_n[0], _n[1])))
        else:
            if self.Verbose:print('Find boundary points at face "N" for {}:{}'.format(self.patchName,self.BoundaryPoints.get('N')))
            self.N_vertices=self.BoundaryPoints.get('N')  
            
        if self.BoundaryPoints.get('S') is None:
            for i in range(np_lines):
                _s = splev(_poloidal_f(i / (np_lines-1)), self.S_spl)
                self.S_vertices.append(Point((_s[0], _s[1])))
        else:
            self.S_vertices=self.BoundariesPoints.get('S')    

        u=[_radial_f(i / (nr_lines-1)) for i in range(nr_lines)]

        Npts=1000
        xy = splev(np.linspace(0,1,Npts), self.E_spl)

        for i in range(nr_lines):
            _e = splev(u[i], self.E_spl)

            E_vertices.append(Point((_e[0], _e[1])))

            _w = splev(u[i], self.W_spl)
            W_vertices.append(Point((_w[0], _w[1])))
        DetailedVertices=False
        if  not isinstance(ShowVertices,bool):
            if ShowVertices=='detailed':
                ShowVertices=True
                DetailedVertices=True
            else:
                raise ValueError('Unknow type for ShowVertices: must be True,False or "detailed"',ShowVertices)

        if visual or ShowVertices:
            if not DetailedVertices:
                markers=['o']*4
                ms=8
                color='black'

            else:
                color=self.color
                ms=6
                markers=['o','X','s','D']
            for vertices,mark in zip([W_vertices, E_vertices, self.N_vertices, self.S_vertices],markers):
                for p in vertices:
                    plt.plot(p.x, p.y, '.', color = color, markersize = ms,marker=mark,markeredgecolor='black')
        # Radial lines of Psi surfaces. Ordered with increasing magnitude, starting with
        # the South boundary of the current Patch, and ending with the North boundary of
        # this current Patch. These will serve as delimiters when constructing cells.
        radial_lines = [self.N]
        radial_vertices = [self.N_vertices]
        dynamic_step = set_dynamic_step()
        if verbose: print('# Interpolate radial lines between North and South patch boundaries..')
        # Interpolate radial lines between North and South patch boundaries.
        self.radial_spl=[]
            
        for i in range(len(W_vertices) - 2):
            #TODO: parallelize tracing of radial lines (draw_line function must be "externalized" in the scope of the script)
            radial_lines.append(grid.eq.draw_line(W_vertices[i + 1], {'line' : self.E}, option = OptionTrace, \
                direction = 'cw', show_plot = 0, dynamic_step = dynamic_step))
            radial_vals = radial_lines[i + 1].fluff(1000)
            Radial_spl, uR = splprep([radial_vals[0], radial_vals[1]], s = 0)
            self.radial_spl.append(Radial_spl)
            vertex_list = []
            for j in range(np_lines):
                u=_poloidal_f(j / (np_lines-1))
                _r = splev(u, self.radial_spl[i])
                Pt=Point((_r[0], _r[1]))
                if self.CorrectDistortion['Active'] and j>0 and j<np_lines-1:
                    Res=self.CorrectDistortion['Resolution']
                    ThetaMin=self.CorrectDistortion['ThetaMin']
                    ThetaMax=self.CorrectDistortion['ThetaMax']
                    umin=_poloidal_f((j-1) / (np_lines-1))
                    umax=_poloidal_f((j+1) / (np_lines-1))
                    Pt1=radial_vertices[i][j]
                    Pt2=radial_vertices[i][j-1]
                    Tag='---- Correcting points: {},{}'.format(i, j)
                    Pt=CorrectDistortion(u,Pt,Pt1,Pt2,self.radial_spl[i],ThetaMin,ThetaMax,umin,umax,Res,visual,Tag,self.Verbose)
                     
                vertex_list.append(Pt)    
                if visual:
                    for p in vertex_list:
                        plt.plot(p.x, p.y, '.', color = 'black', markersize = 8)
                        
            radial_vertices.append(vertex_list)
        radial_lines.append(self.S)
        
        #Correct point on south boundary
        if self.CorrectDistortion['Active']:
            for j in range(1,np_lines-1):
                u=_poloidal_f(j / (np_lines-1))
                Pt=self.S_vertices[j]
                Res=self.CorrectDistortion['Resolution']
                ThetaMin=self.CorrectDistortion['ThetaMin']
                ThetaMax=self.CorrectDistortion['ThetaMax']
                umin=_poloidal_f((j-1) / (np_lines-1))
                umax=_poloidal_f((j+1) / (np_lines-1))
                Pt1=radial_vertices[-1][j]
                Pt2=radial_vertices[-1][j-1]
                Tag='---- Correcting south boundary points:{}'.format(j)
                self.S_vertices[j]=CorrectDistortion(u,Pt,Pt1,Pt2,self.S_spl,ThetaMin,ThetaMax,umin,umax,Res,visual,Tag,self.Verbose)
        radial_vertices.append(self.S_vertices)
        if verbose: print('# Create Cells: South boundary -> North Boundary')
        # Create Cells: South boundary -> North Boundary
        for i in range(len(radial_lines)):
            if radial_lines[i] is self.S:
                break
            cell_grid.append([])
            for j in range(len(radial_vertices[i]) - 1):
                NW = radial_vertices[i][j]
                NE = radial_vertices[i][j+1]
                SW = radial_vertices[i+1][j]
                SE = radial_vertices[i+1][j+1]
                cell_grid[i].append(Cell([Line([NW, NE]), Line([SW, SE]), Line([SE, NE]), Line([SW, NW])]))

        self.cell_grid = cell_grid

    def CheckPatch(self,grid,verbose=False):
        if verbose: print(' # Checking if patch boundaries can be interpolated wiht splines')

        def psi_parameterize(grid, r, z):
            """
            Helper function to be used to generate a
            list of values to parameterize a spline
            in Psi. Returns a list to be used for splprep only
            """
            vmax = grid.psi_norm.get_psi(r[-1], z[-1])
            vmin = grid.psi_norm.get_psi(r[0], z[0])

            vparameterization = np.empty(0)
            for i in range(len(r)):
                vcurr = grid.psi_norm.get_psi(r[i], z[i])
                vparameterization = np.append(vparameterization, abs((vcurr - vmin) / (vmax - vmin)))

            return vparameterization

        def IsMonotonic(psi,x,y,Str):
            if not (non_increasing(psi) or non_decreasing(psi)):
                print(Str +' is not monotonic')
                d=which_non_increasing(psi)
                u=which_increasing(psi)
                plt.plot(x[u],y[u],'o','r')
                plt.plot(x[d],y[d],'o','b')
                raise ValueError(Str+ ' is not monotonic')
            else:
               print(Str +' is monotonic')

       # N_vals = self.N.fluff()
        #S_vals = self.S.reverse_copy().fluff()
        n = 20 if len(self.W.p) > 500 else 100
        W_vals = self.W.reverse_copy().fluff(n)
        n = 20 if len(self.E.p) > 500 else 100
        E_vals = self.E.fluff(num = n)
        if verbose: print(' ## Getting Psi values along the boundaries')
        #PsiN=psi_parameterize(grid, N_vals[0], N_vals[1])
        #PsiS=psi_parameterize(grid, S_vals[0], S_vals[1])
        PsiW=psi_parameterize(grid, W_vals[0], W_vals[1])
        PsiE=psi_parameterize(grid, E_vals[0], E_vals[1])
        if verbose: print(' ## Checking monoticity of Psi values along the boundaries')
        IsMonotonic(PsiW,W_vals[0],W_vals[1],'PsiW')
        IsMonotonic(PsiE,E_vals[0],E_vals[1],'PsiE')
        
def CorrectDistortion(u,Pt,Pt1,Pt2,spl,ThetaMin,ThetaMax,umin,umax,Resolution,visual,Tag,MinTol=1.02,MaxTol=0.98,Verbose=False):
            dumax=(umax-u)/Resolution
            dumin=(u-umin)/Resolution
            Theta=Line([Pt1,Pt]).GetAngle(Line([Pt1,Pt2]))
            if Theta<ThetaMin or Theta>ThetaMax:
                if Verbose: print('{}: u={};Theta={};ThetaMin={};ThetaMax={}'.format(Tag,u,Theta,ThetaMin,ThetaMax))
                if visual:
                        plt.plot(Pt.x, Pt.y, '.', color = 'red', markersize = 8,marker='o')
                        plt.show()
                        plt.draw()
                icount=0
                color='purple'
                while Theta<ThetaMin or Theta>ThetaMax:
                    icount+=1    
                    if Theta<ThetaMin:
                       u=u+dumax
                    elif Theta>ThetaMax:    
                       u=u-dumin
                    if u>umax*MaxTol or u<umin*MinTol:
                       if Verbose: print('[{}]>>>> umax={} umin={};u:{};Theta={}'.format(icount,umax,umin,u,Theta))
                       color='gray'
                       break
                    _r = splev(u, spl)
                    Pt=Point((_r[0], _r[1]))
                    Theta=Line([Pt1,Pt]).GetAngle(Line([Pt1,Pt2])) 
                if Verbose: print('[{}]>>>> u:{};Theta={}'.format(icount,u,Theta))
                if visual:
                    plt.plot(Pt.x, Pt.y, '.', color = color , markersize = 8,marker='s')
            return(Pt)         

  

        # ACCURACY ON PSI_MIN SIDE OF W IDL LINE ISSUE.
class DNL_Patch(Patch):
    def __init__(self, lines, patchName = '', platePatch = False, plateLocation = None):
        super().__init__(lines, patchName, platePatch, plateLocation)

    def make_subgrid(self, grid, np_cells = 2, nr_cells = 2, verbose = False, visual = False,ShowVertices=False):
        """
        Generate a refined grid within a patch.
        This 'refined-grid' within a Patch is a collection
        of num x num Cell objects

        Parameters:
        ----------
        grid : Ingrid object
                To be used for obtaining Efit data and all
                other class information.
        num  : int, optional
                Number to be used to generate num x num
                cells within our Patch.
        """

        def get_func(grid, _func):

            def make_sympy_func(var, expression):
                import sympy as sp
                _f = sp.lambdify(var, expression, 'numpy')
                return _f

            f_str_raw = _func

            f_str_raw = f_str_raw.replace(' ', '')
            delim = f_str_raw.index(',')

            var = f_str_raw[0 : delim]
            expression = f_str_raw[delim + 1 :]

            _func = make_sympy_func(var, expression)

            return _func

        def psi_parameterize(grid, r, z):
            """
            Helper function to be used to generate a
            list of values to parameterize a spline
            in Psi. Returns a list to be used for splprep only
            """
            vmax = grid.psi_norm.get_psi(r[-1], z[-1])
            vmin = grid.psi_norm.get_psi(r[0], z[0])

            vparameterization = np.empty(0)
            for i in range(len(r)):
                vcurr = grid.psi_norm.get_psi(r[i], z[i])
                vparameterization = np.append(vparameterization, abs((vcurr - vmin) / (vmax - vmin)))

            return vparameterization

        def psi_test(leg, grid):
            """
            Determine whether a leg's orientation is increasing in psi
            by checking the endpoints of a given leg.

            Returns True if final Point entry is greater in psi than the
            first Point entry. Returns False otherwise.
            """
            p1 = leg.p[0]
            p2 = leg.p[-1]

            return True if grid.psi_norm.get_psi(p2.x, p2.y) > grid.psi_norm.get_psi(p1.x, p1.y) else False

        def set_dynamic_step(ratio = 0.01):
            """
            Compute a dt value for line_tracing that occurs during subgrid generation of
            this Patch object. The value is taken to be a ratio of the average length of
            a Patch in the radial direction.
            """
            d = 0
            for i in range(nr_lines):
                d += np.sqrt((E_vertices[i].x - W_vertices[i].x)**2 + (E_vertices[i].y - W_vertices[i].y)**2)
            dynamic_step = d / nr_lines
            print('Dynamic-Step value came out to be: {}\n'.format(dynamic_step * ratio))
            return dynamic_step * ratio


        # Allocate space for collection of cell objects.
        # Arbitrary 2D container for now.
        cell_grid = []

        if self.platePatch:
            print('Plate patch!')
            if self.patchName[1] == '1':
                key = 'plate_W1'
            elif self.patchName[1] == '8':
                key = 'plate_E1'
            elif self.patchName[1] == '4':
                key = 'plate_E2'
            elif self.patchName[1] == '5':
                key = 'plate_W2'
            np_cells = grid.settings['target_plates'][key]['np_local']
            _poloidal_func = grid.settings['target_plates'][key]['poloidal_f']
            _poloidal_f = get_func(grid, _poloidal_func)
        else:
            _poloidal_f = lambda x: x

        """
        Organize DNL patches into groupings.
        """

        if self in grid.PRIMARY_SOL:
            try:
                _radial_f = grid.settings['grid_params']['grid_generation']['radial_f_primary_sol']
                valid_function = True
                print('SOL radial transformation: "{}"'.format(_radial_f))
            except KeyError:
                valid_function = False
        elif self in grid.SECONDARY_SOL:
            try:
                _radial_f = grid.settings['grid_params']['grid_generation']['radial_f_secondary_sol']
                valid_function = True
                print('CORE radial transformation: "{}"'.format(_radial_f))
            except KeyError:
                valid_function = False
        elif self in grid.CORE:
            try:
                _radial_f = grid.settings['grid_params']['grid_generation']['radial_f_core']
                valid_function = True
                print('CORE radial transformation: "{}"'.format(_radial_f))
            except KeyError:
                valid_function = False
        elif self in grid.PRIMARY_PF:
            try:
                _radial_f = grid.settings['grid_params']['grid_generation']['radial_f_primary_pf']
                valid_function = True
                print('PF radial transformation: "{}"'.format(_radial_f))
            except KeyError:
                valid_function = False
        elif self in grid.SECONDARY_PF:
            try:
                _radial_f = grid.settings['grid_params']['grid_generation']['radial_f_secondary_pf']
                valid_function = True
                print('PF radial transformation: "{}"'.format(_radial_f))
            except KeyError:
                valid_function = False

        _radial_f = get_func(grid, _radial_f) if valid_function else lambda x: x

        if verbose:
            print('Constructing grid for patch "{}" with dimensions (np, nr) = ({}, {})'.format(self.patchName, np_cells, nr_cells))

        np_lines = np_cells + 1
        nr_lines = nr_cells + 1

        # Create B-Splines along the North and South boundaries.
        N_vals = self.N.fluff()

        N_spl, uN = splprep([N_vals[0], N_vals[1]], s = 0)
        # Reverse the orientation of the South line to line up with the North.

        S_vals = self.S.reverse_copy().fluff()
        S_spl, uS = splprep([S_vals[0], S_vals[1]], s = 0)

        # Create B-Splines along the East and West boundaries.
        # Parameterize EW splines in Psi

        n = 50 if len(self.W.p) > 50 else 100
        W_vals = self.W.reverse_copy().fluff(num = n)
        W_spl, uW = splprep([W_vals[0], W_vals[1]], u = psi_parameterize(grid, W_vals[0], W_vals[1]), s = 0)

        n = 50 if len(self.E.p) > 50 else 100
        E_vals = self.E.fluff(num = n)
        E_spl, uE = splprep([E_vals[0], E_vals[1]], u = psi_parameterize(grid, E_vals[0], E_vals[1]), s = 0)

        # ACCURACY ON PSI_MIN SIDE OF W IDL LINE ISSUE.
        if self.platePatch:

            def f(u, *args):
                _y = splev(u, args[0])
                return grid.psi_norm.get_psi(args[1], args[2]) - grid.psi_norm.get_psi(_y[0], _y[1])

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return array[idx]

            if self.plateLocation == 'W':
                _u = uW
                U_vals = W_vals
                U_spl = W_spl
                plate_north = [N_vals[0][0], N_vals[1][0]]
                plate_south = [S_vals[0][0], S_vals[1][0]]
            elif self.plateLocation == 'E':
                _u = uE
                U_vals = E_vals
                U_spl = E_spl
                plate_north = [N_vals[0][-1], N_vals[1][-1]]
                plate_south = [S_vals[0][-1], S_vals[1][-1]]

            if verbose:
                print('=' * 80 + '\n')
                print('Parameterization of Target Plate in PSI:')
                print(_u)
                print('=' * 80 + '\n')

            lookup = {}
            for i in range(len(_u)):
                lookup[_u[i]] = i
            try:
                plate_north_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args = (U_spl, plate_north[0], plate_north[1])))]
            except ValueError:
                print('NorthIndex: brentq failed. attempting fsolve...')
                try:
                    plate_north_index = lookup[find_nearest(_u, fsolve(f, 0, args = (U_spl, plate_north[0], plate_north[1])))]
                    print('NorthIndex: ...fsolve success!')
                except:
                    print('NorthIndex: ERROR IN PARAMETERIZATION IN PSI')
            try:
                plate_south_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args = (U_spl, plate_south[0], plate_south[1])))]
            except ValueError:
                print('SouthIndex: brentq failed. attempting fsolve...')
                try:
                    plate_south_index = lookup[find_nearest(_u, fsolve(f, 0, args = (U_spl, plate_south[0], plate_south[1])))]
                    print('SouthIndex: ...fsolve success!')
                except:
                    print('SouthIndex: ERROR IN PARAMETERIZATION IN PSI')
            if plate_south_index > plate_north_index:
                # print('WARNING: Caught an index error... Fixing...')
                plate_north_index, plate_south_index = plate_south_index, plate_north_index

            U_vals = [U_vals[0][plate_south_index:plate_north_index+1], U_vals[1][plate_south_index:plate_north_index+1]]
            U_spl, _u = splprep([U_vals[0], U_vals[1]], u = psi_parameterize(grid, U_vals[0], U_vals[1]), s = 0)

            if self.plateLocation == 'W':
                W_vals = U_vals
                W_spl = U_spl
                uW = _u
            elif self.plateLocation == 'E':
                E_vals = U_vals
                E_spl = U_spl
                uE = _u

        # Generate our sub-grid anchor points along the North
        # and South boundaries of our patch.
        N_vertices = []
        S_vertices = []
        E_vertices = []
        W_vertices = []

        for i in range(np_lines):
            _n = splev(_poloidal_f(i / (np_lines-1)), N_spl)
            N_vertices.append(Point((_n[0], _n[1])))

            _s = splev(_poloidal_f(i / (np_lines-1)), S_spl)
            S_vertices.append(Point((_s[0], _s[1])))

        for i in range(nr_lines):
            _e = splev(_radial_f(i / (nr_lines-1)), E_spl)
            E_vertices.append(Point((_e[0], _e[1])))

            _w = splev(_radial_f(i / (nr_lines-1)), W_spl)
            W_vertices.append(Point((_w[0], _w[1])))

        if visual or ShowVertices:
            for (vertices,mark) in zip([W_vertices, E_vertices, N_vertices, S_vertices],['o','+','s','d']):
                for p in vertices:
                    plt.plot(p.x, p.y, '.', color = self.color, markersize = 8,marker=mark)
        # Radial lines of Psi surfaces. Ordered with increasing magnitude, starting with
        # the South boundary of the current Patch, and ending with the North boundary of
        # this current Patch. These will serve as delimiters when constructing cells.
        radial_lines = [self.N]
        radial_vertices = [N_vertices]
        dynamic_step = set_dynamic_step()

        # Interpolate radial lines between North and South patch boundaries.
        for i in range(len(W_vertices) - 2):

            radial_lines.append(grid.eq.draw_line(W_vertices[i + 1], {'line' : self.E}, option = OptionTrace, \
                direction = 'ccw', show_plot = visual, dynamic_step = dynamic_step))

            radial_vals = radial_lines[i + 1].fluff()
            radial_spl, uR = splprep([radial_vals[0], radial_vals[1]], s = 0)
            radial_spline = splev(uR, radial_spl)
            vertex_list = []
            for i in range(np_lines):
                _r = splev(_poloidal_f(i / (np_lines-1)), radial_spl)
                vertex_list.append(Point((_r[0], _r[1])))
                if visual:
                    for p in vertex_list:
                        plt.plot(p.x, p.y, '.', color = 'black', markersize = 8,marker='s')
            radial_vertices.append(vertex_list)
        radial_lines.append(self.S)
        radial_vertices.append(S_vertices)

        # Create Cells: South boundary -> North Boundary
        for i in range(len(radial_lines)):
            if radial_lines[i] is self.S:
                break
            cell_grid.append([])
            for j in range(len(radial_vertices[i]) - 1):
                NW = radial_vertices[i][j]
                NE = radial_vertices[i][j+1]
                SW = radial_vertices[i+1][j]
                SE = radial_vertices[i+1][j+1]
                cell_grid[i].append(Cell([Line([NW, NE]), Line([SW, SE]), Line([SE, NE]), Line([SW, NW])]))

        self.cell_grid = cell_grid

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
    for ind in range(len(line) - 1):
        (x1, y1), (x2, y2) = line[ind], line[ind + 1]
        x, y = p1
        a, b = p2
        d1 = (x-x1)*(y2-y1)-(y-y1)*(x2-x1)
        d2 = (a-x1)*(y2-y1)-(b-y1)*(x2-x1)

        if (np.sign(d1) != np.sign(d2)):
            return True
    return False


def intersect(line1, line2, verbose = False):
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
                         y - line(x, test_line)])

    if isinstance(line2, Line):
        line2 = [(p.x, p.y) for p in line2.p]
    for ind in range(len(line2) - 1):
        # use the mean for initial guess
        test_line = [line2[ind], line2[ind + 1]]
        (a, b), (c, d) = line1
        (i, j), (p, q) = test_line
        guess = (np.mean([a, c, i, p]), np.mean([b, d, j, q]))
        sol, infoDict, ier, mesg = fsolve(f, guess, full_output = True)
        if verbose:
            print('{}: {}'.format(ind, mesg))
        if ier == 1:
            break
    return sol[0], sol[1]

def segment_intersect(line1, line2, verbose = False):
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
    (xa, ya), (xb, yb) = line1
    if isinstance(line2, Line):
        line2 = [(p.x, p.y) for p in line2.p]

    for i in range(len(line2) - 1):
        (xc, yc), (xd, yd) = line2[i], line2[i+1]

        M = np.array([[xb - xa, -xd + xc], [yb - ya, -yd + yc]])
        r = np.array([xc - xa, yc - ya])
        try:
            sol = np.linalg.solve(M, r)
        except np.linalg.LinAlgError:
            continue

        if (sol[0] <= 1) and (sol[1] <= 1) \
            and (sol[0] >= 0) and (sol[1] >= 0):
            return True, [(xc, yc), (xc + sol[1]*(xd-xc), yc + sol[1]*(yd - yc))]
    return False, [(np.nan, np.nan), (np.nan, np.nan)]


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

def UnfoldLabel(Dic:dict,Name:str)->str:
        '''
        Unfold Patch label (e.g. "ICT" -> "Inner Core Top")

        Parameters
        ----------
        Dic : dict
            Dictionnary containing description of acronym characters
        Name : str
            patch label

        Returns
        -------
        str
            Unfolded patch label.

        '''
        Output=[]
        for s in Name:
            if Dic.get(s)!=None:
                Output.append(Dic.get(s)+' ')
            else:
                Output.append(s)
        if len(Output)>0:
            return ''.join(Output)
        else:
            return ''
        
#def CheckRadialVertices(radial_vertices,i):
    # for j in range(len(radial_vertices[i]) - 1):
    #     NW = radial_vertices[i-1][j]
    #     NE = radial_vertices[i-1][j+1]
    #     SW = radial_vertices[i][j]
    #     SE = radial_vertices[i][j+1]
    #     AngleNWE=Getangle(Line([NW, NE]))
        
    #     cell_grid[i].append(Cell([, Line([SW, SE]), ), Line([SW, NW])]))
        