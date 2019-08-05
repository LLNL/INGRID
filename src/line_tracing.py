#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:07:45 2019

@author: watkins35
"""
# want to click on a point close to a zero, adjust to the exact point
# via newton's, and then trace out the contour line.
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from Root_Finder import RootFinder
from time import time
import geometry as geo
from geometry import Point, Line


class LineTracing:
    """ This class traces the polodal and radial lines of a given psi
    function based of the points where the user clicks.
    """
    def __init__(self, grid, params, eps=1e-3, tol=1e-3,
                 numPoints=25, dt=0.02, option='xpt_circ', direction='cw'):
        """ Supports many different types of line drawing,
        parameters
        ----------
        grid - the grid object upon which the lines will be drawn.
                expects one from the Efit_Data class.
        eps - epsilon :: size of the circle drawn around the zero point
        tol - tolerance :: how close to the final point the line must get
                before converging. also a circle
        numPoints - number of points in the circle of radius eps
        option - accepts a string. default is contour
               - 'theta' : draws the contour line where the user clicks
               - 'rho' : draws the line orthogonal to the
                          contour line where the user clicked
               - 'xpt_circ': uses the root finder to find the root
                           closest to where the user clicked. Then
                           finds the points around that circle a distance
                           epsilon away
        """
        self.grid = grid

        self.cid = self.grid.ax.figure.canvas.mpl_connect('button_press_event',
                                                          self)
        self.eps = eps
        self.tol = tol
        self.numPoints = numPoints
        self.dt = dt
        self.option = option
        if self.option == 'xpt_circ':
            self.root = RootFinder(self.grid)
            print('Entering Fuzzy click mode. Click near a zero point.')

        zdim = grid.zmax-grid.zmin
        rdim = grid.rmax-grid.rmin
        self.max_step = params['step_ratio'] * max(rdim, zdim)
        self.first_step = 1e-3  # self.max_step

        self.set_function(option, direction)  # initialize the function

    def differential_theta(self, t, xy):
        """ coupled set of differential equations
        to trace the poloidal lines """
        R, Z = xy
        B_R = (1/R)*self.grid.get_psi(R, Z, tag='vz')
        B_Z = -(1/R)*self.grid.get_psi(R, Z, tag='vr')
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_R/B
        dZ = B_Z/B
        if self.dir == 'cw':
            return np.array([dR, dZ])
        else:
            return -np.array([dR, dZ])

    def differential_rho(self, t, xy):
        """ coupled set of differential equations
        to trace the radial lines """
        R, Z = xy
        B_R = (1/R)*self.grid.get_psi(R, Z, tag='vz')
        B_Z = -(1/R)*self.grid.get_psi(R, Z, tag='vr')
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_Z/B
        dZ = -B_R/B
        if self.dir == 'cw':
            return np.array([dR, dZ])
        else:
            return -np.array([dR, dZ])

    def set_function(self, option, direction):
        self.option = option
        self.dir = direction
        if self.option == 'theta':
            self.function = self.differential_theta
        elif self.option == 'rho':
            self.function = self.differential_rho

    def __call__(self, event):
        """ Activates upon mouse click. Will call the appropriate function
        based on parameters
        """
        if event.button == 3:
            print('...disabling fuzzy click mode...')
            event.canvas.mpl_disconnect(self.cid)
            return
        if event.inaxes != self.grid.ax.axes:
            # safe guards against clicking outside the figure
            return

        x0, y0 = event.xdata, event.ydata
        plt.draw()

        if self.option in ['theta', 'rho']:
            self.draw_line((x0, y0))

        elif self.option == 'xpt_circ':
            # no longer using this... 
            self.root.find_root(x0, y0)
            r, z = self.root.final_root
            self.calc_equal_psi_points(r, z)

            # this part is just plotting
            for key, (x, y) in self.eq_psi.items():
                plt.plot(x, y, 'x', label=key)

            plt.legend()
            plt.draw()
            
    def disconnect(self):
        """ turns of the click functionality """
        self.grid.ax.figure.canvas.mpl_disconnect(self.cid)
        self.root.disconnect()
    
    def calc_equal_psi_points(self, r, z, theta2d=False, err_circles=False,
                              show_eq_psi_points=False):
        """ draws a circle around the xpt, and saves four coordinates for
        each direction the poloidal line will travel.
        """
        def get_circle(eps, x0, y0, psi=False):
            """ traces a circle of given radius (eps) around a point.
            Can be set to return the value of psi at each point
            """
            x, y, z = [], [], []
            for th in self.theta:
                x.append(x0 + eps * np.cos(th))
                y.append(y0 + eps * np.sin(th))
                if psi:
                    z.append(self.grid.get_psi(x[-1], y[-1]))
            if psi:
                return x, y, z
            return x, y

        circ = {}
        zero = {'x': r, 'y': z, 'psi': self.grid.get_psi(r, z)}
        zero_psi_line = np.full(self.numPoints, zero['psi'])
        self.theta = np.linspace(0, 2*np.pi, self.numPoints)

        # get the circle
        circ['x'], circ['y'], circ['psi'] = get_circle(self.eps, r, z,
                                                       psi=True)

        if err_circles:
            # include some other error lines around the zero point
            cx0, cy0 = get_circle(1e-2, r, z)
            self.grid.ax.plot(cx0, cy0, label='1e-2', color='black')
            plt.draw()
            cx1, cy1 = get_circle(1e-3, r, z)
            self.grid.ax.plot(cx1, cy1, label='1e-3', color='green')
            plt.draw()
            cx2, cy2 = get_circle(1e-4, r, z)
            self.grid.ax.plot(cx2, cy2, label='1e-4', color='blue')
            plt.draw()
            plt.legend()

        if theta2d:
            # 2d representation of the circle of psi values
            plt.figure('2d psi circle')
            plt.plot(self.theta, circ['psi'], '.')
            plt.plot(self.theta, zero_psi_line)
            plt.show()

        adjusted_psi = circ['psi'] - zero_psi_line
        sign_changes = []
        zeros = []
        for i in range(len(adjusted_psi)-1):
            if adjusted_psi[i] > 0 and adjusted_psi[i+1] < 0:
                sign_changes.append((i, i+1))
            elif adjusted_psi[i] < 0 and adjusted_psi[i+1] > 0:
                sign_changes.append((i, i+1))
            elif adjusted_psi[i] == 0:
                zeros.append(i)

        root = []

        for i, j in sign_changes:
            # i and j are the idices of our points
            x0, y0 = self.theta[i], adjusted_psi[i]
            x1, y1 = self.theta[j], adjusted_psi[j]
            sol = root_scalar(lambda x: (y1-y0)/(x1-x0) * (x-x0) + y0,
                              bracket=[x0, x1])
            root.append(sol.root)

        equal_psi_coords = []
        for th in root:
            # this is the value for theta
            xr = zero['x'] + self.eps * np.cos(th)
            yr = zero['y'] + self.eps * np.sin(th)
            equal_psi_coords.append((xr, yr))
        for th in zeros:
            # this is the index of the value
            xr = zero['x'] + self.eps * np.cos(self.theta[th])
            yr = zero['y'] + self.eps * np.sin(self.theta[th])
            equal_psi_coords.append((xr, yr))

        self.eq_psi = {'NE': (equal_psi_coords[0]),
                       'NW': (equal_psi_coords[1]),
                       'SW': (equal_psi_coords[2]),
                       'SE': (equal_psi_coords[3])}

        # calculate N, S, E, W
        ne = geo.Vector(self.eq_psi['NE'], (zero['x'], zero['y']))
        nw = geo.Vector(self.eq_psi['NW'], (zero['x'], zero['y']))
        sw = geo.Vector(self.eq_psi['SW'], (zero['x'], zero['y']))
        se = geo.Vector(self.eq_psi['SE'], (zero['x'], zero['y']))

        self.eq_psi['N'] = geo.calc_mid_point(ne, nw)  # NORTH
        self.eq_psi['W'] = geo.calc_mid_point(nw, sw)  # WEST
        self.eq_psi['S'] = geo.calc_mid_point(sw, se)  # SOUTH
        self.eq_psi['E'] = geo.calc_mid_point(se, ne)  # EAST
        
        if show_eq_psi_points:
            # let's us see where the offset points we will trace from are
            for key, (x, y) in self.eq_psi.items():
                plt.plot(x, y, 'x', label=key)
            plt.legend()
            plt.draw()
        

    def draw_line(self, rz_start, rz_end=None, color='green',
                  option=None, direction=None, show_plot=False, text=False):
        """ rz_start and rz_end are defaulted to be the same point.
        Checks the new set of points calculated and if it is near enough to
        the end point, it stops calculating.

        direction :: determines if the function plots clockwise (cw) or
                          counterclockwise (ccw). default is None.
        """

        if option is not None and direction is not None:
            self.set_function(option, direction)

        if rz_end is None:
            print('testing for loop completion')
            test = 'point'
            rz_end = rz_start
            xf, yf = rz_end

        elif isinstance(rz_end, Point):
            print('Testing for point convergence')
            test = 'point'
            xf = rz_end.x
            yf = rz_end.y
            
        elif isinstance(rz_end, Line):
            print('Testing for line convergence')
            test = 'line'
       
        
        elif np.shape(rz_end) == ():
            print('Testing for psi convergence')
            test = 'psi'
            
            
            
            
            

        elif np.shape(rz_end) == (2,):
#            print('Testing for point convergence')
#            test = 'point'
#            xf, yf = rz_end
            print("Point convergence shape error")
            print(rz_end)

        elif np.shape(rz_end) == (2, 2):
            print("Line convergence shape error")
            print(rz_end)


        else:
            test = 'multiple lines'
            print('this could be interesting')

        # reset each iteration
        if isinstance(rz_start, geo.Point):
            ynot = (rz_start.x, rz_start.y)
        else:
            ynot = rz_start
            
        # check to make sure ynot is not an object            
        if isinstance(ynot, Point):
            print('\nynot is a point object')
            print('Object\n')
            
            
        # size for each line segment
        told, tnew = 0, self.dt

        count = 0
        # unpack boundaries
        rmin = self.grid.rmin
        rmax = self.grid.rmax
        zmin = self.grid.zmin
        zmax = self.grid.zmax

        self.time_in_converged = 0
        line = [Point(ynot)]

        def converged(points):
            """ checks for converence of the line in various ways """

            def success(message):
                """ Displays a message so we know the convergence worked. """
                print('Converged via {}.'.format(message))
                print('Iterations: ', count)
                print('Spent {} '.format(self.time_in_converged)
                      + 'seconds checking convergence.')

            def save_line(x, y):
                """ Plots the current line segments and saves
                it for future use
                x: list -- r endpoints
                y: list -- z endpoints
                """
                line.append(Point(x[-1], y[-1]))
                if show_plot:
                    self.grid.ax.plot(x, y, '.-', linewidth='2', color=color)
                    plt.draw()
                    plt.pause(1e-15)

            t1 = time()
            # don't go off the plot
            tol = 1e-3
            if (any(abs(rmin - points[0]) < tol)
                    or any(abs(rmax - points[0]) < tol)
                    or any(abs(zmin - points[1]) < tol)
                    or any(abs(zmax - points[1]) < tol)):
                # this is just here as a safegaurd, none of the lines
                # we care about should go off the grid
                if text: success('edge')
                return True

            # check if any point is close enough to the endpoint
            if test == 'point':
                if (any(abs(points[0]-xf) < self.tol)
                        and any(abs(points[1]-yf) < self.tol)
                        and count > 5):
                    if text: success('endpoint')
                    new_x, new_y = geo.truncate_list(self.x, self.y, xf, yf)
                    save_line([new_x[0], new_x[-1]], [new_y[0], new_y[-1]])
                    return True

            elif test == 'line':
                # endpoints define the latest line segment
                p1 = (points[0][0], points[1][0])
                p2 = (points[0][-1], points[1][-1])
                p1s, p2s = geo.test2points(p1, p2, rz_end)
                if p1s != p2s:
                    if text: success('line crossing')
                    # pass in the line instance
#                    segment = Line([Point(p1), Point(p2)])
                    endLine = rz_end.points()
                    r, z = geo.intersect((p1, p2), endLine)
                    save_line([p1[0], r], [p1[1], z])
                    return True

            elif test == 'psi':
                x1, y1 = points[0][0], points[1][0]
                x2, y2 = points[0][-1], points[1][-1]

                psi1 = self.grid.get_psi(x1, y1)
                psi2 = self.grid.get_psi(x2, y2)

                if (psi1 - rz_end)*(psi2 - rz_end) < 0:
                    if text: success('psi test')
                    # need to find coords for the value of psi that we want
                    
                    def f(x):
                        # must manually calculate y each time we stick to
                        # the line of interest
                        y = (y2-y1)/(x2-x1)*(x-x1)+y1
                        return rz_end - self.grid.get_psi(x, y)

                    sol = root_scalar(f, bracket=[x1, x2])
                    r_psi = sol.root
                    z_psi = (y2-y1)/(x2-x1)*(r_psi-x1)+y1

                    save_line([x1, r_psi], [y1, z_psi])
                    return True

            else:
                print('Error: No termination criteria specified.')

            # if convergence didn't occur
            if count > 0:
                # plot the line like normal
                save_line([self.x[0], self.x[-1]], [self.y[0], self.y[-1]])

            t2 = time()
            self.time_in_converged += t2 - t1

        # keep track of the line segment generated
        points = np.zeros([2, 2])  # initial length is arbitrary
        start = time()
        while not converged(points):
            t_span = (told, tnew)
            # solve the system of differential equations
            sol = solve_ivp(self.function, t_span, ynot, method='LSODA',
                            first_step=self.first_step, max_step=self.max_step,
                            rtol=1e-12, atol=1e-11)
            # unpack
            self.x = sol.y[0]
            self.y = sol.y[1]

            ynot = [self.x[-1], self.y[-1]]


            told, tnew = tnew, tnew + self.dt
            points = (self.x, self.y)

            if count > 350:
                print('did not converge, exiting...')
                print('Iterations: ', count)
                end = time()
                print('Took {} '.format(end-start)
                      + 'seconds trying to converge.')
                break
            count += 1
        end = time()
#        print('Drew for {} seconds\n'.format(end-start))


        return Line(line)


if __name__ == '__main__':
    from Read_Psi_Data import read_psi_data

    plt.close('all')
    grid = read_psi_data()
    grid.Calculate_PDeriv()
    grid.plot_data()
    click = LineTracing(grid, option='rho')
