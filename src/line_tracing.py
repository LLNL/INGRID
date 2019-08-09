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


class LineTracing:
    """
    This class traces the polodal and radial lines of a given psi
    function based of the points where the user clicks.
    
    Parameters
    ----------
    grid : Setup_Grid_Data.Efit_Data
        The grid object upon which the lines will be drawn.
    params : dict
        Dictionary of the grid parameters. Should have been read in
        from a namelist file.
    eps : float, optional
        Short for epsilon. Specifies the size of the circle drawn
        around the zero point.
    tol : float, optional
        Short for tolerance. Specifies how close to the final point the
        line must get before converging. Also defines a circle.
    numPoints : int
        Number of points in the circle of radius eps.
    dt : float, optional
        Specify the size of each line segment that is traced by
        scipy.integrate.solve_ivp.
    option : str, optional
        'theta' draws the poloidal line where the user clicks.
        'rho' draws the radial line where the user clicked.
        'xpt_circ': uses the root finder to find the root closest
        to where the user clicked. Then finds the points around
        that circle a distance epsilon away.
    direction : str, optional
        'cw' or 'ccw'. Specifies clockwise or counterclockwise line
        tracing.
    """

    def __init__(self, grid, params, eps=1e-3, tol=1e-3,
                 numPoints=25, dt=0.02, option='xpt_circ', direction='cw'):
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
        # TODO: set self.first_step as a parameter
        self.first_step = 1e-3  # self.max_step

        # initialize the function
        self._set_function(option, direction)

    def _differential_theta(self, t, xy):
        """
        Coupled set of differential equations
        to trace the poloidal lines.
        """

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

    def _differential_rho(self, t, xy):
        """
        Coupled set of differential equations
        to trace the radial lines.
        """
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

    def _set_function(self, option, direction):
        self.option = option
        self.dir = direction
        if self.option == 'theta':
            self.function = self._differential_theta
        elif self.option == 'rho':
            self.function = self._differential_rho

    def __call__(self, event):
        """
        Activates upon mouse click. Will call the appropriate function
        based on parameters.
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
            self.draw_line((x0, y0), show_plot=True, text=True)

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
        """ Turns off the click functionality """
        self.grid.ax.figure.canvas.mpl_disconnect(self.cid)
        self.root.disconnect()

    def calc_equal_psi_points(self, r, z, theta2d=False, err_circles=False,
                              show_eq_psi_points=False):
        """ Draws a circle around the xpt, and saves four coordinates
        for each direction the poloidal line will travel. Also
        calculates the location of the points that bisect the poloidal
        lines.

        Parameters
        ---------
        r : float
            R coordinate of the point in question.
        z : float
            Z coordinate of the point.
        theta2d : bool, optional
            Displays a plot on a seperate figure of the value of psi on
            circle that is traced versus the value of psi at (r, z).
        err_circles : bool, optional
            Draw some circles around (r, z) to show what the allowed
            error is for the point convergence criteria.
        show_eq_psi_points : bool, optional
            Plots an x on the NE, NW, SE, SW points that are found.
        """

        def get_circle(eps, x0, y0, psi=False):
            """ Traces a circle of given radius (eps) around a point.
            Can be set to return the value of psi at each point.
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
        circ['x'], circ['y'], circ['psi'] = get_circle(self.eps, r, z, psi=True)

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
        """
        Uses scipy.integrate.solve_ivp to trace poloidal or radial
        lines. Uses the LSODA method to solve the differential
        equations. Three options for termination criteria, specified
        by rz_end.

        Parameters
        ----------
        rz_start : array-like or geometry.Point
            Starting location for line tracing.
        rz_end : dict, optional
            Defaults to be rz_start. This is how we specify the
            termination critera. i.e. {'point': Point}, {'line': Line},
            {'psi': Psi}
            Points can be a geometry.Point, or array-like
            i.e. (x, y)
            Lines can be a geometry.Line, or array-like
            i.e. ((x, y), (x, y))
            Psi must be a scalar, i.e. 1.1, and specifies the
            level of psi to stop on.
        color : str, optional
            Specifies the color of the produced grid lines.
        option : str, optional
            Change which differential equation is used in the line
            tracing proccess. 'theta', 'rho'
        direction : str
            determines if the function plots clockwise (cw) or
            counterclockwise (ccw). default is None.
        show_plot : bool, optional
            Show the user real-time tracing and the line tracer works.
        text : bool, optional
            Prints convergence method, number of iterations, and
            time taken to the terminal window.

        Returns
        -------
        line : geometry.Line
            Curved line consisting of the start and end points of each
            segment calculated by solve_ivp. Does not store the
            intermediate points.
        """

        if option is not None and direction is not None:
            self._set_function(option, direction)

        # check rz_start
        if isinstance(rz_start, geo.Point):
            ynot = (rz_start.x, rz_start.y)
        else:
            ynot = rz_start

        # check rz_end:
        if rz_end is None:
            print('testing for loop completion')
            test = 'point'
            rz_end = ynot
            xf, yf = rz_end

        elif rz_end.keys() == ['point']:
            print('Testing for point convergence')
            test = 'point'
            if isinstance(rz_end['point'], geo.Point):
                xf = rz_end['point'].x
                yf = rz_end['point'].y
            else:
                xf = rz_end['point'][0]
                yf = rz_end['point'][1]

        # TODO: generalize the line termination to include "curves", or lines
        # defined by more than one set of points.
        elif rz_end.keys() == ['line']:
            print('Testing for line convergence')
            test = 'line'
            if isinstance(rz_end['line'], geo.Line):
                # extract
                endLine = rz_end['line'].points()
            else:
                # form ((),())
                endLine = rz_end['line']

        elif rz_end.keys() == ['psi']:
            print('Testing for psi convergence')
            test = 'psi'
            psi_test = rz_end['psi']

        else:
            print('rz_end type not recognized')

        # size for each line segment
        told, tnew = 0, self.dt

        count = 0
        # unpack boundaries
        rmin = self.grid.rmin
        rmax = self.grid.rmax
        zmin = self.grid.zmin
        zmax = self.grid.zmax

        self.time_in_converged = 0
        line = [geo.Point(ynot)]

        def converged(points):
            # checks for converence of the line in various ways

            def success(message):
                # Displays a message so we know the convergence worked.
                print('Converged via {}.'.format(message))
                print('Iterations: ', count)
                print('Spent {} '.format(self.time_in_converged)
                      + 'seconds checking convergence.')

            def save_line(x, y):
                # Plots the current line segments and saves
                # it for future use
                # x: list -- r endpoints
                # y: list -- z endpoints
                line.append(geo.Point(x[-1], y[-1]))
                if show_plot:
                    self.grid.ax.plot(x, y, '.-', linewidth='2', color=color)
                    plt.draw()
                    plt.pause(1e-15)

            t1 = time()
            # don't go off the plot
            boundary = [((rmin, zmin), (rmin, zmax)),
                        ((rmin, zmax), (rmax, zmax)),
                        ((rmax, zmax), (rmax, zmin)),
                        ((rmax, zmin), (rmin, zmin))]

            # check for intersections
            for edge in boundary:
                p1 = (points[0][0], points[1][0])
                p2 = (points[0][-1], points[1][-1])
                p1s, p2s = geo.test2points(p1, p2, edge)
                if p1s != p2s:
                    print('The line went over the boundary')
                    if text:
                        success('edge')
                    r, z = geo.intersect((p1, p2), edge)
                    save_line([p1[0], r], [p1[1], z])
                    return True

            # check if any point is close enough to the endpoint
            # this currently works by checking all the intermediate points
            # with the termination point.
            # TODO: develope a method to check if the point is close
            # to a line defined by two points.
            if test == 'point':
                if (any(abs(points[0]-xf) < self.tol)
                        and any(abs(points[1]-yf) < self.tol)
                        and count > 5):
                    if text:
                        success('endpoint')
                    save_line([points[0][0], xf], [points[1][0], yf])
                    return True

            elif test == 'line':
                # endpoints define the latest line segment
                # TODO change the point criteria
                p1 = (points[0][0], points[1][0])
                p2 = (points[0][-1], points[1][-1])
                p1s, p2s = geo.test2points(p1, p2, endLine)
                if p1s != p2s:
                    if text:
                        success('line crossing')
                    r, z = geo.intersect((p1, p2), endLine)
                    save_line([p1[0], r], [p1[1], z])
                    return True

            elif test == 'psi':
                x1, y1 = points[0][0], points[1][0]
                x2, y2 = points[0][-1], points[1][-1]

                psi1 = self.grid.get_psi(x1, y1)
                psi2 = self.grid.get_psi(x2, y2)

                if (psi1 - psi_test)*(psi2 - psi_test) < 0:
                    if text: success('psi test')
                    # need to find coords for the value of psi that we want

                    def f(x):
                        # must manually calculate y each time we stick it into
                        # the line of interest
                        y = (y2-y1)/(x2-x1)*(x-x1)+y1
                        return psi_test - self.grid.get_psi(x, y)

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
                save_line([points[0][0], points[0][-1]], [points[1][0], points[1][-1]])

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

            # TODO: reduce the list passed in to only be the endpoints
            # of the line, and not include all the intermediate points
            # currently there is an error in the calculation of midpoints
            # for the sepratrix, when we try to truncate the points list.
            # Occurs with point converence.
            points = (self.x, self.y)

            if count > 350:
                print('did not converge, exiting...')
                print('Iterations: ', count)
                end = time()
                print('Took {} '.format(end-start)
                      + 'seconds trying Lto converge.')
                break
            count += 1
        end = time()
        if text:
            print('Drew for {} seconds\n'.format(end-start))

        return geo.Line(line)


if __name__ == '__main__':
    from Read_Psi_Data import read_psi_data

    plt.close('all')
    grid = read_psi_data()
    grid.Calculate_PDeriv()
    grid.plot_data()
    click = LineTracing(grid, option='rho')
