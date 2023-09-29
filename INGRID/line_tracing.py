
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, minimize
from time import time
from INGRID.geometry import Point, Line, rotate, find_split_index, \
    segment_intersect, test2points
from INGRID.exceptions import RegionEntered

class LineTracing:
    """
    This class traces the polodal and radial lines of a given psi
    function based of the points where the user clicks.

    Parameters
    ----------
    grid : EfitData.EfitData
        The grid object upon which the lines will be drawn.
    settings : dict
        YAML file containing all INGRID parameters.
    eps : float, optional
        Short for epsilon. Specifies the size of the circle drawn
        around the zero point.
    tol : float, optional
        Short for tolerance. Specifies how close to the final point the
        line must get before converging. Also defines a circle.
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

    def __init__(self, grid, settings, eps=1e-6, tol=5e-5, first_step=1e-5, dt=0.01, option='xpt_circ', direction='cw'):

        self.grid = grid
        self.settings = settings
        self.option = option
        self.NSEW_lookup = {
                            'xpt1': {'coor': {}},
                            'xpt2': {'coor': {}, 'theta': {} }
                            }

        try:
            settings['integrator_settings']['eps']
        except KeyError:
            settings['integrator_settings'].update({'eps': eps})

        self.eps = settings['integrator_settings']['eps']
        print('eps set to {}'.format(self.eps))

        try:
            settings['integrator_settings']['tol']
        except KeyError:
            settings['integrator_settings'].update({'tol': tol})

        self.tol = settings['integrator_settings']['tol']
        print('tol set to {}'.format(self.tol))

        try:
            settings['integrator_settings']['dt']
        except KeyError:
            settings['integrator_settings'].update({'dt': dt})

        self.dt = settings['integrator_settings']['dt']
        print('dt set to {}'.format(self.dt))

        try:
            settings['integrator_settings']['first_step']
        except KeyError:
            settings['integrator_settings'].update({'first_step': first_step})

        self.first_step = settings['integrator_settings']['first_step']
        print('first_step set to {}'.format(self.first_step))

        try:
            settings['integrator_settings']['step_ratio']
        except KeyError:
            settings['integrator_settings'].update({'step_ratio': step_ratio})

        self.step_ratio = settings['integrator_settings']['step_ratio']
        print('step_ratio set to {}'.format(self.step_ratio))

        zdim = grid.zmax - grid.zmin
        rdim = grid.rmax - grid.rmin

        try:
            settings['integrator_settings']['max_step']
        except KeyError:
            settings['integrator_settings'].update({'max_step': self.step_ratio * max(rdim, zdim)})
        self.max_step = settings['integrator_settings']['max_step']
        print('max_step set to {}'.format(self.max_step))

        # initialize the function
        self._set_function(option, direction)

    def _differential_theta(self, t, xy):
        """
        Coupled set of differential equations
        to trace the poloidal lines.
        """

        R, Z = xy
        B_R = (1 / R) * self.grid.get_psi(R, Z, tag='vz')
        B_Z = -(1 / R) * self.grid.get_psi(R, Z, tag='vr')
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_R / B
        dZ = B_Z / B
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
        B_R = (1 / R) * self.grid.get_psi(R, Z, tag='vz')
        B_Z = -(1 / R) * self.grid.get_psi(R, Z, tag='vr')
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_Z / B
        dZ = -B_R / B
        if self.dir == 'cw':
            return np.array([dR, dZ])
        else:
            return -np.array([dR, dZ])

    def _differential_r_const(self, t, xy):
        """
        Coupled set of differential equations
        to trace vertical lines.
        """
        R, Z = xy
        B_R = 0
        B_Z = 1
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_R / B
        dZ = B_Z / B
        if self.dir == 'cw':
            return np.array([dR, dZ])
        else:
            return -np.array([dR, dZ])

    def _differential_z_const(self, t, xy):
        """
        Coupled set of differential equations
        to trace horizontal lines.
        """
        R, Z = xy
        B_R = np.cos(self.tilt_angle)
        B_Z = np.sin(self.tilt_angle)
        B = np.sqrt(B_R**2 + B_Z**2)
        dR = B_R / B
        dZ = B_Z / B
        if self.dir == 'cw':
            return np.array([dR, dZ])
        else:
            return -np.array([dR, dZ])

    def _set_function(self, option, direction):
        self.option = option
        self.dir = direction
        if self.option == 'theta':
            self.function = self._differential_theta
        elif self.option == 'r_const':
            self.function = self._differential_r_const
        elif self.option == 'rho':
            self.function = self._differential_rho
        elif self.option == 'z_const':
            self.function = self._differential_z_const

    def analyze_saddle(self, xpt, xpt_ID):
        """
        Finds theta values to be tested for N and S directions
        """

        rxpt, zxpt = xpt
        hessian = np.zeros((2, 2))
        hessian[0, 0] = self.grid.get_psi(rxpt, zxpt, tag='vrr')
        hessian[1, 1] = self.grid.get_psi(rxpt, zxpt, tag='vzz')
        hessian[0, 1] = self.grid.get_psi(rxpt, zxpt, tag='vrz')
        hessian[1, 0] = self.grid.get_psi(rxpt, zxpt, tag='vrz')

        eigval, eigvect = np.linalg.eig(hessian)

        index = 1
        if np.sign(eigval[0]) == -1: 
            index = 0

        origin = np.array([rxpt, zxpt])

        N = np.array([rxpt + self.first_step * eigvect[0][index], zxpt + self.first_step * eigvect[1][index]])
        W = rotate(N, np.pi / 2, origin)
        S = np.array([rxpt - self.first_step * eigvect[0][index], zxpt - self.first_step * eigvect[1][index]])
        E = rotate(S, np.pi / 2, origin)

        NW = rotate(N, np.pi / 4, origin)
        SW = rotate(W, np.pi / 4, origin)
        SE = rotate(S, np.pi / 4, origin)
        NE = rotate(E, np.pi / 4, origin)

        self.NSEW_lookup[xpt_ID]['coor'] = {'center': xpt, 'N': N, 'S': S, 'E': E, 'W': W,
                                            'NE': NE, 'NW': NW, 'SE': SE, 'SW': SW}

    def rotate_NSEW_lookup(self, xpt_ID, turns=2):
        swap_key = {'N': 'NW', 'S': 'SE',
                    'W': 'SW', 'E': 'NE',
                    'NW': 'W', 'SE': 'E',
                    'NE': 'N', 'SW': 'S',
                    'center': 'center'}
        for i in range(turns):
            temp_dict = {}
            for k, v in self.NSEW_lookup[xpt_ID]['coor'].items():
                temp_dict[swap_key[k]] = v
            self.NSEW_lookup[xpt_ID]['coor'] = temp_dict

    def flip_NSEW_lookup(self, xpt_ID):
        swap_key = {'N': 'S', 'S': 'N',
                    'W': 'E', 'E': 'W',
                    'NW': 'SE', 'SE': 'NW',
                    'NE': 'SW', 'SW': 'NE',
                    'center': 'center'}
        temp_dict = {}
        for k, v in self.NSEW_lookup[xpt_ID]['coor'].items():
            temp_dict[swap_key[k]] = v
        self.NSEW_lookup[xpt_ID]['coor'] = temp_dict

    def map_xpt(self, xpt, magx, xpt_ID='xpt1', visual=False, verbose=False):

        if xpt_ID == 'xpt1':
            self.analyze_saddle(xpt, xpt_ID)
            self._set_function('rho', 'cw')

            NS_buffer = self.NSEW_lookup[xpt_ID]['coor']

            N_sol = solve_ivp(self.function, (0, self.dt), NS_buffer['N'], method='LSODA',
                first_step=self.first_step, max_step=self.max_step, rtol=1e-13, atol=1e-12).y
            S_sol = solve_ivp(self.function, (0, self.dt), NS_buffer['S'], method='LSODA',
                first_step=self.first_step, max_step=self.max_step, rtol=1e-13, atol=1e-12).y

            N_guess = (N_sol[0][-1], N_sol[1][-1])
            S_guess = (S_sol[0][-1], S_sol[1][-1])

            r_bounds = (self.grid.rmin, self.grid.rmax)
            z_bounds = (self.grid.zmin, self.grid.zmax)
            N_minimizer = minimize(self.PsiCostFunc, N_guess, method='L-BFGS-B',
                jac=self.grid.Gradient, bounds=[r_bounds, z_bounds]).x
            S_minimizer = minimize(self.PsiCostFunc, S_guess, method='L-BFGS-B',
                jac=self.grid.Gradient, bounds=[r_bounds, z_bounds]).x

            if (np.linalg.norm(N_minimizer - magx) >= self.eps):
                self.flip_NSEW_lookup(xpt_ID)

            self.config = 'LSN' if self.NSEW_lookup['xpt1']['coor']['N'][1] > xpt[1] else 'USN'

        elif xpt_ID == 'xpt2':

            from matplotlib.patches import Polygon

            def cb_region_check(xk):
                if core_polygon.get_path().contains_point(xk):
                    raise RegionEntered(message='# Entered Core...', region='Core')
                elif pf_polygon.get_path().contains_point(xk):
                    raise RegionEntered(message='# Entered PF...', region='PF')

            def limiter_split(start, end, limiter):
                start_index, sls = find_split_index(start, limiter)
                end_index, sls = find_split_index(end, limiter)
                if end_index <= start_index:
                    limiter.p = limiter.p[start_index:] + limiter.p[:start_index + 1]
                return limiter

            try:
                visual = self.grid.parent.settings['DEBUG']['visual']['SF_analysis']
            except KeyError:
                visual = False
            try:
                verbose = self.grid.parent.settings['DEBUG']['verbose']['SF_analysis']
            except KeyError:
                verbose = False


            limiter = self.grid.parent.LimiterData

            xpt1 = self.NSEW_lookup['xpt1']['coor']
            magx = np.array([self.grid.parent.settings['grid_settings']['rmagx'], self.grid.parent.settings['grid_settings']['zmagx']])

            # Create mid-line
            LHS_Point = Point([magx[0] - 1e6, magx[1]])
            RHS_Point = Point([magx[0] + 1e6, magx[1]])
            midline = Line([LHS_Point, RHS_Point])

            # Generate Vertical Mid-Plane line
            Lower_Point = Point([magx[0], magx[1] - 1e6])
            Upper_Point = Point([magx[0], magx[1] + 1e6])
            topline = Line([Lower_Point, Upper_Point])

            # Drawing Separatrix
            print('# Tracing core region...')
            xptNW_midLine = self.draw_line(xpt1['NW'], {'line': midline}, option='theta', direction='cw', show_plot=visual, text=verbose)
            xptNE_midLine = self.draw_line(xpt1['NE'], {'line': midline}, option='theta', direction='ccw', show_plot=visual, text=verbose)
            midline_topline_west = self.draw_line(xptNW_midLine[-1], {'line': topline}, option='theta', direction='cw', show_plot=visual, text=verbose)
            midline_topline_east = self.draw_line(xptNE_midLine[-1], {'line': topline}, option='theta', direction='ccw', show_plot=visual, text=verbose)

            print('# Tracing PF region...')
            xptSW_limiter = self.draw_line(xpt1['SW'], {'line': limiter}, option='theta', direction='ccw', show_plot=visual, text=verbose)[::-1]
            xptSE_limiter = self.draw_line(xpt1['SE'], {'line': limiter}, option='theta', direction='cw', show_plot=visual, text=verbose)
            core_boundary = xptNW_midLine + midline_topline_west + midline_topline_east[::-1] + xptNE_midLine[::-1]
            core_polygon = Polygon(core_boundary, fill=True, closed=True, color='violet', label='Core')

            pf_boundary = (
                xptSW_limiter + 
                xptSE_limiter + (
                    limiter_split(
                        start=xptSE_limiter[-1], 
                        end=xptSW_limiter[0], 
                        limiter=limiter
                    ).split(xptSE_limiter[-1])[1]).split(xptSW_limiter[0], add_split_point=True)[0]
                )
            pf_polygon = Polygon(pf_boundary, fill=True, closed=True, color='dodgerblue', label='PF_1')

            self.RegionPolygon = {'Core': core_polygon, 'PF_1': pf_polygon}

            if (pf_polygon.get_path().contains_point(xpt)):
                print("# Snowflake-plus...")
                self.analyze_saddle(xpt, xpt_ID='xpt2')
                # Rotate NSEW so that 'W' is set to 'N' and 'E' to 'S'
                self.rotate_NSEW_lookup(xpt_ID=xpt_ID, turns=6)

                self._set_function('rho', 'ccw')
                xpt2 = self.NSEW_lookup['xpt2']['coor']
                N_sol = solve_ivp(self.function, (0, self.dt), xpt2['N'], method='LSODA',
                        first_step=self.first_step, max_step=self.max_step, rtol=1e-13, atol=1e-12).y
                N_guess = (N_sol[0][-1], N_sol[1][-1])

                self.RegionLineCut = self.draw_line(xpt2['N'], {'line_group': [xptSW_limiter, xptSE_limiter, limiter]},
                    option='rho', direction='ccw', show_plot=visual, text=verbose)

                if np.allclose(self.line_group_intersect - limiter, 0.0):
                    self.flip_NSEW_lookup(xpt_ID)
                    xpt2 = self.NSEW_lookup['xpt2']['coor']
                    self.RegionLineCut = self.draw_line(xpt2['N'], {'line_group': [xptSW_limiter, xptSE_limiter, limiter]},
                        option='rho', direction='ccw', show_plot=visual, text=verbose)

                if np.allclose(self.line_group_intersect - xptSE_limiter, 0.0):
                    self.config = 'SF75'
                elif np.allclose(self.line_group_intersect - xptSW_limiter, 0.0):
                    self.config = 'SF105'

            else:
                print("# Snowflake-minus...")

                # Obtain candidate NSEW directions
                self.analyze_saddle(xpt, xpt_ID='xpt2')

                # Prepare minimizer for identification of SF- type
                self._set_function('rho', 'cw')
                xpt2 = self.NSEW_lookup['xpt2']['coor']
                N_sol = solve_ivp(self.function, (0, self.dt), xpt2['N'], method='LSODA',
                        first_step=self.first_step, max_step=self.max_step, rtol=1e-13, atol=1e-12).y
                S_sol = solve_ivp(self.function, (0, self.dt), xpt2['S'], method='LSODA',
                        first_step=self.first_step, max_step=self.max_step, rtol=1e-13, atol=1e-12).y

                N_guess = (N_sol[0][-1], N_sol[1][-1])
                S_guess = (S_sol[0][-1], S_sol[1][-1])

                for guess in [S_guess, N_guess]:
                    try:
                        minimize(self.PsiCostFunc, guess, method='trust-ncg', jac=self.grid.parent.PsiNorm.Gradient, hess=self.grid.parent.PsiNorm.Hessian,
                            options={'initial_trust_radius': self.eps, 'max_trust_radius': self.dt}, callback=cb_region_check)
                    except RegionEntered as e:
                        region = e.region

                        if region == 'Core':
                            # True south should land in region of interest
                            if guess is S_guess:
                                self.flip_NSEW_lookup(xpt_ID)
                                xpt2 = self.NSEW_lookup['xpt2']['coor']

                            # Determine whether SF15 or SF165 based off of intersection with core boundary
                            self.RegionLineCut = self.draw_line(xpt2['N'], {'line_group': [xptNE_midLine, xptNW_midLine,
                                midline_topline_west, midline_topline_east, limiter]},
                                option='rho', direction='cw', show_plot=visual, text=verbose)

                            if np.allclose(self.line_group_intersect - xptNE_midLine, 0.0):
                                self.config = 'SF15'
                            elif np.allclose(self.line_group_intersect - xptNW_midLine, 0.0):
                                self.config = 'SF165'
                            elif np.allclose(self.line_group_intersect - midline_topline_west, 0.0):
                                self.config = 'UDN'
                            elif np.allclose(self.line_group_intersect - midline_topline_east, 0.0):
                                self.config = 'UDN'
                            break

                        elif region == 'PF':
                            # True south should land in region of interest
                            if guess is S_guess:
                                self.flip_NSEW_lookup(xpt_ID)
                                xpt2 = self.NSEW_lookup['xpt2']['coor']

                            # Determine whether SF15 or SF165 based off of intersection with core boundary
                            self.RegionLineCut = self.draw_line(xpt2['N'], {'line_group': [xptSW_limiter, xptSE_limiter]},
                                option='rho', direction='cw', show_plot=visual, text=verbose)

                            if np.allclose(self.line_group_intersect - xptSE_limiter, 0.0):
                                self.config = 'SF45'
                            elif np.allclose(self.line_group_intersect - xptSW_limiter, 0.0):
                                self.config = 'SF135'
                            break
        else:
            raise ValueError(f'# Invalid xpt_ID "{xpt_ID}" in function map_xpt(). Valid IDs are "xpt1" and "xpt2".')

    def SNL_find_NSEW(self, xpt, magx, visual=False):
        """
        Find NSEW based off primary x-point and magnetic axis,

        Parameters
        ----------
        xpt : array-like
            R, Z coordinate of the primary x-point.
        mag : array-like
            R, Z coordinate of the magnetic axis.

        Notes
        -----
        self.LineTracer_psi will contain NSEW information post call.
        """
        self.map_xpt(xpt, magx, xpt_ID='xpt1')

    def DNL_find_NSEW(self, xpt1, xpt2, magx, visual=False):
        """
        Find NSEW based off primary x-point and magnetic axis,

        Parameters
        ----------
        xpt : array-like
            R, Z coordinate of the primary x-point.
        mag : array-like
            R, Z coordinate of the magnetic axis.

        Notes
        -----
        LineTracer_psi will contain NSEW information post call.
        """
        self.map_xpt(xpt1, magx, xpt_ID='xpt1')
        self.map_xpt(xpt2, magx, xpt_ID='xpt2')

    def draw_line(self, rz_start, rz_end=None, color='purple',
                  option=None, direction=None, show_plot=False, text=False, dynamic_step=None, debug=False, Verbose=False):
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
        ynot = rz_start

        if 'point' in rz_end:
            if Verbose: 
                print('# Starting search for Point convergence...')
            test = 'point'
            xf = rz_end['point']
            yf = rz_end['point']

        elif 'line' in rz_end:
            if Verbose: 
                print('# Starting search for Line intersection...')
            test = 'line'
            try:
                self.tilt_angle = rz_end['line'][1]
                endLine = rz_end['line'][0]
            except:
                 endLine = rz_end['line']

        elif 'line_group' in rz_end:
            if Verbose: 
                print('# Starting search for Line (grouping) intersection...')
            test = 'line_group'
            line_group = rz_end['line_group']

        elif 'psi' in rz_end:
            if Verbose: 
                print('# Starting search for Psi value convergence...')
            test = 'psi'
            psi_test = rz_end['psi']

        elif 'psi_horizontal' in rz_end:
            try:
                psi_test = rz_end['psi_horizontal'][0]
                self.tilt_angle = rz_end['psi_horizontal'][1]
                test = 'psi' if self.tilt_angle != 0.0 else 'psi_horizontal'
                if Verbose: 
                    print('# Starting search for Psi value convergence...')
            except:
                psi_test = rz_end['psi_horizontal']
                self.tilt_angle = 0.0
                test = 'psi_horizontal'
                if Verbose: 
                    print('# Starting search for Psi value convergence...')

        elif 'psi_vertical' in rz_end:
            if Verbose: 
                print('# Starting search for Psi value convergence...')
            test = 'psi_vertical'
            psi_test = rz_end['psi_vertical']
        else:
            if Verbose: 
                print('rz_end type not recognized')

        count = 0
        # unpack boundaries
        rmin = self.grid.rmin
        rmax = self.grid.rmax
        zmin = self.grid.zmin
        zmax = self.grid.zmax

        self.time_in_converged = 0
        line = [Point(ynot)]

        def converged(points):
            # checks for converence of the line in various ways

            def success(message):
                # Displays a message so we know the convergence worked.
                if Verbose:
                    print('Converged via {}.'.format(message))
                    print('Iterations: ', count)
                    print('Spent {} '.format(self.time_in_converged)
                          + 'seconds checking convergence.')

            def save_line(x, y, linecolor=color, marker='.-', markersize=1.5):
                # Plots the current line segments and saves
                # it for future use
                # x: list -- r endpoints
                # y: list -- z endpoints
                line.append(Point([x[-1], y[-1]]))
                if show_plot:
                    if hasattr(self.grid, 'ax') is False:
                        self.grid.plot_data()
                    self.grid.ax.plot(x, y, '.-', linewidth=2, color=color, markersize=1.5)
                    plt.draw()
                    plt.pause(np.finfo(float).eps)

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
                intersected, segment = segment_intersect((p1, p2), edge, text)
                if intersected:
                    self.grid.ax.plot(segment[1][0], segment[1][1], marker='s', color='red', markersize=12.5, zorder=20)
                    self.grid.ax.plot(segment[1][0], segment[1][1], marker='o', color='black', markersize=8, zorder=21)
                    raise ValueError(f'# LineTracing Error: The line missed the intended target of type "{test}" and intersected the boundary instead.')

            # check if any point is close enough to the endpoint
            # this currently works by checking all the intermediate points
            # with the termination point.
            # TODO: develope a method to check if the point is close
            # to a line defined by two points.
            if test == 'point':
                # Check if any of our points are withing the given
                # tol value in both the x & y directions, and check if
                # we have at least 5 segments in making up our line.
                if (any(abs(points[0] - xf) < self.tol)
                        and any(abs(points[1] - yf) < self.tol)):
                    success('endpoint')
                    save_line([points[0][0], xf], [points[1][0], yf])
                    return True
            elif test == 'line':
                p1 = (points[0][0], points[1][0])
                p2 = (points[0][-1], points[1][-1])

                intersected, segment = segment_intersect((p1, p2), endLine, text)

                if intersected:
                    success('line crossing')
                    save_line((p1[0], segment[1][0]), (p1[1], segment[1][1]))
                    return True

            elif test == 'line_group':
                # endpoints define the latest line segment
                # TODO change the point criteria
                # p1 = ( x1, y1 )
                p1 = (points[0][0], points[1][0])
                # p2 = ( x2, y2 )
                p2 = (points[0][-1], points[1][-1])

                for L in line_group:
                    result = test2points(p1, p2, L)
                    intersected, segment = segment_intersect((p1, p2), L, text)

                    if result and intersected:
                        success('line crossing')
                        save_line((p1[0], segment[1][0]), (p1[1], segment[1][1]))
                        self.line_group_intersect = L
                        return True

            elif test in ['psi', 'psi_sloped']:
                x1, y1 = points[0][0], points[1][0]
                x2, y2 = points[0][-1], points[1][-1]

                psi1 = self.grid.get_psi(x1, y1)
                psi2 = self.grid.get_psi(x2, y2)

                if (psi1 - psi_test) * (psi2 - psi_test) < 0:
                    success('psi test')
                    def f(x):
                        # must manually calculate y each time we stick it into
                        # the line of interest
                        y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
                        return psi_test - self.grid.get_psi(x, y)

                    sol = root_scalar(f, bracket=[x1, x2])
                    r_psi = sol.root
                    z_psi = (y2 - y1) / (x2 - x1) * (r_psi - x1) + y1
                    if text:
                        if Verbose: print('[x1, y1]: [{}, {}]'.format(x1, y1))
                        if Verbose: print('[x2, y2]: [{}, {}]'.format(x2, y2))

                    save_line([x1, r_psi], [y1, z_psi])
                    return True

            elif test == 'psi_horizontal':
                x1, y1 = points[0][0], points[1][0]
                x2, y2 = points[0][-1], points[1][-1]

                psi1 = self.grid.get_psi(x1, y1)
                psi2 = self.grid.get_psi(x2, y2)

                if (psi1 - psi_test) * (psi2 - psi_test) < 0:
                    success('horizontal psi integration')

                    def fend_R(x):
                        return self.grid.get_psi(x, y2) - psi_test
                    sol = root_scalar(fend_R, bracket=[x1, x2])
                    r_psi = sol.root
                    save_line([x1, r_psi], [y1, y1])
                    return True

            elif test == 'psi_vertical':
                x1, y1 = points[0][0], points[1][0]
                x2, y2 = points[0][-1], points[1][-1]

                psi1 = self.grid.get_psi(x1, y1)
                psi2 = self.grid.get_psi(x2, y2)

                if (psi1 - psi_test) * (psi2 - psi_test) < 0:
                    success('vertical psi integration')

                    def fend_Z(y):
                        return self.grid.get_psi(x2, y) - psi_test
                    sol = root_scalar(fend_Z, bracket=[y1, y2])
                    z_psi = sol.root
                    save_line([x1, x1], [y1, z_psi])
                    return True
            else:
                print('Error: No termination criteria specified.')

            # if convergence didn't occur
            if count > 0:
                # plot the line like normal
                save_line([points[0][0], points[0][-1]], [points[1][0], points[1][-1]])
                # print('Plotting...')

            t2 = time()
            self.time_in_converged += t2 - t1

        # keep track of the line segment generated via 'points' array.
        # This stores two points that give us a single line segment.
        points = np.zeros([2, 2])  # initial length is arbitrary
        start = time()

        # size for each line segment
        told, tnew = 0, self.dt

        # Use minimum of self.dt or dynamic_step value if provided.
        if dynamic_step:
            dt = np.amin([self.dt, dynamic_step])
            if dt < self.dt:
                if Verbose: 
                    print('Using dynamic_step value!\n')
        else:
            dt = self.dt

        if Verbose: 
            print('# Tracing line', end='')

        while not converged(points):
            t_span = (told, tnew)
            # solve the system of differential equations
            sol = solve_ivp(self.function, t_span, ynot, method='LSODA',
                            first_step=self.first_step, max_step=self.max_step,
                            rtol=1e-13, atol=1e-12)
            # unpack
            x0, xf = sol.y[0][0], sol.y[0][-1]
            y0, yf = sol.y[1][0], sol.y[1][-1]

            ynot = [xf, yf]
            told, tnew = tnew, tnew + dt

            # TODO: reduce the list passed in to only be the endpoints
            # of the line, and not include all the intermediate points
            # currently there is an error in the calculation of midpoints
            # for the sepratrix, when we try to truncate the points list.
            # Occurs with point converence.
            points = ((x0, xf), (y0, yf))

            if count > 3500:
                if Verbose:
                    print('did not converge, exiting...')
                    print('Iterations: ', count)
                end = time()
                if Verbose: print('Took {} '.format(end - start)
                      + 'seconds trying Lto converge.')
                break
            count += 1
            if Verbose: print('.', end='', flush=True)
        end = time()

        if Verbose:
            print('Drew for {} seconds\n'.format(end - start))
        print('')
        return Line(line)

    def PsiCostFunc(self, xy):
        x, y = xy
        return self.grid.get_psi(x, y)
