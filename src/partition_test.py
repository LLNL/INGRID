from Ingrid import Ingrid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from geometry import Point, Line, Patch, segment_intersect
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from matplotlib.patches import Polygon
import pathlib

class RegionEntered(Exception):
    def __init__(self, message, region):
        self.message=message
        self.region=region

def partition_domain(Ingrid_obj, xpt2=None, bounds = None):
    """
    partition_domain 

        Description:
            Partition the domain into regions in order to identify whether snowflake plus or minus
            configuration.
    """
    def cb(xk):
        if core_polygon.get_path().contains_point(xk):
            plt.plot(xk[0], xk[1], 'x', color='dodgerblue')
            raise RegionEntered(message='# Entered Core...', region='Core')
        elif pf_polygon.get_path().contains_point(xk):
            plt.plot(xk[0], xk[1], 'x', color='magenta')
            raise RegionEntered(message='# Entered PF...', region='PF')
    def cb_plot1(xk):
        plt.plot(xk[0], xk[1], '.', color='dodgerblue')
    def cb_plot2(xk):
        plt.plot(xk[0], xk[1], '.', color='black')

    def sort_limiter(limiter):
        """
        sort_limiter
             Description:
                Orients limiter data clockwise around the magnetic axis by inspecting 
                NW quadrant of domain.
        """
        def unit_vector(v):
            """
            Returns unit vector
            """
            return v / np.linalg.norm(v)

        def angle_between(u, v):
            """
            Compute angle in radians between vectors u and v
            """
            u_norm = unit_vector(u)
            v_norm = unit_vector(v)
            return np.arctan2(u_norm[0] * v_norm[1] - u_norm[1] * v_norm[0], \
                u_norm[0] * v_norm[0] + u_norm[0] * v_norm[1])

        rmid = (grid.psi_norm.rmin + grid.psi_norm.rmax) / 2
        zmid = (grid.psi_norm.zmin + grid.psi_norm.zmax) / 2

        lookup = {}
        point_list = []
        for p in limiter.p:
            try:
                lookup[(p.x, p.y)]
            except:
                if (p.x >= rmid) and (p.y >= zmid):
                    lookup[(p.x, p.y)] = p
                    point_list.append(p)

        quadrant_boundary = Line(point_list)
        quadrant_boundary.plot('red')

        ordered = False
        start = quadrant_boundary.p[0]
        for p in reversed(quadrant_boundary.p):
            end = p
            if end.coor != start.coor:
                break

        # Endpoints are on same vertical line
        if start.x == end.x:
            # If end point is above start point
            if start.y < end.y: 
                # Return target plate as is
                ordered = True

        # Endpoints are on same horizontal line
        elif start.y == end.y:
            # If start point is left of end point
            if start.x < end.x:
                # Return strike plate as is
                ordered = True

        else:

            # Check the angle between the plate endpoints with
            # tail fixed on the magnetic axis

            grid.parent.AutoRefineMagAxis()
            v_reference = np.array( [ grid.yaml['grid_params']['rmagx'], 
                grid.yaml['grid_params']['zmagx'] ])

            v_start = np.array( [ start.x, start.y ] ) - v_reference
            v_end = np.array( [ end.x, end.y ] ) - v_reference

            if angle_between(v_start, v_end) <= 0:
                ordered = True

        return limiter.copy() if ordered else limiter.reverse_copy()


    grid = Ingrid_obj.current_topology

    try:
        visual = grid.yaml['DEBUG']['visual']['patch_map']
    except KeyError:
        visual = False
    try:
        verbose = grid.yaml['DEBUG']['verbose']['patch_generation']
    except KeyError:
        verbose = False

    WestPlate = Line([Point(i) for i in grid.plate_data['plate_W1']['coordinates']])
    EastPlate = Line([Point(i) for i in grid.plate_data['plate_E1']['coordinates']])

    if bounds:
        rmin, rmax, zmin, zmax = bounds
        grid.parent.OMFIT_psi['RLIM'] = [rmin, rmax, rmax, rmin, rmin]
        grid.parent.OMFIT_psi['ZLIM'] = [zmax, zmax, zmin, zmin, zmax]

    limiter = sort_limiter(Line([Point(p) for p in zip(grid.parent.OMFIT_psi['RLIM'], grid.parent.OMFIT_psi['ZLIM'])]).fluff_copy(100))
    limiter.plot('dodgerblue')

    xpt = grid.eq.NSEW_lookup['xpt1']['coor']
    magx = np.array([grid.yaml['grid_params']['rmagx'], grid.yaml['grid_params']['zmagx']])

    psi_max = grid.yaml['grid_params']['psi_max']
    psi_min_core = grid.yaml['grid_params']['psi_min_core']
    psi_min_pf = grid.yaml['grid_params']['psi_min_pf']

    # Create mid-line
    LHS_Point = Point(magx[0] - 1e6, magx[1])
    RHS_Point = Point(magx[0] + 1e6, magx[1])
    midline = Line([LHS_Point, RHS_Point])

    # Generate Vertical Mid-Plane line
    Lower_Point = Point(magx[0], magx[1] - 1e6)
    Upper_Point = Point(magx[0], magx[1] + 1e6)
    topline = Line([Lower_Point, Upper_Point])

    # Drawing Separatrix
    xptNW_midLine = grid.eq.draw_line(xpt['NW'], {'line' : midline}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
    xptNE_midLine = grid.eq.draw_line(xpt['NE'], {'line' : midline}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
    midline_topline_west = grid.eq.draw_line(xptNW_midLine.p[-1], {'line' : topline}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
    midline_topline_east = grid.eq.draw_line(xptNE_midLine.p[-1], {'line' : topline}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

    xptSW_limiter = grid.eq.draw_line(xpt['SW'], {'line' : limiter}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
    xptSE_limiter = grid.eq.draw_line(xpt['SE'], {'line' : limiter}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

    core_boundary = Line((xptNW_midLine.p + midline_topline_west.p + midline_topline_east.reverse_copy().p + xptNE_midLine.reverse_copy().p))
    core_polygon = Polygon(np.column_stack(core_boundary.points()).T, fill = True, closed = True, color = 'violet', label='Core Region')

    pf_boundary = Line((xptSW_limiter.p + xptSE_limiter.p + (limiter.split(xptSE_limiter.p[-1])[1]).split(xptSW_limiter.p[0], add_split_point = True)[0].p))
    pf_polygon = Polygon(np.column_stack(pf_boundary.points()).T, fill = True, closed = True, color = 'dodgerblue', label='PF Region')

    ax = grid.eq.grid.ax
    ax.add_patch(core_polygon)
    ax.add_patch(pf_polygon)
    ax.legend()
    core_boundary.plot('black')
    xptSW_limiter.plot('black')
    xptSE_limiter.plot('black')

    rxpt2, zxpt2 = xpt2
    grid.parent.FindXPoint2(rxpt2, zxpt2)
    xpt2_coor = (grid.yaml['grid_params']['rxpt2'], grid.yaml['grid_params']['zxpt2'])
    ax.plot(xpt2_coor[0], xpt2_coor[1], 'x', color = 'red')
    
    if (pf_polygon.get_path().contains_point(xpt2_coor)):
        print("# Snowflake-minus")
    else:
        print("# Snowflake-plus")

        # Obtain candidate NSEW directions
        grid.eq.analyze_saddle(xpt2_coor, xpt_ID='xpt2')

        # Prepare minimizer for identification of SF+ type
        grid.eq._set_function('rho', 'cw')
        NS_buffer = grid.eq.NSEW_lookup['xpt2']['coor']
        N_sol = solve_ivp(grid.eq.function, (0, grid.eq.dt), NS_buffer['N'], method='LSODA',\
                first_step=grid.eq.first_step, max_step=grid.eq.max_step, rtol=1e-13, atol=1e-12).y
        S_sol = solve_ivp(grid.eq.function, (0, grid.eq.dt), NS_buffer['S'], method='LSODA',\
                first_step=grid.eq.first_step, max_step=grid.eq.max_step, rtol=1e-13, atol=1e-12).y
        N_guess = (N_sol[0][-1], N_sol[1][-1])
        S_guess = (S_sol[0][-1], S_sol[1][-1])

        minimize(grid.eq.PsiCostFunc, N_guess, method='trust-ncg', jac=grid.psi_norm.Gradient, hess=grid.psi_norm.Hessian,
            options={'initial_trust_radius' : grid.eq.eps, 'max_trust_radius' : grid.eq.dt}, callback=cb_plot1)
        minimize(grid.eq.PsiCostFunc, S_guess, method='trust-ncg', jac=grid.psi_norm.Gradient, hess=grid.psi_norm.Hessian,
            options={'initial_trust_radius' : grid.eq.eps, 'max_trust_radius' : grid.eq.dt}, callback=cb_plot2)

        for guess in [N_guess, S_guess]:
            try:
                minimize(grid.eq.PsiCostFunc, guess, method='trust-ncg', jac=grid.psi_norm.Gradient, hess=grid.psi_norm.Hessian,
                    options={'initial_trust_radius' : grid.eq.eps, 'max_trust_radius' : grid.eq.dt}, callback=cb)
            except RegionEntered as e:
                region = e.region
                if region == 'Core':
                    print('# Identified SF15')
                    # True south should land in region of interest
                    if guess is N_guess: grid.eq.flip_NSEW_lookup()
                    break
                elif region == 'PF':
                    print('# Identified SF45')
                    # True south should land in region of interest
                    if guess is N_guess: grid.eq.flip_NSEW_lookup()
                    break

def ADX_test():
    ADX_case = "../Parameter Files/ADX.yml"
    fname = ADX_case
    ADX = Ingrid(InputFile=fname)
    ADX.Setup()
    xpt2 = (0.95, -0.68)
    partition_domain(ADX, xpt2)
    plt.ioff()
    plt.show()

def SF15_test():
    SF15_case = "../Parameter Files/SF15.yml"
    fname = SF15_case
    SF15 = Ingrid(InputFile=fname)
    SF15.OMFIT_read_psi()
    SF15.calc_efit_derivs()
    SF15.read_target_plates()
    SF15.AutoRefineXPoint()
    SF15.FindXPoint(1.5, -0.62)
    SF15.FindMagAxis(1.71, 0.45)
    SF15.SetMagReference()
    SF15.calc_psinorm()
    SF15.plot_psinorm(interactive=True)
    SF15.init_LineTracing(interactive=True)
    SF15._analyze_topology()
    xpt2 = (2.1, -0.557)
    bounds = [1.0, 2.4, -1, 1.0]
    partition_domain(SF15, xpt2, bounds)
    plt.ioff()
    plt.show()

def SF45_test():
    SF45_case = "../Parameter Files/SF45.yml"
    fname = SF45_case
    SF45 = Ingrid(InputFile=fname)
    SF45.OMFIT_read_psi()
    SF45.calc_efit_derivs()
    SF45.read_target_plates()
    SF45.AutoRefineXPoint()
    SF45.SetMagReference()
    SF45.calc_psinorm()
    SF45.plot_psinorm(interactive=True)
    SF45.init_LineTracing(interactive=True)
    SF45._analyze_topology()
    xpt2 = (2.12, -1.13)
    bounds = [1.0, 2.4, -1, 1.0]
    SF45.FindXPoint(1.5, -0.47)
    SF45.FindMagAxis(1.71, 0.45)
    partition_domain(SF45, xpt2, bounds)
    plt.ioff()
    plt.show()

def SPARC_test():
    SPARC_case = "../Parameter Files/SPARC_SF.yml"
    fname = SPARC_case
    SPARC = Ingrid(InputFile=fname)
    SPARC.Setup()
    xpt2 = (1.7, -1.51)
    partition_domain(SPARC, xpt2)
    plt.ioff()
    plt.show()    

SF15_test()
ADX_test()
SPARC_test()





