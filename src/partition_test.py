from Ingrid import Ingrid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
try:
    style.use('fast')
except:
    pass
from matplotlib.patches import Polygon
from geometry import Point, Line, SNL_Patch, DNL_Patch, segment_intersect
import pathlib

def partition_domain(Ingrid_obj):
    """
    partition_domain 

        Description:
            Partition the domain into regions in order to identify whether snowflake plus or minus
            configuration.
    """

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
                if (p.x <= rmid) and (p.y >= zmid):
                    lookup[(p.x, p.y)] = p
                    point_list.append(p)

        quadrant_boundary = Line(point_list)
        quadrant_boundary.plot('dodgerblue')

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


    limiter = sort_limiter(Line([Point(p) for p in zip(grid.parent.OMFIT_psi['RLIM'], grid.parent.OMFIT_psi['ZLIM'])]).fluff_copy())
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

    grid.parent.FindXPoint2(0.95, -0.68)
    xpt2_coor = (grid.yaml['grid_params']['rxpt2'], grid.yaml['grid_params']['zxpt2'])
    ax.plot(xpt2_coor[0], xpt2_coor[1], 'x', color = 'red')
    
    if (pf_polygon.get_path().contains_point(xpt2_coor)):
        print("# Secondary x-point is contained in primary private flux region.")
    else:
        print("# Secondary x-point is NOT contained in primary private flux region.")



USN_case = "../Parameter Files/USN_YAML_EXAMPLE.yml"
DIIID_case = "../Parameter Files/DIIID_SNL.yml"
SAS_case = "/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS1_modif.yml"
ADX_case = "../Parameter Files/ADX.yml"

fname = ADX_case

ADX = Ingrid(InputFile=fname)
ADX.Setup()
partition_domain(ADX)
plt.ioff()
plt.show()