"""
SNL.py

Description:
    SNL configuration class.

Created: June 18, 2020

"""
import numpy as np
import matplotlib
import pathlib
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt
from TopologyUtils import TopologyUtils
from geometry import Point, Line, Patch, trim_geometry


class SNL(TopologyUtils):
    """
    The SNL (Single-Null) class is the parent class for both upper-single null (USN)
    and lower-single null (LSN) configurations.
    This base class handles the formatting and plotting of data obtained from an LSN or USN
    object.
    Parameter:
        - INGRID_object : Ingrid class object
        All SNL objects are children of the main Ingrid class. INGRID_object provides
        information such as YAML data, efit_psi, psi_norm, and plate_data.
    @author: garcia299
    """

    def __init__(self, Ingrid_obj, config):
        TopologyUtils.__init__(self, Ingrid_obj, config)


    def concat_grid(self, Verbose=False):
        """
        Concatenate all local grids on individual patches into a single
        array with branch cuts
        Parameters:
        ----------
            Verbose : bool
                Verbose flag indicator
        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches.values():
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1


        # Total number of poloidal indices in subgrid.
        np_total = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:-1]])) + 2
        nr_total = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:3]])) + 2

        rm = np.zeros((np_total, nr_total, 5), order = 'F')
        zm = np.zeros((np_total, nr_total, 5), order = 'F')

        ixcell = 0
        jycell = 0

        # Iterate over all the patches in our SNL configuration (we exclude guard cells denoted by '[None]')
        for ixp in range(1, 7):

            nr_sum = 0
            for jyp in range(1, 3):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]
                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch.
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm)):
            rm[i] = rm[i][::-1]
        for i in range(len(zm)):
            zm[i] = zm[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb)
        self.zm = self.add_guardc(zm, ixlb, ixrb)

        try:
            debug = self.settings['DEBUG']['visual']['gridue']
        except:
            debug = False

        if debug:
            self.animate_grid()


    def add_guardc(self, cell_map, ixlb, ixrb, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, boundary):
            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # TODO: Edit the code to reflect this at some point so the next reader is not overwhelmed.
            if boundary == 'left':
                ixn = ix + 1
                iyn = iy
                cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][2])
                cell_map[ix][iy][2] = cell_map[ixn][iyn][1]
                cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][4])
                cell_map[ix][iy][4] = cell_map[ixn][iyn][3]
                cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

            elif boundary == 'right':
                ixn = ix - 1
                iyn = iy
                cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][1])
                cell_map[ix][iy][1] = cell_map[ixn][iyn][2]
                cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][3])
                cell_map[ix][iy][3] = cell_map[ixn][iyn][4]
                cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

            elif boundary == 'bottom':
                ixn = ix
                iyn = iy + 1
                cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][3])
                cell_map[ix][iy][3] = cell_map[ixn][iyn][1]
                cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][4])
                cell_map[ix][iy][4] = cell_map[ixn][iyn][2]
                cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])
            elif boundary == 'top':
                ixn = ix
                iyn = iy - 1
                cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][1])
                cell_map[ix][iy][1] = cell_map[ixn][iyn][3]
                cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][2])
                cell_map[ix][iy][2] = cell_map[ixn][iyn][4]
                cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

            return cell_map

        np = len(cell_map) - 2
        nr = len(cell_map[0]) - 2

        for iy in range(1, nr + 1):
            ix = ixlb
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'top')

        return cell_map


    def AdjustPatch(self,patch):
        primary_xpt = Point(self.eq.NSEW_lookup['xpt1']['coor']['center'])

        tag = patch.get_tag()
        if tag == 'A2':
            patch.adjust_corner(primary_xpt, 'SE')
        elif tag == 'A1':
            patch.adjust_corner(primary_xpt, 'NE')
        elif tag == 'B2':
            patch.adjust_corner(primary_xpt, 'SW')
        elif tag == 'B1':
            patch.adjust_corner(primary_xpt, 'NW')
        elif tag == 'E1':
            patch.adjust_corner(primary_xpt, 'NE')
        elif tag == 'E2':
            patch.adjust_corner(primary_xpt, 'SE')
        elif tag == 'F1':
            patch.adjust_corner(primary_xpt, 'NW')
        elif tag == 'F2':
            patch.adjust_corner(primary_xpt, 'SW')


    # def construct_grid(self, np_cells = 1, nr_cells = 1,Verbose=False,ShowVertices=False,RestartScratch=False,OptionTrace='theta',ExtraSettings={},ListPatches='all', Enforce=True):

    #     # Straighten up East and West segments of our patches,
    #     # Plot borders and fill patches.
        
    #     self.RefreshSettings()

    #     if Verbose: print('Construct Grid')
    #     try:
    #         visual = self.settings['DEBUG']['visual']['subgrid']
    #     except:
    #         visual = False
    #     try:
    #         verbose = self.settings['DEBUG']['verbose']['grid_generation']
    #     except:
    #         verbose = False

    #     verbose=Verbose or verbose
        
            
    #     print('>>> Patches:', [k for k in self.patches.keys()])
    #     if RestartScratch:
    #         self.CurrentListPatch={}
    
    #     for name, patch in self.patches.items():
            
    #         if self.CorrectDistortion.get(name) is not None:
    #            patch.CorrectDistortion=self.CorrectDistortion.get(name)
    #         elif self.CorrectDistortion.get('all') is not None:
    #             patch.CorrectDistortion=self.CorrectDistortion.get('all')
    #         else:
    #             patch.CorrectDistortion={'Active':False}
    #         if (ListPatches=='all' and patch not in self.CurrentListPatch) or (ListPatches!='all' and name in ListPatches):
    #             self.SetPatchBoundaryPoints(patch)
    #             (nr_cells,np_cells)=self.GetNpoints(patch, Enforce=Enforce)
    #             (_radial_f,_poloidal_f)=self.GetFunctions(patch,ExtraSettings=ExtraSettings,Enforce=Enforce)
    #             print('>>> Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(name, np_cells, nr_cells, inspect.getsource(_poloidal_f), inspect.getsource(_radial_f)))
    #             patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f,verbose = verbose, visual = visual,ShowVertices=ShowVertices,OptionTrace=OptionTrace)
    #             self.AdjustPatch(patch)
    #             patch.plot_subgrid()
    #             self.CurrentListPatch[name] = patch


    #     if all(['cell_grid' in patch.__dict__ for patch in self.patches.values()]):
    #         self.concat_grid()
    #         self.set_gridue()
        

    def construct_patches(self):
        """
        Draws lines and creates patches for both USN and LSN configurations.

        Patch Labeling Key:
            I: Inner,
            O: Outer,
            DL: Divertor Leg,
            PF: Private Flux,
            T: Top,
            B: Bottom,
            S: Scrape Off Layer,
            C: Core.
        """
        # TODO: Create a 'lookup' procedure for determining line drawing
        #       orientations and inner-outer locations.

        try:
            visual = self.settings['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False
        try:
            verbose = self.settings['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False
        try:
            inner_tilt = self.settings['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0
        try:
            outer_tilt = self.settings['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0

        self.RefreshSettings()

        if self.settings['limiter']['use_limiter']:
            WestPlate = self.parent.limiter_data.copy()
            EastPlate = self.parent.limiter_data.copy()

        else:
            WestPlate = self.plate_data['plate_W1']
            EastPlate = self.plate_data['plate_E1']

        xpt = self.eq.NSEW_lookup['xpt1']['coor']
        magx = np.array([self.settings['grid_params']['rmagx'] + self.settings['grid_params']['patch_generation']['rmagx_shift'], \
            self.settings['grid_params']['zmagx'] + self.settings['grid_params']['patch_generation']['zmagx_shift']])

        psi_max = self.settings['grid_params']['psi_max']
        psi_min_core = self.settings['grid_params']['psi_min_core']
        psi_min_pf = self.settings['grid_params']['psi_min_pf']

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        inner_midLine = Line([LHS_Point, RHS_Point])
        inner_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        outer_midLine = Line([LHS_Point, RHS_Point])
        outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Drawing Separatrix
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

        # Drawing Lower-SNL region
        if self.settings['grid_params']['patch_generation']['use_NW']:
            tilt = self.settings['grid_params']['patch_generation']['NW_adjust']
            xptW_psiMax = self.eq.draw_line(rotate(xpt['W'], tilt, xpt['center']), {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        else:
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        if self.settings['grid_params']['patch_generation']['use_NE']:
            tilt = self.settings['grid_params']['patch_generation']['NE_adjust']
            xptE_psiMax = self.eq.draw_line(rotate(xpt['E'], tilt, xpt['center']), {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        xpt_WestPlate = self.eq.draw_line(xpt['SW'], {'line' : WestPlate}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xpt_EastPlate = self.eq.draw_line(xpt['SE'], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        iPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : WestPlate}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_WestPlate = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : WestPlate},option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        oPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_EastPlate = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        imidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', \
            direction = 'cw', show_plot = visual, text = verbose)
        
        omidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', \
            direction = 'ccw', show_plot = visual, text = verbose)

        # Integrating horizontally along mid-line towards psiMax and psiMinCore

        imidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_max, inner_tilt)}, option = 'z_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)
        imidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_min_core, inner_tilt)}, option = 'z_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        omidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_max, outer_tilt)}, option = 'z_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        omidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_min_core, outer_tilt)}, option = 'z_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)

        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)

        # A1 Patch
        location = 'W'
        A2_N = iPsiMax_TP.reverse_copy()
        A2_S = xpt_WestPlate
        A2_E = xptW_psiMax.reverse_copy()
        # =====================================================================================
        # Trimming the target_plate to conform to the patch boundary.
        # -------------------------------------------------------------------------------------
        # Recall ITP has a clockwise orientation.
        #
        # The inner 'split' trims all Point objects BEFORE the point of intersection of ITP
        # and A2_S. Call this new Line object Line_A.
        #
        # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
        # and A2_N. This new Line object is the plate facing boundary of the Patch.
        # =====================================================================================
        A2_W = (WestPlate.split(A2_S.p[-1])[1]).split(A2_N.p[0], add_split_point = True)[0]
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # A1 Patch
        location = 'W'
        A1_N = A2_S.reverse_copy()
        
        A1_S = psiMinPF_WestPlate
        A1_E = xptS_psiMinPF
        A1_W = (WestPlate.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # B2 Patch
        
        B2_N = self.eq.draw_line(A2_N.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B2_S = xptNW_midLine.reverse_copy()
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B2_W = xptW_psiMax
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'ISB')

        # B1 Patch
        B1_N = B2_S.reverse_copy()
        B1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B1_W = xptN_psiMinCore.reverse_copy()
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'ICB')

        # C2 Patch
        C2_N = self.eq.draw_line(B2_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C2_S = imidLine_topLine.reverse_copy()
        C2_E = Line([C2_N.p[-1], C2_S.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'IST')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.eq.draw_line(B1_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        C1_E = Line([C1_N.p[-1], C1_S.p[0]])
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'ICT')

        # F2 Patch
        location = 'E'
        F2_N = oPsiMax_TP
        F2_S = xpt_EastPlate.reverse_copy()
        F2_E = (EastPlate.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point = True)[0]
        F2_W = xptE_psiMax
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # F1 Patch
        location = 'E'
        F1_N = F2_S.reverse_copy()
        F1_S = psiMinPF_EastPlate.reverse_copy()
        F1_E = (EastPlate.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point = True)[0]
        F1_W = xptS_psiMinPF.reverse_copy()
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # E2 Patch
        E2_N = self.eq.draw_line(F2_N.p[0], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        E2_S = xptNE_midLine
        E2_E = xptE_psiMax.reverse_copy()
        E2_W = Line([E2_S.p[-1], E2_N.p[0]])
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'OSB')

        # E1 Patch
        E1_N = E2_S.reverse_copy()
        E1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E1_E = xptN_psiMinCore
        E1_W = Line([E1_S.p[-1], E1_N.p[0]])
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'OCB')

        # D2 Patch
        D2_N = self.eq.draw_line(E2_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D2_S = omidLine_topLine
        D2_E = Line([D2_N.p[-1], D2_S.p[0]])
        D2_W = Line([D2_S.p[-1], D2_N.p[0]])
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'OST')

        # D1 Patch
        D1_N = D2_S.reverse_copy()
        D1_S = self.eq.draw_line(E1_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_E = Line([D1_N.p[-1], D1_S.p[0]])
        D1_W = Line([D1_S.p[-1], D1_N.p[0]])
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'OCT')

        patches = [A2, A1, B2, B1, C2, C1, D2, D1, E2, E1, F2, F1]

        self.patches = {}
        for patch in patches:
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch

    def GroupPatches(self):
        p = self.patches
        self.PatchGroup = {'SOL' : (p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL']), 
        'CORE' : (p['ICB'], p['ICT'], p['OCT'], p['OCB']), 
        'PF' : (p['IPF'], p['OPF'])}

    def set_gridue(self):
        """
        Prepare the relevant arrays for writing to GRIDUE.
        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patches['IDL'].npol - 1
        ixpt2 = ixrb - self.patches['ODL'].npol + 1
        iyseparatrix1 = self.patches['IDL'].nrad - 1
        nxm = len(self.rm) - 2
        nym = len(self.rm[0]) - 2

        psi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        br = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bz = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bpol = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bphi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        b = np.zeros((nxm + 2, nym + 2, 5), order = 'F')

        rm = self.rm
        zm = self.zm
        rb_prod = self.efit_psi.rcenter * self.efit_psi.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.efit_psi.get_psi(_r, _z)
                    _br = self.efit_psi.get_psi(_r, _z, tag = 'vz') / _r
                    _bz = -self.efit_psi.get_psi(_r, _z, tag = 'vr') / _r
                    _bpol = np.sqrt(_br ** 2 + _bz ** 2)
                    _bphi = rb_prod / _r
                    _b = np.sqrt(_bpol ** 2 + _bphi ** 2)

                    psi[i][j][k] = _psi
                    br[i][j][k] = _br
                    bz[i][j][k] = _bz
                    bpol[i][j][k] = _bpol
                    bphi[i][j][k] = _bphi
                    b[i][j][k] = _b

        self.gridue_params = {'nxm' : nxm, 'nym' : nym, 'ixpt1' : ixpt1, 'ixpt2' : ixpt2, 'iyseptrx1' : iyseparatrix1, \
            'rm' : self.rm, 'zm' : self.zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b}