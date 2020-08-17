"""
SF135.py

Description:
    SF135 configuration class.

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


class SF135(TopologyUtils):
    def __init__(self, Ingrid_obj, config):
        TopologyUtils.__init__(self, Ingrid_obj, config)

    def concat_grid(self):
        """
        Concatenate all local grids on individual patches into a single
        array with branch cuts
        Parameters:
        ----------
            config : str
                Type of SNL grid to concat.
        """
        patch_matrix = self.patch_matrix

        for patch in self.patches.values():
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

        # Total number of poloidal indices in all subgrids.
        np_total1 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:5]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total1 = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        # Total number of poloidal indices in all subgrids.
        np_total2 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:11]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total2 = int(np.sum([patch[7].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        rm1 = np.zeros((np_total1, nr_total1, 5), order = 'F')
        zm1  = np.zeros((np_total1, nr_total1, 5), order = 'F')
        rm2 = np.zeros((np_total2, nr_total2, 5), order = 'F')
        zm2  = np.zeros((np_total2, nr_total2, 5), order = 'F')

        ixcell = 0
        jycell = 0

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')
        for ixp in range(1, 5):

            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

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
                            rm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')

        ixcell = 0
        jycell = 0

        for ixp in range(7, 11):

            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch.
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm1)):
            rm1[i] = rm1[i][::-1]
        for i in range(len(zm1)):
            zm1[i] = zm1[i][::-1]
        for i in range(len(rm2)):
            rm2[i] = rm2[i][::-1]
        for i in range(len(zm2)):
            zm2[i] = zm2[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb1 = len(rm1) - 2
        ixlb1 = 0
        ixrb2 = len(rm2) - 2
        ixlb2 = 0

        rm1 = self.add_guardc(rm1, ixlb1, ixrb1)
        zm1 = self.add_guardc(zm1, ixlb1, ixrb1)
        rm2 = self.add_guardc(rm2, ixlb2, ixrb2)
        zm2 = self.add_guardc(zm2, ixlb2, ixrb2)

        self.rm = np.concatenate((rm1, rm2))
        self.zm = np.concatenate((zm1, zm2))

        try:
            debug = self.settings['DEBUG']['visual']['gridue']
        except:
            debug = False

        if debug:
            self.animate_grid()

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

        xpt1 = self.eq.NSEW_lookup['xpt1']['coor']
        xpt2 = self.eq.NSEW_lookup['xpt2']['coor']

        magx = np.array([self.settings['grid_params']['rmagx'] + self.settings['grid_params']['patch_generation']['rmagx_shift'], \
            self.settings['grid_params']['zmagx'] + self.settings['grid_params']['patch_generation']['zmagx_shift']])

        psi_max_west = self.settings['grid_params']['psi_max_west']
        psi_max_east = self.settings['grid_params']['psi_max_east']
        psi_core = self.settings['grid_params']['psi_min_core']
        psi_pf_1 = self.settings['grid_params']['psi_min_pf']
        psi_pf_2 = self.settings['grid_params']['psi_pf2']
        psi_separatrix_2 = Point(xpt2['center']).psi(self)

        if self.settings['limiter']['use_limiter']:
            WestPlate1 = self.parent.limiter_data.copy()
            WestPlate2 = self.parent.limiter_data.copy()

            EastPlate1 = self.parent.limiter_data.copy()
            EastPlate2 = self.parent.limiter_data.copy()

        else:
            WestPlate1 = Line([Point(i) for i in self.plate_W1])
            WestPlate2 = Line([Point(i) for i in self.plate_W2])

            EastPlate1 = Line([Point(i) for i in self.plate_E1])
            EastPlate2 = Line([Point(i) for i in self.plate_E2])

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        west_midLine = Line([LHS_Point, RHS_Point])
        # inner_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])
        # outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        F1_E = self.eq.draw_line(xpt1['N'], {'psi' : psi_core}, option='rho', direction='cw',
            show_plot=visual, text=verbose)
        C1_W = F1_E.reverse_copy()

        C1_N = self.eq.draw_line(xpt1['NW'], {'line' : west_midLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        C2_S = C1_N.reverse_copy()

        C1_S = self.eq.draw_line(C1_W.p[0], {'line' : west_midLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()

        D1_N = self.eq.draw_line(C1_N.p[-1], {'line' : topLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        D2_S = D1_N.reverse_copy()

        D1_S = self.eq.draw_line(C1_S.p[0], {'line' : topLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()

        F2_S = self.eq.draw_line(xpt1['NE'], {'line' : east_midLine}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        F1_N = F2_S.reverse_copy()

        F1_S = self.eq.draw_line(F1_E.p[-1], {'line' : east_midLine}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)

        E2_S = self.eq.draw_line(F2_S.p[-1], {'line' : topLine}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        E1_N = E2_S.reverse_copy()

        # Save some key-strokes...
        d = self.eq.draw_line

        E1_S = d(F1_S.p[-1], {'line' : topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)

        B1_E = d(xpt1['S'], {'psi' : psi_pf_1}, option='rho', direction='cw', show_plot=visual, text=verbose)
        G1_W = B1_E.reverse_copy()

        G1_N = d(xpt1['SE'], {'line' : EastPlate1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        G2_S = G1_N.reverse_copy()

        G1_S = d(G1_W.p[0], {'line' : EastPlate1}, option='theta', direction='cw', 
            show_plot=visual, text=verbose).reverse_copy()

        H1_N__B1_N = d(xpt1['SW'], {'line' : WestPlate1}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()

        H1_S__B1_S = d(B1_E.p[-1], {'line' : WestPlate1}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)

        B2_N__C2_N = d(xpt2['NW'], {'line' : west_midLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        D2_N = d(B2_N__C2_N.p[-1], {'line' : topLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        D3_S = D2_N.reverse_copy()

        E2_N = d(D2_N.p[-1], {'line' : east_midLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        E3_S = E2_N.reverse_copy()

        F2_N__G2_N = d(E2_N.p[-1], {'line' : EastPlate1}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        G2_W = d(xpt1['E'], {'line' : F2_N__G2_N}, option='rho', direction='ccw',
            show_plot=visual, text=verbose)
        F2_E = G2_W.reverse_copy()

        F2_N, G2_N = F2_N__G2_N.split(G2_W.p[-1], add_split_point=True)

        F3_S, G3_S = F2_N.reverse_copy(), G2_N.reverse_copy()

        C2_W = d(xpt1['W'], {'line' : B2_N__C2_N}, option='rho', direction='ccw',
            show_plot=visual, text=verbose)
        B2_E = C2_W.reverse_copy()

        B2_N, C2_N = B2_N__C2_N.split(C2_W.p[-1], add_split_point=True)
        B3_S, C3_S = B2_N.reverse_copy(), C2_N.reverse_copy()

        H2_E = d(xpt2['N'], {'line' : H1_N__B1_N}, option='rho', direction='cw',
            show_plot=visual, text=verbose)
        B2_W = H2_E.reverse_copy()

        H1_E = d(H2_E.p[-1], {'line' : H1_S__B1_S}, option='rho', direction='cw',
            show_plot=visual, text=verbose)
        B1_W = H1_E.reverse_copy()

        H1_N, B1_N = H1_N__B1_N.split(H2_E.p[-1], add_split_point=True)
        H2_S, B2_S = H1_N.reverse_copy(), B1_N.reverse_copy()

        B1_S, H1_S = H1_S__B1_S.split(H1_E.p[-1], add_split_point=True)

        H3_S = d(xpt2['NE'], {'line' : WestPlate1}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        H2_N = H3_S.reverse_copy()

        I3_W = d(xpt2['E'], {'psi' : psi_max_east}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        H3_E = I3_W.reverse_copy()

        H3_N = d(H3_E.p[0], {'line' : WestPlate1}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        I3_N = d(H3_E.p[0], {'line' : EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        I2_N = d(xpt2['SE'], {'line' : EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        I3_S = I2_N.reverse_copy()

        A2_E__A1_E = d(xpt2['S'], {'psi' : psi_pf_2}, option='rho', direction='cw',
            show_plot=visual, text=verbose)

        A2_E, A1_E = A2_E__A1_E.split(A2_E__A1_E.p[len(A2_E__A1_E.p)//2], add_split_point=True)
        I1_W, I2_W = A1_E.reverse_copy(), A2_E.reverse_copy()

        A1_S = d(A1_E.p[-1], {'line' : WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        I1_S = d(A1_E.p[-1], {'line' : EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()

        A2_S = d(A1_E.p[0], {'line' : WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        A1_N = A2_S.reverse_copy()

        I1_N = d(A1_E.p[0], {'line' : EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        I2_S = I1_N.reverse_copy()

        A3_S = d(xpt2['SW'], {'line' : WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)
        A2_N = A3_S.reverse_copy()

        B3_W = d(xpt2['W'], {'psi' : psi_max_west}, option='rho', direction='ccw',
            show_plot=visual, text=verbose)
        A3_E = B3_W.reverse_copy()

        A3_N = d(A3_E.p[0], {'line' : WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()

        B3_N__C3_N = d(B3_W.p[-1], {'line' : west_midLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        C3_W = d(B2_N.p[-1], {'line' : B3_N__C3_N}, option='rho', direction='ccw',
            show_plot=visual, text=verbose)
        B3_E = C3_W.reverse_copy()

        B3_N, C3_N = B3_N__C3_N.split(C3_W.p[-1], add_split_point=True)

        D3_N = d(C3_N.p[-1], {'line' : topLine}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        E3_N = d(D3_N.p[-1], {'line' : east_midLine}, option='theta', direction='cw')

        F3_N__G3_N = d(E3_N.p[-1], {'line' : EastPlate1}, option='theta', direction='cw')

        G3_W = d(G2_W.p[-1], {'line' : F3_N__G3_N}, option='rho', direction='ccw',
            show_plot=visual, text=verbose)
        F3_E = G3_W.reverse_copy()

        F3_N, G3_N = F3_N__G3_N.split(G3_W.p[-1], add_split_point=True)

        C3_E = d(C2_N.p[-1], {'psi_horizontal' : psi_max_west}, option='z_const', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        C2_E = d(C2_N.p[-1], {'psi_horizontal' : 1.0}, option='z_const', direction='cw',
            show_plot=visual, text=verbose)
        C1_E = d(C2_E.p[-1], {'psi_horizontal' : psi_core}, option='z_const', direction='cw')

        D3_W, D2_W, D1_W = C3_E.reverse_copy(), C2_E.reverse_copy(), C1_E.reverse_copy()

        E3_E = d(E2_N.p[-1], {'psi_horizontal' : psi_max_west}, option='z_const', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()
        E2_E = d(E2_N.p[-1], {'psi_horizontal' : 1.0}, option='z_const', direction='ccw',
            show_plot=visual, text=verbose)
        E1_E = d(E2_E.p[-1], {'psi_horizontal' : psi_core}, option='z_const', direction='ccw')

        F3_W, F2_W, F1_W = E3_E.reverse_copy(), E2_E.reverse_copy(), E1_E.reverse_copy()

        D3_E = d(D2_N.p[-1], {'psi_vertical' : psi_max_west}, option='r_const', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()
        D2_E = d(D2_N.p[-1], {'psi_vertical' : 1.0}, option='r_const', direction='ccw',
            show_plot=visual, text=verbose)
        D1_E = d(D2_E.p[-1], {'psi_vertical' : psi_core}, option='r_const', direction='ccw')

        E3_W, E2_W, E1_W = D3_E.reverse_copy(), D2_E.reverse_copy(), D1_E.reverse_copy()


        A1_W = trim_geometry(WestPlate1, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate1, A2_S.p[-1], A2_N.p[0])
        A3_W = trim_geometry(WestPlate1, A3_S.p[-1], A3_N.p[0])

        G1_E = trim_geometry(EastPlate1, G1_N.p[-1], G1_S.p[0])
        G2_E = trim_geometry(EastPlate1, G2_N.p[-1], G2_S.p[0])
        G3_E = trim_geometry(EastPlate1, G3_N.p[-1], G3_S.p[0])

        I1_E = trim_geometry(WestPlate2, I1_N.p[-1], I1_S.p[0])
        I2_E = trim_geometry(WestPlate2, I2_N.p[-1], I2_S.p[0])
        I3_E = trim_geometry(WestPlate2, I3_N.p[-1], I3_S.p[0])

        H1_W = trim_geometry(EastPlate2, H1_S.p[-1], H1_N.p[0])
        H2_W = trim_geometry(EastPlate2, H2_S.p[-1], H2_N.p[0])
        H3_W = trim_geometry(EastPlate2, H3_S.p[-1], H3_N.p[0])


        # ============== Patch A1 ==============
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'A1', platePatch = True, plateLocation = 'W')
        # ============== Patch A2 ==============
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'A2', platePatch = True, plateLocation = 'W')
        # ============== Patch A3 ==============
        A3 = Patch([A3_N, A3_E, A3_S, A3_W], patchName = 'A3', platePatch = True, plateLocation = 'W')


        # ============== Patch B1 ==============
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'B1')
        # ============== Patch B2 ==============
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'B2')
        # ============== Patch B3 ==============
        B3 = Patch([B3_N, B3_E, B3_S, B3_W], patchName = 'B3')

        # ============== Patch C1 ==============
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'C1')
        # ============== Patch C2 ==============
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'C2')
        # ============== Patch C3 ==============
        C3 = Patch([C3_N, C3_E, C3_S, C3_W], patchName = 'C3')

        # ============== Patch D1 ==============
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'D1')
        # ============== Patch D2 ==============
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'D2')
        # ============== Patch D3 ==============
        D3 = Patch([D3_N, D3_E, D3_S, D3_W], patchName = 'D3')

        # ============== Patch E1 ==============
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'E1')
        # ============== Patch E2 ==============
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'E2')
        # ============== Patch E3 ==============
        E3 = Patch([E3_N, E3_E, E3_S, E3_W], patchName = 'E3')

        # ============== Patch F1 ==============
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'F1')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'F2')
        # ============== Patch F3 ==============
        F3 = Patch([F3_N, F3_E, F3_S, F3_W], patchName = 'F3')

        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patchName = 'G1', platePatch = True, plateLocation = 'E')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patchName = 'G2', platePatch = True, plateLocation = 'E')
        # ============== Patch G3 ==============
        G3 = Patch([G3_N, G3_E, G3_S, G3_W], patchName = 'G3', platePatch = True, plateLocation = 'E')

        # ============== Patch H1 ==============
        H1 = Patch([H1_N, H1_E, H1_S, H1_W], patchName = 'H1', platePatch = True, plateLocation = 'W')
        # ============== Patch H2 ==============
        H2 = Patch([H2_N, H2_E, H2_S, H2_W], patchName = 'H2', platePatch = True, plateLocation = 'W')
        # ============== Patch H3 ==============
        H3 = Patch([H3_N, H3_E, H3_S, H3_W], patchName = 'H3', platePatch = True, plateLocation = 'W')

        # ============== Patch I1 ==============
        I1 = Patch([I1_N, I1_E, I1_S, I1_W], patchName = 'I1', platePatch = True, plateLocation = 'E')
        # ============== Patch I2 ==============
        I2 = Patch([I2_N, I2_E, I2_S, I2_W], patchName = 'I2', platePatch = True, plateLocation = 'E')
        # ============== Patch I3 ==============
        I3 = Patch([I3_N, I3_E, I3_S, I3_W], patchName = 'I3', platePatch = True, plateLocation = 'E')

        patches = [A3, A2, A1, B3, B2, B1, C3, C2, C1, D3, D2, D1, E3, E2, E1,
                   F3, F2, F1, G3, G2, G1, H3, H2, H1, I3, I2, I1]

        self.patches = {}
        for patch in patches:
            patch.parent = self
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch


    def AdjustPatch(self,patch):
        xpt1 = Point(self.eq.NSEW_lookup['xpt1']['coor']['center'])
        xpt2 = Point(self.eq.NSEW_lookup['xpt2']['coor']['center'])

        tag = patch.get_tag()
        if tag == 'B2':
            patch.adjust_corner(xpt1, 'SE')
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'B1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'C2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'C1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'F2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'F1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'G2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'G1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'A3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'B3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'A2':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'I2':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'I3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'H3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'H2':
            patch.adjust_corner(xpt2, 'NE')


    def GroupPatches(self):
        # p = self.patches
        # self.PatchGroup = {'SOL' : [], 
        # 'CORE' : (p['ICB'], p['ICT'], p['OCT'], p['OCB']), 
        # 'PF' : (p['IPF'], p['OPF'])}
        pass

    def SetupPatchMatrix(self):
        # p = self.patches
        # self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
        #                 [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
        #                 [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
        #                 [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
        #                 ]
        pass
