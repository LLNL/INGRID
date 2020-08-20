"""
SF165.py

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


class SF165(TopologyUtils):
    def __init__(self, Ingrid_obj, config):
        TopologyUtils.__init__(self, Ingrid_obj, config)

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
            west_tilt = self.settings['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            west_tilt = 0.0
        try:
            east_tilt = self.settings['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            east_tilt = 0.0

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
        LHS_Point = Point(magx[0] - 1e6 * np.cos(west_tilt), magx[1] - 1e6 * np.sin(west_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(west_tilt), magx[1] + 1e6 * np.sin(west_tilt))
        west_midLine = Line([LHS_Point, RHS_Point])
        west_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(east_tilt), magx[1] - 1e6 * np.sin(east_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(east_tilt), magx[1] + 1e6 * np.sin(east_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])
        # outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Save some key-strokes...
        d = self.eq.draw_line

        E1_E = d(xpt1['N'], {'psi' : psi_core}, option='rho', direction='cw', show_plot=visual, text=verbose)
        H1_W = E1_E.reverse_copy()

        E1_S = d(E1_E.p[-1], {'line' : east_midLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        D1_S = d(E1_S.p[-1], {'line' : topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)

        B1_S__H1_S = d(E1_E.p[-1], {'line' : west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()
        C1_S = d(B1_S__H1_S.p[0], {'line' : topLine}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()

        H1_N__B1_N = d(xpt1['NW'], {'line' : west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)

        C1_N = d(H1_N__B1_N.p[-1], {'line' : topLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        C2_S = C1_N.reverse_copy()

        E2_S = d(xpt1['NE'], {'line' : east_midLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        E1_N = E2_S.reverse_copy()

        D2_S = d(E2_S.p[-1], {'line' : topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        D1_N = D2_S.reverse_copy()

        H2_E = d(xpt2['N'], {'line' : H1_N__B1_N}, option='rho', direction='cw', show_plot=visual, text=verbose)
        B2_W = H2_E.reverse_copy()

        H1_N, B1_N = H1_N__B1_N.split(H2_E.p[-1], add_split_point=True)
        H2_S, B2_S = H1_N.reverse_copy(), B1_N.reverse_copy()

        H1_E = d(H2_E.p[-1], {'line' : B1_S__H1_S}, option='rho', direction='cw', show_plot=visual, text=verbose)
        B1_W = H1_E.reverse_copy()

        B1_S, H1_S = B1_S__H1_S.split(H1_E.p[-1], add_split_point=True)

        F2_W__F3_W = d(xpt1['E'], {'psi' : psi_max_west}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        B3_W = d(xpt2['W'], {'psi' : psi_max_west}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        A3_E = B3_W.reverse_copy()

        A3_N = d(B3_W.p[-1], {'line' : WestPlate2}, option='theta', direction='ccw', 
            show_plot=visual, text=verbose).reverse_copy()
        B3_N = d(A3_N.p[-1], {'line' : west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        C3_N = d(B3_N.p[-1], {'line' : topLine}, option='theta', direction='cw', show_plot=visual, text=verbose)

        F3_N = d(F2_W__F3_W.p[-1], {'line' : EastPlate1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        E3_N = d(F3_N.p[0], {'line' : east_midLine}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        D3_N = d(E3_N.p[0], {'line' : topLine}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()

        B2_N = d(xpt2['NW'], {'line' : west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        B3_S = B2_N.reverse_copy()

        C2_N = d(B2_N.p[-1], {'line' : topLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        C3_S = C2_N.reverse_copy()

        D2_N = d(C2_N.p[-1], {'line' : east_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        D3_S = D2_N.reverse_copy()

        E2_N = d(D2_N.p[-1], {'line' : F2_W__F3_W}, option='theta', direction='cw', show_plot=visual, text=verbose)
        E3_S = E2_N.reverse_copy()

        F2_W, F3_W = F2_W__F3_W.split(E2_N.p[-1], add_split_point=True)
        E2_E, E3_E = F2_W.reverse_copy(), F3_W.reverse_copy()

        F2_N = d(E2_N.p[-1], {'line' : EastPlate1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        F3_S = F2_N.reverse_copy()

        F1_N = d(xpt1['SE'], {'line' : EastPlate1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        F2_S = F1_N.reverse_copy()

        G1_E = d(xpt1['S'], {'psi' : psi_pf_1}, option='rho', direction='cw', show_plot=visual, text=verbose)
        F1_W = G1_E.reverse_copy()

        G1_S = d(G1_E.p[-1], {'line' : WestPlate1}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        F1_S = d(G1_E.p[-1], {'line' : EastPlate1}, option='theta', direction='cw', 
            show_plot=visual, text=verbose).reverse_copy()

        G2_S = d(xpt1['SW'], {'line' : WestPlate1}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        G1_N = G2_S.reverse_copy()

        H3_S__G3_S = d(xpt2['NE'], {'line' : WestPlate1}, option='theta', direction='ccw', show_plot=visual, text=verbose)

        H2_W = d(xpt1['W'], {'line' : H3_S__G3_S}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        G2_E = H2_W.reverse_copy()

        H3_S, G3_S = H3_S__G3_S.split(H2_W.p[-1], add_split_point=True)
        H2_N, G2_N = H3_S.reverse_copy(), G3_S.reverse_copy()

        I3_W = d(xpt2['E'], {'psi' : psi_max_east}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        H3_E = I3_W.reverse_copy()

        I3_N = d(I3_W.p[-1], {'line' : EastPlate2}, option='theta', direction='cw', show_plot=visual, text=verbose)

        G3_N__H3_N = d(I3_W.p[-1], {'line' : WestPlate1}, option='theta', direction='ccw', 
            show_plot=visual, text=verbose).reverse_copy()

        H3_W = d(G2_N.p[-1], {'line' : G3_N__H3_N}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        G3_E = H3_W.reverse_copy()

        G3_N, H3_N = G3_N__H3_N.split(H3_W.p[-1], add_split_point=True)

        A2_E__A1_E = d(xpt2['S'], {'psi' : psi_pf_2}, option='rho', direction='cw', show_plot=visual, text=verbose)
        A2_E, A1_E = A2_E__A1_E.split(A2_E__A1_E.p[len(A2_E__A1_E.p)//2], add_split_point=True)
        I2_W, I1_W = A2_E.reverse_copy(), A1_E.reverse_copy()

        A1_S = d(A1_E.p[-1], {'line' : WestPlate2}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        I1_S = d(A1_E.p[-1], {'line' : EastPlate2}, option='theta', direction='cw', 
            show_plot=visual, text=verbose).reverse_copy()

        A2_S = d(A2_E.p[-1], {'line' : WestPlate2}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        A1_N = A2_S.reverse_copy()
        I2_S = d(A2_E.p[-1], {'line' : EastPlate2}, option='theta', direction='cw', 
            show_plot=visual, text=verbose).reverse_copy()
        I1_N = I2_S.reverse_copy()

        A3_S = d(xpt2['SW'], {'line' : WestPlate2}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        A2_N = A3_S.reverse_copy()

        I3_S = d(xpt2['SE'], {'line' : EastPlate2}, option='theta', direction='cw', 
            show_plot=visual, text=verbose).reverse_copy()
        I2_N = I3_S.reverse_copy()


        B3_E = self.eq.draw_line(B3_N.p[-1], {'psi_horizontal' : (psi_separatrix_2, west_tilt)}, option = 'z_const', direction = 'cw', 
            show_plot = visual, text = verbose)
        C3_W = B3_E.reverse_copy()

        B2_E = self.eq.draw_line(B3_E.p[-1], {'psi_horizontal' : (1.0, west_tilt)}, option='z_const', direction='cw',
            show_plot=visual, text=verbose)
        C2_W = B2_E.reverse_copy()

        B1_E = self.eq.draw_line(B2_E.p[-1], {'psi_horizontal' : (psi_core, west_tilt)}, option='z_const', direction='cw',
            show_plot=visual, text=verbose)
        C1_W = B1_E.reverse_copy()

        C3_E = self.eq.draw_line(C3_N.p[-1], {'psi_vertical' : psi_separatrix_2}, option = 'r_const', direction = 'ccw', 
            show_plot = visual, text = verbose)
        D3_W = C3_E.reverse_copy()

        C2_E = self.eq.draw_line(C2_N.p[-1], {'psi_vertical' : 1.0}, option = 'r_const', direction = 'ccw', 
            show_plot = visual, text = verbose)
        D2_W = C2_E.reverse_copy()

        C1_E = self.eq.draw_line(C1_N.p[-1], {'psi_vertical' : psi_core}, option = 'r_const', direction = 'ccw', 
            show_plot = visual, text = verbose)
        D1_W = C1_E.reverse_copy()

        D3_E = self.eq.draw_line(D3_N.p[-1], {'psi_horizontal' : (psi_separatrix_2, east_tilt)}, option = 'z_const', direction = 'ccw', 
            show_plot = visual, text = verbose)
        E3_W = D3_E.reverse_copy()

        D2_E = self.eq.draw_line(D3_E.p[-1], {'psi_horizontal' : (1.0, east_tilt)}, option='z_const', direction='ccw',
            show_plot=visual, text=verbose)
        E2_W = D2_E.reverse_copy()

        D1_E = self.eq.draw_line(D2_E.p[-1], {'psi_horizontal' : (psi_core, east_tilt)}, option='z_const', direction='ccw',
            show_plot=visual, text=verbose)
        E1_W = D1_E.reverse_copy()


        A1_W = trim_geometry(WestPlate2, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate2, A2_S.p[-1], A2_N.p[0])
        A3_W = trim_geometry(WestPlate2, A3_S.p[-1], A3_N.p[0])

        F1_E = trim_geometry(EastPlate1, F1_N.p[-1], F1_S.p[0])
        F2_E = trim_geometry(EastPlate1, F2_N.p[-1], F2_S.p[0])
        F3_E = trim_geometry(EastPlate1, F3_N.p[-1], F3_S.p[0])

        I1_E = trim_geometry(WestPlate1, I1_N.p[-1], I1_S.p[0])
        I2_E = trim_geometry(WestPlate1, I2_N.p[-1], I2_S.p[0])
        I3_E = trim_geometry(WestPlate1, I3_N.p[-1], I3_S.p[0])

        G1_W = trim_geometry(EastPlate2, G1_S.p[-1], G1_N.p[0])
        G2_W = trim_geometry(EastPlate2, G2_S.p[-1], G2_N.p[0])
        G3_W = trim_geometry(EastPlate2, G3_S.p[-1], G3_N.p[0])


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
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'F1', platePatch = True, plateLocation = 'E')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'F2', platePatch = True, plateLocation = 'E')
        # ============== Patch F3 ==============
        F3 = Patch([F3_N, F3_E, F3_S, F3_W], patchName = 'F3', platePatch = True, plateLocation = 'E')

        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patchName = 'G1', platePatch = True, plateLocation = 'W')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patchName = 'G2', platePatch = True, plateLocation = 'W')
        # ============== Patch G3 ==============
        G3 = Patch([G3_N, G3_E, G3_S, G3_W], patchName = 'G3', platePatch = True, plateLocation = 'W')

        # ============== Patch H1 ==============
        H1 = Patch([H1_N, H1_E, H1_S, H1_W], patchName = 'H1')
        # ============== Patch H2 ==============
        H2 = Patch([H2_N, H2_E, H2_S, H2_W], patchName = 'H2')
        # ============== Patch H3 ==============
        H3 = Patch([H3_N, H3_E, H3_S, H3_W], patchName = 'H3')

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
        if tag == 'A3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'A2':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'B3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'B2':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'E2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'E1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'F2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'F1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'G1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'G2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'H1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'H2':
            patch.adjust_corner(xpt1, 'SW')
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'H3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'I3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'I2':
            patch.adjust_corner(xpt2, 'NW')

    def GroupPatches(self):
        # p = self.patches
        # self.PatchGroup = {'SOL' : [], 
        # 'CORE' : (p['ICB'], p['ICT'], p['OCT'], p['OCB']), 
        # 'PF' : (p['IPF'], p['OPF'])}
        pass

    def set_gridue(self):
        """
        set_gridue:
            Prepares 'self.gridue_params' dictionary with required data.
            The self.gridue_params attribute is used to write a gridue
            formatted file
        Parameters:
            N/A
        Return:
            N/A
        """

        ixlb = 0
        ixrb = len(self.rm) - 2

        nxm = len(self.rm) - 2
        nym = len(self.rm[0]) - 2
        iyseparatrix1 = self.patches['A1'].nrad + self.patches['A2'].nrad - 2
        iyseparatrix2 = self.patches['E1'].nrad - 1
        iyseparatrix3 = iyseparatrix2
        iyseparatrix4 = iyseparatrix1

        ix_plate1 = 0
        ix_cut1 = self.patches['A1'].npol - 1

        ix_cut2=0
        for alpha in ['A', 'B', 'C', 'D', 'E']:
            ix_cut2 += self.patches[alpha+'1'].npol - 1

        ix_plate2=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F']:
            ix_plate2 += self.patches[alpha+'3'].npol - 1

        ix_plate3 = ix_plate2 + 2

        ix_cut3=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            ix_cut3 += self.patches[alpha+'1'].npol - 1
        ix_cut3+=2

        ix_cut4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            ix_cut4 += self.patches[alpha+'2'].npol - 1
        ix_cut4+=2

        ix_plate4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
            ix_plate4 += self.patches[alpha+'1'].npol - 1
        ix_plate4+=2

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

        self.gridue_params = {'nxm' : nxm, 'nym' : nym, 'iyseparatrix1' : iyseparatrix1, 'iyseparatrix2' : iyseparatrix2, \
                'ix_plate1' : ix_plate1, 'ix_cut1' : ix_cut1, 'ix_cut2' : ix_cut2, 'ix_plate2' : ix_plate2, 'iyseparatrix3' : iyseparatrix3, \
                'iyseparatrix4' : iyseparatrix4, 'ix_plate3' : ix_plate3, 'ix_cut3' : ix_cut3, 'ix_cut4' : ix_cut4, 'ix_plate4' : ix_plate4, \
                'rm' : self.rm, 'zm' : self.zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b, '_FILLER_' : -1}

