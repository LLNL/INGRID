"""
UDN.py

Description:
    UDN (unbalanced double-null) configuration class.

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


class UDN(TopologyUtils):
    def __init__(self, Ingrid_obj, config):
        TopologyUtils.__init__(self, Ingrid_obj, config)

    def construct_patches(self):
        """
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
        xpt2_psi = self.psi_norm.get_psi(xpt2['center'][0], xpt2['center'][1])

        magx = np.array([self.settings['grid_params']['rmagx'] + self.settings['grid_params']['patch_generation']['rmagx_shift'], \
            self.settings['grid_params']['zmagx'] + self.settings['grid_params']['patch_generation']['zmagx_shift']])

        psi_max_west = self.settings['grid_params']['psi_max_west']
        psi_max_east = self.settings['grid_params']['psi_max_east']
        psi_min_core = self.settings['grid_params']['psi_min_core']
        psi_min_pf = self.settings['grid_params']['psi_min_pf']
        psi_min_pf_2 = self.settings['grid_params']['psi_pf2']

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
        # topLine.plot()


        # Tracing primary-separatrix: core-boundary

        xpt1N__psiMinCore = self.eq.draw_line(xpt1['N'], {'psi' : psi_min_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        B1_W = xpt1N__psiMinCore
        G1_E = B1_W.reverse_copy()


        xpt1NW__west_midLine = self.eq.draw_line(xpt1['NW'], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_N = xpt1NW__west_midLine
        B2_S = xpt1NW__west_midLine.reverse_copy()


        xpt1NE__east_midLine = self.eq.draw_line(xpt1['NE'], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        G1_N = xpt1NE__east_midLine.reverse_copy()
        G2_S = xpt1NE__east_midLine


        # Tracing core: psi-min-core region (rho = 1)
        

        psiMinCore__west_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_S = psiMinCore__west_midLine_core.reverse_copy()


        psiMinCore__east_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        G1_S = psiMinCore__east_midLine_core


        xpt1__psiMinPF1 = self.eq.draw_line(xpt1['S'], {'psi' : psi_min_pf}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        A1_E = xpt1__psiMinPF1
        H1_W = A1_E.reverse_copy()


        psiMinPF1__WestPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A1_S = psiMinPF1__WestPlate1


        psiMinPF1__EastPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        H1_S = psiMinPF1__EastPlate1.reverse_copy()


        xpt1__WestPlate1 = self.eq.draw_line(xpt1['SW'], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A2_S = xpt1__WestPlate1
        A1_N = A2_S.reverse_copy()

        # I1_N / I2_S
        xpt1__EastPlate1 = self.eq.draw_line(xpt1['SE'], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        H1_N = xpt1__EastPlate1
        H2_S = H1_N.reverse_copy()

        # Drawing the portion of the separatrix found in the double-null configuration

        xpt2N__outer_core = self.eq.draw_line(xpt2['N'], {'psi' : 1.0}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        C2_E = xpt2N__outer_core
        F2_W = C2_E.reverse_copy()


        xpt2N__outer_core__inner_core = self.eq.draw_line(xpt2N__outer_core.p[-1], \
            {'psi' : psi_min_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        C1_E = xpt2N__outer_core__inner_core
        F1_W = C1_E.reverse_copy()


        C1_S = self.eq.draw_line(C1_E.p[-1], {'line' : west_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

        F1_S = self.eq.draw_line(C1_E.p[-1], {'line' : west_midLine}, option = 'theta',
            direction = 'cw', show_plot = visual, text = verbose).reverse_copy()


        C1_N = self.eq.draw_line(xpt2N__outer_core.p[-1], {'line' : west_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        C2_S = C1_N.reverse_copy()

        F1_N = self.eq.draw_line(xpt2N__outer_core.p[-1], {'line' : east_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        F2_S = F1_N.reverse_copy()


        xpt2NW__east_midLine = self.eq.draw_line(xpt2['NW'], {'line' : east_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        F2_N = xpt2NW__east_midLine
        F3_S = F2_N.reverse_copy()


        xpt2NE__west_midLine = self.eq.draw_line(xpt2['NE'], {'line' : west_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        C2_N = xpt2NE__west_midLine.reverse_copy()
        C3_S = C2_N.reverse_copy()


        xpt2__east_midLine__EastPlate1 = self.eq.draw_line(xpt2NW__east_midLine.p[-1], {'line' : EastPlate1}, 
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        xpt2__west_midLine__WestPlate1 = self.eq.draw_line(xpt2NE__west_midLine.p[-1], {'line' : WestPlate1}, 
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)


        H2_W = self.eq.draw_line(xpt1['E'], {'line' : xpt2__east_midLine__EastPlate1}, option = 'rho', direction = 'ccw', 
            show_plot = visual, text = verbose)
        G2_E = H2_W.reverse_copy()


        B2_W = self.eq.draw_line(xpt1['W'], {'line' : xpt2__west_midLine__WestPlate1}, option = 'rho', direction = 'ccw', 
            show_plot = visual, text = verbose)
        A2_E = B2_W.reverse_copy()


        G2_N, H2_N = xpt2__east_midLine__EastPlate1.split(H2_W.p[-1], add_split_point=True)
        G3_S, H3_S = G2_N.reverse_copy(), H2_N.reverse_copy()
        

        A2_N, B2_N = xpt2__west_midLine__WestPlate1.reverse_copy().split(B2_W.p[-1], add_split_point=True)
        A3_S, B3_S = A2_N.reverse_copy(), B2_N.reverse_copy()


        xpt2__psiMinPF2 = self.eq.draw_line(xpt2['S'], {'psi' : psi_min_pf_2}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        E2_E, E1_E = xpt2__psiMinPF2.split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p)//2], add_split_point = True)
        D1_W = E1_E.reverse_copy()
        D2_W = E2_E.reverse_copy()


        E1_S = self.eq.draw_line(E1_E.p[-1], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose)

        D1_S = self.eq.draw_line(E1_E.p[-1], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose).reverse_copy()

        
        E1_N = self.eq.draw_line(E2_E.p[-1], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        E2_S = E1_N.reverse_copy()


        D1_N = self.eq.draw_line(E2_E.p[-1], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        D2_S = D1_N.reverse_copy()


        E2_N = self.eq.draw_line(xpt2['SW'], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        E3_S = E2_N.reverse_copy()


        D2_N = self.eq.draw_line(xpt2['SE'], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        D3_S = D2_N.reverse_copy()


        E3_E = self.eq.draw_line(xpt2['W'], {'psi' : psi_max_east}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        F3_W = E3_E.reverse_copy()


        D3_W = self.eq.draw_line(xpt2['E'], {'psi' : psi_max_west}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        C3_E = D3_W.reverse_copy()


        F3_N = self.eq.draw_line(F3_W.p[-1], {'line' : east_midLine}, option = 'theta',
            direction = 'cw', show_plot = visual, text = verbose)


        C3_N = self.eq.draw_line(D3_W.p[-1], {'line' : west_midLine}, option = 'theta',
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()


        D3_N = self.eq.draw_line(D3_W.p[-1], {'line' : EastPlate2}, option = 'theta', direction = 'cw', 
            show_plot=visual, text=verbose)

        E3_N = self.eq.draw_line(E3_E.p[0], {'line' : WestPlate2}, option = 'theta', direction = 'ccw', 
            show_plot=visual, text=verbose).reverse_copy()


        east_midLine__EastPlate1 = self.eq.draw_line(F3_N.p[-1], {'line' : EastPlate1}, 
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        west_midLine__WestPlate1 = self.eq.draw_line(C3_N.p[0], {'line' : WestPlate1},
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)


        H3_W = self.eq.draw_line(H2_W.p[-1], {'line' : east_midLine__EastPlate1},
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        G3_E = H3_W.reverse_copy()

        B3_W = self.eq.draw_line(B2_W.p[-1], {'line' : west_midLine__WestPlate1},
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        A3_E = B3_W.reverse_copy()

        A3_N, B3_N = west_midLine__WestPlate1.reverse_copy().split(B3_W.p[-1], add_split_point=True)
        G3_N, H3_N = east_midLine__EastPlate1.split(H3_W.p[-1], add_split_point=True)


        B1_E = self.eq.draw_line(B1_N.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C1_W = B1_E.reverse_copy()

        B2_E = self.eq.draw_line(B2_N.p[-1], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C2_W = B2_E.reverse_copy()

        B3_E = self.eq.draw_line(B3_N.p[-1], {'psi_horizontal' : xpt2_psi}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C3_W = B3_E.reverse_copy()


        F1_E = self.eq.draw_line(F1_N.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        G1_W = F1_E.reverse_copy()

        F2_E = self.eq.draw_line(F2_N.p[-1], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        G2_W = F2_E.reverse_copy()

        F3_E = self.eq.draw_line(F3_N.p[-1], {'psi_horizontal' : xpt2_psi}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        G3_W = F3_E.reverse_copy()



        A1_W = trim_geometry(WestPlate1, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate1, A2_S.p[-1], A2_N.p[0])
        A3_W = trim_geometry(WestPlate1, A3_S.p[-1], A3_N.p[0])

        D1_E = trim_geometry(WestPlate2, D1_N.p[-1], D1_S.p[0])
        D2_E = trim_geometry(WestPlate2, D2_N.p[-1], D2_S.p[0])
        D3_E = trim_geometry(WestPlate2, D3_N.p[-1], D3_S.p[0])

        E1_W = trim_geometry(EastPlate2, E1_S.p[-1], E1_N.p[0])
        E2_W = trim_geometry(EastPlate2, E2_S.p[-1], E2_N.p[0])
        E3_W = trim_geometry(EastPlate2, E3_S.p[-1], E3_N.p[0])

        H1_E = trim_geometry(EastPlate1, H1_N.p[-1], H1_S.p[0])
        H2_E = trim_geometry(EastPlate1, H2_N.p[-1], H2_S.p[0])
        H3_E = trim_geometry(EastPlate1, H3_N.p[-1], H3_S.p[0])

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
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'D1', platePatch = True, plateLocation = 'E')
        # ============== Patch D2 ==============
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'D2', platePatch = True, plateLocation = 'E')
        # ============== Patch D3 ==============
        D3 = Patch([D3_N, D3_E, D3_S, D3_W], patchName = 'D3', platePatch = True, plateLocation = 'E')

        # ============== Patch E1 ==============
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'E1', platePatch = True, plateLocation = 'W')
        # ============== Patch E2 ==============
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'E2', platePatch = True, plateLocation = 'W')
        # ============== Patch E3 ==============
        E3 = Patch([E3_N, E3_E, E3_S, E3_W], patchName = 'E3', platePatch = True, plateLocation = 'W')

        # ============== Patch F1 ==============
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'F1')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'F2')
        # ============== Patch F3 ==============
        F3 = Patch([F3_N, F3_E, F3_S, F3_W], patchName = 'F3')

        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patchName = 'G1')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patchName = 'G2')
        # ============== Patch G3 ==============
        G3 = Patch([G3_N, G3_E, G3_S, G3_W], patchName = 'G3')

        # ============== Patch H1 ==============
        H1 = Patch([H1_N, H1_E, H1_S, H1_W], patchName = 'H1', platePatch = True, plateLocation = 'E')
        # ============== Patch H2 ==============
        H2 = Patch([H2_N, H2_E, H2_S, H2_W], patchName = 'H2', platePatch = True, plateLocation = 'E')
        # ============== Patch H3 ==============
        H3 = Patch([H3_N, H3_E, H3_S, H3_W], patchName = 'H3', platePatch = True, plateLocation = 'E')

        patches = [A3, A2, A1, B3, B2, B1, C3, C2, C1, D3, D2, D1, E3, E2, E1,
                   F3, F2, F1, G3, G2, G1, H3, H2, H1]

        self.patches = {}
        for patch in patches:
            patch.parent = self
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch

    def AdjustPatch(self,patch):
        xpt1 = Point(self.eq.NSEW_lookup['xpt1']['coor']['center'])
        xpt2 = Point(self.eq.NSEW_lookup['xpt2']['coor']['center'])

        tag = patch.get_tag()
        if tag == 'A2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'B2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'B1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'A1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'H1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'H2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'G2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'G1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'C3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'C2':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'D2':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'D3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'E2':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'E3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'F3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'F2':
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
        iyseparatrix1 = self.patches['A1'].nrad - 1
        iyseparatrix2 = self.patches['C1'].nrad + self.patches['C2'].nrad - 2
        iyseparatrix3 = iyseparatrix2
        iyseparatrix4 = iyseparatrix1
        
        ix_plate1 = 0
        ix_cut1 = self.patches['A1'].npol - 1

        ix_cut2=0
        for alpha in ['A', 'B', 'C']:
            ix_cut2 += self.patches[alpha+'1'].npol - 1

        ix_plate2=0
        for alpha in ['A', 'B', 'C', 'D']:
            ix_plate2 += self.patches[alpha+'3'].npol - 1


        ix_plate3 = ix_plate2 + 2

        ix_cut3=0
        for alpha in ['A', 'B', 'C', 'D', 'E']:
            ix_cut3 += self.patches[alpha+'2'].npol - 1
        ix_cut3+=2

        ix_cut4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            ix_cut4 += self.patches[alpha+'1'].npol - 1
        ix_cut4+=2

        ix_plate4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
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
