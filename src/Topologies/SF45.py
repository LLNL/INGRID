"""
SF45.py

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


class SF45(TopologyUtils):
    def __init__(self, Ingrid_obj, config):
        TopologyUtils.__init__(self, Ingrid_obj, config)

    def construct_patches(self):

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

        # H1_E / B1_W 
        xpt1N__psiMinCore = self.eq.draw_line(xpt1['N'], {'psi' : psi_min_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        E1_E = xpt1N__psiMinCore
        B1_W = E1_E.reverse_copy()

        # B1_N / B2_S
        xpt1NW__west_midLine = self.eq.draw_line(xpt1['NW'], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_N = xpt1NW__west_midLine
        B2_S = B1_N.reverse_copy()

        # H2_S / H1_N
        xpt1NE__east_midLine = self.eq.draw_line(xpt1['NE'], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E1_N = xpt1NE__east_midLine.reverse_copy()
        E2_S = E1_N.reverse_copy()

        # C1_N / C2_S
        west_midLine__topLine = self.eq.draw_line(xpt1NW__west_midLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C1_N = west_midLine__topLine
        C2_S = C1_N.reverse_copy()

        # D1_N / D2_S
        east_midLine__topLine = self.eq.draw_line(xpt1NE__east_midLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_N = east_midLine__topLine.reverse_copy()
        D2_S = D1_N.reverse_copy()

        # Tracing core: psi-min-core region (rho = 1)
        
        # / B1_S
        psiMinCore__west_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_S = psiMinCore__west_midLine_core.reverse_copy()

        # H1_E1
        psiMinCore__east_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E1_S = psiMinCore__east_midLine_core

        # / C1_S
        west_midLine_core__topLine = self.eq.draw_line(psiMinCore__west_midLine_core.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C1_S = west_midLine_core__topLine.reverse_copy()

        # D1_S
        east_midLine_core__topLine = self.eq.draw_line(psiMinCore__east_midLine_core.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_S = east_midLine_core__topLine

        # A1_E / I1_W
        xpt1__psiMinPF1 = self.eq.draw_line(xpt1['S'], {'psi' : psi_min_pf}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        A1_E = xpt1__psiMinPF1
        F1_W = A1_E.reverse_copy()

        # A1_S
        psiMinPF1__WestPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A1_S = psiMinPF1__WestPlate1

        # / I1_S
        psiMinPF1__EastPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        
        # !!! I1_S = psiMinPF1__EastPlate1.reverse_copy()

        # A2_S / A1_N
        xpt1__WestPlate1 = self.eq.draw_line(xpt1['SW'], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A2_S = xpt1__WestPlate1
        A1_N = A2_S.reverse_copy()

        # I1_N / I2_S
        xpt1__EastPlate1 = self.eq.draw_line(xpt1['SE'], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        # Drawing the portion of the separatrix found in the double-null configuration

        # E2_E / H2_W
        F2_E = self.eq.draw_line(xpt2['N'], {'line' : xpt1__EastPlate1}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        F1_E = self.eq.draw_line(F2_E.p[-1], {'line' : psiMinPF1__EastPlate1}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)

        I1_W = F1_E.reverse_copy()
        I2_W = F2_E.reverse_copy()

        F1_N, I1_N = xpt1__EastPlate1.split(I2_W.p[0], add_split_point=True)
        I2_S = I1_N.reverse_copy()
        F2_S = F1_N.reverse_copy()
        I1_S, F1_S = psiMinPF1__EastPlate1.reverse_copy().split(I1_W.p[0], add_split_point=True)
        
        # H2_N__I2_N
        I2_N = self.eq.draw_line(xpt2['NW'], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        I3_S = I2_N.reverse_copy()

        # E3_S / E2_N
        xpt2NE__east_midLine = self.eq.draw_line(xpt2['NE'], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        F2_W = self.eq.draw_line(xpt1['E'], {'line' : xpt2NE__east_midLine}, \
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        E2_E = F2_W.reverse_copy()

        E2_N, F2_N = xpt2NE__east_midLine.reverse_copy().split(F2_W.p[-1], add_split_point=True)
        F3_S = F2_N.reverse_copy()
        E3_S = E2_N.reverse_copy()

        # D3_S / D2_N
        xpt2__east_midLine__topLine = self.eq.draw_line(xpt2NE__east_midLine.p[-1], {'line' : topLine}, 
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D3_S = xpt2__east_midLine__topLine
        D2_N = D3_S.reverse_copy()

        # C3_S / C2_N
        xpt2__topLine__west_midLine = self.eq.draw_line(xpt2__east_midLine__topLine.p[-1], {'line' : west_midLine},
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        C3_S = xpt2__topLine__west_midLine
        C2_N = C3_S.reverse_copy()

        xpt2__west_midLine__WestPlate1 = self.eq.draw_line(xpt2__topLine__west_midLine.p[-1], {'line' : WestPlate1},
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)


        B2_W = self.eq.draw_line(xpt1['W'], {'line' : xpt2__west_midLine__WestPlate1},
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        A2_E = B2_W.reverse_copy()

        # / A3_S, B3_S
        A2_N, B2_N = xpt2__west_midLine__WestPlate1.reverse_copy().split(B2_W.p[-1], add_split_point=True)
        A3_S = A2_N.reverse_copy()
        B3_S = B2_N.reverse_copy()

        xpt2__psiMinPF2 = self.eq.draw_line(xpt2['S'], {'psi' : psi_min_pf_2}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)

        G1_W, G2_W = xpt2__psiMinPF2.reverse_copy().split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p)//2], add_split_point = True)
        H1_E = G1_W.reverse_copy()
        H2_E = G2_W.reverse_copy()

        G1_S = self.eq.draw_line(G1_W.p[0], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose).reverse_copy()
        H1_S = self.eq.draw_line(G1_W.p[0], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose)
        
        G1_N = self.eq.draw_line(G2_W.p[0], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        G2_S = G1_N.reverse_copy()

        H1_N = self.eq.draw_line(G2_W.p[0], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        H2_S = H1_N.reverse_copy()

        G2_N = self.eq.draw_line(xpt2['SE'], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        G3_S = G2_N.reverse_copy()

        H2_N = self.eq.draw_line(xpt2['SW'], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        H3_S = H2_N.reverse_copy()


        H3_E = self.eq.draw_line(xpt2['W'], {'psi' : psi_max_east}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        I3_W = H3_E.reverse_copy()


        H3_N = self.eq.draw_line(H3_E.p[0], {'line' : EastPlate1}, option = 'theta',
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()

        I3_N = self.eq.draw_line(H3_E.p[0], {'line' : EastPlate1}, option = 'theta',
            direction = 'cw', show_plot = visual, text = verbose)

        G3_W = self.eq.draw_line(xpt2['E'], {'psi' : psi_max_west}, option = 'rho',
            direction = 'ccw', show_plot = visual, text = verbose)
        F3_E = G3_W.reverse_copy()


        G3_N = self.eq.draw_line(G3_W.p[-1], {'line' : EastPlate2}, option = 'theta', 
            direction = 'cw', show_plot = visual, text = verbose)

        G3_W__east_midLine = self.eq.draw_line(G3_W.p[-1], {'line' : east_midLine}, option = 'theta',
            direction = 'ccw', show_plot = visual, text = verbose)

        F3_W = self.eq.draw_line(F2_W.p[-1], {'line' : G3_W__east_midLine}, option = 'rho',
            direction = 'ccw', show_plot = visual, text = verbose)
        E3_E = F3_W.reverse_copy()
        
        E3_N, F3_N = G3_W__east_midLine.reverse_copy().split(F3_W.p[-1], add_split_point=True)

        D3_N = self.eq.draw_line(E3_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', 
            show_plot = visual, text = verbose).reverse_copy()
        C3_N = self.eq.draw_line(D3_N.p[0], {'line' : west_midLine}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()

        west_midLine__WestPlate1 = self.eq.draw_line(C3_N.p[0], {'line' : WestPlate1}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose)
        B3_W = self.eq.draw_line(B2_W.p[-1], {'line' : west_midLine__WestPlate1}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        A3_E = B3_W.reverse_copy()

        A3_N, B3_N = west_midLine__WestPlate1.reverse_copy().split(B3_W.p[-1], add_split_point=True)

        B1_E = self.eq.draw_line(B1_N.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C1_W = B1_E.reverse_copy()

        B2_E = self.eq.draw_line(B2_N.p[-1], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C2_W = B2_E.reverse_copy()

        B3_E = self.eq.draw_line(B3_N.p[-1], {'psi_horizontal' : Point(xpt2['center']).psi(self)}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C3_W = B3_E.reverse_copy()

        C1_E = self.eq.draw_line(C1_N.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D1_W = C1_E.reverse_copy()

        C2_E = self.eq.draw_line(C2_N.p[-1], {'psi_vertical' : 1.0}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D2_W = C2_E.reverse_copy()

        C3_E = self.eq.draw_line(C3_N.p[-1], {'psi_vertical' : Point(xpt2['center']).psi(self)}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D3_W = C3_E.reverse_copy()

        E1_W = self.eq.draw_line(E1_N.p[0], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D1_E = E1_W.reverse_copy()

        E2_W = self.eq.draw_line(E2_N.p[0], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D2_E = E2_W.reverse_copy()

        E3_W = self.eq.draw_line(E3_N.p[0], {'psi_horizontal' : Point(xpt2['center']).psi(self)}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D3_E = E3_W.reverse_copy()

        A1_W = trim_geometry(WestPlate1, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate1, A2_S.p[-1], A2_N.p[0])
        A3_W = trim_geometry(WestPlate1, A3_S.p[-1], A3_N.p[0])

        H1_W = trim_geometry(WestPlate2, H1_S.p[-1], H1_N.p[0])
        H2_W = trim_geometry(WestPlate2, H2_S.p[-1], H2_N.p[0])
        H3_W = trim_geometry(WestPlate2, H3_S.p[-1], H3_N.p[0])

        I1_E = trim_geometry(EastPlate1, I1_N.p[-1], I1_S.p[0])
        I2_E = trim_geometry(EastPlate1, I2_N.p[-1], I2_S.p[0])
        I3_E = trim_geometry(EastPlate1, I3_N.p[-1], I3_S.p[0])

        G1_E = trim_geometry(EastPlate2, G1_N.p[-1], G1_S.p[0])
        G2_E = trim_geometry(EastPlate2, G2_N.p[-1], G2_S.p[0])
        G3_E = trim_geometry(EastPlate2, G3_N.p[-1], G3_S.p[0])

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
        if tag == 'A2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'A1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'B2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'B1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'E2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'E1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'F3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'F2':
            patch.adjust_corner(xpt1, 'SW')
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'F1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'G3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'G2':
            patch.adjust_corner(xpt2, 'NW')
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
        iyseparatrix2 = iyseparatrix1
        iyseparatrix3 = self.patches['A1'].nrad + self.patches['A2'].nrad - 2
        iyseparatrix4 = iyseparatrix3

        ix_plate1 = 0
        ix_cut1 = self.patches['A1'].npol - 1

        ix_cut2=0
        for alpha in ['A', 'B', 'C', 'D', 'E']:
            ix_cut2 += self.patches[alpha+'1'].npol - 1

        ix_plate2=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            ix_plate2 += self.patches[alpha+'3'].npol - 1

        ix_plate3 = ix_plate2 + 2

        ix_cut3=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F']:
            ix_cut3 += self.patches[alpha+'2'].npol - 1

        ix_cut4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            ix_cut4 += self.patches[alpha+'1'].npol - 1
        ix_cut4 +=2


        ix_plate4=0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
            ix_plate4 += self.patches[alpha+'1'].npol - 1
        ix_plate4 +=2

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

