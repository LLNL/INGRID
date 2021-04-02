"""
The ``sf15h`` module contains :class:`SF15H` for representing a Snowflake-15
half domain configuration.

Child of base :class:`utils.TopologyUtils`.

"""
import numpy as np
import matplotlib
import pathlib
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt
from INGRID.utils import TopologyUtils
from INGRID.geometry import Point, Line, Patch, trim_geometry
from collections import OrderedDict


class SF15H(TopologyUtils):
    """
    The `SF15H` class for handling `Snowflake-15` half-domain configurations within a tokamak.

    Parameters
    ----------
    Ingrid_obj : Ingrid
        Ingrid object the SF15H object is being managed by.

    config : str, optional
        String code representing the configuration.

    Attributes
    ----------
    ConnexionMap : dict
        A mapping defining dependencies between Patch objects for grid generation.

    patches : dict
        The collection of Patch objects representing the topology.

    Attributes
    ----------
    ConnexionMap : dict
        A mapping defining dependencies between Patch objects for grid generation.

    patches : dict
        The collection of Patch objects representing the topology.
    """
    def __init__(self, Ingrid_obj: 'ingrid.Ingrid', config: str = 'SF15H'):
        TopologyUtils.__init__(self, Ingrid_obj, config)

        self.ConnexionMap = {
            'A1': {'N': ('A2', 'S')},
            'A2': {'N': ('A3', 'S')},
            'A3': None,

            'B1': {'N': ('B2', 'S')},
            'B2': {'N': ('B3', 'S'), 'W': ('A2', 'E')},
            'B3': {'W': ('A3', 'E')},

            'C1': {'N': ('C2', 'S')},
            'C2': {'N': ('C3', 'S')},
            'C3': None,

            'D1': {'N': ('D2', 'S'), 'W': ('E1', 'E')},
            'D2': {'N': ('D3', 'S'), 'W': ('E2', 'E')},
            'D3': {'W': ('C3', 'E')},

            'E1': {'N': ('E2', 'S')},
            'E2': {'N': ('E3', 'S')},
            'E3': None,

            'F1': {'N': ('F2', 'S'), 'W': ('C1', 'E')},
            'F2': {'N': ('F3', 'S'), 'W': ('C2', 'E')},
            'F3': {'W': ('E3', 'E')},

            'G1': {'N': ('G2', 'S'), 'W': ('F2', 'E')},
            'G2': {'N': ('G3', 'S'), 'W': ('F3', 'E')},
            'G3': {'W': ('F3', 'E')},

        }
        self.gridue_settings = None

    def AdjustGrid(self) -> None:
        """
        Adjust the grid so that no holes occur at x-points, and cell grid
        faces are alligned

        A small epsilon radius is swept out around x-points during Patch
        line tracing. This simple tidies up a grid.

        Parameters
        ----------
        patch : Patch
            The patch to tidy up (will only adjust if next to x-point).
        """

        for patch in self.patches.values():
            # Adjust cell to any adjacent x-point
            self.AdjustPatch(patch)

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
            magx_tilt_1 = self.settings['grid_settings']['patch_generation']['magx_tilt_1']
        except KeyError:
            magx_tilt_1 = 0.0
        try:
            magx_tilt_2 = self.settings['grid_settings']['patch_generation']['magx_tilt_2']
        except KeyError:
            magx_tilt_2 = 0.0

        xpt1 = self.LineTracer.NSEW_lookup['xpt1']['coor']
        xpt2 = self.LineTracer.NSEW_lookup['xpt2']['coor']

        magx = np.array([self.settings['grid_settings']['rmagx'] + self.settings['grid_settings']['patch_generation']['rmagx_shift'],
            self.settings['grid_settings']['zmagx'] + self.settings['grid_settings']['patch_generation']['zmagx_shift']])

        xpt2_psi = self.PsiNorm.get_psi(xpt2['center'][0], xpt2['center'][1])
        psi_1 = self.settings['grid_settings']['psi_1']
        psi_2 = self.settings['grid_settings']['psi_2']
        psi_core = self.settings['grid_settings']['psi_core']
        psi_pf_1 = self.settings['grid_settings']['psi_pf_1']
        psi_pf_2 = self.settings['grid_settings']['psi_pf_2']

        if self.settings['grid_settings']['patch_generation']['strike_pt_loc'] == 'limiter':
            WestPlate1 = self.parent.LimiterData.copy()
            WestPlate2 = self.parent.LimiterData.copy()

            EastPlate1 = self.parent.LimiterData.copy()
            EastPlate2 = self.parent.LimiterData.copy()

        else:
            WestPlate1 = self.PlateData['plate_W1']
            WestPlate2 = self.PlateData['plate_W2']

            EastPlate1 = self.PlateData['plate_E1']
            EastPlate2 = self.PlateData['plate_E2']

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(magx_tilt_1), magx[1] - 1e6 * np.sin(magx_tilt_1))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(magx_tilt_1), magx[1] + 1e6 * np.sin(magx_tilt_1))
        midline_1 = Line([LHS_Point, RHS_Point])

        LHS_Point = Point(magx[0] - 1e6 * np.cos(magx_tilt_2), magx[1] - 1e6 * np.sin(magx_tilt_2))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(magx_tilt_2), magx[1] + 1e6 * np.sin(magx_tilt_2))
        midline_2 = Line([LHS_Point, RHS_Point])

        # Tracing primary-separatrix: core-boundary

        xpt1N__psiMinCore = self.LineTracer.draw_line(xpt1['N'], {'psi': psi_core},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        F1_E = xpt1N__psiMinCore
        B1_W = F1_E.reverse_copy()

        # B1_N / B2_S
        xpt1NW__midline_1 = self.LineTracer.draw_line(xpt1['NW'], {'line': midline_1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        B1_N = xpt1NW__midline_1
        B2_S = B1_N.reverse_copy()

        # Split into C1_N and F1_N
        xpt1NE__midline_2 = self.LineTracer.draw_line(xpt1['NE'], {'line': midline_2},
            option='theta', direction='ccw', show_plot=visual, text=verbose)

        # Tracing core: psi-min-core region (rho = 1)

        # / B1_S
        psiMinCore__midline_1_core = self.LineTracer.draw_line(xpt1N__psiMinCore.p[-1], {'line': midline_1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        B1_S = psiMinCore__midline_1_core.reverse_copy()

        # Split into C1_S and F1_S
        psiMinCore__midline_2_core = self.LineTracer.draw_line(xpt1N__psiMinCore.p[-1], {'line': midline_2},
            option='theta', direction='ccw', show_plot=visual, text=verbose)

        # A1_E / G1_W
        xpt1__psiMinPF1 = self.LineTracer.draw_line(xpt1['S'], {'psi': psi_pf_1},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        A1_E = xpt1__psiMinPF1
        G1_W = A1_E.reverse_copy()

        # A1_S
        psiMinPF1__WestPlate1 = self.LineTracer.draw_line(xpt1__psiMinPF1.p[-1], {'line': WestPlate1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        A1_S = psiMinPF1__WestPlate1

        # / G1_S
        psiMinPF1__EastPlate1 = self.LineTracer.draw_line(xpt1__psiMinPF1.p[-1], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        G1_S = psiMinPF1__EastPlate1.reverse_copy()

        # A2_S / A1_N
        xpt1__WestPlate1 = self.LineTracer.draw_line(xpt1['SW'], {'line': WestPlate1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        A2_S = xpt1__WestPlate1
        A1_N = A2_S.reverse_copy()

        # I1_N / I2_S
        xpt1__EastPlate1 = self.LineTracer.draw_line(xpt1['SE'], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        G1_N = xpt1__EastPlate1
        G2_S = G1_N.reverse_copy()

        # Drawing the portion of the separatrix found in the double-null configuration

        xpt2N__outer_core = self.LineTracer.draw_line(xpt2['N'], {'line': xpt1NE__midline_2},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        C2_E = xpt2N__outer_core
        F2_W = C2_E.reverse_copy()

        xpt2N__outer_core__inner_core = self.LineTracer.draw_line(xpt2N__outer_core.p[-1],
            {'line': psiMinCore__midline_2_core},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        C1_E = xpt2N__outer_core__inner_core
        F1_W = C1_E.reverse_copy()

        C1_N, F1_N = xpt1NE__midline_2.reverse_copy().split(xpt2N__outer_core.p[-1], add_split_point=True)
        C2_S = C1_N.reverse_copy()
        F2_S = F1_N.reverse_copy()

        F1_S, C1_S = psiMinCore__midline_2_core.split(xpt2N__outer_core__inner_core.p[-1], add_split_point=True)

        xpt2NW__EastPlate1 = self.LineTracer.draw_line(xpt2['NW'], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)

        xpt2NE__midline_2 = self.LineTracer.draw_line(xpt2['NE'], {'line': midline_2},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        C3_S = xpt2NE__midline_2
        C2_N = C3_S.reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt1_W']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_W_tilt']
            B2_W = self.LineTracer.draw_line(xpt1['W'], {'psi_horizontal':(xpt2_psi, tilt)}, option='z_const', direction='ccw', show_plot=visual, text=verbose)
        else:
            B2_W = self.LineTracer.draw_line(xpt1['W'], {'psi': xpt2_psi}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        A2_E = B2_W.reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt1_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_E_tilt']
            G2_W = self.LineTracer.draw_line(xpt1['E'], {'line': (xpt2NW__EastPlate1, tilt)}, option='z_const', direction='cw', show_plot=visual, text=verbose)
        else:
            G2_W = self.LineTracer.draw_line(xpt1['E'], {'line': xpt2NW__EastPlate1}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        F2_E = G2_W.reverse_copy()

        B2_N = self.LineTracer.draw_line(B2_W.p[-1], {'line': midline_1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        A2_N = self.LineTracer.draw_line(B2_W.p[-1], {'line': WestPlate1}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        A3_S = A2_N.reverse_copy()
        B3_S = B2_N.reverse_copy()

        F2_N, G2_N = xpt2NW__EastPlate1.split(G2_W.p[-1], add_split_point=True)
        F3_S = F2_N.reverse_copy()
        G3_S = G2_N.reverse_copy()

        xpt2__psiMinPF2 = self.LineTracer.draw_line(xpt2['S'], {'psi': psi_pf_2},
            option='rho', direction='cw', show_plot=visual, text=verbose)

        D1_W, D2_W = xpt2__psiMinPF2.reverse_copy().split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p) // 2], add_split_point=True)
        E1_E = D1_W.reverse_copy()
        E2_E = D2_W.reverse_copy()

        D1_S = self.LineTracer.draw_line(D1_W.p[0], {'line': EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()
        E1_S = self.LineTracer.draw_line(D1_W.p[0], {'line': WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)

        D1_N = self.LineTracer.draw_line(D2_W.p[0], {'line': EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        D2_S = D1_N.reverse_copy()

        E1_N = self.LineTracer.draw_line(D2_W.p[0], {'line': WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        E2_S = E1_N.reverse_copy()

        D2_N = self.LineTracer.draw_line(xpt2['SE'], {'line': EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        D3_S = D2_N.reverse_copy()

        E2_N = self.LineTracer.draw_line(xpt2['SW'], {'line': WestPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        E3_S = E2_N.reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt2_W']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt2_W_tilt']
            F3_W = self.LineTracer.draw_line(xpt2['W'],
                {'psi_horizontal': (psi_2, tilt)},
                option='z_const', direction='ccw',
                show_plot=visual, text=verbose)
        else:
            F3_W = self.LineTracer.draw_line(xpt2['W'], {'psi': psi_2}, option='rho', direction='ccw',
                show_plot=visual, text=verbose)
        E3_E = F3_W.reverse_copy()

        F3_W__EastPlate1 = self.LineTracer.draw_line(F3_W.p[-1], {'line': EastPlate1}, option='theta',
            direction='cw', show_plot=visual, text=verbose)

        if self.settings['grid_settings']['patch_generation']['use_xpt1_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_E_tilt']
            G3_W = self.LineTracer.draw_line(G2_W.p[-1], {'line': (F3_W__EastPlate1, tilt)},
                option='z_const', direction='cw', show_plot=visual, text=verbose)
        else:
            G3_W = self.LineTracer.draw_line(G2_W.p[-1], {'line': F3_W__EastPlate1},
                option='rho', direction='ccw', show_plot=visual, text=verbose)

        F3_E = G3_W.reverse_copy()

        F3_N, G3_N = F3_W__EastPlate1.split(G3_W.p[-1], add_split_point=True)

        E3_N = self.LineTracer.draw_line(F3_W.p[-1], {'line': WestPlate2}, option='theta',
            direction='ccw', show_plot=visual, text=verbose).reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt2_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt2_E_tilt']
            D3_W = self.LineTracer.draw_line(xpt2['E'],
                {'psi_horizontal': (psi_1, tilt)},
                option='z_const', direction='cw',
                show_plot=visual, text=verbose)
        else:
            D3_W = self.LineTracer.draw_line(xpt2['E'], {'psi': psi_1}, option='rho', direction='ccw',
                show_plot=visual, text=verbose)

        C3_E = D3_W.reverse_copy()

        D3_N = self.LineTracer.draw_line(D3_W.p[-1], {'line': EastPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        C3_N = self.LineTracer.draw_line(D3_W.p[-1], {'line': midline_2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt1_W']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_W_tilt']
            B3_W = self.LineTracer.draw_line(B2_W.p[-1], {'psi_horizontal': (psi_1, tilt)}, option='z_const', direction='ccw',
                show_plot=visual, text=verbose)
        else:
            B3_W = self.LineTracer.draw_line(B2_W.p[-1], {'psi': psi_1}, option='rho', direction='ccw',
                show_plot=visual, text=verbose)

        A3_E = B3_W.reverse_copy()

        B3_N = self.LineTracer.draw_line(B3_W.p[-1], {'line': midline_1}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        A3_N = self.LineTracer.draw_line(B3_W.p[-1], {'line': WestPlate1}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()

        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B3_E = Line([B3_N.p[-1], B3_S.p[0]])

        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C3_W = Line([C3_S.p[-1], C3_N.p[0]])

        A1_W = trim_geometry(WestPlate1, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate1, A2_S.p[-1], A2_N.p[0])
        A3_W = trim_geometry(WestPlate1, A3_S.p[-1], A3_N.p[0])

        E1_W = trim_geometry(WestPlate2, E1_S.p[-1], E1_N.p[0])
        E2_W = trim_geometry(WestPlate2, E2_S.p[-1], E2_N.p[0])
        E3_W = trim_geometry(WestPlate2, E3_S.p[-1], E3_N.p[0])

        G1_E = trim_geometry(EastPlate1, G1_N.p[-1], G1_S.p[0])
        G2_E = trim_geometry(EastPlate1, G2_N.p[-1], G2_S.p[0])
        G3_E = trim_geometry(EastPlate1, G3_N.p[-1], G3_S.p[0])

        D1_E = trim_geometry(EastPlate2, D1_N.p[-1], D1_S.p[0])
        D2_E = trim_geometry(EastPlate2, D2_N.p[-1], D2_S.p[0])
        D3_E = trim_geometry(EastPlate2, D3_N.p[-1], D3_S.p[0])

        # ============== Patch A1 ==============
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patch_name='A1', plate_patch=True, plate_location='W')
        # ============== Patch A2 ==============
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patch_name='A2', plate_patch=True, plate_location='W')
        # ============== Patch A3 ==============
        A3 = Patch([A3_N, A3_E, A3_S, A3_W], patch_name='A3', plate_patch=True, plate_location='W')

        # ============== Patch B1 ==============
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patch_name='B1')
        # ============== Patch B2 ==============
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patch_name='B2')
        # ============== Patch B3 ==============
        B3 = Patch([B3_N, B3_E, B3_S, B3_W], patch_name='B3')

        # ============== Patch C1 ==============
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patch_name='C1')
        # ============== Patch C2 ==============
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patch_name='C2')
        # ============== Patch C3 ==============
        C3 = Patch([C3_N, C3_E, C3_S, C3_W], patch_name='C3')

        # ============== Patch D1 ==============
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patch_name='D1', plate_patch=True, plate_location='E')
        # ============== Patch D2 ==============
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patch_name='D2', plate_patch=True, plate_location='E')
        # ============== Patch D3 ==============
        D3 = Patch([D3_N, D3_E, D3_S, D3_W], patch_name='D3', plate_patch=True, plate_location='E')

        # ============== Patch E1 ==============
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patch_name='E1', plate_patch=True, plate_location='W')
        # ============== Patch E2 ==============
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patch_name='E2', plate_patch=True, plate_location='W')
        # ============== Patch E3 ==============
        E3 = Patch([E3_N, E3_E, E3_S, E3_W], patch_name='E3', plate_patch=True, plate_location='W')

        # ============== Patch F1 ==============
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patch_name='F1')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patch_name='F2')
        # ============== Patch F3 ==============
        F3 = Patch([F3_N, F3_E, F3_S, F3_W], patch_name='F3')

        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patch_name='G1', plate_patch=True, plate_location='E')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patch_name='G2', plate_patch=True, plate_location='E')
        # ============== Patch G3 ==============
        G3 = Patch([G3_N, G3_E, G3_S, G3_W], patch_name='G3', plate_patch=True, plate_location='E')

        patches = [A3, B3, C3, D3, E3, F3, G3, A2, B2, C2, F2, G2,
                   E2, D2, A1, G1, E1, D1, C1, F1, B1]

        self.patches = {}
        for patch in patches:
            patch.parent = self
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patch_name] = patch
        self.OrderPatches()

    def OrderPatches(self):
        pass

    def AdjustPatch(self, patch):
        xpt1 = Point(self.LineTracer.NSEW_lookup['xpt1']['coor']['center'])
        xpt2 = Point(self.LineTracer.NSEW_lookup['xpt2']['coor']['center'])

        tag = patch.get_tag()
        if tag == 'A2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'A1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'B2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'B1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'F2':
            patch.adjust_corner(xpt1, 'SE')
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'F1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'G2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'G1':
            patch.adjust_corner(xpt1, 'NW')
        elif tag == 'F3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'E3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'E2':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'D2':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'D3':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'C3':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'C2':
            patch.adjust_corner(xpt2, 'NE')

    def GroupPatches(self):
        pass

    def set_gridue(self):
        """
        Prepares a ``gridue_settings`` dictionary with required data
        for writing a gridue file.

        Parameters
        ----------

        Returns
        -------

        """

        ixlb = 0
        ixrb = len(self.rm) - 2

        nxm = len(self.rm) - 2
        nym = len(self.rm[0]) - 2
        iyseparatrix1 = self.patches['A1'].nrad - 1
        iyseparatrix2 = self.patches['E1'].nrad + self.patches['E2'].nrad - 2
        iyseparatrix3 = iyseparatrix2
        iyseparatrix4 = iyseparatrix1

        ix_plate1 = 0
        ix_cut1 = self.patches['A1'].npol - 1

        ix_cut2 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E']:
            ix_cut2 += self.patches[alpha + '1'].npol - 1

        ix_plate2 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F']:
            ix_plate2 += self.patches[alpha + '3'].npol - 1

        ix_plate3 = ix_plate2 + 2

        ix_cut3 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            ix_cut3 += self.patches[alpha + '2'].npol - 1
        ix_cut3 += 2

        ix_cut4 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            ix_cut4 += self.patches[alpha + '1'].npol - 1
        ix_cut4 += 2

        ix_plate4 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
            ix_plate4 += self.patches[alpha + '1'].npol - 1
        ix_plate4 += 2

        psi = np.zeros((nxm + 2, nym + 2, 5), order='F')
        br = np.zeros((nxm + 2, nym + 2, 5), order='F')
        bz = np.zeros((nxm + 2, nym + 2, 5), order='F')
        bpol = np.zeros((nxm + 2, nym + 2, 5), order='F')
        bphi = np.zeros((nxm + 2, nym + 2, 5), order='F')
        b = np.zeros((nxm + 2, nym + 2, 5), order='F')

        rm = self.rm
        zm = self.zm
        rb_prod = self.PsiUNorm.rcenter * self.PsiUNorm.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.PsiUNorm.get_psi(_r, _z)
                    _br = self.PsiUNorm.get_psi(_r, _z, tag='vz') / _r
                    _bz = -self.PsiUNorm.get_psi(_r, _z, tag='vr') / _r
                    _bpol = np.sqrt(_br ** 2 + _bz ** 2)
                    _bphi = rb_prod / _r
                    _b = np.sqrt(_bpol ** 2 + _bphi ** 2)

                    psi[i][j][k] = _psi
                    br[i][j][k] = _br
                    bz[i][j][k] = _bz
                    bpol[i][j][k] = _bpol
                    bphi[i][j][k] = _bphi
                    b[i][j][k] = _b

        self.gridue_settings = {
            'nxm': nxm, 'nym': nym, 'iyseparatrix1': iyseparatrix1, 'iyseparatrix2': iyseparatrix2,
            'ix_plate1': ix_plate1, 'ix_cut1': ix_cut1, 'ix_cut2': ix_cut2, 'ix_plate2': ix_plate2, 'iyseparatrix3': iyseparatrix3,
            'iyseparatrix4': iyseparatrix4, 'ix_plate3': ix_plate3, 'ix_cut3': ix_cut3, 'ix_cut4': ix_cut4, 'ix_plate4': ix_plate4,
            'rm': self.rm, 'zm': self.zm, 'psi': psi, 'br': br, 'bz': bz, 'bpol': bpol, 'bphi': bphi, 'b': b, '_FILLER_': -1
        }

        return self.gridue_settings
