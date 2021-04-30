"""
The ``snl`` module contains :class:`SNL` for representing a single-null
topology/configuration.

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
from INGRID.geometry import Point, Line, Patch, trim_geometry, rotate
from collections import OrderedDict


class SNL(TopologyUtils):
    """
    The SNL class for handling `Lower Single Null` (LSN) and
    `Upper Single Null` (USN) configurations within a tokamak.

    Parameters
    ----------
    Ingrid_obj : Ingrid
        Ingrid object the SNL object is being managed by.

    config : str
        String code representing the configuration (for SNL it can be 'LSN' or 'USN').

    Attributes
    ----------
    ConnexionMap : dict
        A mapping defining dependencies between Patch objects for grid generation.

    patches : dict
        The collection of Patch objects representing the topology.
    """

    def __init__(self, Ingrid_obj: 'ingrid.Ingrid', config: str):
        TopologyUtils.__init__(self, Ingrid_obj, config)

        self.ConnexionMap = {
            'A1': {'N': ('A2', 'S')},
            'B1': {'N': ('B2', 'S')},
            'C1': {'N': ('C2', 'S'), 'W': ('B1', 'E')},
            'D1': {'N': ('D2', 'S'), 'W': ('C1', 'E')},
            'E1': {'N': ('E2', 'S'), 'W': ('D1', 'E'), 'E': ('B1', 'W')},
            'F1': {'N': ('F2', 'S'), 'W': ('A1', 'E')},
            'A2': None,
            'B2': {'W': ('A2', 'E')},
            'C2': {'W': ('B2', 'E')},
            'D2': {'W': ('C2', 'E')},
            'E2': {'W': ('D2', 'E')},
            'F2': {'W': ('E2', 'E')},
        }

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

            # Adjust cell grid face along vertical plane
            poloidal_tag, radial_tag = patch.get_tag()
            if poloidal_tag == 'C':
                patch.AdjustBorder('E', [p for p in self.patches.values() if p.get_tag() == 'D' + radial_tag][0])

            # Circular patch configuration requires adjustment of border to close loop.
            # Convention chosen: 'E' indicates closed loop

            try:
                if patch.TerminatesLoop:
                    # Get patch name of adjacent patch for linking boundary points
                    pname = self.PatchTagMap[self.ConnexionMap.get(patch.get_tag())['E'][0]]
                    patch.AdjustBorder('E', self.patches[pname])
            except:
                pass

    def AdjustPatch(self, patch):
        primary_xpt = Point(self.LineTracer.NSEW_lookup['xpt1']['coor']['center'])

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

    def construct_patches(self):
        """
        Create the Patch map with :class:`LineTracing`.

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

        if self.settings['grid_settings']['patch_generation']['strike_pt_loc'] == 'limiter':
            WestPlate = self.parent.LimiterData.copy()
            EastPlate = self.parent.LimiterData.copy()

        else:
            WestPlate = self.PlateData['plate_W1']
            EastPlate = self.PlateData['plate_E1']

        xpt = self.LineTracer.NSEW_lookup['xpt1']['coor']
        magx = np.array([self.settings['grid_settings']['rmagx'] + self.settings['grid_settings']['patch_generation']['rmagx_shift'],
            self.settings['grid_settings']['zmagx'] + self.settings['grid_settings']['patch_generation']['zmagx_shift']])

        psi_1 = self.settings['grid_settings']['psi_1']
        psi_core = self.settings['grid_settings']['psi_core']
        psi_pf_1 = self.settings['grid_settings']['psi_pf_1']

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(magx_tilt_1), magx[1] - 1e6 * np.sin(magx_tilt_1))
        # RHS_Point = Point(magx[0] + 1e6 * np.cos(magx_tilt_1), magx[1] + 1e6 * np.sin(magx_tilt_1))
        RHS_Point = Point(magx[0], magx[1])
        midline_1 = Line([LHS_Point, RHS_Point])

        # LHS_Point = Point(magx[0] - 1e6 * np.cos(magx_tilt_2), magx[1] - 1e6 * np.sin(magx_tilt_2))
        LHS_Point = Point(magx[0], magx[1])
        RHS_Point = Point(magx[0] + 1e6 * np.cos(magx_tilt_2), magx[1] + 1e6 * np.sin(magx_tilt_2))
        midline_2 = Line([LHS_Point, RHS_Point])

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # If USN, we swap east and west lines
        if self.config == 'USN':
            temp = midline_1.copy()
            midline_1 = midline_2
            midline_2 = temp

        # Drawing Separatrix
        xptNW_midLine = self.LineTracer.draw_line(xpt['NW'], {'line': midline_1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        xptN_psiMinCore = self.LineTracer.draw_line(xpt['N'], {'psi': psi_core}, option='rho', direction='cw', show_plot=visual, text=verbose)
        xptNE_midLine = self.LineTracer.draw_line(xpt['NE'], {'line': midline_2}, option='theta', direction='ccw', show_plot=visual, text=verbose)

        # Drawing Lower-SNL region
        if self.settings['grid_settings']['patch_generation']['use_xpt1_W']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_W_tilt']
            d = 'ccw'
            xptW_psiMax = self.LineTracer.draw_line(xpt['W'], {'psi_horizontal': (psi_1, tilt)}, option='z_const', direction=d, show_plot=visual, text=verbose)
        else:
            xptW_psiMax = self.LineTracer.draw_line(xpt['W'], {'psi': psi_1}, option='rho', direction='ccw', show_plot=visual, text=verbose)

        if self.settings['grid_settings']['patch_generation']['use_xpt1_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_E_tilt']
            d = 'cw'
            xptE_psiMax = self.LineTracer.draw_line(xpt['E'], {'psi_horizontal': (psi_1, tilt)}, option='z_const', direction=d, show_plot=visual, text=verbose)
        else:
            xptE_psiMax = self.LineTracer.draw_line(xpt['E'], {'psi': psi_1}, option='rho', direction='ccw', show_plot=visual, text=verbose)

        xpt_WestPlate = self.LineTracer.draw_line(xpt['SW'], {'line': WestPlate}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        xptS_psiMinPF = self.LineTracer.draw_line(xpt['S'], {'psi': psi_pf_1}, option='rho', direction='cw', show_plot=visual, text=verbose)
        xpt_EastPlate = self.LineTracer.draw_line(xpt['SE'], {'line': EastPlate}, option='theta', direction='cw', show_plot=visual, text=verbose)
        iPsiMax_TP = self.LineTracer.draw_line(xptW_psiMax.p[-1], {'line': WestPlate}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        psiMinPF_WestPlate = self.LineTracer.draw_line(xptS_psiMinPF.p[-1], {'line': WestPlate}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        oPsiMax_TP = self.LineTracer.draw_line(xptE_psiMax.p[-1], {'line': EastPlate}, option='theta', direction='cw', show_plot=visual, text=verbose)
        psiMinPF_EastPlate = self.LineTracer.draw_line(xptS_psiMinPF.p[-1], {'line': EastPlate}, option='theta', direction='cw', show_plot=visual, text=verbose)

        imidLine_topLine = self.LineTracer.draw_line(xptNW_midLine.p[-1], {'line': topLine}, option='theta',
            direction='cw', show_plot=visual, text=verbose)

        omidLine_topLine = self.LineTracer.draw_line(xptNE_midLine.p[-1], {'line': topLine}, option='theta',
            direction='ccw', show_plot=visual, text=verbose)

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
        A2_W = (WestPlate.split(A2_S.p[-1])[1]).split(A2_N.p[0], add_split_point=True)[0]
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patch_name='A2', plate_patch=True, plate_location=location)

        # A1 Patch
        location = 'W'
        A1_N = A2_S.reverse_copy()

        A1_S = psiMinPF_WestPlate
        A1_E = xptS_psiMinPF
        A1_W = (WestPlate.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point=True)[0]
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patch_name='A1', plate_patch=True, plate_location=location)

        # B2 Patch

        B2_N = self.LineTracer.draw_line(A2_N.p[-1], {'line': midline_1}, option='theta', direction='cw', show_plot=visual, text=verbose)
        B2_S = xptNW_midLine.reverse_copy()
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B2_W = xptW_psiMax
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patch_name='B2')

        # B1 Patch
        B1_N = B2_S.reverse_copy()
        B1_S = self.LineTracer.draw_line(xptN_psiMinCore.p[-1], {'line': midline_1}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()
        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B1_W = xptN_psiMinCore.reverse_copy()
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patch_name='B1')

        # C2 Patch
        C2_N = self.LineTracer.draw_line(B2_N.p[-1], {'line': topLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        C2_S = imidLine_topLine.reverse_copy()
        C2_E = Line([C2_N.p[-1], C2_S.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patch_name='C2')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.LineTracer.draw_line(B1_S.p[0], {'line': topLine}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()
        C1_E = Line([C1_N.p[-1], C1_S.p[0]])
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patch_name='C1')

        # F2 Patch
        location = 'E'
        F2_N = oPsiMax_TP
        F2_S = xpt_EastPlate.reverse_copy()
        F2_E = (EastPlate.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point=True)[0]
        F2_W = xptE_psiMax
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patch_name='F2', plate_patch=True, plate_location=location)

        # F1 Patch
        location = 'E'
        F1_N = F2_S.reverse_copy()
        F1_S = psiMinPF_EastPlate.reverse_copy()
        F1_E = (EastPlate.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point=True)[0]
        F1_W = xptS_psiMinPF.reverse_copy()
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patch_name='F1', plate_patch=True, plate_location=location)

        # E2 Patch
        E2_N = self.LineTracer.draw_line(F2_N.p[0], {'line': midline_2}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        E2_S = xptNE_midLine
        E2_E = xptE_psiMax.reverse_copy()
        E2_W = Line([E2_S.p[-1], E2_N.p[0]])
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patch_name='E2')

        # E1 Patch
        E1_N = E2_S.reverse_copy()
        E1_S = self.LineTracer.draw_line(xptN_psiMinCore.p[-1], {'line': midline_2}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        E1_E = xptN_psiMinCore
        E1_W = Line([E1_S.p[-1], E1_N.p[0]])
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patch_name='E1')

        # D2 Patch
        D2_N = self.LineTracer.draw_line(E2_N.p[0], {'line': topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        D2_S = omidLine_topLine
        D2_E = Line([D2_N.p[-1], D2_S.p[0]])
        D2_W = Line([D2_S.p[-1], D2_N.p[0]])
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patch_name='D2')

        # D1 Patch
        D1_N = D2_S.reverse_copy()
        D1_S = self.LineTracer.draw_line(E1_S.p[-1], {'line': topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        D1_E = Line([D1_N.p[-1], D1_S.p[0]])
        D1_W = Line([D1_S.p[-1], D1_N.p[0]])
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patch_name='D1')

        patches = [A2, B2, C2, D2, E2, F2, A1, F1, B1, C1, D1, E1]

        self.patches = {}
        for patch in patches:
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patch_name] = patch
        self.OrderPatches()

    def GroupPatches(self):
        p = self.patches
        self.PatchGroup = {'SOL': (p['A2'], p['B2'], p['C2'], p['D2'], p['E2'], p['F2']),
                           'CORE': (p['B1'], p['C1'], p['D1'], p['E1']),
                           'PF': (p['A1'], p['F1'])}

    def OrderPatches(self):
        patches = ['A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'A1', 'F1', 'B1', 'C1', 'D1', 'E1']
        self.patches = OrderedDict([(pname, self.patches[pname]) for pname in patches])

    def set_gridue(self):
        """
        Prepares a ``gridue_settings`` dictionary with required data
        for writing a gridue file.

        Parameters
        ----------

        Returns
        -------

        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patches['A2'].npol - 1
        ixpt2 = ixrb - self.patches['F2'].npol + 1
        iyseparatrix1 = self.patches['A1'].nrad - 1
        nxm = len(self.rm) - 2
        nym = len(self.rm[0]) - 2

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
            'nxm': nxm, 'nym': nym, 'ixpt1': ixpt1, 'ixpt2': ixpt2, 'iyseptrx1': iyseparatrix1,
            'rm': self.rm, 'zm': self.zm, 'psi': psi, 'br': br, 'bz': bz, 'bpol': bpol, 'bphi': bphi, 'b': b
        }

        return self.gridue_settings
