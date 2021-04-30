"""
The ``snlh`` module contains :class:`SNLH` for representing a half domain single-null
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


class SNLH(TopologyUtils):
    """
    The SNLH class for handling half domain `Lower Single Null` (LSNH) and
    `Upper Single Null` (USNH) configurations within a tokamak.

    Parameters
    ----------
    Ingrid_obj : Ingrid
        Ingrid object the SNL object is being managed by.

    config : str
        String code representing the configuration (for SNL it can be 'LSNH' or 'USNH').

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
            'C1': {'N': ('C2', 'S'), 'E': ('B1', 'W')},
            'D1': {'N': ('D2', 'S'), 'W': ('A1', 'E')},
            'A2': None,
            'B2': {'W': ('A2', 'E')},
            'C2': None,
            'D2': {'W': ('C2', 'E')},
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
        elif tag == 'C1':
            patch.adjust_corner(primary_xpt, 'NE')
        elif tag == 'C2':
            patch.adjust_corner(primary_xpt, 'SE')
        elif tag == 'D1':
            patch.adjust_corner(primary_xpt, 'NW')
        elif tag == 'D2':
            patch.adjust_corner(primary_xpt, 'SW')

    def construct_patches(self):
        """
        Create the Patch map with :class:`LineTracing`.

        Patch Labeling Key:

        - I: Inner,
        - O: Outer,
        - DL: Divertor Leg,
        - PF: Private Flux,
        - T: Top,
        - B: Bottom,
        - S: Scrape Off Layer,
        - C: Core.

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

        # If USN, we swap east and west lines
        if self.config == 'USNH':
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

        # D2 Patch
        location = 'E'
        D2_N = oPsiMax_TP
        D2_S = xpt_EastPlate.reverse_copy()
        D2_E = (EastPlate.split(D2_N.p[-1])[1]).split(D2_S.p[0], add_split_point=True)[0]
        D2_W = xptE_psiMax
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patch_name='D2', plate_patch=True, plate_location=location)

        # D1 Patch
        location = 'E'
        D1_N = D2_S.reverse_copy()
        D1_S = psiMinPF_EastPlate.reverse_copy()
        D1_E = (EastPlate.split(D1_N.p[-1])[1]).split(D1_S.p[0], add_split_point=True)[0]
        D1_W = xptS_psiMinPF.reverse_copy()
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patch_name='D1', plate_patch=True, plate_location=location)

        # C2 Patch
        C2_N = self.LineTracer.draw_line(D2_N.p[0], {'line': midline_2}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        C2_S = xptNE_midLine
        C2_E = xptE_psiMax.reverse_copy()
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patch_name='C2')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.LineTracer.draw_line(xptN_psiMinCore.p[-1], {'line': midline_2}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        C1_E = xptN_psiMinCore
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patch_name='C1')

        patches = [A2, B2, C2, D2, A1, D1, B1, C1]

        self.patches = {}
        for patch in patches:
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patch_name] = patch
        self.OrderPatches()

    def GroupPatches(self):
        p = self.patches
        self.PatchGroup = {'SOL': (p['A2'], p['B2'], p['C2'], p['D2']),
                           'CORE': (p['B1'], p['C1']),
                           'PF': (p['A1'], p['D1'])}

    def OrderPatches(self):
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

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patches['A2'].npol - 1

        ixpt2 = 0
        for tag in ['A', 'B', 'C']:
            ixpt2 += self.patches[tag + '2'].npol - 1
        ixpt2 += 2  # Account for guard cells

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
