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
from topology_utils import Topology
from geometry import Point, Line, Patch, trim_geometry
from collections import OrderedDict


class SNL(Topology):
    """
    The SNL (Single-Null) class is the parent class for both upper-single null (USN)
    and lower-single null (LSN) configurations.
    This base class handles the formatting and plotting of data obtained from an LSN or USN
    object.
    Parameter:
        - INGRID_object : Ingrid class object
        All SNL objects are children of the main Ingrid class. INGRID_object provides
        information such as YAML data, efit_psi, psi_norm, and PlateData.
    @author: garcia299
    """

    def __init__(self, Ingrid_obj, config):
        Topology.__init__(self, Ingrid_obj, config)

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
            west_tilt = self.settings['grid_settings']['patch_generation']['west_tilt']
        except KeyError:
            west_tilt = 0.0
        try:
            east_tilt = self.settings['grid_settings']['patch_generation']['east_tilt']
        except KeyError:
            east_tilt = 0.0

        self.RefreshSettings()

        if self.settings['limiter']['use_limiter']:
            WestPlate = self.parent.LimiterData.copy()
            EastPlate = self.parent.LimiterData.copy()

        else:
            WestPlate = self.PlateData['plate_W1']
            EastPlate = self.PlateData['plate_E1']

        xpt = self.LineTracer.NSEW_lookup['xpt1']['coor']
        magx = np.array([self.settings['grid_settings']['rmagx'] + self.settings['grid_settings']['patch_generation']['rmagx_shift'],
            self.settings['grid_settings']['zmagx'] + self.settings['grid_settings']['patch_generation']['zmagx_shift']])

        psi_max = self.settings['grid_settings']['psi_max']
        psi_core = self.settings['grid_settings']['psi_core']
        psi_pf_1 = self.settings['grid_settings']['psi_pf_1']

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(west_tilt), magx[1] - 1e6 * np.sin(west_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(west_tilt), magx[1] + 1e6 * np.sin(west_tilt))
        west_midLine = Line([LHS_Point, RHS_Point])

        LHS_Point = Point(magx[0] - 1e6 * np.cos(east_tilt), magx[1] - 1e6 * np.sin(east_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(east_tilt), magx[1] + 1e6 * np.sin(east_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # If USN, we swap east and west lines
        if self.config == 'USN':
            temp = west_midLine.copy()
            west_midLine = east_midLine
            east_midLine = temp

        # Drawing Separatrix
        xptNW_midLine = self.LineTracer.draw_line(xpt['NW'], {'line': west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        xptN_psiMinCore = self.LineTracer.draw_line(xpt['N'], {'psi': psi_core}, option='rho', direction='cw', show_plot=visual, text=verbose)
        xptNE_midLine = self.LineTracer.draw_line(xpt['NE'], {'line': east_midLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)

        # Drawing Lower-SNL region
        if self.settings['grid_settings']['patch_generation']['use_NW']:
            tilt = self.settings['grid_settings']['patch_generation']['NW_adjust']
            xptW_psiMax = self.LineTracer.draw_line(rotate(xpt['W'], tilt, xpt['center']), {'psi_horizontal': (psi_max, tilt)}, option='z_const', direction='ccw', show_plot=visual, text=verbose)
        else:
            xptW_psiMax = self.LineTracer.draw_line(xpt['W'], {'psi': psi_max}, option='rho', direction='ccw', show_plot=visual, text=verbose)

        if self.settings['grid_settings']['patch_generation']['use_NE']:
            tilt = self.settings['grid_settings']['patch_generation']['NE_adjust']
            xptE_psiMax = self.LineTracer.draw_line(rotate(xpt['E'], tilt, xpt['center']), {'psi_horizontal': (psi_max, tilt)}, option='z_const', direction='cw', show_plot=visual, text=verbose)
        else:
            xptE_psiMax = self.LineTracer.draw_line(xpt['E'], {'psi': psi_max}, option='rho', direction='ccw', show_plot=visual, text=verbose)

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

        # Integrating horizontally along mid-line towards psiMax and psiMinCore

        imidLine_psiMax = self.LineTracer.draw_line(xptNW_midLine.p[-1], {'psi_horizontal': (psi_max, west_tilt)}, option='z_const',
                direction='ccw' if self.config == 'LSN' else 'cw', show_plot=visual, text=verbose)
        imidLine_psiMinCore = self.LineTracer.draw_line(xptNW_midLine.p[-1], {'psi_horizontal': (psi_core, west_tilt)}, option='z_const',
                direction='cw' if self.config == 'LSN' else 'ccw', show_plot=visual, text=verbose)
        omidLine_psiMax = self.LineTracer.draw_line(xptNE_midLine.p[-1], {'psi_horizontal': (psi_max, east_tilt)}, option='z_const',
                direction='cw' if self.config == 'LSN' else 'ccw', show_plot=visual, text=verbose)
        omidLine_psiMinCore = self.LineTracer.draw_line(xptNE_midLine.p[-1], {'psi_horizontal': (psi_core, east_tilt)}, option='z_const',
                direction='ccw' if self.config == 'LSN' else 'cw', show_plot=visual, text=verbose)

        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.LineTracer.draw_line(omidLine_topLine.p[-1], {'psi_vertical': psi_max}, option='r_const',
                direction='cw' if self.config == 'LSN' else 'ccw', show_plot=visual, text=verbose)
        topLine_psiMinCore = self.LineTracer.draw_line(omidLine_topLine.p[-1], {'psi_vertical': psi_core}, option='r_const',
                direction='ccw' if self.config == 'LSN' else 'cw', show_plot=visual, text=verbose)

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
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patchName='IDL' if self.config == 'LSN' else 'ODL',
            platePatch=True, plateLocation=location)

        # A1 Patch
        location = 'W'
        A1_N = A2_S.reverse_copy()

        A1_S = psiMinPF_WestPlate
        A1_E = xptS_psiMinPF
        A1_W = (WestPlate.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point=True)[0]
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patchName='IPF' if self.config == 'LSN' else 'OPF',
            platePatch=True, plateLocation=location)

        # B2 Patch

        B2_N = self.LineTracer.draw_line(A2_N.p[-1], {'line': west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        B2_S = xptNW_midLine.reverse_copy()
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B2_W = xptW_psiMax
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patchName='ISB' if self.config == 'LSN' else 'OSB')

        # B1 Patch
        B1_N = B2_S.reverse_copy()
        B1_S = self.LineTracer.draw_line(xptN_psiMinCore.p[-1], {'line': west_midLine}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()
        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B1_W = xptN_psiMinCore.reverse_copy()
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patchName='ICB' if self.config == 'LSN' else 'OCB')

        # C2 Patch
        C2_N = self.LineTracer.draw_line(B2_N.p[-1], {'line': topLine}, option='theta', direction='cw', show_plot=visual, text=verbose)
        C2_S = imidLine_topLine.reverse_copy()
        C2_E = Line([C2_N.p[-1], C2_S.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patchName='IST' if self.config == 'LSN' else 'OST')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.LineTracer.draw_line(B1_S.p[0], {'line': topLine}, option='theta', direction='cw', show_plot=visual, text=verbose).reverse_copy()
        C1_E = Line([C1_N.p[-1], C1_S.p[0]])
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patchName='ICT' if self.config == 'LSN' else 'OCT')

        # F2 Patch
        location = 'E'
        F2_N = oPsiMax_TP
        F2_S = xpt_EastPlate.reverse_copy()
        F2_E = (EastPlate.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point=True)[0]
        F2_W = xptE_psiMax
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName='ODL' if self.config == 'LSN' else 'IDL',
            platePatch=True, plateLocation=location)

        # F1 Patch
        location = 'E'
        F1_N = F2_S.reverse_copy()
        F1_S = psiMinPF_EastPlate.reverse_copy()
        F1_E = (EastPlate.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point=True)[0]
        F1_W = xptS_psiMinPF.reverse_copy()
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName='OPF' if self.config == 'LSN' else 'IPF',
            platePatch=True, plateLocation=location)

        # E2 Patch
        E2_N = self.LineTracer.draw_line(F2_N.p[0], {'line': east_midLine}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        E2_S = xptNE_midLine
        E2_E = xptE_psiMax.reverse_copy()
        E2_W = Line([E2_S.p[-1], E2_N.p[0]])
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName='OSB' if self.config == 'LSN' else 'ISB')

        # E1 Patch
        E1_N = E2_S.reverse_copy()
        E1_S = self.LineTracer.draw_line(xptN_psiMinCore.p[-1], {'line': east_midLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        E1_E = xptN_psiMinCore
        E1_W = Line([E1_S.p[-1], E1_N.p[0]])
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName='OCB' if self.config == 'LSN' else 'ICB')

        # D2 Patch
        D2_N = self.LineTracer.draw_line(E2_N.p[0], {'line': topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        D2_S = omidLine_topLine
        D2_E = Line([D2_N.p[-1], D2_S.p[0]])
        D2_W = Line([D2_S.p[-1], D2_N.p[0]])
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName='OST' if self.config == 'LSN' else 'IST')

        # D1 Patch
        D1_N = D2_S.reverse_copy()
        D1_S = self.LineTracer.draw_line(E1_S.p[-1], {'line': topLine}, option='theta', direction='ccw', show_plot=visual, text=verbose)
        D1_E = Line([D1_N.p[-1], D1_S.p[0]])
        D1_W = Line([D1_S.p[-1], D1_N.p[0]])
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName='OCT' if self.config == 'LSN' else 'ICT')

        patches = [A2, B2, C2, D2, E2, F2, A1, F1, B1, C1, D1, E1]

        self.patches = {}
        for patch in patches:
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch
        self.OrderPatches()

    def GroupPatches(self):
        p = self.patches
        self.PatchGroup = {'SOL': (p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL']),
        'CORE': (p['ICB'], p['ICT'], p['OCT'], p['OCB']),
        'PF': (p['IPF'], p['OPF'])}

    def OrderPatches(self):
        if self.config == 'LSN':
            patches = ['IDL', 'ISB', 'IST', 'OST', 'OSB', 'ODL', 'IPF', 'OPF', 'ICB', 'ICT', 'OCT', 'OCB']
        else:
            patches = ['ODL', 'OSB', 'OST', 'IST', 'ISB', 'IDL', 'OPF', 'IPF', 'OCB', 'OCT', 'ICT', 'ICB']

        self.patches = OrderedDict([(pname, self.patches[pname]) for pname in patches])

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

        self.gridue_settings = {'nxm': nxm, 'nym': nym, 'ixpt1': ixpt1, 'ixpt2': ixpt2, 'iyseptrx1': iyseparatrix1,
            'rm': self.rm, 'zm': self.zm, 'psi': psi, 'br': br, 'bz': bz, 'bpol': bpol, 'bphi': bphi, 'b': b}
