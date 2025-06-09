"""
The ``udn`` module contains :class:`CDN` for representing a connected
double-null topology/configuration.

Child of base :class:`utils.TopologyUtils`.

"""
import numpy as np
import matplotlib
try:
    matplotlib.use("TkAgg")
except:
    pass
from INGRID.utils import TopologyUtils
from INGRID.geometry import Point, Line, Patch, trim_geometry
from collections import OrderedDict


class CDN(TopologyUtils):
    """
    The `CDN` class for handling `Connected Double Null` configurations
    within a tokamak.

    Parameters
    ----------
    Ingrid_obj : Ingrid
        Ingrid object the CDN object is being managed by.

    config : str, optional
        String code representing the configuration.

    Attributes
    ----------
    ConnexionMap : dict
        A mapping defining dependencies between Patch objects for grid generation.

    patches : dict
        The collection of Patch objects representing the topology.
    """
    def __init__(self, Ingrid_obj: 'ingrid.Ingrid', config: 'CDN'):
        TopologyUtils.__init__(self, Ingrid_obj, config)

        self.ConnexionMap = {

            'B1': None,
            'C1': {'W': ('B1', 'E')},
            'F1': {'W': ('C1', 'E')},
            'G1': {'W': ('F1', 'E'), 'E': ('B1', 'W')},

            'A1': None,
            'D1': {'W': ('E1', 'E')},
            'E1': None,
            'H1': {'W': ('A1', 'E')},

            'A2': {'S': ('A1', 'N')},
            'B2': {'S': ('B1', 'N'), 'W': ('A2', 'E')},
            'C2': {'S': ('C1', 'N'), 'W': ('B2', 'E')},
            'D2': {'S': ('D1', 'N'), 'W': ('C2', 'E')},
            'E2': {'S': ('E1', 'N')},
            'F2': {'S': ('F1', 'N'), 'W': ('E2', 'E')},
            'G2': {'S': ('G1', 'N'), 'W': ('F2', 'E')},
            'H2': {'S': ('H1', 'N'), 'W': ('G2', 'E')},
        }

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
        xpt2_psi = self.PsiNorm.get_psi(xpt2['center'][0], xpt2['center'][1])

        magx = np.array([self.settings['grid_settings']['rmagx'] + self.settings['grid_settings']['patch_generation']['rmagx_shift'],
            self.settings['grid_settings']['zmagx'] + self.settings['grid_settings']['patch_generation']['zmagx_shift']])

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

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Tracing primary-separatrix: core-boundary
        xpt1N__psiMinCore = self.LineTracer.draw_line(xpt1['N'], {'psi': psi_core},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        B1_W = xpt1N__psiMinCore.reverse_copy()
        G1_E = B1_W.reverse_copy()

        xpt1NW__midline_1 = self.LineTracer.draw_line(xpt1['NW'], {'line': midline_1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        B1_N = xpt1NW__midline_1
        B2_S = xpt1NW__midline_1.reverse_copy()

        xpt1NE__midline_2 = self.LineTracer.draw_line(xpt1['NE'], {'line': midline_2},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        G1_N = xpt1NE__midline_2.reverse_copy()
        G2_S = xpt1NE__midline_2

        # Tracing core: psi-min-core region (rho = 1)

        psiMinCore__midline_1_core = self.LineTracer.draw_line(xpt1N__psiMinCore.p[-1], {'line': midline_1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        B1_S = psiMinCore__midline_1_core.reverse_copy()

        psiMinCore__midline_2_core = self.LineTracer.draw_line(xpt1N__psiMinCore.p[-1], {'line': midline_2},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        G1_S = psiMinCore__midline_2_core
        xpt1__psiMinPF1 = self.LineTracer.draw_line(xpt1['S'], {'psi': psi_pf_1},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        A1_E = xpt1__psiMinPF1
        H1_W = A1_E.reverse_copy()
        psiMinPF1__WestPlate1 = self.LineTracer.draw_line(xpt1__psiMinPF1.p[-1], {'line': WestPlate1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        A1_S = psiMinPF1__WestPlate1

        psiMinPF1__EastPlate1 = self.LineTracer.draw_line(xpt1__psiMinPF1.p[-1], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        H1_S = psiMinPF1__EastPlate1.reverse_copy()

        xpt1__WestPlate1 = self.LineTracer.draw_line(xpt1['SW'], {'line': WestPlate1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        A2_S = xpt1__WestPlate1
        A1_N = A2_S.reverse_copy()

        # I1_N / I2_S
        xpt1__EastPlate1 = self.LineTracer.draw_line(xpt1['SE'], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        H1_N = xpt1__EastPlate1
        H2_S = H1_N.reverse_copy()

        # Drawing the portion of the separatrix found in the double-null configuration
        xpt2N__outer_core = self.LineTracer.draw_line(xpt2['N'], {'psi': psi_core},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        C1_E = xpt2N__outer_core
        F1_W = C1_E.reverse_copy()

        C1_S = self.LineTracer.draw_line(C1_E.p[-1], {'line': midline_1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)

        F1_S = self.LineTracer.draw_line(C1_E.p[-1], {'line': midline_1}, option='theta',
            direction='cw', show_plot=visual, text=verbose).reverse_copy()

        C1_N = self.LineTracer.draw_line(xpt2['N'], {'line': midline_1},
            option='theta', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        C2_S = C1_N.reverse_copy()

        F1_N = self.LineTracer.draw_line(xpt2['N'], {'line': midline_2},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        F2_S = F1_N.reverse_copy()

        xpt2__psiMinPF2 = self.LineTracer.draw_line(xpt2['E'], {'psi': psi_pf_1},
            option='rho', direction='cw', show_plot=visual, text=verbose)
        E1_E = xpt2__psiMinPF2
        D1_W = E1_E.reverse_copy()

        E1_S = self.LineTracer.draw_line(E1_E.p[-1], {'line': EastPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose)

        D1_S = self.LineTracer.draw_line(E1_E.p[-1], {'line': WestPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose).reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt2_W']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt2_W_tilt']
            E2_E = self.LineTracer.draw_line(xpt2['W'], {'psi_horizontal': (psi_1, tilt)}, option='z_const', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        else:
            E2_E = self.LineTracer.draw_line(xpt2['W'], {'psi': psi_1}, option='rho', direction='ccw', show_plot=visual, text=verbose).reverse_copy()
        F2_W = E2_E.reverse_copy()
        if self.settings['grid_settings']['patch_generation']['use_xpt2_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt2_E_tilt']
            D2_W = self.LineTracer.draw_line(xpt2['E'], {'psi_horizontal': (psi_2, tilt)}, option='z_const', direction='ccw', show_plot=visual, text=verbose)
        else:
            D2_W = self.LineTracer.draw_line(xpt2['E'], {'psi': psi_2}, option='rho', direction='ccw', show_plot=visual, text=verbose)
        C2_E = D2_W.reverse_copy()

        E1_N = self.LineTracer.draw_line(xpt2['S'], {'line': EastPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()
        E2_S = E1_N.reverse_copy()

        D1_N = self.LineTracer.draw_line(xpt2['S'], {'line': WestPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)
        D2_S = D1_N.reverse_copy()

        F2_N = self.LineTracer.draw_line(F2_W.p[-1], {'line': midline_2}, option='theta',
            direction='cw', show_plot=visual, text=verbose)

        C2_N = self.LineTracer.draw_line(C2_E.p[0], {'line': midline_1}, option='theta',
            direction='ccw', show_plot=visual, text=verbose).reverse_copy()

        D2_N = self.LineTracer.draw_line(D2_W.p[-1], {'line': WestPlate2}, option='theta', direction='cw',
            show_plot=visual, text=verbose)

        E2_N = self.LineTracer.draw_line(F2_W.p[-1], {'line': EastPlate2}, option='theta', direction='ccw',
            show_plot=visual, text=verbose).reverse_copy()


        midline_2__EastPlate1 = self.LineTracer.draw_line(F2_N.p[-1], {'line': EastPlate1},
            option='theta', direction='cw', show_plot=visual, text=verbose)
        midline_1__WestPlate1 = self.LineTracer.draw_line(C2_N.p[0], {'line': WestPlate1},
            option='theta', direction='ccw', show_plot=visual, text=verbose)
        if self.settings['grid_settings']['patch_generation']['use_xpt1_E']:
            tilt = self.settings['grid_settings']['patch_generation']['xpt1_E_tilt']
            H2_W = self.LineTracer.draw_line(xpt1['E'], {'line': (midline_2__EastPlate1, tilt)}, option='z_const', direction='cw', show_plot=visual, text=verbose)
        else:
            H2_W = self.LineTracer.draw_line(xpt1['E'], {'line': midline_2__EastPlate1},
            option='rho', direction='ccw', show_plot=visual, text=verbose)
        G2_E = H2_W.reverse_copy()

        if self.settings['grid_settings']['patch_generation']['use_xpt1_W']:
            tilt = -self.settings['grid_settings']['patch_generation']['xpt1_W_tilt']
            B2_W = self.LineTracer.draw_line(xpt1['W'], {'line': (midline_1__WestPlate1, tilt)}, option='z_const', direction='ccw', show_plot=visual, text=verbose)
        else:
            B2_W = self.LineTracer.draw_line(xpt1['W'], {'line': midline_1__WestPlate1},
            option='rho', direction='ccw', show_plot=visual, text=verbose)

        A2_E = B2_W.reverse_copy()

        A2_N, B2_N = midline_1__WestPlate1.reverse_copy().split(B2_W.p[-1], add_split_point=True)
        G2_N, H2_N = midline_2__EastPlate1.split(H2_W.p[-1], add_split_point=True)
        B1_E = self.LineTracer.draw_line(B1_N.p[-1], {'psi_horizontal': psi_core}, option='z_const',
            direction='ccw', show_plot=visual, text=verbose)
        C1_W = B1_E.reverse_copy()

        B2_E = self.LineTracer.draw_line(B2_N.p[-1], {'psi_horizontal': 1.0}, option='z_const',
            direction='ccw', show_plot=visual, text=verbose)
        C2_W = B2_E.reverse_copy()
        F1_E = self.LineTracer.draw_line(F1_N.p[-1], {'psi_horizontal': psi_core}, option='z_const',
            direction='cw', show_plot=visual, text=verbose)
        G1_W = F1_E.reverse_copy()

        F2_E = self.LineTracer.draw_line(F2_N.p[-1], {'psi_horizontal': 1.0}, option='z_const',
            direction='cw', show_plot=visual, text=verbose)
        G2_W = F2_E.reverse_copy()

        visual=True



        A1_W = trim_geometry(WestPlate1, A1_S.p[-1], A1_N.p[0])
        A2_W = trim_geometry(WestPlate1, A2_S.p[-1], A2_N.p[0])

        D1_E = trim_geometry(WestPlate2, D1_N.p[-1], D1_S.p[0])
        D2_E = trim_geometry(WestPlate2, D2_N.p[-1], D2_S.p[0])

        E1_W = trim_geometry(EastPlate2, E1_S.p[-1], E1_N.p[0])
        E2_W = trim_geometry(EastPlate2, E2_S.p[-1], E2_N.p[0])

        H1_E = trim_geometry(EastPlate1, H1_N.p[-1], H1_S.p[0])
        H2_E = trim_geometry(EastPlate1, H2_N.p[-1], H2_S.p[0])

        # ============== Patch A1 ==============
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patch_name='A1', plate_patch=True, plate_location='W')
        # ============== Patch A2 ==============
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patch_name='A2', plate_patch=True, plate_location='W')
        
        # ============== Patch B1 ==============
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patch_name='B1')
        # ============== Patch B2 ==============
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patch_name='B2')
        
        # ============== Patch C1 ==============
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patch_name='C1')
        # ============== Patch C2 ==============
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patch_name='C2')
        
        # ============== Patch D1 ==============
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patch_name='D1', plate_patch=True, plate_location='E')
        # ============== Patch D2 ==============
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patch_name='D2', plate_patch=True, plate_location='E')
        
        # ============== Patch E1 ==============
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patch_name='E1', plate_patch=True, plate_location='W')
        # ============== Patch E2 ==============
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patch_name='E2', plate_patch=True, plate_location='W')
       
        # ============== Patch F1 ==============
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patch_name='F1')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patch_name='F2')
       
        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patch_name='G1')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patch_name='G2')
        
        # ============== Patch H1 ==============
        H1 = Patch([H1_N, H1_E, H1_S, H1_W], patch_name='H1', plate_patch=True, plate_location='E')
        # ============== Patch H2 ==============
        H2 = Patch([H2_N, H2_E, H2_S, H2_W], patch_name='H2', plate_patch=True, plate_location='E')
        
        patches = [
                   B1, C1, F1, G1, 
                   A1, D1, E1, H1,
                   A2, B2, C2, D2, E2, F2, G2, H2,
                ]
        self.patches = {}
        for patch in patches:
            patch.parent = self
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patch_name] = patch
        self.OrderPatches()

    def OrderPatches(self):

        patches = [
                   'B1', 'C1', 'F1', 'G1', 
                   'A1', 'E1', 'D1', 'H1',
                   'A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2',
        ]

        self.patches = OrderedDict([(pname, self.patches[pname]) for pname in patches])

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
            #poloidal_tag, radial_tag = patch.get_tag()
            #if poloidal_tag == 'C' and radial_tag in ['1', '2']:
            #    patch.AdjustBorder('E', self.patches['F' + radial_tag])

            # Circular patch configuration requires adjustment of border to close loop.
            # Convention chosen: 'E' indicates closed loop
            try :
              pname = self.PatchTagMap[self.ConnexionMap.get(patch.get_tag())['E'][0]]
              patch.AdjustBorder('E', self.patches[pname])
            except :
              pass
            #try:
            #    if patch.TerminatesLoop:
            #        # Get patch name of adjacent patch for linking boundary points
            #        pname = self.PatchTagMap[self.ConnexionMap.get(patch.get_tag())['E'][0]]
            #        patch.AdjustBorder('E', self.patches[pname])
            #except:
            #    pass

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
        elif tag == 'C2':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'C1':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'D2':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'D1':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'E2':
            patch.adjust_corner(xpt2, 'SE')
        elif tag == 'E1':
            patch.adjust_corner(xpt2, 'NE')
        elif tag == 'F2':
            patch.adjust_corner(xpt2, 'SW')
        elif tag == 'F1':
            patch.adjust_corner(xpt2, 'NW')
        elif tag == 'G2':
            patch.adjust_corner(xpt1, 'SE')
        elif tag == 'G1':
            patch.adjust_corner(xpt1, 'NE')
        elif tag == 'H2':
            patch.adjust_corner(xpt1, 'SW')
        elif tag == 'H1':
            patch.adjust_corner(xpt1, 'NW')

    def GroupPatches(self):
        # p = self.patches
        # self.PatchGroup = {'SOL' : [],
        # 'CORE' : (p['B1'], p['C1'], p['D1'], p['E1']),
        # 'PF' : (p['A1'], p['F1'])}
        pass

    def set_gridue(self) -> dict:
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
        iyseparatrix2 = self.patches['C1'].nrad - 1
        iyseparatrix3 = iyseparatrix2
        iyseparatrix4 = iyseparatrix1

        ix_plate1 = 0
        ix_cut1 = self.patches['A1'].npol - 1

        ix_cut2 = 0
        for alpha in ['A', 'B', 'C']:
            ix_cut2 += self.patches[alpha + '1'].npol - 1

        ix_plate2 = 0
        for alpha in ['A', 'B', 'C', 'D']:
            ix_plate2 += self.patches[alpha + '2'].npol - 1

        ix_plate3 = ix_plate2 + 2

        ix_cut3 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E']:
            ix_cut3 += self.patches[alpha + '2'].npol - 1
        ix_cut3 += 2

        ix_cut4 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
            ix_cut4 += self.patches[alpha + '1'].npol - 1
        ix_cut4 += 2

        ix_plate4 = 0
        for alpha in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
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
