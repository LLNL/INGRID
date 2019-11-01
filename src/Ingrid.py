#!/usr/bin/env pythonu
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:17:21 2019

@author: watkins35, garcia299
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import yaml as _yaml_
import GUI.IngridApp as IA
import Uegrid.uegrid as uegrid

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder
from geometry import Point, Line, SNL_Patch


class Ingrid:
    """ An interactive grid generator for edge plasmas in a tokamak
    Accepts a dictionary generated from a namelist file that contains
    the parameters.
    
    Parameters
    ----------
    yaml : dict
        Params dictionary object contains two dictionaries
        First is 'files' which contains the keys: geqdsk, itp, otp.
        The second is 'grid params' which has: psi_max, psi_min_core,
        psi_min_pf, Rmagx, Zmagx, Rxpt, Zxpt
    
    """

    def __init__(self, params = {}):

        self.default_grid_params = { \
            'config' : '', 'num_xpt' : 1, \
            'np_global' : 3, 'nr_global' : 2, \
            'psi_max' : 0.0, 'psi_max_r' : 0.0, 'psi_max_z' : 0.0, \
            'psi_min_core' : 0.0, 'psi_min_core_r' : 0.0, 'psi_min_core_z' : 0.0, \
            'psi_min_pf' : 0.0, 'psi_min_pf_r' : 0.0, 'psi_min_pf_z' : 0.0, \
            'rmagx' : 0.0, 'zmagx' : 0.0, \
            'rxpt' : 0.0, 'zxpt' : 0.0 \
        }

        self.default_integrator_params = { \
            'dt' : 0.01, 'eps' : 5e-5, \
            'first_step' : 1e-5, 'step_ratio' : 0.02, \
            'tol' : 5e-3 \
        }

        self.default_target_plates_params = { \
            'plate_E1' : {'file' : '', 'name' : '', 'np_local' : 3, 'nr_local' : 2, 'poloidal_f' : 'x, x'}, \
            'plate_W1' : {'file' : '', 'name' : '', 'np_local' : 3, 'nr_local' : 2, 'poloidal_f' : 'x, x'}
        }
        # Process params as a YAML file. Ensure ends up as a raw dictionary.s
        self.process_yaml(params)
        self.get_yaml()

    def process_yaml(self, params):

            yaml_lookup = ['eqdsk', 'grid_params', 'integrator_params', 'target_plates']

            for key in yaml_lookup:
                try:
                    # Check YAML file (if passed into INGRID)
                    params[key]
                except KeyError:
                    # set clean
                    print('Could not find "{}" in YAML file.'.format(key))
                    params[key] = self.get_default_values(key)
                
                if key == 'eqdsk':
                    self.eqdsk = params[key]
                elif key == 'grid_params':
                    self.grid_params = params[key]
                elif key == 'integrator_params':
                    self.integrator_params = params[key]
                elif key == 'target_plates':
                    self.target_plates = params[key]

        # Store YAML data as dictionary object
            self.yaml = params

    def get_yaml(self):
        for item in self.yaml.keys():
            print('=' * 80)
            print('INGRID {}: {}'.format(item, self.yaml[item]))
            print('=' * 80 + '\n')

    def get_default_values(self, key):
        if key == 'eqdsk':
            return None
        elif key == 'grid_params':
            return self.default_grid_params
        elif key == 'integrator_params':
            return self.default_integrator_params
        elif key == 'target_plates':
            return self.default_target_plates_params
        else:
            print('Key not recognized... Add default values to source code for support.')

    def OMFIT_read_psi(self):
        """
        Python class to read the psi data in from an ascii file.
        Saves the boundary information and generated efit_data instance
        """

        g = OMFITgeqdsk(self.yaml['eqdsk'])

        nxefit = g['NW']
        nyefit = g['NH']
        rdim = g['RDIM']
        zdim = g['ZDIM']
        zmid = g['ZMID']
        rgrid1 = g['RLEFT']

        rcenter = g['RCENTR']
        bcenter = g['BCENTR']

        psi = g['PSIRZ'].T

        # TODO: possibly use the limiters to determine where the strke plates
        # are located. rlim and zlim contain lists of the r and z coordinates
        # for the limiter surrounding the plasma
        # rlim = g['RLIM']  # limiter - similar to strike plates
        # zlim = g['ZLIM']


        # calc array for r and z
        rmin = rgrid1
        rmax = rmin + rdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        # reproduce efit grid
        self.efit_psi = Efit_Data(rmin, rmax, nxefit,
                                  zmin, zmax, nyefit,
                                  rcenter, bcenter, name='Efit Data')
        self.efit_psi.set_v(psi)

    def read_target_plates(self):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        """

        target_plates = self.yaml['target_plates']
        plate_data = {}

        for plate in target_plates:

            plate_data[plate] = []

            if not 'name' in target_plates[plate].keys():
                target_plates[plate].update({'name' : plate})
            
            elif target_plates[plate]['name'] == '':
                target_plates[plate]['name'] = plate

            try:
                with open(target_plates[plate]['file']) as f:
                    for line in f:
                        point = line.strip()
                        if point.startswith('#'):
                            # move along to the next iteration
                            # this simulates a comment
                            continue
                        x = float(point.split(',')[0])
                        y = float(point.split(',')[1])
                        plate_data[plate].append((x, y))
                print('Using target plate "{}": {}'.format(target_plates[plate]['name'], plate_data[plate]))
            except KeyError:
                print('No inner target plate file.')
            except FileNotFoundError:
                pass

        self.plate_data = plate_data

    def plot_target_plates(self):
        """ Plots all target plates on the current figure """
        
        try:
            for plate in self.plate_data:
                if self.plate_data[plate]:
                    coor = np.array(self.plate_data[plate])
                    plt.plot(coor[:, 0], coor[:, 1], label=plate)
                    plt.draw()
        except AttributeError:
            pass

    def calc_efit_derivs(self):
        """ Calculates the partial derivatives using finite differences.
        Wrapper for the member function in the efit class.
        """
        self.efit_psi.Calculate_PDeriv()

    def plot_efit_data(self):
        """ Generates the plot that we will be able to manipulate
        using the root finder """
        self.efit_psi.clear_plot()
        self.efit_psi.plot_data()
    
    def find_roots(self, tk_controller = None):
        """ Displays a plot, and has the user click on an approximate
        zero point. Uses a root finder to adjust to the more exact point.
        Right click to disable.
        """
        # All we need to do is pass in the crude grid with its derivatives.
        # The interpolation for derivatives is handled by our newton functon
        # make plot - click on point - refine null point
        # the roots are saved at self.root_finder.roots

        self.root_finder = RootFinder(self.efit_psi, controller = tk_controller)

    def toggle_root_finder(self):
        """ Activates or deactivates the root finder ability. Enables
        the user to save the location where they last clicked.
        """
        self.root_finder.toggle_root_finding()

    def calc_psinorm(self):
        """ Uses magx and xpt1 to normalize the psi data. Furthur calculations
        will use this information """
        from Interpol.Setup_Grid_Data import Efit_Data

        # use the same bounds as the efit data
        self.psi_norm = Efit_Data(self.efit_psi.rmin, self.efit_psi.rmax,
                                  self.efit_psi.nr, self.efit_psi.zmin,
                                  self.efit_psi.zmax, self.efit_psi.nz,
                                  self.efit_psi.rcenter, self.efit_psi.bcenter,
                                  name='psi norm')
        psi = self.efit_psi.v
        psi_magx = self.efit_psi.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.efit_psi.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx))/(psi_xpt1 - psi_magx)
        self.psi_norm.set_v(psinorm)
        self.psi_norm.Calculate_PDeriv()
        self.psi_norm.plot_data()

    def init_LineTracing(self):
        """ 
        Initializes the line tracing class for the construction
        of the grid.
        """
        try:
            # Check if the LineTracing class has been initialized.
            # Exception will be thrown if it no initialization has occured.
            if not self.eq:
                raise AttributeError

        except AttributeError:
            self.eq = LineTracing(self.psi_norm, self.yaml)

        self.eq.disconnect()

    def _get_configuration(self):
        return self.eq.config

    def _classify_gridtype(self):
        if self.yaml['grid_params']['num_xpt'] == 1:
            self.eq.SNL_find_NSEW(self.xpt1, self.magx)
        elif self.yaml['grid_params']['num_xpt'] == 2:
            print('Double null configurations not yet supported...')

        self.yaml['grid_params']['config'] = self.eq.config

        return self.yaml['grid_params']['config']

    def _analyze_topology(self):
        
        config = self._classify_gridtype()
        if config == 'LSN':
            ingrid_topology = LSN(self)
        elif config == 'USN':
            ingrid_topology = USN(self)

        self.current_topology = ingrid_topology

class SNL(Ingrid):

    def __init__(self, INGRID_object):
        super().__init__(params = INGRID_object.yaml)
        self.efit_psi = INGRID_object.efit_psi
        self.psi_norm = INGRID_object.psi_norm
        self.eq = INGRID_object.eq

        self.plate_data = INGRID_object.plate_data

    def construct_patches(self):
        print('SNL patches')

    def grid_diagram(self):
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
          'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
          'seagreen', 'firebrick', 'saddlebrown']

        try:
            plt.close('INGRID: Subgrid Diagram')
        except:
            pass

        plt.figure('INGRID: Subgrid Diagram', figsize=(6,10))
        for patch in self.patches:
            patch.plot_subgrid()
            print('patch completed...')
        
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('INGRID SNL Subgrid')
        plt.show()

    def patch_diagram(self):
        """ Generates the patch diagram for a given configuration. """
        
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']

        try:
            plt.close('INGRID: Patch Diagram')
        except:
            pass

        plt.figure('INGRID: Patch Diagram', figsize=(6, 10))
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('SNL Patch Diagram')

        for i in range(len(self.patches)):
            self.patches[i].plot_border('green')
            self.patches[i].fill(colors[i])

        plt.show() 

    def get_configuration(self):
        """ 
        Returns a string indicating whether the 
        """
        return self.config

    def set_gridue(self):
        """
        Prepare the relevant arrays for writing to GRIDUE.
        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patch_lookup['IDL'].npol - 1
        ixpt2 = ixrb - self.patch_lookup['ODL'].npol + 1
        iyseparatrix1 = self.patch_lookup['IDL'].nrad - 1
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
        rb_prod = self.psi_norm.rcenter * self.psi_norm.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.psi_norm.get_psi(_r, _z)
                    _br = self.psi_norm.get_psi(_r, _z, tag = 'vz') / _r
                    _bz = self.psi_norm.get_psi(_r, _z, tag = 'vr') / _r
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

    def write_gridue(self):
        g = self.gridue_params
        status = uegrid.write_gridue(np.int32(g['ixpt1']), np.int32(g['ixpt2']), np.int32(g['iyseptrx1']),\
                (g['rm'].astype(np.double)), g['zm'].astype(np.double), g['psi'].astype(np.double), g['br'].astype(np.double), g['bz'].astype(np.double),\
                g['bpol'].astype(np.double), g['bphi'].astype(np.double), g['b'].astype(np.double))

    def export(self):
        """ Saves the grid as an ascii file """
        # TODO export the points the patches contain, but don't overlap
        # any points
        self.write_gridue()

    def animate_grid(self):

        try:
            plt.close('Test Plot')
        except:
            pass
        plt.figure('Test Plot', figsize=(6, 10))
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('SNL Patch Diagram')

        k = [1,2,4,3,1]

        for i in range(len(self.rm)):
            for j in range(len(self.rm[0])):
                plt.plot(self.rm[i][j][k], self.zm[i][j][k])
                plt.pause(0.05)

class LSN(SNL, Ingrid):

    def __init__(self, INGRID_object):
        super().__init__(INGRID_object)
        self.config = 'LSN'
        print('=' * 80)
        print('LSN Object!')
        print('=' * 80 + '\n')
    
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
    
        debug = False

        self.itp = self.plate_data['plate_W1']
        self.otp = self.plate_data['plate_E1']

        xpt = self.eq.eq_psi
        magx = np.array([self.grid_params['rmagx'], self.grid_params['zmagx']])
        psi_max = self.grid_params['psi_max']
        psi_min_core = self.grid_params['psi_min_core']
        psi_min_pf = self.grid_params['psi_min_pf']

        ITP = Line([p for p in [Point(i) for i in self.itp]])
        OTP = Line([p for p in [Point(i) for i in self.otp]])

        # Generate Horizontal Mid-Plane line
        LHS_Point = Point(magx[0] - 1e6, magx[1])
        RHS_Point = Point(magx[0] + 1e6, magx[1])
        midLine = Line([LHS_Point, RHS_Point])

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Drawing Separatrix
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = debug)
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        
        # Drawing Lower-SNL region
        xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = debug)        


        iPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : ITP}, option = 'theta', direction = 'ccw', show_plot = debug)
        oPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : OTP}, option = 'theta', direction = 'cw', show_plot = debug)

        imidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)
        omidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)

        xpt_ITP = self.eq.draw_line(xpt['SW'], {'line' : ITP}, option = 'theta', direction = 'ccw', show_plot = debug)
        xpt_OTP = self.eq.draw_line(xpt['SE'], {'line' : self.otp}, option = 'theta', direction = 'cw', show_plot = debug)
        
        # Integrating horizontally along mid-line towards psiMax and psiMinCore
        imidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : psi_max}, option = 'z_const', \
                direction = 'ccw', show_plot = debug)
        imidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', \
                direction = 'cw', show_plot = debug)
        omidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : psi_max}, option = 'z_const', \
                direction = 'cw', show_plot = debug)
        omidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', \
                direction = 'ccw', show_plot = debug)
        
        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'cw', show_plot = debug)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'ccw', show_plot = debug)

        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : ITP},option = 'theta', direction = 'ccw', show_plot = debug)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.otp}, option = 'theta', direction = 'cw', show_plot = debug)
        
        # IDL Patch
        IDL_N = iPsiMax_TP.reverse_copy()
        IDL_S = xpt_ITP
        IDL_E = xptW_psiMax.reverse_copy()
        IDL_W = ITP.reverse_copy()
        location = 'W'
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP
        IPF_E = xptS_psiMinPF
        IPF_W = ITP.reverse_copy()
        location = 'W'
        IPF = SNL_Patch([IPF_N, IPF_E, IPF_S, IPF_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # ISB Patch
        ISB_N = self.eq.draw_line(IDL_N.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
        ISB_S = xptNW_midLine.reverse_copy()
        ISB_E = Line([ISB_N.p[-1], ISB_S.p[0]])
        ISB_W = xptW_psiMax
        ISB = SNL_Patch([ISB_N, ISB_E, ISB_S, ISB_W], patchName = 'ISB')

        # ICB Patch
        ICB_N = ISB_S.reverse_copy()
        ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
        ICB_E = Line([ICB_N.p[-1], ICB_S.p[0]])
        ICB_W = xptN_psiMinCore.reverse_copy()
        ICB = SNL_Patch([ICB_N, ICB_E, ICB_S, ICB_W], patchName = 'ICB')

        # IST Patch
        IST_N = self.eq.draw_line(ISB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)
        IST_S = imidLine_topLine.reverse_copy()
        IST_E = Line([IST_N.p[-1], IST_S.p[0]])
        IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        IST = SNL_Patch([IST_N, IST_E, IST_S, IST_W], patchName = 'IST')

        # ICT Patch
        ICT_N = IST_S.reverse_copy()
        ICT_S = self.eq.draw_line(ICB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
        ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
        ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        ICT = SNL_Patch([ICT_N, ICT_E, ICT_S, ICT_W], patchName = 'ICT')

        # ODL Patch 
        ODL_N = oPsiMax_TP
        ODL_S = xpt_OTP.reverse_copy()
        ODL_E = OTP.reverse_copy()
        ODL_W = xptE_psiMax
        location = 'E'
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP.reverse_copy()
        OPF_E = OTP.reverse_copy()
        OPF_W = xptS_psiMinPF.reverse_copy()
        location = 'E'
        OPF = SNL_Patch([OPF_N, OPF_E, OPF_S, OPF_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # OSB Patch 
        OSB_N = self.eq.draw_line(ODL_N.p[0], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
        OSB_S = xptNE_midLine
        OSB_E = xptE_psiMax.reverse_copy()
        OSB_W = Line([OSB_S.p[-1], OSB_N.p[0]])
        OSB = SNL_Patch([OSB_N, OSB_E, OSB_S, OSB_W], patchName = 'OSB')

        # OCB Patch
        OCB_N = OSB_S.reverse_copy()
        OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        OCB_E = xptN_psiMinCore
        OCB_W = Line([OCB_S.p[-1], OCB_N.p[0]])
        OCB = SNL_Patch([OCB_N, OCB_E, OCB_S, OCB_W], patchName = 'OCB')

        # OST Patch
        OST_N = self.eq.draw_line(OSB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
        OST_S = omidLine_topLine
        OST_E = Line([OST_N.p[-1], OST_S.p[0]])
        OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        OST = SNL_Patch([OST_N, OST_E, OST_S, OST_W], patchName = 'OST')

        # OCT Patch
        OCT_N = OST_S.reverse_copy()
        OCT_S = self.eq.draw_line(OCB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        OCT_E = Line([OCT_N.p[-1], OCT_S.p[0]])
        OCT_W = Line([OCT_S.p[-1], OCT_N.p[0]])
        OCT = SNL_Patch([OCT_N, OCT_E, OCT_S, OCT_W], patchName = 'OCT')

        self.patches = [IDL, IPF, ISB, ICB, IST, ICT, OST, OCT, OSB, OCB, ODL, OPF]
        names = ['IDL', 'IPF', 'ISB', 'ICB', 'IST', 'ICT', 'OST', 'OCT', 'OSB', 'OCB', 'ODL', 'OPF']
        patch_lookup = {}
        for patch in self.patches:
            patch_lookup[patch.patchName] = patch
            patch.plot_border()
            patch.fill()
        self.patch_lookup = patch_lookup

        p = self.patch_lookup
        self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                        [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
                        [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
                        [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                        ]

    def construct_grid(self, np_cells = 1, nr_cells = 1):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.
        primary_xpt = Point([self.grid_params['rxpt'], self.grid_params['zxpt']])
        for patch in self.patches:
            patch.make_subgrid(self, np_cells, nr_cells)

            if patch.patchName == 'IDL':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'IPF':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'ISB':
                patch.adjust_corner(primary_xpt, 'SW')
            elif patch.patchName == 'ICB':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'OCB':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'OSB':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'OPF':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'ODL':
                patch.adjust_corner(primary_xpt, 'SW')

        self.concat_grid(self.get_configuration())
        self.set_gridue()

    def add_guardc(self, cell_map, ixlb, ixrb, config, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, config, boundary):

            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # Edit the code to reflect this at some point so the next reader is not overwhelmed.
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
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'top')

        return cell_map

    def concat_grid(self, config):
        """
        Concatenate all local grids on individual patches into a single 
        array with branch cuts

        Parameters:
        ----------
            config : str
                Type of SNL grid to concat.

        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).


        p = self.patch_lookup
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the 
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches:
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

        # Number of Poloidal patch indices with guard patches
        np_patch = len(patch_matrix[0])
        # Number of Radial patch indices with guard patches
        nr_patch = len(patch_matrix)

        # Total number of poloidal indices in subgrid.
        np_total = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:-1]])) + 2
        nr_total = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:3]])) + 2

        rm = np.zeros((np_total, nr_total, 5), order = 'F')
        zm = np.zeros((np_total, nr_total, 5), order = 'F')

        ixcell = 0
        jycell = 0

        pol_const = patch_matrix[1][1].npol - 1
        rad_const = patch_matrix[1][1].nrad - 1

        for ixp in range(1, np_patch - 1):
            for jyp in range(1, nr_patch - 1):

                for ixl in range(patch_matrix[jyp][ixp].npol - 1):
                    for jyl in range(patch_matrix[jyp][ixp].nrad - 1):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) - pol_const + ixl + 1
                        jycell = int(np.sum([patch.nrad - 1 for patch in patch_matrix[1][1:jyp+1]])) - rad_const + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].x
                            zm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb, config)
        self.zm = self.add_guardc(zm, ixlb, ixrb, config)

        self.animate_grid()

class USN(SNL, Ingrid):
    def __init__(self, INGRID_object):
        self.config = 'USN'
        print('=' * 80)
        print('USN Object!')
        print('=' * 80 + '\n')
        super().__init__(INGRID_object)
    
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
    


        debug = False

        self.itp = self.plate_data['plate_W1']
        self.otp = self.plate_data['plate_E1']

        xpt = self.eq.eq_psi
        magx = np.array([self.grid_params['rmagx'], self.grid_params['zmagx']])
        psi_max = self.grid_params['psi_max']
        psi_min_core = self.grid_params['psi_min_core']
        psi_min_pf = self.grid_params['psi_min_pf']

        ITP = Line([p for p in [Point(i) for i in self.itp]])
        OTP = Line([p for p in [Point(i) for i in self.otp]])

        # Generate Horizontal Mid-Plane line
        LHS_Point = Point(magx[0] - 1e6, magx[1])
        RHS_Point = Point(magx[0] + 1e6, magx[1])
        midLine = Line([LHS_Point, RHS_Point])

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Drawing Separatrix
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = debug)
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        
        # Drawing Lower-SNL region
        xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = debug)        

        iPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : self.itp}, option = 'theta', direction = 'cw', show_plot = debug)
        oPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = debug)

        imidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        omidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)

        xpt_ITP = self.eq.draw_line(xpt['SE'], {'line' : self.itp}, option = 'theta', direction = 'cw', show_plot = debug)
        xpt_OTP = self.eq.draw_line(xpt['SW'], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = debug)

        # Integrating horizontally along mid-line towards psiMax and psiMinCore
        omidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : psi_max}, option = 'z_const', \
                direction = 'cw', show_plot = debug)
        omidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', \
                direction = 'ccw', show_plot = debug)
        imidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : psi_max}, option = 'z_const', \
                direction = 'ccw', show_plot = debug)
        imidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', \
                direction = 'cw', show_plot = debug)
        
        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'ccw', show_plot = debug)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'cw', show_plot = debug)
        
        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.itp},option = 'theta', direction = 'cw', show_plot = debug)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = debug)
        
        # IDL Patch
        IDL_N = iPsiMax_TP
        IDL_S = xpt_ITP.reverse_copy()
        IDL_E = Line([IDL_N.p[-1], IDL_S.p[0]])
        IDL_W = xptE_psiMax
        location = 'E'
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP.reverse_copy()
        IPF_E = Line([IPF_N.p[-1], IPF_S.p[0]])
        IPF_W = xptS_psiMinPF.reverse_copy()
        location = 'E'
        IPF = SNL_Patch([IPF_N, IPF_E, IPF_S, IPF_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # ISB Patch
        ISB_N = self.eq.draw_line(IDL_N.p[0], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
        ISB_S = xptNE_midLine
        ISB_E = xptE_psiMax.reverse_copy()
        ISB_W = Line([ISB_S.p[-1], ISB_N.p[0]])
        ISB = SNL_Patch([ISB_N, ISB_E, ISB_S, ISB_W], patchName = 'ISB')

        # ICB Patch
        ICB_N = ISB_S.reverse_copy()
        ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        ICB_E = xptN_psiMinCore
        ICB_W = Line([ICB_S.p[-1], ICB_N.p[0]])
        ICB = SNL_Patch([ICB_N, ICB_E, ICB_S, ICB_W], patchName = 'ICB')

        # IST Patch
        IST_N = self.eq.draw_line(ISB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
        IST_S = imidLine_topLine
        IST_E = Line([IST_N.p[-1], IST_S.p[0]])
        IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        IST = SNL_Patch([IST_N, IST_E, IST_S, IST_W], patchName = 'IST')

        # ICT Patch
        ICT_N = IST_S.reverse_copy()
        ICT_S = self.eq.draw_line(ICB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
        ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        ICT = SNL_Patch([ICT_N, ICT_E, ICT_S, ICT_W], patchName = 'ICT')

        # ODL Patch
        ODL_N = oPsiMax_TP.reverse_copy()
        ODL_S = xpt_OTP
        ODL_E = xptW_psiMax.reverse_copy()
        ODL_W = Line([ODL_S.p[-1], ODL_N.p[0]])
        location = 'W'
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP
        OPF_E = xptS_psiMinPF
        OPF_W = Line([OPF_S.p[-1], OPF_N.p[0]])
        location = 'W'
        OPF = SNL_Patch([OPF_N, OPF_E, OPF_S, OPF_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # OSB Patch
        OSB_N = self.eq.draw_line(ODL_N.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
        OSB_S = xptNW_midLine.reverse_copy()
        OSB_E = Line([OSB_N.p[-1], OSB_S.p[0]])
        OSB_W = xptW_psiMax
        OSB = SNL_Patch([OSB_N, OSB_E, OSB_S, OSB_W], patchName = 'OSB')

        # OCB Patch
        OCB_N = OSB_S.reverse_copy()
        OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
        OCB_E = Line([OCB_N.p[-1], OCB_S.p[0]])
        OCB_W = xptN_psiMinCore.reverse_copy()
        OCB = SNL_Patch([OCB_N, OCB_E, OCB_S, OCB_W], patchName = 'OCB')

        # OST Patch
        OST_N = self.eq.draw_line(OSB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)
        OST_S = omidLine_topLine.reverse_copy()
        OST_E = Line([OST_N.p[-1], OST_S.p[0]])
        OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        OST = SNL_Patch([OST_N, OST_E, OST_S, OST_W], patchName = 'OST')

        # OCT Patch
        OCT_N = OST_S.reverse_copy()
        OCT_S = self.eq.draw_line(OCB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
        OCT_E = Line([OCT_N.p[-1], OCT_S.p[0]])
        OCT_W = Line([OCT_S.p[-1], OCT_N.p[0]])
        OCT = SNL_Patch([OCT_N, OCT_E, OCT_S, OCT_W], patchName = 'OCT')

        self.patches = [OPF, ODL, OCB, OSB, OCT, OST, ICT, IST, ICB, ISB, IPF, IDL]
        patch_lookup = {}
        for patch in self.patches:
            patch_lookup[patch.patchName] = patch
            patch.plot_border()
            patch.fill()
        self.patch_lookup = patch_lookup

        p = self.patch_lookup
        self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                        [[None], p['ODL'], p['OSB'], p['OST'], p['IST'], p['ISB'], p['IDL'], [None]], \
                        [[None], p['OPF'], p['OCB'], p['OCT'], p['ICT'], p['ICB'], p['IPF'], [None]], \
                        [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                        ]

    def construct_grid(self, np_cells = 1, nr_cells = 1):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.

        primary_xpt = Point([self.grid_params['rxpt'], self.grid_params['zxpt']])
        for patch in self.patches:
            print(patch.patchName)
            patch.make_subgrid(self, np_cells, nr_cells)

            if patch.patchName == 'IDL':
                patch.adjust_corner(primary_xpt, 'SW')
            elif patch.patchName == 'IPF':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'ISB':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'ICB':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'OCB':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'OSB':
                patch.adjust_corner(primary_xpt, 'SW')
            elif patch.patchName == 'OPF':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'ODL':
                patch.adjust_corner(primary_xpt, 'SE')

        self.concat_grid(self.get_configuration())
        self.set_gridue()


    def add_guardc(self, cell_map, ixlb, ixrb, config, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, config, boundary):

            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # Edit the code to reflect this at some point so the next reader is not overwhelmed.
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
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'top')

        return cell_map

    def _add_guardc(self, cell_map, ixlb, ixrb, config, nxpt = 1, eps = 1e-3):

        def _set_guard(cell_map, ix, iy, eps, config, boundary):

            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # Edit the code to reflect this at some point so the next reader is not overwhelmed.
            if boundary == 'left':
                ixn = ix - 1
                iyn = iy
                cell_map[ix][iy][2] = cell_map[ixn][iyn][1]
                cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][2])
                cell_map[ix][iy][4] = cell_map[ixn][iyn][3]
                cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][4])
                cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

            elif boundary == 'right':
                ixn = ix + 1
                iyn = iy
                cell_map[ix][iy][1] = cell_map[ixn][iyn][2]
                cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][1])
                cell_map[ix][iy][3] = cell_map[ixn][iyn][4]
                cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][3])
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
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'top')

        return cell_map    

    def concat_grid(self, config):
        """
        Concatenate all local grids on individual patches into a single 
        array with branch cuts

        Parameters:
        ----------
            config : str
                Type of SNL grid to concat.

        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).


        p = self.patch_lookup
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the 
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches:
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

        # Number of Poloidal patch indices with guard patches
        np_patch = len(patch_matrix[0])
        # Number of Radial patch indices with guard patches
        nr_patch = len(patch_matrix)

        # Total number of poloidal indices in subgrid.
        np_total = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:-1]])) + 2
        nr_total = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:3]])) + 2

        rm = np.zeros((np_total, nr_total, 5), order = 'F')
        zm = np.zeros((np_total, nr_total, 5), order = 'F')

        ixcell = 0
        jycell = 0

        pol_const = patch_matrix[1][1].npol - 1
        rad_const = patch_matrix[1][1].nrad - 1

        for ixp in range(1, np_patch - 1):
            for jyp in range(1, nr_patch - 1):

                for ixl in range(patch_matrix[jyp][ixp].npol - 1):
                    for jyl in range(patch_matrix[jyp][ixp].nrad - 1):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) - pol_const + ixl + 1
                        jycell = int(np.sum([patch.nrad - 1 for patch in patch_matrix[1][1:jyp+1]])) - rad_const + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].x
                            zm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb, config)
        self.zm = self.add_guardc(zm, ixlb, ixrb, config)

        self.animate_grid()

    def _concat_grid(self, config):
        """
        Concatenate all local grids on individual patches into a single 
        array with branch cuts

        Parameters:
        ----------
            config : str
                Type of SNL grid to concat.

        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).


        p = self.patch_lookup
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the 
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches:
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

        # Number of Poloidal patch indices with guard patches
        np_patch = len(patch_matrix[0])
        # Number of Radial patch indices with guard patches
        nr_patch = len(patch_matrix)

        # Total number of poloidal indices in subgrid.
        np_total = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:-1]])) + 2
        nr_total = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:3]])) + 2

        rm = np.zeros((np_total, nr_total, 5), order = 'F')
        zm = np.zeros((np_total, nr_total, 5), order = 'F')

        ixcell = 0
        jycell = 0

        pol_const = patch_matrix[1][1].npol - 1
        rad_const = patch_matrix[1][1].nrad - 1

        for ixp in range(1, np_patch - 1):
            for jyp in range(1, nr_patch - 1):

                for ixl in range(patch_matrix[jyp][ixp].npol - 1):
                    for jyl in range(patch_matrix[jyp][ixp].nrad - 1):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) - ixl
                        jycell = int(np.sum([patch.nrad - 1 for patch in patch_matrix[1][1:jyp+1]])) - rad_const + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].x
                            zm[ixcell][nr_total - jycell - 1][ind] = patch_matrix[jyp][ixp].cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb, config)
        self.zm = self.add_guardc(zm, ixlb, ixrb, config)
        

        self.animate_grid()


def set_params_GUI():
    """
    
    GUI version for Ingrid. Gets the files used, and opens a plot
    of the Efit data. Allows user to input values for magx, xpt, 
    and psi levels.

    """
    try:
        import tkinter as tk
    except:
        import Tkinter as tk
    try:
        import tkinter.messagebox as messagebox
    except:
        import tkMessageBox as messagebox 
    def on_closing():
        if messagebox.askyesno('', 'Are you sure you want to quit?'):
            plt.close('all')
            IngridWindow.destroy()
    IngridWindow = IA.IngridApp()
    IngridWindow.geometry("830x490")
    IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
    IngridWindow.mainloop()
