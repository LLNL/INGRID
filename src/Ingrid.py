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
import f90nml
import IngridApp as IA


class Ingrid:
    """ An interactive grid generator for edge plasmas in a tokamak
    Accepts a dictionary generated from a namelist file that contains
    the parameters.
    
    Parameters
    ----------
    nml : dict
        Params dictionary object contains two dictionaries
        First is 'files' which contains the keys: geqdsk, itp, otp.
        The second is 'grid params' which has: psi_max, psi_min_core,
        psi_min_pf, Rmagx, Zmagx, Rxpt, Zxpt
    
    """

    def __init__(self, nml = { 'files' : {}, 'grid_params' : {} } ):
        self.files = nml['files']
        self.grid_params = nml['grid_params']

        print('Welcome to Ingrid!\n')
        
    def setup(self):
        """ Add the magnetic axis and the x-point """
        self.add_magx(self.grid_params['rmagx'], self.grid_params['zmagx'])
        self.add_xpt1(self.grid_params['rxpt'], self.grid_params['zxpt'])

    # TODO update this method to be able to change parameters.
    def set_param(self, key=None, value=None):
        """ Sets the parameters that we want, without
        adding new and unwanted options.
        User can pass in new values directly in,
        or can be prompted.
        
        Parameters
        ----------
        key : str, optional
            Keyword for the dictionary object. Specify which value we
            wish to change. Will be prompted for it if is None.
        value : varies, optional
            Value associated with a key. Must match the type of the
            value currently tied to the key.
        
        
        """
        def check_key(self, key=None, first_time=True):
            while True:
                try:
                    if not first_time:
                        # defaults to a string
                        key = raw_input('Enter Key: ')
                    if key not in self._grid_params:
                        raise ValueError
                    print('Key: "{}" is in Params'.format(key))
                    break

                except ValueError:
                    print('That key is not in Grid Parameters.')
                    first_time = False
            return key

        def check_value(self, key, value=None, first_time=True):
            while True:
                try:
                    if not first_time or value is None:
                        # Need to preserve the type of the value entered
                        value = input('Enter Value: ')
                    if not isinstance(self._grid_params[key],
                                      type(value)):
                        raise TypeError
                    self._grid_params[key] = value
                    print('New value recorded as "{}".'.format(value))
                    break

                except TypeError:
                    print('That value is the wrong type.')
                    first_time = False

        # in case the user didn't enter the key
        if key is None:
            key = check_key(first_time=False)
            check_value(key, first_time=False)

        # check if the key the user entered is in grid params
        else:
            key = check_key(key)
            check_value(key, value)

    # TODO: update this method.
    def get_param(self, key=None):
        """ Returns the value associated with a key. """
        first_time = True
        while True:
            try:
                if not first_time or key is None:
                    key = raw_input('Enter Key: ')
                if key not in self._grid_params:
                    raise ValueError
                return self._grid_params[key]

            except ValueError:
                print('Key not a grid parameter.')
                first_time = False

    def print_params(self):
        """ Lists the parameters that the user should know about. """
        for key, value in self._grid_params.items():
            print(key, value)


    def OMFIT_read_psi(self):
        """
        Python class to read the psi data in from an ascii file.
        Saves the boundary information and generated efit_data instance
        """
        
        from OMFITgeqdsk import OMFITgeqdsk
        from Interpol.Setup_Grid_Data import Efit_Data

        g = OMFITgeqdsk(self.files['geqdsk'])

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


    def import_psi_data(self):
        """ Same as OMFIT_read_psi, but calls a fortran function to
        read in the psi data. Saves and creates the same instance.
        """
        import Efit.efit as efit
        from Interpol.Setup_Grid_Data import Efit_Data
        efit.readg(self.files['geqdsk'])  # defines the variables
        # extract dimensions and boundaries
        nxefit, nyefit, nbdry, nlim = efit.get_nxy()

        fold, fpol, pres, qpsi, \
            rbdry, zbdry, xlim, ylim, xdim, zdim, \
            rcentr, rgrid1, zmid, Rmagx, Zmagx, \
            simagx, sibdry, bcentr = efit.get_psi(nxefit, nyefit,
                                                  nbdry, nlim)
        # calc array for r and z
        rmin = rgrid1
        rmax = rmin + xdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        # reproduce efit grid
        self.efit_psi = Efit_Data(rmin, rmax, nxefit,
                                  zmin, zmax, nyefit,
                                  rcentr, bcentr,
                                  name='Efit Data')
        self.efit_psi.set_v(fold)  # this is psi


    def read_target_plate(self):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        """

        self.itp = []
        try:
            with open(self.files['itp']) as f:
                for line in f:
                    point = line.strip()
                    if point.startswith('#'):
                        # move along to the next iteration
                        # this simulates a comment
                        continue
                    x = float(point.split(',')[0])
                    y = float(point.split(',')[1])
                    self.itp.append((x, y))
            print('Using inner target plate', self.itp)
        except KeyError:
            print('No inner target plate file.')

        self.otp = []
        try:
            with open(self.files['otp']) as f:
                for line in f:
                    point = line.strip()
                    if point.startswith('#'):
                        # move along to the next iteration
                        continue
                    x = float(point.split(',')[0])
                    y = float(point.split(',')[1])
                    self.otp.append((x, y))
            print('Using outer target plate', self.otp)
        except KeyError:
            print('No outer target plate file.')


    def plot_target_plate(self):
        """ Plots the inner and outer target plates on the current figure """
        
        itp = np.array(self.itp)
        otp = np.array(self.otp)

        plt.plot(itp[:, 0], itp[:, 1], label='itp')
        plt.plot(otp[:, 0], otp[:, 1], label='otp')
        plt.draw()

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
        from Root_Finder import RootFinder
        self.root_finder = RootFinder(self.efit_psi, controller = tk_controller)

    def save_root(self):
        """ Returns the value of the root finder. This will either be the
        root, or the click, depending on the root finding settings.
        """
        return self.root_finder.final_root

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


    def compute_eq_psi(self):
        """ Initializes the line tracing class for the construction
        of the grid. Necessary to be separate so we can wiat until the
        user clicks on the root.
        """
        from line_tracing import LineTracing
        self.eq = LineTracing(self.psi_norm, self.grid_params,
                              option='xpt_circ', eps = 1e-3)
        self.eq.find_NSEW(self.xpt1, self.magx)
        self.eq.disconnect()

    def get_SNL_configuration(self):
        """ 
        Returns a string indicating whether the 
        """
        return self.eq.SNL_CONFIG

    def construct_SNL_patches(self):
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
    
        from geometry import Point, Patch, Line

        debug = False

        # Get starting directions from primary x-point
        self.compute_eq_psi()
        # Sign check to determine if Upper Single Null or Lower Single Null
        # These coordinates tell us which quadrant in a 2D cartesian plane
        # our 'N' unit vector points towards.
        SNL_CONFIG = self.get_SNL_configuration()

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
        midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])
        topLine.plot()

        # Drawing Separatrix
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = debug)
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
        
        # Drawing Lower-SNL region
        xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = debug)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = debug)        

        if SNL_CONFIG == 'USN':
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

        elif SNL_CONFIG == 'LSN':
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
        if SNL_CONFIG == 'USN':
            IDL_N = iPsiMax_TP
            IDL_S = xpt_ITP.reverse_copy()
            IDL_E = Line([IDL_N.p[-1], IDL_S.p[0]])
            IDL_W = xptE_psiMax
            location = 'E'
        elif SNL_CONFIG == 'LSN':
            IDL_N = iPsiMax_TP.reverse_copy()
            IDL_S = xpt_ITP
            IDL_E = xptW_psiMax.reverse_copy()
            IDL_W = ITP.reverse_copy()
            location = 'W'
        IDL = Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        if SNL_CONFIG == 'USN':
            IPF_N = IDL_S.reverse_copy()
            IPF_S = psiMinPF_ITP.reverse_copy()
            IPF_E = Line([IPF_N.p[-1], IPF_S.p[0]])
            IPF_W = xptS_psiMinPF.reverse_copy()
            location = 'E'
        elif SNL_CONFIG == 'LSN':
            IPF_N = IDL_S.reverse_copy()
            IPF_S = psiMinPF_ITP
            IPF_E = xptS_psiMinPF
            IPF_W = ITP.reverse_copy()
            location = 'W'
        IPF = Patch([IPF_N, IPF_E, IPF_S, IPF_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # ISB Patch
        if SNL_CONFIG == 'USN':
            ISB_N = self.eq.draw_line(IDL_N.p[0], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
            ISB_S = xptNE_midLine
            ISB_E = xptE_psiMax.reverse_copy()
            ISB_W = Line([ISB_S.p[-1], ISB_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            ISB_N = self.eq.draw_line(IDL_N.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
            ISB_S = xptNW_midLine.reverse_copy()
            ISB_E = Line([ISB_N.p[-1], ISB_S.p[0]])
            ISB_W = xptW_psiMax
        ISB = Patch([ISB_N, ISB_E, ISB_S, ISB_W], patchName = 'ISB')

        # ICB Patch
        if SNL_CONFIG == 'USN':
            ICB_N = ISB_S.reverse_copy()
            ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
            ICB_E = xptN_psiMinCore
            ICB_W = Line([ICB_S.p[-1], ICB_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            ICB_N = ISB_S.reverse_copy()
            ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
            ICB_E = Line([ICB_N.p[-1], ICB_S.p[0]])
            ICB_W = xptN_psiMinCore.reverse_copy()
        ICB = Patch([ICB_N, ICB_E, ICB_S, ICB_W], patchName = 'ICB')

        # IST Patch
        if SNL_CONFIG == 'USN':
            IST_N = self.eq.draw_line(ISB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
            IST_S = imidLine_topLine
            IST_E = Line([IST_N.p[-1], IST_S.p[0]])
            IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            IST_N = self.eq.draw_line(ISB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)
            IST_S = imidLine_topLine.reverse_copy()
            IST_E = Line([IST_N.p[-1], IST_S.p[0]])
            IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        IST = Patch([IST_N, IST_E, IST_S, IST_W], patchName = 'IST')

        # ICT Patch
        if SNL_CONFIG == 'USN':
            ICT_N = IST_S.reverse_copy()
            ICT_S = self.eq.draw_line(ICB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)
            ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
            ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            ICT_N = IST_S.reverse_copy()
            ICT_S = self.eq.draw_line(ICB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
            ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
            ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        ICT = Patch([ICT_N, ICT_E, ICT_S, ICT_W], patchName = 'ICT')

        # ODL Patch
        if SNL_CONFIG == 'USN':
            ODL_N = oPsiMax_TP.reverse_copy()
            ODL_S = xpt_OTP
            ODL_E = xptW_psiMax.reverse_copy()
            ODL_W = Line([ODL_S.p[-1], ODL_N.p[0]])
            location = 'W'
        elif SNL_CONFIG == 'LSN':
            ODL_N = oPsiMax_TP
            ODL_S = xpt_OTP.reverse_copy()
            ODL_E = OTP.reverse_copy()
            ODL_W = xptE_psiMax
            location = 'E'
        ODL = Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        if SNL_CONFIG == 'USN':
            OPF_N = ODL_S.reverse_copy()
            OPF_S = psiMinPF_OTP
            OPF_E = xptS_psiMinPF
            OPF_W = Line([OPF_S.p[-1], OPF_N.p[0]])
            location = 'W'
        elif SNL_CONFIG == 'LSN':
            OPF_N = ODL_S.reverse_copy()
            OPF_S = psiMinPF_OTP.reverse_copy()
            OPF_E = OTP.reverse_copy()
            OPF_W = xptS_psiMinPF.reverse_copy()
            location = 'E'
        OPF = Patch([OPF_N, OPF_E, OPF_S, OPF_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # OSB Patch
        if SNL_CONFIG == 'USN':
            OSB_N = self.eq.draw_line(ODL_N.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug)
            OSB_S = xptNW_midLine.reverse_copy()
            OSB_E = Line([OSB_N.p[-1], OSB_S.p[0]])
            OSB_W = xptW_psiMax
        elif SNL_CONFIG == 'LSN':
            OSB_N = self.eq.draw_line(ODL_N.p[0], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
            OSB_S = xptNE_midLine
            OSB_E = xptE_psiMax.reverse_copy()
            OSB_W = Line([OSB_S.p[-1], OSB_N.p[0]])
        OSB = Patch([OSB_N, OSB_E, OSB_S, OSB_W], patchName = 'OSB')

        # OCB Patch
        if SNL_CONFIG == 'USN':
            OCB_N = OSB_S.reverse_copy()
            OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
            OCB_E = Line([OCB_N.p[-1], OCB_S.p[0]])
            OCB_W = xptN_psiMinCore.reverse_copy()
        elif SNL_CONFIG == 'LSN':
            OCB_N = OSB_S.reverse_copy()
            OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : midLine}, option = 'theta', direction = 'ccw', show_plot = debug)
            OCB_E = xptN_psiMinCore
            OCB_W = Line([OCB_S.p[-1], OCB_N.p[0]])
        OCB = Patch([OCB_N, OCB_E, OCB_S, OCB_W], patchName = 'OCB')

        # OST Patch
        if SNL_CONFIG == 'USN':
            OST_N = self.eq.draw_line(OSB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug)
            OST_S = omidLine_topLine.reverse_copy()
            OST_E = Line([OST_N.p[-1], OST_S.p[0]])
            OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            OST_N = self.eq.draw_line(OSB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug).reverse_copy()
            OST_S = omidLine_topLine
            OST_E = Line([OST_N.p[-1], OST_S.p[0]])
            OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        OST = Patch([OST_N, OST_E, OST_S, OST_W], patchName = 'OST')

        # OCT Patch
        if SNL_CONFIG == 'USN':
            OCT_N = OST_S.reverse_copy()
            OCT_S = self.eq.draw_line(OCB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = debug).reverse_copy()
            OCT_E = Line([OCT_N.p[-1], OCT_S.p[0]])
            OCT_W = Line([OCT_S.p[-1], OCT_N.p[0]])
        elif SNL_CONFIG == 'LSN':
            OCT_N = OST_S.reverse_copy()
            OCT_S = self.eq.draw_line(OCB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = debug)
            OCT_E = Line([OCT_N.p[-1], OCT_S.p[0]])
            OCT_W = Line([OCT_S.p[-1], OCT_N.p[0]])
        OCT = Patch([OCT_N, OCT_E, OCT_S, OCT_W], patchName = 'OCT')

        self.patches = [IDL, IPF, ISB, ICB, IST, ICT, OST, OCT, OSB, OCB, ODL, OPF]
        names = ['IDL', 'IPF', 'ISB', 'ICB', 'IST', 'ICT', 'OST', 'OCT', 'OSB', 'OCB', 'ODL', 'OPF']
        patch_lookup = {}
        for i in range(len(names)):
            patch_lookup[names[i]] = self.patches[i]
        self.patch_lookup = patch_lookup


        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.
        from timeit import default_timer as timer
        subgrid_start = timer()

        primary_xpt = Point([self.grid_params['rxpt'], self.grid_params['zxpt']])

        i = 0
        for patch in self.patches:
            print(names[i])
            i += 1
            patch.make_subgrid(self, num = 4)
            patch.plot_border()
            patch.fill()
            #TODO: Make this it's own function? It's a bit cumbersome looking...
            if SNL_CONFIG == 'USN':
                if patch is IDL:
                    patch.adjust_corner(primary_xpt, 'SW')
                elif patch is IPF:
                    patch.adjust_corner(primary_xpt, 'NW')
                elif patch is ISB:
                    patch.adjust_corner(primary_xpt, 'SE')
                elif patch is ICB:
                    patch.adjust_corner(primary_xpt, 'NE')
                elif patch is OCB:
                    patch.adjust_corner(primary_xpt, 'NW')
                elif patch is OSB:
                    patch.adjust_corner(primary_xpt, 'SW')
                elif patch is OPF:
                    patch.adjust_corner(primary_xpt, 'NE')
                elif patch is ODL:
                    patch.adjust_corner(primary_xpt, 'SE')
            elif SNL_CONFIG == 'LSN':
                if patch is IDL:
                    patch.adjust_corner(primary_xpt, 'SE')
                elif patch is IPF:
                    patch.adjust_corner(primary_xpt, 'NE')
                elif patch is ISB:
                    patch.adjust_corner(primary_xpt, 'SW')
                elif patch is ICB:
                    patch.adjust_corner(primary_xpt, 'NW')
                elif patch is OCB:
                    patch.adjust_corner(primary_xpt, 'NE')
                elif patch is OSB:
                    patch.adjust_corner(primary_xpt, 'SE')
                elif patch is OPF:
                    patch.adjust_corner(primary_xpt, 'NW')
                elif patch is ODL:
                    patch.adjust_corner(primary_xpt, 'SW')
        subgrid_end = timer()
        print('Subgrid Generation took: {}s'.format(subgrid_end - subgrid_start))

    def grid_diagram(self):
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
          'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
          'seagreen', 'firebrick', 'saddlebrown']
        plt.figure('INGRID', figsize=(6,10))
        for patch in self.patches:
            patch.plot_subgrid()
            print('patch completed...')
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('INGRID Subgrid')
        plt.show()

    def grid_diagram_debug():
        plt.figure('DEBUG')
        import pdb
        pdb.set_trace()

    def patch_diagram(self):
        """ Generates the patch diagram for a given configuration. """
        
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']

        plt.figure('INGRID', figsize=(6, 10))
        for i in range(len(self.patches)):
            self.patches[i].plot_border('green')
            self.patches[i].fill(colors[i])
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('Example Patch')
        plt.show() 


    def concat_grid(self, config):
        """
        Concatenate all local grids on individual patches into a single 
        array with branch cuts

        Parameters:
        ----------
            config : str
                Type of grid we are dealing with.
                Default: 'SNL'

        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).



        p = self.patch_lookup
        patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                        [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
                        [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
                        [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                        ]

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

        if config == 'LSN':
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

        elif config == 'USN':
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

    def add_guardc(self, cell_map, ixlb, ixrb, config, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, config, boundary):

            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # Edit the code to reflect this at some point so the next reader is not overwhelmed.
            if boundary == 'left':

                if config == 'LSN':
                    ixn = ix + 1
                    iyn = iy
                    cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][2])
                    cell_map[ix][iy][2] = cell_map[ixn][iyn][1]
                    cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][4])
                    cell_map[ix][iy][4] = cell_map[ixn][iyn][3]
                    cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])
                elif config == 'USN':
                    ixn = ix - 1
                    iyn = iy
                    cell_map[ix][iy][2] = cell_map[ixn][iyn][1]
                    cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][2])
                    cell_map[ix][iy][4] = cell_map[ixn][iyn][3]
                    cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][4])
                    cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

            elif boundary == 'right':

                if config == 'LSN':
                    ixn = ix - 1
                    iyn = iy
                    cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][1])
                    cell_map[ix][iy][1] = cell_map[ixn][iyn][2]
                    cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][3])
                    cell_map[ix][iy][3] = cell_map[ixn][iyn][4]
                    cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])
                elif config == 'USN':
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
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'left' if config == 'LSN' else 'right')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'right' if config == 'LSN' else 'left')

            if nxpt > 1:
                print('Only SNL is currently supported...')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, config = config, boundary = 'top')

        return cell_map

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

        self.gridue = {'nxm' : nxm, 'nym' : nym, 'ixpt1' : ixpt1, 'ixpt2' : ixpt2, 'iyseptrx1' : iyseparatrix1, \
            'rm' : self.rm, 'zm' : self.zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b}

    def write_gridue(self):
        import Uegrid.uegrid as uegrid
        g = self.gridue
        status = uegrid.write_gridue(np.int32(g['ixpt1']), np.int32(g['ixpt2']), np.int32(g['iyseptrx1']),\
                (g['rm'].astype(np.double)), g['zm'].astype(np.double), g['psi'].astype(np.double), g['br'].astype(np.double), g['bz'].astype(np.double),\
                g['bpol'].astype(np.double), g['bphi'].astype(np.double), g['b'].astype(np.double))

    def export(self):
        """ Saves the grid as an ascii file """
        # TODO export the points the patches contain, but don't overlap
        # any points
        self.concat_grid(self.get_SNL_configuration())
        self.set_gridue()
        self.write_gridue()


    def test_interpol(self, option=2, nfine=100, ncrude=10, tag='v'):
        """ Provides a demonstration and test of the bicubic interpolation
        methods used. Wrapper for the real code from another module.
        
        Parameters
        ----------
        option : int, optional
            Specify which of a set of function to use as the test.
            Accepts 1, 11, 12, 13, 2, 3, 4, 5, 51, 52
            See Test_Functions.get_f for more detail.
        nfine : int, optional
            Density of the fine grid we will interpolate onto.
        ncrude : int, optional
            Density of the crude grid the sample data will be generated
            onto.
        tag : str, optional
            Specify if it is wanted to test the derivative interpolation
            methods. Accepts 'v', 'vr', 'vz', 'vrz'.
        
        """
        from Interpol.Test_Interpol import test_interpol
        test_interpol(option, nfine, ncrude, tag)

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
    IngridWindow.geometry("805x485")
    IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
    IngridWindow.mainloop()

def set_params():
    """ Interactive mode for Ingrid. Gets the files used, and opens
    a plot of the Efit data. Prompts the user for magx, xpt, and psi levels.
    Saves the data in a namelist file."""
    
    def paws():
        # Helps the code to wait for the user to select points
        raw_input("Press the <ENTER> key to continue...")

    nml = {'files': {}, 'grid_params': {}}

    # get files from the user for data, inner and outer strike plates
    # This is an acceptable form - also takes full path names
    #    Geqdsk_file = "../data/SNL/neqdsk"
    #    Inner_plate_file  = "../data/SNL/itp1"
    #    Outer_plate_file  = "../data/SNL/itp1"
    from pathlib2 import Path
    while True:
        geqsdk_path = Path(raw_input("Enter the geqsdk filename: ").strip())  # remove any whitespace
        if geqsdk_path.is_file():
            break
        print('Provided geqsdk filename could not be found...')
    nml['files']['geqdsk'] = str(geqsdk_path)
    del geqsdk_path

    while True:
        itp_path = Path(raw_input("Enter the inner strike plate filename: ").strip())  # remove any whitespace
        if itp_path.is_file():
            break
        print('Provided inner strike plate filename could not be found...')
    nml['files']['itp'] = str(itp_path)
    del itp_path

    while True:
        otp_path = Path(raw_input("Enter the outer strike plate filename: ").strip())  # remove any whitespace
        if otp_path.is_file():
            break
        print('Provided outer strike plate filename could not be found...')
    nml['files']['otp'] = str(otp_path)
    del otp_path

    # now we can read the data and fine tune our parameters
    grid = Ingrid(nml)
    grid.OMFIT_read_psi()
    grid.read_target_plate()
    grid.calc_efit_derivs()
    grid.plot_efit_data()
    grid.plot_target_plate()
    grid.find_roots()

    # find the null-point (magnetic axis)
    print("Click on the magnetic axis") ##-need blocking here!
    paws()
    magx = grid.save_root()
    nml['grid_params']['rmagx'] = magx[0]
    nml['grid_params']['zmagx'] = magx[1]

    # find the x-point
    print("Click on the x-point")
    paws()
    xpt = grid.save_root()
    nml['grid_params']['rxpt'] = xpt[0]
    nml['grid_params']['zxpt'] = xpt[1]

    # enter or click on the location of the psi values
    psi_magx = grid.efit_psi.get_psi(magx[0], magx[1])
    psi_xpt1 = grid.efit_psi.get_psi(xpt[0], xpt[1])

    # get max psi level
    max_psi = raw_input("Enter max psi level, or press <ENTER> to click on the location: ")
    if max_psi == '':
        grid.toggle_root_finder()
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(max_psi)
    nml['grid_params']['psi_max'] = psi

    # get min psi value near the core
    core_psi = raw_input("Enter min psi level near the core plasma, or press <ENTER> to click on the location: ")
    if core_psi == '':
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(core_psi)
    nml['grid_params']['psi_min_core'] = psi

    # get min psi value near pf region
    pf_psi = raw_input("Enter min psi level near the private flux region, or press <ENTER> to click on the location: ")
    if pf_psi == '':
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(pf_psi)
    nml['grid_params']['psi_min_pf'] = psi

    outFile = raw_input("Enter output filename [default is 'grid_params.nml']: ").strip()
    if outFile == '':
        outFile = 'grid_params.nml'

    f90nml.write(nml, outFile, force=True)  # force tag is to overwrite the previous file
    print("Saved parameters to '{}'.".format(outFile))
    plt.close('Efit Data') # finish with the raw Psi data


def run(param_file = None):
    """ Reads a namelist file containing the parameters for Ingrid, and
    runs the grid generator for a single null configuration.

    Parameters
    ----------
    param_file : str, optional
        String containing path to a user provided parameter *.nml file.
        If not provided, user will be prompted to manually enter a path
        to a parameter file.
    """
    from pathlib2 import Path
    def paws():
        # Helps the code to wait for the user to select points
        raw_input("Press the <ENTER> key to continue...")

    test = IA.IngridApp()
    test.mainloop()
    if param_file:
        path_name = Path(param_file)
        if path_name.is_file() and path_name.suffix == '.nml':
            nml_file = param_file
        del path_name
    else:
        nml_file = raw_input("Enter params filename [params.nml]: ")
        if nml_file == '':
            # default is an example case for snl
            nml_file = 'params.nml'
    
    nml = f90nml.read(nml_file)

    grid = Ingrid(nml)
    grid.setup()
    grid.OMFIT_read_psi()
    grid.read_target_plate()
    grid.calc_psinorm()
    grid.plot_target_plate()
    grid.compute_eq_psi()

    #-construct the patch-map for SNL
    from time import time
    start = time()
    grid.construct_SNL_patches()
    end = time()
    
    # TODO: finish writing the refine patches method.
    # grid.refine_patches()
    
    paws()
    grid.patch_diagram()
    paws()

    print("Time for grid: {} seconds.".format(end-start))


if __name__ == '__main__':
    prep_input()

