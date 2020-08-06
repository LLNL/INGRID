#!/usr/bin/env pythonu
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:17:21 2019

@author: watkins35, garcia299
"""
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
import pathlib
import inspect
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt
import yaml as yml
import os
from pathlib import Path
from GUI import IngridApp as IA
from GUI import SimpleGUI as sg

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder

from Topologies.SNL import SNL
from Topologies.SF15 import SF15
from Topologies.SF45 import SF45
from Topologies.SF75 import SF75
from Topologies.UDN import UDN

from geometry import Point, Line, Patch, segment_intersect, orientation_between
from scipy.optimize import root, minimize

class Ingrid:
    """
    The Ingrid class manages all processes pertaining to patch map and grid generation.
    File I/O, Geometry file prep, parameter instantiation, and configuration classification
    are driven by Ingrid.

    Parameters
    ----------
    params : dict (optional)
        Dictionary object containing the following YAML entries:
        'eqdsk'
        'grid_params'
        'integrator_params'
        'target_plates'
        'DEBUG'

        A copy of 'params' is stored within this Ingrid object.
        Should the user not provide 'param', or entries are missing within the YAML file,
        default values will populate the YAML copy that is within the Ingrid object.
    Description of YAML file and 'param' dictionary:
    # =================================================================================================
    # grid_params:
    # -------------------------------------------------------------------------------------------------
    #     Topologically significant parameters for INGRID to
    #     utilize during runtime.
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    #     config - (str) Topology associated with current parameters.
    #               Currently supported: LSN and USN
    # -------------------------------------------------------------------------------------------------
    #     num_xpt - (int) Number of x-points
    # -------------------------------------------------------------------------------------------------
    #     psi_* - (double) Bounds of grid domain determined
    #              by normalized psi efit data.
    # -------------------------------------------------------------------------------------------------
    #     patch_generation - (dict) User settings related to patch map generation.
    #
    #          rmagx_shift - (double) Translate the r-coordinate of the magnetic axis
    #                         used for patch bounds.
    #          zmagx_shift - (double) Translate the z-coordinate of the magnetic axis
    #                         for the patch bounds.
    #          inner_tilt - (double) RADIAN value to tilt the inner horizontal-line patch
    #                        boundary that intersects the magnetic-axis
    #          outer_tilt - (double) RADIAN value to tilt the outer horizontal-line patch
    #                        boundary that intersects the magnetic-axis
    #          use_NW - (bool/logical-int) Trace linearly from primary-xpt in the 'NW' direction
    #                    rather than tracing 'W' in rho.
    #          NW_adjust - (double) RADIAN value to adjust the direction of the 'NW' linear trace.
    #
    #          use_NE - (bool/logical-int) Trace linearly from primary-xpt in the 'NE' direction
    #                    rather than tracing 'E' in rho.
    #          NE_adjust - (double) RADIAN value to adjust the direction of the 'NE' linear trace.
    # -------------------------------------------------------------------------------------------------
    #     grid_generation - (dict) User settings related to grid generation.
    #
    #          np_global - (int) Default number of POLOIDAL cells to generate for all
    #                       patches not adjacent to a target plate.
    #          nr_global - (int) Default number of RADIAL cells to generate for all
    #                       patches not adjacent to a target plate.
    #
    #          np_core - (int) Number of POLOIDAL cells in the core plasma region
    #          nr_core - (int) Number of RADIAL cells in the core plasma region
    #          radial_f_core - (str) User defined cell "distortion" function for
    #                           radial cells in the core plasma region.
    #
    #          np_sol - (int) Number of POLOIDAL cells in the scrape-off layer
    #          nr_sol - (int) Number of RADIAL cells in the scrape-off layer
    #          radial_f_sol - (str) User defined cell "distortion" function for
    #                           radial cells in the scrape off layer.
    #
    #          np_pf - (int) Number of POLOIDAL cells in the private-flux region
    #          nr_pf - (int) Number of RADIAL cells in the private-flux region
    #          radial_f_pf - (str) User defined cell "distortion" function for
    #                           radial cells in the private-flux region.
    # -------------------------------------------------------------------------------------------------
    #     rmagx - (double) r coordinate of magnetic axis
    # -------------------------------------------------------------------------------------------------
    #     rxpt - (double) r coordinate of primary x-point
    # -------------------------------------------------------------------------------------------------
    #     zmagx - (double) z coordinate of magnetic axis
    # -------------------------------------------------------------------------------------------------
    #     zxpt - (double) z coordinate of primary x-point
    # =================================================================================================
    # =================================================================================================
    # integrator_params:
    # -------------------------------------------------------------------------------------------------
    #     Integrator and line_tracing class settings.
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    #     dt - (double) Integration time-step for line_tracing class.
    # -------------------------------------------------------------------------------------------------
    #     eps - (double) Radius of epsilon neighborhood around x-point
    #            containing seed-point directions associated with NSEW.
    # -------------------------------------------------------------------------------------------------
    #     first_step - (double) Initial step size for LSODA integrator.
    # -------------------------------------------------------------------------------------------------
    #     step_ratio - (double) A ratio of RZ-domain dimensions to
    #                   determine max_step parameter for LSODA integrator.
    # -------------------------------------------------------------------------------------------------
    #     tol - (double) Tolerance value defining convergence criterion
    #            for integration with line_tracing class.
    # =================================================================================================
    # =================================================================================================
    # target_plates:
    # -------------------------------------------------------------------------------------------------
    #     Settings for all target plates to be
    #     considered during INGRID runtime.
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    #     plate_E1: (dict) Plate EAST of primary x-point (where NORTH
    #                is defined as the path towards magnetic axis).
    # -------------------------------------------------------------------------------------------------
    #     plate_W1: (dict) Plate WEST of primary x-point (where NORTH
    #                is defined as the path towards magnetic axis).
    # -------------------------------------------------------------------------------------------------
    #     file: (str) Text file containing coordinates defining a target
    #            plate.
    # -------------------------------------------------------------------------------------------------
    #     name: (str) User-defined name associated with plate.
    # -------------------------------------------------------------------------------------------------
    #     np_local: (int) Number of POLOIDAL cells to generate within
    #                patches adjacent to the respective target plate.
    # -------------------------------------------------------------------------------------------------
    #     nr_local: (int) Number of RADIAL cells to generate within
    #                patches adjacent to the respective target plate.
    # -------------------------------------------------------------------------------------------------
    #     poloidal_f: (str) User defined function for 'distortion' of
    #                  cell placement within a target plate.
    # -------------------------------------------------------------------------------------------------
    #     USAGE OF "poloidal_f":
    #                 - The sequence of characters before ',' must
    #                   denote the variable name VAR to operate upon.
    #
    #                 - Place a user-defined mathematical expression
    #                   after the ',' delimiter utilizing the variable
    #                   VAR defined before-hand.
    #
    #                 - The package SymPy is used to convert user-input
    #                   to a mathematical expression. This requires the
    #                   input to be a valid Python expression.
    #
    #                 - For detailed information on available functions,
    #                   see Sympy docs pertaining to "Functions" module.
    # =================================================================================================
    # =================================================================================================
    # DEBUG:
    # -------------------------------------------------------------------------------------------------
    #     Controls for DEBUG mode in INGRID.
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    #     visual - (dict) Dictionary of bool/logical-int values to
    #               activate live-plotting during INGRID operations.
    #
    #         find_NSEW - (bool/logical-int) Show line_tracing during
    #                      search for North/South from primary x-point
    #         patch_map - (bool/logical-int) Show line_tracing during
    #                      generation of patch-map
    #         subgrid - (bool/logical-int) Show line_tracing during
    #                      generation of subgrids within each patch
    #         gridue - (bool/logical-int) Plot gridue data with
    #                      guard cells in order of increasing poloidal
    #                      and radial indices
    # -------------------------------------------------------------------------------------------------
    #     verbose - (dict) Dictionary of bool/logical-int values to
    #                activate verbose output during INGRID operations.
    #
    #         target_plates - (bool/logical-int) Print all target plate
    #                          coordinates to terminal.
    #         patch_generation - (bool/logical-int) Print all patch names
    #                             and intersection/convergence events
    #                             that occur during patch generation
    #         grid_generation - (bool/logical-int) Print cell information
    #                            of a patch, and the populating of arrays
    #                            containing gridue data during meshgrid
    #                            generation.
    # =================================================================================================
    """

    def __init__(self, params = {},**kwargs):
        self.InputFile=None

        self.settings_lookup = ['grid_params', 'integrator_params', \
                            'target_plates', 'DEBUG']

        self.SetDefaultParams()
        self.process_yaml(params)
        self.settings['eqdsk']=None
        #Process input
        self.ProcessKeywords(**kwargs)

        self.PrintSummaryInput()

    def StartGUI(self):
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
                self.IngridWindow.destroy()
        self.IngridWindow = sg.IngridApp(IngridSession=self)
        self.IngridWindow.title('INGRID v0.0')
        self.IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
        self.IngridWindow.mainloop()

    def SimpleGUI(self):
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
                self.IngridWindow.destroy()
        self.IngridWindow = sg.Ingrid_GUI(IngridSession=self)
        self.IngridWindow.title('INGRID v0.1')
        self.IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
        self.IngridWindow.mainloop()

    def SetDefaultParams(self):

        # TODO: Establish parameter conventions consistent for both SNL and DNL topologies

        # TODO: Provide default values for GUI version of SetValues function (see IngridApp.py)
        self.default_grid_params = { \
            'config' : '', 'num_xpt' : 1, 'nlevs' : 30, 'full_domain' : True, \
            'psi_max' : 0.0, 'psi_max_r' : 0.0, 'psi_max_z' : 0.0, \
            'psi_min_core' : 0.0, 'psi_min_core_r' : 0.0, 'psi_min_core_z' : 0.0, \
            'psi_min_pf' : 0.0, 'psi_min_pf_r' : 0.0, 'psi_min_pf_z' : 0.0, \
            'psi_pf2' : 0.0, 'psi_pf2_r' : 0.0, 'psi_pf2_z' : 0.0, \
            'psi_max_inner' : 0.0, 'psi_max_r_inner' : 0.0, 'psi_max_z_inner' : 0.0, \
            'psi_max_outer' : 0.0, 'psi_max_r_outer' : 0.0, 'psi_max_z_outer' : 0.0, \
            'rmagx' : None, 'zmagx' : None, \
            'rxpt' : 0.0, 'zxpt' : 0.0, 'rxpt2' : 0.0, 'zxpt2' : 0.0, \
            'grid_generation' : { \
                'np_global' : 3, 'nr_global' : 2, \
                'np_primary_sol' : 2, 'np_core' : 2, 'np_primary_pf' : 2, \
                'np_secondary_sol' : 2, 'np_secondary_pf' : 2, \
                'nr_primary_sol' : 2, 'nr_core' : 2, 'nr_primary_pf' : 2, \
                'nr_secondary_sol' : 2, 'nr_secondary_pf' : 2, \
                'radial_f_primary_sol' : 'x, x', 'radial_f_secondary_sol' : 'x, x', \
                'radial_f_primary_pf' : 'x, x', 'radial_f_secondary_pf' : 'x, x', \
                'radial_f_core' : 'x, x', 'poloidal_f' : 'x, x',\

                # Temporary attributes
                'np_sol' : 2, 'nr_sol' : 2, 'np_pf' : 2, 'nr_pf' : 2
            }, \
            'patch_generation' : { \
                'rmagx_shift' : 0.0, 'zmagx_shift' : 0.0, \
                'inner_tilt' : 0.0, 'outer_tilt' : 0.0, \
                'rxpt' : 0.0, 'zxpt' : 0.0, \
                'use_NW' : False, 'use_NE' : False, \
                'use_secondary_NW' : False, 'use_secondary_NE' : False, \
                'use_SW' : False, 'use_SE' : False, \
                'NW_adjust' : -0.785398, 'NE_adjust' : 0.785398, \
                'secondary_NW_adjust' : -0.785398, 'secondary_NE_adjust' : 0.785398, \
            } \
        }

        self.default_integrator_params = { \
            'dt' : 0.01, 'eps' : 5e-5, \
            'first_step' : 1e-5, 'step_ratio' : 0.02, \
            'tol' : 5e-3 \
        }

        self.default_target_plates_params = { \
            'plate_E1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0, 'data' : {}}, \
            'plate_E2' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0, 'data' : {}}, \
            'plate_W1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0, 'data' : {}},
            'plate_W2' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0, 'data' : {}}
        }

        self.default_limiter_params = {
            'file' : '', 'use_limiter' : False, 'use_efit_bounds' : False, 'rshift' : 0.0, 'zshift' : 0.0, 'data' : {}
        }

        self.default_DEBUG_params = { \
            'visual' : {'find_NSEW' : False, 'patch_map' : False, 'subgrid' : False, 'gridue' : False, 'SF_analysis' : False}, \
            'verbose' : {'target_plates' : False, 'patch_generation' : False, 'grid_generation' : False, 'SF_analysis' : False}
        }

        self.default_values_lookup = {'grid_params' : self.default_grid_params, 'integrator_params' : self.default_integrator_params, \
            'target_plates' : self.default_target_plates_params, 'limiter' : self. default_limiter_params, 
            'DEBUG' : self.default_DEBUG_params}

        self.plate_data = {'plate_W1' : {}, 'plate_E1' : {}, 'plate_W2' : {}, 'plate_E2' : {}}


    def ProcessKeywords(self,**kwargs):
        for k, v in kwargs.items():
            if k=='InputFile' or k=='yaml':
                print('# Processing Input File:',v)
                self.InputFile=v
                self.process_yaml(self.ReadyamlFile(v))
                continue
            if k=='W1TargetFile' or k=='w1':
                self.settings['target_plates']['plate_W1']['file']=v
                continue
            if k=='E1TargetFile' or k=='e1':
                self.settings['target_plates']['plate_E1']['file']=v
                continue
            if k=='E2TargetFile' or k=='e2':
                self.settings['target_plates']['plate_E2']['file']=v
                continue
            if k=='W2TargetFile' or k=='w2':
                self.settings['target_plates']['plate_W2']['file']=v
                continue
            if k in ['LimiterFile', 'limiter', 'wall']:
                self.settings['limiter']['file']=v
            if k=='EqFile' or k=='eq':
                self.settings['eqdsk']=v
                continue
            print('Keyword "'+k +'" unknown and ignored...')

    def PrintSummaryInput(self):
        print('')
        print('Equilibrium File:',self.settings['eqdsk'])

        if self.settings['limiter']['use_limiter']:
            print('Limiter File:')
            if self.settings['limiter']['file'] == '':
                print( ' # Limiter:','Using default eqdsk (RLIM, ZLIM) coordinates')
            else:
                print( ' # Limiter',self.settings['limiter']['file'])
        else:
            print('Target Files:')
            #try:
            print( ' # W1:',self.settings['target_plates']['plate_W1']['file'])
            print( ' # E1:',self.settings['target_plates']['plate_E1']['file'])
            print( ' # E2:',self.settings['target_plates']['plate_E2']['file'])
            print( ' # W2:',self.settings['target_plates']['plate_W2']['file'])
        print('')

    def PrintSummaryParams(self):
        print('')
        print('Summary:')
        print( ' # Number of x-points:',self.settings['grid_params']['num_xpt'])
        print( ' # Use full domain:',self.settings['grid_params']['full_domain'])
        print( ' # Use limiter:',self.settings['limiter']['use_limiter'])
        print('')
        self.PrintSummaryInput()

    def process_yaml(self, params,verbose=False):
        """
        Parse the contents of a YAML dump (dictionary object).
        Parameter:
            - params : dict
            Dictionary object conforming to INGRID YAML requirements.
        @author: garcia299
        """

        def get_default_values(item, sub_item = None, attribute = None):
            """
            Helper function for processing the YAML file dump.
            Determines if the entry within the YAML file dump is valid and currently
            recognized by INGRID.
            @author: garcia299
            """
            try:
                default_values = self.default_values_lookup[item]
            except KeyError:
                print('Key not recognized... Add default values to source code for support.')
                return None

            if item and sub_item and attribute:
                return self.default_values_lookup[item][sub_item][attribute]
            elif item and sub_item:
                return self.default_values_lookup[item][sub_item]
            elif item:
                return self.default_values_lookup[item]

        # First level entries within YAML file dump.
        for item in self.default_values_lookup.keys():

            try:
                # Check to see if the item is in the YAML dump.
                params[item]
            except KeyError:
                if verbose: print('Could not find "{}" in YAML file.'.format(item))
                params[item] = get_default_values(item)
                if verbose: print('Populated "{}" with default value of "{}".\n'.format(item, params[item]))
                continue

            # Second level entries within YAML file dump.
            for sub_item in self.default_values_lookup[item].keys():
                try:
                    params[item][sub_item]
                except KeyError:
                    if verbose: print('Could not find "{}/{}" in YAML file.'.format(item, sub_item))
                    params[item][sub_item] = get_default_values(item, sub_item)
                    if verbose: print('Populated "{}/{}" with default value of "{}".\n'.format(item, sub_item, params[item][sub_item]))
                    continue
                if item in ['grid_params', 'target_plates'] \
                and sub_item in ['patch_generation', 'grid_generation', 'plate_E1', 'plate_W1', 'plate_E2', 'plate_W2']:
                    for plate_attribute in self.default_values_lookup[item][sub_item].keys():
                        try:
                            params[item][sub_item][plate_attribute]
                        except:
                            params[item][sub_item][plate_attribute] = get_default_values(item, sub_item, plate_attribute)

        # Update references to YAML entries.
        self.settings = params
        self.grid_params = params['grid_params']
        self.integrator_params = params['integrator_params']
        self.target_plates = params['target_plates']
        self.DEBUG = params['DEBUG']


    def OMFIT_read_psi(self):
        """
        Python class to read the psi data in from an ascii file.
        Saves the boundary information and generated efit_data instance
        @author: watkins35
        """

        g = OMFITgeqdsk(self.settings['eqdsk'])

        nxefit = g['NW']
        nyefit = g['NH']
        rdim = g['RDIM']
        zdim = g['ZDIM']
        zmid = g['ZMID']
        rgrid1 = g['RLEFT']

        rcenter = g['RCENTR']
        bcenter = g['BCENTR']

        rmagx = g['RMAXIS']
        zmagx = g['ZMAXIS']

        rlimiter = g['RLIM']
        zlimiter = g['ZLIM']

        psi = g['PSIRZ'].T

        # calc array for r and z
        rmin = rgrid1
        rmax = rmin + rdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        # reproduce efit grid
        self.efit_psi = Efit_Data(rmin, rmax, nxefit,
                                  zmin, zmax, nyefit,
                                  rcenter, bcenter,
                                  rlimiter, zlimiter,
                                  rmagx, zmagx, 
                                  name='Efit Data', parent=self)
        self.efit_psi.set_v(psi)

        if self.settings['grid_params']['rmagx'] == None or self.settings['grid_params']['zmagx'] == None:
            self.settings['grid_params']['rmagx'], self.settings['grid_params']['zmagx'] = (self.efit_psi.rmagx, self.efit_psi.zmagx)

        self.OMFIT_psi = g

    def ParseFileCoordinates(self, fpath, rshift=0, zshift=0):

        R, Z = [], []

        if Path(fpath).is_file() and Path(fpath).suffix in ['.txt']:
            try:
                with open(fpath) as f:

                    for line in f:

                        if line.startswith('#'):
                            # move along to the next iteration
                            # this simulates a comment
                            continue
                        # we deal with separator ',' or ' '
                        point = line.replace('D', 'e').replace('(', '').replace(')', '').strip()
                        if point.count(',')>0:
                            x = float(point.split(',')[0])
                            y = float(point.split(',')[1])
                        else:

                            x = float(point.split(' ')[0])
                            y = float(point.split(' ')[1])

                        R.append(x - rshift)
                        Z.append(y - zshift)
            except:
                raise ValueError(f"# Error occured when reading data from file {fpath}:\t 'open(fpath)' error")

        else:
            raise ValueError(f"# Error occur when reading data from file {fpath}:\t file does not exist or is not of extension '*.txt'")

        return R, Z


    def SetGeometry(self, settings):
        """
        SetGeometry

            Allows the user to pass input data for setting target plates or limiter
            coordinates
        """

        for k, v in settings.items():

            TargetPlate = False
            Limiter = False
            DefaultEQ = False

            if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                TargetPlate = True
            elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                TargetPlate = True
            elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                TargetPlate = True
            elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                TargetPlate = True
            elif k.lower() in ['limiter', 'wall']:
                Limiter = True

            if type(v) == str:
                if Limiter and v.lower() == 'eq':
                    DefaultEQ = True
                else:
                    x, y = self.ParseFileCoordinates(v)
            elif type(v) in [list, tuple]:
                x, y = v
            elif type(v) == dict:
                xk, yk = v
                x, y = v[xk], v[yk]
            else:
                raise ValueError(f"# Invalid data type {type(v)} was used for setting geometry.")

            if settings.get('rshift'):
                rshift = settings['rshift']
            else:
                rshift = 0

            if settings.get('zshift'):
                zshift = settings['zshift']
            else:
                zshift = 0

            if TargetPlate: 
                self.SetTargetPlate({k : [x, y]}, rshift=rshift, zshift=zshift)
            elif Limiter and not DefaultEQ:
                self.SetLimiter([x, y], rshift=rshift, zshift=zshift)
            elif Limiter and DefaultEQ:
                self.SetLimiter()
            else:
                raise ValueError(f"# Invalid key '{k}' was used for setting geometry.")

    def SaveTargetData(self):
        with open(self.InputFile, 'r') as f:
            settings = yml.load(f, Loader=yml.Loader)

        plate_data = {}

        for k, v in self.plate_data.items():
            if v:
                plate_data[k] = {'R' : None, 'Z' : None}
                plate_data[k]['R'] = [R[0] for R in v.points()]
                plate_data[k]['Z'] = [Z[1] for Z in v.points()]
            else:
                plate_data[k] = v
        
            settings['target_plates'][k].update({'data' : plate_data[k]})

        with open(self.InputFile, 'w') as f:
            yml.dump(settings, f)


    def LoadTargetData(self):
        with open(self.InputFile, 'r') as f:
            settings = yml.load(f, Loader=yml.Loader)

        for k in settings['target_plates'].keys():
            if settings['target_plates'][k]['data']:
                self.plate_data[k] = Line([Point(p) 
                    for p in zip(settings['target_plates'][k]['data']['R'], settings['target_plates'][k]['data']['Z'])])
            else:
                self.plate_data[k] = None

    def SaveLimiterData(self):
        with open(self.InputFile, 'r') as f:
            settings = yml.load(f, Loader=yml.Loader)

        limiter_data = {'R' : [R[0] for R in self.limiter_data.points()],
                        'Z' : [Z[1] for Z in self.limiter_data.points()]}

        settings['limiter'].update({'data' : limiter_data})

        with open(self.InputFile, 'w') as f:
            yml.dump(settings, f)

    def LoadLimiterData(self):
        with open(self.InputFile, 'r') as f:
            settings = yml.load(f, Loader=yml.Loader)

        self.limiter_data = Line([Point(P) for P in zip(settings['limiter']['data']['R'], settings['limiter']['data']['Z'])])


    def SetLimiter(self, coordinates=None, fpath=None, rshift=None, zshift=None):
        
        if rshift is None:
            rshift = self.settings['limiter']['rshift']
        if zshift is None:
            zshift = self.settings['limiter']['zshift']

        use_efit_bounds = self.settings['limiter']['use_efit_bounds']

        if fpath:
            try:
                print('# Processing file for limiter data : {}'.format(self.settings['limiter']['file']))
                self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = self.ParseFileCoordinates(self.settings['limiter']['file'])
            except:
                raise ValueError(f"# Error in method 'SetLimiter' with fpath={fpath}")

        elif coordinates is not None:
            RLIM, ZLIM = coordinates
            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = [r - rshift for r in RLIM], [z - zshift for z in ZLIM]

        elif coordinates is None:
            g=OMFITgeqdsk(self.settings['eqdsk'])
            RLIM, ZLIM = g['RLIM'], g['ZLIM']
            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = RLIM - rshift, ZLIM - zshift

        if use_efit_bounds:
            coordinates = [(self.efit_psi.rmin + 1e-2 - rshift, self.efit_psi.zmin + 1e-2 - zshift),
                               (self.efit_psi.rmax - 1e-2 - rshift, self.efit_psi.zmin + 1e-2 - zshift),
                               (self.efit_psi.rmax - 1e-2 - rshift, self.efit_psi.zmax - 1e-2 - zshift),
                               (self.efit_psi.rmin + 1e-2 - rshift, self.efit_psi.zmax - 1e-2 - zshift),
                               (self.efit_psi.rmin + 1e-2 - rshift, self.efit_psi.zmin + 1e-2 - zshift)]

            RLIM, ZLIM = [], []
            for c in coordinates:
                r, z = c
                RLIM.append(r - rshift)
                ZLIM.append(z - zshift)

            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = np.array(RLIM), np.array(ZLIM)

        self.OrderLimiter()

    def SetTargetPlate(self, settings, rshift=0, zshift=0):

        for _k, _v in settings.items():
            k = _k
            v = _v
            break

        if k.lower() in ['w1', 'westplate1', 'plate_w1']:
            k='plate_W1'
        elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
            k='plate_E1'
        elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
            k='plate_W2'
        elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
            k='plate_E2'
        else:
            raise ValueError(f"# Invalid key '{k}' provided for 'SetTargetPlate'")

        print(f"# Setting coordinates for '{k}'")

        R, Z = v

        # Make sure coordinates are unique
        data = np.array([c for c in zip(R, Z)])
        a,index=np.unique(data,return_index=True,axis=0)
        index.sort()
        self.plate_data[k]=Line([Point(x,y) for x,y in data[index]])
        self.OrderTargetPlate(k)

    def read_target_plates(self):
        for plate in self.settings['target_plates']:
            try:
                self.SetGeometry({plate : self.settings['target_plates'][plate]['file']})
            except:
                pass

    def OrderTargetPlate(self, plate_key):
        """
        Ensures the target plate points are oriented clockwise around the
        magnetic axis.

        @author: garcia299
        """

        k = plate_key

        if not isinstance(self.plate_data.get(k), Line):
            raise ValueError(f"# Plate '{k}' is not loaded as a Line object. Check that 'SetGeometry({{'{k}' : ... }})' has been run.")

        plate = self.plate_data[k]

        start = plate.p[0]
        end = plate.p[-1]

        # Endpoints are on same vertical line
        if start.x == end.x:
            # If end point is above start point
            if start.y < end.y: 
                # Return target plate as is
                ordered_plate = plate.copy()
            else:
                # Flip plate orientation
                ordered_plate = plate.reverse_copy()

        # Endpoints are on same horizontal line
        elif start.y == end.y:
            # If start point is right of end point
            if start.x > end.x:
                # Return strike plate as is
                ordered_plate = plate.copy()
            else:
                # Flip plate orientation
                ordered_plate = plate.reverse_copy()

        else:

            # Check the angle between the plate endpoints with
            # tail fixed on the magnetic axis
            v_reference = np.array( [ self.settings['grid_params']['rmagx'], 
                self.settings['grid_params']['zmagx'] ])

            v_start = np.array( [ start.x, start.y ] )
            v_end = np.array( [ end.x, end.y ] )

            if orientation_between(v_start, v_end, v_reference) <= 0:
                ordered_plate = plate.copy()
            else:
                ordered_plate = plate.reverse_copy()

            # Gather raw data
            self.plate_data[k] = ordered_plate

    def OrderLimiter(self):
        """
        Ensures limiter geometry is oriented clockwise around magnetic axis
        """

        limiter = Line([Point(p) for p in zip(self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'])]).fluff_copy(100)

        rmid = (self.efit_psi.rmin + self.efit_psi.rmax) / 2
        zmid = (self.efit_psi.zmin + self.efit_psi.zmax) / 2

        lookup = {}
        point_list = []
        for p in limiter.p:
            try:
                lookup[(p.x, p.y)]
            except:
                if (p.x >= rmid) and (p.y >= zmid):
                    lookup[(p.x, p.y)] = p
                    point_list.append(p)

        quadrant_boundary = Line(point_list)

        ordered = False
        start = quadrant_boundary.p[0]
        for p in reversed(quadrant_boundary.p):
            end = p
            if end.coor != start.coor:
                break

        # Endpoints are on same vertical line
        if start.x == end.x:
            # If end point is above start point
            if start.y < end.y: 
                # Return target plate as is
                ordered = True

        # Endpoints are on same horizontal line
        elif start.y == end.y:
            # If start point is left of end point
            if start.x < end.x:
                # Return strike plate as is
                ordered = True

        else:

            # Check the angle between the plate endpoints with
            # tail fixed on the magnetic axis
            v_reference = np.array([self.settings['grid_params']['rmagx'], 
                self.settings['grid_params']['zmagx']])

            v_start = np.array( [ start.x, start.y ] )
            v_end = np.array( [ end.x, end.y ] )

            if orientation_between(v_start, v_end, v_reference) <= 0:
                ordered = True

        self.limiter_data = limiter.copy() if ordered else limiter.reverse_copy()

    def plot_target_plates(self):
        """
        Plot the plate_data stored within the Ingrid object to the current figure.

        @author: watkins35, garcia299
        """
        ind = -1
        colorbank = {0 : 'blue', 1 : 'orange', 2 : 'firebrick', 3 : 'green'}
        for k, v in self.plate_data.items():
            ind += 1
            if type(v) == Line:
                v.plot(label=k, color=colorbank[ind])

    def plot_limiter_data(self):
        self.limiter_data.plot(color='dodgerblue', label='Limiter')

    def plot_strike_geometry(self):
        try:
            if self.settings['grid_params']['num_xpt'] == 2 or self.settings['limiter']['use_limiter']:
                self.plot_limiter_data()
                if self.settings['limiter']['use_limiter'] == False: 
                    self.plot_target_plates()
            else:
                self.plot_target_plates()
        except:
            self.plot_target_plates()

    def calc_efit_derivs(self):
        """
        Calculates the partial derivatives using finite differences.
        Wrapper for the member function in the efit class.
        @author: watkins35
        """
        self.efit_psi.Calculate_PDeriv()

    def plot_efit_data(self):
        """
        Generates the plot that we will be able to manipulate
        using the root finder
        @author: watkins35
        """
        self.efit_psi.plot_data(self.settings['grid_params']['nlevs'])

    def FindMagAxis(self,x:float,y:float)->None:
        # r_bounds = (self.efit_psi.rmin, self.efit_psi.rmax)
        # z_bounds = (self.efit_psi.zmin, self.efit_psi.zmax)
        # sol = minimize(fun=self.efit_psi.PsiFunction, x0=np.array([x,y]),
        #     method='L-BFGS-B', jac=self.efit_psi.Gradient, bounds=[r_bounds, z_bounds])
        sol = root(self.efit_psi.Gradient,[x,y])
        self.settings['grid_params']['rmagx']=sol.x[0]
        self.settings['grid_params']['zmagx']=sol.x[1]
        self.magx = (sol.x[0], sol.x[1])

    def FindXPoint(self,x:float,y:float)->None:
        sol = root(self.efit_psi.Gradient,[x,y])
        self.settings['grid_params']['rxpt']=sol.x[0]
        self.settings['grid_params']['zxpt']=sol.x[1]
        self.xpt1 = (sol.x[0], sol.x[1])

    def FindXPoint2(self,x:float,y:float)->None:
        sol = root(self.efit_psi.Gradient,[x,y])
        self.settings['grid_params']['rxpt2']=sol.x[0]
        self.settings['grid_params']['zxpt2']=sol.x[1]
        self.xpt2 = (sol.x[0], sol.x[1])

    def AutoRefineMagAxis(self):
        self.FindMagAxis(self.settings['grid_params']['rmagx'],self.settings['grid_params']['zmagx'])

    def AutoRefineXPoint(self):
        self.FindXPoint(self.settings['grid_params']['rxpt'],self.settings['grid_params']['zxpt'])

    def AutoRefineXPoint2(self):
        self.FindXPoint2(self.settings['grid_params']['rxpt2'],self.settings['grid_params']['zxpt2'])

    def find_roots(self, tk_controller = None):
        """ Displays a plot, and has the user click on an approximate
        zero point. Uses a root finder to adjust to the more exact point.
        Right click to disable.
        Parameter:
            - tk_controller : Tk object
            A reference to the root TK object (for usage in GUI mode).
        @author: watkins35, garcia299
        """

        self.root_finder = RootFinder(self.efit_psi, controller = tk_controller)

    def toggle_root_finder(self):
        """
        Activates or deactivates the root finder ability. Enables
        the user to save the location where they last clicked.
        @author: watkins35
        """
        self.root_finder.toggle_root_finding()

    def calc_psinorm(self):
        """ Uses magx and xpt1 to normalize the psi data. Furthur calculations
        will use this information
        @author: watkins35, garcia299
        """
        from Interpol.Setup_Grid_Data import Efit_Data

        # use the same bounds as the efit data
        self.psi_norm = Efit_Data(self.efit_psi.rmin, self.efit_psi.rmax,
                                  self.efit_psi.nr, self.efit_psi.zmin,
                                  self.efit_psi.zmax, self.efit_psi.nz,
                                  self.efit_psi.rcenter, self.efit_psi.bcenter,
                                  self.efit_psi.rlimiter, self.efit_psi.zlimiter,
                                  self.efit_psi.rmagx, self.efit_psi.zmagx, 
                                  name='Normalized Efit Data', parent = self)
        psi = self.efit_psi.v
        psi_magx = self.efit_psi.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.efit_psi.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx))/(psi_xpt1 - psi_magx)
        self.psi_norm.set_v(psinorm)
        self.psi_norm.Calculate_PDeriv()

    def plot_psinorm(self,interactive=True, target_plates=False):
        """
        Plot the psi_norm data stored in the Ingrid object.
        """
        self.psi_norm.plot_data(self.settings['grid_params']['nlevs'],interactive)

    def find_psi_lines(self, tk_controller = None):
        self.psi_finder = RootFinder(self.psi_norm, mode = 'psi_finder', controller = tk_controller)

    def PrepLineTracing(self, interactive=True):
        """
        Initializes the line tracing class for the construction
        of the grid.
        Parameter:
            - refresh : boolean
            Re-initialize the LineTracing class.
        @author: watkins35, garcia299
        """
        
        self.eq = LineTracing(self.psi_norm, self.settings, interactive=interactive)
        if interactive:
            self.eq.disconnect()

    def ClassifyTopology(self, visual=False):

        print('')
        print("# Begin classification....")
        print('')

        if self.settings['grid_params']['num_xpt'] == 1:
            self.eq.SNL_find_NSEW(self.xpt1, self.magx, visual)

        elif self.settings['grid_params']['num_xpt'] == 2:
            self.eq.DNL_find_NSEW(self.xpt1, self.xpt2, self.magx, visual)

        self.settings['grid_params']['config'] = self.eq.config

    def AnalyzeTopology(self,verbose=False):

        try:
            visual = self.settings['DEBUG']['visual']['find_NSEW']
        except:
            visual = False

        self.PrepLineTracing(interactive=False)
        self.ClassifyTopology(visual=visual)
        
        # Updated after call to ClassifyTopology --> self.settings['grid_params']['config']

        print('')
        print('# Identified {} configuration.'.format(self.settings['grid_params']['config']))
        print('')

        self.SetTopology(self.settings['grid_params']['config'])

    def SetTopology(self, topology:str):
        if topology in ['LSN', 'USN']:
            ingrid_topology = SNL(self, topology)

        elif topology == 'UDN':
            ingrid_topology = UDN(self, topology)

        elif topology == 'SF15':
            ingrid_topology = SF15(self, topology)

        elif topology == 'SF45':
            ingrid_topology = SF45(self, topology)

        elif topology == 'SF75':
            ingrid_topology = SF75(self, topology)

        self.current_topology = ingrid_topology

    def get_config(self):
        return self.eq.config

    def export(self, fname = 'gridue'):
        """ Saves the grid as an ascii file """
        if isinstance(self.current_topology, SNL):
            if self.WritegridueSNL(self.current_topology.gridue_params, fname):
                print(f"# Saved gridue file as '{fname}'.")
        elif isinstance(self.current_topology, DNL):
            self.WritegridueDNL(self.current_topology.gridue_params, fname)

    def WritegridueSNL(self, gridue_params, fname = 'gridue'):

        def format_header(gridue):
            header_items = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iyseptrx1']
            header = ''
            for item in header_items:
                header += '{}'.format(gridue[item]).rjust(4)

            header += '\n'
            return header

        def format_body(data):

            delim_val = 0
            delim_char = ''
            body = ''

            for n in range(5):
                for j in range(len(data[0])):
                    for i in range(len(data)):
                        delim_val += 1
                        val = np.format_float_scientific(data[i][j][n], precision = 15, unique = False).rjust(23).replace('e', 'D')
                        if delim_val == 3:
                            delim_val = 0
                            delim_char = '\n'
                        body += val + delim_char
                        delim_char = ''

            if delim_val % 3 != 0:
                body += '\n'

            return body


        f = open(fname, mode = 'w')
        f.write(format_header(gridue_params) + '\n')

        body_items = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
        for item in body_items:
            f.write(format_body(gridue_params[item]) + '\n')

        runidg = 'iogridue'
        f.write(runidg + '\n')

        f.close()

        return True

    def WriteGridueDNL(self, gridue_params, fname = 'gridue'):

        def format_header(gridue):
            header_rows = [['nxm', 'nym'], \
                    ['iyseparatrix1', 'iyseparatrix2'], \
                    ['ix_plate1', 'ix_cut1', '_FILLER_', 'ix_cut2', 'ix_plate2'],\
                    ['iyseparatrix3', 'iyseparatrix4'], \
                    ['ix_plate3', 'ix_cut3', '_FILLER_', 'ix_cut4', 'ix_plate4']]

            header = ''
            for header_items in header_rows:
                for item in header_items:
                    header += '{}'.format(gridue[item]).rjust(4)
                header += '\n'
            return header

        def format_body(data):

            delim_val = 0
            delim_char = ''
            body = ''

            for n in range(5):
                for j in range(len(data[0])):
                    for i in range(len(data)):
                        delim_val += 1
                        val = np.format_float_scientific(data[i][j][n], precision = 15, unique = False).rjust(23).replace('e', 'D')
                        if delim_val == 3:
                            delim_val = 0
                            delim_char = '\n'
                        body += val + delim_char
                        delim_char = ''

            if delim_val % 3 != 0:
                body += '\n'

            return body


        f = open(fname, mode = 'w')
        f.write(format_header(gridue_params) + '\n')

        body_items = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
        for item in body_items:
            f.write(format_body(gridue_params[item]) + '\n')

        runidg = 'iogridue'
        f.write(runidg + '\n')

        f.close()

    def CreateSubgrid(self,ShowVertices=False,color=None,RestartScratch=False,NewFig=True,OptionTrace='theta',ExtraSettings={},ListPatches='all',Enforce=False):

        try:
            np_cells = self.settings['grid_params']['grid_generation']['np_global']
        except KeyError:
            np_cells = 2
            print('yaml file did not contain parameter np_global. Set to default value of 2...')

        try:
            nr_cells = self.settings['grid_params']['grid_generation']['nr_global']
        except KeyError:
            nr_cells = 2
            print('yaml file did not contain parameter nr_global. Set to default value of 2...')

        print('Value Check for local plates:')
        _i = self.current_topology.settings['target_plates']
        for plate in _i:
            print('Name: {}\n np_local: {}\n'.format(_i[plate]['name'], _i[plate]['np_local']))
        if NewFig:
            self.FigGrid=plt.figure('INGRID: Grid', figsize=(6, 10))
            ax=self.FigGrid.gca()
            ax.set_title(label=f'{self.current_topology.config} Grid Diagram')
            ax.set_aspect('equal', adjustable='box')

        self.GetConnexionMap()

        self.current_topology.construct_grid(np_cells, nr_cells,ShowVertices=ShowVertices,RestartScratch=RestartScratch,OptionTrace=OptionTrace,ExtraSettings=ExtraSettings,ListPatches=ListPatches,Enforce=Enforce)
        return True

    def SetMagReference(self:object,topology:str='SNL')->None:
        self.magx= (self.settings['grid_params']['rmagx'], self.settings['grid_params']['zmagx'])
        self.Xpt = (self.settings['grid_params']['rxpt'],self.settings['grid_params']['zxpt'])
        self.xpt1 = self.Xpt

        if topology == 'DNL':
            self.xpt2 = (self.settings['grid_params']['rxpt2'], self.settings['grid_params']['zxpt2'])

    def PlotPsiNormMagReference(self):
        (x,y)=self.magx
        self.psi_norm.ax.plot(x,y,'+',color='yellow',ms=15,linewidth=5)
        (x,y)=self.xpt1
        self.psi_norm.ax.plot(x,y,'+',color='white',ms=15,linewidth=5)
        if self.settings['grid_params']['num_xpt'] == 2:
            try:
                (x,y)=self.xpt2
                self.psi_norm.ax.plot(x,y,'+',color='white',ms=15,linewidth=5)
            except:
                pass

    def PlotPsiNormBounds(self:object)->None:
        nxpt = self.settings['grid_params']['num_xpt']
        if nxpt == 1:
            Dic={'psi_max':'lime',
                  'psi_min_core':'cyan',
                  'psi_min_pf':'white'}
        elif nxpt == 2:
            Dic={ 'psi_min_core':'cyan',
                  'psi_max_west':'magenta',
                  'psi_max_east':'lime',
                  'psi_min_pf':'white',
                  'psi_pf2':'yellow'}
        for k,c in Dic.items():
                self.psi_norm.PlotLevel(self.settings['grid_params'][k],color=Dic[k],label=k)

        self.psi_norm.PlotLevel(1.0,color='red',label='Primary Separatrix')
        if nxpt == 2: 
            self.psi_norm.PlotLevel(
            self.psi_norm.get_psi(self.xpt2[0], self.xpt2[1]), color='blue',label = 'Secondary Separatrix'
            )

        self.psi_norm.ax.legend(bbox_to_anchor=(0.5, -0.15), loc='lower center', ncol=len(Dic) // 2)

    def PlotTopologyAnalysis(self):
        self.psi_norm.ax.add_patch(self.eq.RegionPolygon['core'])
        self.psi_norm.ax.add_patch(self.eq.RegionPolygon['pf'])
        self.eq.RegionLineCut.plot(color='white')

    def PlotPsiLevel(efit_psi:object,level:float,Label:str='')->None:
        plt.contour(efit_psi.r, efit_psi.z, efit_psi.v, [level], colors = 'red', label = Label)

    def Setup(self,interactive=True, **kwargs):
        # TODO have an option to not connect event on figure in command-line mode (so no tk needed)

        topology = 'SNL' if self.settings['grid_params']['num_xpt'] == 1 else 'DNL'

        self.OMFIT_read_psi()
        self.calc_efit_derivs()
        self.AutoRefineMagAxis()
        self.AutoRefineXPoint()
        if topology == 'DNL':
            self.AutoRefineXPoint2()
        self.read_target_plates()
        self.SetLimiter()
        self.SetMagReference(topology)
        self.calc_psinorm()
        self.plot_psinorm(interactive=interactive)
        self.plot_strike_geometry()
        self.PrepLineTracing(interactive=interactive)
        self.AnalyzeTopology()

    def ShowSetup(self):
        self.plot_psinorm()
        self.PlotPsiNormMagReference()
        self.PlotPsiNormBounds()
        self.plot_strike_geometry()

    def ConstructPatches(self):
        self.GetPatchTagMap(self.current_topology.get_config())
        self.PrepLineTracing(interactive=False)
        self.current_topology.construct_patches()
        self.current_topology.CheckPatches()

    def ShowPatchMap(self):
        self.current_topology.patch_diagram()

    def GetConnexionMap(self):
        if isinstance(self.current_topology, SNL):
            # SNL connexion map based off patch tag
            ConnexionMap = {}
            for patch in self.current_topology.patches.values():
                tag = patch.get_tag()
                if tag[1] == '1':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '2', 'S')}
            self.current_topology.ConnexionMap = ConnexionMap

        elif isinstance(self.current_topology, SF45):
            # DNL connexion map based off patch tag
            ConnexionMap = {}
            for patch in self.current_topology.patches.values():
                tag = patch.get_tag()
                if tag[1] == '1':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '2', 'S')}
                elif tag[1] == '2':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '3', 'S')}
            self.current_topology.ConnexionMap = ConnexionMap
    
    def GetPatchTagMap(self, config):
        if config in ['SNL', 'LSN', 'USN']:
            PatchTagMap = {
            'A1' : 'IPF', 'A2' : 'IDL',
            'B1' : 'ICB', 'B2' : 'ISB',
            'C1' : 'ICT', 'C2' : 'IST',
            'D1' : 'OCT', 'D2' : 'OST',
            'E1' : 'OCB', 'E2' : 'OSB',
            'F1' : 'OPF', 'F2' : 'ODL',
            }
        else:
            PatchTagMap = {}
            TempLabels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
            for l in TempLabels:
                for i in range(1,4):
                    PatchTagMap[l+str(i)] = l+str(i)

        # Make it bijective.
        PatchNameMap = {}
        for tag, name in PatchTagMap.items():
            PatchNameMap[name] = tag
        return {**PatchTagMap, **PatchNameMap}

    def Importgridue(self, fname:str = 'gridue')->dict:
        """Import UEDGE grid file as dictionary."""
        try:
            f= open(fname, mode = 'r')
            Values = [int(x) for x in next(f).split()]
            HeaderItems = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iyseptrx1']
            gridue_params=dict(zip(HeaderItems, Values))
            next(f)
            BodyItems = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
            Str={ i : [] for i in BodyItems }
            k=iter(Str.keys())
            Key=next(k)
            for line in f:
                if line=='iogridue\n':
                    continue
                if line=='\n':
                    try:
                        Key=next(k)
                    except:
                        continue
                    print(Key)
                else:
                    Str[Key].append(line)
            f.close()
            nx=gridue_params['nxm']+2
            ny=gridue_params['nym']+2
            for k,v in Str.items():
                L=(''.join(v).replace('\n','').replace('D','e')).split()
                l=iter(L)
                vv= next(l)

                data_=np.zeros((nx,ny,5))
                for n in range(5):
                        for j in range(ny):
                            for i in range(nx):

                                data_[i][j][n]=float(vv)

                                try:
                                    vv=next(l)
                                except:
                                    continue
                gridue_params[k]=data_
            return gridue_params
        except Exception as e:
            print(repr(e))

    @classmethod
    def CheckOverlapCells(Grid,Verbose=False):
        from shapely.geometry import Polygon
        r=Grid['rm']
        z=Grid['zm']
        idx=[1,2,4,3]
        p=[]                  
        pinfo=[]    
        Nx=len(r)
        Ny=len(r[0])
        # polygon = [[Polygon([(x,y) for (x,y) in zip(r[i,j,idx],z[i,j,idx]]) for y in range(0,Ny)] for i range(0,Nx)]
                          
        for i in range(Nx):
          for j in range(Ny):
              c=[(r[i,j,idxx],z[i,j,idxx]) for idxx in idx]
              if Verbose: print(c)
              p.append(Polygon(c))
              pinfo.append((i,j))
        ListIntersect=[]  
        for p1,pinfo1 in zip(p,pinfo):
            for p2,pinfo2 in zip(p,pinfo):
              if p1.intersects(p2) and np.sum(abs(np.array(pinfo1)-np.array(pinfo2)))>2:
                  ListIntersect.append((pinfo1,pinfo2))
                  print('p1:{} and p2:{} intersect!'.format(pinfo1,pinfo2))
        return ListIntersect
  
                
    @staticmethod    
    def Plotgridue(GridueParams,Fill=False,ax=None,Verbose=False,facecolor=None,edgecolor='black'):
        """Plot UEDGE grid."""
        if Verbose:
            print('PlotGridue')
        r=GridueParams['rm']
        z=GridueParams['zm']
        Nx=len(r)
        Ny=len(r[0])
        patches=[]
        if ax is None:
            ax=plt.gca()
        idx=[np.array([1,2,4,3,1])]
        for i in range(Nx):
            for j in range(Ny):
                    p=matplotlib.patches.Polygon(np.concatenate((r[i][j][idx],z[i][j][idx])).reshape(2,5).T,closed=True,fill=False,edgecolor=edgecolor)
                    if not Fill:
                        ax.add_patch(p) #draw the contours
                    else:
                            patches.append(p)

        if Fill:
            Collec=matplotlib.collections.PatchCollection(patches,match_original=True,facecolors=None,edgecolors=edgecolor)
            #Collec.set_facecolor(None)
            ax.add_collection(Collec)
        plt.ylim(z.min(),z.max())
        plt.xlim(r.min(),r.max())
        plt.show()

    @staticmethod
    def ReadyamlFile(FileName:str)->dict:
        '''
        Read a yaml file and return a dictionnary

        Parameters
        ----------
        FileName : str


        Raises
        ------
        IOError


        Returns
        -------
        dict
            Params dictionnary

        '''
        File=pathlib.Path(FileName).expanduser()
        if File.exists() and File.is_file():
            try:
                with open(FileName, 'r') as f:
                    ymlDict=yml.load(f, Loader=yml.Loader)
                return ymlDict
            except:
                raise IOError('Cannot load the yml file: '+FileName)

        else:
            print('Current Working Directory:'+os.getcwd())
            raise IOError('Cannot find the file: ' +FileName)


def QuickStart():
    QuickGrid = Ingrid()
    QuickGrid.SimpleGUI()




