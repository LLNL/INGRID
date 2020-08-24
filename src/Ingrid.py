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
from time import time
from GUI import SimpleGUI as sg

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder

from Topologies.SNL import SNL
from Topologies.SF15 import SF15
from Topologies.SF45 import SF45
from Topologies.SF75 import SF75
from Topologies.SF105 import SF105
from Topologies.SF135 import SF135
from Topologies.SF165 import SF165
from Topologies.UDN import UDN

from geometry import Point, Line, Patch, segment_intersect, orientation_between
from scipy.optimize import root, minimize
from collections import OrderedDict

class Ingrid:
    """
    The Ingrid class manages all processes pertaining to patch map and grid generation.
    File I/O, Geometry file prep, parameter instantiation, and configuration classification
    are driven by Ingrid.
    """

    def __init__(self, params = {},**kwargs):
        self.InputFile=None
        self.config=None

        self.settings_lookup = ['grid_params', 'integrator_params', \
                            'target_plates', 'DEBUG']

        self.SetDefaultParams()
        self.process_yaml(params)
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
        self.IngridWindow.title('INGRID')
        self.IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
        self.IngridWindow.mainloop()

    def SaveSettingsFile(self, fname=None, settings={}):
        self.process_yaml(settings)
        if Path(fname).suffix == '.yml':
            with open(fname, 'w') as f:
                yml.dump(self.settings, f)

    def SetDefaultParams(self):

        self.default_grid_params = { \
            'num_xpt' : 1, 'nlevs' : 30, 'full_domain' : True, \
            'psi_max' : 0.0, \
            'psi_core' : 0.0, \
            'psi_pf_1' : 0.0,  \
            'psi_pf_2' : 0.0, \
            'psi_max_west' : 0.0, \
            'psi_max_east' : 0.0, \
            'rmagx' : None, 'zmagx' : None, \
            'rxpt' : 0.0, 'zxpt' : 0.0, 'rxpt2' : 0.0, 'zxpt2' : 0.0, \
            'grid_generation' : { \
                'np_default' : 2, 'nr_default' : 2,
                'poloidal_f_default' : 'x, x', 'radial_f_default' : 'x, x'\
            }, \
            'patch_generation' : { \
                'rmagx_shift' : 0.0, 'zmagx_shift' : 0.0, \
                'west_tilt' : 0.0, 'east_tilt' : 0.0, \
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
            'plate_E1' : {'file' : '', 'name' : '', 'rshift': 0.0, 'zshift' : 0.0}, \
            'plate_E2' : {'file' : '', 'name' : '', 'rshift': 0.0, 'zshift' : 0.0}, \
            'plate_W1' : {'file' : '', 'name' : '', 'rshift': 0.0, 'zshift' : 0.0}, \
            'plate_W2' : {'file' : '', 'name' : '', 'rshift': 0.0, 'zshift' : 0.0}, \
        }

        self.default_limiter_params = {
            'file' : '', 'use_limiter' : False, 'use_efit_bounds' : False, 'rshift' : 0.0, 'zshift' : 0.0
        }

        self.default_patch_data_params = {
            'file' : '', 'use_file' : False,
            'preferences' : {'new_file' : False, 'new_fname' : ''}
        }

        self.default_DEBUG_params = { \
            'visual' : {'find_NSEW' : False, 'patch_map' : False, 'subgrid' : False, 'gridue' : False, 'SF_analysis' : False}, \
            'verbose' : {'target_plates' : False, 'patch_generation' : False, 'grid_generation' : False, 'SF_analysis' : False}
        }

        self.default_values_lookup = {
        'eqdsk' : '', \
        'grid_params' : self.default_grid_params, 'integrator_params' : self.default_integrator_params, \
            'target_plates' : self.default_target_plates_params, 'limiter' : self. default_limiter_params, 
            'patch_data' : self.default_patch_data_params,
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

        if (self.settings['patch_data']['use_file']) \
        and (Path(self.settings['patch_data']['file']).is_file()) \
        and (Path(self.settings['patch_data']['file']).suffix == '.npy'):
            print('Patch data-file:')
            print( ' # Patch data-file:',self.settings['patch_data']['file'])

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
        print( ' # Use patch data-file:',self.settings['patch_data']['use_file'])
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
            if item == 'eqdsk':
                continue

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


    def SetGeometry(self, settings, rshift=None, zshift=None):
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

            if rshift == None:
                rshift = 0.0
            if zshift == None:
                zshift = 0.0

            if type(v) == str:
                if Limiter and v.lower() == 'eq':
                    DefaultEQ = True
                elif Path(v).suffix == '.txt':
                    x, y = self.ParseFileCoordinates(v)
                elif Path(v).suffix == '.npy':
                    x, y = np.load(v)
            elif type(v) in [list, tuple]:
                x, y = v
            elif type(v) == dict:
                if TargetPlate:
                    if v.keys() == self.default_target_plate_params[k].keys():
                        if Path(v['file']).suffix == '.txt':
                            x, y = self.ParseFileCoordinates(v)
                        elif Path(v['file']).suffix == '.npy':
                            x, y = np.load(v)
                        else:
                            raise ValueError(f"# No file provided for target plate '{k}'")
                        rshift = v['rshift']
                        zshift = v['zshift']
                elif Limiter:
                    if v.keys() == self.default_limiter_params.keys():
                        if Path(v['file']).suffix == '.txt':
                            x, y = self.ParseFileCoordinates(v)
                        elif Path(v['file']).suffix == '.npy':
                            x, y = np.load(v)
                        else:
                            DefaultEQ = True
                        rshift = v['rshift']
                        zshift = v['zshift']
                else:
                    for _k, _v in v.items():
                        if _k.lower() in ['x', 'r']:
                            x = _v
                        elif _k.lower() in ['y', 'z']:
                            y = _v
                        elif _k.lower() == 'file':
                            if Path(_v).suffix == '.txt':
                                x, y = self.ParseFileCoordinates(_v)
                            elif Path(_v).suffix == '.npy':
                                x, y = np.load(_v)
                            else:
                                raise ValueError(f"# Invalid file extension '{Path(_v).suffix}' was used within dic in SetGeometry.")
                        elif _k.lower() == 'rshift':
                            rshift = _v
                        elif _k.lower() == 'zshift':
                            zshift = _v
                        else:
                            raise ValueError(f"# Invalid key '{_k}' was used within dic in SetGeometry.")
            else:
                raise ValueError(f"# Invalid data type {type(v)} was used for setting geometry.")

            if TargetPlate: 
                self.SetTargetPlate({k : [x, y]}, rshift=rshift, zshift=zshift)
            elif Limiter and not DefaultEQ:
                self.SetLimiter([x, y], rshift=rshift, zshift=zshift)
            elif Limiter and DefaultEQ:
                self.SetLimiter(rshift=rshift, zshift=zshift)
            else:
                raise ValueError(f"# Invalid key '{k}' was used for setting geometry.")

    def SaveTargetData(self, settings=None, timestamp=False):

        if type(settings) == dict:

            for k, fname in settings.items():
                if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                    k='plate_W1'
                elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                    k='plate_E1'
                elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                    k='plate_W2'
                elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                    k='plate_E2'
                else:
                    raise ValueError(f"# Invalid key '{k}' provided for 'SaveTargetDatafile'")

                v = self.plate_data[k]

                if v:
                    R = np.array(v.xval)
                    Z = np.array(v.yval)
                else:
                    R = np.array([None])
                    Z = np.array([None])

                if timestamp == False:
                    f = fname+'_'+k
                else:
                    # Save to current working directory with timestamp
                    f = fname+'_'+k+'_'+str(int(time()))
                np.save(f, np.array([R,Z]))
                print("# Saved target plate data to file "+f+'.npy')

        elif settings == None:
            for k, v in self.plate_data.items():
                if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                    k='plate_W1'
                elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                    k='plate_E1'
                elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                    k='plate_W2'
                elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                    k='plate_E2'
                else:
                    raise ValueError(f"# Invalid key '{k}' provided for 'SaveTargetDatafile'")

                if v:
                    R = np.array(v.xval)
                    Z = np.array(v.yval)
                else:
                    R = np.array([None])
                    Z = np.array([None])

                if timestamp == False:
                    f = k+'_data'
                else:
                    # Save to current working directory with timestamp
                    from time import time
                    f = k+'_data_'+str(int(time()))
                np.save(f, np.array([R,Z]))
                print("# Saved target plate data to file "+f+'.npy')
        else:
            raise ValueError(f"# Error in SaveTargetDatafile handling parameter of type {type(settings)}. Must be dict or None.")



    def LoadTargetData(self, settings):
        if type(settings) == dict:

            for k, fname in settings.items():
                if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                    k='plate_W1'
                elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                    k='plate_E1'
                elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                    k='plate_W2'
                elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                    k='plate_E2'
                else:
                    raise ValueError(f"# Invalid key '{k}' provided for 'LoadTargetDatafile'")

                try:
                    # Protect against None type
                    Path(fname).is_file()
                except:
                    raise ValueError(f"# Error in 'LoadTargetDatafile()': File '{fname}' does not exist")

                suffix = Path(fname).suffix

                if suffix == '.npy':
                    R, Z = np.load(fname)
                elif suffix == '.txt':
                    R, Z = self.ParseFileCoordinates(fname)
                else:
                    raise ValueError(f"# Error in 'LoadTargetDatafile()': file type '{suffix}' is not supported (requires '.txt', '.npy')")

                if R[0] is None or Z[0] is None:
                    self.plate_data[k] = None
                else:
                    self.plate_data[k] = Line([Point(P) for P in zip(R, Z)])

                print(f"# Loaded target plate '{k}'' data from file "+fname)

        else:
            raise ValueError(f"# Error in LoadTargetDatafile handling parameter of type {type(settings)}. Must be dict with format {{PLATE_KEY : FPATH}}")

    def SaveLimiterData(self, fpath=None):

        try:
            R = np.array(self.limiter_data.xval)
            Z = np.array(self.limiter_data.yval)
        except:
            raise ValueError("# INGRID object contains no limiter data to save!")

        if type(fpath) == str:
            f = fpath
        else:
            # Save to current working directory with timestamp
            from time import time
            f = str(Path('').cwd())+'/'+'LimiterData_'+str(int(time()))
        np.save(f, np.array([R,Z]))
        print("# Saved limiter data to file "+f+'.npy')

    def LoadLimiterData(self, fpath=None):
        
        if fpath:
            f = fpath
        else:
            fpath = self.settings['limiter']['file']

        if f and Path(f).is_file():
            suffix = Path(f).suffix
            if suffix == '.npy':
                R, Z = np.load(f)
            elif suffix == '.txt':
                R, Z = self.ParseFileCoordinates(f)
            else:
                raise ValueError(f"# Error in 'LoadLimiterDatafile()': file type '{suffix}' is not supported (requires '.txt', '.npy')")

            self.limiter_data = Line([Point(P) for P in zip(R, Z)])
            print(f"# Loaded limiter data from file "+f)
        else:
            raise ValueError(f"# Error in 'LoadLimiterDatafile()': File '{f}' does not exist")
                


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

        if len(RLIM) == 0 or len(ZLIM) == 0:
            use_efit_bounds=True

        if use_efit_bounds:
            coordinates = [(self.efit_psi.rmin + 1e-2, self.efit_psi.zmin + 1e-2),
                               (self.efit_psi.rmax - 1e-2, self.efit_psi.zmin + 1e-2),
                               (self.efit_psi.rmax - 1e-2, self.efit_psi.zmax - 1e-2),
                               (self.efit_psi.rmin + 1e-2, self.efit_psi.zmax - 1e-2),
                               (self.efit_psi.rmin + 1e-2, self.efit_psi.zmin + 1e-2)]

            RLIM, ZLIM = [], []
            for c in coordinates:
                r, z = c
                RLIM.append(r)
                ZLIM.append(z)

            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = np.array(RLIM) - rshift, np.array(ZLIM) - zshift

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
        self.plate_data[k]=Line([Point(x-rshift,y-zshift) for x,y in data[index]])
        self.OrderTargetPlate(k)

    def read_target_plates(self):
        for plate in self.settings['target_plates']:
            try:
                self.SetGeometry({plate : self.settings['target_plates'][plate]['file']}, 
                    rshift=self.settings['target_plates'][plate]['zshift'],
                    zshift = self.settings['target_plates'][plate]['zshift'])
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

        self.config = self.eq.config

    def AnalyzeTopology(self,verbose=False):

        try:
            visual = self.settings['DEBUG']['visual']['find_NSEW']
        except:
            visual = False

        self.PrepLineTracing(interactive=False)
        self.ClassifyTopology(visual=visual)
        
        # Updated after call to ClassifyTopology --> self.config

        print('')
        print('# Identified {} configuration.'.format(self.config))
        print('')

        self.SetTopology(self.config)

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

        elif topology == 'SF105':
            ingrid_topology = SF105(self, topology)

        elif topology == 'SF135':
            ingrid_topology = SF135(self, topology)

        elif topology == 'SF165':
            ingrid_topology = SF165(self, topology)

        self.current_topology = ingrid_topology

    def get_config(self):
        return self.eq.config

    def ExportGridue(self, fname = 'gridue'):
        """ Saves the grid as an ascii file """
        self.PrepGridue()
        if type(self.current_topology) in [SNL]:
            if self.WriteGridueSNL(self.current_topology.gridue_params, fname):
                print(f"# Successfully saved gridue file as '{fname}'")
        elif type(self.current_topology) in [SF15, SF45, SF75, SF105, SF135, SF165, UDN]:
            if self.WriteGridueDNL(self.current_topology.gridue_params, fname):
                print(f"# Successfully saved gridue file as '{fname}'")

    def PrepGridue(self):
        self.current_topology.SetupPatchMatrix()
        self.current_topology.concat_grid()
        self.current_topology.set_gridue()

    def WriteGridueSNL(self, gridue_params, fname = 'gridue'):

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
        return True

    def CreateSubgrid(self,ShowVertices=False,color=None,RestartScratch=False,NewFig=True,OptionTrace='theta',ExtraSettings={},ListPatches='all',Enforce=False):

        self.current_topology.RefreshSettings()
        print('Value Check for local plates:')
        _i = self.current_topology.settings['target_plates']
        if NewFig:
            try:
                plt.close(self.FigGrid)
            except:
                pass
            self.FigGrid=plt.figure(f'INGRID: {self.current_topology.config} Grid', figsize=(6, 10))
            ax=self.FigGrid.gca()
            ax.set_xlim([self.efit_psi.rmin, self.efit_psi.rmax])
            ax.set_ylim([self.efit_psi.zmin, self.efit_psi.zmax])
            ax.set_aspect('equal', adjustable='box')

            ax.set_xlabel('R')
            ax.set_ylabel('Z')
            ax.set_title(f'{self.current_topology.config} Grid Diagram')
            ax.set_aspect('equal', adjustable='box')

        self.GetConnexionMap()

        self.current_topology.construct_grid(ShowVertices=ShowVertices,RestartScratch=RestartScratch,OptionTrace=OptionTrace,ExtraSettings=ExtraSettings,ListPatches=ListPatches,Enforce=Enforce)
        return True

    def SetMagReference(self:object,topology:str='SNL')->None:
        self.magx= (self.settings['grid_params']['rmagx'], self.settings['grid_params']['zmagx'])
        self.Xpt = (self.settings['grid_params']['rxpt'],self.settings['grid_params']['zxpt'])
        self.xpt1 = self.Xpt

        if self.settings['grid_params']['num_xpt']==2:
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
                  'psi_core':'cyan',
                  'psi_pf_1':'white'}
        elif nxpt == 2:
            Dic={ 'psi_core':'cyan',
                  'psi_max_west':'magenta',
                  'psi_max_east':'lime',
                  'psi_pf_1':'white',
                  'psi_pf_2':'yellow'}
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

    def PlotMidplane(self):
        try:
            west_tilt = self.settings['grid_params']['patch_generation']['west_tilt']
        except KeyError:
            west_tilt = 0.0
        try:
            east_tilt = self.settings['grid_params']['patch_generation']['east_tilt']
        except KeyError:
            east_tilt = 0.0

        magx = Point(np.array([self.settings['grid_params']['rmagx'] + self.settings['grid_params']['patch_generation']['rmagx_shift'], \
            self.settings['grid_params']['zmagx'] + self.settings['grid_params']['patch_generation']['zmagx_shift']]))

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx.x - 1e6 * np.cos(west_tilt), magx.y - 1e6 * np.sin(west_tilt))
        RHS_Point = Point(magx.x, magx.y)
        west_midLine = Line([LHS_Point, RHS_Point])
        west_midLine.plot('white')

        LHS_Point = Point(magx.x, magx.y)
        RHS_Point = Point(magx.x + 1e6 * np.cos(east_tilt), magx.y + 1e6 * np.sin(east_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])
        east_midLine.plot('white')

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
        self.PrepLineTracing()

    def ShowSetup(self):
        self.plot_psinorm()
        self.PlotPsiNormMagReference()
        self.PlotPsiNormBounds()
        self.plot_strike_geometry()
        self.PlotMidplane()

    def CheckPatches(self,verbose=False):
        self.current_topology.CheckPatches(verbose=verbose)

    def ConstructPatches(self):
        self.current_topology.RefreshSettings()
        self.GetPatchTagMap(self.current_topology.get_config())
        self.PrepLineTracing(interactive=False)
        self.current_topology.construct_patches()
        self.current_topology.GroupPatches()
        self.CheckPatches()

        if self.settings['patch_data']['preferences']['new_file']:
            self.SavePatches(self.settings['patch_data']['preferences']['new_fname'])

    def SavePatches(self, fname=''):
        if fname in ['', None]:
            fname = self.current_topology.config + '_patches_' + str(int(time()))

        patch_data = [patch.as_np() for patch in self.current_topology.patches.values()]
        patch_data.insert(0, self.current_topology.config)
        patch_data.insert(1, self.current_topology.eq.NSEW_lookup)
        np.save(fname, np.array(patch_data))
        print(f"# Saved patch data for file {fname}.npy")

    def LoadPatches(self, fname=None):
        if fname == None:
            fname = self.settings['patch_data']['file']
        config, xpt_data, patches = self.ReconstructPatches(fname)
        self.process_yaml(Ingrid.ReadyamlFile(self.InputFile))
        self.Setup()
        self.eq.NSEW_lookup = xpt_data
        self.SetTopology(config)
        self.current_topology.patches = patches
        self.current_topology.OrderPatches()
        self.current_topology.SetupPatchMatrix()
        self.CheckPatches()

    def ReconstructPatches(self, fname):
        '''
        Reconstruct a patch from Patch.as_np() formatted npy file
        '''
        data = np.load(fname, allow_pickle=True)
        config = data[0]
        xpt_data = data[1]

        patches = {}

        for raw_patch in data[2:]:
            patch_data, cell_data, patch_settings = raw_patch
            NR = patch_data[0]
            NZ = patch_data[1]
            ER = patch_data[2]
            EZ = patch_data[3]
            SR = patch_data[4]
            SZ = patch_data[5]
            WR = patch_data[6]
            WZ = patch_data[7]

            N = Line([Point(p) for p in zip(NR, NZ)])
            E = Line([Point(p) for p in zip(ER, EZ)])
            S = Line([Point(p) for p in zip(SR, SZ)])
            W = Line([Point(p) for p in zip(WR, WZ)])

            patch = Patch([N, E, S, W], patchName=patch_settings['patchName'], platePatch=patch_settings['platePatch'],
                    plateLocation=patch_settings['plateLocation'], PatchTagMap=patch_settings['PatchTagMap'])
            #patch.cell_grid = cell_grid
            patches[patch.patchName] = patch

        patches = OrderedDict([(k, v) for k,v in patches.items()])
        
        return config, xpt_data, patches


    def ShowPatchMap(self):
        try:
            plt.close(self.PatchFig)
        except:
            pass
        self.PatchFig = plt.figure('INGRID: ' + self.current_topology.config + ' Patches', figsize=(6, 10))
        self.PatchAx = self.PatchFig.add_subplot(111)
        self.current_topology.patch_diagram(fig=self.PatchFig, ax=self.PatchAx)
        self.plot_strike_geometry()

    def GetConnexionMap(self):
        alpha = {
            'A' : 'A',
            'B' : 'A', 
            'C' : 'B', 
            'D' : 'C',
            'E' : 'D',
            'F' : 'E', 
            'G' : 'F',
            'H' : 'G',
            'I' : 'H'
            }
        if isinstance(self.current_topology, SNL):
            # SNL connexion map based off patch tag
            ConnexionMap = {}
            for patch in self.current_topology.patches.values():
                tag = patch.get_tag()
                if tag[1] == '1':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '2', 'S'),
                        'W' : (alpha[tag[0]] + '1', 'E')}
            self.current_topology.ConnexionMap = ConnexionMap

        elif type(self.current_topology) in [SF15, SF45, SF75, SF105, SF135, SF165, UDN]:
            # DNL connexion map based off patch tag
            ConnexionMap = {}
            for patch in self.current_topology.patches.values():
                tag = patch.get_tag()
                if tag[1] == '1':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '2', 'S'),
                    'E' : (alpha[tag[0]] + '1', 'W')}
                elif tag[1] == '2':
                    ConnexionMap[patch.get_tag()]={'N' : (tag[0] + '3', 'S'),
                    'E' : (alpha[tag[0]] + '1', 'W')}
            self.current_topology.ConnexionMap = ConnexionMap
    
    def GetPatchTagMap(self, config):
        if config == 'LSN':
            PatchTagMap = {
            'A1' : 'IPF', 'A2' : 'IDL',
            'B1' : 'ICB', 'B2' : 'ISB',
            'C1' : 'ICT', 'C2' : 'IST',
            'D1' : 'OCT', 'D2' : 'OST',
            'E1' : 'OCB', 'E2' : 'OSB',
            'F1' : 'OPF', 'F2' : 'ODL',
            }
        elif config == 'USN':
            PatchTagMap = {
            'A1' : 'OPF', 'A2' : 'ODL',
            'B1' : 'OCB', 'B2' : 'OSB',
            'C1' : 'OCT', 'C2' : 'OST',
            'D1' : 'ICT', 'D2' : 'IST',
            'E1' : 'ICB', 'E2' : 'ISB',
            'F1' : 'IPF', 'F2' : 'IDL',
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




