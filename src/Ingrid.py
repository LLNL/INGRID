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
import yaml as _yaml_
import os
from GUI import IngridApp as IA

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder
from Topologies import SNL
from geometry import Point, Line, SNL_Patch, DNL_Patch, segment_intersect
from scipy.optimize import root

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

        self.yaml_lookup = ['grid_params', 'integrator_params', \
                            'target_plates', 'DEBUG']

        self.SetDefaultParams()
        self.process_yaml(params)
        self.yaml['eqdsk']=None
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
        self.IngridWindow = IA.IngridApp(IngridSession=self)
        self.IngridWindow.title('INGRID v0.0')
        self.IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
        self.IngridWindow.mainloop()

    def SetDefaultParams(self):

        # TODO: Establish parameter conventions consistent for both SNL and DNL topologies

        # TODO: Provide default values for GUI version of SetValues function (see IngridApp.py)
        self.default_grid_params = { \
            'config' : '', 'num_xpt' : 1, 'nlevs' : 30, \
            'psi_max' : 0.0, 'psi_max_r' : 0.0, 'psi_max_z' : 0.0, \
            'psi_min_core' : 0.0, 'psi_min_core_r' : 0.0, 'psi_min_core_z' : 0.0, \
            'psi_min_pf' : 0.0, 'psi_min_pf_r' : 0.0, 'psi_min_pf_z' : 0.0, \
            'psi_pf2' : 0.0, 'psi_pf2_r' : 0.0, 'psi_pf2_z' : 0.0, \
            'psi_max_inner' : 0.0, 'psi_max_r_inner' : 0.0, 'psi_max_z_inner' : 0.0, \
            'psi_max_outer' : 0.0, 'psi_max_r_outer' : 0.0, 'psi_max_z_outer' : 0.0, \
            'rmagx' : 0.0, 'zmagx' : 0.0, \
            'rxpt' : 0.0, 'zxpt' : 0.0, 'rxpt2' : 0.0, 'zxpt2' : 0.0, \
            'grid_generation' : { \
                'np_global' : 3, 'nr_global' : 2, \
                'np_primary_sol' : 2, 'np_core' : 2, 'np_primary_pf' : 2, \
                'np_secondary_sol' : 2, 'np_secondary_pf' : 2, \
                'nr_primary_sol' : 2, 'nr_core' : 2, 'nr_primary_pf' : 2, \
                'nr_secondary_sol' : 2, 'nr_secondary_pf' : 2, \
                'radial_f_primary_sol' : 'x, x', 'radial_f_secondary_sol' : 'x, x', \
                'radial_f_primary_pf' : 'x, x', 'radial_f_secondary_pf' : 'x, x', \
                'radial_f_core' : 'x, x', \

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
                'NW_adjust' : 0.0, 'NE_adjust' : 0.0, \
                'secondary_NW_adjust' : 0.0, 'secondary_NE_adjust' : 0.0, \
            } \
        }

        self.default_integrator_params = { \
            'dt' : 0.01, 'eps' : 5e-5, \
            'first_step' : 1e-5, 'step_ratio' : 0.02, \
            'tol' : 5e-3 \
        }

        self.default_target_plates_params = { \
            'plate_E1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0}, \
            'plate_E2' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0}, \
            'plate_W1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0}, \
            'plate_W2' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0} \
        }

        self.default_DEBUG_params = { \
            'visual' : {'find_NSEW' : False, 'patch_map' : False, 'subgrid' : False, 'gridue' : False}, \
            'verbose' : {'target_plates' : False, 'patch_generation' : False, 'grid_generation' : False}
        }

        self.default_values_lookup = {'grid_params' : self.default_grid_params, 'integrator_params' : self.default_integrator_params, \
            'target_plates' : self.default_target_plates_params, 'DEBUG' : self.default_DEBUG_params}

    def ProcessKeywords(self,**kwargs):
        for k, v in kwargs.items():
            if k=='InputFile' or k=='yaml':
                print('# Processing Input File:',v)
                self.InputFile=v
                self.process_yaml(Ingrid.ReadyamlFile(v))
                continue
            if k=='W1TargetFile' or k=='w1':
                self.yaml['target_plates']['plate_W1']['file']=v
                continue
            if k=='E1TargetFile' or k=='e1':
                self.yaml['target_plates']['plate_E1']['file']=v
                continue
            if k=='E2TargetFile' or k=='e2':
                self.yaml['target_plates']['plate_E2']['file']=v
                continue
            if k=='W2TargetFile' or k=='w2':
                self.yaml['target_plates']['plate_W2']['file']=v
                continue
            if k=='EqFile' or k=='eq':
                self.yaml['eqdsk']=v
                continue
            print('Keyword "'+k +'" unknown and ignored...')
    def PrintSummaryInput(self):
        print('Equilibrium File:',self.yaml['eqdsk'])
        print('Target Files:')
        #try:
        print( ' # W1:',self.yaml['target_plates']['plate_W1']['file'])
        print( ' # E1:',self.yaml['target_plates']['plate_E1']['file'])
        print( ' # E2:',self.yaml['target_plates']['plate_E2']['file'])
        print( ' # W2:',self.yaml['target_plates']['plate_W2']['file'])
       # except:
          #  print('error')
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
        self.yaml = params
        self.grid_params = params['grid_params']
        self.integrator_params = params['integrator_params']
        self.target_plates = params['target_plates']
        self.DEBUG = params['DEBUG']

    def print_yaml(self):
        """
        Print the YAML file copy stored within the Ingrid object.
        @author: garcia299
        """
        for item in self.yaml.keys():
            print('=' * 80)
            print('INGRID {}: {}'.format(item, self.yaml[item]))
            print('=' * 80 + '\n')

    def OMFIT_read_psi(self):
        """
        Python class to read the psi data in from an ascii file.
        Saves the boundary information and generated efit_data instance
        @author: watkins35
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

        rmagx = g['RMAXIS']
        zmagx = g['ZMAXIS']

        rlimiter = g['RLIM']
        zlimiter = g['ZLIM']

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
                                  rcenter, bcenter,
                                  rlimiter, zlimiter,
                                  rmagx, zmagx, name='Efit Data')
        self.efit_psi.set_v(psi)

        self.yaml['grid_params']['rmagx'], self.yaml['grid_params']['zmagx'] = (self.efit_psi.rmagx, self.efit_psi.zmagx)

        self.OMFIT_psi = g

    def read_target_plates(self):
        """
        Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        @author: watkins35, garcia299
        """

        # Determine if in DEBUG verbose mode.
        try:
            debug = self.yaml['DEBUG']['verbose']['target_plates']
        except:
            debug = False

        # Create references and plate_data container.
        target_plates = self.yaml['target_plates']
        plate_data = {}

        for plate in target_plates:

            # Check for valid file path.
            #JG: not good because it does not check if file exists
            # That only makes sense when used with the GUI where file existence is automatic with the FilePicker method get_file()
            # We assume that not relevant plate files are either set to None
            # Could be managed through knowing the topology since we know which plate are needed
            print('# Processing file for plate {} : {}'.format(plate,target_plates[plate]['file']))




            # Some target plate geometries have an associated 'z-shift' that
            # translates the geometry. Application of 'z-shift' occurs during reading.

            zshift = target_plates[plate]['zshift']

            # plate_data contains a list of geometric coordinates, and the min/max psi values
            # that lie on the plate geometry. This can be used (and is) for checking for valid
            # psi entries during patch map generation.
            plate_data[plate] = {'coordinates' : [], 'psi_min_rz' : (), 'psi_max_rz' : ()}

            # Associated label with a target_plate (remove this feature?)
            if not 'name' in target_plates[plate].keys():
                target_plates[plate].update({'name' : plate})
            elif target_plates[plate]['name'] == '':
                target_plates[plate]['name'] = plate

            if target_plates[plate]['file'] is not None:
                try:
                    with open(target_plates[plate]['file']) as f:
                        #First we check if zshift is given
                        for line in f:
                            if line.count('zshift')>0:
                                zshift=float(line.split('=')[-1])
                                if debug: print('zshift=',zshift)
                        f.seek(0)
                        for line in f:

                            point = line.strip()
                            if point.startswith('#'):
                                # move along to the next iteration
                                # this simulates a comment
                                continue

                            if line.count('zshift')>0:
                                continue
                            # we deal with separator ',' or ' '

                            if point.count(',')>0:
                                x = float(point.split(',')[0])
                                y = float(point.split(',')[1])
                            else:

                                x = float(point.split(' ')[0])
                                y = float(point.split(' ')[1])

                            plate_data[plate]['coordinates'].append((x, y - zshift))

                    # Check that data points are unique.
                    # If not, the fitting with splines during the grid generation will fail
                    data=np.array(plate_data[plate]['coordinates'])
                    a,index=np.unique(data,return_index=True,axis=0)
                    index.sort()
                    plate_data[plate]['coordinates']=[(x,y) for x,y in data[index]]

                    temp_max = temp_min = plate_data[plate]['coordinates'][0]



                    # Determine the min/max psi values that lie on the current target plate.
                    for (r, z) in plate_data[plate]['coordinates']:
                        curr_v = self.efit_psi.get_psi(r, z)
                        min_v_test = self.efit_psi.get_psi(temp_min[0], temp_min[1])
                        max_v_test = self.efit_psi.get_psi(temp_max[0], temp_max[1])
                        temp_min = (r, z) if curr_v <= min_v_test else temp_min
                        temp_max = (r, z) if curr_v >= max_v_test else temp_max

                    plate_data[plate]['psi_min_rz'] = temp_min
                    plate_data[plate]['psi_max_rz'] = temp_max

                except Exception as e:
                    print(repr(e))
                    raise ValueError("Something went wrong when processing the target plate file:{}. \
                                     No entry is expected in the yaml file when no file is required".format(target_plates[plate]['file']))

            else:
                print('No file provided for plate:{} '.format(plate))

        # Save plate_data dictionary within the Ingrid object.
        self.plate_data = plate_data

    def order_target_plates(self):
        """
        Ensures the target plate points are oriented clockwise around the
        magnetic axis.

        @author: garcia299
        """
        def unit_vector(v):
            """
            Returns unit vector
            """
            return v / np.linalg.norm(v)

        def angle_between(u, v):
            """
            Compute angle in radians between vectors u and v
            """
            u_norm = unit_vector(u)
            v_norm = unit_vector(v)
            return np.arctan2(u_norm[0] * v_norm[1] - u_norm[1] * v_norm[0], \
                u_norm[0] * v_norm[0] + u_norm[0] * v_norm[1])


        for k in self.plate_data.keys():

            try:
                if self.plate_data[k]['coordinates'] == []:
                    continue
            except:
                continue

            plate = Line([Point(p) for p in self.plate_data[k]['coordinates']])

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

                self.AutoRefineMagAxis()
                v_reference = np.array( [ self.yaml['grid_params']['rmagx'], 
                    self.yaml['grid_params']['zmagx'] ])

                v_start = np.array( [ start.x, start.y ] ) - v_reference
                v_end = np.array( [ end.x, end.y ] ) - v_reference

                if angle_between(v_start, v_end) <= 0:
                    ordered_plate = plate.copy()
                else:
                    ordered_plate = plate.reverse_copy()

            # Gather raw data
            self.plate_data[k]['coordinates'] = ordered_plate.points()

    def plot_target_plates(self):
        """
        Plot the plate_data stored within the Ingrid object to the current figure.

        @author: watkins35, garcia299
        """

        try:
            for plate in self.plate_data:
                try:
                    if self.plate_data[plate]:
                        coor = np.array(self.plate_data[plate]['coordinates'])
                        plt.plot(coor[:, 0], coor[:, 1], label=plate)
                        plt.draw()
                except:
                    continue
        except:
            pass

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
        self.efit_psi.clear_plot()
        self.efit_psi.plot_data(self.yaml['grid_params']['nlevs'])

    def FindMagAxis(self,x:float,y:float)->None:
        sol = root(self.efit_psi.Gradient,[x,y])
        self.yaml['grid_params']['rmagx']=sol.x[0]
        self.yaml['grid_params']['zmagx']=sol.x[1]

    def FindXPoint(self,x:float,y:float)->None:
        sol = root(self.efit_psi.Gradient,[x,y])
        self.yaml['grid_params']['rxpt']=sol.x[0]
        self.yaml['grid_params']['zxpt']=sol.x[1]

    def FindXPoint2(self,x:float,y:float)->None:
        sol = root(self.efit_psi.Gradient,[x,y])
        self.yaml['grid_params']['rxpt2']=sol.x[0]
        self.yaml['grid_params']['zxpt2']=sol.x[1]

    def AutoRefineMagAxis(self):
        self.FindMagAxis(self.yaml['grid_params']['rmagx'],self.yaml['grid_params']['zmagx'])

    def AutoRefineXPoint(self):
        self.FindXPoint(self.yaml['grid_params']['rxpt'],self.yaml['grid_params']['zxpt'])

    def AutoRefineXPoint2(self):
        self.FindXPoint2(self.yaml['grid_params']['rxpt2'],self.yaml['grid_params']['zxpt2'])

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
                                  self.efit_psi.rmagx, self.efit_psi.zmagx, name='psi norm')
        psi = self.efit_psi.v
        psi_magx = self.efit_psi.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.efit_psi.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx))/(psi_xpt1 - psi_magx)
        self.psi_norm.set_v(psinorm)
        self.psi_norm.Calculate_PDeriv()

    def plot_psinorm(self,interactive=True):
        """
        Plot the psi_norm data stored in the Ingrid object.
        """
        self.psi_norm.plot_data(self.yaml['grid_params']['nlevs'],interactive)
        self.plot_target_plates()

    def find_psi_lines(self, tk_controller = None):
        self.psi_finder = RootFinder(self.psi_norm, mode = 'psi_finder', controller = tk_controller)

    def init_LineTracing(self, refresh = False,interactive=True):
        """
        Initializes the line tracing class for the construction
        of the grid.
        Parameter:
            - refresh : boolean
            Re-initialize the LineTracing class.
        @author: watkins35, garcia299
        """
        try:
            # Check if the LineTracing class has been initialized.
            # Exception will be thrown if it no initialization has occured.
            if not self.eq:
                raise AttributeError

        except AttributeError:
            self.eq = LineTracing(self.psi_norm, self.yaml,interactive=interactive)
        if refresh:
            self.eq = LineTracing(self.psi_norm, self.yaml,interactive=interactive)
        self.eq.disconnect()

    def refresh_LineTracing(self):
        self.init_LineTracing(refresh = True)

    def _classify_gridtype(self):
        """
        Analyze the topology around the primary x-point in order to determine the configuration.
        This is a wrapper for the LineTracing class method SNL_find_NSEW.
        @author: garcia299
        """
        print(" # Analyzing topology....")
        try:
            debug = self.yaml['DEBUG']['visual']['find_NSEW']
        except:
            debug = False

        if self.yaml['grid_params']['num_xpt'] == 1:
            self.eq.SNL_find_NSEW(self.xpt1, self.magx, debug)

        elif self.yaml['grid_params']['num_xpt'] == 2:
            self.eq.DNL_find_NSEW(self.xpt1, self.xpt2, self.magx, debug)

        self.yaml['grid_params']['config'] = self.eq.config
        return self.eq.config

    def _analyze_topology(self,verbose=False):
        if verbose: print(" # Analyzing topology....")
        config = self._classify_gridtype()
        if config in ['LSN', 'USN']:
            ingrid_topology = SNL(self, config)
        elif config in ['DNL']:
            ingrid_topology = DNL(self, config)

        self.current_topology = ingrid_topology

    def _get_configuration(self):
        return self.eq.config

    def export(self, fname = 'gridue'):
        """ Saves the grid as an ascii file """
        if isinstance(self.current_topology, SNL):
            self.WritegridueSNL(self.current_topology.gridue_params, fname)
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
                ymlDict=_yaml_.full_load(File.read_text())
                return ymlDict
            except:
                raise IOError('Cannot load the yml file: '+FileName)

        else:
            print('Current Working Directory:'+os.getcwd())
            raise IOError('Cannot find the file: ' +FileName)

    def CreateSubgrid(self,ShowVertices=False,color=None,RestartScratch=False,NewFig=True,OptionTrace='theta',ExtraSettings={},ListPatches='all'):


        try:
            np_cells = self.yaml['grid_params']['grid_generation']['np_global']
        except KeyError:
            np_cells = 2
            print('yaml file did not contain parameter np_global. Set to default value of 2...')

        try:
            nr_cells = self.yaml['grid_params']['grid_generation']['nr_global']
        except KeyError:
            nr_cells = 2
            print('yaml file did not contain parameter nr_global. Set to default value of 2...')


        print('Value Check for local plates:')
        _i = self.current_topology.yaml['target_plates']
        for plate in _i:
            print('Name: {}\n np_local: {}\n'.format(_i[plate]['name'], _i[plate]['np_local']))
        if NewFig:
            self.FigGrid=plt.figure('INGRID: Grid', figsize=(10,10))
            ax=self.FigGrid.gca()
        self.current_topology.construct_grid(np_cells, nr_cells,ShowVertices=ShowVertices,RestartScratch=RestartScratch,OptionTrace=OptionTrace,ExtraSettings=ExtraSettings,ListPatches=ListPatches)


            #Plotgridue(self.current_topology.gridue_params,ax=ax,Verbose=True,facecolor=None,edgecolor=color)
            #self.current_topology.grid_diagram()
            #self.current_topology.grid_diagram()
        #     return True
        # except Exception as e:
        #     print('=' * 80)
        #     print('ERROR DURING GRID GENERATION: {}'.format(repr(e)))
        #     print('=' * 80 + '\n')

    def SetMagReference(self:object,topology:str='SNL')->None:
        self.magx= (self.yaml['grid_params']['rmagx'], self.yaml['grid_params']['zmagx'])
        self.Xpt = (self.yaml['grid_params']['rxpt'],self.yaml['grid_params']['zxpt'])
        self.xpt1 = self.Xpt

        if topology == 'DNL':
            self.xpt2 = (self.yaml['grid_params']['rxpt2'], self.yaml['grid_params']['zxpt2'])

    def PlotMagReference(self):
        (x,y)=self.xpt1
        plt.plot(x,y,'+',color='red',ms=15)
        (x,y)=self.magx
        plt.plot(x,y,'+',color='magenta',ms=15)
        plt.show()

    def PlotPsiNormBounds(self:object)->None:
        Dic={'psi_max':'blue',
              'psi_min_core':'magenta',
              'psi_max_outer':'blue',
              'psi_max_inner':'blue',
              'psi_min_pf':'green',
              'psi_pf2':'yellow'}
        for k,c in Dic.items():
                print(k)
                print(self.yaml['grid_params'][k])
                self.psi_norm.PlotLevel(self.yaml['grid_params'][k],color=Dic[k],label=k)



        self.psi_norm.PlotLevel(1.0,color='red',label='Separatrix')
        plt.legend()


    def PlotPsiLevel(efit_psi:object,level:float,Label:str='')->None:
        plt.contour(efit_psi.r, efit_psi.z, efit_psi.v, [level], colors = 'red', label = Label)

    def Setup(self,interactive=True):
        # TODO have an option to not connect event on figure in command-line mode (so no tk needed)
        print('### Setup INGRID')
        self.OMFIT_read_psi()
        self.calc_efit_derivs()
        self.read_target_plates()
        self.AutoRefineMagAxis()
        self.AutoRefineXPoint()
        self.SetMagReference()
        self.calc_psinorm()
        self.plot_psinorm(interactive=interactive)
        self.init_LineTracing(interactive=interactive)
        self._analyze_topology()

    def ShowSetup(self):
        self.PlotPsiNormBounds()
        self.plot_target_plates()
        self.PlotMagReference()

    def ConstructPatches(self):
        self.current_topology.construct_patches()
        self.current_topology.patch_diagram()
        self.current_topology.CheckPatches()

def Importgridue(fname:str = 'gridue')->dict:
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



