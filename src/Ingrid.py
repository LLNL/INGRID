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

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder
from geometry import Point, Line, SNL_Patch, DNL_Patch, segment_intersect


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

    def __init__(self, params = {}):

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

        self.yaml_lookup = ['grid_params', 'integrator_params', \
                            'target_plates', 'DEBUG']

        self.default_values_lookup = {'grid_params' : self.default_grid_params, 'integrator_params' : self.default_integrator_params, \
            'target_plates' : self.default_target_plates_params, 'DEBUG' : self.default_DEBUG_params}

        self.process_yaml(params)
        self.print_yaml()

    def process_yaml(self, params):
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
                print('Could not find "{}" in YAML file.'.format(item))
                params[item] = get_default_values(item)
                print('Populated "{}" with default value of "{}".\n'.format(item, params[item]))
                continue

            # Second level entries within YAML file dump.
            for sub_item in self.default_values_lookup[item].keys():
                try:
                    params[item][sub_item]
                except KeyError:
                    print('Could not find "{}/{}" in YAML file.'.format(item, sub_item))
                    params[item][sub_item] = get_default_values(item, sub_item)
                    print('Populated "{}/{}" with default value of "{}".\n'.format(item, sub_item, params[item][sub_item]))
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
            try:
                open(target_plates[plate]['file'])
            except FileNotFoundError:
                continue

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

            with open(target_plates[plate]['file']) as f:
                for line in f:
                    point = line.strip()
                    if point.startswith('#'):
                        # move along to the next iteration
                        # this simulates a comment
                        continue
                    x = float(point.split(',')[0])
                    y = float(point.split(',')[1])
                    plate_data[plate]['coordinates'].append((x, y - zshift))

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

            if debug:
                print('Using target plate "{}": {}'.format(target_plates[plate]['name'], plate_data[plate]))

        # Save plate_data dictionary within the Ingrid object.
        self.plate_data = plate_data

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

    def order_plate_points(self, plate, location = 'LITP', orientation = 'cw'):
        """
        Order the plate points within a Line object in a clockwise or counter-clockwise orientation.
        This is used for patch generation and determining boundaries for Patch objects adjacent to
        target plates.
        Parameters:
            - plate : Line object
            - location: string
            This string value should start with an 'L' or 'U'. These values correspond to the 
            target plate names: 
                'LITP' : (lower-inner target-plate)
                'LOTP' : (lower-outer target-plate)
                'UITP' : (upper-inner target-plate)
                'UOTP' : (upper-outer target-plate)
            - orientation: string
            String values should be 'cw' or 'ccw'.
        @author: garcia299
        """

        # The hard-coded logic is for a plate located above the magnetic axis.
        #
        # loc_sgn = 'location sign'.
        #
        # This sign-coefficient adapts the code to a plate below the magnetic axis.

        loc_sgn = 1 if location[0] == 'U' else -1

        # Assume the plate is already ordered.
        start = plate.p[0]
        end = plate.p[-1]

        # Endpoints on same vertical line.
        if start.x == end.x:
            # If 'end' point above 'start' point.
            if end.y - start.y > 0:
                # Return target plate as is
                return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
            else:
                # Else flip plate orientation.
                return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p

        # Endpoints on same horizontal line.
        elif start.y == end.y:
            # If 'end' point to the right of 'start' point.
            if start.x - end.x > 0:
                # Return target plate as is
                return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
            else:
                # Else flip plate orientation
                return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p

        # Endpoints are on sloped line.
        # Check if 'end' point to the right of 'start' point.
        elif loc_sgn * (end.x - start.x) > 0:
            return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
        else:
            return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p

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
                                  self.efit_psi.rmagx, self.efit_psi.zmagx, name='psi norm')
        psi = self.efit_psi.v
        psi_magx = self.efit_psi.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.efit_psi.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx))/(psi_xpt1 - psi_magx)
        self.psi_norm.set_v(psinorm)
        self.psi_norm.Calculate_PDeriv()

    def plot_psinorm(self):
        """
        Plot the psi_norm data stored in the Ingrid object.
        """
        self.psi_norm.plot_data(self.yaml['grid_params']['nlevs'])
        self.plot_target_plates()

    def find_psi_lines(self, tk_controller = None):
        self.psi_finder = RootFinder(self.psi_norm, mode = 'psi_finder', controller = tk_controller)

    def init_LineTracing(self, refresh = False):
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
            self.eq = LineTracing(self.psi_norm, self.yaml)
        if refresh:
            self.eq = LineTracing(self.psi_norm, self.yaml)
        self.eq.disconnect()

    def refresh_LineTracing(self):
        self.init_LineTracing(refresh = True)

    def catagorize_patches(self):
        """
        Group Patch objects into their relevant topological regions.
        TODO: Create a class in geometry.py for grouping Patches.
        @author: garcia299
        """
        m = self.patch_matrix
        if isinstance(self, SNL):
            self.SOL = m[1][1:-1]
            self.CORE = m[2][2:-2]
            self.PF = [m[2][1], m[2][-2]]
        elif isinstance(self, DNL):
            self.PRIMARY_SOL = m[2][1:4] + m[2][8:11]
            self.SECONDARY_SOL = m[1][1:5] + m[1][7:11]
            self.CORE = m[3][2:4] + m[3][8:10]
            self.PRIMARY_PF = [m[3][1], m[3][-2]]
            self.SECONDARY_PF = [m[3][4], m[2][4], m[3][7], m[2][7]]

    def _classify_gridtype(self):
        """
        Analyze the topology around the primary x-point in order to determine the configuration.
        This is a wrapper for the LineTracing class method SNL_find_NSEW.
        @author: garcia299
        """
        try:
            debug = self.yaml['DEBUG']['visual']['find_NSEW']
        except:
            debug = False
        if self.yaml['grid_params']['num_xpt'] == 1:
            self.eq.SNL_find_NSEW(self.xpt1, self.magx, debug)
        elif self.yaml['grid_params']['num_xpt'] == 2:
            self.eq.DNL_find_NSEW(self.xpt1, self.xpt2, self.magx, debug)
        
        self.yaml['grid_params']['config'] = self.eq.config
        return self.yaml['grid_params']['config']

    def _analyze_topology(self):
        config = self._classify_gridtype()
        if config == 'LSN':
            ingrid_topology = LSN(self)
        elif config == 'USN':
            ingrid_topology = USN(self)
        elif config == 'DNL':
            ingrid_topology = DNL(self)

        self.current_topology = ingrid_topology

    def _get_configuration(self):
        return self.eq.config

    def export(self, fname = 'gridue'):
        """ Saves the grid as an ascii file """
        if isinstance(self.current_topology, SNL):
            self.write_gridue_SNL(self.current_topology.gridue_params, fname)
        elif isinstance(self.current_topology, DNL):
            self.write_gridue_DNL(self.current_topology.gridue_params, fname)

    def write_gridue_SNL(self, gridue_params, fname = 'gridue'):
        
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

    def write_gridue_DNL(self, gridue_params, fname = 'gridue'):
        
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


class DNL(Ingrid):
    def __init__(self, INGRID_object):
        super().__init__(params = INGRID_object.yaml)
        self.efit_psi = INGRID_object.efit_psi
        self.psi_norm = INGRID_object.psi_norm
        self.eq = INGRID_object.eq
        self.plate_data = INGRID_object.plate_data
        self.config = 'DNL'
        print('=' * 80)
        print('DNL Object!')
        print('=' * 80 + '\n')

    def patch_diagram(self):
        """ Generates the patch diagram for a given configuration. """
        
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown', 'c', 
                  'm', 'dodgerblue', 'darkorchid', 'crimson', 
                  'darkorange', 'lightgreen', 'lightseagreen', 'indigo', 
                  'mediumvioletred', 'mistyrose', 'darkolivegreen', 'rebeccapurple']

        try:
            plt.close('INGRID: Patch Map')
        except:
            pass

        plt.figure('INGRID: Patch Map', figsize=(6, 10))
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('DNL Patch Diagram')

        for i in range(len(self.patches)):
            self.patches[i].plot_border('green')
            self.patches[i].fill(colors[i])

        plt.show() 

    def grid_diagram(self):

        try:
            plt.close('INGRID: Grid')
        except:
            pass

        plt.figure('INGRID: Grid', figsize=(6,10))
        for patch in self.patches:
            patch.plot_subgrid()
        
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('INGRID DNL Subgrid')
        plt.show()

    def animate_grid(self):
        """
        animate_grid:
            Sequentially draw each cell contained withing the 
            finalized RM and ZM arrays. Used for verification
            of poloidal/radial index space ordering of cell data.

        Parameters:
            N/A

        Return:
            N/A

        @author: garcia299
        """

        try:
            plt.close('INGRID: Debug')
        except:
            pass
        plt.figure('INGRID: Debug', figsize=(6, 10))
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('visualize gridue')

        k = [1,2,4,3,1]

        for i in range(len(self.rm)):
            for j in range(len(self.rm[0])):
                plt.plot(self.rm[i][j][k], self.zm[i][j][k])
                plt.pause(0.01)

    def concat_grid(self):
        """
        concat_grid:
            Concatenate each refined grid within a Patch object into
            UEDGE formatted arrays RM and ZM.

        Parameters:
            N/A

        Return:
            N/A

        @author: garcia299
        """

        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the 
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches:
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

            print('"{}" has npol = {} and nrad = {}'.format(patch.patchName, patch.npol, patch.nrad))

        # Total number of poloidal indices in all subgrids.
        np_total1 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:5]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total1 = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        # Total number of poloidal indices in all subgrids.
        np_total2 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:11]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total2 = int(np.sum([patch[7].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        rm1 = np.zeros((np_total1, nr_total1, 5), order = 'F')
        zm1  = np.zeros((np_total1, nr_total1, 5), order = 'F')
        rm2 = np.zeros((np_total2, nr_total2, 5), order = 'F')
        zm2  = np.zeros((np_total2, nr_total2, 5), order = 'F')

        ixcell = 0
        jycell = 0

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')
        for ixp in range(1, 5):
            
            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch. 
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                        print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')

        ixcell = 0
        jycell = 0

        for ixp in range(7, 11):
            
            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch. 
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                        print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm1)):
            rm1[i] = rm1[i][::-1]
        for i in range(len(zm1)):
            zm1[i] = zm1[i][::-1]
        for i in range(len(rm2)):
            rm2[i] = rm2[i][::-1]
        for i in range(len(zm2)):
            zm2[i] = zm2[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb1 = len(rm1) - 2
        ixlb1 = 0
        ixrb2 = len(rm2) - 2
        ixlb2 = 0

        rm1 = self.add_guardc(rm1, ixlb1, ixrb1)
        zm1 = self.add_guardc(zm1, ixlb1, ixrb1)
        rm2 = self.add_guardc(rm2, ixlb2, ixrb2)
        zm2 = self.add_guardc(zm2, ixlb2, ixrb2)

        self.rm = np.concatenate((rm1, rm2))
        self.zm = np.concatenate((zm1, zm2))

        try:
            debug = self.yaml['DEBUG']['visual']['gridue']
        except:
            debug = False
        
        if debug:
            self.animate_grid()

    def add_guardc(self, cell_map, ixlb, ixrb, eps = 1e-3):
        """
        add_guardc:
            Method for adding guard cells to a refined patch map.

        Parameters:
            cell_map : list-like
                - The refined patch map with shape = (np, nr, 5).
            ixlb : int
                - Index corresponding to the left boundary of the
                  cell_map object.
            ixrb : int
                - Index corresponding to the right boundary of the
                  cell_map object.
            eps : float
                - Epsilon value for guard cell size. For usage in
                  linear interpolation of non-guard cells.

        Returns:
            cell_map : list-like
                - The modified cell_map parameter with guard cells
                  included.
        """

        def set_guard(cell_map, ix, iy, eps, boundary):
            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # TODO: Edit the code to reflect this at some point so the next reader is not overwhelmed.
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
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'top')

        return cell_map

    def set_gridue(self):
        """
        set_gridue:
            Prepares 'self.gridue_params' dictionary with required data.
            The self.gridue_params attribute is used to write a gridue
            formatted file

        Parameters:
            N/A
        Return:
            N/A
        """

        ixlb = 0
        ixrb = len(self.rm) - 2

        nxm = len(self.rm) - 2
        nym = len(self.rm[0]) - 2
        iyseparatrix1 = self.patch_lookup['C2'].nrad - 1
        iyseparatrix2 = self.patch_lookup['B5'].nrad + self.patch_lookup['C5'].nrad - 2
        ix_plate1 = 0
        ix_cut1 = self.patch_lookup['A1'].npol - 1
        ix_cut2 = self.patch_lookup['A1'].npol + self.patch_lookup['A2'].npol + self.patch_lookup['A3'].npol - 3
        ix_plate2 = ix_cut2 + self.patch_lookup['A4'].npol - 1
        iyseparatrix3 = iyseparatrix2
        iyseparatrix4 = iyseparatrix1
        ix_plate3 = ix_plate2 + 2
        ix_cut3 = ix_plate3 + self.patch_lookup['A5'].npol - 1
        ix_cut4 = ix_cut3 + self.patch_lookup['A6'].npol + self.patch_lookup['A7'].npol - 2
        ix_plate4 = ix_cut4 + self.patch_lookup['A8'].npol - 1

        psi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        br = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bz = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bpol = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bphi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        b = np.zeros((nxm + 2, nym + 2, 5), order = 'F')

        rm = self.rm
        zm = self.zm
        rb_prod = self.efit_psi.rcenter * self.efit_psi.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.efit_psi.get_psi(_r, _z)
                    _br = self.efit_psi.get_psi(_r, _z, tag = 'vz') / _r
                    _bz = -self.efit_psi.get_psi(_r, _z, tag = 'vr') / _r
                    _bpol = np.sqrt(_br ** 2 + _bz ** 2)
                    _bphi = rb_prod / _r
                    _b = np.sqrt(_bpol ** 2 + _bphi ** 2)

                    psi[i][j][k] = _psi
                    br[i][j][k] = _br
                    bz[i][j][k] = _bz
                    bpol[i][j][k] = _bpol
                    bphi[i][j][k] = _bphi
                    b[i][j][k] = _b

        self.gridue_params = {'nxm' : nxm, 'nym' : nym, 'iyseparatrix1' : iyseparatrix1, 'iyseparatrix2' : iyseparatrix2, \
                'ix_plate1' : ix_plate1, 'ix_cut1' : ix_cut1, 'ix_cut2' : ix_cut2, 'ix_plate2' : ix_plate2, 'iyseparatrix3' : iyseparatrix3, \
                'iyseparatrix4' : iyseparatrix4, 'ix_plate3' : ix_plate3, 'ix_cut3' : ix_cut3, 'ix_cut4' : ix_cut4, 'ix_plate4' : ix_plate4, \
                'rm' : self.rm, 'zm' : self.zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b, '_FILLER_' : -1}


    def construct_grid(self, np_cells = 3, nr_cells = 3):
        """
        construct_grid:
            Method for facilitating the refining of Patches within the Ingrid object.
            Grid construction as well as grid concatenation and gridue
            file prep is completed here.

        Parameters:
            np_cells : int
                - Default number of poloidal cells to generate in a Patch
            nr_cells : int
                - Default number of radial cells to generate in a Patch

        # TODO: Remove parameters and rely on np_global and nr_global in YAML
        """

        # Check for visual mode in Patch refinement.
        try:
            visual = self.yaml['DEBUG']['visual']['subgrid']
        except:
            visual = False

        # Check for verbose output in Patch refinement.
        try:
            verbose = self.yaml['DEBUG']['verbose']['grid_generation']
        except:
            verbose = False

        # Check for specified number of poloidal cells for patches in Primary SOL.
        # Uses np_global if no specification.
        try:
            np_primary_sol = self.yaml['grid_params']['grid_generation']['np_primary_sol']
        except:
            np_primary_sol = self.yaml['grid_params']['grid_generation']['np_global']

        # Check for specified number of poloidal cells for patches in Secondary SOL.
        # Uses np_global if no specification.
        try:
            np_secondary_sol = self.yaml['grid_params']['grid_generation']['np_secondary_sol']
        except:
            np_secondary_sol = self.yaml['grid_params']['grid_generation']['np_global']

        # Check for specified number of poloidal cells for patches in Core.
        # Uses np_global if no specification.
        try:
            np_core = self.yaml['grid_params']['grid_generation']['np_core']
        except:
            np_core = self.yaml['grid_params']['grid_generation']['np_global']

        # Check for specified number of poloidal cells for patches in Primary PF.
        # Uses np_global if no specification.
        try:
            np_primary_pf = self.yaml['grid_params']['grid_generation']['np_primary_pf']
        except:
            np_primary_pf = self.yaml['grid_params']['grid_generation']['np_global']

        # Check for specified number of poloidal cells for patches in Secondary PF.
        # Uses np_global if no specification.
        try:
            np_secondary_pf = self.yaml['grid_params']['grid_generation']['np_secondary_pf']
        except:
            np_secondary_pf = self.yaml['grid_params']['grid_generation']['np_global']

        # Check for specified number of radial cells for patches in Primary SOL.
        # Uses nr_global if no specification.
        try:
            nr_primary_sol = self.yaml['grid_params']['grid_generation']['nr_primary_sol']
        except:
            nr_primary_sol = self.yaml['grid_params']['grid_generation']['nr_global']

        # Check for specified number of radial cells for patches in Secondary SOL.
        # Uses nr_global if no specification.
        try:
            nr_secondary_sol = self.yaml['grid_params']['grid_generation']['nr_secondary_sol']
        except:
            nr_secondary_sol = self.yaml['grid_params']['grid_generation']['nr_global']

        # Check for specified number of radial cells for patches in Core.
        # Uses nr_global if no specification.
        try:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_core']
        except:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_global']

        # Check for specified number of radial cells for patches in Primary PF.
        # Uses nr_global if no specification.
        try:
            nr_primary_pf = self.yaml['grid_params']['grid_generation']['nr_primary_pf']
        except:
            nr_primary_pf = self.yaml['grid_params']['grid_generation']['nr_global']

        # Check for specified number of radial cells for patches in Secondary PF.
        # Uses nr_global if no specification.
        try:
            nr_secondary_pf = self.yaml['grid_params']['grid_generation']['nr_secondary_pf']
        except:
            nr_secondary_pf = self.yaml['grid_params']['grid_generation']['nr_global']


        # Local lookup table helper.
        lookup = {self.PRIMARY_SOL : [nr_primary_sol, np_primary_sol, 'PRIMARY_SOL'], \
                    self.SECONDARY_SOL : [nr_secondary_sol, np_secondary_sol, 'SECONDARY_SOL'], \
                    self.CORE : [nr_core, np_core, 'CORE'], \
                    self.PRIMARY_PF : [nr_primary_pf, np_primary_pf, 'PRIMARY_PF'], \
                    self.SECONDARY_PF : [nr_secondary_pf, np_secondary_pf, 'SECONDARY_PF']}

        # Iterate through patches and refine into grids.
        for patch in self.patches:
            for key in lookup.keys():
                if patch not in key:
                    continue
                nr_cells = lookup[key][0]
                np_cells = lookup[key][1]
                topos_group = lookup[key][2]
            print('=====================================')
            print('=====================================')
            print('CONSTRUCTING GRID FOR PATCH: {}'.format(patch.patchName))
            print('-------------------------------------')
            print('Patch "{}" is located in "{}"'.format(patch.patchName, topos_group))
            print('-------------------------------------')
            print('Number of poloidal cells to generate: {}'.format(np_cells))
            print('-------------------------------------')
            print('Number of radial cells to generate: {}'.format(nr_cells))
            print('=====================================')
            print('=====================================')
            patch.make_subgrid(self, np_cells, nr_cells, verbose = verbose, visual = visual)

            # Tidy up primary x-point
            primary_xpt = Point(self.xpt1)
            secondary_xpt = Point(self.xpt2)

            if patch.patchName == 'B1':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'C1':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'B2':
                patch.adjust_corner(primary_xpt, 'SW')
            elif patch.patchName == 'C2':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'C7':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'B7':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'C8':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'B8':
                patch.adjust_corner(primary_xpt, 'SW')

            # Tidy up secondary x-point
            elif patch.patchName == 'A3':
                patch.adjust_corner(secondary_xpt, 'SE')
            elif patch.patchName == 'B3':
                patch.adjust_corner(secondary_xpt, 'NE')
            elif patch.patchName == 'A4':
                patch.adjust_corner(secondary_xpt, 'SW')
            elif patch.patchName == 'B4':
                patch.adjust_corner(secondary_xpt, 'NW')
            elif patch.patchName == 'B5':
                patch.adjust_corner(secondary_xpt, 'NE')
            elif patch.patchName == 'A5':
                patch.adjust_corner(secondary_xpt, 'SE')
            elif patch.patchName == 'B6':
                patch.adjust_corner(secondary_xpt, 'NW')
            elif patch.patchName == 'A6':
                patch.adjust_corner(secondary_xpt, 'SW')

        self.concat_grid()
        self.set_gridue()

    def construct_patches(self):
        """
        construct_patches:
            Draws lines and creates patches based off the X-points, Magnetic Axis, 
            and specified Psi surfaces for the DNL configuration.

            Patch Labeling Key:
                A-C : Radial location of a Patch. Convention here has A being the outermost 
                      and C being the innermost.
                1-8 : Poloidal location of a Patch. Convention here follows the clockwise 
                      orientation.

        Parameters: 
            No user specified parameters.

        Returns:
            No return value.
        """

        def order_plate_points(plate, location = 'UITP', orientation = 'cw'):
            """
            order_plate_points:
                Sets the points in the target plate to have an orientation
                increasing in R and Z. Checks the endpoints to satisfy this criteria.

            Parameters:
                plate : Line object
                    - Line object representing the target plate geometry to order.
                location : str
                    - String indicating the target plate location. 
                      Must be:
                        'UITP' : Upper Inner Target Plate
                        'UOTP' : Upper Outer Target Plate
                        'LITP' : Lower Inner Target Plate
                        'LOTP' : Lower Outer Target Plate
                orientation : str
                    - Orientation to order plate points in. 
                      Must be:
                        'cw' : Clockwise
                        'ccw' : Counter Clockwise
            Returns:
                List
                    - List of Point objects ordered in the orientation specified.
            """

            loc_sgn = 1 if location[0] == 'U' else -1

            start = plate.p[0]
            end = plate.p[-1]
            # Endpoints on same vertical line.
            if start.x == end.x:
                # If rhs endpoint above lhs endpoint.
                if end.y - start.y > 0:
                    # Return target plate as is
                    return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
                else:
                    # Else flip plate orientation.
                    return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p
            # Endpoints on same horizontal line.
            elif start.y == end.y:
                # If lhs endpoint to the left of rhs endpoint.
                if end.x - start.x > 0:
                    # Return target plate as is
                    return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
                else:
                    # Else flip plate orientation
                    return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p
            # Endpoints are on sloped line.
            # Check if lhs endpoint is on the left of rhs endpoint
            elif loc_sgn * (end.x - start.x) > 0:
                return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
            else:
                return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p


        # Check for visual mode in Patch Map generation.
        try:
            visual = self.yaml['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False

        # Check for verbose mode in Patch Map generation.
        try:
            verbose = self.yaml['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False

        # Check for specified radian adjustment for INNER-most horizontal line
        # through magnetic axis.
        try:
            inner_tilt = self.yaml['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0

        # Check for specified radian adjustment for OUTER-most horizontal line
        # through magnetic axis.
        try:
            outer_tilt = self.yaml['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0

        # Ingrid plate attributes as references to geometric plate data.
        """
        Plate Labeling:
            W : Inner
            E : Outer
            1 : Closest to primary X-point
            2 : Closest to secondary X-point
        """
        self.plate_W1 = self.plate_data['plate_W1']['coordinates']
        self.plate_E1 = self.plate_data['plate_E1']['coordinates']
        self.plate_E2 = self.plate_data['plate_E2']['coordinates']
        self.plate_W2 = self.plate_data['plate_W2']['coordinates']

        # References to NSEW rz coordinates for xpt1 and xpt2.
        xpt1_dict = self.eq.eq_psi['xpt1']
        xpt2_dict = self.eq.eq_psi['xpt2']

        # References to NSEW theta values for xpt1 and xpt2.
        xpt1_theta = self.eq.eq_psi_theta['xpt1']
        xpt1_theta = self.eq.eq_psi_theta['xpt2']

        # Psi value corresponding to separatrix through xpt1
        sptrx1_v = self.psi_norm.get_psi(self.xpt1[0], self.xpt1[1])
        # Psi value corresponding to separatrix through xpt2
        sptrx2_v = self.psi_norm.get_psi(self.xpt2[0], self.xpt2[1])

        # RZ oordinates for magnetic axis with the corresponding r-shift and z-shift applied.
        magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
            self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

        # Get relevant psi surface values.
        psi_max = self.yaml['grid_params']['psi_max']
        psi_min_core = self.yaml['grid_params']['psi_min_core']
        psi_min_pf = self.yaml['grid_params']['psi_min_pf']
        psi_max_outer = self.yaml['grid_params']['psi_max_outer']
        psi_max_inner = self.yaml['grid_params']['psi_max_inner']
        psi_min_pf_2 = self.yaml['grid_params']['psi_pf2']

        # Process raw RZ coordinates of target plate data into Line objects.
        LITP = Line(order_plate_points(Line([Point(i) for i in self.plate_W1]), location = 'LITP'))
        LOTP = Line(order_plate_points(Line([Point(i) for i in self.plate_E1]), location = 'LOTP'))
        UITP = Line(order_plate_points(Line([Point(i) for i in self.plate_E2]), location = 'UITP'))
        UOTP = Line(order_plate_points(Line([Point(i) for i in self.plate_W2]), location = 'UOTP'))

        # Generate INNER Horizontal Mid-Plane line and apply inner-tilt radian adjustment (if any)
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        inner_midLine = Line([LHS_Point, RHS_Point])

        # Generate OUTER Horizontal Mid-Plane line and apply outer-tilt radian adjustment (if any)
        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        outer_midLine = Line([LHS_Point, RHS_Point])

        # Generate Vertical Mid-Plane line that intersects the secondary x-pt and the magnetic axis
        # TODO : Maybe create new criteria for this 'topLine'?
        slp = [self.xpt2[0] - magx[0], self.xpt2[1] - magx[1]]
        slp = slp[1] / slp[0]
        topLine_tilt = np.arctan(slp)
        Upper_Point = Point(-1e6, slp * (-1e6 - self.xpt2[0]) + self.xpt2[1])
        Lower_Point = Point(1e6, slp * (1e6 - self.xpt2[0]) + self.xpt2[1])
        topLine = Line([Lower_Point, Upper_Point])

        # ====================================================================================
        # Beginning of line-tracing for primary separatrix and core region.
        # ====================================================================================
        """
        Draw line from primary x-point towards the core. 
        Stops at psi surface corresponding to psiMinCore value.

            Starting tag : xpt1N
            Ending Tag   : psiMinCore
        """
        xpt1N__psiMinCore = self.eq.draw_line(xpt1_dict['N'], {'psi' : psi_min_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from primary x-point towards the inner mid-line. 
        Stops when intersection occurs with inner_midLine object.

            Starting tag : xpt1NW
            Ending tag   : sptrx1imidLine
        """
        xpt1NW__sptrx1imidLine = self.eq.draw_line(xpt1_dict['NW'], {'line' : inner_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from primary x-point towards the outer mid-line. 
        Stops when intersection occurs with outer_midLine object.

            Starting tag : xpt1NE
            Ending tag   : sptrx1omidLine
        """
        xpt1NE__sptrx1omidLine = self.eq.draw_line(xpt1_dict['NE'], {'line' : outer_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from inner_midLine on primary separatrix to topLine. 
        Stops when intersection occurs with topLine object.

            Starting tag : sptrx1imidLine
            Ending tag   : topLine
        """
        sptrx1imidLine__topLine = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from outer_midLine on primary separatrix to topLine.
        Stops when intersection occurs with topLine object.

            Starting tag : sptrx1omidLine
            Ending tag   : topLine
        """
        sptrx1omidLine__topLine = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line in the core region from psiMinCore to inner_midLine.
        Stops when intersection occurs with inner_midLine object.

            Starting tag : psiMinCore
            Ending tag   : imidLineCore
        """
        psiMinCore__imidLineCore = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : inner_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line in the core region from psiMinCore to outer_midLine.
        Stops when intersection occurs with outer_midLine object.

            Starting tag : psiMinCore
            Ending tag   : omidLineCore
        """
        psiMinCore__omidLineCore = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : outer_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from inner_midLine in core region to topLine.
        Stops when intersection occurs with topLine object.

            Starting tag : imidLineCore
            Ending tag   : topLine
        """
        imidLineCore__topLine = self.eq.draw_line(psiMinCore__imidLineCore.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from outer_midLine in core region to topLine.
        Stops when intersection occurs with topLine object.

            Starting tag : omidLineCore
            Ending tag   : topLine
        """
        omidLineCore__topLine = self.eq.draw_line(psiMinCore__omidLineCore.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        # ====================================================================================
        # End of line-tracing for primary separatrix and core region.
        # ====================================================================================
        # ####################################################################################

        # Create line objects representing the primary separatrix and the core region boundaries.
        sptrx1 = Line([xpt1NW__sptrx1imidLine.p + sptrx1imidLine__topLine.p + sptrx1omidLine__topLine.p[-2::-1] + xpt1NE__sptrx1omidLine.p[-2::-1]][0])
        core = Line([psiMinCore__imidLineCore.p + imidLineCore__topLine.p + omidLineCore__topLine.p[-2::-1] + psiMinCore__omidLineCore.p[-2::-1]][0])
        # Note: These will handy for termination criteria in later line tracing.

        # ####################################################################################
        # ====================================================================================
        # Beginning of line-tracing for primary private-flux region.
        # ====================================================================================
        """
        Draw line from primary x-point towards the primary private-flux region. 
        Stops at psi surface corresponding to psiMinPF1 value.

            Starting tag : xpt1
            Ending tag   : psiMinPF1
        """
        xpt1__psiMinPF1 = self.eq.draw_line(xpt1_dict['S'], {'psi' : psi_min_pf}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from primary private-flux surface to lower inner target-plate. 
        Stops upon intersection with LITP Line object.

            Starting tag : psiMinPF1
            Ending tag   : LITP
        """
        psiMinPF1__LITP = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : LITP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from primary private-flux surface to lower outer target-plate. 
        Stops upon intersection with LOTP Line object.

            Starting tag : psiMinPF1
            Ending tag   : LOTP
        """
        psiMinPF1__LOTP = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : LOTP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from primary x-point to lower inner target-plate. 
        Stops upon intersection with LITP Line object

            Starting tag : xpt1
            Ending tag   : LITP
        """
        xpt1__LITP = self.eq.draw_line(xpt1_dict['SW'], {'line' : LITP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from primary x-point to lower outer target-plate. 
        Stops upon intersection with LOTP Line object

            Starting tag : xpt1
            Ending tag   : LOTP
        """
        xpt1__LOTP = self.eq.draw_line(xpt1_dict['SE'], {'line' : LOTP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        # ====================================================================================
        # End of line-tracing for primary private-flux region.
        # ====================================================================================
        # ####################################################################################

        # ####################################################################################
        # ====================================================================================
        # Beginning of line-tracing for secondary separatrix.
        # ====================================================================================
        """
        Draw line from secondary x-point to primary separatrix.
        Stops upon intersection with sptrx1 Line object.

            Starting tag : xpt2
            Ending tag   : sptrx1
        """
        xpt2__sptrx1 = self.eq.draw_line(xpt2_dict['N'], {'line' : (sptrx1, topLine_tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from secondary x-point to upper inner target-plate. 
        Stops when intersection occurs with UITP Line object.

            Starting tag : xpt2
            Ending tag   : UITP
        """
        xpt2__UITP = self.eq.draw_line(xpt2_dict['SE'], {'line' : UITP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from secondary x-point to upper outer target-plate.
        Stops when intersection occurs with UOTP Line oject.

            Starting tag : xpt2
            Ending tag   : UOTP
        """
        xpt2__UOTP = self.eq.draw_line(xpt2_dict['SW'], {'line' : UOTP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from secondary x-point to the inner mid-line.
        Stops upon intersection with inner_midLine Line object.

            Starting tag : xpt2NE
            Ending tag   : sptrx2imidLine
        """
        xpt2NE__sptrx2imidLine = self.eq.draw_line(xpt2_dict['NE'], {'line' : inner_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from secondary x-point to the outer mid-line.
        Stops upon intersection with outer_midLine Line object.

            Starting tag : xpt2NW
            Ending tag   : sptrx2omidLine
        """
        xpt2NW__sptrx2omidLine = self.eq.draw_line(xpt2_dict['NW'], {'line' : outer_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from inner_midLine on secondary separatrix to lower inner target-plate. 
        Stops when intersection occurs with LITP Line object.

            Starting tag : sptrx2imidLine
            Ending tag   : LITP
        """
        sptrx2imidLine__LITP = self.eq.draw_line(xpt2NE__sptrx2imidLine.p[-1], {'line' : LITP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from outer_midLine on secondary separatrix to lower outer target-plate. 
        Stops when intersection occurs with LOTP Line object.

            Starting tag : sptrx2omidLine
            Ending tag   : LOTP
        """
        sptrx2omidLine__LOTP = self.eq.draw_line(xpt2NW__sptrx2omidLine.p[-1], {'line' : LOTP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        # ====================================================================================
        # End of line-tracing for secondary separatrix.
        # ====================================================================================
        # ####################################################################################

        # Saving portion of secondary separatrix on the inner and outer sides as Line objects
        # Note: These will be handy for termination criteria in later line tracing.
        sptrx2_inner = Line([xpt2NE__sptrx2imidLine.p + sptrx2imidLine__LITP.p][0])
        sptrx2_outer = Line([xpt2NW__sptrx2omidLine.p + sptrx2omidLine__LOTP.p][0])

        # ####################################################################################
        # ====================================================================================
        # Beginning of line-tracing for secondary private-flux region.
        # ====================================================================================
        """
        Draw line from secondary x-point towards secondary private-flux region.
        Stops upon reaching psi-surface corresponding to psi_min_pf_2 psi value.

            Starting tag : xpt2
            Ending tag   : psiMinPF2
        """
        xpt2__psiMinPF2 = self.eq.draw_line(xpt2_dict['S'], {'psi' : psi_min_pf_2}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        """
        Split the above Line object into two separate Line objects.
        Note: Splitting of Line object is required to conform with final index space mappings.
        
        Labeling key:
            xpt2__psiMinPF2_A : The first half of the split line
                i.e. Starting at xpt2__psiMinPF2.p[0] and ending at mid-point.
            xpt2__psiMinPF2_B : The second half of the split line
                i.e. Starting at mid-point and ending at xpt2__psiMinPF2.p[-1]
        """
        xpt2__psiMinPF2_A, xpt2__psiMinPF2_B = xpt2__psiMinPF2.split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p)//2], add_split_point = True)

        """
        Draw line from psiMinPF2_A to upper inner target-plate. 
        Stops when intersection occurs with UITP Line object.

            Starting tag : psiMinPF2_A
            Ending tag   : UITP
        """
        psiMinPF2_A__UITP = self.eq.draw_line(xpt2__psiMinPF2_B.p[0], {'line' : UITP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from psiMinPF2_A to upper outer target-plate. 
        Stops when intersection occurs with UOTP Line object.

            Starting tag : psiMinPF2_A
            Ending tag   : UOTP
        """
        psiMinPF2_A__UOTP = self.eq.draw_line(xpt2__psiMinPF2_B.p[0], {'line' : UOTP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from psiMinPF2_B to upper inner target-plate. 
        Stops when intersection occurs with UITP Line object.

            Starting tag : psiMinPF2_B
            Ending tag   : UITP
        """
        psiMinPF2_B__UITP = self.eq.draw_line(xpt2__psiMinPF2_B.p[-1], {'line' : UITP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose, debug = True)
        """
        Draw line from psiMinPF2_B to upper outer target-plate. 
        Stops when intersection occurs with UOTP Line object.

            Starting tag : psiMinPF2_B
            Ending tag   : UOTP
        """
        psiMinPF2_B__UOTP = self.eq.draw_line(xpt2__psiMinPF2_B.p[-1], {'line' : UOTP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        # ====================================================================================
        # End of line-tracing for secondary private-flux region.
        # ====================================================================================
        # ####################################################################################

        # ####################################################################################
        # ====================================================================================
        # Beginning of line-tracing for outer.
        # ====================================================================================
        """
        Draw line from secondary x-point to psi-surface corresponding to value of psi_max_inner.
        Stops line tracing when integrator reaches value of psi_max_inner.

            Starting tag : xpt2
            Ending tag   : psiMaxInner
        
        Check if 'use_secondary_NE' is set to True. This allows for NE line to be drawn tangent
        to secondary separatrix. The radian value of 'secondary_NE_adjust' allows for theta
        adjustment of the generated line.
        """
        if self.yaml['grid_params']['patch_generation']['use_secondary_NE']:
            v = np.array([xpt2NE__sptrx2imidLine.p[1].x - self.xpt2[0], xpt2NE__sptrx2imidLine.p[1].y - self.xpt2[1]])
            tilt = np.arccos(np.dot(v, np.array([-1, 0])) / np.linalg.norm(v)) \
                    + self.yaml['grid_params']['patch_generation']['secondary_NE_adjust']

            xpt2__psiMaxInner = self.eq.draw_line(xpt2_dict['E'], {'psi_horizontal' : (psi_max_inner, tilt)}, \
                option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        # If option is not set (default behavior), line trace orthogonal to psi surface.
        else:
            xpt2__psiMaxInner = self.eq.draw_line(xpt2_dict['E'], {'psi' : psi_max_inner}, \
                option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from secondary x-point to psi-surface corresponding to value of psi_max_outer.
        Stops line tracing when integrator reaches value of psi_max_outer.

            Starting tag : xpt2
            Ending tag   : psiMaxOuter
        
        Check if 'use_secondary_NW' is set to True. This allows for NW line to be drawn tangent
        to secondary separatrix. The radian value of 'secondary_NW_adjust' allows for theta
        adjustment of the generated line.
        """
        if self.yaml['grid_params']['patch_generation']['use_secondary_NW']:
            v = np.array([xpt2NW__sptrx2omidLine.p[1].x - self.xpt2[0], xpt2NW__sptrx2omidLine.p[1].y - self.xpt2[1]])
            tilt = -np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) \
                    + self.yaml['grid_params']['patch_generation']['secondary_NW_adjust']

            xpt2__psiMaxOuter = self.eq.draw_line(xpt2_dict['W'], {'psi_horizontal' : (psi_max_outer, tilt)}, \
                option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        # If option is not set (default behavior), line trace orthogonal to psi surface.
        else:
            xpt2__psiMaxOuter = self.eq.draw_line(xpt2_dict['W'], {'psi' : psi_max_outer}, \
                option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line along psi-surface corresponding to value of 'psi_max_inner' towards upper
        inner target-plate.
        Stops line-tracing when integrator intersects UITP Line object.

            Starting tag : psiMaxInner
            Ending tag   : UITP
        """
        psiMaxInner__UITP = self.eq.draw_line(xpt2__psiMaxInner.p[-1], {'line' : UITP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line along psi-surface corresponding to value of 'psi_max_outer' towards upper
        outer target-plate.
        Stops line-tracing when integrator intersects UOTP Line object.

            Starting tag : psiMaxOuter
            Ending tag   : UOTP
        """
        psiMaxOuter__UOTP = self.eq.draw_line(xpt2__psiMaxOuter.p[-1], {'line' : UOTP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line along psi-surface corresponding to value of 'psi_max_inner' towards
        inner mid-line.
        Stops line-tracing when integrator intersects inner_midLine Line object.

            Starting tag : psiMaxInner
            Ending tag   : imidLine
        """
        psiMaxInner__imidLine = self.eq.draw_line(xpt2__psiMaxInner.p[-1], {'line' : inner_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line along psi-surface corresponding to value of 'psi_max_outer' towards
        outer mid-line.
        Stops line-tracing when integrator intersects outer_midLine Line object.

            Starting tag : psiMaxOuter
            Ending tag   : omidLine
        """
        psiMaxOuter__omidLine = self.eq.draw_line(xpt2__psiMaxOuter.p[-1], {'line' : outer_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        """
        Draw line from inner_midLine on psi-surface of 'psi_max_inner' towards lower
        inner target-plate.
        Stops when intersection occurs with LITP Line object.

            Starting tag : imidLine
            Ending tag   : LITP
        """
        imidLine__LITP = self.eq.draw_line(psiMaxInner__imidLine.p[-1], {'line' : LITP}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        """
        Draw line from outer_midLine on psi-surface of 'psi_max_outer' towards lower
        outer target-plate.
        Stops when intersection occurs with LOTP Line object.

            Starting tag : omidLine
            Ending tag   : LOTP
        """
        omidLine__LOTP = self.eq.draw_line(psiMaxOuter__omidLine.p[-1], {'line' : LOTP}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        sptrx3_inner = Line([psiMaxInner__imidLine.p + imidLine__LITP.p][0])
        sptrx3_outer = Line([psiMaxOuter__omidLine.p + omidLine__LOTP.p][0])

        # Computing angle between line tangent to sptrx1 in NW direction and unit vector < 1, 0 >
        v = np.array([xpt1NW__sptrx1imidLine.p[1].x - self.xpt1[0], xpt1NW__sptrx1imidLine.p[1].y - self.xpt1[1]])
        tilt = np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) + np.pi/6
        xpt1W__sptrx2Lower = self.eq.draw_line(xpt1_dict['W'], {'line' : (sptrx2imidLine__LITP, tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        sptrx2__LowerPsiMaxInner = self.eq.draw_line(xpt1W__sptrx2Lower.p[-1], {'line' : (imidLine__LITP, tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

        # Computing angle between line tangent to sptrx1 in NE direction and unit vector < 1, 0 >
        v = np.array([xpt1NE__sptrx1omidLine.p[1].x - self.xpt1[0], xpt1NE__sptrx1omidLine.p[1].y - self.xpt1[1]])
        tilt = np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) - np.pi/6
        xpt1E__sptrx2Lower = self.eq.draw_line(xpt1_dict['E'], {'line': (sptrx2omidLine__LOTP, tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        sptrx2__LowerPsiMaxOuter = self.eq.draw_line(xpt1E__sptrx2Lower.p[-1], {'line' : (omidLine__LOTP, tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

        sptrx1__core_top = self.eq.draw_line(xpt2__sptrx1.p[-1], {'line' : (core, topLine_tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        sptrx1__core_inner = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : (core, inner_tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        sptrx1__core_outer = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : (core, outer_tilt)}, \
            option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        sptrx1__sptrx2_inner = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : (sptrx2_inner, inner_tilt)}, \
            option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        sptrx1__sptrx2_outer = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : (sptrx2_outer, outer_tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        sptrx2_inner__sptrx3_inner = self.eq.draw_line(xpt2NE__sptrx2imidLine.p[-1], {'line' : (sptrx3_inner, inner_tilt)}, \
            option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        sptrx2_outer__sptrx3_outer = self.eq.draw_line(xpt2NW__sptrx2omidLine.p[-1], {'line' : (sptrx3_outer, outer_tilt)}, \
            option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

        imidLine__LowerPsiMaxInner, LowerPsiMaxInner__LITP = imidLine__LITP.split(sptrx2__LowerPsiMaxInner.p[-1], add_split_point = True)
        sptrx2imidLine__sptrx2Lower, sptrx2Lower__LITP = sptrx2imidLine__LITP.split(xpt1W__sptrx2Lower.p[-1], add_split_point = True)

        # ============== Patch A1 ==============
        location = 'W'
        A1_N = LowerPsiMaxInner__LITP.reverse_copy()
        A1_S = sptrx2Lower__LITP
        A1_E = sptrx2__LowerPsiMaxInner.reverse_copy()
        # =====================================================================================
        # Trimming the target_plate to conform to the patch boundary.
        # -------------------------------------------------------------------------------------
        # Recall LITP has a clockwise orientation.
        #
        # The inner 'split' trims all Point objects BEFORE the point of intersection of LITP 
        # and A1_S. Call this new Line object Line_A.
        #
        # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
        # and A1_N. This new Line object is the plate facing boundary of the Patch.
        # =====================================================================================
        A1_W = (LITP.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
        A1 = DNL_Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'A1', platePatch = True, plateLocation = location)

        # ============== Patch B1 ==============
        location = 'W'
        B1_N = A1_S.reverse_copy()
        B1_E = xpt1W__sptrx2Lower.reverse_copy()
        B1_S = xpt1__LITP
        B1_W = (LITP.split(B1_S.p[-1])[1]).split(B1_N.p[0], add_split_point = True)[0]
        B1 = DNL_Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'B1', platePatch = True, plateLocation = location)

        # ============== Patch C1 ==============
        location = 'W'
        C1_N = B1_S.reverse_copy()
        C1_E = xpt1__psiMinPF1
        C1_S = psiMinPF1__LITP
        C1_W = (LITP.split(C1_S.p[-1])[1]).split(C1_N.p[0], add_split_point = True)[0]
        C1 = DNL_Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'C1', platePatch = True, plateLocation = location)

        # ============== Patch A2 ==============
        A2_N = imidLine__LowerPsiMaxInner.reverse_copy()
        A2_E = sptrx2_inner__sptrx3_inner.reverse_copy()
        A2_S = sptrx2imidLine__sptrx2Lower
        A2_W = A1_E.reverse_copy()
        A2 = DNL_Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'A2')

        # ============== Patch B2 ==============
        B2_N = A2_S.reverse_copy()
        B2_E = sptrx1__sptrx2_inner.reverse_copy()
        B2_S = xpt1NW__sptrx1imidLine.reverse_copy()
        B2_W = B1_E.reverse_copy()
        B2 = DNL_Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'B2')

        # ============== Patch C2 ==============
        C2_N = B2_S.reverse_copy()
        C2_E = sptrx1__core_inner
        C2_S = psiMinCore__imidLineCore.reverse_copy()
        C2_W = xpt1N__psiMinCore.reverse_copy()
        C2 = DNL_Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'C2')

        # ============== Patch A3 ==============
        A3_N = psiMaxInner__imidLine.reverse_copy()
        A3_E = xpt2__psiMaxInner.reverse_copy()
        A3_S = xpt2NE__sptrx2imidLine
        A3_W = A2_E.reverse_copy()
        A3 = DNL_Patch([A3_N, A3_E, A3_S, A3_W], patchName = 'A3')

        # ============== Patch B3 ==============
        B3_N = A3_S.reverse_copy()
        B3_E = xpt2__sptrx1
        B3_S = sptrx1imidLine__topLine.reverse_copy()
        B3_W = B2_E.reverse_copy()
        B3 = DNL_Patch([B3_N, B3_E, B3_S, B3_W], patchName = 'B3')

        # ============== Patch C3 ==============
        C3_N = B3_S.reverse_copy()
        C3_E = sptrx1__core_top
        C3_S = imidLineCore__topLine.reverse_copy()
        C3_W = C2_E.reverse_copy()
        C3 = DNL_Patch([C3_N, C3_E, C3_S, C3_W], patchName = 'C3')

        # ============== Patch A4 ==============
        location = 'E'
        A4_N = psiMaxInner__UITP
        A4_S = xpt2__UITP.reverse_copy()
        A4_E = (UITP.split(A4_N.p[-1])[1]).split(A4_S.p[0], add_split_point = True)[0]
        A4_W = A3_E.reverse()
        A4 = DNL_Patch([A4_N, A4_E, A4_S, A4_W], patchName = 'A4', plateLocation = location)

        # ============== Patch B4 ==============
        location = 'E'
        B4_N = A4_S.reverse_copy()
        B4_S = psiMinPF2_A__UITP.reverse_copy()
        B4_E = (UITP.split(B4_N.p[-1]))[1].split(B4_S.p[0], add_split_point = True)[0]
        B4_W = xpt2__psiMinPF2_A.reverse_copy()
        B4 = DNL_Patch([B4_N, B4_E, B4_S, B4_W], patchName = 'B4', plateLocation = location)
        # ============== Patch C4 ==============
        location = 'E'
        C4_N = B4_S.reverse_copy()
        C4_S = psiMinPF2_B__UITP.reverse_copy()
        C4_E = (UITP.split(C4_N.p[-1])[1]).split(C4_S.p[0], add_split_point = True)[0]
        C4_W = xpt2__psiMinPF2_B.reverse_copy()
        C4 = DNL_Patch([C4_N, C4_E, C4_S, C4_W], patchName = 'C4', plateLocation = location)
        # ============== Patch A5 ==============
        location = 'W'
        A5_N = psiMaxOuter__UOTP.reverse_copy()
        A5_S = xpt2__UOTP
        A5_E = xpt2__psiMaxOuter.reverse_copy()
        A5_W = (UOTP.split(A5_S.p[-1])[1]).split(A5_N.p[0], add_split_point = True)[0]
        A5 = DNL_Patch([A5_N, A5_E, A5_S, A5_W], patchName = 'A5', plateLocation = location)

        # ============== Patch B5 ==============
        location = 'W'
        B5_N = A5_S.reverse_copy()
        B5_S = psiMinPF2_A__UOTP
        B5_E = xpt2__psiMinPF2_A
        B5_W = (UOTP.split(B5_S.p[-1])[1]).split(B5_N.p[0], add_split_point = True)[0]
        B5 = DNL_Patch([B5_N, B5_E, B5_S, B5_W], patchName = 'B5', plateLocation = location)

        # ============== Patch C5 ==============
        location = 'W'
        C5_N = B5_S.reverse_copy()
        C5_S = psiMinPF2_B__UOTP
        C5_E = xpt2__psiMinPF2_B
        C5_W = (UOTP.split(C5_S.p[-1])[1]).split(C5_N.p[0], add_split_point = True)[0]
        C5 = DNL_Patch([C5_N, C5_E, C5_S, C5_W], patchName = 'C5', plateLocation = location)

        # ============== Patch A6 ==============
        A6_N = psiMaxOuter__omidLine
        A6_S = xpt2NW__sptrx2omidLine.reverse_copy()
        A6_E = sptrx2_outer__sptrx3_outer.reverse_copy()
        A6_W = A5_E.reverse_copy()
        A6 = DNL_Patch([A6_N, A6_E, A6_S, A6_W], patchName = 'A6')

        # ============== Patch B6 ==============
        B6_N = A6_S.reverse_copy()
        B6_S = sptrx1omidLine__topLine
        B6_E = sptrx1__sptrx2_outer.reverse_copy()
        B6_W = B3_E.reverse_copy()
        B6 = DNL_Patch([B6_N, B6_E, B6_S, B6_W], patchName = 'B6')

        # ============== Patch C6 ==============
        C6_N = B6_S.reverse_copy()
        C6_S = omidLineCore__topLine
        C6_E = sptrx1__core_outer
        C6_W = C3_E.reverse_copy()
        C6 = DNL_Patch([C6_N, C6_E, C6_S, C6_W], patchName = 'C6')

        # ============== Patch A7 ==============
        omidLine__LowerPsiMaxOuter, LowerPsiMaxOuter__LOTP = omidLine__LOTP.split(sptrx2__LowerPsiMaxOuter.p[-1], add_split_point = True)
        sptrx2imidLine__sptrx2Lower, sptrx2Lower__LOTP = sptrx2omidLine__LOTP.split(xpt1E__sptrx2Lower.p[-1], add_split_point = True)
        A7_N = omidLine__LowerPsiMaxOuter
        A7_S = sptrx2imidLine__sptrx2Lower.reverse_copy()
        A7_E = sptrx2__LowerPsiMaxOuter.reverse_copy()
        A7_W = A6_E.reverse_copy()
        A7 = DNL_Patch([A7_N, A7_E, A7_S, A7_W], patchName = 'A7')

        # ============== Patch B7 ==============
        B7_N = A7_S.reverse_copy()
        B7_S = xpt1NE__sptrx1omidLine
        B7_E = xpt1E__sptrx2Lower.reverse_copy()
        B7_W = B6_E.reverse_copy()
        B7 = DNL_Patch([B7_N, B7_E, B7_S, B7_W], patchName = 'B7')

        # ============== Patch C7 ==============
        C7_N = B7_S.reverse_copy()
        C7_S = psiMinCore__omidLineCore
        C7_E = C2_W.reverse_copy()
        C7_W = C6_E.reverse_copy()
        C7 = DNL_Patch([C7_N, C7_E, C7_S, C7_W], patchName = 'C7')

        # ============== Patch A8 ==============
        location = 'E'
        A8_N = LowerPsiMaxOuter__LOTP
        A8_S = sptrx2Lower__LOTP.reverse_copy()
        A8_E = (LOTP.split(A8_N.p[-1])[1]).split(A8_S.p[0], add_split_point = True)[0]
        A8_W = A7_E.reverse_copy()
        A8 = DNL_Patch([A8_N, A8_E, A8_S, A8_W], patchName = 'A8', plateLocation = location)

        # ============== Patch B8 ==============
        location = 'E'
        B8_N = A8_S.reverse_copy()
        B8_S = xpt1__LOTP.reverse_copy()
        B8_E = (LOTP.split(B8_N.p[-1])[1]).split(B8_S.p[0], add_split_point = True)[0]
        B8_W = B7_E.reverse_copy()
        B8 = DNL_Patch([B8_N, B8_E, B8_S, B8_W], patchName = 'B8', plateLocation = location)

        # ============== Patch C8 ==============
        location = 'E'
        C8_N = B8_S.reverse_copy()
        C8_S = psiMinPF1__LOTP.reverse_copy()
        C8_E = (LOTP.split(C8_N.p[-1])[1]).split(C8_S.p[0], add_split_point = True)[0]
        C8_W = C1_E.reverse_copy()
        C8 = DNL_Patch([C8_N, C8_E, C8_S, C8_W], patchName = 'C8', plateLocation = location)

        self.patches = [A1, B1, C1, A2, B2, C2, A3, B3, C3, A4, B4, C4, A5, B5, C5, A6, B6, C6, A7, B7, C7, A8, B8, C8]
        patch_lookup = {}
        for patch in self.patches:
            patch_lookup[patch.patchName] = patch
            patch.plot_border()
            patch.fill()
        self.patch_lookup = patch_lookup

        p = self.patch_lookup
        self.patch_matrix = [[[None],  [None],  [None],  [None],  [None], [None], [None], [None], [None], [None], [None], [None]], \
                            [[None], p['A1'], p['A2'], p['A3'], p['A4'], [None], [None], p['A5'], p['A6'], p['A7'], p['A8'], [None]],  \
                            [[None], p['B1'], p['B2'], p['B3'], p['B4'], [None], [None], p['B5'], p['B6'], p['B7'], p['B8'], [None]],  \
                            [[None], p['C1'], p['C2'], p['C3'], p['C4'], [None], [None], p['C5'], p['C6'], p['C7'], p['C8'], [None]], \
                            [[None],  [None],  [None],  [None],  [None], [None], [None], [None], [None], [None], [None], [None]]  \
                            ]

        self.catagorize_patches()

class SNL(Ingrid):
    """
    The SNL (Single-Null) class is the parent class for both upper-single null (USN)
    and lower-single null (LSN) configurations. 
    This base class handles the formatting and plotting of data obtained from an LSN or USN
    object.
    Parameter:
        - INGRID_object : Ingrid class object
        All SNL objects are children of the main Ingrid class. INGRID_object provides
        information such as YAML data, efit_psi, psi_norm, and plate_data.
    @author: garcia299
    """

    def __init__(self, INGRID_object):

        super().__init__(params = INGRID_object.yaml)
        self.efit_psi = INGRID_object.efit_psi
        self.psi_norm = INGRID_object.psi_norm
        self.eq = INGRID_object.eq

        self.plate_data = INGRID_object.plate_data

    def grid_diagram(self):
        """
        Create Grid matplotlib figure for an SNL object.
        @author: watkins35, garcia299
        """
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
          'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
          'seagreen', 'firebrick', 'saddlebrown']

        try:
            plt.close('INGRID: Grid')
        except:
            pass

        plt.figure('INGRID: Grid', figsize=(6,10))
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
        """ 
        Generates the patch diagram for a given configuration. 
        @author: watkins35, garcia299
        """
        
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']

        try:
            plt.close('INGRID: Patch Map')
        except:
            pass

        plt.figure('INGRID: Patch Map', figsize=(6, 10))
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
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the 
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches:
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1



        # Total number of poloidal indices in subgrid.
        np_total = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:-1]])) + 2
        nr_total = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:3]])) + 2

        rm = np.zeros((np_total, nr_total, 5), order = 'F')
        zm = np.zeros((np_total, nr_total, 5), order = 'F')

        ixcell = 0
        jycell = 0

        # Iterate over all the patches in our SNL configuration (we exclude guard cells denoted by '[None]')
        for ixp in range(1, 7):
            
            nr_sum = 0
            for jyp in range(1, 3):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]
                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch. 
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                        print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm)):
            rm[i] = rm[i][::-1]
        for i in range(len(zm)):
            zm[i] = zm[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb, config)
        self.zm = self.add_guardc(zm, ixlb, ixrb, config)

        try:
            debug = self.yaml['DEBUG']['visual']['gridue']
        except:
            debug = False
        
        if debug:
            self.animate_grid()


    def add_guardc(self, cell_map, ixlb, ixrb, config, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, config, boundary):
            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # TODO: Edit the code to reflect this at some point so the next reader is not overwhelmed.
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


    def animate_grid(self):

        try:
            plt.close('INGRID: Debug')
        except:
            pass
        plt.figure('INGRID: Debug', figsize=(6, 10))
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('visualize gridue')

        k = [1,2,4,3,1]

        for i in range(len(self.rm)):
            for j in range(len(self.rm[0])):
                plt.plot(self.rm[i][j][k], self.zm[i][j][k])
                plt.pause(0.01)

class LSN(SNL, Ingrid):

    def __init__(self, INGRID_object):
        super().__init__(INGRID_object)
        self.config = 'LSN'
        print('=' * 80)
        print('LSN Object!')
        print('=' * 80 + '\n')

    def construct_grid(self, np_cells = 1, nr_cells = 1):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.

        try:
            visual = self.yaml['DEBUG']['visual']['subgrid']
        except:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['grid_generation']
        except:
            verbose = False

        try:
            np_sol = self.yaml['grid_params']['grid_generation']['np_sol']
        except:
            np_sol = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_core = self.yaml['grid_params']['grid_generation']['np_core']
        except:
            np_core = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_pf = self.yaml['grid_params']['grid_generation']['np_pf']
        except:
            np_pf = self.yaml['grid_params']['grid_generation']['np_global']

        if np_sol != np_core:
            print('WARNING: SOL and CORE must have equal POLOIDAL np values!\nSetting np values' \
                + ' to the minimum of np_sol and np_core.\n')
            np_sol = np_core = np.amin([np_sol, np_core])

        try:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_sol']
        except:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_global']

        try:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_core']
        except:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_global']
        try:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_pf']
        except:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_global']

        if nr_pf != nr_core:
            print('WARNING: PF and CORE must have equal RADIAL nr values!\nSetting nr values' \
                + ' to the minimum of nr_pf and nr_core.\n')
            nr_pf = nr_core = np.amin([nr_pf, nr_core])

        primary_xpt = Point([self.yaml['grid_params']['rxpt'], self.yaml['grid_params']['zxpt']])
        for patch in self.patches:
            if patch in self.SOL:
                nr_cells = nr_sol
                np_cells = np_sol
                print('Patch "{}" is in SOL'.format(patch.patchName))
            elif patch in self.CORE:
                nr_cells = nr_core
                np_cells = np_core
                print('Patch "{}" is in CORE'.format(patch.patchName))
            elif patch in self.PF:
                nr_cells = nr_pf
                np_cells = np_pf
                print('Patch "{}" is in PF'.format(patch.patchName))

            patch.make_subgrid(self, np_cells, nr_cells, verbose = verbose, visual = visual)

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
            visual = self.yaml['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False
        try:
            inner_tilt = self.yaml['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0
        try:
            outer_tilt = self.yaml['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0

        self.itp = self.plate_data['plate_W1']['coordinates']
        self.otp = self.plate_data['plate_E1']['coordinates']

        xpt = self.eq.eq_psi
        magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
            self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

        psi_max = self.yaml['grid_params']['psi_max']
        psi_min_core = self.yaml['grid_params']['psi_min_core']
        psi_min_pf = self.yaml['grid_params']['psi_min_pf']

        ITP = Line(self.order_plate_points(Line([Point(i) for i in self.itp]), location = 'LITP', orientation = 'cw'))
        OTP = Line(self.order_plate_points(Line([Point(i) for i in self.otp]), location = 'LOTP', orientation = 'cw'))

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        inner_midLine = Line([LHS_Point, RHS_Point])
        inner_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        outer_midLine = Line([LHS_Point, RHS_Point])
        outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Drawing Separatrix
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        
        # Drawing Lower-SNL region
        if self.yaml['grid_params']['patch_generation']['use_NW']:
            tilt = self.eq.eq_psi_theta['NW'] + self.yaml['grid_params']['patch_generation']['NW_adjust']
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        if self.yaml['grid_params']['patch_generation']['use_NE']:
            tilt = self.eq.eq_psi_theta['NE'] + self.yaml['grid_params']['patch_generation']['NE_adjust']
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        
        xpt_ITP = self.eq.draw_line(xpt['SW'], {'line' : ITP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)        
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xpt_OTP = self.eq.draw_line(xpt['SE'], {'line' : OTP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        
        iPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : ITP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : ITP},option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        oPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : OTP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : OTP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        imidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        omidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        
        # Integrating horizontally along mid-line towards psiMax and psiMinCore
        imidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_max, inner_tilt)}, option = 'z_const', \
                direction = 'ccw', show_plot = visual, text = verbose)
        imidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_min_core, inner_tilt)}, option = 'z_const', \
                direction = 'cw', show_plot = visual, text = verbose)
        omidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_max, outer_tilt)}, option = 'z_const', \
                direction = 'cw', show_plot = visual, text = verbose)
        omidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_min_core, outer_tilt)}, option = 'z_const', \
                direction = 'ccw', show_plot = visual, text = verbose)
        
        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'cw', show_plot = visual, text = verbose)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'ccw', show_plot = visual, text = verbose)

        # IDL Patch
        location = 'W'
        IDL_N = iPsiMax_TP.reverse_copy()
        IDL_S = xpt_ITP
        IDL_E = xptW_psiMax.reverse_copy()
        # =====================================================================================
        # Trimming the target_plate to conform to the patch boundary.
        # -------------------------------------------------------------------------------------
        # Recall ITP has a clockwise orientation.
        #
        # The inner 'split' trims all Point objects BEFORE the point of intersection of ITP 
        # and IDL_S. Call this new Line object Line_A.
        #
        # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
        # and IDL_N. This new Line object is the plate facing boundary of the Patch.
        # =====================================================================================
        IDL_W = (ITP.split(IDL_S.p[-1])[1]).split(IDL_N.p[0], add_split_point = True)[0]
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        location = 'W'
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP
        IPF_E = xptS_psiMinPF
        IPF_W = (ITP.split(IPF_S.p[-1])[1]).split(IPF_N.p[0], add_split_point = True)[0]
        IPF = SNL_Patch([IPF_N, IPF_E, IPF_S, IPF_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # ISB Patch
        ISB_N = self.eq.draw_line(IDL_N.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        ISB_S = xptNW_midLine.reverse_copy()
        ISB_E = Line([ISB_N.p[-1], ISB_S.p[0]])
        ISB_W = xptW_psiMax
        ISB = SNL_Patch([ISB_N, ISB_E, ISB_S, ISB_W], patchName = 'ISB')

        # ICB Patch
        ICB_N = ISB_S.reverse_copy()
        ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        ICB_E = Line([ICB_N.p[-1], ICB_S.p[0]])
        ICB_W = xptN_psiMinCore.reverse_copy()
        ICB = SNL_Patch([ICB_N, ICB_E, ICB_S, ICB_W], patchName = 'ICB')

        # IST Patch
        IST_N = self.eq.draw_line(ISB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        IST_S = imidLine_topLine.reverse_copy()
        IST_E = Line([IST_N.p[-1], IST_S.p[0]])
        IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        IST = SNL_Patch([IST_N, IST_E, IST_S, IST_W], patchName = 'IST')

        # ICT Patch
        ICT_N = IST_S.reverse_copy()
        ICT_S = self.eq.draw_line(ICB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
        ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        ICT = SNL_Patch([ICT_N, ICT_E, ICT_S, ICT_W], patchName = 'ICT')

        # ODL Patch 
        location = 'E'
        ODL_N = oPsiMax_TP
        ODL_S = xpt_OTP.reverse_copy()
        ODL_E = (OTP.split(ODL_N.p[-1])[1]).split(ODL_S.p[0], add_split_point = True)[0]
        ODL_W = xptE_psiMax
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        location = 'E'
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP.reverse_copy()
        OPF_E = (OTP.split(OPF_N.p[-1])[1]).split(OPF_S.p[0], add_split_point = True)[0]
        OPF_W = xptS_psiMinPF.reverse_copy()
        OPF = SNL_Patch([OPF_N, OPF_E, OPF_S, OPF_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # OSB Patch 
        OSB_N = self.eq.draw_line(ODL_N.p[0], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        OSB_S = xptNE_midLine
        OSB_E = xptE_psiMax.reverse_copy()
        OSB_W = Line([OSB_S.p[-1], OSB_N.p[0]])
        OSB = SNL_Patch([OSB_N, OSB_E, OSB_S, OSB_W], patchName = 'OSB')

        # OCB Patch
        OCB_N = OSB_S.reverse_copy()
        OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        OCB_E = xptN_psiMinCore
        OCB_W = Line([OCB_S.p[-1], OCB_N.p[0]])
        OCB = SNL_Patch([OCB_N, OCB_E, OCB_S, OCB_W], patchName = 'OCB')

        # OST Patch
        OST_N = self.eq.draw_line(OSB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        OST_S = omidLine_topLine
        OST_E = Line([OST_N.p[-1], OST_S.p[0]])
        OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        OST = SNL_Patch([OST_N, OST_E, OST_S, OST_W], patchName = 'OST')

        # OCT Patch
        OCT_N = OST_S.reverse_copy()
        OCT_S = self.eq.draw_line(OCB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
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

        self.catagorize_patches()


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
        rb_prod = self.efit_psi.rcenter * self.efit_psi.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.efit_psi.get_psi(_r, _z)
                    _br = self.efit_psi.get_psi(_r, _z, tag = 'vz') / _r
                    _bz = -self.efit_psi.get_psi(_r, _z, tag = 'vr') / _r
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

class USN(SNL, Ingrid):
    def __init__(self, INGRID_object):
        self.config = 'USN'
        print('=' * 80)
        print('USN Object!')
        print('=' * 80 + '\n')
        super().__init__(INGRID_object)


    def construct_grid(self, np_cells = 1, nr_cells = 1):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.

        try:
            visual = self.yaml['DEBUG']['visual']['subgrid']
        except:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['grid_generation']
        except:
            verbose = False

        try:
            np_sol = self.yaml['grid_params']['grid_generation']['np_sol']
        except:
            np_sol = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_core = self.yaml['grid_params']['grid_generation']['np_core']
        except:
            np_core = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_pf = self.yaml['grid_params']['grid_generation']['np_pf']
        except:
            np_pf = self.yaml['grid_params']['grid_generation']['np_global']

        if np_sol != np_core:
            print('WARNING: SOL and CORE must have equal POLOIDAL np values!\nSetting np values' \
                + ' to the minimum of np_sol and np_core.\n')
            np_sol = np_core = np.amin([np_sol, np_core])

        try:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_sol']
        except:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_global']

        try:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_core']
        except:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_global']
        try:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_pf']
        except:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_global']

        if nr_pf != nr_core:
            print('WARNING: PF and CORE must have equal RADIAL nr values!\nSetting nr values' \
                + ' to the minimum of nr_pf and nr_core.\n')
            nr_pf = nr_core = np.amin([nr_pf, nr_core])

        primary_xpt = Point([self.yaml['grid_params']['rxpt'], self.yaml['grid_params']['zxpt']])
        for patch in self.patches:
            if patch in self.SOL:
                nr_cells = nr_sol
                np_cells = np_sol
                print('Patch "{}" is in SOL'.format(patch.patchName))
            elif patch in self.CORE:
                nr_cells = nr_core
                np_cells = np_core
                print('Patch "{}" is in CORE'.format(patch.patchName))
            elif patch in self.PF:
                nr_cells = nr_pf
                np_cells = np_pf
                print('Patch "{}" is in PF'.format(patch.patchName))

            patch.make_subgrid(self, np_cells, nr_cells, verbose = verbose, visual = visual)

            if patch.patchName == 'ODL':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'OPF':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'OSB':
                patch.adjust_corner(primary_xpt, 'SW')
            elif patch.patchName == 'OCB':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'ICB':
                patch.adjust_corner(primary_xpt, 'NE')
            elif patch.patchName == 'ISB':
                patch.adjust_corner(primary_xpt, 'SE')
            elif patch.patchName == 'IPF':
                patch.adjust_corner(primary_xpt, 'NW')
            elif patch.patchName == 'IDL':
                patch.adjust_corner(primary_xpt, 'SW')

        self.concat_grid(self.get_configuration())
        self.set_gridue()

    def construct_patches(self):

        try:
            visual = self.yaml['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False
        try:
            inner_tilt = self.yaml['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0
        try:
            outer_tilt = self.yaml['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0

        self.itp = self.plate_data['plate_E1']['coordinates']
        self.otp = self.plate_data['plate_W1']['coordinates']

        xpt = self.eq.eq_psi
        magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
            self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

        psi_max = self.yaml['grid_params']['psi_max']
        psi_min_core = self.yaml['grid_params']['psi_min_core']
        psi_min_pf = self.yaml['grid_params']['psi_min_pf']

        ITP = Line(self.order_plate_points(Line([Point(i) for i in self.itp]), location = 'UITP', orientation = 'cw'))
        OTP = Line(self.order_plate_points(Line([Point(i) for i in self.otp]), location = 'UOTP', orientation = 'cw'))

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        inner_midLine = Line([LHS_Point, RHS_Point])
        inner_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        outer_midLine = Line([LHS_Point, RHS_Point])
        outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])

        # Drawing Separatrix
        xptNW_midLine = self.eq.draw_line(xpt['NW'], {'line' : outer_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        xptN_psiMinCore = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xptNE_midLine = self.eq.draw_line(xpt['NE'], {'line' : inner_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        
        # Drawing Lower-SNL region
        if self.yaml['grid_params']['patch_generation']['use_NW']:
            tilt = self.eq.eq_psi_theta['NW'] + self.yaml['grid_params']['patch_generation']['NW_adjust']
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        if self.yaml['grid_params']['patch_generation']['use_NE']:
            tilt = self.eq.eq_psi_theta['NE'] + self.yaml['grid_params']['patch_generation']['NE_adjust']
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        
        xpt_OTP = self.eq.draw_line(xpt['SW'], {'line' : OTP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xpt_ITP = self.eq.draw_line(xpt['SE'], {'line' : ITP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        
        oPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : OTP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : OTP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        iPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : ITP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : ITP},option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        

        imidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        omidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        # Integrating horizontally along mid-line towards psiMax and psiMinCore
        omidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_max, outer_tilt)}, option = 'z_const', \
                direction = 'cw', show_plot = visual, text = verbose, color = 'red')
        omidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_min_core, outer_tilt)}, option = 'z_const', \
                direction = 'ccw', show_plot = visual, text = verbose, color = 'blue')
        imidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_max, inner_tilt)}, option = 'z_const', \
                direction = 'ccw', show_plot = visual, text = verbose, color = 'green')
        imidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_min_core, inner_tilt)}, option = 'z_const', \
                direction = 'cw', show_plot = visual, text = verbose, color = 'black')
        
        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'ccw', show_plot = visual, text = verbose)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'cw', show_plot = visual, text = verbose)

        # IDL Patch
        location = 'E'   
        IDL_N = iPsiMax_TP
        IDL_S = xpt_ITP.reverse_copy()  
        # =====================================================================================
        # The inner 'split' trims all Point objects BEFORE the point of intersection of ITP 
        # and IDL_N. Call this new Line object Line_A.
        #
        # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
        # and IDL_S. This new Line object is the plate facing boundary of the Patch.
        # =====================================================================================
        IDL_E = (ITP.split(IDL_N.p[-1])[1]).split(IDL_S.p[0], add_split_point = True)[0]
        IDL_W = xptE_psiMax   
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        location = 'E'
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP.reverse_copy()
        IPF_E = (ITP.split(IPF_N.p[-1])[1]).split(IPF_S.p[0], add_split_point = True)[0]
        IPF_W = xptS_psiMinPF.reverse_copy()
        IPF = SNL_Patch([IPF_N, IPF_E, IPF_S, IPF_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # ISB Patch
        ISB_N = self.eq.draw_line(IDL_N.p[0], {'line' : inner_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        ISB_S = xptNE_midLine
        ISB_E = xptE_psiMax.reverse_copy()
        ISB_W = Line([ISB_S.p[-1], ISB_N.p[0]])
        ISB = SNL_Patch([ISB_N, ISB_E, ISB_S, ISB_W], patchName = 'ISB')

        # ICB Patch
        ICB_N = ISB_S.reverse_copy()
        ICB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        ICB_E = xptN_psiMinCore
        ICB_W = Line([ICB_S.p[-1], ICB_N.p[0]])
        ICB = SNL_Patch([ICB_N, ICB_E, ICB_S, ICB_W], patchName = 'ICB')
        
        # IST Patch
        IST_N = self.eq.draw_line(ISB_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        IST_S = imidLine_topLine
        IST_E = Line([IST_N.p[-1], IST_S.p[0]])
        IST_W = Line([IST_S.p[-1], IST_N.p[0]])
        IST = SNL_Patch([IST_N, IST_E, IST_S, IST_W], patchName = 'IST')

        # ICT Patch
        ICT_N = IST_S.reverse_copy()
        ICT_S = self.eq.draw_line(ICB_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        ICT_E = Line([ICT_N.p[-1], ICT_S.p[0]])
        ICT_W = Line([ICT_S.p[-1], ICT_N.p[0]])
        ICT = SNL_Patch([ICT_N, ICT_E, ICT_S, ICT_W], patchName = 'ICT')
        
        # ODL Patch
        location = 'W'
        ODL_N = oPsiMax_TP.reverse_copy()
        ODL_S = xpt_OTP
        ODL_E = xptW_psiMax.reverse_copy()
        ODL_W = (OTP.split(ODL_S.p[-1])[1]).split(ODL_N.p[0], add_split_point = True)[0]
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        location = 'W'
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP
        OPF_E = xptS_psiMinPF
        OPF_W = (OTP.split(OPF_S.p[-1])[1]).split(OPF_N.p[0], add_split_point = True)[0]
        OPF = SNL_Patch([OPF_N, OPF_E, OPF_S, OPF_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # OSB Patch
        OSB_N = self.eq.draw_line(ODL_N.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        OSB_S = xptNW_midLine.reverse_copy()
        OSB_E = Line([OSB_N.p[-1], OSB_S.p[0]])
        OSB_W = xptW_psiMax
        OSB = SNL_Patch([OSB_N, OSB_E, OSB_S, OSB_W], patchName = 'OSB')

        # OCB Patch
        OCB_N = OSB_S.reverse_copy()
        OCB_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        OCB_E = Line([OCB_N.p[-1], OCB_S.p[0]])
        OCB_W = xptN_psiMinCore.reverse_copy()
        OCB = SNL_Patch([OCB_N, OCB_E, OCB_S, OCB_W], patchName = 'OCB')

        # OST Patch
        OST_N = self.eq.draw_line(OSB_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        OST_S = omidLine_topLine.reverse_copy()
        OST_E = Line([OST_N.p[-1], OST_S.p[0]])
        OST_W = Line([OST_S.p[-1], OST_N.p[0]])
        OST = SNL_Patch([OST_N, OST_E, OST_S, OST_W], patchName = 'OST')

        # OCT Patch
        OCT_N = OST_S.reverse_copy()
        OCT_S = self.eq.draw_line(OCB_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        OCT_E = Line([OCT_N.p[-1], OCT_S.p[0]])
        OCT_W = Line([OCT_S.p[-1], OCT_N.p[0]])
        OCT = SNL_Patch([OCT_N, OCT_E, OCT_S, OCT_W], patchName = 'OCT')

        self.patches = [ODL, OPF, OSB, OCB, OST, OCT, IST, ICT, ISB, ICB, IDL, IPF]
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
        self.catagorize_patches()

    def set_gridue(self):
        """
        Prepare the relevant arrays for writing to GRIDUE.
        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patch_lookup['ODL'].npol - 1
        ixpt2 = ixrb - self.patch_lookup['IDL'].npol + 1
        iyseparatrix1 = self.patch_lookup['ODL'].nrad - 1
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
        rb_prod = self.efit_psi.rcenter * self.efit_psi.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = self.efit_psi.get_psi(_r, _z)
                    _br = self.efit_psi.get_psi(_r, _z, tag = 'vz') / _r
                    _bz = -self.efit_psi.get_psi(_r, _z, tag = 'vr') / _r
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
    IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
    IngridWindow.mainloop()
