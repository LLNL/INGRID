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
import GUI.IngridApp as IA

from OMFITgeqdsk import OMFITgeqdsk
from Interpol.Setup_Grid_Data import Efit_Data
from line_tracing import LineTracing
from Root_Finder import RootFinder
from geometry import Point, Line, SNL_Patch, segment_intersect


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
            'config' : '', 'num_xpt' : 1, \
            'psi_max' : 0.0, 'psi_max_r' : 0.0, 'psi_max_z' : 0.0, \
            'psi_min_core' : 0.0, 'psi_min_core_r' : 0.0, 'psi_min_core_z' : 0.0, \
            'psi_min_pf' : 0.0, 'psi_min_pf_r' : 0.0, 'psi_min_pf_z' : 0.0, \
            'rmagx' : 0.0, 'zmagx' : 0.0, \
            'grid_generation' : { \
                'np_global' : 3, 'nr_global' : 2, \
                'nr_sol' : 2, 'nr_core' : 2, 'nr_pf' : 2, \
                'radial_f_sol' : 'x, x', 'radial_f_core' : 'x, x', \
                'radial_f_pf' : 'x, x', \
            }, \
            'patch_generation' : { \
                'rmagx_shift' : 0.0, 'zmagx_shift' : 0.0, \
                'inner_tilt' : 0.0, 'outer_tilt' : 0.0, \
                'rxpt' : 0.0, 'zxpt' : 0.0, \
                'use_NW' : False, 'use_NE' : False, \
                'use_SW' : False, 'use_SE' : False,
                'NW_adjust' : 0.0, 'NE_adjust' : 0.0, \
            } \
        }

        self.default_integrator_params = { \
            'dt' : 0.01, 'eps' : 5e-5, \
            'first_step' : 1e-5, 'step_ratio' : 0.02, \
            'tol' : 5e-3 \
        }

        self.default_target_plates_params = { \
            'plate_E1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0}, \
            'plate_W1' : {'file' : '', 'name' : '', 'np_local' : 3, 'poloidal_f' : 'x, x', 'zshift' : 0.0}
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
                continue

            # Second level entries within YAML file dump.
            for sub_item in self.default_values_lookup[item].keys():
                try:
                    params[item][sub_item]
                except KeyError:
                    print('Could not find "{}/{}" in YAML file.'.format(item, sub_item))
                    params[item][sub_item] = get_default_values(item, sub_item)
                    continue
                if item in ['grid_params', 'target_plates'] and sub_item in ['patch_generation', 'grid_generation', 'plate_E1', 'plate_W1']:
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
            if end.x - start.x > 0:
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
        self.efit_psi.plot_data()
    
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
                                  name='psi norm')
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
        self.psi_norm.plot_data()
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
        self.SOL = m[1][1:-1]
        self.CORE = m[2][2:-2]
        self.PF = [m[2][1], m[2][-2]]

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
            print('Double null configurations not yet supported...')
        
        self.yaml['grid_params']['config'] = self.eq.config
        return self.yaml['grid_params']['config']

    def _analyze_topology(self):
        """
        Analyze the topological structure of the domain in order to determine the configuration.
        Once analysis is complete, the relevant SNL or DNL Object will be instantiated.

        @author: garcia299
        """
        config = self._classify_gridtype()
        if config == 'LSN':
            ingrid_topology = LSN(self)
        elif config == 'USN':
            ingrid_topology = USN(self)

        self.current_topology = ingrid_topology

    def _get_configuration(self):
        return self.eq.config

    def export(self, fname = 'gridue'):
        """ Saves the grid as an ascii file """
        self.write_gridue(self.current_topology.gridue_params, fname)

    def write_gridue(self, gridue_params, fname = 'gridue'):
        """
        Write the INGRID data to a gridue formatted file.
        """
        
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
        xpt_OTP = self.eq.draw_line(xpt['SE'], {'line' : self.otp}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        
        iPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : ITP}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : ITP},option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        oPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : OTP}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.otp}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

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


    def set_gridue_manual(self,rm,zm,nx,ny,ixpt1,ixpt2,iysptrx):
        """
        Prepare the relevant arrays for writing to GRIDUE.
        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(rm) - 2
        nxm = len(rm) - 2
        nym = len(rm[0]) - 2

        psi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        br = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bz = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bpol = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        bphi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
        b = np.zeros((nxm + 2, nym + 2, 5), order = 'F')

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

        self.gridue_params = {'nxm' : nx, 'nym' : ny, 'ixpt1' : ixpt1, 'ixpt2' : ixpt2, 'iyseptrx1' : iysptrx, \
            'rm' : rm, 'zm' : zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b}


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
        
        xpt_OTP = self.eq.draw_line(xpt['SW'], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xpt_ITP = self.eq.draw_line(xpt['SE'], {'line' : self.itp}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        
        oPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_OTP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.otp}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        iPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : self.itp}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_ITP = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : self.itp},option = 'theta', direction = 'cw', show_plot = visual, text = verbose)        

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

        import pdb
        pdb.set_trace()

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
    IngridWindow.geometry("830x490")
    IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
    IngridWindow.mainloop()
