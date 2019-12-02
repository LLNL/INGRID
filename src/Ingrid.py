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
from geometry import Point, Line, SNL_Patch, segment_intersect


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

        # Process params as a YAML file. Ensure ends up as a raw dictionary.s
        self.process_yaml(params)
        self.get_yaml()

    def process_yaml(self, params):

        for item in self.default_values_lookup.keys():

            try:
                params[item]
            except KeyError:
                print('Could not find "{}" in YAML file.'.format(item))
                params[item] = self.get_default_values(item)
                continue

            for sub_item in self.default_values_lookup[item].keys():
                try:
                    params[item][sub_item]
                except KeyError:
                    print('Could not find "{}/{}" in YAML file.'.format(item, sub_item))
                    params[item][sub_item] = self.get_default_values(item, sub_item)
                    continue
                if item in ['grid_params', 'target_plates'] and sub_item in ['patch_generation', 'grid_generation', 'plate_E1', 'plate_W1']:
                    for plate_attribute in self.default_values_lookup[item][sub_item].keys():
                        try:
                            params[item][sub_item][plate_attribute]
                        except:
                            params[item][sub_item][plate_attribute] = self.get_default_values(item, sub_item, plate_attribute)

        # Store YAML data as dictionary object
        self.yaml = params
        self.grid_params = params['grid_params']
        self.integrator_params = params['integrator_params']
        self.target_plates = params['target_plates']
        self.DEBUG = params['DEBUG']

    def get_yaml(self):
        for item in self.yaml.keys():
            print('=' * 80)
            print('INGRID {}: {}'.format(item, self.yaml[item]))
            print('=' * 80 + '\n')

    def get_default_values(self, item, sub_item = None, attribute = None):
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

        self.OMFIT_psi = g

    def read_target_plates(self):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        """

        try:
            debug = self.yaml['DEBUG']['verbose']['target_plates']
        except:
            debug = False

        target_plates = self.yaml['target_plates']
        plate_data = {}

        for plate in target_plates:

            zshift = target_plates[plate]['zshift']

            plate_data[plate] = {'coordinates' : [], 'psi_min_rz' : (), 'psi_max_rz' : ()}

            if not 'name' in target_plates[plate].keys():
                target_plates[plate].update({'name' : plate})
            
            elif target_plates[plate]['name'] == '':
                target_plates[plate]['name'] = plate

            try:
                open(target_plates[plate]['file'])
            except FileNotFoundError:
                continue


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

            for i in range(len(plate_data[plate]['coordinates'])):
                curr_v = self.efit_psi.get_psi(plate_data[plate]['coordinates'][i][0], plate_data[plate]['coordinates'][i][1])
                min_v_test = self.efit_psi.get_psi(temp_min[0], temp_min[1])
                max_v_test = self.efit_psi.get_psi(temp_max[0], temp_max[1])
                temp_min = plate_data[plate]['coordinates'][i] if curr_v <= min_v_test else temp_min
                temp_max = plate_data[plate]['coordinates'][i] if curr_v >= max_v_test else temp_max
            
            plate_data[plate]['psi_min_rz'] = temp_min
            plate_data[plate]['psi_max_rz'] = temp_max

            if debug:
                print('Using target plate "{}": {}'.format(target_plates[plate]['name'], plate_data[plate]))

         
        self.plate_data = plate_data

    def plot_target_plates(self):
        """ Plots all target plates on the current figure """
        
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

    def plot_psinorm(self):
        self.psi_norm.plot_data()
        self.plot_target_plates()

    def find_psi_lines(self, tk_controller = None):
        self.psi_finder = RootFinder(self.psi_norm, mode = 'psi_finder', controller = tk_controller)

    def init_LineTracing(self, refresh = False):
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
        if refresh:
            self.eq = LineTracing(self.psi_norm, self.yaml)
        self.eq.disconnect()

    def refresh_LineTracing(self):
        self.init_LineTracing(refresh = True)

    def catagorize_patches(self):
        m = self.patch_matrix
        self.SOL = m[1][1:-1]
        self.CORE = m[2][2:-2]
        self.PF = [m[2][1], m[2][-2]]

    def _classify_gridtype(self):
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

    def __init__(self, INGRID_object):
        super().__init__(params = INGRID_object.yaml)
        self.efit_psi = INGRID_object.efit_psi
        self.psi_norm = INGRID_object.psi_norm
        self.eq = INGRID_object.eq

        self.plate_data = INGRID_object.plate_data

    def grid_diagram(self):
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
        """ Generates the patch diagram for a given configuration. """
        
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

        def order_plate_points(plate):
            """
            Sets the points in the target plate to have an orientation
            increasing in R and Z. Checks the endpoints to satisfy this criteria.
            """
            left = plate.p[0]
            right = plate.p[-1]
            # Endpoints on same vertical line.
            if left.x == right.x:
                # If rhs endpoint above lhs endpoint.
                if right.y > left.y:
                    # Return target plate as is
                    return plate.copy().p[::]
                else:
                    # Else flip plate orientation.
                    return plate.copy().p[::-1]
            # Endpoints on same horizontal line.
            elif left.y == right.y:
                # If lhs endpoint to the left of rhs endpoint.
                if left.x < right.x:
                    # Return target plate as is
                    return plate.copy().p
                else:
                    # Else flip plate orientation
                    return plate.copy().p[::-1]
            # Endpoints are on sloped line.
            # Check if lhs endpoint is on the left of rhs endpoint
            elif left.x < right.x:
                return plate.copy().p
            else:
                return plate.copy().p[::-1]

        def get_insert_index(arr, val):

            # Iterate over all points in a given target plate
            
            # Check the current point and the following point
            # since we are looking for where the given value will
            # need to be inserted in an already populated list.

            # 
            index_found = True
            for i in range(len(arr) - 1):
                p1, p2 = arr[i], arr[i+1]
                if (p1.x, p1.y) == (p2.x, p2.y):
                    continue
                # vertical segment
                if p1.x == p2.x:
                    # value in-between vertically increasing segment
                    if (val.y >= p1.y and val.y <= p2.y):
                        break
                # horizontal segment
                elif p1.y == p2.y:
                    # value in-between horizontally increasing segment
                    if (val.x >= p1.x and val.x <= p2.x):
                        break
                # sloped segment
                elif p1.x < p2.x:
                    if (val.x >= p1.x and val.x <= p2.x) \
                        and (val.y <= p1.y and val.y >= p2.y):
                        break
                    if (val.x >= p1.x and val.x <= p2.x) \
                        and (val.y >= p1.y and val.y <= p2.y):
                        break
                elif p1.x > p2.x:
                    if (val.x >= p2.x and val.x <= p1.x)\
                        and (val.y <= p2.y and val.y >= p1.y):
                        break
                elif i == len(arr.p) - 2:
                    index_found = False
            
            return i+1 if index_found else null

        def set_face(plate, south_point, north_point, location):
            """
            Trim a Line object adjacent to a divertor plate.
            The Line object will be defined by the interval of points
            [min_point, max_point] that form a subset of a target plate.
            """

            new_leg = order_plate_points(plate)
            i = get_insert_index(new_leg, south_point)
            j = get_insert_index(new_leg, north_point)

            if location == 'W':
                new_leg.insert(i, south_point)
                new_leg.insert(j, north_point)
            elif location == 'E':
                new_leg.insert(j, north_point)
                new_leg.insert(i, south_point)

            lookup = {}
            for i in range(len(new_leg)):
                lookup[new_leg[i]] = i

            start = lookup[north_point]
            end = lookup[south_point]

            if start > end:
                start, end = end, start

            # Return a reverse copy to maintain clockwise orientation within the LSN SNL patch.
            return Line([p for p in new_leg[start:end+1]]).reverse_copy()

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

        ITP = Line(order_plate_points(Line([Point(i) for i in self.itp])))
        OTP = Line(order_plate_points(Line([Point(i) for i in self.otp])))

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
        IDL_N = iPsiMax_TP.reverse_copy()
        IDL_S = xpt_ITP
        IDL_E = xptW_psiMax.reverse_copy()
        location = 'W'
        IDL_W = set_face(ITP, IDL_S.p[-1], IDL_N.p[0], location = location)
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP
        IPF_E = xptS_psiMinPF
        location = 'W'
        IPF_W = set_face(ITP, IPF_S.p[-1], IPF_N.p[0], location = location)
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
        ODL_N = oPsiMax_TP
        ODL_S = xpt_OTP.reverse_copy()
        location = 'E'
        ODL_E = set_face(OTP, ODL_S.p[0], ODL_N.p[-1], location = location)
        ODL_W = xptE_psiMax
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP.reverse_copy()
        location = 'E'
        OPF_E = set_face(OTP, OPF_S.p[0], OPF_N.p[-1], location = location)
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

        def order_plate_points(plate):
            """
            Sets the points in the target plate to have an orientation
            increasing in R and Z. Checks the endpoints to satisfy this criteria.
            """
            start = plate.p[0]
            end = plate.p[-1]
            # Endpoints on same vertical line.
            if start.x == end.x:
                # If rhs endpoint above lhs endpoint.
                if end.y < start.y:
                    # Return target plate as is
                    return plate.copy().p
                else:
                    # Else flip plate orientation.
                    return plate.copy().p[::-1]
            # Endpoints on same horizontal line.
            elif start.y == start.y:
                # If lhs endpoint to the left of rhs endpoint.
                if end.x > start.x:
                    # Return target plate as is
                    return plate.copy().p
                else:
                    # Else flip plate orientation
                    return plate.copy().p[::-1]
            # Endpoints are on sloped line.
            # Check if lhs endpoint is on the left of rhs endpoint
            elif end.x < start.x:
                return plate.copy().p
            else:
                return plate.copy().p[::-1]

        def get_insert_index(arr, val):

            # Iterate over all points in a given target plate
            
            # Check the current point and the following point
            # since we are looking for where the given value will
            # need to be inserted in an already populated list.

            index_found = True
            for i in range(len(arr) - 1):
                p1, p2 = arr[i], arr[i+1]
                if (p1.x, p1.y) == (p2.x, p2.y):
                    continue
                # vertical segment
                if p1.x == p2.x:
                    # value in-between vertically increasing segment
                    if (val.y <= p1.y and val.y >= p2.y):
                        break
                # horizontal segment
                elif p1.y == p2.y:
                    # value in-between horizontally increasing segment
                    if (val.x <= p1.x and val.x >= p2.x):
                        break
                # sloped segment
                elif p1.x > p2.x:
                    if (val.x <= p1.x and val.x >= p2.x) \
                        and (val.y <= p1.y and val.y >= p2.y):
                        break
                    if (val.x <= p1.x and val.x >= p2.x) \
                        and (val.y >= p1.y and val.y <= p2.y):
                        break
                elif p1.x < p2.x:
                    if (val.x <= p2.x and val.x >= p1.x)\
                        and (val.y >= p2.y and val.y <= p1.y):
                        break
                elif i == len(arr.p) - 2:
                    index_found = False
            
            return i+1 if index_found else null

        def set_face(plate, south_point, north_point, location):
            """
            Trim a Line object adjacent to a divertor plate.
            The Line object will be defined by the interval of points
            [min_point, max_point] that form a subset of a target plate.
            """

            new_leg = order_plate_points(plate)
            i = get_insert_index(new_leg, south_point)
            j = get_insert_index(new_leg, north_point)

            if location == 'W':
                new_leg.insert(i, south_point)
                new_leg.insert(j, north_point)
            elif location == 'E':
                new_leg.insert(j, north_point)
                new_leg.insert(i, south_point)

            lookup = {}
            for i in range(len(new_leg)):
                lookup[new_leg[i]] = i

            start = lookup[north_point]
            end = lookup[south_point]

            if start > end:
                start, end = end, start

            # Return a reverse copy to maintain clockwise orientation within the LSN SNL patch.
            return Line([p for p in new_leg[start:end+1]]).reverse_copy()

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

        ITP = Line(order_plate_points(Line([Point(i) for i in self.itp])))
        OTP = Line(order_plate_points(Line([Point(i) for i in self.otp])))

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
        IDL_N = iPsiMax_TP
        IDL_S = xpt_ITP.reverse_copy()
        IDL_W = xptE_psiMax
        location = 'E'        
        IDL_E = set_face(ITP, IDL_S.p[0], IDL_N.p[-1], location = location)
        IDL = SNL_Patch([IDL_N, IDL_E, IDL_S, IDL_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # IPF Patch
        IPF_N = IDL_S.reverse_copy()
        IPF_S = psiMinPF_ITP.reverse_copy()
        IPF_W = xptS_psiMinPF.reverse_copy()
        location = 'E'
        IPF_E = set_face(ITP, IPF_S.p[0], IPF_N.p[-1], location = location)
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
        ODL_N = oPsiMax_TP.reverse_copy()
        ODL_S = xpt_OTP
        ODL_E = xptW_psiMax.reverse_copy()
        location = 'W'
        ODL_W = set_face(OTP, ODL_S.p[-1], ODL_N.p[0], location = location)
        ODL = SNL_Patch([ODL_N, ODL_E, ODL_S, ODL_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # OPF Patch
        OPF_N = ODL_S.reverse_copy()
        OPF_S = psiMinPF_OTP
        OPF_E = xptS_psiMinPF
        location = 'W'
        OPF_W = set_face(OTP, OPF_S.p[-1], OPF_N.p[0], location = location)
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
