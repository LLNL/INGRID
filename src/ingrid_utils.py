from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt
import pathlib
import inspect
from scipy.optimize import root, minimize

import yaml as yml
import os
from pathlib import Path
from time import time
from collections import OrderedDict

from OMFITgeqdsk import OMFITgeqdsk
from interpol import EfitData
from Root_Finder import RootFinder
from line_tracing import LineTracing
from geometry import Point, Line, Patch, segment_intersect, orientation_between


class IngridUtils():

    def __init__(self, settings: dict = {}, **kwargs):
        self.InputFile = None
        self.config = None
        self.settings_lookup = [
            'grid_settings',
            'integrator_settings',
            'target_plates',
            'DEBUG'
        ]
        self.GeometryAx = {
            'plate_W1': None,
            'plate_E1': None,
            'plate_W2': None,
            'plate_E2': None,
            'limiter': None
        }
        self.SetDefaultParams()
        self.process_yaml(settings)
        self.ProcessKeywords(**kwargs)

    def get_config(self):
        return self.LineTracer.config

    def SetDefaultParams(self):

        self.default_grid_settings = {
            'num_xpt': 1,
            'nlevs': 30,
            'full_domain': True,
            'psi_max': 0.0,
            'psi_core': 0.0,
            'psi_pf_1': 0.0,
            'psi_pf_2': 0.0,
            'psi_max_west': 0.0,
            'psi_max_east': 0.0,
            'rmagx': None,
            'zmagx': None,
            'rxpt': 0.0,
            'zxpt': 0.0,
            'rxpt2': 0.0,
            'zxpt2': 0.0,
            'grid_generation': {
                'np_default': 2,
                'nr_default': 2,
                'poloidal_f_default': 'x, x',
                'radial_f_default': 'x, x'
            },
            'patch_generation': {
                'rmagx_shift': 0.0,
                'zmagx_shift': 0.0,
                'west_tilt': 0.0,
                'east_tilt': 0.0,
                'use_NW': False,
                'use_NE': False,
                'use_secondary_NW': False,
                'use_secondary_NE': False,
                'use_SW': False,
                'use_SE': False,
                'NW_adjust': -0.785398,  # All values of pi / 4 radians.
                'NE_adjust': 0.785398,
                'secondary_NW_adjust': -0.785398,
                'secondary_NE_adjust': 0.785398,
            }
        }

        self.default_integrator_settings = {
            'dt': 0.01,
            'eps': 5e-5,
            'first_step': 1e-5,
            'step_ratio': 0.02,
            'tol': 5e-3
        }

        self.default_target_plate_settings = {
            'plate_E1': {
                'file': '',
                'name': '',
                'rshift': 0.0,
                'zshift': 0.0
            },

            'plate_E2': {
                'file': '',
                'name': '',
                'rshift': 0.0,
                'zshift': 0.0
            },

            'plate_W1': {
                'file': '',
                'name': '',
                'rshift': 0.0,
                'zshift': 0.0
            },

            'plate_W2': {
                'file': '',
                'name': '',
                'rshift': 0.0,
                'zshift': 0.0
            },
        }

        self.default_limiter_settings = {
            'file': '',
            'use_limiter': False,
            'use_efit_bounds': False,
            'use_default_data': True,
            'rshift': 0.0,
            'zshift': 0.0
        }

        self.default_patch_data_settings = {
            'file': '',
            'use_file': False,
            'preferences': {
                'new_file': False,
                'new_fname': ''
            }
        }

        self.default_DEBUG_settings = {
            'visual': {
                'find_NSEW': False,
                'patch_map': False,
                'subgrid': False,
                'gridue': False,
                'SF_analysis': False
            },

            'verbose': {
                'target_plates': False,
                'patch_generation': False,
                'grid_generation': False,
                'SF_analysis': False
            }
        }

        self.default_values_lookup = {
            'eqdsk': '',
            'grid_settings': self.default_grid_settings,
            'integrator_settings': self.default_integrator_settings,
            'target_plates': self.default_target_plate_settings,
            'limiter': self. default_limiter_settings,
            'patch_data': self.default_patch_data_settings,
            'DEBUG': self.default_DEBUG_settings
        }

        self.PlateData = {
            'plate_W1': {},
            'plate_E1': {},
            'plate_W2': {},
            'plate_E2': {}
        }

    def ProcessKeywords(self, **kwargs: str):
        for k, v in kwargs.items():
            if k == 'InputFile' or k == 'yaml':
                print('# Processing Input File:', v)
                self.InputFile = v
                self.process_yaml(self.ReadyamlFile(v))
                continue
            if k == 'W1TargetFile' or k == 'w1':
                self.settings['target_plates']['plate_W1']['file'] = v
                continue
            if k == 'E1TargetFile' or k == 'e1':
                self.settings['target_plates']['plate_E1']['file'] = v
                continue
            if k == 'E2TargetFile' or k == 'e2':
                self.settings['target_plates']['plate_E2']['file'] = v
                continue
            if k == 'W2TargetFile' or k == 'w2':
                self.settings['target_plates']['plate_W2']['file'] = v
                continue
            if k in ['LimiterFile', 'limiter', 'wall']:
                self.settings['limiter']['file'] = v
            if k == 'EqFile' or k == 'eq':
                self.settings['eqdsk'] = v
                continue
            print('Keyword "' + k + '" unknown and ignored...')

    def process_yaml(self, settings, verbose=False):
        """
        Parse the contents of a YAML dump (dictionary object).
        Parameter:
            - settings : dict
            Dictionary object conforming to INGRID YAML requirements.
        @author: garcia299
        """

        def get_default_values(item, sub_item=None, attribute=None):
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
                settings[item]
            except KeyError:
                settings[item] = get_default_values(item)
                if verbose:
                    print('Could not find "{}" in YAML file.'.format(item))
                    print('Populated "{}" with default value of "{}".\n'.format(item, settings[item]))
                continue

            # Second level entries within YAML file dump.
            if item == 'eqdsk':
                continue

            for sub_item in self.default_values_lookup[item].keys():
                try:
                    settings[item][sub_item]
                except KeyError:
                    settings[item][sub_item] = get_default_values(item, sub_item)
                    if verbose:
                        print('Could not find "{}/{}" in YAML file.'.format(item, sub_item))
                        print('Populated "{}/{}" with default value of "{}".\n'.format(item, sub_item, settings[item][sub_item]))
                    continue

                if item in ['grid_settings', 'target_plates'] \
                        and sub_item in ['patch_generation', 'grid_generation', 'plate_E1', 'plate_W1', 'plate_E2', 'plate_W2']:
                    for plate_attribute in self.default_values_lookup[item][sub_item].keys():
                        try:
                            settings[item][sub_item][plate_attribute]
                        except:
                            settings[item][sub_item][plate_attribute] = get_default_values(item, sub_item, plate_attribute)

        # Update references to YAML entries.
        self.settings = settings
        self.grid_settings = settings['grid_settings']
        self.integrator_settings = settings['integrator_settings']
        self.target_plates = settings['target_plates']
        self.DEBUG = settings['DEBUG']

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
        self.PsiUNorm = EfitData(rmin, rmax, nxefit,
                                 zmin, zmax, nyefit,
                                 rcenter, bcenter,
                                 rlimiter, zlimiter,
                                 rmagx, zmagx,
                                 name='Efit Data', parent=self)
        self.PsiUNorm.set_v(psi)

        if self.settings['grid_settings']['rmagx'] is None or self.settings['grid_settings']['zmagx'] is None:
            self.settings['grid_settings']['rmagx'], self.settings['grid_settings']['zmagx'] = (self.PsiUNorm.rmagx, self.PsiUNorm.zmagx)

        self.OMFIT_psi = g

    def calc_efit_derivs(self) -> None:
        """ Calculate the partial derivatives using finite differences
        for the un

        Wrapper for the member function in the EfitData class.
        """
        self.PsiUNorm.Calculate_PDeriv()

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
                        if point.count(',') > 0:
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

    def SetLimiter(self, coordinates=None, fpath=None, rshift=None, zshift=None):

        if rshift is None:
            rshift = self.settings['limiter']['rshift']
        if zshift is None:
            zshift = self.settings['limiter']['zshift']

        use_efit_bounds = self.settings['limiter']['use_efit_bounds']

        if fpath is not None:
            try:
                print('# Processing file for limiter data : {}'.format(self.settings['limiter']['file']))
                self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = self.ParseFileCoordinates(self.settings['limiter']['file'])
            except:
                raise ValueError(f"# Error in method 'SetLimiter' with fpath={fpath}")

        elif coordinates is not None:
            RLIM, ZLIM = coordinates
            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = [r - rshift for r in RLIM], [z - zshift for z in ZLIM]

        elif coordinates is None:
            g = OMFITgeqdsk(self.settings['eqdsk'])
            RLIM, ZLIM = g['RLIM'], g['ZLIM']
            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = RLIM - rshift, ZLIM - zshift

        if len(RLIM) == 0 or len(ZLIM) == 0:
            use_efit_bounds = True

        if use_efit_bounds:
            coordinates = [(self.PsiUNorm.rmin + 1e-2, self.PsiUNorm.zmin + 1e-2),
                           (self.PsiUNorm.rmax - 1e-2, self.PsiUNorm.zmin + 1e-2),
                           (self.PsiUNorm.rmax - 1e-2, self.PsiUNorm.zmax - 1e-2),
                           (self.PsiUNorm.rmin + 1e-2, self.PsiUNorm.zmax - 1e-2),
                           (self.PsiUNorm.rmin + 1e-2, self.PsiUNorm.zmin + 1e-2)]

            RLIM, ZLIM = [], []
            for c in coordinates:
                r, z = c
                RLIM.append(r)
                ZLIM.append(z)

            self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'] = np.array(RLIM) - rshift, np.array(ZLIM) - zshift

        self.LimiterData = Line([Point(p) for p in zip(self.OMFIT_psi['RLIM'], self.OMFIT_psi['ZLIM'])])

    def SetTargetPlate(self, settings, rshift=0, zshift=0):

        for _k, _v in settings.items():
            k = _k
            v = _v
            break

        if k.lower() in ['w1', 'westplate1', 'plate_w1']:
            k = 'plate_W1'
        elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
            k = 'plate_E1'
        elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
            k = 'plate_W2'
        elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
            k = 'plate_E2'
        else:
            raise ValueError(f"# Invalid key '{k}' provided for 'SetTargetPlate'")

        print(f"# Setting coordinates for '{k}'")

        R, Z = v

        # Make sure coordinates are unique
        data = np.array([c for c in zip(R, Z)])
        a, index = np.unique(data, return_index=True, axis=0)
        index.sort()
        self.PlateData[k] = Line([Point(x - rshift, y - zshift) for x, y in data[index]])

    def SetTargetPlates(self):
        for plate in self.settings['target_plates']:
            try:
                self.SetGeometry({plate: self.settings['target_plates'][plate]['file']},
                                 rshift=self.settings['target_plates'][plate]['rshift'],
                                 zshift=self.settings['target_plates'][plate]['zshift'])
            except:
                continue

    def OrderTargetPlates(self):
        for k in self.PlateData.keys():
            if type(self.PlateData[k]) == Line:
                self.OrderTargetPlate(k)

    def OrderTargetPlate(self, plate_key):
        """
        Ensures the target plate points are oriented clockwise around the
        magnetic axis.

        @author: garcia299
        """

        k = plate_key

        if not isinstance(self.PlateData.get(k), Line):
            raise ValueError(f"# Plate '{k}' is not loaded as a Line object. Check that 'SetGeometry({{'{k}' : ... }})' has been run.")

        plate = self.PlateData[k]

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
            v_reference = np.array([self.settings['grid_settings']['rmagx'],
                                    self.settings['grid_settings']['zmagx']])

            v_start = np.array([start.x, start.y])
            v_end = np.array([end.x, end.y])

            if orientation_between(v_start, v_end, v_reference) <= 0:
                ordered_plate = plate.copy()
            else:
                ordered_plate = plate.reverse_copy()

            # Gather raw data
            self.PlateData[k] = ordered_plate

    def OrderLimiter(self):
        """
        Ensures limiter geometry is oriented clockwise around magnetic axis
        """

        limiter = self.LimiterData.fluff_copy(100)

        rmid = (self.PsiUNorm.rmin + self.PsiUNorm.rmax) / 2
        zmid = (self.PsiUNorm.zmin + self.PsiUNorm.zmax) / 2

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
            v_reference = np.array([self.settings['grid_settings']['rmagx'],
                                    self.settings['grid_settings']['zmagx']])

            v_start = np.array([start.x, start.y])
            v_end = np.array([end.x, end.y])

            if orientation_between(v_start, v_end, v_reference) <= 0:
                ordered = True

        if ordered is False:
            self.LimiterData = self.LimiterData.reverse_copy()

    def FindMagAxis(self, x: float, y: float) -> None:
        # r_bounds = (self.PsiUNorm.rmin, self.PsiUNorm.rmax)
        # z_bounds = (self.PsiUNorm.zmin, self.PsiUNorm.zmax)
        # sol = minimize(fun=self.PsiUNorm.PsiFunction, x0=np.array([x,y]),
        #     method='L-BFGS-B', jac=self.PsiUNorm.Gradient, bounds=[r_bounds, z_bounds])
        sol = root(self.PsiUNorm.Gradient, [x, y])
        self.settings['grid_settings']['rmagx'] = sol.x[0]
        self.settings['grid_settings']['zmagx'] = sol.x[1]
        self.magx = (sol.x[0], sol.x[1])

    def FindXPoint(self, x: float, y: float) -> None:
        sol = root(self.PsiUNorm.Gradient, [x, y])
        self.settings['grid_settings']['rxpt'] = sol.x[0]
        self.settings['grid_settings']['zxpt'] = sol.x[1]
        self.xpt1 = (sol.x[0], sol.x[1])

    def FindXPoint2(self, x: float, y: float) -> None:
        sol = root(self.PsiUNorm.Gradient, [x, y])
        self.settings['grid_settings']['rxpt2'] = sol.x[0]
        self.settings['grid_settings']['zxpt2'] = sol.x[1]
        self.xpt2 = (sol.x[0], sol.x[1])

    def find_roots(self, tk_controller=None):
        """ Displays a plot, and has the user click on an approximate
        zero point. Uses a root finder to adjust to the more exact point.
        Right click to disable.
        Parameter:
            - tk_controller : Tk object
            A reference to the root TK object (for usage in gui mode).
        @author: watkins35, garcia299
        """

        self.root_finder = RootFinder(self.PsiUNorm, controller=tk_controller)

    def toggle_root_finder(self):
        """
        Activates or deactivates the root finder ability. Enables
        the user to save the location where they last clicked.
        @author: watkins35
        """
        self.root_finder.toggle_root_finding()

    def find_psi_lines(self, tk_controller=None):
        self.psi_finder = RootFinder(self.PsiNorm, mode='psi_finder', controller=tk_controller)

    def PrepLineTracing(self, interactive=True):
        """
        Initializes the line tracing class for the construction
        of the grid.
        Parameter:
            - refresh : boolean
            Re-initialize the LineTracing class.
        @author: watkins35, garcia299
        """

        self.LineTracer = LineTracing(self.PsiNorm, self.settings, interactive=interactive)
        if interactive:
            self.LineTracer.disconnect()

    def ClassifyTopology(self, visual=False):

        print('')
        print("# Begin classification....")
        print('')

        if self.settings['grid_settings']['num_xpt'] == 1:
            self.LineTracer.SNL_find_NSEW(self.xpt1, self.magx, visual)

        elif self.settings['grid_settings']['num_xpt'] == 2:
            self.LineTracer.DNL_find_NSEW(self.xpt1, self.xpt2, self.magx, visual)

        self.config = self.LineTracer.config

    def WriteGridueSNL(self, gridue_settings, fname='gridue'):

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
                        val = np.format_float_scientific(data[i][j][n], precision=15, unique=False).rjust(23).replace('e', 'D')
                        if delim_val == 3:
                            delim_val = 0
                            delim_char = '\n'
                        body += val + delim_char
                        delim_char = ''

            if delim_val % 3 != 0:
                body += '\n'

            return body

        f = open(fname, mode='w')
        f.write(format_header(gridue_settings) + '\n')

        body_items = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
        for item in body_items:
            f.write(format_body(gridue_settings[item]) + '\n')

        runidg = 'iogridue'
        f.write(runidg + '\n')

        f.close()

        return True

    def WriteGridueDNL(self, gridue_settings, fname='gridue'):

        def format_header(gridue):
            header_rows = [
                ['nxm', 'nym'],
                ['iyseparatrix1', 'iyseparatrix2'],
                ['ix_plate1', 'ix_cut1', '_FILLER_', 'ix_cut2', 'ix_plate2'],
                ['iyseparatrix3', 'iyseparatrix4'],
                ['ix_plate3', 'ix_cut3', '_FILLER_', 'ix_cut4', 'ix_plate4']
            ]

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
                        val = np.format_float_scientific(data[i][j][n], precision=15, unique=False).rjust(23).replace('e', 'D')
                        if delim_val == 3:
                            delim_val = 0
                            delim_char = '\n'
                        body += val + delim_char
                        delim_char = ''

            if delim_val % 3 != 0:
                body += '\n'

            return body

        f = open(fname, mode='w')
        f.write(format_header(gridue_settings) + '\n')

        body_items = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
        for item in body_items:
            f.write(format_body(gridue_settings[item]) + '\n')

        runidg = 'iogridue'
        f.write(runidg + '\n')

        f.close()
        return True

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

            patch = Patch([N, E, S, W], patchName=patch_settings['patchName'],
                          platePatch=patch_settings['platePatch'],
                          plateLocation=patch_settings['plateLocation'],
                          PatchTagMap=patch_settings['PatchTagMap'])

            patches[patch.patchName] = patch

        patches = OrderedDict([(k, v) for k, v in patches.items()])

        return config, xpt_data, patches

    def CheckPatches(self, verbose=False):
        self.CurrentTopology.CheckPatches(verbose=verbose)

    def PrepGridue(self):
        self.CurrentTopology.SetupPatchMatrix()
        self.CurrentTopology.concat_grid()
        self.CurrentTopology.set_gridue()
