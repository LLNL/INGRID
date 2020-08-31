#!/usr/bin/env pythonu
# -*- coding: utf-8 -*-
"""Ingrid module for interfacing all grid generator capabilities.

This module contains the Ingrid class that drives all code functionality.

The Ingrid class is to encapsulate all functionality and allow the user to
easily take advantage of advanced features of the code.

"""
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
from interpol import Bicubic
from ingrid_utils import IngridUtils
from topologies.snl import SNL
from topologies.sf15 import SF15
from topologies.sf45 import SF45
from topologies.sf75 import SF75
from topologies.sf105 import SF105
from topologies.sf135 import SF135
from topologies.sf165 import SF165
from topologies.udn import UDN
from Root_Finder import RootFinder
from line_tracing import LineTracing
from geometry import Point, Line


def QuickStart() -> None:
    """ Start Ingrid in gui mode with default settings.

    Rather than providing data to the class constructor, a user can opt to
    start Ingrid immediately in it's gui form. This is useful for users who
    are not familiar with the Ingrid class and it's capabiliites. Advanced
    users may still find QuickStart useful, but would also use the code in
    self-authored scripts.

    """
    QuickGrid = Ingrid()
    QuickGrid.SimpleGUI()


class Ingrid(IngridUtils):
    """Ingrid class for managing the grid generator code.

    Child of class 'IngridUtils'. Class 'Ingrid' is to contain attributes
    and methods utilized by beginner code users.

    Args:
        settings (optional): Dictionary representation of the settings file.

            The dictionary is obtained via a YAML dump while operating in gui
            mode. Providing no dictionary will populate the Ingrid object
            with a default value settings attribute. Any missing entries
            provided by a user will be populated with default values by Ingrid.

            (Refer to the Ingrid "Settings File Format" for a complete descri-
            ption of the settings dictionary)

        **kwargs: Keyword arguments for processing of input data.
            (See method 'ProcessKeywords' in class 'IngridUtils' for a
            list of supported keywords)

    Attributes:
        settings (dict): Dictionary representation of the settings file.

        PlateData (dict): Dictionary mapping target plates to their
            corresponding class 'Line' objects.

        LimiterData (:obj:`Line`): 'Line' class object representing
            tokamak limiter.

        magx (tuple): (R, Z) coordinates of magnetic-axis (float entries).

        xpt1 (tuple): (R, Z) coordinates of primary x-point (float entries).

        xpt2 (tuple): (R, Z) coordinates of secondary x-point (float entries).

        PsiUNorm (:obj:'Efit'): 'EfitData' class object containing unormalized
            psi data from provided neqdsk file in settings.

        PsiNorm (:obj:'Efit'): 'EfitData' class object containing normalized
            psi data calculated from class attribute 'PsiUNorm'.

        CurrentTopology (:obj:'Topology'): Tokamak magnetic configuration
            type currently being operated on. Class 'Topology' acts as the
            base class for all magnetic configurations.

    """

    def __init__(self, settings: dict = {}, **kwargs):
        IngridUtils.__init__(self, settings, **kwargs)
        self.PrintSummaryInput()

    def LoadEFIT(self, fpath: str) -> None:
        self.settings['eqdsk'] = fpath
        self.OMFIT_read_psi()
        self.calc_efit_derivs()

    def SimpleGUI(self) -> None:
        """Start GUI for Ingrid.

        Will assume usage on a machine with tk GUI capabilities.
        No prior settings file is required as the user will be
        prompted with an option to generate a new file.
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
                self._IngridWindow.destroy()

        from gui.simple_gui import IngridGUI
        self.IngridWindow = IngridGUI(IngridSession=self)
        self.IngridWindow.title('INGRID')
        self.IngridWindow.protocol('WM_DELETE_WINDOW', on_closing)
        self.IngridWindow.mainloop()

    def SaveSettingsFile(self, fname: str = '', settings: dict = {}) -> 'Path':
        """Save a new settings .yml file.

        Args:
            fname (optional): Name of new settings '.yml' file.
                If default value of '', then Ingrid will generate
                a '.yml' file named 'INGRID_Session' appended with
                a timestamp.

            settings (optional): Ingrid settings dictionary to dump into
                the '.yml' file. Defaults to empty dict which produces a
                template settings file.

        Returns:
            Path (:obj: 'Path'): Path object of saved file.
        """

        if fname == '':
            fname = 'INGRID_Session' + str(int(time())) + '.yml'

        fpath = Path(fname)

        if fpath.suffix == '':
            fname += '.yml'
        elif fpath.suffix != '.yml':
            fname = fname[:-4] + '.yml'

        self.process_yaml(settings)

        with open(fname, 'w') as f:
            yml.dump(self.settings, f)

        return fpath

    def PrintSummaryInput(self) -> None:
        """ Print a summary of the currently loaded data files

        Will print relevant EQDSK, patch data files, target plate files,
        and limiter files.
        """
        print('')
        print('Equilibrium File:', self.settings['eqdsk'])

        if (self.settings['patch_data']['use_file']) \
            and (Path(self.settings['patch_data']['file']).is_file()) \
                and (Path(self.settings['patch_data']['file']).suffix == '.npy'):

            print('Patch data-file:')
            print(' # Patch data-file:', self.settings['patch_data']['file'])

        if self.settings['limiter']['use_limiter']:
            print('Limiter File:')

            if self.settings['limiter']['file'] == '':
                print(' # Limiter:', 'Using default eqdsk (RLIM, ZLIM) coordinates')
            else:
                print(' # Limiter', self.settings['limiter']['file'])
        else:
            print('Target Files:')
            print(' # W1:', self.settings['target_plates']['plate_W1']['file'])
            print(' # E1:', self.settings['target_plates']['plate_E1']['file'])
            print(' # E2:', self.settings['target_plates']['plate_E2']['file'])
            print(' # W2:', self.settings['target_plates']['plate_W2']['file'])
        print('')

    def PrintSummaryParams(self) -> None:
        """ Print a summary of key settings values.
        """
        print('')
        print('Summary:')
        print(' # Number of x-points:', self.settings['grid_settings']['num_xpt'])
        print(' # Use full domain:', self.settings['grid_settings']['full_domain'])
        print(' # Use limiter:', self.settings['limiter']['use_limiter'])
        print(' # Use patch data-file:', self.settings['patch_data']['use_file'])
        print('')
        self.PrintSummaryInput()

    def SetGeometry(self, geo_items: dict, rshift: float = 0.0, zshift: float = 0.0) -> None:
        """ Define and load the tokamak strike geometry into the Ingrid object.

        Allows the user to set target plate data and/or limiter data that
        will be used to generate a patch map.

        Args:
            geo_items: Argument dict specifying which item(s) to set and by
                what way. Multiple geometry items can be set at once, but
                the following key-value format must be obeyed for any number
                of entries.

                All keys for 'geo_items' are type 'str'. Accepted key
                values are as follows:

                Geometry:  Key value:
                Plate W1 - 'plate_W1', 'W1'
                Plate E1 - 'plate_E1', 'E1'
                Plate W2 - 'plate_W2', 'W2'
                Plate E2 - 'plate_E2', 'E2'
                Limiter  - 'limiter', 'wall'

                Note: the above keys are NOT case sensitive.


                Corresponding key values can vary in data type depending
                on means of setting geometry. The types and their usage
                are as follows:


                'str': Path to '.txt' file containing coordinate data or
                to Ingrid generated '.npy' file (obtained via methods
                'SaveLimiterData' and 'SaveTargetPlateData').

                Note: When setting Limiter geometry, the user can provide the
                value 'default' to set limiter data to that which is contained
                in the loaded neqdsk file.


                'list', 'tuple': Iterables provided as values must be
                of length == 2. The first entry corresponds to R coordinate
                information, and the second entry corresponds to Z coordinate
                information. This information can be in the form of a list or
                NumPy array with shape == (N,).


                'dict': One can map to a dictionary taking on a variety of
                formats depending on need.

                Setting geometry with explicit RZ coordinates:
                {'R': R_data, 'Z': z_data}
                Where R_data and Z_data are in the form of a list or NumPy
                array with shape == (N,) as above.

                Setting geometry with data from external file:
                {'file': str}

                Notes: The above dict option support keys 'rshift' and 'zshift'
                for the user to provide transformations to individual geometry
                items (see examples).

                Because the core 'settings' attribute contained by the
                Ingrid class contains dict structures itself, the user can also
                provide settings['limiter'] and settings['target_plates'][k]
                to method 'SetGeometry' (where k corresponds to a plate key).


        rshift (optional): Translate 'geo_items' to coordinate R',
            with R' = R - rshift. Will override all 'rshift' provided entries
            in 'geo_items'.

        zshift (optional): Translate 'geo_items' to coordinate Z',
            with Z' = Z - zshift. Will override all 'zshift' provided entries
            in 'geo_items'.


        Examples:
            Setting default limiter data contained in loaded neqdsk
            file:

            >>> MyIG = Ingrid()
            >>> MyIG.SetGeometry({'limiter': 'default'})


            Setting target plate 'E1' with numpy array:

            >>> MyIG = Ingrid()
            >>> MyIG.SetGeometry({'E1': 'E1_data.npy'})


            Setting limiter data with Ingrid '.npy' file:

            >>> MyIG = Ingrid()
            >>> MyIG.SetGeometry({'limiter': 'LimiterData.npy'})


            Setting both target plates 'W1' and 'E1' with '.txt'
            and '.npy' files while only applying rshift and zshift
            to target plate 'E1'
            (Note the dict structure used for specifying 'E1'):

            >>> MyIG = Ingrid()
            >>> geometry_dict = {
            ... 'W1': 'W1_data.txt',
            ... 'E1': {
            ... 'file': 'E1_data.npy',
            ... 'rshift': 1.23, 'zshift': 3.14
            ... }
            ... }
            >>> MyIG.SetGeometry(geometry_dict)


            Setting plate 'W2' with user provided NumPy array
            for RZ coordinates:
            >>> MyIG = Ingrid()
            >>> R_data = my_numpy_array_1
            >>> Z_data = my_numpy_array_2
            >>> RZ_dict = {'R': my_numpy_array_1, 'Z': my_numpy_array_2}
            >>> MyIG.SetGeometry({'plate_W1': RZ_dict})

        Raises:
            ValueError: If file path provided does not lead to actual file.
            ValueError: If file provided is not of suffix '.txt' or '.npy'.
            ValueError: If invalid 'geo_items' key provided.
            ValueError: If 'geo_items' dict value contains invalid key.
            ValueError: If value associated with 'geo_items' key is not of
                supported data type.
        """

        def _parse_v_for_file_values(k_str, v_dict):
            if v_dict.get('file') is not None:

                if Path(v['file']).is_file() is False:
                    raise ValueError(f"# File '{v_dict['file']}'' provided in settings dict corresponding to '{k_str}' does not exist.")

                if Path(v_dict['file']).suffix == '.txt':
                    x_ret_val, y_ret_val = self.ParseFileCoordinates(v_dict['file'])
                elif Path(v_dict['file']).suffix == '.npy':
                    x_ret_val, y_ret_val = np.load(v_dict['file'])
                else:
                    raise ValueError(f"# Invalid file extension '{Path(v_dict['file']).suffix}' used in SetGeometry.")

                return (x_ret_val, y_ret_val)

        def _parse_v_for_xy_coordinates(v_dict):
            x_flag = False
            y_flag = False
            if v_dict.get('x') is not None:
                x_ret_val = v_dict['x']
                x_flag = True
            if v_dict.get('X') is not None:
                x_ret_val = v_dict['X']
                x_flag = True
            if v_dict.get('r') is not None:
                x_ret_val = v_dict['r']
                x_flag = True
            if v_dict.get('R') is not None:
                x_ret_val = v_dict['R']
                x_flag = True
            if v_dict.get('y') is not None:
                y_ret_val = v_dict['y']
                y_flag = True
            if v_dict.get('Y') is not None:
                y_ret_val = v_dict['Y']
                y_flag = True
            if v_dict.get('z') is not None:
                y_ret_val = v_dict['z']
                y_flag = True
            if v_dict.get('Z') is not None:
                y_ret_val = v_dict['Z']
                y_flag = True

            if x_flag is True and y_flag is True:
                return (x_ret_val, y_ret_val)
            else:
                return None

        for k, v in geo_items.items():

            TargetPlate = False
            Limiter = False
            DefaultEQ = False

            if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                TargetPlate = True
                k = 'plate_W1'
            elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                TargetPlate = True
                k = 'plate_E1'
            elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                TargetPlate = True
                k = 'plate_W2'
            elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                TargetPlate = True
                k = 'plate_E2'
            elif k.lower() in ['limiter', 'wall']:
                Limiter = True
                k = 'limiter'
            else:
                raise ValueError(f"# Invalid key '{k}' used in SetGeometry.")

            if type(v) == str:  # Path to file.
                if Limiter and v.lower() in ['default', 'eq']:
                    DefaultEQ = True
                elif Path(v).suffix == '.txt':
                    x, y = self.ParseFileCoordinates(v)
                elif Path(v).suffix == '.npy':
                    x, y = np.load(v)
                else:
                    raise ValueError(f"# Invalid file extension '{Path(v).suffix}' used in SetGeometry.")

            elif type(v) == dict:  # Dict of settings.

                xy_result = _parse_v_for_xy_coordinates(v)

                if xy_result is None:  # Dict structure with settings file keys.
                    if TargetPlate:
                        if v.get('file') not in ['', None]:  # Attempt to read file.
                            x, y = _parse_v_for_file_values(k, v)
                        else:
                            raise ValueError(f"# No data source for target plate '{k}' identified.")
                    elif Limiter:
                        if v.get('file') not in ['', None]:  # Attempt to read file.
                            x, y = _parse_v_for_file_values(k, v)
                        else:
                            DefaultEQ = True
                elif type(xy_result) == tuple:
                    x, y = xy_result

                else:
                    raise ValueError(f"# Invalid key '{k}' was used within dic in SetGeometry.")

                if v.get('rshift') is not None and rshift == 0.0:
                    rshift = v['rshift']
                if v.get('zshift') is not None and zshift == 0.0:
                    zshift = v['zshift']

            elif type(v) in [list, tuple]:
                x, y = v

            else:
                raise ValueError(f"# Invalid data type {type(v)} was used for setting geometry.")

            if TargetPlate:
                self.SetTargetPlate({k: [x, y]}, rshift=rshift, zshift=zshift)
            elif Limiter and not DefaultEQ:
                self.SetLimiter([x, y], rshift=rshift, zshift=zshift)
            elif Limiter and DefaultEQ:
                self.SetLimiter(rshift=rshift, zshift=zshift)
            else:
                raise ValueError(f"# Invalid key '{k}' was used for setting geometry.")

    def SaveGeometryData(self, geo_items: dict, timestamp: bool = False) -> None:
        """ Save strike geometry Line object RZ coordinates as '.npy' file.

        Args:
            geo_items: Argument dict specifying which target plate to save
                and the file name/path

                {geo_name: data_fname, ... }

                Where both geo_name and data_fname are of type str. Multiple
                data_fname files can be saved at once in a manner similar to
                method 'SetGeometry' (see documentation for 'SetGeometry').

            timestamp (optional): Append a time stamp to the end of the files.

        Raises:
            ValueError: If 'geo_items' is not type 'dict'.
            ValueError: If invalid 'geo_items' key provided.
            ValueError: If data_fname entry is not of type 'str'.
            ValueError: If data_fname entry is an empty string.
            ValueError: If requested strike geometry Line to save has no data.
        """

        if type(settings) is dict:

            if timestamp is True:  # Get same timestamp for all files.
                f_tail = '_' + str(int(time()))
            else:
                f_tail = ''

            for k, v in settings.items():
                if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                    k = 'plate_W1'
                elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                    k = 'plate_E1'
                elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                    k = 'plate_W2'
                elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                    k = 'plate_E2'
                elif k.lower() in ['limiter', 'wall']:
                    k = 'limiter'
                else:
                    raise ValueError(f"# Invalid key '{k}' provided for 'SaveGeometryData'")

                if type(v) is not str:
                    raise ValueError(f"# Value {v} associated with key '{k}' must be of type str")

                fpath = Path(v)

                if fpath.name.strip() == '':
                    raise ValueError(f"# Cannot save '{k}' to empty file name.")

                fname = fpath.fname

                if k == 'limiter':
                    geo = self.LimiterData
                else:
                    geo = self.PlateData[k]

                if type(geo) is Line:
                    R = np.array(geo.xval)
                    Z = np.array(geo.yval)
                else:
                    raise ValueError(f"# No Line object has been instantiated for {k}")

                fname += f_tail  # Add on f_tail value from start.

                np.save(fname, np.array([R, Z]))
                print(f"# Saved '{k}' plate data to file " + {fname} + ".npy")

        else:
            raise ValueError(f"# Argument 'geo_items' must be of type dict.")

    def LoadGeometryData(self, geo_items: dict) -> None:
        """ Load strike geometry RZ coordinates from external file.

        Args:
            geo_items: Argument dict specifying which target plate to load
                data into. Format is as follows:

                {geo_name: data_fname, ... }

                Where both geo_name and data_fname are of type str. Multiple
                data_fname files can be loaded at once in a manner similar to
                method 'SetGeometry' (see documentation for 'SetGeometry').

        Raises:
            ValueError: If 'geo_items' is not type 'dict'.
            ValueError: If invalid 'geo_items' key provided.
            ValueError: If data_fname entry is not of type 'str'.
            ValueError: If data_fname is not a file.
            ValueError: If data_fname file is not of format '.txt' or '.npy'
        """
        if type(geo_items) is dict:

            for k, v in geo_items.keys():
                if k.lower() in ['w1', 'westplate1', 'plate_w1']:
                    k = 'plate_W1'
                elif k.lower() in ['e1', 'eastplate1', 'plate_e1']:
                    k = 'plate_E1'
                elif k.lower() in ['w2', 'westplate2', 'plate_w2']:
                    k = 'plate_W2'
                elif k.lower() in ['e2', 'eastplate2', 'plate_e2']:
                    k = 'plate_E2'
                elif k.lower() in ['limiter', 'wall']:
                    k = 'limiter'
                else:
                    raise ValueError(f"# Invalid key '{k}' provided for 'LoadGeometryData'")

                if type(v) is not str:
                    raise ValueError(f"# fname associated with key '{k}' must be of type {type('')}. Received type {type(v)}")

                fpath = Path(v)

                if fpath.is_file() is False:
                    raise ValueError(f"# fname associated with key '{k}' does not exist.")

                suffix = fpath.suffix
                fname = fpath.name

                if suffix == '.npy':
                    R, Z = np.load(fname)
                elif suffix == '.txt':
                    R, Z = self.ParseFileCoordinates(fname)
                else:
                    raise ValueError(f"# File type '{suffix}' is not supported (requires '.txt', '.npy')")

                if k == 'limiter':
                    destination = self.LimiterData
                else:
                    destination = self.PlateData[k]

                if R[0] is None or Z[0] is None:
                    destination = None
                else:
                    destination = Line([Point(P) for P in zip(R, Z)])

                print(f"# Loaded '{k}' data from file " + fname)

        else:
            raise ValueError(f"# LoadGeometryData input must be of type dict with format {{GEO_KEY : FNAME}}")

    def PlotStrikeGeometry(self) -> None:
        """ Plot all strike geometry to be used for drawing the Patch Map.
        """
        try:
            [ax_line.remove() for StrikeGeometry in self.GeometryAx.values() for ax_line in StrikeGeometry.lines]
        except:
            pass

        try:
            if self.settings['grid_settings']['num_xpt'] == 2 or self.settings['limiter']['use_limiter']:
                self.PlotLimiter()
                if self.settings['limiter']['use_limiter'] is False:
                    self.PlotTargetPlates()
            else:
                self.PlotTargetPlates()
        except:
            self.PlotTargetPlates()

    def PlotTargetPlates(self) -> None:
        """ Plot all PlateData and remove outdated plate line artists.
        """
        plate_colors = ['blue', 'orange', 'firebrick', 'green']
        for k, c in zip([key for key in self.PlateData.keys()], plate_colors):
            self.PlotTargetPlate(plate_key=k, color=c)

    def PlotTargetPlate(self, plate_key: str, color: str = 'red') -> None:
        """ Plot a target plate corresponding to a plate key.

        Args:
            plate_key: An Ingrid supported target plate key (see method
                SetGeometry for supported plate keys)

            color (optional): Color to provide to matplotlib
        """

        if plate_key in [k for k in self.PlateData.keys()]:
            try:
                [ax_line.remove()
                    for ax_line in self.GeometryAx[plate_key].lines
                    if ax_line.get_label() == plate_key]
            except:
                pass
            if type(self.PlateData[plate_key]) is Line:
                self.GeometryAx[plate_key] = self.PlateData[plate_key].plot(label=plate_key, color=color)
            else:
                pass
        else:
            raise ValueError(f"# Value '{plate_key}' is not an Ingrid supported target plate key.")

    def PlotLimiter(self) -> None:
        """ Plot limiter geometry.
        """

        try:
            [ax_line.remove() for ax_line in self.GeometryAx['limiter'].lines if ax_line.get_label() == 'limiter']
        except:
            pass
        self.GeometryAx['limiter'] = self.LimiterData.plot(color='dodgerblue', label='limiter')

    def PlotPsiUNorm(self) -> None:
        """ Plot unnormalized psi data.
        """
        try:
            if self.PsiUNorm.ax.__dict__['collections'] is not None:
                self.PsiUNorm.ax.collections = []
                plt.draw()
        except:
            pass
        self.PsiUNorm.plot_data(self.settings['grid_settings']['nlevs'])

    def PlotPsiNorm(self) -> None:
        """ Plot normalized psi data.
        """
        try:
            if self.PsiNorm.ax.__dict__['collections'] is not None:
                self.PsiNorm.ax.collections = []
                plt.draw()
        except:
            pass
        self.PsiNorm.plot_data(self.settings['grid_settings']['nlevs'])

    def PlotPsiNormBounds(self) -> None:
        """ Plot contour lines associated with psi boundary values provided.

        This method extracts psi values from the 'settings' dict and plots the
        psi level. In addition to the psi values in 'settings', the primary
        and, if applicable, secondary separatrix are plotted as well.

        If the user is operating on a single null configuration, the psi values
        plotted are 'psi_max', 'psi_core', 'psi_pf_1'.

        If the user is operating on a case with two x-points, the psi values
        are the same as above but with 'psi_max_west', 'psi_max_east', and
        'psi_pf_2' also included in the plot.
        """
        nxpt = self.settings['grid_settings']['num_xpt']
        if nxpt == 1:
            Dic = {'psi_max': 'lime',
                   'psi_core': 'cyan',
                   'psi_pf_1': 'white'}
        elif nxpt == 2:
            Dic = {'psi_core': 'cyan',
                   'psi_max_west': 'magenta',
                   'psi_max_east': 'lime',
                   'psi_pf_1': 'white',
                   'psi_pf_2': 'yellow'}

        for k, c in Dic.items():
            self.PsiNorm.PlotLevel(self.settings['grid_settings'][k], color=Dic[k], label=k)

        self.PsiNorm.PlotLevel(1.0, color='red', label='Primary Separatrix')
        if nxpt == 2:
            self.PsiNorm.PlotLevel(
                self.PsiNorm.get_psi(self.xpt2[0], self.xpt2[1]), color='blue', label='Secondary Separatrix'
            )

        self.PsiNorm.ax.legend(bbox_to_anchor=(0.5, -0.15), loc='lower center', ncol=len(Dic) // 2)

    def PlotPsiNormMagReference(self) -> None:
        """ Plot a marker on the magnetic axis and all x-points of interest.
        """
        (x, y) = self.magx
        self.PsiNorm.ax.plot(x, y, '+', color='yellow', ms=15, linewidth=5)
        (x, y) = self.xpt1
        self.PsiNorm.ax.plot(x, y, '+', color='white', ms=15, linewidth=5)
        if self.settings['grid_settings']['num_xpt'] == 2:
            try:
                (x, y) = self.xpt2
                self.PsiNorm.ax.plot(x, y, '+', color='white', ms=15, linewidth=5)
            except:
                pass

    def PlotTopologyAnalysis(self) -> None:
        """ Shade the private flux, core, and show where the secondary x-point
        travels.

        This method can be used to interpret which type of configuration the
        user is handling.
        """
        self.PsiNorm.ax.add_patch(self.LineTracer.RegionPolygon['core'])
        self.PsiNorm.ax.add_patch(self.LineTracer.RegionPolygon['pf'])
        self.LineTracer.RegionLineCut.plot(color='white')

    def PlotPsiLevel(self, efit_psi: object, level: float, Label: str = '') -> None:
        """ Plot a contour corresponding to a psi level.

        Args:
            efit_psi: The 'EfitData' object to get the psi data from.

            level: The psi value to plot.

            Label (optional): Label to provide to matplotlib.pyplot.contour.
        """
        plt.contour(efit_psi.r, efit_psi.z, efit_psi.v, [level], colors='red', label=Label)

    def PlotMidplane(self) -> None:
        """ Plot midplane line through magnetic axis with any applied
        transformations specified in settings.

        This method can be used to inspect the effects of 'west_tilt',
        'east_tilt', 'rmagx_shift', and 'zmagx_shift'.
        """
        try:
            west_tilt = self.settings['grid_settings']['patch_generation']['west_tilt']
        except KeyError:
            west_tilt = 0.0
        try:
            east_tilt = self.settings['grid_settings']['patch_generation']['east_tilt']
        except KeyError:
            east_tilt = 0.0

        R = self.settings['grid_settings']['rmagx'] + self.settings['grid_settings']['patch_generation']['rmagx_shift']
        Z = self.settings['grid_settings']['zmagx'] + self.settings['grid_settings']['patch_generation']['zmagx_shift']
        magx = Point(np.array([R, Z]))

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx.x - 1e6 * np.cos(west_tilt), magx.y - 1e6 * np.sin(west_tilt))
        RHS_Point = Point(magx.x, magx.y)
        west_midLine = Line([LHS_Point, RHS_Point])
        west_midLine.plot('white')

        LHS_Point = Point(magx.x, magx.y)
        RHS_Point = Point(magx.x + 1e6 * np.cos(east_tilt), magx.y + 1e6 * np.sin(east_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])
        east_midLine.plot('white')

    def PlotPatches(self) -> None:
        """ Plot the patch map that was generated with method 'CreatePatches'
        """
        try:
            plt.close(self._PatchFig)
        except:
            pass
        self._PatchFig = plt.figure('INGRID: ' + self.CurrentTopology.config + ' Patches', figsize=(6, 10))
        self.PatchAx = self._PatchFig.add_subplot(111)
        self.CurrentTopology.patch_diagram(fig=self._PatchFig, ax=self.PatchAx)
        self.PlotStrikeGeometry()

    def CloseOpenPlots(self):
        """ Close all open matplotlib figures.
        """
        try:
            plt.close('all')
        except:
            pass

    def AutoRefineMagAxis(self) -> None:
        """ Refine magnetic axis RZ coordinates.

        Will apply a root_finder method to 'rmagx' and 'zmagx' float values
        stored in attribute "settings['grid_settings']".
        """
        self.FindMagAxis(self.settings['grid_settings']['rmagx'], self.settings['grid_settings']['zmagx'])

    def AutoRefineXPoint(self) -> None:
        """ Refine primary x-point RZ coordinates.

        Will apply a root_finder method to 'rxpt' and 'zxpt' float values
        stored in attribute "settings['grid_settings']".
        """
        self.FindXPoint(self.settings['grid_settings']['rxpt'], self.settings['grid_settings']['zxpt'])

    def AutoRefineXPoint2(self) -> None:
        """ Refine secondary x-point RZ coordinates.

        Will apply a root_finder method to 'rxpt2' and 'zxpt2' float values
        stored in attribute "settings['grid_settings']".
        """
        self.FindXPoint2(self.settings['grid_settings']['rxpt2'], self.settings['grid_settings']['zxpt2'])

    def SetMagReference(self: object, topology: str = 'SNL') -> None:
        self.magx = (self.settings['grid_settings']['rmagx'], self.settings['grid_settings']['zmagx'])
        self.xpt1 = (self.settings['grid_settings']['rxpt'], self.settings['grid_settings']['zxpt'])

        if self.settings['grid_settings']['num_xpt'] == 2:
            self.xpt2 = (self.settings['grid_settings']['rxpt2'], self.settings['grid_settings']['zxpt2'])

    def CalcPsiNorm(self) -> None:
        """ Normalize psi data to a refined magnetic axis and primary x-point
        """

        # use the same bounds as the efit data
        self.PsiNorm = EfitData(self.PsiUNorm.rmin, self.PsiUNorm.rmax,
                                self.PsiUNorm.nr, self.PsiUNorm.zmin,
                                self.PsiUNorm.zmax, self.PsiUNorm.nz,
                                self.PsiUNorm.rcenter, self.PsiUNorm.bcenter,
                                self.PsiUNorm.rlimiter, self.PsiUNorm.zlimiter,
                                self.PsiUNorm.rmagx, self.PsiUNorm.zmagx,
                                name='Normalized Efit Data', parent=self)
        psi = self.PsiUNorm.v
        psi_magx = self.PsiUNorm.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.PsiUNorm.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx)) / (psi_xpt1 - psi_magx)
        self.PsiNorm.set_v(psinorm)
        self.PsiNorm.Calculate_PDeriv()

    def AnalyzeTopology(self) -> None:
        """ Perform analysis on normalized psi to determine magnetic topolgy.
        """

        try:
            visual = self.settings['DEBUG']['visual']['find_NSEW']
        except:
            visual = False

        self.PrepLineTracing(interactive=False)
        self.ClassifyTopology(visual=visual)

        print('')
        print('# Identified {} configuration.'.format(self.config))
        print('')

        self.SetTopology(self.config)

    def SetTopology(self, topology: str) -> None:
        """ Initialize the current topology to a particular magnetic topology.

        Args:
            topology: String literal corresponding to which magnetic topology
                to initialize.

                Values can be:
                    'LSN': Lower Single Null
                    'USN': Upper Single Null
                    'UDN': Unbalanced Double Null
                    'SF15': Snowflake-15
                    'SF45': Snowflake-45
                    'SF75': Snowflake-75
                    'SF105': Snowflake-105
                    'SF135': Snowflake-135
                    'SF165': Snowflake-165

        Raises:
            ValueError: If user provided unrecognized string entry.
        """
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

        else:
            raise ValueError(f"# Invalid entry: no Ingrid topology {topology} exists.")

        self.CurrentTopology = ingrid_topology

    def StartSetup(self, **kwargs) -> None:
        """ A collection of essential tasks before generating a patch map from
        scratch.

        The user should ensure that the 'settings' dict is populated with
        the correct paths to relevant neqdsk data and geometry files.

        The user should ensure 'settings[''grid_settings''][''num_xpt'] is
        set with the correct integer value.
        """

        if self.settings['grid_settings']['num_xpt'] == 1:
            topology = 'SNL'
        elif self.settings['grid_settings']['num_xpt'] == 2:
            topology = 'DNL'
        else:
            v_error_str = '# Invalid number of x-points.' \
                + f'(User provided {self.settings["grid_settings"]["num_xpt"]}.' \
                + 'Must be <= 2).'
            raise ValueError(v_error_str)

        self.LoadEFIT(self.settings['eqdsk'])
        self.AutoRefineMagAxis()
        self.AutoRefineXPoint()
        if topology == 'DNL':
            self.AutoRefineXPoint2()
        if self.settings['limiter']['use_limiter']:
            self.SetGeometry({'limiter': self.settings['limiter']})
        self.SetTargetPlates()
        self.SetMagReference(topology)
        self.CalcPsiNorm()
        self.PrepLineTracing()

    def ShowSetup(self) -> None:
        """ Show Ingrid setup that a patch map will be generated from.

        This method plots normalized psi data, psi boundaries, strike geometry,
        and midplane lines through the magnetic axis.
        """
        self.CloseOpenPlots()
        self.PlotPsiNorm()
        self.PlotPsiNormMagReference()
        self.PlotStrikeGeometry()
        self.PlotPsiNormBounds()
        self.PlotMidplane()
        self.PrintSummaryParams()

    def ConstructPatches(self) -> None:
        """ Create a patch map that can be refined into a grid.

        This method assumes the user has either loaded patch data
        from a previous Ingrid session (see 'LoadPatches'), or that
        the user has already successfully called method 'AnalyzeTopology'.

        Should the user want to automatically enable patch saving, the user
        should set the entry

        'settings[''patch_data''][''preferences''][''new_file'']'

        with a value of 'True'.

        """
        self.CurrentTopology.RefreshSettings()
        self.OrderLimiter()
        self.OrderTargetPlates()
        self.PrepLineTracing(interactive=False)
        self.CurrentTopology.construct_patches()
        self.CurrentTopology.GroupPatches()
        self.CheckPatches()

        if self.settings['patch_data']['preferences']['new_file']:
            self.SavePatches(self.settings['patch_data']['preferences']['new_fname'])

    def SavePatches(self, fname: str = '') -> None:
        """ Save patches as '.npy' file for later reconstruction in Ingrid.

        This file contains the information required to reconstruct patches
        at a later time and bypass the line_tracing.

        Args:
            fname (optional): Name of file/location for patch data.
        """
        if fname in ['', None]:
            fname = self.CurrentTopology.config + '_patches_' + str(int(time()))

        patch_data = [patch.as_np() for patch in self.CurrentTopology.patches.values()]
        patch_data.insert(0, self.CurrentTopology.config)
        patch_data.insert(1, self.CurrentTopology.LineTracer.NSEW_lookup)
        np.save(fname, np.array(patch_data))
        print(f"# Saved patch data for file {fname}.npy")

    def LoadPatches(self, fname: str = '') -> None:
        """ Load patches stored in an Ingrid generated '.npy' file.

        Args:
            fname (optional): Path to patch data. If no fname is provided to
                method 'LoadPatches', Ingrid code will check the settings
                'dict' for a file under entry

                'settings[''patch_data''][''file'']'

        """

        if type(fname) is not str:
            raise ValueError('# User did not provide a string to patch data.')

        if fname.strip() == '':  # Check if settings contains patch data.
            fname = self.settings['patch_data']['file']
        config, xpt_data, patches = self.ReconstructPatches(fname)
        self.process_yaml(Ingrid.ReadyamlFile(self.InputFile))
        self.StartSetup()
        self.LineTracer.NSEW_lookup = xpt_data
        self.SetTopology(config)
        self.CurrentTopology.patches = patches
        self.CurrentTopology.OrderPatches()
        self.CurrentTopology.SetupPatchMatrix()
        self.CheckPatches()

    def CreateSubgrid(self, NewFig: bool = True, ShowVertices: bool = False) -> None:
        """ Refine a generated patch map into a grid for exporting.

        Args:
            NewFig (optional): Plot the created grid on a new figure.

            ShowVertices (optional): Plot vertices in refined grid with
                bolded markers.

            ListPatches (optional): Generate a grid for a particular patch.
                Requires the correct patch name associated with the 'Patch'
                object.

                Warning: Grid generation is order dependent. Specifying a
                particular patch to generate a grid would only be done in
                rare cases and require the user to know dependencies for
                the particular patch.

        """

        if NewFig:
            try:
                plt.close(self._GridFig)
            except:
                pass
            self._GridFig = plt.figure(f'INGRID: {self.CurrentTopology.config} Grid', figsize=(6, 10))
            ax = self._GridFig.gca()
            ax.set_xlim([self.PsiUNorm.rmin, self.PsiUNorm.rmax])
            ax.set_ylim([self.PsiUNorm.zmin, self.PsiUNorm.zmax])
            ax.set_aspect('equal', adjustable='box')

            ax.set_xlabel('R')
            ax.set_ylabel('Z')
            ax.set_title(f'{self.CurrentTopology.config} Grid Diagram')
            ax.set_aspect('equal', adjustable='box')

        self.CurrentTopology.construct_grid(ShowVertices=ShowVertices, RestartScratch=RestartScratch,
                                            ListPatches=ListPatches)

    def ExportGridue(self, fname: str = 'gridue') -> None:
        """ Export a gridue file for the created grid.

        Args:
            fname (optional): Name of gridue file to save.

        """
        self.PrepGridue()
        if type(self.CurrentTopology) in [SNL]:
            if self.WriteGridueSNL(self.CurrentTopology.gridue_settings, fname):
                print(f"# Successfully saved gridue file as '{fname}'")
        elif type(self.CurrentTopology) in [SF15, SF45, SF75, SF105, SF135, SF165, UDN]:
            if self.WriteGridueDNL(self.CurrentTopology.gridue_settings, fname):
                print(f"# Successfully saved gridue file as '{fname}'")

    @staticmethod
    def ImportGridue(self, fname: str = 'gridue') -> dict:
        """Import UEDGE grid file as dictionary.

        Args:
            fname (optional): Path/file name to gridue formatted file.

        Returns:
            A dict containing header and body information from the gridue file.

        """
        try:
            f = open(fname, mode='r')
            Values = [int(x) for x in next(f).split()]
            HeaderItems = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iyseptrx1']
            gridue_settings = dict(zip(HeaderItems, Values))
            next(f)
            BodyItems = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
            Str = {i: [] for i in BodyItems}
            k = iter(Str.keys())
            Key = next(k)
            for line in f:
                if line == 'iogridue\n':
                    continue
                if line == '\n':
                    try:
                        Key = next(k)
                    except:
                        continue
                    print(Key)
                else:
                    Str[Key].append(line)
            f.close()
            nx = gridue_settings['nxm'] + 2
            ny = gridue_settings['nym'] + 2
            for k, v in Str.items():
                L = (''.join(v).replace('\n', '').replace('D', 'e')).split()
                _l = iter(L)
                vv = next(_l)

                data_ = np.zeros((nx, ny, 5))
                for n in range(5):
                    for j in range(ny):
                        for i in range(nx):

                            data_[i][j][n] = float(vv)

                            try:
                                vv = next(_l)
                            except:
                                continue
                gridue_settings[k] = data_
            return gridue_settings
        except Exception as e:
            print(repr(e))

    @staticmethod
    def Plotgridue(GridueParams: dict, edgecolor='black', ax: object = None):
        """Plot UEDGE grid from 'dict' obtained from method 'ImportGridue'

        Args:
            GridueParams: 'dict' of gridue header and body information.

            edgecolor (optional): Color of grid.

            ax (optional): Matplotlib axes to plot on.

        """
        r = GridueParams['rm']
        z = GridueParams['zm']
        Nx = len(r)
        Ny = len(r[0])
        patches = []
        if ax is None:
            ax = plt.gca()
        idx = [np.array([1, 2, 4, 3, 1])]
        for i in range(Nx):
            for j in range(Ny):
                p = matplotlib.patches.Polygon(np.concatenate((r[i][j][idx], z[i][j][idx])).reshape(2, 5).T, closed=True, fill=False, edgecolor=edgecolor)
                if not Fill:
                    ax.add_patch(p)  # draw the contours
                else:
                    patches.append(p)

        plt.ylim(z.min(), z.max())
        plt.xlim(r.min(), r.max())
        plt.show()

    @staticmethod
    def ReadyamlFile(FileName: str) -> dict:
        """Read a yaml file and return a dictionary

        Args:
            FileName: Path/file name of '.yml' parameter file to
                represent as a dictionary.

        Returns:
            Settings file represented as a dictionary

        Raises:
            IOError: If error occurs while loading yml file.

        """
        File = pathlib.Path(FileName).expanduser()
        if File.exists() and File.is_file():
            try:
                with open(FileName, 'r') as f:
                    ymlDict = yml.load(f, Loader=yml.Loader)
                return ymlDict
            except:
                raise IOError('Cannot load the yml file: ' + FileName)

        else:
            print('Current Working Directory:' + os.getcwd())
            raise IOError('Cannot find the file: ' + FileName)
