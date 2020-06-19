#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 08:58:43 2019

@author: garcia299
"""

from __future__ import print_function

from sys import platform as sys_pf
import matplotlib
try:
    matplotlib.use("TkAgg")
except:
    pass
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.pyplot import close
from matplotlib.figure import Figure
from scipy.optimize import root
from numpy import full_like
import yaml

try:
    from pathlib import Path
except:
    from pathlib2 import Path
try:
    import tkinter as tk
except:
    import Tkinter as tk
try:
    import tkinter.filedialog as tkFD
except:
    import tkFileDialog as tkFD
try:
    from tkinter import messagebox as tkMB
except:
    import tkMessageBox as tkMB

import Ingrid
import Topologies
import Root_Finder

import pdb

helv_large = 'Helvetica 13 bold'
helv_medium = 'Helvetica 9 bold'

class IngridApp(tk.Tk):
    """
    IngridApp:
    ----------
        - Primary controller/hub for the GUI application.
        - Contains the 'FilePicker' and 'ParamPicker' pages,
          along with our INGRID data corresponding to the
          current session.

    TODO:
    -----
        - Determine if helper functions can be established
          to streamline the other classes.
    """

    def __init__(self, master = None,IngridSession=None,*args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand = True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        self.gui_mode = 'SNL'
        self.update_gui_dimensions()

        # attached to the parent Ingrid instance instead of new instance
        if IngridSession is None:
            self.IngridSession = Ingrid.Ingrid(params = {'grid_params' : {'num_xpt' : 1}})
        else:
            self.IngridSession = IngridSession
        self.populate_GUI()

        #automatic loading of input files
        if self.IngridSession.InputFile is not None:
            self.frames[FilePicker].load_param_file(ExistingParamFile=self.IngridSession.yaml)
            self.frames[ParamPicker].load_files()




    def populate_GUI(self):
        self.frames = {}

        for F in (FilePicker, ParamPicker):
            frame = F(self.container, self)
            frame.IngridSession=self.IngridSession
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.FilePickerFrame=self.frames[FilePicker]
        self.ParamPickerFrame=self.frames[ParamPicker]
        self.set_global_lookup()
        self.show_frame(ParamPicker)
        self.IngridMenubar = MenuBarControl(self)

    def show_frame(self, item):
        frame = self.frames[item]
        frame.tkraise()
        self.geometry(self.gui_dimensions[item])

    def change_numxpt(self, config = 'SNL'):
        self.gui_mode = config
        self.update_gui_dimensions()
        return self.reset_data(message = 'Changing number of x-points requires a reset.\nContinue?')

    def set_global_lookup(self):
        # TODO collect data from entry. Why not doing it straight into yaml file?
        _ref = self.frames[ParamPicker]
        if self.gui_mode == 'SNL':
            self.global_lookup = {'rmagx' : _ref.MagFrame.R_EntryText, 'zmagx' : _ref.MagFrame.Z_EntryText,\
                      'rxpt'  : _ref.XptFrame1.R_EntryText, 'zxpt'  : _ref.XptFrame1.Z_EntryText,\
                      'psi_max_r' : _ref.PsiMaxFrame.R_EntryText, 'psi_max_z' : _ref.PsiMaxFrame.Z_EntryText,\
                      'psi_min_core_r' : _ref.PsiMinFrame.R_EntryText, 'psi_min_core_z' : _ref.PsiMinFrame.Z_EntryText,\
                      'psi_min_pf_r' : _ref.PsiPrivateFrame.R_EntryText, 'psi_min_pf_z' : _ref.PsiPrivateFrame.Z_EntryText,\
                      'psi_max' : _ref.PsiMaxFrame.Psi_EntryText, 'psi_min_core' : _ref.PsiMinFrame.Psi_EntryText,\
                      'psi_min_pf' : _ref.PsiPrivateFrame.Psi_EntryText

            }

        elif self.gui_mode == 'DNL':
            self.global_lookup = {'rmagx' : _ref.MagFrame.R_EntryText, 'zmagx' : _ref.MagFrame.Z_EntryText,\
                      'rxpt'  : _ref.XptFrame1.R_EntryText, 'zxpt'  : _ref.XptFrame1.Z_EntryText,\
                      'rxpt2' : _ref.XptFrame2.R_EntryText, 'zxpt2' : _ref.XptFrame2.Z_EntryText, \
                      'psi_max_r' : _ref.PsiMaxFrame.R_EntryText, 'psi_max_z' : _ref.PsiMaxFrame.Z_EntryText,\
                      'psi_min_core_r' : _ref.PsiMinFrame.R_EntryText, 'psi_min_core_z' : _ref.PsiMinFrame.Z_EntryText,\
                      'psi_min_pf_r' : _ref.PsiPrivateFrame.R_EntryText, 'psi_min_pf_z' : _ref.PsiPrivateFrame.Z_EntryText,\
                      'psi_max' : _ref.PsiMaxFrame.Psi_EntryText, 'psi_min_core' : _ref.PsiMinFrame.Psi_EntryText,\
                      'psi_min_pf' : _ref.PsiPrivateFrame.Psi_EntryText,\
                      'psi_max_r_outer' : _ref.PsiMaxOuterFrame.R_EntryText, 'psi_max_z_outer' : _ref.PsiMaxOuterFrame.Z_EntryText, \
                      'psi_max_outer' : _ref.PsiMaxOuterFrame.Psi_EntryText, \
                      'psi_max_r_inner' : _ref.PsiMaxInnerFrame.R_EntryText, 'psi_max_z_inner' : _ref.PsiMaxInnerFrame.Z_EntryText, \
                      'psi_max_inner' : _ref.PsiMaxInnerFrame.Psi_EntryText, \
                      'psi_pf2_r' : _ref.PsiPrivate2Frame.R_EntryText, 'psi_pf2_z' : _ref.PsiPrivate2Frame.Z_EntryText, \
                      'psi_pf2' : _ref.PsiPrivate2Frame.Psi_EntryText,
            }

    def update_gui_dimensions(self):
        self.gui_dimensions = {}
        if self.gui_mode == 'SNL':
            self.gui_dimensions = {FilePicker : "550x270", ParamPicker : "1300x475"}
        elif self.gui_mode == 'DNL':
            self.gui_dimensions = {FilePicker : "550x380", ParamPicker : "1155x475"}

    def reset_data(self, message = 'Are you sure you want to reset?'):
        if tkMB.askyesno('', message):
            self.frames = {}

            try:
                # Close all open figures
                close('all')
            except:
                pass

            self.populate_GUI()

            if self.gui_mode == 'SNL':
                nxpt = 1
            elif self.gui_mode == 'DNL':
                nxpt = 2
            self.IngridSession = Ingrid.Ingrid({'grid_params' : {'num_xpt' : nxpt}})

            self.show_frame(ParamPicker)
            self.geometry(self.gui_dimensions[ParamPicker])
            return True

        else:
            return False

    def exit_app(self):
        if tkMB.askyesno('', 'Are you sure you want to quit?'):
            try:
                plt.close('all')
            except:
                pass
            self.destroy()



class MenuBarControl(tk.Tk):

    def __init__(self, controller):
        self.controller = controller
        self.menubar = tk.Menu(controller)

        controller.config(menu = self.menubar)

        ingridTab = tk.Menu(self.menubar)
        ingridTab.add_command(label='Preferences...', command=self.open_preferences)
        ingridTab.add_separator()
        ingridTab.add_command(label='Select Data', command=controller.frames[ParamPicker].load_files)
        ingridTab.add_command(label='Load YAML Parameter File', command = self.load_yaml)
        ingridTab.add_command(label='Save Parameters', command=controller.frames[ParamPicker].saveParameters)
        ingridTab.add_separator()
        ingridTab.add_command(label='Restart INGRID session', command=controller.reset_data)
        ingridTab.add_command(label="Quit INGRID", command=controller.exit_app)
        patchTab = tk.Menu(self.menubar)
        patchTab.add_command(label='Generate Patches', command = controller.ParamPickerFrame.CreatePatches)
        gridTab = tk.Menu(self.menubar)
        gridTab.add_command(label = 'Generate Grid', command = controller.ParamPickerFrame.GenerateGrid)
        gridTab.add_command(label = 'Export gridue', command = controller.ParamPickerFrame.write_gridue)

        self.menubar.add_cascade(label = 'INGRID', menu = ingridTab)
        self.menubar.add_cascade(label = 'Patch Controls', menu = patchTab)
        self.menubar.add_cascade(label = 'Grid Controls', menu = gridTab)

        patchTab.entryconfig('Generate Patches', state = 'disabled')
        gridTab.entryconfig('Generate Grid', state = 'disabled')
        gridTab.entryconfig('Export gridue', state = 'disabled')

        self.ingridTab = ingridTab
        self.patchTab = patchTab
        self.gridTab = gridTab

        self.preferences_window = None



    def load_yaml(self):
        self.controller.frames[ParamPicker].load_files()
        self.controller.frames[FilePicker].load_param_file()

    def enable_patch_generation(self):
        self.patchTab.entryconfig('Generate Patches', state = 'normal')

    def enable_grid_generation(self):
        self.gridTab.entryconfig('Generate Grid', state = 'normal')

    def enable_gridue_export(self):
        self.gridTab.entryconfig('Export gridue', state = 'normal')

    def open_preferences(self):

        def apply_settings():
            lookup = {'1' : 'SNL', '2' : 'DNL'}
            if lookup[nxpt_spinbox.get()] != self.controller.gui_mode:
                if not self.controller.change_numxpt(lookup[nxpt_spinbox.get()]):
                    pass

        def close_settings():
            self.preferences_window.destroy()
            self.preferences_window = None

        if self.preferences_window:
            self.preferences_window.destroy()
            self.preferences_window = None

        self.preferences_window = tk.Toplevel(self.controller)

        entry_container = tk.Frame(self.preferences_window)
        entry_container.grid(row = 0, column = 0, columnspan = 2, padx = 10, pady = 10, sticky = 'nsew')

        nxpt_label = tk.Label(entry_container, text = 'Number of x-points:')
        nxpt_label.grid(row = 0, column = 0, padx = 10, pady = 5, sticky = 'nsew')

        nxpt_text = tk.StringVar()
        if self.controller.gui_mode == 'SNL':
            txt = '1'
        elif self.controller.gui_mode == 'DNL':
            txt = '2'
        nxpt_text.set(txt)
        nxpt_spinbox = tk.Spinbox(entry_container, from_ = 1, to = 2, width = 5, state = 'readonly', textvariable=nxpt_text.get())
        nxpt_spinbox.grid(row = 0, column = 1, padx = 10, pady = 5, sticky = 'nsew')


        button_container = tk.Frame(self.preferences_window)
        button_container.grid(row = 1, column = 0, columnspan = 2, padx = 10, pady = 10, sticky = 'nsew')
        confirm = tk.Button(button_container, text = 'Apply settings', command = apply_settings)
        confirm.grid(row = 1, column = 0, padx = 10, pady = 5, sticky = 'nsew')
        close = tk.Button(button_container, text = 'Close', command = close_settings)
        close.grid(row = 1, column = 1, padx = 10, pady = 5, sticky = 'nsew')

class FilePickerFrame(tk.Frame):
    """
    Helper Class that represents a widget in the FilePicker page.
    Each object represents a specific file to be loaded.

    Parameters:
    ----------

    parent        :    tk.Frame

    controller    :    tk.Frame

    params        :    dict

    Dictionary containing tk.Entry text, tk.Button text,
    and a command for the particular tk.Button.

        Dictionary Format:
        {
        'ButtonText' : str,
        'EntryText' : str,
        'ButtonCommand' : function
        }

    """
    def __init__(self, parent, controller, params):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        self.FP_ButtonText = tk.StringVar()
        self.FP_ButtonText.set(params['ButtonText'])

        self.FP_EntryText = tk.StringVar()
        self.FP_EntryText.set(params['EntryText'])

        self.FP_Entry = tk.Entry(self, text = self.FP_EntryText, width = 40, \
                                 disabledbackground = '#f8f8ff', state = 'normal')
        self.FP_Button = tk.Button(self, text = self.FP_ButtonText.get(), width = 20, \
                                   font = helv_medium, command = params['ButtonCommand'])

        self.FP_Entry.grid(row = 0, column = 0, columnspan = 2, padx = 10, pady = 10, sticky = 'NSEW')
        self.FP_Button.grid(row = 0, column = 2, padx = 10, pady = 10, sticky = 'EW')

        self.isLoaded = False
        self._handle = None

    def fileLoaded(self, Path):
        print('File loaded:',Path)
        self.FP_EntryText.set(str(Path))
        self.FP_Entry.delete(0,tk.END)
        self.FP_Entry.insert(0,str(Path))
        self._handle = Path
        self.FP_Button.config(fg = 'lime green')
        self.isLoaded = True

    def invalidPath(self, Path):
        message = "Path is invalid: '{}'".format(str(Path))
        self._handle = None
        self.FP_Entry.delete(0,tk.END)
        self.FP_Entry.insert(0,str(Path))
        self.FP_Button.config(fg = 'red')
        self.isLoaded = False

    def invalidType(self, Path):
        message = "File is not *.txt: '{}'".format(str(Path.name))
        self._handle = None
        self.FP_EntryText.set(message)
        self.FP_Button.config(fg = 'red')
        self.isLoaded = False

class FilePicker(tk.Frame):
    """
    FilePicker:
    -----------
        - Class representing the Frame that contains all file selection
          functionality.

    TODO:
    -----
        - Transition from tk.pack geometry manager to tk.grid
          geometry manager.
        - Create Classes representing objects currently populating
          the FilePicker Frame.
        - General code-cleanup...

    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.IngridSession=controller.IngridSession
        self.parent = parent

        self.preview_loaded = False

        Title = tk.Label(self, text = "Ingrid.", font=helv_large)
        Title.grid(row = 0, column = 0, columnspan = 3, \
            padx = 10, pady = 5, sticky = 'NSEW')
        Title.grid(row = 0, column = 0, padx = 10, pady = 10, sticky = 'EW')

        self.eqdskFrame = FilePickerFrame(self, controller, \
            {'ButtonText' : 'Select EQDSK File', \
             'EntryText'  : 'Select an EQDSK File', \
             'ButtonCommand' : self.load_eqdsk_file \
            })
        self.itpFrame = FilePickerFrame(self, controller, \
            {'ButtonText' : 'Select Inner Plate', \
             'EntryText'  : 'Select an inner strike plate configuration file.', \
             'ButtonCommand' : self.load_itp_file \
            })
        self.otpFrame = FilePickerFrame(self, controller, \
            {'ButtonText' : 'Select Outer Plate', \
             'EntryText'  : 'Select an outer strike plate configuration file.', \
             'ButtonCommand' : self.load_otp_file \
            })

        self.itp2Frame = FilePickerFrame(self, controller, \
            {'ButtonText' : 'Select Inner Plate 2', \
             'EntryText'  : 'Select inner strike plate 2 configuration file.', \
             'ButtonCommand' : self.load_itp2_file \
            })
        self.otp2Frame = FilePickerFrame(self, controller, \
            {'ButtonText' : 'Select Outer Plate 2', \
             'EntryText'  : 'Select outer strike plate 2 configuration file.', \
             'ButtonCommand' : self.load_otp2_file \
            })

        self.paramFrame = FilePickerFrame(self, controller,\
            {'ButtonText' : 'Load Parameter File', \
             'EntryText'  : 'Select a pre-existing YAML format file.', \
             'ButtonCommand' : self.load_param_file \
             })

        if self.controller.gui_mode == 'DNL':
            self.FP_Frames = [self.eqdskFrame, self.itpFrame, self.otpFrame,\
                              self.itp2Frame, self.otp2Frame]
        elif self.controller.gui_mode == 'SNL':
            self.FP_Frames = [self.eqdskFrame, self.itpFrame, self.otpFrame]

        self.FP_handles = [f._handle for f in self.FP_Frames]

        for i in range(len(self.FP_Frames)):
            self.FP_Frames[i].grid(row = i + 1, column = 0, \
                padx = 10, pady = 5, sticky = 'NSEW')

        self.ControlPanel = tk.Frame(self)
        self.ControlPanel.grid(row = len(self.FP_Frames) + 1, column = 0, \
                padx = 10, pady = 5, sticky = 'NSEW')
        self.previewButton = tk.Button(self.ControlPanel, text = 'Preview Loaded Data', \
            font = helv_medium, state = 'disabled', command = self.preview_data)
        self.confirmButton = tk.Button(self.ControlPanel, text = 'Confirm', \
            font = helv_medium, state = 'disabled', command = self.confirm_data)
        self.resetButton = tk.Button(self.ControlPanel, text = 'Reset', \
            font = helv_medium, command = self.controller.reset_data)

        self.previewButton.grid(row = 0, column = 0, columnspan = 1, padx = 10, pady = 10, \
            sticky = 'NSEW')
        self.confirmButton.grid(row = 0, column = 1, padx = 10, pady = 10, sticky = 'NSEW')
        self.resetButton.grid(row = 0, column = 2, padx = 10, pady = 10, sticky = 'NSEW')


    def get_file(self):
        """
        General method to allow user to select a file via GUI.

        Return vals:
        -----------
            fpath : PathLib object
                    PathLib2 instance if file-exists.
                    Empty string if no file selected.

            True/False  : bool
                    Return whether a file was selected or not.

        """
        try:
            fpath = Path(tkFD.askopenfilename(initialdir = '..', title = 'Select File'))
            if fpath.is_file():
                return fpath, True
            else:
                raise TypeError
        except TypeError:
            return '', False

    def load_eqdsk_file(self, showFileDialog = True, toLoad = None):
        """
        Loads EQDSK file.

        Parameters:
        -----------
            toLoad : string, optional
                    Path to an EQDSK file that we
                    will attempt to load.
        Return Vals:
        ------------
            True/False : bool
                    Return whether or not the file
                    was successfully loaded.

        Post-Call:
        ----------
            Updates eqdsk_loaded flag to True if
            an EQDSK file was loaded successfully.

            Updates eqdsk_loaded flag to False if
            wrong file-type selected.

        """
        if showFileDialog:
            eqdsk_file, valid_path = self.get_file()
            if valid_path:
                self.controller.IngridSession.yaml['eqdsk'] = str(eqdsk_file)
                self.eqdskFrame.fileLoaded(eqdsk_file)
                self.update_frame_state()
        elif not showFileDialog:
            eqdsk_file = Path(toLoad)
            valid_path = eqdsk_file.is_file()
            if valid_path:
                print('In load_eqdsk_file:',eqdsk_file)

                self.controller.IngridSession.yaml['eqdsk'] = str(eqdsk_file)
                #print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                self.eqdskFrame.fileLoaded(eqdsk_file)
                self.update_frame_state()
            elif not valid_path:
                self.eqdskFrame.invalidPath(eqdsk_file)


    def load_itp_file(self, showFileDialog = True, toLoad = None):
        """
        Loads *.txt file containing ITP geometry.

        Post-Call:
        ----------
            Updates itp_loaded flag to True if
            a *.txt file was loaded successfully.

        """
        if showFileDialog:
            itp_file, valid_path = self.get_file()
            if valid_path and itp_file.suffix == '.txt':
                self.controller.IngridSession
                self.controller.IngridSession.yaml['target_plates']['plate_E1']['file'] = str(itp_file)
                print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                self.itpFrame.fileLoaded(itp_file)
                self.update_frame_state()
            elif not valid_path:
                self.itpFrame.invalidPath(itp_file)
            else:
                self.itpFrame.invalidType(itp_file)
        elif not showFileDialog:
            itp_file = Path(toLoad)
            valid_path = itp_file.is_file()
            if valid_path and itp_file.suffix == '.txt':
                print('In load_itp:',itp_file)

                self.controller.IngridSession.yaml['target_plates']['plate_E1']['file'] = str(itp_file)
                #print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                self.itpFrame.fileLoaded(itp_file)
                self.update_frame_state()

    def load_itp2_file(self, showFileDialog = True, toLoad = None):
        """
        Loads *.txt file containing ITP geometry.

        Post-Call:
        ----------
            Updates itp_loaded flag to True if
            a *.txt file was loaded successfully.

        """
        if showFileDialog:
            itp_file, valid_path = self.get_file()
            if valid_path and itp_file.suffix == '.txt':
                self.controller.IngridSession
                self.controller.IngridSession.yaml['target_plates']['plate_E2']['file'] = str(itp_file)
                print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                self.itp2Frame.fileLoaded(itp_file)
                self.update_frame_state()
            elif not valid_path:
                self.itp2Frame.invalidPath(itp_file)
            else:
                self.itp2Frame.invalidType(itp_file)
        elif not showFileDialog:
            if toLoad is not None:
                itp_file = Path(toLoad)
                valid_path = itp_file.is_file()
                if valid_path and itp_file.suffix == '.txt':
                    self.controller.IngridSession.yaml['target_plates']['plate_E2']['file'] = str(itp_file)
                    #print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                    self.itp2Frame.fileLoaded(itp_file)
                    self.update_frame_state()

    def load_otp_file(self, showFileDialog = True, toLoad=None):
        """
        Loads *.txt file containing OTP geometry.

        Post-Call:
        ----------
            Updates otp_loaded flag to True if
            a *.txt file was loaded successfully.

        """
        if showFileDialog:
            otp_file, valid_path = self.get_file()
            if valid_path and otp_file.suffix == '.txt':
                print(self.controller.IngridSession.yaml)
                self.controller.IngridSession.yaml['target_plates']['plate_W1']['file'] = str(otp_file)
                self.otpFrame.fileLoaded(otp_file)
                self.update_frame_state()
            elif not valid_path:
                self.otpFrame.invalidPath(otp_file)
            else:
                self.otpFrame.invalidType(otp_file)
        elif not showFileDialog:
            otp_file = Path(toLoad)
            valid_path = otp_file.is_file()
            if valid_path and otp_file.suffix == '.txt':
                print('In load_otp:', otp_file)
                self.controller.IngridSession.yaml['target_plates']['plate_W1']['file'] = str(otp_file)
                #print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                self.otpFrame.fileLoaded(otp_file)
                self.update_frame_state()

    def load_otp2_file(self, showFileDialog = True, toLoad=None):
        """
        Loads *.txt file containing OTP geometry.

        Post-Call:
        ----------
            Updates otp_loaded flag to True if
            a *.txt file was loaded successfully.

        """
        if showFileDialog:
            otp_file, valid_path = self.get_file()
            if valid_path and otp_file.suffix == '.txt':
                print(self.controller.IngridSession.yaml)
                self.controller.IngridSession.yaml['target_plates']['plate_W2']['file'] = str(otp_file)
                self.otp2Frame.fileLoaded(otp_file)
                self.update_frame_state()
            elif not valid_path:
                self.otp2Frame.invalidPath(otp_file)
            else:
                self.otp2Frame.invalidType(otp_file)
        elif not showFileDialog:
            if toLoad is not None:
                otp_file = Path(toLoad)
                valid_path = otp_file.is_file()
                if valid_path and otp_file.suffix == '.txt':
                    self.controller.IngridSession.yaml['target_plates']['plate_W2']['file'] = str(otp_file)
                    #print(yaml.dump(self.controller.IngridSession.yaml, indent = 4))
                    self.otp2Frame.fileLoaded(otp_file)
                    self.update_frame_state()


    def load_param_file(self,ExistingParamFile=None):
        """
        Loads *.yaml file containing INGRID parameters.
        Said *.yaml contains session settings.

        Post-Call:
        ----------
            Updates param_loaded flag to True if
            a *.yaml file was loaded successfully.

        TODO:
        -----
            Check whether the *.yaml file actually
            contained all the required files. This
            can be done by inspecting the yaml object.

        """
        if ExistingParamFile is None:
            param_file, valid_path = self.get_file()
            suffix_list = ['.txt',  '.yml', '.yaml', '']
            if valid_path and param_file.suffix in suffix_list:
                _yaml = yaml.load(param_file.read_text())
                print(yaml.dump(_yaml, indent = 4))
                self.parse_yaml(_yaml)
                self.paramFrame.fileLoaded(param_file)
                self.update_frame_state()
                #print(self.controller.IngridSession.yaml)
            elif not valid_path:
                return False
        else:

            self.CheckFile(ExistingParamFile)
            self.paramFrame.fileLoaded('')
            self.update_frame_state()
            return True
    def parse_yaml(self, _yaml):
            lookup = { \
            'eqdsk' : self.load_eqdsk_file, \
            'target_plates' : {'plate_E1' : self.load_itp_file, 'plate_W1' : self.load_otp_file,
                               'plate_W2' : self.load_otp2_file, 'plate_E2' : self.load_itp2_file} \
            }
            for item in _yaml:
                print(item)
                if item == 'eqdsk':
                    lookup[item](showFileDialog = False, toLoad = _yaml[item])
                elif item == 'target_plates':
                    for sub_item in _yaml[item]:
                        lookup[item][sub_item](showFileDialog = False, toLoad = _yaml[item][sub_item]['file'])
            self.controller.IngridSession.process_yaml(_yaml)

    def CheckFile(self,yaml):
            lookup = { \
            'eqdsk' : self.load_eqdsk_file, \
            'target_plates' : {'plate_E1' : self.load_itp_file, 'plate_W1' : self.load_otp_file,
                               'plate_W2' : self.load_otp2_file, 'plate_E2' : self.load_itp2_file} \
            }
            for item in yaml:
                if item == 'eqdsk':
                    lookup[item](showFileDialog = False, toLoad = yaml[item])
                elif item == 'target_plates':
                    for sub_item in yaml[item]:
                        lookup[item][sub_item](showFileDialog = False, toLoad = yaml[item][sub_item]['file'])

    def update_frame_state(self):
        """
        Keeps track of which files are loaded via
        our file-flags.

        When all files are loaded successfully, the
        user can proceed to the next Frame (ParamPicker)

        """
        # Button state bool flag would be nice here.
        # We'd be able to call a 'flip_state' function instead of
        # manually setting the config.
        if self.eqdskFrame.isLoaded:
            self.previewButton.config(state = 'normal')
        if not self.eqdskFrame.isLoaded:
            self.previewButton.config(state = 'disabled')
        if self.files_ready():
            self.confirmButton.config(state = 'normal')

    def update_plate_handles(self):
        self.FP_handles = [f._handle for f in self.FP_Frames]

    def files_ready(self):
        """
        Helper that determines if all the required files
        have been loaded or not.

        Return Vals:
            True/False : bool
                        Return whether or not all files
                        have been loaded.

        """
        files_ready = True
        for item in self.FP_Frames:
            if item is self.paramFrame:
                pass
            elif not item.isLoaded:
                files_ready = False
        return files_ready

    def preview_data(self):

        # self.controller.IngridSession = Ingrid.Ingrid(params = self.controller.IngridSession.yaml)
        self.controller.IngridSession.OMFIT_read_psi()

        if self.itpFrame.isLoaded:
            self.controller.IngridSession.read_target_plates()
            self.update_plate_handles()
        if self.otpFrame.isLoaded:
            self.controller.IngridSession.read_target_plates()
            self.update_plate_handles()
        self.controller.IngridSession.calc_efit_derivs()
        self.controller.IngridSession.plot_efit_data()
        self.controller.IngridSession.plot_target_plates()


    def confirm_data(self):
        for file in self.FP_Frames:
            if not file._handle in self.FP_handles:
                self.preview_data()
        for item in self.controller.frames[ParamPicker].RF_Frames:
            item.Edit_Button.config(state = 'normal')

        self.controller.IngridSession.find_roots(tk_controller = self.controller)
        try:
            self.controller.IngridSession.AutoRefineMagAxis()
            self.controller.IngridSession.AutoRefineXpoint()
        except: pass

        self.controller.frames[ParamPicker].load_frame_entries()
        self.controller.frames[ParamPicker].acceptRF_Button.config(state = 'normal')
        self.controller.show_frame(ParamPicker)
        self.controller.geometry(self.controller.gui_dimensions[ParamPicker])
        self.controller.IngridSession.root_finder.disconnect()

class ParamPicker(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.IngridSession = controller.IngridSession
        print(controller)
        print(parent)
        # TODO: Create a ParameterFrame Class. This will be useful when we have
        #       multiple X-points to deal with.
        #define and pack container

        self.CreateContainers()
        self.CreateButtons()
        self.SetupChildrenFrames()







        self.ActiveFrame = []
        self.curr_click = ()

        self.efit_loaded = False
        self.psinorm_loaded = False

    def SetupChildrenFrames(self):
        self.MagFrame = RootFinderFrame(self.RF_Container, self, title = 'Magnetic-Axis', style = helv_medium)
        self.XptFrame1 = RootFinderFrame(self.RF_Container, self, title = 'Primary X-Point', style = helv_medium)
        self.XptFrame2 = RootFinderFrame(self.RF_Container, self, title = 'Secondary X-Point', style = helv_medium)

        self.PsiMaxFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Max', style = helv_medium)
        self.PsiMinFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Core', style = helv_medium)
        self.PsiPrivateFrame = PsiFinderFrame(self.PF_Container, self, title = 'Primary PF', style = helv_medium)

        self.PsiMaxInnerFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Max Inner', style = helv_medium)
        self.PsiMaxOuterFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Max Outer', style = helv_medium)
        self.PsiPrivate2Frame = PsiFinderFrame(self.PF_Container, self, title = 'Secondary PF', style = helv_medium)
        self.GridParameterFrame = GridParameterFrame(self.GridParam_Container, self, title = 'Primary SOL', style = helv_medium)
        if self.controller.gui_mode == 'DNL':
            self.RF_Frames = [self.MagFrame, self.XptFrame1, self.XptFrame2]
            self.PF_Frames = [self.PsiMaxInnerFrame, self.PsiMaxOuterFrame, self.PsiMinFrame, self.PsiPrivateFrame, self.PsiPrivate2Frame]
        elif self.controller.gui_mode == 'SNL':
            self.RF_Frames = [self.MagFrame, self.XptFrame1]
            self.PF_Frames = [self.PsiMaxFrame, self.PsiMinFrame, self.PsiPrivateFrame]
        self.GridParams_Frames=[self.GridParameterFrame]
        self.All_Frames = []

        for F in self.RF_Frames:
            self.All_Frames.append(F)
        for F in self.PF_Frames:
            self.All_Frames.append(F)
        #for F in self.GridParams_Frames: # no need for that
            #self.All_Frames.append(F)

        for i, Frame in enumerate(self.RF_Frames):
            Frame.grid(row = i, column = 0, padx = 5, pady = 10, sticky = 'nsew')

        i = 0
        j = 0
        for Frame in self.PF_Frames:
            j += 1 if i == 0 else 0
            Frame.grid(row =  i, column = j, padx = 5, pady = 10, sticky = 'nsew')
            i += 1 if i < 2 else -i

        i = 0
        j = 0
        for Frame in self.GridParams_Frames:
            j += 1 if i == 0 else 0
            Frame.grid(row =  i, column = j, padx = 5, pady = 10, sticky = 'nsew')
            i += 1 if i < 2 else -i

    def CreateContainers(self):
        self.RF_Container = tk.LabelFrame(self, text = '1.', font = helv_medium, relief = 'groove')
        self.PF_Container = tk.LabelFrame(self, text = '2.', font = helv_medium)
        self.GridParam_Container = tk.LabelFrame(self, text = '.3', font = helv_medium)
        self.Controls_Container = tk.Frame(self)
        self.Action_Container = tk.LabelFrame(self.Controls_Container, text = 'INGRID', font = helv_medium)
        self.Settings_Container = tk.LabelFrame(self.Controls_Container, text = 'Settings', font = helv_medium)
        self.RF_Container.grid(row = 0, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PF_Container.grid(row = 0, column = 1, padx = 10, pady = 10, sticky = 'nsew')
        self.GridParam_Container.grid(row = 0, column = 2, padx = 10, pady = 10, sticky = 'nsew')
        self.GridButtons_Container=tk.LabelFrame(self)
        self.GridButtons_Container.grid(row = 2, column = 2, padx = 10, pady = 2, sticky = 'nsew')

    def CreateButtons(self):
        self.GenerateGrid_Button = tk.Button(self.GridButtons_Container, text = 'Generate Grid', font = helv_medium, \
                                         state = 'disabled', command = self.GenerateGrid)
        self.CreatePatches_Button = tk.Button(self.GridButtons_Container, text = 'Create Patches', font = helv_medium, \
                                         state = 'disabled', command = self.CreatePatches)
        self.CreatePatches_Button.grid(row = 0, column = 0, padx = 0, pady = 0, sticky = 'n')
        self.GenerateGrid_Button.grid(row = 0, column = 1, padx = 0, pady = 0, sticky = 'n')



        self.save_Button = tk.Button(self.Settings_Container, text = 'Save Parameters', \
                                     font = helv_medium, state = 'disabled', width = 30, \
                                     command = self.saveParameters)
        self.load_Button = tk.Button(self.Settings_Container, text = 'Select Files', font = helv_medium, \
                                     width = 30, command = self.load_files)
        self.settings_Button = tk.Button(self.Action_Container, text = 'Grid/Integrator Params', font = helv_medium,\
                                    width = 30, command = self.settings_window, state = 'disabled')
        self.export_Button = tk.Button(self.Action_Container, text = 'Export GRIDUE', font = helv_medium, \
                                     width = 30, state = 'disabled', command = self.write_gridue)
        self.acceptRF_Button = tk.Button(self, text = 'Confirm geometry', font = helv_medium, \
                                         state = 'disabled', command = self.set_RFValues)
        self.acceptRF_Button.grid(row = 2, column = 0, columnspan = 1, padx = 0, pady = 0)
        self.acceptPF_Button = tk.Button(self, text = 'Confirm psi', font = helv_medium, \
                                         state = 'disabled', command = self.set_PsiValues)
        self.acceptPF_Button.grid(row = 2, column = 1, columnspan = 1, padx = 0, pady = 0)

    def EnableGridGeneration(self):
        self.GridParameterFrame.enable_frame()
        self.GenerateGrid_Button.config(state = 'normal')


    def GenerateGrid(self):
        if self.GridParameterFrame.GetValues():
            if self.IngridSession.CreateSubgrid():
                self.controller.IngridMenubar.enable_gridue_export()



    def refine(self):
        sol = root(self.controller.IngridSession.root_finder.func, \
                [float(self.ActiveFrame[0].R_EntryText.get()), \
                 float(self.ActiveFrame[0].Z_EntryText.get())] )
        self.ActiveFrame[0].R_EntryText.set('{0:.12f}'.format(sol.x[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.12f}'.format(sol.x[1]))


    def activate_frame(self, calling_frame):
        if self.ActiveFrame:
            #self.ActiveFrame[0].disable_frame()
            self.ActiveFrame[0].toggle_mouse()
            self.ActiveFrame[0].Edit_Button.config(state = 'normal')
            self.ActiveFrame.pop()
        calling_frame.enable_frame()
        self.ActiveFrame.append(calling_frame)
        self.ActiveFrame[0].toggle_mouse()
        if calling_frame in self.RF_Frames:
            self.efit_loaded = True
        if calling_frame in self.PF_Frames:
            self.psinorm_loaded = True

    def unlock_PF_Frame(self):
        try:
            close(self.controller.IngridSession.efit_psi.name)
        except:
            pass

        for F in self.PF_Frames:
            F.Edit_Button.config(state = 'normal')

        self.controller.IngridSession.calc_psinorm()
        self.controller.IngridSession.plot_psinorm()
        self.controller.IngridSession.PlotPsiNormBounds()
        self.controller.IngridSession.find_psi_lines(tk_controller = self.controller)
        self.controller.IngridSession.psi_finder.disconnect()
        self.controller.IngridSession.init_LineTracing()
        self.controller.frames[ParamPicker].acceptPF_Button.config(state = 'normal')

    def unlock_controls(self):
        self.settings_Button.config(state = 'normal')
        self.CreatePatches_Button.config(state = 'normal')
        self.save_Button.config(state = 'normal')

    def load_files(self):
        self.controller.frames[FilePicker].preview_loaded = False
        self.controller.show_frame(FilePicker)
        self.controller.geometry(self.controller.gui_dimensions[FilePicker])

    def update_root_finder(self):
        self.ActiveFrame[0].R_EntryText.set('{0:.12f}'.format(self.curr_click[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.12f}'.format(self.curr_click[1]))
        if self.ActiveFrame[0] in self.PF_Frames:
            psi = self.controller.IngridSession.psi_norm.get_psi(self.curr_click[0],self.curr_click[1])
            self.ActiveFrame[0].Psi_EntryText.set('{0:.12f}'.format(psi))

    def set_RFValues(self):

        lookup = self.controller.global_lookup
        for key in lookup.keys():
            self.controller.IngridSession.yaml['grid_params'][key] = float(lookup[key].get())

        self.MagAxis = (self.controller.IngridSession.yaml['grid_params']['rmagx'], self.controller.IngridSession.yaml['grid_params']['zmagx'])
        self.controller.IngridSession.magx = self.MagAxis
        self.Xpt = (self.controller.IngridSession.yaml['grid_params']['rxpt'], self.controller.IngridSession.yaml['grid_params']['zxpt'])
        self.controller.IngridSession.xpt1 = self.Xpt

        if self.controller.gui_mode == 'DNL':
            self.Xpt2 = (self.controller.IngridSession.yaml['grid_params']['rxpt2'], self.controller.IngridSession.yaml['grid_params']['zxpt2'])
            self.controller.IngridSession.xpt2 = self.Xpt2

        # TODO: Exception handling behind the scenes for ensuring PF_Frame is indeed ready.

        self.acceptRF_Button.config(text = 'Entries Saved!', fg = 'lime green')
        # plt.close('all')
        self.controller.IngridSession.root_finder.disconnect()
        for F in self.RF_Frames:
            if F.active_frame:
                F.update_frame()
                F.config(text = F.title, fg = 'black')
                F.Edit_Button.config(text = 'Edit ' + F.title, fg = 'black')
        self.UnlockPsiValues()

    def set_PsiValues(self):
        #Add toggle options
        if not self.IsRFLocked:
            self.LockPsiValues()
        else:
            self.UnlockPsiValues()

    def UnlockPsiValues(self):
        self.IsRFLocked=False
        self.unlock_PF_Frame()
        self.acceptPF_Button.config(state = 'normal')

    def LockPsiValues(self):
        self.IsRFLocked=True
        for F in self.PF_Frames:
            if F.R_EntryText.get() == '':
                F.R_EntryText.set('0.0')
            if F.Z_EntryText.get() == '':
                F.Z_EntryText.set('0.0')

        lookup = self.controller.global_lookup

        for key in lookup.keys():
            self.controller.IngridSession.yaml['grid_params'][key] = float(lookup[key].get())

        self.unlock_controls()
        self.controller.IngridSession.psi_finder.disconnect()
        for F in self.PF_Frames:
            if F.active_frame:
                F.update_frame()
                F.config(text = F.title, fg = 'black')
                F.Edit_Button.config(text = 'Edit ' + F.title, fg = 'black')
        self.acceptPF_Button.config(text = 'Entries Saved!', fg = 'lime green')
        self.controller.IngridSession.grid_params = self.controller.IngridSession.yaml['grid_params']
        self.controller.IngridSession.integrator_params = self.controller.IngridSession.yaml['integrator_params']
        self.controller.IngridSession.target_plates = self.controller.IngridSession.yaml['target_plates']

        self.controller.IngridMenubar.enable_patch_generation()

    def set_INGRID_params(self):

        self.set_RFValues()
        self.set_PsiValues()
        self.GridParameterFrame.GetValues()
        # TODO: why not keep the parameters in the yaml dic?
        self.controller.IngridSession.grid_params = self.controller.IngridSession.yaml['grid_params']
        self.controller.IngridSession.integrator_params = self.controller.IngridSession.yaml['integrator_params']
        self.controller.IngridSession.target_plates = self.controller.IngridSession.yaml['target_plates']

    def refresh_yaml(self):
        if self.controller.frames[FilePicker].paramFrame.isLoaded:
            f = Path(self.controller.frames[FilePicker].paramFrame.FP_EntryText.get())
            self.controller.IngridSession.process_yaml(yaml.load(f.read_text()))

    def CreatePatches(self):

        #self.refresh_yaml()
        self.set_INGRID_params()

        result = self.check_psi_validity()

        if not result['psi_valid']:
            tkMB.showerror("Psi Value Error", result['message'])
            return

        self.controller.IngridSession._analyze_topology()
        self.controller.IngridSession.init_LineTracing(refresh = True)

        self.Grid = self.controller.IngridSession.current_topology
        self.Grid.yaml = self.controller.IngridSession.yaml


        if isinstance(self.Grid, Topologies.SNL):
            self.Grid.magx = self.MagAxis
            self.Grid.xpt1 = self.Xpt
        elif isinstance(self.Grid, Topologies.DNL):
            self.Grid.magx = self.MagAxis
            self.Grid.xpt1 = self.Xpt
            self.Grid.xpt2 = self.Xpt2

        for F in self.PF_Frames:
            if F.active_mouse.get() == False:
                F.toggle_mouse()
            F.disable_frame()

        try:
            self.Grid.construct_patches()
            self.Grid.patch_diagram()
            self.controller.IngridMenubar.enable_grid_generation()
            self.EnableGridGeneration()
            self.CreatePatches_Button.config(fg='lime green')
        except Exception as e:
            print('=' * 80)
            print('ERROR DURING PATCH CONSTRUCTION: {}'.format(repr(e)))
            print('=' * 80 + '\n')

    def check_psi_validity(self):
        result = {'psi_valid' : True, 'message' : 'Success'}
        return result

        grid = self.controller.IngridSession

        (r,z) = grid.plate_data['plate_E1']['psi_max_rz']
        E1_max = grid.psi_norm.get_psi(r,z)
        print('E1_max at {} with psi value of {}'.format((r,z), E1_max))
        (r,z) = grid.plate_data['plate_E1']['psi_min_rz']
        E1_min = grid.psi_norm.get_psi(r,z)
        print('E1_min at {} with psi value of {}'.format((r,z), E1_min))

        (r,z) = grid.plate_data['plate_W1']['psi_max_rz']
        W1_max = grid.psi_norm.get_psi(r,z)
        print('W1_max at {} with psi value of {}'.format((r,z), W1_max))
        (r,z) = grid.plate_data['plate_W1']['psi_min_rz']
        W1_min = grid.psi_norm.get_psi(r,z)
        print('W1_min at {} with psi value of {}'.format((r,z), W1_min))

        psi_list = ['psi_max', 'psi_min_pf'] if self.controller.gui_mode == 'SNL' else ['psi_min_pf']
        psi_error = {'psi_max' : False, 'psi_min_pf' : False}
        E1_error = False
        W1_error = False

        error_message = 'The following do not intersect a target plate: \n'
        for v in psi_list:
            if grid.yaml['grid_params'][v] > E1_max or \
                grid.yaml['grid_params'][v] < E1_min:
                psi_error[v] = True
                E1_error = True
                error_message += '"' + v + '"' + ' at plate_E1' + '\n'
            if grid.yaml['grid_params'][v] > W1_max or \
                grid.yaml['grid_params'][v] < W1_min:
                psi_error[v] = True
                W1_error = True
                error_message += '"' + v + '"' + ' at plate_W1' + '\n'
        if E1_error or W1_error:
            result['message'] = error_message
            result['psi_valid'] = False

        return result


    def write_gridue(self):
        fname = tkFD.asksaveasfilename(initialdir = '.', title = 'Save File', defaultextension ='', initialfile = 'gridue')
        if fname != '':
            self.controller.IngridSession.export(fname)
            print('Saved gridue file as "{}"'.format(fname))
        else:
            print('Cancelling export...')

    def saveParameters(self):
        self.set_INGRID_params()
        fname = Path(tkFD.asksaveasfilename(initialdir = '.', title = 'Save File', defaultextension ='.yml'))
        try:
            documentation = Path('../Docs/INGRID_YAML_CONTROLS.txt')
            fname.write_text(documentation.resolve().read_text() + '\n' + yaml.dump(self.controller.IngridSession.yaml, indent = 4))
        except:
            fname.write_text(yaml.dump(self.controller.IngridSession.yaml, indent = 4))  # force tag is to overwrite the previous file
        print("Saved parameters to '{}'.".format(fname))
    #TODO: why adding a layer on top of the yaml dic?
    def load_frame_entries(self):

        ignore = ['config', 'num_xpt', 'np_global', 'nr_global']

        for param in self.controller.IngridSession.yaml['grid_params']:
            if param in ignore:
                continue
            try:
                self.controller.global_lookup[param].set('{0:.12f}'.format(self.controller.IngridSession.yaml['grid_params'][param]))
                print('Set "{}" to {}'.format(param, self.controller.IngridSession.yaml['grid_params'][param]))
            except:
                print('Did not find: ' + str(param))

    def settings_window(self):

        def confirm_settings():
            self.preferences_window.withdraw()

        confirm = tk.Button(container4, text = 'Confirm', command = confirm_settings)
        confirm.grid(row = i + 1, column = 0, padx = 10, pady = 5, sticky = 'nsew')

        close = tk.Button(container4, text = 'Close')
        close.grid(row = i + 1, column = 1, padx = 10, pady = 5, sticky = 'nsew')


class SettingsFrame(tk.Frame):
    def __init__(self, parent, controller, params):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        self.Save_ButtonText = tk.StringVar()
        self.Save_ButtonText.set('Save')

        self.Edit_ButtonText = tk.StringVar()
        self.Edit_ButtonText.set('Edit')

        self.Setting_LabelText = tk.StringVar()
        self.Setting_LabelText.set(params['LabelText'])
        self.Setting_EntryText = tk.StringVar()
        self.Setting_EntryText.set(params['EntryText'])

        self.Setting_Label = tk.Label(self, text = self.Setting_LabelText.get())

        self.Setting_Entry = tk.Entry(self, text = self.Setting_EntryText, width = 40, \
                                 disabledbackground = '#f8f8ff', state = 'disabled')
        self.Save_Button = tk.Button(self, text = self.Save_ButtonText.get(), width = 20, \
                                   font = helv_medium)
        self.Edit_Button = tk.Button(self, text = self.Edit_ButtonText.get(), width = 20, \
                                   font = helv_medium)

        self.Setting_Label.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = 'NSEW')
        self.Setting_Entry.grid(row = 0, column = 1, columnspan = 2, padx = 10, pady = 10, sticky = 'NSEW')
        self.Save_Button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = 'EW')
        self.Edit_Button.grid(row = 0, column = 4, padx = 10, pady = 10, sticky = 'EW')

        self.isLoaded = False

class RootFinderFrame(tk.LabelFrame):

    def __init__(self, parent, controller, title, style):
        tk.LabelFrame.__init__(self, parent, text = title, font = style, relief = 'flat')

        rbWrapper = tk.LabelFrame(self, text = 'Mouse Entry:', font = helv_medium)
        R_Label = tk.Label(self, text = 'R:', font = helv_medium)
        Z_Label = tk.Label(self, text = 'Z:', font = helv_medium)
        self.title = title
        self.controller = controller
        self.R_EntryText = tk.StringVar()
        self.Z_EntryText = tk.StringVar()
        self.All_EntryText = (self.R_EntryText, self.Z_EntryText)

        self.R_Entry = tk.Entry(self, text = self.R_EntryText, disabledbackground = '#eee9e9')
        self.Z_Entry = tk.Entry(self, text = self.Z_EntryText, disabledbackground = '#eee9e9')
        self.All_Entries = [self.R_Entry, self.Z_Entry]

        self.Refine_Button = tk.Button(self, text = 'Refine ' + self.title, font = helv_medium, \
                                       command = controller.refine, width = 20)
        self.Edit_Button = tk.Button(self, text = 'Edit ' + self.title, font = helv_medium, \
                command = lambda : self.toggle_edit(parent, controller), width = 20)
        self.All_Buttons = (self.Refine_Button, self.Edit_Button)

        self.active_mouse = tk.BooleanVar()
        self.active_mouse.set(0)

        self.active_frame = False

        self.On_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Enabled', font = helv_medium, \
                                            variable = self.active_mouse, value = 1, \
                                            command = self.toggle_mouse)
        self.Off_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Disabled', font = helv_medium, \
                                            variable = self.active_mouse, value = 0, \
                                            command = self.toggle_mouse)
        self.All_Radiobuttons = (self.On_Radiobutton, self.Off_Radiobutton)

        R_Label.grid(row = 0, column = 0, sticky = 'w', padx = 2, pady = 2)
        Z_Label.grid(row = 1, column = 0, sticky = 'w', padx = 2, pady = 2)

        self.R_Entry.grid(row = 0, column = 1, sticky = 'NSew', padx = 2, pady = 2)
        self.Z_Entry.grid(row = 1, column = 1, sticky = 'NSew', padx = 2, pady = 2)

        self.Refine_Button.grid(row = 2, column = 1, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        self.Edit_Button.grid(row = 2, column = 2, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        rbWrapper.grid(row = 0, rowspan = 2, column = 2, columnspan = 1, sticky = 'ns', padx = 10, pady = 5)
        self.On_Radiobutton.grid(row = 1, column = 2, sticky = 'we', padx = 2, pady = 2)
        self.Off_Radiobutton.grid(row = 2, column = 2, sticky = 'we', padx = 2, pady = 2)

        self.disable_frame()

    def toggle_edit(self, parent, controller):
        self.config(text = 'Editing ' + self.title, fg = 'blue')
        self.Edit_Button.config(text = 'Editing ' + self.title, fg = 'blue')
        controller.activate_frame(self)
        for F in controller.All_Frames:
            if F != self:

                F.config(text = F.title, fg = 'black')
                F.Edit_Button.config(text = 'Edit ' + F.title, fg = 'black')

    def toggle_mouse(self):
        if self.active_mouse.get() is True:
            self.controller.controller.IngridSession.root_finder.connect()
            print('Connected')
        else:
            self.controller.controller.IngridSession.root_finder.disconnect()
            print('Disconnected')

    def update_frame(self):
        if self.active_frame == False:
            self.enable_frame()
        else:
            self.disable_frame()

    def disable_frame(self):
        self.active_frame = False
        self.active_mouse.set(0)
        self.Refine_Button.config(state = 'disabled')
        self.Edit_Button.config(state = 'disabled')
        for entry in self.All_Entries:
            entry.config(state = 'disabled')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'disabled')

    def enable_frame(self):
        self.active_frame = True
        self.active_mouse.set(1)
        self.Refine_Button.config(state = 'normal')
        for entry in self.All_Entries:
            entry.config(state = 'normal')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'normal')

class PsiFinderFrame(tk.LabelFrame):

    def __init__(self, parent, controller, title, style):
        self.title = title
        self.controller = controller
        tk.LabelFrame.__init__(self, parent, text = title, font = style, relief = 'flat')
        R_Label = tk.Label(self, text = 'R:', font = helv_medium)
        Z_Label = tk.Label(self, text = 'Z:', font = helv_medium)
        Psi_Label = tk.Label(self, text = 'Psi:', font = helv_medium)

        rbWrapper = tk.LabelFrame(self, text = 'Mouse Entry:', font = helv_medium)
        self.R_EntryText = tk.StringVar()
        self.Z_EntryText = tk.StringVar()
        self.Psi_EntryText = tk.StringVar()
        self.All_EntryText = [self.R_EntryText, self.Z_EntryText, self.Psi_EntryText]

        self.R_Entry = tk.Entry(self, text = self.R_EntryText, disabledbackground = '#eee9e9')
        self.Z_Entry = tk.Entry(self, text = self.Z_EntryText, disabledbackground = '#eee9e9')
        self.Psi_Entry = tk.Entry(self, text = self.Psi_EntryText, disabledbackground = '#eee9e9')
        self.All_Entries = [self.R_Entry, self.Z_Entry, self.Psi_Entry]

        self.Edit_Button = tk.Button(self, text = 'Edit ' + self.title, font = helv_medium, \
                command = self.toggle_edit, width = 20)

        self.All_Buttons = [self.Edit_Button]

        self.active_mouse = tk.BooleanVar()
        self.active_mouse.set(0)
        self.active_frame = False

        self.On_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Enabled', font = helv_medium, \
                                            variable = self.active_mouse, value = 1, \
                                            command = self.toggle_mouse)
        self.Off_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Disabled', font = helv_medium, \
                                            variable = self.active_mouse, value = 0, \
                                            command = self.toggle_mouse)
        self.All_Radiobuttons = (self.On_Radiobutton, self.Off_Radiobutton)

        R_Label.grid(row = 0, column = 0, sticky = 'w', padx = 2, pady = 2)
        Z_Label.grid(row = 1, column = 0, sticky = 'w', padx = 2, pady = 2)
        Psi_Label.grid(row = 2, column = 0, sticky = 'w', padx = 2, pady = 2)

        self.R_Entry.grid(row = 0, column = 1, sticky = 'nsew', padx = 2, pady = 2)
        self.Z_Entry.grid(row = 1, column = 1, sticky = 'nsew', padx = 2, pady = 2)
        self.Psi_Entry.grid(row = 2, column = 1, sticky = 'nswe', padx = 2, pady = 2)
        self.Psi_Entry.bind('<Return>', self.plot_level)

        self.Edit_Button.grid(row = 2, column = 2, columnspan = 1, sticky = 'we', padx = 2, pady = 2)


        rbWrapper.grid(row = 0, rowspan = 2, column = 2, columnspan = 1, sticky = 'ns', padx = 10, pady = 5)
        self.On_Radiobutton.grid(row = 1, column = 2, sticky = 'we', padx = 2, pady = 2)
        self.Off_Radiobutton.grid(row = 2, column = 2, sticky = 'we', padx = 2, pady = 2)

        self.disable_frame()

    def toggle_mouse(self):
        if self.active_mouse.get() is True:
            self.controller.controller.IngridSession.psi_finder.connect()
            print('psi_finder connected')
        else:
            self.controller.controller.IngridSession.psi_finder.disconnect()
            print('psi_finder disconnected')

    def toggle_edit(self):
        pass
        self.config(text = 'Editing ' + self.title, fg = 'blue')
        self.Edit_Button.config(text = 'Editing ' + self.title, fg = 'blue')
        self.controller.activate_frame(self)
        for F in self.controller.All_Frames:
            if F != self:
                F.config(text = F.title, fg = 'black')
                F.Edit_Button.config(text = 'Edit ' + F.title, fg = 'black')

    def disable_frame(self):
        self.active_frame = False
        self.active_mouse.set(0)
        self.Edit_Button.config(state = 'disabled')
        for entry in self.All_Entries:
            entry.config(state = 'disabled')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'disabled')

    def update_frame(self):
        if self.active_frame == False:
            self.enable_frame()
        else:
            self.disable_frame()

    def enable_frame(self):
        self.active_frame = True
        self.active_mouse.set(1)
        for entry in self.All_Entries:
            entry.config(state = 'normal')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'normal')

    def plot_level(self, event = None):
        self.controller.controller.IngridSession.psi_finder._level()


class GridParameterFrame(tk.LabelFrame):

    def __init__(self, parent, controller, title, style,Verbose=False):
        Verbose=True
        self.controller = controller
        self.IngridSession = controller.IngridSession
        self.Verbose=Verbose
        self.title = title
        self.controller=controller
        tk.LabelFrame.__init__(self, parent, text = title, font = style, relief = 'flat')




        #set label and entry for Npol/Nread for each region (#TODO: must add handling of DNL)

        #This is the key setting: set label/emtry for function ( #TODO:must add handling of DNL)
        # You only need to modify this line to add whatever entry and field
        #TODO Add a label entry to have nicer description of each item.

        self.ListNpts={}
        KeysNpts=['grid_params','grid_generation']

        self.ListNpts['SOLNpol']=('Npol SOL',KeysNpts+['np_sol'])
        self.ListNpts['SOLNrad']=('Nrad SOL',KeysNpts+['nr_sol'])
        self.ListNpts['CoreNpol']=('Npol Core',KeysNpts+['np_core'])
        self.ListNpts['CoreNrad']=('Nrad  Core',KeysNpts+['nr_core'])
        self.ListNpts['PFRNpol']=('Npol PFR',KeysNpts+['np_pf'])
        self.ListNpts['PFRNrad']=('Nrad PFR',KeysNpts+['nr_pf'])

        self.ListFunc={}
        KeysFunc=['grid_params','grid_generation']
        KeysTarget=['target_plates']
        self.ListFunc['rad_f_pf']=('f(x) rad PFR',KeysFunc+['radial_f_primary_pf'])
        self.ListFunc['rad_f_sol']=('f(x) rad SOL',KeysFunc+['radial_f_primary_sol'])
        self.ListFunc['rad_f_core']=('f(x) rad SOL',KeysFunc+['radial_f_core'])
        self.ListFunc['pol_f_W1']=('f(x) pol W1',KeysTarget+['plate_W1','poloidal_f'])
        self.ListFunc['pol_f_E1']=('f(x) pol E1',KeysTarget+['plate_E1','poloidal_f'])

        self.Functions = {'x,x','x,x**2','x,exp(x)','x,x*exp(x)','x,x**3'}
        #Create tk objects
        self.All_Entries=[]
        Nelem=0
        Nelem=self.CreateNptsEntry(self.ListNpts,Nelem)
        Nelem=self.CreateFuncEntry(self.ListFunc,Nelem)
        self.SetValues()



        self.active_mouse = tk.BooleanVar()
        self.active_mouse.set(0)

    def CreateFuncEntry(self,List,Nelem):
        count=0
        for Name,Properties in List.items():
            exec("self.{}_Label = tk.Label(self, text = '{}:', font = helv_medium)".format(Name,Properties[0]),globals(),locals())
            exec("self.{}_Label.grid(row = {}, column = 0, sticky = 'w', padx = 2, pady = 2)".format(Name,count+Nelem),globals(),locals())
            exec("self.{}= tk.StringVar()".format(Name),globals(),locals())
            com="self.{}_Entry = tk.Entry(self, text = self.{}, disabledbackground = '#eee9e9')".format(Name,Name)
            exec(com,globals(),locals())
            com="self.All_Entries.append(self.{}_Entry)".format(Name)
            exec(com,globals(),locals())
            com="self.{}.trace_add('write', lambda *_,self=self:self.CheckFunc('{}'))".format(Name,Name)
            exec(com,globals(),locals())
            exec("self.{}_Entry.grid(row = {}, column = 1, sticky = 'nsew', padx = 2, pady = 2)".format(Name,count+Nelem),globals(),locals())
            exec("self.{}_MenuVar= tk.StringVar()".format(Name),globals(),locals())
            com="self.{}_Menu = tk.OptionMenu(self, self.{}_MenuVar, *self.Functions)".format(Name,Name)
            exec(com,globals(),locals())
            com="self.{}_MenuVar.set('x')".format(Name)
            exec(com,globals(),locals())
            com="self.{}.set('x')".format(Name)
            exec(com,globals(),locals())
            com="self.{}_MenuVar.trace_add('write', lambda *_,self=self:self.SetFunc('{}'))".format(Name,Name)
            exec(com,globals(),locals())
            exec("self.{}_Menu.grid(row = {}, column = 2, sticky = 'nsew', padx = 2, pady = 2)".format(Name,count+Nelem),globals(),locals())
            exec("self.{}_Menu.config(width=15)".format(Name),globals(),locals())
            exec("self.All_Entries.append(self.{}_Entry)".format(Name),globals(),locals())
            count=count+1
        return count

    def CreateNptsEntry(self,List,Nelem):
        count=0
        for Name,Properties in List.items():
            exec("self.{}_Label = tk.Label(self, text = '{}:', font = helv_medium)".format(Name,Properties[0]),globals(),locals())
            exec("self.{}_Label.grid(row = {}, column = 0, sticky = 'w', padx = 2, pady = 2)".format(Name,count+Nelem),globals(),locals())
            exec("self.{}= tk.StringVar()".format(Name),globals(),locals())
            com="self.{}_Entry = tk.Entry(self, text = self.{}, disabledbackground = '#eee9e9')".format(Name,Name)
            exec(com,globals(),locals())
            com="self.All_Entries.append(self.{}_Entry)".format(Name)
            exec(com,globals(),locals())
            com="self.{}.trace_add('write', lambda *_,self=self:self.CheckNpts('{}'))".format(Name,Name,Name)
            exec(com,globals(),locals())
            exec("self.{}_Entry.grid(row = {}, column = 1, sticky = 'nsew', padx = 2, pady = 2)".format(Name,count+Nelem),globals(),locals())
            exec("self.All_Entries.append(self.{}_Entry)".format(Name),globals(),locals())
            count=count+1
        return count

    def toggle_mouse(self):
        pass
    def toggle_edit(self):
        self.config(text = 'Editing ' + self.title, fg = 'blue')
        self.Edit_Button.config(text = 'Editing ' + self.title, fg = 'blue')
        self.controller.activate_frame(self)
        for F in self.controller.All_Frames:
            if F != self:
                F.config(text = F.title, fg = 'black')
                F.Edit_Button.config(text = 'Edit ' + F.title, fg = 'black')


    def disable_frame(self):
        self.active_frame = False
        self.active_mouse.set(0)
        self.Edit_Button.config(state = 'disabled')
        for entry in self.All_Entries:
            entry.config(state = 'disabled')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'disabled')

    def update_frame(self):
        if self.active_frame == False:
            self.enable_frame()
        else:
            self.disable_frame()

    def enable_frame(self):
        self.active_frame = True
        self.active_mouse.set(1)
        for entry in self.All_Entries:
            entry.config(state = 'normal')

    def SetFunc(self,Name):
        if self.Verbose: print('SetFunc:',Name)
        exec("self.{}.set(self.{}_MenuVar.get())".format(Name,Name))

    def CheckNpts(self,Name):
        self.val=''
        exec("self.val=self.{}.get()".format(Name,Name),globals(),locals())
        try:
            if int(self.val)<1:
                ExpValid=False
            else:
                ExpValid=True
        except:
            ExpValid=False

        if ExpValid:
                exec("self.{}_Entry.config(bg = 'GREEN')".format(Name),globals(),locals())
        else:
                exec("self.{}_Entry.config(bg = 'RED')".format(Name),globals(),locals())
        return ExpValid

    def CheckFunc(self,Name):
        self.val=''
        exec("self.val=self.{}.get()".format(Name,Name),globals(),locals())

        try:
            com='lambda {}:{}'.format(self.val.split(',')[0],self.val.split(',')[1])
            if self.Verbose: print(com)
            eval(com)
            ExpValid=True
        except:
            print('Unable to parse function {} for entry {}. Lambda function expected:"x,f(x)"'.format(self.val,Name))
            ExpValid=False

        if ExpValid:
                exec("self.{}_Entry.config(bg = 'GREEN')".format(Name),globals(),locals())
        else:
                exec("self.{}_Entry.config(bg = 'RED')".format(Name),globals(),locals())
        return ExpValid

    def SetValues(self):
        # TODO: Provide default values when starting GUI with no provided parameter file.
        for Name,Properties in self.ListNpts.items():
            com="self.{}.set(str(self.IngridSession.yaml['{}']['{}']['{}']))".format(Name,Properties[1][0],Properties[1][1],Properties[1][2])
            if self.Verbose: print(com)
            exec(com,globals(),locals())
        for Name,Properties in self.ListFunc.items():
            com="self.{}.set(str(self.IngridSession.yaml['{}']['{}']['{}']))".format(Name,Properties[1][0],Properties[1][1],Properties[1][2])
            if self.Verbose: print(com)
            exec(com,globals(),locals())


    def GetValues(self):
        for Name,Properties in self.ListNpts.items():
            if not self.CheckNpts(Name): return False
            com="self.IngridSession.yaml['{}']['{}']['{}']=int(self.{}.get())".format(Properties[1][0],Properties[1][1],Properties[1][2],Name)
            if self.Verbose: print(com)
            exec(com,globals(),locals())
        for Name,Properties in self.ListFunc.items():
            if not self.CheckFunc(Name): return False
            com="self.IngridSession.yaml['{}']['{}']['{}']=self.{}.get()".format(Properties[1][0],Properties[1][1],Properties[1][2],Name)
            if self.Verbose: print(com)
            exec(com,globals(),locals())
        return True

