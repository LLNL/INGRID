#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 08:58:43 2019

@author: garcia299
"""

from __future__ import print_function

from sys import platform as sys_pf
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from scipy.optimize import root 
from numpy import full_like
import f90nml

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
import Root_Finder

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

    def __init__(self, master = None, *args, **kwargs):
        
        self.nml = { 'files' : {} , 'grid_params' : {} }
        self.IngridSession = Ingrid.Ingrid(self.nml)

        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (FilePicker, ParamPicker):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(ParamPicker)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

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
                                 disabledbackground = '#f8f8ff', state = 'disabled')
        self.FP_Button = tk.Button(self, text = self.FP_ButtonText.get(), width = 20, \
                                   font = helv_medium, command = params['ButtonCommand'])

        self.FP_Entry.grid(row = 0, column = 0, columnspan = 2, padx = 10, pady = 10, sticky = 'NSEW')
        self.FP_Button.grid(row = 0, column = 2, padx = 10, pady = 10, sticky = 'EW')

        self.isLoaded = False

    def fileLoaded(self, Path):
        self.FP_EntryText.set(str(Path))
        self.FP_Button.config(fg = 'lime green')
        self.isLoaded = True

    def invalidPath(self, Path):
        message = "Path is invalid: '{}'".format(str(Path))
        self.FP_EntryText.set(message)
        self.FP_Button.config(fg = 'red')
        self.isLoaded = False
    def invalidType(self, Path):
        message = "File is not *.txt: '{}'".format(str(Path.name))
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
        
        # Quick access to top-level window.
        # This ROOT is our IngridApp instance.
        self.ROOT = self.winfo_toplevel()

        # TODO: Use controller over self.ROOT. 
        #       Work towards code uniformity.
        self.controller = controller

        self.preview_loaded = False
        
        Title = tk.Label(self, text = "Ingrid.", font=helv_large)
        Title.grid(row = 0, column = 0, columnspan = 3, \
            padx = 10, pady = 10, sticky = 'NSEW')
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
        self.paramFrame = FilePickerFrame(self, controller,\
            {'ButtonText' : 'Load Parameter File', \
             'EntryText'  : 'Select a pre-existing *.nml file.', \
             'ButtonCommand' : self.load_param_file \
             })

        self.FP_Frames = [self.eqdskFrame, self.itpFrame, \
                          self.otpFrame, self.paramFrame]
        
        for i in range(len(self.FP_Frames)):
            self.FP_Frames[i-1].grid(row = i + 1, column = 0, \
                padx = 10, pady = 5, sticky = 'NSEW')

        self.ControlPanel = tk.Frame(self)
        self.ControlPanel.grid(row = len(self.FP_Frames) + 1, column = 0, \
                padx = 10, pady = 5, sticky = 'NSEW')
        self.previewButton = tk.Button(self.ControlPanel, text = 'Preview Loaded Files', \
            font = helv_medium, state = 'disabled', command = self.preview_data)
        self.confirmButton = tk.Button(self.ControlPanel, text = 'Confirm', \
            font = helv_medium, state = 'disabled', command = self.confirm_data)

        self.previewButton.grid(row = 0, column = 0, columnspan = 1, padx = 10, pady = 10, \
            sticky = 'NSEW')
        self.confirmButton.grid(row = 0, column = 1, padx = 10, pady = 10, sticky = 'NSEW')
         

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
                print(self.ROOT.nml)
                self.ROOT.nml['files']['geqdsk'] = str(eqdsk_file)
                self.eqdskFrame.fileLoaded(eqdsk_file)
                self.update_frame_state()
        elif not showFileDialog:
            eqdsk_file = Path(toLoad)
            valid_path = eqdsk_file.is_file()
            if valid_path:
                print(self.ROOT.nml)
                self.ROOT.nml['files']['geqdsk'] = str(eqdsk_file)
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
                self.ROOT.nml['files']['itp'] = str(itp_file)
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
                self.ROOT.nml['files']['itp'] = str(itp_file)
                self.itpFrame.fileLoaded(itp_file)
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
                self.ROOT.nml['files']['otp'] = str(otp_file)
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
                self.ROOT.nml['files']['otp'] = str(otp_file)
                self.otpFrame.fileLoaded(otp_file)
                self.update_frame_state()  

    def load_param_file(self):
        """
        Loads *.nml file containing INGRID parameters.
        Said *.nml contains session settings.

        Post-Call:
        ----------
            Updates param_loaded flag to True if
            a *.nml file was loaded successfully.

        TODO:
        -----
            Check whether the *.nml file actually
            contained all the required files. This
            can be done by inspecting the nml object.

        """
        param_file, valid_path = self.get_file()
        if valid_path and param_file.suffix == '.nml':
            nml = f90nml.read(str(param_file))
            self.ROOT.nml = nml
            self.check_nml(nml)
            self.paramFrame.fileLoaded(param_file)
            self.update_frame_state()

    def check_nml(self, nml):
            lookup = { \
                'geqdsk' : self.load_eqdsk_file, \
                'itp' : self.load_itp_file, \
                'otp' : self.load_otp_file \
                }
            try:
                for file in nml['files']:
                    lookup[file](showFileDialog = False, toLoad = nml['files'][file])
            except:
                print('Did not find file: ' + file)

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
        
        self.ROOT.IngridSession = Ingrid.Ingrid(nml = self.ROOT.nml)
        self.ROOT.IngridSession.OMFIT_read_psi()

        if self.itpFrame.isLoaded:
            self.ROOT.IngridSession.read_target_plate()
        if self.otpFrame.isLoaded:
            self.ROOT.IngridSession.read_target_plate()
        self.ROOT.IngridSession.calc_efit_derivs()
        self.ROOT.IngridSession.plot_efit_data()
        self.ROOT.IngridSession.plot_target_plate()

        self.preview_loaded = True

    def confirm_data(self):
        for item in self.controller.frames[ParamPicker].RF_Frames:
            item.Edit_Button.config(state = 'normal')
        self.controller.frames[ParamPicker].loadParameters()
        self.ROOT.show_frame(ParamPicker)
        self.ROOT.geometry('805x485')
        if not self.preview_loaded:
            self.preview_data()
        self.ROOT.IngridSession.find_roots(tk_controller = self.controller)

class ParamPicker(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.ROOT = self.winfo_toplevel()
        self.controller = controller

        # TODO: Create a ParameterFrame Class. This will be useful when we have
        #       multiple X-points to deal with.
        self.RF_Container = tk.LabelFrame(self, text = '1.', font = helv_medium, relief = 'groove')
        self.PF_Container = tk.LabelFrame(self, text = '2.', font = helv_medium)
        self.Controls_Container = tk.LabelFrame(self, text = 'Settings', font = helv_medium)

        self.RF_Container.grid(row = 1, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PF_Container.grid(row = 0, rowspan = 2, column = 1, padx = 10, pady = 10, sticky = 'nsew')
        self.Controls_Container.grid(row = 0, column = 0, padx = 10, pady = 10, sticky = 'NSew')

        self.createPatches_Button = tk.Button(self.Controls_Container, text = 'Generate Patches', \
                                              font = helv_medium, state = 'disabled', \
                                              command = self.createPatches)
        self.save_Button = tk.Button(self.Controls_Container, text = 'Save Parameters', \
                                     font = helv_medium, state = 'disabled', \
                                     command = self.saveParameters)
        self.load_Button = tk.Button(self.Controls_Container, text = 'Load Files', font = helv_medium, \
                                     command = self.load_files)
        self.createPatches_Button.grid(row = 0, column = 1, padx = 2, pady = 2)
        self.save_Button.grid(row = 0, column = 2, padx = 2, pady = 2)
        self.load_Button.grid(row = 0, column = 0, padx = 2, pady = 2)

        self.MagFrame = RootFinderFrame(self.RF_Container, self, title = 'Magnetic-Axis', style = helv_medium)
        self.XptFrame = RootFinderFrame(self.RF_Container, self, title = 'Primary X-Point', style = helv_medium)
        self.PsiMaxFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Max', style = helv_medium)        
        self.PsiMinFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Min', style = helv_medium)
        self.PsiPrivateFrame = PsiFinderFrame(self.PF_Container, self, title = 'Psi-Private', style = helv_medium)
        self.RF_Frames = [self.MagFrame, self.XptFrame]
        self.PF_Frames = [self.PsiMaxFrame, self.PsiMinFrame, self.PsiPrivateFrame]

        self.All_Frames = []

        for F in self.RF_Frames:
            self.All_Frames.append(F)
        for F in self.PF_Frames:
            self.All_Frames.append(F)
        
        self.acceptRF_Button = tk.Button(self.RF_Container, text = 'Confirm Entries', font = helv_medium, \
                                         command = self.set_RFValues)
        self.acceptRF_Button.grid(row = 3, column = 0, columnspan = 1, padx = 10, pady = 10)

        self.acceptPF_Button = tk.Button(self.PF_Container, text = 'Confirm Entries', font = helv_medium, \
                                         command = self.set_PsiValues)
        self.acceptPF_Button.grid(row = 6, column = 0, columnspan = 1, padx = 10, pady = 10)

        self.MagFrame.grid(row = 1, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.XptFrame.grid(row = 2, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PsiMaxFrame.grid(row = 3, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PsiMinFrame.grid(row = 4, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PsiPrivateFrame.grid(row = 5, column = 0, padx = 10, pady = 10, sticky = 'nsew')

        # Use push/pop system for keeping track. Initially all inactive
        self.ActiveFrame = []
        self.curr_click = ()
    
    def refine(self):
        sol = root(self.controller.IngridSession.root_finder.func, \
                [float(self.ActiveFrame[0].R_EntryText.get()), \
                 float(self.ActiveFrame[0].Z_EntryText.get())] )
        self.ActiveFrame[0].R_EntryText.set('{0:.12f}'.format(sol.x[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.12f}'.format(sol.x[1]))
        

    def activate_frame(self, calling_frame):
        if self.ActiveFrame:
            self.ActiveFrame[0].disable_frame()
            self.ActiveFrame[0].toggle_mouse()
            self.ActiveFrame[0].Edit_Button.config(state = 'normal')
            self.ActiveFrame.pop()
        calling_frame.enable_frame()
        self.ActiveFrame.append(calling_frame)
        self.ActiveFrame[0].toggle_mouse()

    def unlock_PF_Frame(self):
        for F in self.PF_Frames:
            F.Edit_Button.config(state = 'normal')
        self.controller.IngridSession.calc_psinorm()

    def unlock_controls(self):
        self.createPatches_Button.config(state = 'normal')
        self.save_Button.config(state = 'normal') 

    def load_files(self):
        self.controller.frames[FilePicker].preview_loaded = False
        self.controller.show_frame(FilePicker)
        self.controller.geometry("550x340")

    def update_root_finder(self):
        self.ActiveFrame[0].R_EntryText.set('{0:.12f}'.format(self.curr_click[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.12f}'.format(self.curr_click[1]))
        if self.ActiveFrame[0] in self.PF_Frames:
            psi_magx = self.controller.IngridSession.efit_psi.get_psi(self.MagAxis[0], self.MagAxis[1])
            psi_xpt = self.controller.IngridSession.efit_psi.get_psi(self.Xpt[0], self.Xpt[1])
            psi_efit = self.controller.IngridSession.efit_psi.get_psi(self.curr_click[0],self.curr_click[1])
            psi = (psi_efit - full_like(psi_efit, psi_magx))/(psi_xpt - psi_magx)
            self.ActiveFrame[0].Psi_EntryText.set('{0:.12f}'.format(psi))

    def set_RFValues(self):
        self.controller.nml['grid_params']['rmagx'] = float(self.MagFrame.R_EntryText.get())
        self.controller.nml['grid_params']['zmagx'] = float(self.MagFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['rxpt'] = float(self.XptFrame.R_EntryText.get())
        self.controller.nml['grid_params']['zxpt'] = float(self.XptFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['step_ratio'] = 0.02

        self.MagAxis = (self.controller.nml['grid_params']['rmagx'], self.controller.nml['grid_params']['zmagx'])
        self.Xpt = (self.controller.nml['grid_params']['rxpt'], self.controller.nml['grid_params']['zxpt'])

        self.controller.IngridSession.magx = self.MagAxis
        self.controller.IngridSession.xpt1 = self.Xpt
        
        # TODO: Exception handling behind the scenes for ensuring PF_Frame is indeed ready.
        
        self.acceptRF_Button.config(text = 'Entries Saved!', fg = 'lime green')
        self.unlock_PF_Frame()
    
    def set_PsiValues(self):
        
        for F in self.PF_Frames:
            if F.R_EntryText.get() == '':
                F.R_EntryText.set('0.0')
            if F.Z_EntryText.get() == '':
                F.Z_EntryText.set('0.0')

        self.controller.nml['grid_params']['psi_max_r'] = float(self.PsiMaxFrame.R_EntryText.get())
        self.controller.nml['grid_params']['psi_max_z'] = float(self.PsiMaxFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['psi_max'] = float(self.PsiMaxFrame.Psi_EntryText.get())

        self.controller.nml['grid_params']['psi_min_core_r'] = float(self.PsiMinFrame.R_EntryText.get())
        self.controller.nml['grid_params']['psi_min_core_z'] = float(self.PsiMinFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['psi_min_core'] = float(self.PsiMinFrame.Psi_EntryText.get())

        self.controller.nml['grid_params']['psi_min_pf_r'] = float(self.PsiPrivateFrame.R_EntryText.get())
        self.controller.nml['grid_params']['psi_min_pf_z'] = float(self.PsiPrivateFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['psi_min_pf'] = float(self.PsiPrivateFrame.Psi_EntryText.get())
        self.controller.nml['grid_params']['step_ratio'] = 0.02
        
        self.unlock_controls()
        self.acceptPF_Button.config(text = 'Entries Saved!', fg = 'lime green')

    def createPatches(self):

        debug = False

        import time as ti

        print(self.ROOT.nml)
        self.set_RFValues()
        self.set_PsiValues()
        self.controller.IngridSession.magx = self.MagAxis
        self.controller.IngridSession.xpt1 = self.Xpt
        self.controller.IngridSession.calc_psinorm()
        self.winfo_toplevel().nml = f90nml.read(str(self.winfo_toplevel().frames[FilePicker].eqdskFrame.FP_EntryText.get()))
        
        patch_time = ti.time()
        self.controller.IngridSession.construct_SNL_patches()
        patch_time = ti.time() - patch_time

        subgrid_time = ti.time()
        self.controller.IngridSession.construct_SNL_subgrid()
        subgrid_time = ti.time() - subgrid_time

        patchDiagram_time = ti.time()
        self.controller.IngridSession.patch_diagram()
        patchDiagram_time = ti.time() - patchDiagram_time

        gridDiagram_time = ti.time()
        self.controller.IngridSession.grid_diagram()
        gridDiagram_time = ti.time() - gridDiagram_time

        export_time = ti.time()
        self.controller.IngridSession.export()
        export_time = ti.time() - export_time

        print('========================================================')
        print('SUMMARY:')
        print("'construct_SNL_patches' took: {}s".format(patch_time))
        print("'construct_SNL_subgrid' took: {}s".format(subgrid_time))
        print("'patch_diagram' took: {}s".format(patchDiagram_time))
        print("'grid_diagram' took: {}s".format(gridDiagram_time))
        print("'export' took: {}s".format(export_time))
        print('========================================================')


    def saveParameters(self):
        self.set_RFValues()
        self.set_PsiValues()
        fname = tkFD.asksaveasfilename(initialdir = '.', title = 'Save File',defaultextension ='.nml')
        print(fname)
        f90nml.write(self.controller.nml, fname, force=True)  # force tag is to overwrite the previous file
        print("Saved parameters to '{}'.".format(fname))

    def loadParameters(self):
        lookup = {'rmagx' : self.MagFrame.R_EntryText, 'zmagx' : self.MagFrame.Z_EntryText,\
                  'rxpt'  : self.XptFrame.R_EntryText, 'zxpt'  : self.XptFrame.Z_EntryText,\
                  'psi_max_r' : self.PsiMaxFrame.R_EntryText, 'psi_max_z' : self.PsiMaxFrame.Z_EntryText,\
                  'psi_min_core_r' : self.PsiMinFrame.R_EntryText, 'psi_min_core_z' : self.PsiMinFrame.Z_EntryText,\
                  'psi_min_pf_r' : self.PsiPrivateFrame.R_EntryText, 'psi_min_pf_z' : self.PsiPrivateFrame.Z_EntryText,\
                  'psi_max' : self.PsiMaxFrame.Psi_EntryText, 'psi_min_core' : self.PsiMinFrame.Psi_EntryText,\
                  'psi_min_pf' : self.PsiPrivateFrame.Psi_EntryText\
                  }
        for param in self.controller.nml['grid_params']:
            try:
                lookup[param].set('{0:.12f}'.format(self.controller.nml['grid_params'][param]))
            except:
                print('Did not find: ' + str(param))


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
        self.Edit_Button = tk.Button(self, text = 'Set ' + self.title, font = helv_medium, \
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
        controller.activate_frame(self)
        for F in controller.All_Frames:
            if F != self:
                F.config(text = F.title, fg = 'black')

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

        self.Edit_Button = tk.Button(self, text = 'Set ' + self.title, font = helv_medium, \
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

        self.Edit_Button.grid(row = 2, column = 2, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        rbWrapper.grid(row = 0, rowspan = 2, column = 2, columnspan = 1, sticky = 'ns', padx = 10, pady = 5)
        self.On_Radiobutton.grid(row = 1, column = 2, sticky = 'we', padx = 2, pady = 2)
        self.Off_Radiobutton.grid(row = 2, column = 2, sticky = 'we', padx = 2, pady = 2)

        self.disable_frame()

    def toggle_mouse(self):
        if self.active_mouse.get() is True:
            self.controller.controller.IngridSession.root_finder.connect()
            print('Connected')
        else:
            self.controller.controller.IngridSession.root_finder.disconnect()
            print('Disconnected')

    def toggle_edit(self):
        self.config(text = 'Editing ' + self.title, fg = 'blue')
        self.controller.activate_frame(self)
        for F in self.controller.All_Frames:
            if F != self:
                F.config(text = F.title, fg = 'black')

    def disable_frame(self):
        self.active_frame = False
        self.active_mouse.set(0)
        self.Edit_Button.config(state = 'disabled')
        for entry in self.All_Entries:
            entry.config(state = 'disabled')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'disabled')

    def enable_frame(self):
        self.active_frame = True
        self.active_mouse.set(1)
        for entry in self.All_Entries:
            entry.config(state = 'normal')
        for rb in self.All_Radiobuttons:
            rb.config(state = 'normal')

