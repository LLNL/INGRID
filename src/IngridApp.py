#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 08:58:43 2019

@author: garcia299
"""

from __future__ import print_function
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from scipy.optimize import root 
from numpy import full_like
import f90nml

from pathlib2 import Path

import Tkinter as tk
import tkFileDialog as tkFD
import tkMessageBox as tkMB

import Ingrid
import Root_Finder

LARGE_FONT = 'Helvetica 13 bold'
MEDIUM_FONT = 'Helvetica 9 bold'

class IngridApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        
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

class FilePicker(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent) 
        
        # Quick access to top-level window.
        # This ROOT is our IngridApp instance.
        self.ROOT = self.winfo_toplevel()

        # TODO: Use controller over self.ROOT. 
        #       Work towards code uniformity.
        self.controller = controller

        # Bool flags for keeping track of which files are loaded.
        self.eqdsk_loaded = False
        self.itp_loaded = False
        self.otp_loaded = False
        self.param_loaded = False

        self.preview_loaded = False
        
        # TODO: Gather all buttons in a list or dictionary, and
        #       populate the 'StartPage' with a method? 
        #       Grid-like appearance would be a good idea...
        label = tk.Label(self, text = "Ingrid.", font=LARGE_FONT)
        label.pack(side = 'top', padx = 5, pady = 3)

        self.container1 = tk.Frame(self)
        self.eqdsk_button = tk.Button(self.container1, anchor = 'w', width = 16, \
                                    text = "Select EQDSK File", command=self.load_eqdsk_file, \
                                    font = MEDIUM_FONT)
        self.eqdsk_text = tk.StringVar()
        self.eqdsk_text.set('Select an EQDSK file.')
        self.eqdsk_entry = tk.Entry(self.container1, text = self.eqdsk_text, width = 40, \
                                  state = 'disabled', disabledbackground = '#f8f8ff')

        self.container1.pack(side = 'top', fill = tk.X, padx = 5, pady = 2)
        self.eqdsk_button.pack(side = 'right', padx = 5,   pady = 2)
        self.eqdsk_entry.pack(side = 'left', fill = tk.X)

        self.container2 = tk.Frame(self)
        self.itp_button = tk.Button(self.container2, anchor = 'w', width = 16, \
                                    text = 'Select Inner Plate', command = self.load_itp_file, \
                                    font = MEDIUM_FONT)
        self.itp_text = tk.StringVar()
        self.itp_text.set('Select an inner strike plate configuration file.')
        self.container2.pack(side = 'top', fill = tk.X, padx = 5,   pady = 2)
        self.itp_entry = tk.Entry(self.container2, text = self.itp_text, width = 40, \
                                  state = 'disabled', disabledbackground = '#f8f8ff')
        self.itp_button.pack(side = 'right', padx = 5,   pady = 2)
        self.itp_entry.pack(side = 'left', fill = tk.X)

        self.container3 = tk.Frame(self)
        self.otp_button = tk.Button(self.container3, anchor = 'w', width = 16, \
                                    text = 'Select Outer Plate', command = self.load_otp_file, \
                                    font = MEDIUM_FONT)
        self.otp_text = tk.StringVar()
        self.otp_text.set('Select an outer strike plate configuration file.')
        self.container3.pack(side = 'top', fill = tk.X, padx = 5,   pady = 2)
        self.otp_entry = tk.Entry(self.container3, text = self.otp_text, width = 40, \
                                  state = 'disabled', disabledbackground = '#f8f8ff') 
        self.otp_button.pack(side = 'right', padx = 5,   pady = 2)
        self.otp_entry.pack(side = 'left', fill = tk.X)
        
        self.container4 = tk.Frame(self)
        self.param_button = tk.Button(self.container4, anchor = 'w', width = 16, \
                                    text = 'Load Parameter File', command = self.load_param_file, \
                                    font = MEDIUM_FONT)
        self.param_text = tk.StringVar()
        self.param_text.set('Load a pre-existing parameter *.nml file.')
        self.container4.pack(side = 'top', fill = tk.X, padx = 5,   pady = 2)
        self.param_entry = tk.Entry(self.container4, text = self.param_text, width = 40, \
                                    state = 'disabled', disabledbackground = '#f8f8ff')
        self.param_button.pack(side = 'right', padx = 5,   pady = 2)
        self.param_entry.pack(side = 'left', fill = tk.X)

        self.container5 = tk.Frame(self)
        self.container5.pack(side = 'top', padx = 5,   pady = 2)
        self.preview_button = tk.Button(self.container5, text = 'Preview Data', command = self.preview_data, \
                               font = MEDIUM_FONT, state = 'disabled')
        self.preview_button.pack( pady = 2, padx = 5, side = 'left')

        self.proceed_button = tk.Button(self.container5, text = 'Confirm', state = 'disabled', \
                                command = self.next_page)
        self.proceed_button.pack(side = 'right')
        print(self.master)

    def get_file(self, ftypes = (("All Files", "*"), ("Text Files", "*.txt"))):
        try:
            fpath = Path(tkFD.askopenfilename(initialdir = '.', title = 'Select File', \
                         filetypes = ftypes ))
            if fpath.is_file():
                return fpath, True
            else:
                raise TypeError
        except TypeError:
            return '', False

    def load_eqdsk_file(self): 
        eqdsk_file, valid_path = self.get_file()
        if valid_path and eqdsk_file.suffix == '':
            self.ROOT.nml['files']['geqdsk'] = str(eqdsk_file)
            print(self.ROOT.nml)
            self.eqdsk_button.config(fg = 'lime green')
            self.eqdsk_text.set(str(eqdsk_file))
            self.eqdsk_loaded = True
            self.update_frame_state()
        elif eqdsk_file.suffix != '':
            self.eqdsk_button.config(fg = 'red')
            self.eqdsk_text.set("Invalid File Type '" + eqdsk_file.suffix + "' selected.")
            self.eqdsk_loaded = False
            self.update_frame_state()

    def load_itp_file(self):
        itp_file, valid_path = self.get_file()
        if valid_path and itp_file.suffix == '.txt':
            self.ROOT.nml['files']['itp'] = str(itp_file)
            print(self.ROOT.nml)
            self.itp_button.config(fg = 'lime green')
            self.itp_text.set(str(itp_file))
            self.itp_loaded = True
            self.update_frame_state()

    def load_otp_file(self):
        otp_file, valid_path = self.get_file()
        if valid_path and otp_file.suffix == '.txt':
            self.ROOT.nml['files']['otp'] = str(otp_file)
            print(self.ROOT.nml)
            self.otp_button.config(fg = 'lime green')
            self.otp_text.set(str(otp_file))
            self.otp_loaded = True
            self.update_frame_state()

    def load_param_file(self):
        param_file, valid_path = self.get_file( (('NML Files', '*.nml'), ('All Files', '*')) )
        if valid_path and param_file.suffix == '.nml':
            self.ROOT.nml = f90nml.read(str(param_file))
            print(self.ROOT.nml)
            self.param_button.config(fg = 'lime green')
            self.param_text.set(str(param_file))
            self.param_loaded = True
            self.update_all_widgets()
            self.update_frame_state()

    def update_frame_state(self):
        # Button state bool flag would be nice here.
        # We'd be able to call a 'flip_state' function instead of
        # manually setting the config.
        if self.eqdsk_loaded:
            self.preview_button.config(state = 'normal')
        if not self.eqdsk_loaded:
            self.preview_button.config(state = 'disabled')
        if self.files_ready():
            self.proceed_button.config(state = 'normal')
    
    def update_all_widgets(self):
        self.eqdsk_button.config(fg = 'lime green')
        self.eqdsk_text.set(self.ROOT.nml['files']['geqdsk'])
        self.eqdsk_loaded = True

        self.itp_button.config(fg = 'lime green')
        self.itp_text.set(self.ROOT.nml['files']['itp'])
        self.itp_loaded = True

        self.otp_button.config(fg = 'lime green')
        self.otp_text.set(self.ROOT.nml['files']['otp'])
        self.otp_loaded = True

    def files_ready(self):
        if self.param_loaded and self.preview_loaded:
            return True
        if self.eqdsk_loaded and self.itp_loaded and self.otp_loaded and self.preview_loaded:
            return True
        return False

    def preview_data(self):
        # TODO: Make sure ONLY if self.ROOT.nml good to go.. execute this code. 
         
        # self.ROOT.IngridSession.setup()
        print(self.ROOT.nml)
        self.ROOT.IngridSession = Ingrid.Ingrid(nml = self.ROOT.nml)
        self.ROOT.IngridSession.OMFIT_read_psi()
        print(self.ROOT.IngridSession.efit_psi)
        self.ROOT.IngridSession.read_target_plate()
        self.ROOT.IngridSession.calc_efit_derivs()
        self.ROOT.IngridSession.plot_efit_data()
        self.ROOT.IngridSession.plot_target_plate()

        self.preview_loaded = True
        self.update_frame_state()

    def next_page(self):
        for item in self.controller.frames[ParamPicker].All_RefineFrames:
            item.Edit_Button.config(state = 'normal')

        
        self.ROOT.IngridSession.find_roots(tk_controller = self.controller)
        self.ROOT.show_frame(ParamPicker)
        self.ROOT.geometry('760x530')

class ParamPicker(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.ROOT = self.winfo_toplevel()
        self.controller = controller

        # TODO: Create a ParameterFrame Class. This will be useful when we have
        #       multiple X-points to deal with.
        self.RF_Container = tk.LabelFrame(self, text = '1.', font = MEDIUM_FONT, relief = 'groove')
        self.PF_Container = tk.LabelFrame(self, text = '2.', font = MEDIUM_FONT)
        self.Controls_Container = tk.Frame(self)

        self.RF_Container.grid(row = 1, column = 0, padx = 10, pady = 10, sticky = 'nsew')
        self.PF_Container.grid(row = 0, rowspan = 2, column = 1, padx = 10, pady = 10, sticky = 'nsew')
        self.Controls_Container.grid(row = 0, column = 0, padx = 10, pady = 20, stick = 'sew')

        self.createPatches_Button = tk.Button(self.Controls_Container, text = 'Generate Patches', \
                                              font = MEDIUM_FONT, state = 'disabled', \
                                              command = self.createPatches)
        self.save_Button = tk.Button(self.Controls_Container, text = 'Save Parameters', \
                                     font = MEDIUM_FONT, state = 'disabled', \
                                     command = self.saveParameters)
        self.load_Button = tk.Button(self.Controls_Container, text = 'Load Files', font = MEDIUM_FONT, \
                                     command = self.load_files)
        self.createPatches_Button.grid(row = 1, column = 1, padx = 2, pady = 2, sticky = 'nsew')
        self.save_Button.grid(row = 1, column = 2, padx = 2, pady = 2, sticky = 'nsew')
        self.load_Button.grid(row = 1, column = 0, padx = 2, pady = 2, sticky = 'nsew')

        self.MagFrame = RefineFrame(self.RF_Container, self, title = 'Magnetic-Axis', style = MEDIUM_FONT)
        self.XptFrame = RefineFrame(self.RF_Container, self, title = 'Primary X-Point', style = MEDIUM_FONT)
        self.PsiMaxFrame = PsiFrame(self.PF_Container, self, title = 'Psi-Max', style = MEDIUM_FONT)        
        self.PsiMinFrame = PsiFrame(self.PF_Container, self, title = 'Psi-Min', style = MEDIUM_FONT)
        self.PsiPrivateFrame = PsiFrame(self.PF_Container, self, title = 'Psi-Private', style = MEDIUM_FONT)
        self.All_RefineFrames = [self.MagFrame, self.XptFrame]
        self.All_PsiFrames = [self.PsiMaxFrame, self.PsiMinFrame, self.PsiPrivateFrame]

        self.All_Frames = []

        for F in self.All_RefineFrames:
            self.All_Frames.append(F)
        for F in self.All_PsiFrames:
            self.All_Frames.append(F)
        
        self.acceptRF_Button = tk.Button(self.RF_Container, text = 'Confirm Entries', font = MEDIUM_FONT, \
                                         command = self.set_RFValues)
        self.acceptRF_Button.grid(row = 3, column = 0, columnspan = 1, padx = 10, pady = 10)

        self.acceptPF_Button = tk.Button(self.PF_Container, text = 'Confirm Entries', font = MEDIUM_FONT, \
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
        self.ActiveFrame[0].R_EntryText.set('{0:.7f}'.format(sol.x[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.7f}'.format(sol.x[1]))
        

    def activate_frame(self, calling_frame):
        if self.ActiveFrame:
            self.ActiveFrame[0].disable_frame()
            self.ActiveFrame[0].toggle_mouse()
            self.ActiveFrame[0].Edit_Button.config(state = 'normal')
            self.ActiveFrame.pop()
        calling_frame.enable_frame()
        self.ActiveFrame.append(calling_frame)
        self.ActiveFrame[0].toggle_mouse()

    def unlock_PsiFrame(self):
        for F in self.All_PsiFrames:
            F.Edit_Button.config(state = 'normal')
        self.controller.IngridSession.calc_psinorm()

    def unlock_controls(self):
        self.createPatches_Button.config(state = 'normal')
        self.save_Button.config(state = 'normal') 

    def load_files(self):
        self.controller.frames[FilePicker].preview_loaded = False
        self.controller.show_frame(FilePicker)
        self.controller.geometry('450x215')

    def update_root_finder(self):
        self.ActiveFrame[0].R_EntryText.set('{0:.7f}'.format(self.curr_click[0]))
        self.ActiveFrame[0].Z_EntryText.set('{0:.7f}'.format(self.curr_click[1]))
        if self.ActiveFrame[0] in self.All_PsiFrames:
            psi_magx = self.controller.IngridSession.efit_psi.get_psi(self.MagAxis[0], self.MagAxis[1])
            psi_xpt = self.controller.IngridSession.efit_psi.get_psi(self.Xpt[0], self.Xpt[1])
            psi_efit = self.controller.IngridSession.efit_psi.get_psi(self.curr_click[0],self.curr_click[1])
            psi = (psi_efit - full_like(psi_efit, psi_magx))/(psi_xpt - psi_magx)
            self.ActiveFrame[0].Psi_EntryText.set('{0:.7f}'.format(psi))
            self.controller.IngridSession.compute_eq_psi()

    def set_RFValues(self):
        self.controller.nml['grid_params']['Rmagx'] = float(self.MagFrame.R_EntryText.get())
        self.controller.nml['grid_params']['Zmagx'] = float(self.MagFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['Rxpt'] = float(self.XptFrame.R_EntryText.get())
        self.controller.nml['grid_params']['Zxpt'] = float(self.XptFrame.Z_EntryText.get())
        self.controller.nml['grid_params']['step_ratio'] = 0.02

        self.MagAxis = (self.controller.nml['grid_params']['Rmagx'], self.controller.nml['grid_params']['Zmagx'])
        self.Xpt = (self.controller.nml['grid_params']['Rxpt'], self.controller.nml['grid_params']['Zxpt'])

        self.controller.IngridSession.magx = self.MagAxis
        self.controller.IngridSession.xpt1 = self.Xpt
        
        # TODO: Exception handling behind the scenes for ensuring PsiFrame is indeed ready.
        
        self.acceptRF_Button.config(text = 'Entries Saved!', fg = 'lime green')
        self.unlock_PsiFrame()
    
    def set_PsiValues(self):
        self.controller.nml['grid_params']['psi_max'] = float(self.PsiMaxFrame.Psi_EntryText.get())
        self.controller.nml['grid_params']['psi_min_core'] = float(self.PsiMinFrame.Psi_EntryText.get())
        self.controller.nml['grid_params']['psi_min_pf'] = float(self.PsiPrivateFrame.Psi_EntryText.get())
        self.controller.nml['grid_params']['step_ratio'] = 0.02

        self.unlock_controls()
        self.acceptPF_Button.config(text = 'Entries Saved!', fg = 'lime green')

    def createPatches(self):
        self.set_RFValues()
        self.set_PsiValues()
        self.controller.IngridSession.magx = self.MagAxis
        print('Magnetic Axis: ' + str(self.MagAxis))
        self.controller.IngridSession.xpt1 = self.Xpt
        self.controller.IngridSession.calc_psinorm()
        self.controller.IngridSession.compute_eq_psi()
        self.controller.IngridSession.construct_SNL_patches()
        self.controller.IngridSession.patch_diagram()

    def saveParameters(self):
        self.set_RFValues()
        self.set_PsiValues()
        fname = tkFD.asksaveasfilename(initialdir = '.', title = 'Save File',defaultextension ='.nml')
        print(fname)
        f90nml.write(self.controller.nml, fname, force=True)  # force tag is to overwrite the previous file
        print("Saved parameters to '{}'.".format(fname))

class RefineFrame(tk.LabelFrame):

    def __init__(self, parent, controller, title, style):
        tk.LabelFrame.__init__(self, parent, text = title, font = style, relief = 'flat')

        rbWrapper = tk.LabelFrame(self, text = 'Mouse Entry:', font = MEDIUM_FONT)
        R_Label = tk.Label(self, text = 'R:', font = MEDIUM_FONT)
        Z_Label = tk.Label(self, text = 'Z:', font = MEDIUM_FONT)
        self.title = title
        self.controller = controller
        self.R_EntryText = tk.StringVar()
        self.Z_EntryText = tk.StringVar()
        self.All_EntryText = (self.R_EntryText, self.Z_EntryText)

        self.R_Entry = tk.Entry(self, text = self.R_EntryText, disabledbackground = '#eee9e9')
        self.Z_Entry = tk.Entry(self, text = self.Z_EntryText, disabledbackground = '#eee9e9')
        self.All_Entries = [self.R_Entry, self.Z_Entry]

        self.Refine_Button = tk.Button(self, text = 'Refine ' + self.title, font = MEDIUM_FONT, \
                                       command = controller.refine, width = 15)
        self.Edit_Button = tk.Button(self, text = 'Set ' + self.title, font = MEDIUM_FONT, \
                command = lambda : self.toggle_edit(parent, controller), width = 15)
        self.All_Buttons = (self.Refine_Button, self.Edit_Button)

        self.active_mouse = tk.BooleanVar()
        self.active_mouse.set(0)

        self.active_frame = False

        self.On_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Enabled', font = MEDIUM_FONT, \
                                            variable = self.active_mouse, value = 1, \
                                            command = self.toggle_mouse)
        self.Off_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Disabled', font = MEDIUM_FONT, \
                                            variable = self.active_mouse, value = 0, \
                                            command = self.toggle_mouse)
        self.All_Radiobuttons = (self.On_Radiobutton, self.Off_Radiobutton)

        R_Label.grid(row = 0, column = 0, sticky = 'w', padx = 2, pady = 2)
        Z_Label.grid(row = 1, column = 0, sticky = 'w', padx = 2, pady = 2)
        
        self.R_Entry.grid(row = 0, column = 1, sticky = 'ew', padx = 2, pady = 2)
        self.Z_Entry.grid(row = 1, column = 1, sticky = 'ew', padx = 2, pady = 2)
        
        self.Refine_Button.grid(row = 2, column = 1, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        self.Edit_Button.grid(row = 2, column = 2, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        rbWrapper.grid(row = 0, rowspan = 2, column = 2, columnspan = 1, sticky = 'new', padx = 10, pady = 5)
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

class PsiFrame(tk.LabelFrame):

    def __init__(self, parent, controller, title, style):
        self.title = title
        self.controller = controller
        tk.LabelFrame.__init__(self, parent, text = title, font = style, relief = 'flat')
        R_Label = tk.Label(self, text = 'R:', font = MEDIUM_FONT)
        Z_Label = tk.Label(self, text = 'Z:', font = MEDIUM_FONT)
        Psi_Label = tk.Label(self, text = 'Psi:', font = MEDIUM_FONT)

        rbWrapper = tk.LabelFrame(self, text = 'Mouse Entry:', font = MEDIUM_FONT)
        self.R_EntryText = tk.StringVar()
        self.Z_EntryText = tk.StringVar()
        self.Psi_EntryText = tk.StringVar()
        self.All_EntryText = [self.R_EntryText, self.Z_EntryText, self.Psi_EntryText]

        self.R_Entry = tk.Entry(self, text = self.R_EntryText, disabledbackground = '#eee9e9')
        self.Z_Entry = tk.Entry(self, text = self.Z_EntryText, disabledbackground = '#eee9e9')
        self.Psi_Entry = tk.Entry(self, text = self.Psi_EntryText, disabledbackground = '#eee9e9')
        self.All_Entries = [self.R_Entry, self.Z_Entry, self.Psi_Entry]

        self.Edit_Button = tk.Button(self, text = 'Set ' + self.title, font = MEDIUM_FONT, \
                command = self.toggle_edit, width = 15)

        self.All_Buttons = [self.Edit_Button]

        self.active_mouse = tk.BooleanVar()
        self.active_mouse.set(0)
        self.active_frame = False

        self.On_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Enabled', font = MEDIUM_FONT, \
                                            variable = self.active_mouse, value = 1, \
                                            command = self.toggle_mouse)
        self.Off_Radiobutton = tk.Radiobutton(rbWrapper, text = 'Disabled', font = MEDIUM_FONT, \
                                            variable = self.active_mouse, value = 0, \
                                            command = self.toggle_mouse)
        self.All_Radiobuttons = (self.On_Radiobutton, self.Off_Radiobutton)

        R_Label.grid(row = 0, column = 0, sticky = 'w', padx = 2, pady = 2)
        Z_Label.grid(row = 1, column = 0, sticky = 'w', padx = 2, pady = 2)
        Psi_Label.grid(row = 2, column = 0, sticky = 'w', padx = 2, pady = 2)

        self.R_Entry.grid(row = 0, column = 1, sticky = 'ew', padx = 2, pady = 2)
        self.Z_Entry.grid(row = 1, column = 1, sticky = 'ew', padx = 2, pady = 2)
        self.Psi_Entry.grid(row = 2, column = 1, sticky = 'we', padx = 2, pady = 2)

        self.Edit_Button.grid(row = 2, column = 2, columnspan = 1, sticky = 'we', padx = 2, pady = 2)

        rbWrapper.grid(row = 0, rowspan = 2, column = 2, columnspan = 1, sticky = 'new', padx = 10, pady = 5)
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

