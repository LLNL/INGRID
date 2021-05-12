from __future__ import print_function

from sys import platform as sys_pf
from sys import path

from os.path import getmtime
import matplotlib

try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt

try:
    from pathlib import Path
except:
    from pathlib2 import Path
try:
    import tkinter as tk
except:
    import Tkinter as tk
try:
    import tkinter.filedialog as filedialog
except:
    import tkFileDialog as tkFD
try:
    from tkinter import messagebox as tkMB
except:
    import tkMessageBox as tkMB

import yaml
from INGRID.ingrid import Ingrid
from time import time


class IngridGUI(tk.Tk):
    def __init__(self, master=None, IngridSession=None, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        self.gui_dimensions = {FilePicker: "550x270"}

        # attached to the parent Ingrid instance instead of new instance
        if IngridSession is None:
            self.NewIG()
        else:
            self.Ingrid = IngridSession

        # self.NewSettingsPrompt()
        self.PopulateGUI()

        #automatic loading of input files
        if self.Ingrid.InputFile is not None:
            self.frames[FilePicker].load_param_file(
                ExistingParamFile=self.Ingrid.yaml)

    def ShowFrame(self, item):
        frame = self.frames[item]
        frame.tkraise()
        self.geometry(self.gui_dimensions[FilePicker])

    def NewSettingsPrompt(self):
        fname = ''
        while fname in ['', None]:
            if tkMB.askyesno('', 'Create a new settings file for this INGRID session?'):
                fname = 'INGRID_Session' + str(int(time())) + '.yml'
                fname = filedialog.asksaveasfilename(initialdir='.', title='New YAML',
                    defaultextension='.yml', initialfile=fname)
                if fname not in ['', None]:
                    self.Ingrid.SaveSettingsFile(fname=fname)
                    break
            else:
                break

    def Reset(self, message='Are you sure you want to reset?'):
        if tkMB.askyesno('', message):

            try:
                # Close all open figures
                plt.close('all')
            except:
                pass

            self.PopulateGUI()
            self.ResetIG()
            self.ShowFrame(FilePicker)
            self.geometry(self.gui_dimensions[ParamPicker])
            return True

        else:
            return False

    def Exit(self):
        if tkMB.askyesno('', 'Are you sure you want to quit?'):
            try:
                plt.close('all')
            except:
                pass
            self.destroy()

    def NewIG(self):
        self.Ingrid = Ingrid()

    def ResetIG(self):
        self.Ingrid = None

    def PopulateGUI(self):
        self.frames = {}

        for F in [FilePicker]:
            frame = F(self.container, self)
            frame.Ingrid = self.Ingrid
            self.frames[F] = frame
            #frame.grid(row=0, column=0, sticky="nsew")


class FilePicker(tk.Tk):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller
        self.PreviewPlot = None

        self.EntryFrame = tk.Frame(parent)
        self.ControlFrame = tk.Frame(parent)

        self.EntryFrame.grid(row=0, column=0, padx=10, pady=10)
        self.ControlFrame.grid(row=1, column=0, padx=10, pady=10)

        self.ParamFileLabel_String = tk.StringVar()
        self.ParamFileEntry_String = tk.StringVar()
        self.ParamFileLabel_String.set('Parameter File Path:')

        self.ParamFileLabel = tk.Label(self.EntryFrame,
            text=self.ParamFileLabel_String.get())
        self.ParamFileEntry = tk.Entry(self.EntryFrame,
            text=self.ParamFileEntry_String, width=50,
            disabledbackground='#f8f8ff', state='disabled')
        self.ParamFileButton = tk.Button(self.EntryFrame,
            text='Select Parameter File', command=self.LoadParameterFile)

        self.ParamFileLabel.grid(row=0, column=0, sticky='nsew')
        self.ParamFileEntry.grid(row=0, column=1, sticky='nsew')
        self.ParamFileButton.grid(row=0, column=2, sticky='nsew')

        self.ViewDataButton = tk.Button(self.ControlFrame,
            text='View Loaded File', command=self.ViewData)
        self.AnalyzeTopologyButton = tk.Button(self.ControlFrame,
            text='Analyze Topology', command=self.AnalyzeTopology)
        self.CreatePatchesButton = tk.Button(self.ControlFrame,
            text='Create Patches', command=self.CreatePatches)
        self.CreateSubgridButton = tk.Button(self.ControlFrame,
            text='Create Grid', command=self.CreateSubgrid)
        self.ExportGridueButton = tk.Button(self.ControlFrame,
            text='Export gridue', command=self.ExportGridue)
        self.QuitButton = tk.Button(self.ControlFrame,
            text='Quit', command=self.Quit)

        self.ViewDataButton.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
        self.CreatePatchesButton.grid(row=0, column=2, padx=10, pady=10, sticky='nsew')
        self.CreateSubgridButton.grid(row=0, column=3, padx=10, pady=10, sticky='nsew')
        self.ExportGridueButton.grid(row=0, column=4, padx=10, pady=10, sticky='nsew')
        self.QuitButton.grid(row=0, column=5, padx=10, pady=10, sticky='nsew')

    def NewCase(self):
        self.controller.NewIG()
        self.Ingrid = self.controller.Ingrid

    def ProcessParameterFile(self, fname):
        self.NewCase()
        self.Ingrid.InputFile = fname
        self.Ingrid.PopulateSettings(Ingrid.ReadYamlFile(fname))
        self.ParamFileMtime = getmtime(fname)
        self.ParamFileName = fname

    def LoadParameterFile(self):
        fname = filedialog.askopenfilename(title='Select YAML File')
        fpath = Path(fname)
        if fname == '':
            pass
        elif fpath.is_file() and fpath.suffix in ['.yml', '.yaml']:
            self.ParamFileEntry_String.set(f'"{fpath}"')
            self.ProcessParameterFile(fname)
            self.ReadyIngridData()
        else:
            self.ParamFileEntry_String.set(f'Invalid file extension "{fpath.suffix}"')

    def ReadyIngridData(self):
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
        self.Ingrid.StartSetup()

    def AnalyzeTopology(self):
        IG = self.Ingrid
        self.ViewData()
        IG.AnalyzeTopology()

        if IG.settings['grid_settings']['num_xpt'] == 2:
            IG.PlotTopologyAnalysis()

    def CreatePatches(self):
        IG = self.Ingrid
        patch_data_available = False

        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
            IG.RefreshSettings()

        try:
            if IG.settings['patch_data']['use_file'] is True:
                if Path(IG.settings['patch_data']['file']).is_file():
                    patch_data_available = True
        except:
            pass

        if patch_data_available is True:
            IG.LoadPatches()
        else:
            self.AnalyzeTopology()
            IG.ConstructPatches()

        IG.PlotPatches()

    def CreateSubgrid(self):
        IG = self.Ingrid
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
            IG.RefreshSettings()
        IG.ConstructGrid()
        IG.PlotGrid()

    def ExportGridue(self):
        IG = self.Ingrid
        fname = filedialog.asksaveasfilename(initialdir='.', title='Save File', defaultextension='', initialfile='gridue')
        if fname != '':
            IG.ExportGridue(fname)

    def ViewData(self):
        IG = self.Ingrid
        self.ReadyIngridData()
        IG.ShowSetup(view_mode=IG.settings['grid_settings']['view_mode'])

    def Quit(self):
        self.controller.Exit()
