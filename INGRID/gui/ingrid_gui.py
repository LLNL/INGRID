import matplotlib
try:
    matplotlib.use("TkAgg")
except:
    pass
try:
    import tkinter as tk
except:
    import Tkinter as tk

import matplotlib.pyplot as plt
import tkinter.filedialog as filedialog
from tkinter import messagebox as tkMB
from INGRID.ingrid import Ingrid
from time import time
from os.path import getmtime
from pathlib import Path

class IngridGUI:
        
    def __init__(self, IngridSession=None) -> None:

        #
        # Track the Ingrid session or init a new one
        #
        self.IngridSession = IngridSession if IngridSession else Ingrid()

        #
        # Starts a tk session
        #
        self.tk_session = tk.Tk()
        self.containter = tk.Frame(master=self.tk_session)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)

        #
        # Default GUI dimensions for each frame type
        #
        self.gui_dimensions = {
            FilePicker: "550x270"
        }

        #
        # Populate the GUI with frames
        #
        self.PopulateGUI()

        #
        # Automatic loading of input files using the Ingrid
        # session
        #
        if self.IngridSession.InputFile:
            self.frames[FilePicker].load_param_file(
                ExistingParamFile=self.IngridSession.yaml
            )

    def NewSettingsPrompt(self):
        fname = ''
        while fname in ['', None]:
            if tkMB.askyesno('', 'Create a new settings file for this INGRID session?'):
                fname = 'INGRID_Session' + str(int(time())) + '.yml'
                fname = filedialog.asksaveasfilename(initialdir='.', title='New YAML',
                    defaultextension='.yml', initialfile=fname)
                if fname not in ['', None]:
                    self.IngridSession.SaveSettingsFile(fname=fname)
                    break
            else:
                break

    def Reset(self, message='Are you sure you want to reset?'):
        """
        Reset the GUI and Ingrid state

        Parameters
        ----------
        message: str, optional
            Box message to prompt users with

        Returns
        -------
            Bool indicating prompt response
        """

        prompt_response = tkMB.askyesno('', message)

        if prompt_response:
            #
            # Close all open matplotlib figures
            #
            try:
                plt.close('all')
            except:
                pass

            #
            # Initialize the GUI frames on hand
            #
            self.InitializeFrames()

            #
            # Reset the Ingrid session state to None
            #
            self.IngridSession = None
            #
            # Reset the frame view to the FilePicker
            #
            frame = self.frames[FilePicker]
            frame.tkraise()
            self.tk_session.geometry(self.gui_dimensions[FilePicker])
            #
            # Set FilePicker gui dimensions
            #
            self.geometry(self.gui_dimensions[FilePicker])

        return prompt_response

    def Exit(self):
        if tkMB.askyesno('', 'Are you sure you want to quit?'):
            try:
                plt.close('all')
            except:
                pass
            self.destroy()

    def InitializeFrames(self):
        frames = {}
        for FrameType in [FilePicker]:
            frame = FrameType(master=self.tk_session)
            frame.IngridSession = self.IngridSession
            frames[FrameType] = frame
        self.frames = frames
        return frames

class FilePicker(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.master      = master
        self.PreviewPlot = None

        self.EntryFrame   = tk.Frame(master=master)
        self.ControlFrame = tk.Frame(master=master)

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
        self.controller.Ingrid = Ingrid()
        self.IngridSession = self.controller.Ingrid

    def ProcessParameterFile(self, fname):
        self.NewCase()
        self.IngridSession.InputFile = fname
        self.IngridSession.PopulateSettings(Ingrid.ReadYamlFile(fname))
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
        self.IngridSession.StartSetup()

    def AnalyzeTopology(self):
        IG = self.IngridSession
        self.ViewData()
        IG.AnalyzeTopology()

        if IG.settings['grid_settings']['num_xpt'] == 2:
            IG.PlotTopologyAnalysis()

    def CreatePatches(self):
        IG = self.IngridSession
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
        IG = self.IngridSession
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
            IG.RefreshSettings()
        IG.ConstructGrid()
        IG.PlotGrid()

    def ExportGridue(self):
        IG = self.IngridSession
        fname = filedialog.asksaveasfilename(initialdir='.', title='Save File', defaultextension='', initialfile='gridue')
        if fname != '':
            IG.ExportGridue(fname)

    def ViewData(self):
        IG = self.IngridSession
        self.ReadyIngridData()
        IG.ShowSetup(view_mode=IG.settings['grid_settings']['view_mode'])

    def Quit(self):
        self.master.Exit()
