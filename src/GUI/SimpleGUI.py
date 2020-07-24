from __future__ import print_function

from sys import platform as sys_pf
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
import Ingrid
import Root_Finder

#helv_large = 'Helvetica 13 bold'
#helv_medium = 'Helvetica 9 bold'

class Ingrid_GUI(tk.Tk):
    def __init__(self, master = None,IngridSession=None,*args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand = True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        self.gui_dimensions = {FilePicker : "550x270"}

        # attached to the parent Ingrid instance instead of new instance
        if IngridSession is None:
            self.ResetIngrid()
        else:
            self.Ingrid = IngridSession
        
        self.PopulateGUI()

        #automatic loading of input files
        if self.Ingrid.InputFile is not None:
            self.frames[FilePicker].load_param_file(
                ExistingParamFile=self.Ingrid.yaml)

    def ShowFrame(self, item):
        frame = self.frames[item]
        frame.tkraise()
        self.geometry(self.gui_dimensions[FilePicker])


    def Reset(self, message = 'Are you sure you want to reset?'):
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
        self.Ingrid = Ingrid.Ingrid()

    def ResetIG(self):
        self.Ingrid = None

    def PopulateGUI(self):
        self.frames = {}

        for F in [FilePicker]:
            frame = F(self.container, self)
            frame.Ingrid=self.Ingrid
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

        self.EntryFrame.grid(row=0,column=0,padx=10,pady=10)
        self.ControlFrame.grid(row=1,column=0,padx=10,pady=10)

        self.ParamFileLabel_String = tk.StringVar()
        self.ParamFileEntry_String = tk.StringVar()
        self.ParamFileLabel_String.set('Parameter File Path:')
        
        self.ParamFileLabel = tk.Label(self.EntryFrame,
            text=self.ParamFileLabel_String.get())
        self.ParamFileEntry = tk.Entry(self.EntryFrame, 
            text=self.ParamFileEntry_String, width=40,
            disabledbackground='#f8f8ff', state='disabled')
        self.ParamFileButton = tk.Button(self.EntryFrame,
            text='Select Parameter File', command=self.LoadParameterFile)

        self.ParamFileLabel.grid(row=0,column=0,sticky='nsew')
        self.ParamFileEntry.grid(row=0,column=1,sticky='nsew')
        self.ParamFileButton.grid(row=0,column=2,sticky='nsew')


        self.ViewDataButton = tk.Button(self.ControlFrame,
            text='View Loaded File', command=self.ViewData)
        self.AnalyzeTopologyButton = tk.Button(self.ControlFrame,
            text='Analyze Topology', command=self.AnalyzeTopology)
        self.CreatePatchesButton = tk.Button(self.ControlFrame,
            text='Create Patches', command=self.CreatePatches)
        self.QuitButton = tk.Button(self.ControlFrame,
            text='Quit', command=self.Quit)

        self.ViewDataButton.grid(row=0,column=0,padx=10,pady=10,sticky='nsew')
        self.AnalyzeTopologyButton.grid(row=0,column=1,padx=10,pady=10,sticky='nsew')
        self.CreatePatchesButton.grid(row=0,column=2,padx=10,pady=10,sticky='nsew')
        self.QuitButton.grid(row=0,column=3,padx=10,pady=10,sticky='nsew')

    def ProcessParameterFile(self, fname):
        self.Ingrid.process_yaml(Ingrid.ReadyamlFile(fname))
        self.ParamFileMtime = getmtime(fname)
        self.ParamFileName = fname

    def LoadParameterFile(self):
        fname = filedialog.askopenfilename(title='Select YAML File')
        fpath = Path(fname)
        if fpath.is_file() and fpath.suffix in ['.yml', '.yaml']:
            self.controller.NewIG()
            self.ProcessParameterFile(fname)
            self.ReadyIngridData()
        else:
            self.ParamFileEntry_String.set('Invalid file extension "{fpath.suffix}"')

    def ReadyIngridData(self):
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
        IG = self.Ingrid
        topology = 'SNL' if IG.yaml['grid_params']['num_xpt'] == 1 else 'DNL'
        IG.OMFIT_read_psi()
        IG.calc_efit_derivs()
        IG.AutoRefineMagAxis()
        IG.AutoRefineXPoint()
        if topology == 'DNL':
            IG.AutoRefineXPoint2()
        IG.read_target_plates()
        IG.set_limiter()
        IG.SetMagReference(topology)
        IG.calc_psinorm()

    def CloseOpenViews(self):
        try:
            plt.close('all')
        except:
            pass

    def PreviewIngridData(self):
        self.CloseOpenViews()
        IG = self.Ingrid
        IG.plot_psinorm()
        IG.PlotPsiNormMagReference()
        IG.plot_strike_geometry()
        IG.PlotPsiNormBounds()
        IG.PrintSummaryParams()

    def AnalyzeTopology(self):
        IG = self.Ingrid
        self.ViewData()
        IG.AnalyzeTopology()

        if IG.current_topology.config in ['SF15', 'SF45']:
            IG.PlotTopologyAnalysis()

    def CreatePatches(self):
        IG = self.Ingrid
        IG.ConstructPatches()

    def ViewData(self):
        self.ReadyIngridData()
        self.PreviewIngridData()

    def Quit(self):
        self.controller.Exit()











