import matplotlib
import os

try:
    if os.environ.get('DISPLAY','') == '':
        print('no display found. Using non-interactive Agg backend')
        matplotlib.use('Agg')
    else:
        matplotlib.use("TkAgg")
except:
   active_backend = matplotlib.get_backend()
   msg = 'Warning: Could not set matplotlib backend to "TkAgg". '
   msg = f'The active backend is "{active_backend}"'
   print(msg)

import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.messagebox as tk_messagebox
import tkinter.filedialog as tk_filedialog
from time import time
from os.path import getmtime
from pathlib import Path
from INGRID.ingrid import Ingrid
from INGRID.exceptions import TkInitializationSuccess

class IngridGUI:
    """
    Class for supervising tkinter GUI and Ingrid interactions

    Parameters
    ----------
    IngridSession: INGRID obj, default=None
        An active INGRID session to connect to the GUI.
        If None, we assume the user would like a new session.

    """
        
    def __init__(self, IngridSession=None) -> None:
        #
        # Track the Ingrid session or init a new one
        #
        self.IngridSession = IngridSession if IngridSession else Ingrid()

        #
        # Starts a tk session
        #
        self.tk_session = tk.Tk()
        self.main_frame = tk.Frame(master=self.tk_session)
        self.main_frame.grid_rowconfigure(0, weight=1)
        self.main_frame.grid_columnconfigure(0, weight=1)

        #
        # Default GUI dimensions for each frame type
        #
        self.gui_dimensions = {
            FilePicker: "550x270"
        }

        #
        # Populate the GUI with frames
        #
        self.InitializeFrames()

        #
        # Automatic loading of input files using the Ingrid
        # session
        #
        if self.IngridSession.InputFile:
            self.frames[FilePicker].load_param_file(
                ExistingParamFile=self.IngridSession.yaml
            )

    def NewSettingsPrompt(self, message='Create a new settings file for this INGRID session?'):
        """
        Prompt user with dialog for new settings/config file

        Parameters
        ----------
        message: str, optional
            String to forward to the Yes/No dialog box
        """

        fname = ''

        while fname in ['', None]:
            #
            # Break out of box prompt if "no" selected
            #
            if not tk_messagebox.askyesno('', message):
                break

            #
            # Otherwise, proceed to settings dialog prompt
            #
            fname = 'INGRID_Session' + str(int(time())) + '.yml'
            fname = tk_filedialog.asksaveasfilename(
                        initialfile      = fname,
                        initialdir       = '.', 
                        title            = 'New YAML',
                        defaultextension = '.yml', 
                    )

            #
            # If we get an empty string or None for the filename,
            # we need to re-prompt the user until a valid name is
            # present, or they exit the dialog properly.
            #
            if fname in ['', None]:
                continue

            #
            # Save an Ingrid settings file when we have a valid filename
            #
            self.IngridSession.SaveSettingsFile(fname=fname)
            break

    def Run(self, test_initialization: bool = False):
        """
        Initiate the tk.mainloop() call

        Parameters
        ----------
        test_initialization : optional
            Flag to trigger a TkInitializationSuccess exception when
            entering the method calling tk.mainloop(). Suggests a successful
            initialization of the IngridGUI class on the tk side of things.
        """
        def on_closing():
            self.Exit()
        self.tk_session.title('INGRID')
        self.tk_session.protocol('WM_DELETE_WINDOW', on_closing)
        if test_initialization:
            self.tk_session.after(1, self.tk_session.destroy())
        self.tk_session.mainloop()
        #
        # If test_initialization flag is set, the gui mainloop
        # will auto-destroy after initialization. We raise this
        # exception to be caught in the test suite
        #
        if test_initialization:
            raise TkInitializationSuccess

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

        prompt_response = tk_messagebox.askyesno('', message)

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

    def Exit(self, message='Are you sure you want to quit?'):
        """
        Prompt the user with the Exit menu dialog

        Parameters
        ----------
        message: str, optional
            String to forward to message box
        """
        if tk_messagebox.askyesno('', message):
            try:
                plt.close('all')
            except:
                pass
            self.tk_session.destroy()

    def InitializeFrames(self):
        """
        Initialize the tkinter frames of interest for the GUI
        session.

        We currently only use the FilePicker frame class, but this
        pattern could support multiple "pages"/frames per GUI session.

        Returns
        -------
            A dict mapping Frame class to initialized Frame object
        """
        frames = {}
        for FrameType in [FilePicker]:
            frame = FrameType(controller=self, tk_session=self.tk_session)
            frames[FrameType] = frame
        self.frames = frames
        return frames

class FilePicker(tk.Frame):
    """
    FilePicker represents the tkinter Frame for selecting parameter
    files and running the workflow.

    Parameters
    ----------
    controller: IngridGUI
        The IngridGUI instance that is controlling the FilePicker.
        Allows us to reference the active Ingrid session from the Frame.
    tk_session: tk.Tk object
        The tk session the frame is controlled by.

    """
    def __init__(self, controller, tk_session):
        super().__init__(master=tk_session)
        self.controller  = controller
        self.tk_session  = tk_session
        self.PreviewPlot = None

        #
        # Initialize frames that drive the FilePicker
        #
        self.EntryFrame   = tk.Frame(master=tk_session)
        self.ControlFrame = tk.Frame(master=tk_session)

        #
        # Place the frames in appropriate locations
        #
        self.EntryFrame.grid(
            row    = 0, 
            column = 0, 
            padx   = 10, 
            pady   = 10
        )

        self.ControlFrame.grid(
            row    = 1, 
            column = 0, 
            padx   = 10, 
            pady   = 10
        )

        #
        # Initialize parameter file label settings
        #
        self.ParamFileLabel_String = tk.StringVar()
        self.ParamFileLabel_String.set('Parameter File Path:')
        self.ParamFileLabel = tk.Label(
            master = self.EntryFrame,
            text   = self.ParamFileLabel_String.get()
        )

        #
        # Initialize parameter file entry settings
        #
        self.ParamFileEntry_String = tk.StringVar()
        self.ParamFileEntry = tk.Entry(
            master             = self.EntryFrame,
            text               = self.ParamFileEntry_String,
            width              = 50,
            state              = 'disabled',
            disabledbackground = '#f8f8ff'
        )

        #
        # Initialize parameter file button settings
        #
        self.ParamFileButton = tk.Button(
            master  = self.EntryFrame,
            text    = 'Select Parameter File', 
            command = self.LoadParameterFile
        )

        #
        # Place the intialized items in the tk session
        #
        self.ParamFileLabel.grid(
            row    = 0, 
            column = 0, 
            sticky = 'nsew'
        )

        self.ParamFileEntry.grid(
            row    = 0, 
            column = 1, 
            sticky = 'nsew'
        )

        self.ParamFileButton.grid(
            row    = 0, 
            column = 2, 
            sticky = 'nsew'
        )

        #
        # Initialize the remaining control buttons for the FilePicker
        #
        self.ViewDataButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'View Loaded File',
            command = self.ViewData
        )
        self.AnalyzeTopologyButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'Analyze Topology', 
            command = self.AnalyzeTopology
        )
        self.CreatePatchesButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'Create Patches', 
            command = self.CreatePatches
        )
        self.CreateSubgridButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'Create Grid', 
            command = self.CreateSubgrid
        )
        self.ExportGridueButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'Export gridue', 
            command = self.ExportGridue
        )
        self.QuitButton = tk.Button(
            master  = self.ControlFrame,
            text    = 'Quit', 
            command = self.Quit
        )

        #
        # Place the buttons in the tk session
        #
        self.ViewDataButton.grid(
            row    = 0, 
            column = 0, 
            padx   = 10, 
            pady   = 10, 
            sticky = 'nsew'
        )
        self.CreatePatchesButton.grid(
            row    = 0, 
            column = 2, 
            padx   = 10, 
            pady   = 10, 
            sticky = 'nsew'
        )
        self.CreateSubgridButton.grid(
            row    = 0, 
            column = 3, 
            padx   = 10,
            pady   = 10, 
            sticky = 'nsew'
        )
        self.ExportGridueButton.grid(
            row    = 0, 
            column = 4, 
            padx   = 10, 
            pady   = 10, 
            sticky = 'nsew'
        )
        self.QuitButton.grid(
            row    = 0, 
            column = 5, 
            padx   = 10, 
            pady   = 10, 
            sticky = 'nsew'
        )

    def NewCase(self):
        """
        Helper method for initializing a new Ingrid session/case
        """
        self.controller.IngridSession = Ingrid()

    def ProcessParameterFile(self, fname):
        """
        Helper method for processing an Ingrid parameter/config file
        
        Parameters
        ----------
        fname: str
            Absolute path to the Ingrid parameter file to process
        """
        self.NewCase()
        self.controller.IngridSession.InputFile = fname
        self.controller.IngridSession.PopulateSettings(Ingrid.ReadYamlFile(fname))
        self.ParamFileMtime = getmtime(fname)
        self.ParamFileName  = fname

    def LoadParameterFile(self):
        """
        Helper method for loading an Ingrid parameter/config file from
        the user system
        """
        fname = tk_filedialog.askopenfilename(title='Select YAML File')
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
        """
        Helper method for booting up Ingrid. 
        
        Processes data in parameter/config file and initializes an Ingrid session
        if appropriate.
        """
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
        self.controller.IngridSession.StartSetup()

    def AnalyzeTopology(self):
        """
        Method for analyzing data topology. 

        To be used by tkinter control button and interfaces with 
        Ingrid engine.
        """
        self.ViewData()
        self.controller.IngridSession.AnalyzeTopology()

        if self.controller.IngridSession.settings['grid_settings']['num_xpt'] == 2:
            self.controller.IngridSession.PlotTopologyAnalysis()

    def CreatePatches(self):
        """
        Method for creating Ingrid Patches. 

        To be used by tkinter control button and interfaces with 
        Ingrid engine.
        """
        patch_data_available = False

        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
            self.controller.IngridSession.RefreshSettings()

        try:
            if self.controller.IngridSession.settings['patch_data']['use_file'] is True:
                if Path(self.controller.IngridSession.settings['patch_data']['file']).is_file():
                    patch_data_available = True
        except:
            pass

        if patch_data_available is True:
            self.controller.IngridSession.LoadPatches()
        else:
            self.AnalyzeTopology()
            self.controller.IngridSession.ConstructPatches()

        self.controller.IngridSession.PlotPatches()

    def CreateSubgrid(self):
        """
        Method for creating Ingrid local meshes. 

        To be used by tkinter control button and interfaces with 
        Ingrid engine.
        """
        if self.ParamFileMtime != getmtime(self.ParamFileName):
            self.ProcessParameterFile(self.ParamFileName)
            self.controller.IngridSession.RefreshSettings()
        self.controller.IngridSession.ConstructGrid()
        self.controller.IngridSession.PlotGrid()

    def ExportGridue(self):
        """
        Method for exporting gridue files. 

        To be used by tkinter control button and interfaces with 
        Ingrid engine.
        """
        fname = tk_filedialog.asksaveasfilename(
                    initialdir       = '.', 
                    title            = 'Save File', 
                    defaultextension = '', 
                    initialfile      = 'gridue'
                )

        if fname != '':
            self.controller.IngridSession.ExportGridue(fname)

    def ViewData(self):
        """
        Method for visualizing loaded data.

        To be used by tkinter control button and interfaces with 
        Ingrid engine.
        """
        self.ReadyIngridData()
        view_mode = self.controller.IngridSession.settings['grid_settings']['view_mode']
        self.controller.IngridSession.ShowSetup(view_mode=view_mode)

    def Quit(self):
        """
        Helper method for exiting the GUI
        """
        self.controller.Exit()
