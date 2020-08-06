"""
SNL.py

Description:
    SNL configuration class.

Created: June 18, 2020

"""
import numpy as np
import matplotlib
import pathlib
import inspect
import yaml as _yaml_
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry import Point, Line, Patch, segment_intersect, angle_between, rotate, find_split_index


class SNL():
    """
    The SNL (Single-Null) class is the parent class for both upper-single null (USN)
    and lower-single null (LSN) configurations.
    This base class handles the formatting and plotting of data obtained from an LSN or USN
    object.
    Parameter:
        - INGRID_object : Ingrid class object
        All SNL objects are children of the main Ingrid class. INGRID_object provides
        information such as YAML data, efit_psi, psi_norm, and plate_data.
    @author: garcia299
    """

    def __init__(self, Ingrid_obj, config):

        self.parent = Ingrid_obj
        self.config = config
        self.settings = Ingrid_obj.settings
        self.plate_data = Ingrid_obj.plate_data

        self.parent.OrderTargetPlates()
        self.PatchTagMap = self.parent.GetPatchTagMap(config='SNL')

        self.eq = Ingrid_obj.eq
        self.efit_psi = Ingrid_obj.efit_psi
        self.psi_norm = Ingrid_obj.psi_norm
        self.eq = Ingrid_obj.eq
        self.CurrentListPatch={}
        self.Verbose=False
        self.plate_data = Ingrid_obj.plate_data
        self.CorrectDistortion={}

    def grid_diagram(self,ax=None):
        """
        Create Grid matplotlib figure for an SNL object.
        @author: watkins35, garcia299
        """
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
          'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
          'seagreen', 'firebrick', 'saddlebrown']

        fig=plt.figure('INGRID: Grid', figsize=(6,10))
        if ax is None:
            ax=fig.gca()
        for patch in self.patches:
            patch.plot_subgrid(ax)
            print('patch completed...')

        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('INGRID SNL Subgrid')
        plt.show()

    def patch_diagram(self):
        """
        Generates the patch diagram for a given configuration.
        @author: watkins35, garcia299
        """

        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']


        self.FigPatch=plt.figure('INGRID: Patch Map', figsize=(6, 10))
        self.FigPatch.clf()
        ax=self.FigPatch.subplots(1,1)
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        ax.set_aspect('equal', adjustable='box')

        ax.set_xlabel('R')
        ax.set_ylabel('Z')
        self.FigPatch.suptitle('SNL Patch Diagram')

        for i, patch in enumerate(self.patches.values()):
            patch.plot_border('green')
            patch.fill(colors[i])
            patch.color=colors[i]
        # ax.legend()
        plt.show()

    def get_config(self):
        """
        Returns a string indicating whether the
        """
        return self.config

    def concat_grid(self, Verbose=False):
        """
        Concatenate all local grids on individual patches into a single
        array with branch cuts
        Parameters:
        ----------
            Verbose : bool
                Verbose flag indicator
        """
        # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).
        patch_matrix = self.patch_matrix

        # Get some poloidal and radial information from each patch to attribute to the
        # local subgrid.
        # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
        #       the value of 1 to the cell number to get the accurate number of lines.


        for patch in self.patches.values():
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
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm)):
            rm[i] = rm[i][::-1]
        for i in range(len(zm)):
            zm[i] = zm[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb = len(rm) - 2
        ixlb = 0
        self.rm = self.add_guardc(rm, ixlb, ixrb)
        self.zm = self.add_guardc(zm, ixlb, ixrb)

        try:
            debug = self.settings['DEBUG']['visual']['gridue']
        except:
            debug = False

        if debug:
            self.animate_grid()


    def add_guardc(self, cell_map, ixlb, ixrb, nxpt = 1, eps = 1e-3):

        def set_guard(cell_map, ix, iy, eps, boundary):
            # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
            # TODO: Edit the code to reflect this at some point so the next reader is not overwhelmed.
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
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'left')
            ix = ixrb + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'right')

        for ix in range(np + 2):
            iy = 0
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'bottom')
            iy = nr + 1
            cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'top')

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

    def get_func(self, _func):

            def make_sympy_func(var, expression):
                import sympy as sp
                _f = sp.lambdify(var, expression, 'numpy')
                return _f

            f_str_raw = _func

            f_str_raw = f_str_raw.replace(' ', '')
            delim = f_str_raw.index(',')

            var = f_str_raw[0 : delim]
            expression = f_str_raw[delim + 1 :]

            _func = make_sympy_func(var, expression)
            # TODO: Check Normalization of the function to 1
            return _func

    def CheckFunction(self,Str,Verbose=False):
        # can use a fancier test of tehe lambda function like in get_func above
        try:
            com='lambda {}:{}'.format(Str.split(',')[0],Str.split(',')[1])
            if Verbose: print(com)
            eval(com)
            ExpValid=True
        except:
            raise ValueError('Unable to parse function {} for entry {}. Lambda function expected:"x,f(x)"'.format(self.val,Name))
            ExpValid=False
        return ExpValid

    def GetFunctions(self,Patch,Enforce=True,Verbose=False,ExtraSettings={}):
        if Patch in self.SOL:
            try:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_sol']
            except:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_primary_sol']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('SOL radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.CORE:
            _radial_f = self.settings['grid_params']['grid_generation']['radial_f_core']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('Core radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.PF:
            try:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_pf']
            except:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_primary_pf']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('Core radial transformation: "{}"'.format(_radial_f))
        if valid_function:
            _radial_f = self.get_func( _radial_f)
        else:
            _radial_f=lambda x:x

        if Patch in self.SOL:
            _poloidal_f = self.settings['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('SOL poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.CORE:
            _poloidal_f = self.settings['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('Core poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.PF:
            _poloidal_f = self.settings['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('Core poloidal transformation: "{}"'.format(_poloidal_f))
        if valid_function:
            _poloidal_f = self.get_func( _poloidal_f)
        else:
            _poloidal_f=lambda x:x


        if Patch.platePatch:
            if self.config == 'LSN':
                if Patch.plateLocation == 'W':
                    _poloidal_f = self.settings['target_plates']['plate_E1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
                elif Patch.plateLocation == 'E':
                    _poloidal_f = self.settings['target_plates']['plate_W1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if self.config == 'USN':
                if Patch.plateLocation == 'E':
                    _poloidal_f = self.settings['target_plates']['plate_E1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
                elif Patch.plateLocation == 'W':
                    _poloidal_f = self.settings['target_plates']['plate_W1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if valid_function:
                _poloidal_f = self.get_func( _poloidal_f)
            else:
                _poloidal_f = lambda x: x
        print('Extra Settings:',ExtraSettings)
        
        if ExtraSettings.get(Patch.patchName) is not None:
            Extra=ExtraSettings.get(Patch.patchName)
            print('Extra:',Extra)
            if Extra.get('poloidal_f') is not None:
                if self.CheckFunction(Extra.get('poloidal_f'),Enforce):
                    print('Overriding poloidal function for patch {} with f={}'.format(Patch.patchName,Extra.get('poloidal_f')))
                    _poloidal_f=self.get_func(Extra.get('poloidal_f'))
            if Extra.get('radial_f') is not None:
                if self.CheckFunction(Extra.get('radial_f'),Enforce):
                    print('Overriding radial function for patch {} with f={}'.format(Patch.patchName,Extra.get('radial_f')))
                    _radial_f=self.get_func(Extra.get('radial_f'))
            
                 
        

        return (_radial_f,_poloidal_f)


    def GetNpoints(self,Patch,Enforce=True,Verbose=False):
        try:
            np_sol = self.settings['grid_params']['grid_generation']['np_sol']
        except:
            np_sol = self.settings['grid_params']['grid_generation']['np_global']
        try:
            np_core = self.settings['grid_params']['grid_generation']['np_core']
        except:
            np_core = self.settings['grid_params']['grid_generation']['np_global']
        try:
            np_pf = self.settings['grid_params']['grid_generation']['np_pf']
        except:
            np_pf = self.settings['grid_params']['grid_generation']['np_global']

        if np_sol != np_core:
            if Enforce:
                raise ValueError('SOL and CORE must have equal POLOIDAL np values')
            else:
                print('WARNING: SOL and CORE must have equal POLOIDAL np values!\nSetting np values' \
                + ' to the minimum of np_sol={} and np_core={}.\n'.format(np_sol,np_core))
            np_sol = np_core = np.amin([np_sol, np_core])

        try:
            nr_sol = self.settings['grid_params']['grid_generation']['nr_sol']
        except:
            nr_sol = self.settings['grid_params']['grid_generation']['nr_global']

        try:
            nr_core = self.settings['grid_params']['grid_generation']['nr_core']
        except:
            nr_core = self.settings['grid_params']['grid_generation']['nr_global']
        try:
            nr_pf = self.settings['grid_params']['grid_generation']['nr_pf']
        except:
            nr_pf = self.settings['grid_params']['grid_generation']['nr_global']

        if nr_pf != nr_core:
            if Enforce:
                raise ValueError('PF and CORE must have equal RADIAL nr values')
            else:
                print('WARNING: PF and CORE must have equal RADIAL nr values!\nSetting nr values' \
                + ' to the minimum of nr_pf={} and nr_core={}.\n'.format(nr_pf,nr_core))
                nr_pf = nr_core = np.amin([nr_pf, nr_core])

        if Verbose:
            print("nr_sol={}; nr_pf={};nr_core={}".format(nr_sol,nr_pf,nr_core))
            print("np_sol={}; np_pf={};np_core={}".format(np_sol,np_pf,np_core))

        if Patch in self.SOL:
            nr_cells = nr_sol
            np_cells = np_sol
            if Verbose: print('Patch "{}" is in SOL'.format(patch.patchName))
        elif Patch in self.CORE:
            nr_cells = nr_core
            np_cells = np_core
            if Verbose: print('Patch "{}" is in CORE'.format(patch.patchName))
        elif Patch in self.PF:
            nr_cells = nr_pf
            np_cells = np_pf
            if Verbose: print('Patch "{}" is in PF with nr={} and np={}'.format(patch.patchName,np_cells,nr_cells))

        # overwrite poloidal values for target patches
        if Patch.platePatch:
            if self.config == 'LSN':
                if Patch.plateLocation == 'W':
                    np_cells = self.settings['target_plates']['plate_W1']['np_local']
                elif Patch.plateLocation == 'E':
                    np_cells = self.settings['target_plates']['plate_E1']['np_local']
            if self.config == 'USN':
                if Patch.plateLocation == 'E':
                    np_cells = self.settings['target_plates']['plate_W1']['np_local']
                elif Patch.plateLocation == 'W':
                    np_cells = self.settings['target_plates']['plate_E1']['np_local']

        return (nr_cells,np_cells)

    def AdjustPatch(self,patch):
        primary_xpt = Point(self.eq.NSEW_lookup['xpt1']['coor']['center'])

        tag = patch.get_tag()
        if tag == 'A2':
            patch.adjust_corner(primary_xpt, 'SE')
        elif tag == 'A1':
            patch.adjust_corner(primary_xpt, 'NE')
        elif tag == 'B2':
            patch.adjust_corner(primary_xpt, 'SW')
        elif tag == 'B1':
            patch.adjust_corner(primary_xpt, 'NW')
        elif tag == 'E1':
            patch.adjust_corner(primary_xpt, 'NE')
        elif tag == 'E2':
            patch.adjust_corner(primary_xpt, 'SE')
        elif tag == 'F1':
            patch.adjust_corner(primary_xpt, 'NW')
        elif tag == 'F2':
            patch.adjust_corner(primary_xpt, 'SW')
    
    def AdjustPatches(self):
        for patch in self.patches.values():
            self.AdjustPatch(patch)
                
    def CheckPatches(self,verbose=False):
        for name, patch in self.patches.items():
            if patch.platePatch:
                print(' # Checking patch: ', name)
                patch.CheckPatch(self)

    def SavePatches(self,FileName):
        with open(FileName,'w') as File:
            _yaml_.dump(self.patches,File)

    def LoadPatches(self,FileName):
        with open(FileName,'r') as File:
            self.patches = _yaml_.load(File,Loader=_yaml_.Loader)
        print(self.patches)
        if isinstance(self.patches, list):
            patch_dict = {}
            for patch in self.patches:
                patch_dict[patch.patchName] = patch
            self.patches = patch_dict
        for patch in self.patches.values():
            patch.plot_border()
            patch.fill()
        p = self.patches
        print('>>> Loaded Patches:', [k for k in self.patches.keys()])

        self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                    [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
                    [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
                    [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                    ]
        self.categorize_patches()

    def construct_grid(self, np_cells = 1, nr_cells = 1,Verbose=False,ShowVertices=False,RestartScratch=False,OptionTrace='theta',ExtraSettings={},ListPatches='all', Enforce=True):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.
        if Verbose: print('Construct Grid')
        try:
            visual = self.settings['DEBUG']['visual']['subgrid']
        except:
            visual = False
        try:
            verbose = self.settings['DEBUG']['verbose']['grid_generation']
        except:
            verbose = False

        verbose=Verbose or verbose
        
            
        print('>>> Patches:', [k for k in self.patches.keys()])
        if RestartScratch:
            self.CurrentListPatch={}
    
        for name, patch in self.patches.items():
            
            if self.CorrectDistortion.get(name) is not None:
               patch.CorrectDistortion=self.CorrectDistortion.get(name)
            elif self.CorrectDistortion.get('all') is not None:
                patch.CorrectDistortion=self.CorrectDistortion.get('all')
            else:
                patch.CorrectDistortion={'Active':False}
            if (ListPatches=='all' and patch not in self.CurrentListPatch) or (ListPatches!='all' and name in ListPatches):
                self.SetPatchBoundaryPoints(patch)
                (nr_cells,np_cells)=self.GetNpoints(patch, Enforce=Enforce)
                (_radial_f,_poloidal_f)=self.GetFunctions(patch,ExtraSettings=ExtraSettings,Enforce=Enforce)
                print('>>> Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(name, np_cells, nr_cells, inspect.getsource(_poloidal_f), inspect.getsource(_radial_f)))
                patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f,verbose = verbose, visual = visual,ShowVertices=ShowVertices,OptionTrace=OptionTrace)
                self.AdjustPatch(patch)
                patch.plot_subgrid()
                self.CurrentListPatch[name] = patch



        if all(['cell_grid' in patch.__dict__ for patch in self.patches.values()]):
            self.concat_grid()
            self.set_gridue()
        
    def SetPatchBoundaryPoints(self,Patch):
            if self.ConnexionMap.get(Patch.get_tag()) is not None:
                if self.Verbose: print('Find connexion map for patch {}'.format(Patch.patchName))
                for Boundary,AdjacentPatch in self.ConnexionMap.get(Patch.get_tag()).items():
                    Patch.BoundaryPoints[Boundary]=self.GetBoundaryPoints(AdjacentPatch)
                    if self.Verbose: print('Find Boundaries points for {}'.format(Patch.patchName))
                
    def GetBoundaryPoints(self,AdjacentPatch):
        if AdjacentPatch is not None:
            PatchName=AdjacentPatch[0]
            Boundary=AdjacentPatch[1]
            for name, patch in self.patches.items():
                if name == PatchName:
                   if Boundary=='S':
                       return patch.S_vertices
                   elif Boundary=='N':
                       return patch.N_vertices
                   elif Boundary=='E':
                       return patch.W_vertices
                   elif Boundary=='W':
                       return patch.W_vertices 
        return None

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

        try:
            visual = self.settings['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False
        try:
            verbose = self.settings['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False
        try:
            inner_tilt = self.settings['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0
        try:
            outer_tilt = self.settings['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0


        WestPlate = self.plate_data['plate_W1']
        EastPlate = self.plate_data['plate_E1']

        xpt = self.eq.NSEW_lookup['xpt1']['coor']
        magx = np.array([self.settings['grid_params']['rmagx'] + self.settings['grid_params']['patch_generation']['rmagx_shift'], \
            self.settings['grid_params']['zmagx'] + self.settings['grid_params']['patch_generation']['zmagx_shift']])

        psi_max = self.settings['grid_params']['psi_max']
        psi_min_core = self.settings['grid_params']['psi_min_core']
        psi_min_pf = self.settings['grid_params']['psi_min_pf']

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
        if self.settings['grid_params']['patch_generation']['use_NW']:
            tilt = self.settings['grid_params']['patch_generation']['NW_adjust']
            xptW_psiMax = self.eq.draw_line(rotate(xpt['W'], tilt, xpt['center']), {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
        else:
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        if self.settings['grid_params']['patch_generation']['use_NE']:
            tilt = self.settings['grid_params']['patch_generation']['NE_adjust']
            xptE_psiMax = self.eq.draw_line(rotate(xpt['E'], tilt, xpt['center']), {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        xpt_WestPlate = self.eq.draw_line(xpt['SW'], {'line' : WestPlate}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        xptS_psiMinPF = self.eq.draw_line(xpt['S'], {'psi' : psi_min_pf}, option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        xpt_EastPlate = self.eq.draw_line(xpt['SE'], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        iPsiMax_TP = self.eq.draw_line(xptW_psiMax.p[-1], {'line' : WestPlate}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        psiMinPF_WestPlate = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : WestPlate},option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        oPsiMax_TP = self.eq.draw_line(xptE_psiMax.p[-1], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        psiMinPF_EastPlate = self.eq.draw_line(xptS_psiMinPF.p[-1], {'line' : EastPlate}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

        imidLine_topLine = self.eq.draw_line(xptNW_midLine.p[-1], {'line' : topLine}, option = 'theta', \
            direction = 'cw', show_plot = visual, text = verbose)
        
        omidLine_topLine = self.eq.draw_line(xptNE_midLine.p[-1], {'line' : topLine}, option = 'theta', \
            direction = 'ccw', show_plot = visual, text = verbose)

        # Integrating horizontally along mid-line towards psiMax and psiMinCore

        imidLine_psiMax = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_max, inner_tilt)}, option = 'z_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)
        imidLine_psiMinCore = self.eq.draw_line(xptNW_midLine.p[-1], {'psi_horizontal' : (psi_min_core, inner_tilt)}, option = 'z_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        omidLine_psiMax = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_max, outer_tilt)}, option = 'z_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        omidLine_psiMinCore = self.eq.draw_line(xptNE_midLine.p[-1], {'psi_horizontal' : (psi_min_core, outer_tilt)}, option = 'z_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)

        # Integrating vertically along top-line towards psiMax and psiMinCore
        topLine_psiMax = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_max}, option = 'r_const', \
                direction = 'cw' if self.config == 'LSN' else 'ccw', show_plot = visual, text = verbose)
        topLine_psiMinCore = self.eq.draw_line(omidLine_topLine.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', \
                direction = 'ccw' if self.config == 'LSN' else 'cw', show_plot = visual, text = verbose)

        # A1 Patch
        location = 'W'
        A2_N = iPsiMax_TP.reverse_copy()
        A2_S = xpt_WestPlate
        A2_E = xptW_psiMax.reverse_copy()
        # =====================================================================================
        # Trimming the target_plate to conform to the patch boundary.
        # -------------------------------------------------------------------------------------
        # Recall ITP has a clockwise orientation.
        #
        # The inner 'split' trims all Point objects BEFORE the point of intersection of ITP
        # and A2_S. Call this new Line object Line_A.
        #
        # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
        # and A2_N. This new Line object is the plate facing boundary of the Patch.
        # =====================================================================================
        A2_W = (WestPlate.split(A2_S.p[-1])[1]).split(A2_N.p[0], add_split_point = True)[0]
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # A1 Patch
        location = 'W'
        A1_N = A2_S.reverse_copy()
        
        A1_S = psiMinPF_WestPlate
        A1_E = xptS_psiMinPF
        A1_W = (WestPlate.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # B2 Patch
        
        B2_N = self.eq.draw_line(A2_N.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B2_S = xptNW_midLine.reverse_copy()
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B2_W = xptW_psiMax
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'ISB')

        # B1 Patch
        B1_N = B2_S.reverse_copy()
        B1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B1_W = xptN_psiMinCore.reverse_copy()
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'ICB')

        # C2 Patch
        C2_N = self.eq.draw_line(B2_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C2_S = imidLine_topLine.reverse_copy()
        C2_E = Line([C2_N.p[-1], C2_S.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'IST')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.eq.draw_line(B1_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        C1_E = Line([C1_N.p[-1], C1_S.p[0]])
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'ICT')

        # F2 Patch
        location = 'E'
        F2_N = oPsiMax_TP
        F2_S = xpt_EastPlate.reverse_copy()
        F2_E = (EastPlate.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point = True)[0]
        F2_W = xptE_psiMax
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # F1 Patch
        location = 'E'
        F1_N = F2_S.reverse_copy()
        F1_S = psiMinPF_EastPlate.reverse_copy()
        F1_E = (EastPlate.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point = True)[0]
        F1_W = xptS_psiMinPF.reverse_copy()
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # E2 Patch
        E2_N = self.eq.draw_line(F2_N.p[0], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        E2_S = xptNE_midLine
        E2_E = xptE_psiMax.reverse_copy()
        E2_W = Line([E2_S.p[-1], E2_N.p[0]])
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'OSB')

        # E1 Patch
        E1_N = E2_S.reverse_copy()
        E1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E1_E = xptN_psiMinCore
        E1_W = Line([E1_S.p[-1], E1_N.p[0]])
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'OCB')

        # D2 Patch
        D2_N = self.eq.draw_line(E2_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D2_S = omidLine_topLine
        D2_E = Line([D2_N.p[-1], D2_S.p[0]])
        D2_W = Line([D2_S.p[-1], D2_N.p[0]])
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'OST')

        # D1 Patch
        D1_N = D2_S.reverse_copy()
        D1_S = self.eq.draw_line(E1_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_E = Line([D1_N.p[-1], D1_S.p[0]])
        D1_W = Line([D1_S.p[-1], D1_N.p[0]])
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'OCT')

        patches = [A2, A1, B2, B1, C2, C1, D2, D1, E2, E1, F2, F1]

        self.patches = {}
        for patch in patches:
            patch.PatchTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch

        p = self.patches
        self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                        [[None], p[A2.patchName], p[B2.patchName], p[C2.patchName], p[D2.patchName], p[E2.patchName], p[F2.patchName], [None]], \
                        [[None], p[A1.patchName], p[B1.patchName], p[C1.patchName], p[D1.patchName], p[E1.patchName], p[F1.patchName], [None]], \
                        [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                        ]

        self.categorize_patches()

    def categorize_patches(self):
        m = self.patch_matrix
        self.SOL = m[1][1:-1]
        self.CORE = m[2][2:-2]
        self.PF = [m[2][1], m[2][-2]]

    def set_gridue(self):
        """
        Prepare the relevant arrays for writing to GRIDUE.
        """

        # RECALL: self.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(self.rm) - 2
        ixpt1 = self.patches['IDL'].npol - 1
        ixpt2 = ixrb - self.patches['ODL'].npol + 1
        iyseparatrix1 = self.patches['IDL'].nrad - 1
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