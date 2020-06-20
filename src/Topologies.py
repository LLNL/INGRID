"""
Topology.py

Description:
    Collection of available Ingrid patch map topologies.

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
from geometry import Point, Line, SNL_Patch, DNL_Patch, segment_intersect


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
        self.yaml = Ingrid_obj.yaml
        self.plate_data = Ingrid_obj.plate_data

        self.parent.order_target_plates()

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

        #try:
            #plt.close('INGRID: Grid')
        #except:
            #pass
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
        ax.legend()
        plt.show()

    def get_configuration(self):
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
            debug = self.yaml['DEBUG']['visual']['gridue']
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
            _radial_f = self.yaml['grid_params']['grid_generation']['radial_f_sol']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('SOL radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.CORE:
            _radial_f = self.yaml['grid_params']['grid_generation']['radial_f_core']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('Core radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.PF:
            _radial_f = self.yaml['grid_params']['grid_generation']['radial_f_pf']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('Core radial transformation: "{}"'.format(_radial_f))
        if valid_function:
            _radial_f = self.get_func( _radial_f)
        else:
            _radial_f=lambda x:x

        if Patch in self.SOL:
            _poloidal_f = self.yaml['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('SOL poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.CORE:
            _poloidal_f = self.yaml['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('Core poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.PF:
            _poloidal_f = self.yaml['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('Core poloidal transformation: "{}"'.format(_poloidal_f))
        if valid_function:
            _poloidal_f = self.get_func( _poloidal_f)
        else:
            _poloidal_f=lambda x:x


        if Patch.platePatch:
            if self.config == 'LSN':
                if Patch.plateLocation == 'W':
                    _poloidal_f = self.yaml['target_plates']['plate_E1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
                elif Patch.plateLocation == 'E':
                    _poloidal_f = self.yaml['target_plates']['plate_W1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if self.config == 'USN':
                if Patch.plateLocation == 'E':
                    _poloidal_f = self.yaml['target_plates']['plate_E1']['poloidal_f']
                    valid_function=self.CheckFunction(_poloidal_f,Enforce)
                elif Patch.plateLocation == 'W':
                    _poloidal_f = self.yaml['target_plates']['plate_W1']['poloidal_f']
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
            np_sol = self.yaml['grid_params']['grid_generation']['np_sol']
        except:
            np_sol = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_core = self.yaml['grid_params']['grid_generation']['np_core']
        except:
            np_core = self.yaml['grid_params']['grid_generation']['np_global']
        try:
            np_pf = self.yaml['grid_params']['grid_generation']['np_pf']
        except:
            np_pf = self.yaml['grid_params']['grid_generation']['np_global']

        if np_sol != np_core:
            if Enforce:
                raise ValueError('SOL and CORE must have equal POLOIDAL np values')
            else:
                print('WARNING: SOL and CORE must have equal POLOIDAL np values!\nSetting np values' \
                + ' to the minimum of np_sol={} and np_core={}.\n'.format(np_sol,np_core))
            np_sol = np_core = np.amin([np_sol, np_core])

        try:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_sol']
        except:
            nr_sol = self.yaml['grid_params']['grid_generation']['nr_global']

        try:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_core']
        except:
            nr_core = self.yaml['grid_params']['grid_generation']['nr_global']
        try:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_pf']
        except:
            nr_pf = self.yaml['grid_params']['grid_generation']['nr_global']

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
                    np_cells = self.yaml['target_plates']['plate_W1']['np_local']
                elif Patch.plateLocation == 'E':
                    np_cells = self.yaml['target_plates']['plate_E1']['np_local']
            if self.config == 'USN':
                if Patch.plateLocation == 'E':
                    np_cells = self.yaml['target_plates']['plate_W1']['np_local']
                elif Patch.plateLocation == 'W':
                    np_cells = self.yaml['target_plates']['plate_E1']['np_local']

        return (nr_cells,np_cells)
    def AdjustPatch(self,patch):
        primary_xpt = Point([self.yaml['grid_params']['rxpt'], self.yaml['grid_params']['zxpt']])
        # if self.config == 'USN':
        #     if patch.patchName == 'F2':
        #         patch.adjust_corner(primary_xpt, 'SE')
        #     elif patch.patchName == 'F1':
        #         patch.adjust_corner(primary_xpt, 'NE')
        #     elif patch.patchName == 'E2':
        #         patch.adjust_corner(primary_xpt, 'SW')
        #     elif patch.patchName == 'E1':
        #         patch.adjust_corner(primary_xpt, 'NW')
        #     elif patch.patchName == 'B1':
        #         patch.adjust_corner(primary_xpt, 'NE')
        #     elif patch.patchName == 'B2':
        #         patch.adjust_corner(primary_xpt, 'SE')
        #     elif patch.patchName == 'A1':
        #         patch.adjust_corner(primary_xpt, 'NW')
        #     elif patch.patchName == 'A2':
        #         patch.adjust_corner(primary_xpt, 'SW')
        # if self.config == 'LSN':

        tag = patch.name2tag()
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
                
                
                
    def CheckPatches(self,verbose=False):
        for name, patch in self.patches.items():
            if patch.platePatch:
                print(' # Checking patch: ', name)
                patch.CheckPatch(self)

    def SavePatches(self,FileName):
        with open(FileName,'w') as File:
            _yaml_.dump(self.patches,File)

            # self.platePatch = platePatch
            # self.plateLocation = plateLocation
            # self.patchName = patchName
            # self.lines
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


        # if self.config == 'USN':
        #     self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
        #                 [[None], p['F2'], p['E2'], p['D2'], p['C2'], p['B2'], p['A2'], [None]], \
        #                 [[None], p['F1'], p['E1'], p['D1'], p['C1'], p['B1'], p['A1'], [None]], \
        #                 [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
        #                 ]

        # if self.config == 'LSN':
        self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
                    [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
                    [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
                    [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
                    ]
        self.categorize_patches()

    def construct_grid(self, np_cells = 1, nr_cells = 1,Verbose=False,ShowVertices=False,RestartScratch=False,OptionTrace='theta',ExtraSettings={},ListPatches='all'):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.
        if Verbose: print('Construct Grid')
        try:
            visual = self.yaml['DEBUG']['visual']['subgrid']
        except:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['grid_generation']
        except:
            verbose = False

        verbose=Verbose or verbose
        
            
        print('>>> Patches:', [k for k in self.patches.keys()])
        if RestartScratch:
            self.CurrentListPatch={}
        self.GetConnexionMap() 
    
        for name, patch in self.patches.items():
            
            if self.CorrectDistortion.get(name) is not None:
               patch.CorrectDistortion=self.CorrectDistortion.get(name)
            elif self.CorrectDistortion.get('all') is not None:
                patch.CorrectDistortion=self.CorrectDistortion.get('all')
            else:
                patch.CorrectDistortion={'Active':False}
            if (ListPatches=='all' and patch not in self.CurrentListPatch) or (ListPatches!='all' and name in ListPatches):
                self.SetPatchBoundaryPoints(patch)
                (nr_cells,np_cells)=self.GetNpoints(patch)
                (_radial_f,_poloidal_f)=self.GetFunctions(patch,ExtraSettings=ExtraSettings)
                print('>>> Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(name, np_cells, nr_cells, inspect.getsource(_poloidal_f), inspect.getsource(_radial_f)))
                patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f,verbose = verbose, visual = visual,ShowVertices=ShowVertices,OptionTrace=OptionTrace)
                self.AdjustPatch(patch)
                patch.plot_subgrid()
                plt.pause(1)
                plt.show()
                plt.ion()
                self.CurrentListPatch[name] = patch



        if all(['cell_grid' in patch.__dict__ for patch in self.patches.values()]):
            self.concat_grid()
            self.set_gridue()
        
    def SetPatchBoundaryPoints(self,Patch):
            if self.ConnexionMap.get(Patch.patchName) is not None:
                if self.Verbose: print('Find connexion map for patch {}'.format(Patch.patchName))
                for Boundary,AdjacentPatch in self.ConnexionMap.get(Patch.patchName).items():
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
            
    def GetConnexionMap(self):
        # Connexion MAP for LSN
        self.ConnexionMap={}
        for name, patch in self.patches.items():
            if name[1]=='C':
                self.ConnexionMap[name]={'N':(name[0]+'S'+name[2],'S')}
        self.ConnexionMap['IPF']={'N':('IDL','S')}
        self.ConnexionMap['OPF']={'N':('ODL','S')}
        if self.Verbose: print('ConnexionMap:',self.ConnexionMap)

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
            visual = self.yaml['DEBUG']['visual']['patch_map']
        except KeyError:
            visual = False
        try:
            verbose = self.yaml['DEBUG']['verbose']['patch_generation']
        except KeyError:
            verbose = False
        try:
            inner_tilt = self.yaml['grid_params']['patch_generation']['inner_tilt']
        except KeyError:
            inner_tilt = 0.0
        try:
            outer_tilt = self.yaml['grid_params']['patch_generation']['outer_tilt']
        except KeyError:
            outer_tilt = 0.0


        WestPlate = Line([Point(i) for i in self.plate_data['plate_W1']['coordinates']])
        EastPlate = Line([Point(i) for i in self.plate_data['plate_E1']['coordinates']])

        xpt = self.eq.NSEW_lookup['xpt1']['coor']
        magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
            self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

        psi_max = self.yaml['grid_params']['psi_max']
        psi_min_core = self.yaml['grid_params']['psi_min_core']
        psi_min_pf = self.yaml['grid_params']['psi_min_pf']

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
        if self.yaml['grid_params']['patch_generation']['use_NW']:
            tilt = self.eq.NSEW_lookup['xpt1']['theta']['NW'] + self.yaml['grid_params']['patch_generation']['NW_adjust']
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
        else:
            xptW_psiMax = self.eq.draw_line(xpt['W'], {'psi' : psi_max}, option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

        if self.yaml['grid_params']['patch_generation']['use_NE']:
            tilt = self.eq.NSEW_lookup['xpt1']['theta']['NE'] + self.yaml['grid_params']['patch_generation']['NE_adjust']
            xptE_psiMax = self.eq.draw_line(xpt['E'], {'psi_horizontal' : (psi_max, tilt)}, option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
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
        A2 = SNL_Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'IDL', platePatch = True, plateLocation = location)

        # A1 Patch
        location = 'W'
        A1_N = A2_S.reverse_copy()
        
        A1_S = psiMinPF_WestPlate
        A1_E = xptS_psiMinPF
        A1_W = (WestPlate.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
        A1 = SNL_Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'IPF', platePatch = True, plateLocation = location)

        # B2 Patch
        
        B2_N = self.eq.draw_line(A2_N.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B2_S = xptNW_midLine.reverse_copy()
        B2_E = Line([B2_N.p[-1], B2_S.p[0]])
        B2_W = xptW_psiMax
        B2 = SNL_Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'ISB')

        # B1 Patch
        B1_N = B2_S.reverse_copy()
        B1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : inner_midLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        B1_E = Line([B1_N.p[-1], B1_S.p[0]])
        B1_W = xptN_psiMinCore.reverse_copy()
        B1 = SNL_Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'ICB')

        # C2 Patch
        C2_N = self.eq.draw_line(B2_N.p[-1], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C2_S = imidLine_topLine.reverse_copy()
        C2_E = Line([C2_N.p[-1], C2_S.p[0]])
        C2_W = Line([C2_S.p[-1], C2_N.p[0]])
        C2 = SNL_Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'IST')

        # C1 Patch
        C1_N = C2_S.reverse_copy()
        C1_S = self.eq.draw_line(B1_S.p[0], {'line' : topLine}, option = 'theta', direction = 'cw', show_plot = visual, text = verbose).reverse_copy()
        C1_E = Line([C1_N.p[-1], C1_S.p[0]])
        C1_W = Line([C1_S.p[-1], C1_N.p[0]])
        C1 = SNL_Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'ICT')

        # F2 Patch
        location = 'E'
        F2_N = oPsiMax_TP
        F2_S = xpt_EastPlate.reverse_copy()
        F2_E = (EastPlate.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point = True)[0]
        F2_W = xptE_psiMax
        F2 = SNL_Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'ODL', platePatch = True, plateLocation = location)

        # F1 Patch
        location = 'E'
        F1_N = F2_S.reverse_copy()
        F1_S = psiMinPF_EastPlate.reverse_copy()
        F1_E = (EastPlate.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point = True)[0]
        F1_W = xptS_psiMinPF.reverse_copy()
        F1 = SNL_Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'OPF', platePatch = True, plateLocation = location)

        # E2 Patch
        E2_N = self.eq.draw_line(F2_N.p[0], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        E2_S = xptNE_midLine
        E2_E = xptE_psiMax.reverse_copy()
        E2_W = Line([E2_S.p[-1], E2_N.p[0]])
        E2 = SNL_Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'OSB')

        # E1 Patch
        E1_N = E2_S.reverse_copy()
        E1_S = self.eq.draw_line(xptN_psiMinCore.p[-1], {'line' : outer_midLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E1_E = xptN_psiMinCore
        E1_W = Line([E1_S.p[-1], E1_N.p[0]])
        E1 = SNL_Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'OCB')

        # D2 Patch
        D2_N = self.eq.draw_line(E2_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D2_S = omidLine_topLine
        D2_E = Line([D2_N.p[-1], D2_S.p[0]])
        D2_W = Line([D2_S.p[-1], D2_N.p[0]])
        D2 = SNL_Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'OST')

        # D1 Patch
        D1_N = D2_S.reverse_copy()
        D1_S = self.eq.draw_line(E1_S.p[-1], {'line' : topLine}, option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_E = Line([D1_N.p[-1], D1_S.p[0]])
        D1_W = Line([D1_S.p[-1], D1_N.p[0]])
        D1 = SNL_Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'OCT')

        patches = [A2, A1, B2, B1, C2, C1, D2, D1, E2, E1, F2, F1]

        self.patches = {}
        for patch in patches:
            self.patches[patch.patchName] = patch
            patch.plot_border()
            patch.fill()

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


# class DNL(Ingrid.Ingrid):
#     def __init__(self, INGRID_object):
#         super().__init__(params = INGRID_object.yaml)
#         self.efit_psi = INGRID_object.efit_psi
#         self.psi_norm = INGRID_object.psi_norm
#         self.eq = INGRID_object.eq
#         self.plate_data = INGRID_object.plate_data
#         self.config = 'DNL'
#         print('=' * 80)
#         print('DNL Object!')
#         print('=' * 80 + '\n')

#     def patch_diagram(self):
#         """ Generates the patch diagram for a given configuration. """

#         colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
#                   'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
#                   'seagreen', 'firebrick', 'saddlebrown', 'c',
#                   'm', 'dodgerblue', 'darkorchid', 'crimson',
#                   'darkorange', 'lightgreen', 'lightseagreen', 'indigo',
#                   'mediumvioletred', 'mistyrose', 'darkolivegreen', 'rebeccapurple']

#         try:
#             plt.close('INGRID: Patch Map')
#         except:
#             pass

#         plt.figure('INGRID: Patch Map', figsize=(6, 10))
#         plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
#         plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
#         plt.gca().set_aspect('equal', adjustable='box')
#         plt.xlabel('R')
#         plt.ylabel('Z')
#         plt.title('DNL Patch Diagram')

#         for i in range(len(self.patches)):
#             self.patches[i].plot_border('green')
#             self.patches[i].fill(colors[i])

#         plt.show()

#     def grid_diagram(self):

#         try:
#             plt.close('INGRID: Grid')
#         except:
#             pass

#         plt.figure('INGRID: Grid', figsize=(6,10))
#         for patch in self.patches:
#             patch.plot_subgrid()

#         plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
#         plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
#         plt.gca().set_aspect('equal', adjustable='box')
#         plt.xlabel('R')
#         plt.ylabel('Z')
#         plt.title('INGRID DNL Subgrid')
#         plt.show()

#     def animate_grid(self):

#         try:
#             plt.close('INGRID: Debug')
#         except:
#             pass
#         plt.figure('INGRID: Debug', figsize=(6, 10))
#         plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
#         plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
#         plt.gca().set_aspect('equal', adjustable='box')
#         plt.xlabel('R')
#         plt.ylabel('Z')
#         plt.title('visualize gridue')

#         k = [1,2,4,3,1]

#         for i in range(len(self.rm)):
#             for j in range(len(self.rm[0])):
#                 plt.plot(self.rm[i][j][k], self.zm[i][j][k])
#                 plt.pause(0.01)

#     def concat_grid(self):
#         """
#         Concatenate all local grids on individual patches into a single
#         array with branch cuts
#         Parameters:
#         ----------
#             config : str
#                 Type of SNL grid to concat.
#         """
#         # Patch Matrix corresponds to the SNL Patch Map (see GINGRED paper).
#         patch_matrix = self.patch_matrix

#         # Get some poloidal and radial information from each patch to attribute to the
#         # local subgrid.
#         # NOTE: npol and nrad refer to the actual lines in the subgrid. Because of this, we must add
#         #       the value of 1 to the cell number to get the accurate number of lines.


#         for patch in self.patches:
#             patch.npol = len(patch.cell_grid[0]) + 1
#             patch.nrad = len(patch.cell_grid) + 1

#             print('"{}" has npol = {} and nrad = {}'.format(patch.patchName, patch.npol, patch.nrad))

#         # Total number of poloidal indices in all subgrids.
#         np_total1 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:5]])) + 2

#         # Total number of radial indices in all subgrids.
#         nr_total1 = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:4]])) + 2

#         # Total number of poloidal indices in all subgrids.
#         np_total2 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:11]])) + 2

#         # Total number of radial indices in all subgrids.
#         nr_total2 = int(np.sum([patch[7].nrad - 1 for patch in patch_matrix[1:4]])) + 2

#         rm1 = np.zeros((np_total1, nr_total1, 5), order = 'F')
#         zm1  = np.zeros((np_total1, nr_total1, 5), order = 'F')
#         rm2 = np.zeros((np_total2, nr_total2, 5), order = 'F')
#         zm2  = np.zeros((np_total2, nr_total2, 5), order = 'F')

#         ixcell = 0
#         jycell = 0

#         # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')
#         for ixp in range(1, 5):

#             nr_sum = 0
#             for jyp in range(1, 4):
#                 # Point to the current patch we are operating on.
#                 local_patch = patch_matrix[jyp][ixp]

#                 if local_patch == [None]:
#                     continue

#                 nr_sum += local_patch.nrad - 1

#                 # Access the grid that is contained within this local_patch.
#                 # ixl - number of poloidal cells in the patch.
#                 for ixl in range(len(local_patch.cell_grid[0])):
#                     # jyl - number of radial cells in the patch
#                     for jyl in range(len(local_patch.cell_grid)):

#                         ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:ixp+1]])) \
#                                 - len(local_patch.cell_grid[0]) + ixl + 1

#                         jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

#                         ind = 0
#                         for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
#                             rm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
#                             zm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
#                             ind += 1
#                             if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

#         # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')

#         ixcell = 0
#         jycell = 0

#         for ixp in range(7, 11):

#             nr_sum = 0
#             for jyp in range(1, 4):
#                 # Point to the current patch we are operating on.
#                 local_patch = patch_matrix[jyp][ixp]

#                 if local_patch == [None]:
#                     continue

#                 nr_sum += local_patch.nrad - 1

#                 # Access the grid that is contained within this local_patch.
#                 # ixl - number of poloidal cells in the patch.
#                 for ixl in range(len(local_patch.cell_grid[0])):
#                     # jyl - number of radial cells in the patch
#                     for jyl in range(len(local_patch.cell_grid)):

#                         ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:ixp+1]])) \
#                                 - len(local_patch.cell_grid[0]) + ixl + 1

#                         jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

#                         ind = 0
#                         for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
#                             rm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
#                             zm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
#                             ind += 1
#                             if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

#         # Flip indices into gridue format.
#         for i in range(len(rm1)):
#             rm1[i] = rm1[i][::-1]
#         for i in range(len(zm1)):
#             zm1[i] = zm1[i][::-1]
#         for i in range(len(rm2)):
#             rm2[i] = rm2[i][::-1]
#         for i in range(len(zm2)):
#             zm2[i] = zm2[i][::-1]

#         # Add guard cells to the concatenated grid.
#         ixrb1 = len(rm1) - 2
#         ixlb1 = 0
#         ixrb2 = len(rm2) - 2
#         ixlb2 = 0

#         rm1 = self.add_guardc(rm1, ixlb1, ixrb1)
#         zm1 = self.add_guardc(zm1, ixlb1, ixrb1)
#         rm2 = self.add_guardc(rm2, ixlb2, ixrb2)
#         zm2 = self.add_guardc(zm2, ixlb2, ixrb2)

#         self.rm = np.concatenate((rm1, rm2))
#         self.zm = np.concatenate((zm1, zm2))

#         try:
#             debug = self.yaml['DEBUG']['visual']['gridue']
#         except:
#             debug = False

#         if debug:
#             self.animate_grid()

#     def add_guardc(self, cell_map, ixlb, ixrb, eps = 1e-3):

#         def set_guard(cell_map, ix, iy, eps, boundary):
#             # Note: 'USN' and 'right' is really just 'LSN' and 'left' settings.
#             # TODO: Edit the code to reflect this at some point so the next reader is not overwhelmed.
#             if boundary == 'left':
#                 ixn = ix + 1
#                 iyn = iy
#                 cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][2])
#                 cell_map[ix][iy][2] = cell_map[ixn][iyn][1]
#                 cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][4])
#                 cell_map[ix][iy][4] = cell_map[ixn][iyn][3]
#                 cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

#             elif boundary == 'right':
#                 ixn = ix - 1
#                 iyn = iy
#                 cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][1])
#                 cell_map[ix][iy][1] = cell_map[ixn][iyn][2]
#                 cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][3])
#                 cell_map[ix][iy][3] = cell_map[ixn][iyn][4]
#                 cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

#             elif boundary == 'bottom':
#                 ixn = ix
#                 iyn = iy + 1
#                 cell_map[ix][iy][1] = cell_map[ixn][iyn][1] + eps * (cell_map[ixn][iyn][1] - cell_map[ixn][iyn][3])
#                 cell_map[ix][iy][3] = cell_map[ixn][iyn][1]
#                 cell_map[ix][iy][2] = cell_map[ixn][iyn][2] + eps * (cell_map[ixn][iyn][2] - cell_map[ixn][iyn][4])
#                 cell_map[ix][iy][4] = cell_map[ixn][iyn][2]
#                 cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])
#             elif boundary == 'top':
#                 ixn = ix
#                 iyn = iy - 1
#                 cell_map[ix][iy][3] = cell_map[ixn][iyn][3] + eps * (cell_map[ixn][iyn][3] - cell_map[ixn][iyn][1])
#                 cell_map[ix][iy][1] = cell_map[ixn][iyn][3]
#                 cell_map[ix][iy][4] = cell_map[ixn][iyn][4] + eps * (cell_map[ixn][iyn][4] - cell_map[ixn][iyn][2])
#                 cell_map[ix][iy][2] = cell_map[ixn][iyn][4]
#                 cell_map[ix][iy][0] = 0.25 * (cell_map[ix][iy][1] + cell_map[ix][iy][2] + cell_map[ix][iy][3] + cell_map[ix][iy][4])

#             return cell_map

#         np = len(cell_map) - 2
#         nr = len(cell_map[0]) - 2

#         for iy in range(1, nr + 1):
#             ix = ixlb
#             cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'left')
#             ix = ixrb + 1
#             cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'right')

#         for ix in range(np + 2):
#             iy = 0
#             cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'bottom')
#             iy = nr + 1
#             cell_map = set_guard(cell_map, ix, iy, eps, boundary = 'top')

#         return cell_map

#     def set_gridue(self):
#         """
#         Prepare the relevant arrays for writing to GRIDUE.
#         """

#         ixlb = 0
#         ixrb = len(self.rm) - 2

#         nxm = len(self.rm) - 2
#         nym = len(self.rm[0]) - 2
#         iyseparatrix1 = self.patch_lookup['C2'].nrad - 1
#         iyseparatrix2 = self.patch_lookup['B5'].nrad + self.patch_lookup['C5'].nrad - 2
#         ix_plate1 = 0
#         ix_cut1 = self.patch_lookup['A1'].npol - 1
#         ix_cut2 = self.patch_lookup['A1'].npol + self.patch_lookup['A2'].npol + self.patch_lookup['A3'].npol - 3
#         ix_plate2 = ix_cut2 + self.patch_lookup['A4'].npol - 1
#         iyseparatrix3 = iyseparatrix2
#         iyseparatrix4 = iyseparatrix1
#         ix_plate3 = ix_plate2 + 2
#         ix_cut3 = ix_plate3 + self.patch_lookup['A5'].npol - 1
#         ix_cut4 = ix_cut3 + self.patch_lookup['A6'].npol + self.patch_lookup['A7'].npol - 2
#         ix_plate4 = ix_cut4 + self.patch_lookup['A8'].npol - 1

#         psi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
#         br = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
#         bz = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
#         bpol = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
#         bphi = np.zeros((nxm + 2, nym + 2, 5), order = 'F')
#         b = np.zeros((nxm + 2, nym + 2, 5), order = 'F')

#         rm = self.rm
#         zm = self.zm
#         rb_prod = self.efit_psi.rcenter * self.efit_psi.bcenter

#         for i in range(len(b)):
#             for j in range(len(b[0])):
#                 for k in range(5):
#                     _r = rm[i][j][k]
#                     _z = zm[i][j][k]

#                     _psi = self.efit_psi.get_psi(_r, _z)
#                     _br = self.efit_psi.get_psi(_r, _z, tag = 'vz') / _r
#                     _bz = -self.efit_psi.get_psi(_r, _z, tag = 'vr') / _r
#                     _bpol = np.sqrt(_br ** 2 + _bz ** 2)
#                     _bphi = rb_prod / _r
#                     _b = np.sqrt(_bpol ** 2 + _bphi ** 2)

#                     psi[i][j][k] = _psi
#                     br[i][j][k] = _br
#                     bz[i][j][k] = _bz
#                     bpol[i][j][k] = _bpol
#                     bphi[i][j][k] = _bphi
#                     b[i][j][k] = _b

#         self.gridue_params = {'nxm' : nxm, 'nym' : nym, 'iyseparatrix1' : iyseparatrix1, 'iyseparatrix2' : iyseparatrix2, \
#                 'ix_plate1' : ix_plate1, 'ix_cut1' : ix_cut1, 'ix_cut2' : ix_cut2, 'ix_plate2' : ix_plate2, 'iyseparatrix3' : iyseparatrix3, \
#                 'iyseparatrix4' : iyseparatrix4, 'ix_plate3' : ix_plate3, 'ix_cut3' : ix_cut3, 'ix_cut4' : ix_cut4, 'ix_plate4' : ix_plate4, \
#                 'rm' : self.rm, 'zm' : self.zm, 'psi' : psi, 'br' : br, 'bz' : bz, 'bpol' : bpol, 'bphi' : bphi, 'b' : b, '_FILLER_' : -1}


#     def construct_grid(self, np_cells = 3, nr_cells = 3,Verbose=False):

#         # Straighten up East and West segments of our patches,
#         # Plot borders and fill patches.

#         try:
#             visual = self.yaml['DEBUG']['visual']['subgrid']
#         except:
#             visual = False
#         try:
#             verbose = self.yaml['DEBUG']['verbose']['grid_generation']
#         except:
#             verbose = False

#         verbose=verbose or Verbose

#         try:
#             np_primary_sol = self.yaml['grid_params']['grid_generation']['np_primary_sol']
#         except:
#             np_primary_sol = self.yaml['grid_params']['grid_generation']['np_global']
#         try:
#             np_secondary_sol = self.yaml['grid_params']['grid_generation']['np_secondary_sol']
#         except:
#             np_secondary_sol = self.yaml['grid_params']['grid_generation']['np_global']
#         try:
#             np_core = self.yaml['grid_params']['grid_generation']['np_core']
#         except:
#             np_core = self.yaml['grid_params']['grid_generation']['np_global']
#         try:
#             np_primary_pf = self.yaml['grid_params']['grid_generation']['np_primary_pf']
#         except:
#             np_primary_pf = self.yaml['grid_params']['grid_generation']['np_global']
#         try:
#             np_secondary_pf = self.yaml['grid_params']['grid_generation']['np_secondary_pf']
#         except:
#             np_secondary_pf = self.yaml['grid_params']['grid_generation']['np_global']
#         """
#         if np_sol != np_core:
#             print('WARNING: SOL and CORE must have equal POLOIDAL np values!\nSetting np values' \
#                 + ' to the minimum of np_sol and np_core.\n')
#             np_sol = np_core = np.amin([np_sol, np_core])
#         """
#         try:
#             nr_primary_sol = self.yaml['grid_params']['grid_generation']['nr_primary_sol']
#         except:
#             nr_primary_sol = self.yaml['grid_params']['grid_generation']['nr_global']

#         try:
#             nr_secondary_sol = self.yaml['grid_params']['grid_generation']['nr_secondary_sol']
#         except:
#             nr_secondary_sol = self.yaml['grid_params']['grid_generation']['nr_global']

#         try:
#             nr_core = self.yaml['grid_params']['grid_generation']['nr_core']
#         except:
#             nr_core = self.yaml['grid_params']['grid_generation']['nr_global']
#         try:
#             nr_primary_pf = self.yaml['grid_params']['grid_generation']['nr_primary_pf']
#         except:
#             nr_primary_pf = self.yaml['grid_params']['grid_generation']['nr_global']
#         try:
#             nr_secondary_pf = self.yaml['grid_params']['grid_generation']['nr_secondary_pf']
#         except:
#             nr_secondary_pf = self.yaml['grid_params']['grid_generation']['nr_global']
#         """
#         if nr_pf != nr_core:
#             print('WARNING: PF and CORE must have equal RADIAL nr values!\nSetting nr values' \
#                 + ' to the minimum of nr_pf and nr_core.\n')
#             nr_pf = nr_core = np.amin([nr_pf, nr_core])
#         """
#         for patch in self.patches:
#             if patch in self.PRIMARY_SOL:
#                 nr_cells = nr_primary_sol
#                 np_cells = np_primary_sol
#                 print('Patch "{}" is in PRIMARY_SOL'.format(patch.patchName))
#             elif patch in self.SECONDARY_SOL:
#                 nr_cells = nr_secondary_sol
#                 np_cells = np_secondary_sol
#                 print('Patch "{}" is in SECONDARY_SOL'.format(patch.patchName))
#             elif patch in self.CORE:
#                 nr_cells = nr_core
#                 np_cells = np_core
#                 print('Patch "{}" is in CORE'.format(patch.patchName))
#             elif patch in self.PRIMARY_PF:
#                 nr_cells = nr_primary_pf
#                 np_cells = np_primary_pf
#                 print('Patch "{}" is in PRIMARY_PF'.format(patch.patchName))
#             elif patch in self.SECONDARY_PF:
#                 nr_cells = nr_secondary_pf
#                 np_cells = np_secondary_pf
#                 print('Patch "{}" is in SECONDARY_PF'.format(patch.patchName))

#             print('CONSTRUCTING GRID FOR PATCH: {}'.format(patch.patchName))
#             patch.make_subgrid(self, np_cells, nr_cells, verbose = verbose, visual = visual)

#             # Tidy up primary x-point
#             primary_xpt = Point(self.xpt1)
#             secondary_xpt = Point(self.xpt2)

#             if patch.patchName == 'B1':
#                 patch.adjust_corner(primary_xpt, 'SE')
#             elif patch.patchName == 'C1':
#                 patch.adjust_corner(primary_xpt, 'NE')
#             elif patch.patchName == 'B2':
#                 patch.adjust_corner(primary_xpt, 'SW')
#             elif patch.patchName == 'C2':
#                 patch.adjust_corner(primary_xpt, 'NW')
#             elif patch.patchName == 'C7':
#                 patch.adjust_corner(primary_xpt, 'NE')
#             elif patch.patchName == 'B7':
#                 patch.adjust_corner(primary_xpt, 'SE')
#             elif patch.patchName == 'C8':
#                 patch.adjust_corner(primary_xpt, 'NW')
#             elif patch.patchName == 'B8':
#                 patch.adjust_corner(primary_xpt, 'SW')

#             # Tidy up secondary x-point
#             elif patch.patchName == 'A3':
#                 patch.adjust_corner(secondary_xpt, 'SE')
#             elif patch.patchName == 'B3':
#                 patch.adjust_corner(secondary_xpt, 'NE')
#             elif patch.patchName == 'A4':
#                 patch.adjust_corner(secondary_xpt, 'SW')
#             elif patch.patchName == 'B4':
#                 patch.adjust_corner(secondary_xpt, 'NW')
#             elif patch.patchName == 'B5':
#                 patch.adjust_corner(secondary_xpt, 'NE')
#             elif patch.patchName == 'A5':
#                 patch.adjust_corner(secondary_xpt, 'SE')
#             elif patch.patchName == 'B6':
#                 patch.adjust_corner(secondary_xpt, 'NW')
#             elif patch.patchName == 'A6':
#                 patch.adjust_corner(secondary_xpt, 'SW')

#         self.concat_grid()
#         self.set_gridue()

#     def construct_patches(self):
#         """
#         Draws lines and creates patches for both USN and LSN configurations.

#         Patch Labeling Key:
#             I: Inner,
#             O: Outer,
#             DL: Divertor Leg,
#             PF: Private Flux,
#             T: Top,
#             B: Bottom,
#             S: Scrape Off Layer,
#             C: Core.
#         """
#         # TODO: Create a 'lookup' procedure for determining line drawing
#         #       orientations and inner-outer locations.

#         def order_plate_points(plate, location = 'UITP', orientation = 'cw'):
#             """
#             Sets the points in the target plate to have an orientation
#             increasing in R and Z. Checks the endpoints to satisfy this criteria.
#             """

#             loc_sgn = 1 if location[0] == 'U' else -1

#             start = plate.p[0]
#             end = plate.p[-1]
#             # Endpoints on same vertical line.
#             if start.x == end.x:
#                 # If rhs endpoint above lhs endpoint.
#                 if end.y - start.y > 0:
#                     # Return target plate as is
#                     return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
#                 else:
#                     # Else flip plate orientation.
#                     return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p
#             # Endpoints on same horizontal line.
#             elif start.y == end.y:
#                 # If lhs endpoint to the left of rhs endpoint.
#                 if end.x - start.x > 0:
#                     # Return target plate as is
#                     return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
#                 else:
#                     # Else flip plate orientation
#                     return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p
#             # Endpoints are on sloped line.
#             # Check if lhs endpoint is on the left of rhs endpoint
#             elif loc_sgn * (end.x - start.x) > 0:
#                 return plate.copy().p if orientation == 'cw' else plate.copy().p[::-1]
#             else:
#                 return plate.copy().p[::-1] if orientation == 'cw' else plate.copy().p


#         try:
#             visual = self.yaml['DEBUG']['visual']['patch_map']
#         except KeyError:
#             visual = False
#         try:
#             verbose = self.yaml['DEBUG']['verbose']['patch_generation']
#         except KeyError:
#             verbose = False
#         try:
#             inner_tilt = self.yaml['grid_params']['patch_generation']['inner_tilt']
#         except KeyError:
#             inner_tilt = 0.0
#         try:
#             outer_tilt = self.yaml['grid_params']['patch_generation']['outer_tilt']
#         except KeyError:
#             outer_tilt = 0.0

#         self.plate_W1 = self.plate_data['plate_W1']['coordinates']
#         self.plate_E1 = self.plate_data['plate_E1']['coordinates']
#         self.plate_E2 = self.plate_data['plate_E2']['coordinates']
#         self.plate_W2 = self.plate_data['plate_W2']['coordinates']

#         xpt1_dict = self.eq.NSEW_lookup['xpt1']['coor']
#         xpt1_theta = self.eq.NSEW_lookup['xpt1']['theta']
#         xpt2_dict = self.eq.NSEW_lookup['xpt2']['coor']
#         xpt1_theta = self.eq.NSEW_lookup['xpt2']['theta']

#         sptrx1_v = self.psi_norm.get_psi(self.xpt1[0], self.xpt1[1])
#         sptrx2_v = self.psi_norm.get_psi(self.xpt2[0], self.xpt2[1])

#         magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
#             self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

#         psi_max = self.yaml['grid_params']['psi_max']
#         psi_min_core = self.yaml['grid_params']['psi_min_core']
#         psi_min_pf = self.yaml['grid_params']['psi_min_pf']

#         psi_max_outer = self.yaml['grid_params']['psi_max_outer']
#         psi_max_inner = self.yaml['grid_params']['psi_max_inner']
#         psi_min_pf_2 = self.yaml['grid_params']['psi_pf2']

#         LITP = Line(order_plate_points(Line([Point(i) for i in self.plate_W1]), location = 'LITP'))
#         LOTP = Line(order_plate_points(Line([Point(i) for i in self.plate_E1]), location = 'LOTP'))

#         UITP = Line(order_plate_points(Line([Point(i) for i in self.plate_E2]), location = 'UITP'))
#         UOTP = Line(order_plate_points(Line([Point(i) for i in self.plate_W2]), location = 'UOTP'))

#         # Generate Horizontal Mid-Plane lines
#         LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
#         RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
#         inner_midLine = Line([LHS_Point, RHS_Point])
#         # inner_midLine.plot()

#         LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
#         RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
#         outer_midLine = Line([LHS_Point, RHS_Point])
#         # outer_midLine.plot()

#         # Generate Vertical Mid-Plane line that intersects the secondary x-pt and the magnetic axis
#         slp = [self.xpt2[0] - magx[0], self.xpt2[1] - magx[1]]
#         slp = slp[1] / slp[0]
#         topLine_tilt = np.arctan(slp)

#         Upper_Point = Point(-1e6, slp * (-1e6 - self.xpt2[0]) + self.xpt2[1])
#         Lower_Point = Point(1e6, slp * (1e6 - self.xpt2[0]) + self.xpt2[1])
#         topLine = Line([Lower_Point, Upper_Point])
#         # topLine.plot()

#         # Drawing the portion of separatrix similar to single-null configuration.
#         xpt1N__psiMinCore = self.eq.draw_line(xpt1_dict['N'], {'psi' : psi_min_core}, \
#             option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
#         xpt1NW__sptrx1imidLine = self.eq.draw_line(xpt1_dict['NW'], {'line' : inner_midLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         xpt1NE__sptrx1omidLine = self.eq.draw_line(xpt1_dict['NE'], {'line' : outer_midLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

#         sptrx1imidLine__topLine = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : topLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx1omidLine__topLine = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : topLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         psiMinCore__imidLineCore = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : inner_midLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         psiMinCore__omidLineCore = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : outer_midLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         imidLineCore__topLine = self.eq.draw_line(psiMinCore__imidLineCore.p[-1], {'line' : topLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         omidLineCore__topLine = self.eq.draw_line(psiMinCore__omidLineCore.p[-1], {'line' : topLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         sptrx1 = Line([xpt1NW__sptrx1imidLine.p + sptrx1imidLine__topLine.p + sptrx1omidLine__topLine.p[-2::-1] + xpt1NE__sptrx1omidLine.p[-2::-1]][0])
#         core = Line([psiMinCore__imidLineCore.p + imidLineCore__topLine.p + omidLineCore__topLine.p[-2::-1] + psiMinCore__omidLineCore.p[-2::-1]][0])

#         xpt1__psiMinPF1 = self.eq.draw_line(xpt1_dict['S'], {'psi' : psi_min_pf}, \
#             option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
#         psiMinPF1__LITP = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : LITP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         psiMinPF1__LOTP = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : LOTP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         xpt1__LITP = self.eq.draw_line(xpt1_dict['SW'], {'line' : LITP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         xpt1__LOTP = self.eq.draw_line(xpt1_dict['SE'], {'line' : LOTP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)

#         # Drawing the portion of the separatrix found in the double-null configuration
#         xpt2__sptrx1 = self.eq.draw_line(xpt2_dict['N'], {'line' : (sptrx1, topLine_tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         xpt2NE__sptrx2imidLine = self.eq.draw_line(xpt2_dict['NE'], {'line' : inner_midLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         xpt2NW__sptrx2omidLine = self.eq.draw_line(xpt2_dict['NW'], {'line' : outer_midLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx2imidLine__LITP = self.eq.draw_line(xpt2NE__sptrx2imidLine.p[-1], {'line' : LITP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         sptrx2omidLine__LOTP = self.eq.draw_line(xpt2NW__sptrx2omidLine.p[-1], {'line' : LOTP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx2_inner = Line([xpt2NE__sptrx2imidLine.p + sptrx2imidLine__LITP.p][0])
#         sptrx2_outer = Line([xpt2NW__sptrx2omidLine.p + sptrx2omidLine__LOTP.p][0])

#         xpt2__psiMinPF2 = self.eq.draw_line(xpt2_dict['S'], {'psi' : psi_min_pf_2}, \
#             option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
#         xpt2__psiMinPF2_A, xpt2__psiMinPF2_B = xpt2__psiMinPF2.split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p)//2], add_split_point = True)

#         psiMinPF2_A__UITP = self.eq.draw_line(xpt2__psiMinPF2_B.p[0], {'line' : UITP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         psiMinPF2_A__UOTP = self.eq.draw_line(xpt2__psiMinPF2_B.p[0], {'line' : UOTP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         psiMinPF2_B__UITP = self.eq.draw_line(xpt2__psiMinPF2_B.p[-1], {'line' : UITP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose, debug = True)
#         psiMinPF2_B__UOTP = self.eq.draw_line(xpt2__psiMinPF2_B.p[-1], {'line' : UOTP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         xpt2__UITP = self.eq.draw_line(xpt2_dict['SE'], {'line' : UITP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         xpt2__UOTP = self.eq.draw_line(xpt2_dict['SW'], {'line' : UOTP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

#         if self.yaml['grid_params']['patch_generation']['use_secondary_NE']:
#             v = np.array([xpt2NE__sptrx2imidLine.p[1].x - self.xpt2[0], xpt2NE__sptrx2imidLine.p[1].y - self.xpt2[1]])
#             tilt = np.arccos(np.dot(v, np.array([-1, 0])) / np.linalg.norm(v)) + self.yaml['grid_params']['patch_generation']['secondary_NE_adjust']
#             xpt2__psiMaxInner = self.eq.draw_line(xpt2_dict['E'], {'psi_horizontal' : (psi_max_inner, tilt)}, \
#                 option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
#         else:
#             xpt2__psiMaxInner = self.eq.draw_line(xpt2_dict['E'], {'psi' : psi_max_inner}, \
#                 option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
#         if self.yaml['grid_params']['patch_generation']['use_secondary_NW']:
#             v = np.array([xpt2NW__sptrx2omidLine.p[1].x - self.xpt2[0], xpt2NW__sptrx2omidLine.p[1].y - self.xpt2[1]])
#             tilt = -np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) + self.yaml['grid_params']['patch_generation']['secondary_NW_adjust']
#             xpt2__psiMaxOuter = self.eq.draw_line(xpt2_dict['W'], {'psi_horizontal' : (psi_max_outer, tilt)}, \
#                 option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         else:
#             xpt2__psiMaxOuter = self.eq.draw_line(xpt2_dict['W'], {'psi' : psi_max_outer}, \
#                 option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)

#         psiMaxInner__UITP = self.eq.draw_line(xpt2__psiMaxInner.p[-1], {'line' : UITP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         psiMaxOuter__UOTP = self.eq.draw_line(xpt2__psiMaxOuter.p[-1], {'line' : UOTP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         psiMaxInner__imidLine = self.eq.draw_line(xpt2__psiMaxInner.p[-1], {'line' : inner_midLine}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         psiMaxOuter__omidLine = self.eq.draw_line(xpt2__psiMaxOuter.p[-1], {'line' : outer_midLine}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         imidLine__LITP = self.eq.draw_line(psiMaxInner__imidLine.p[-1], {'line' : LITP}, \
#             option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
#         omidLine__LOTP = self.eq.draw_line(psiMaxOuter__omidLine.p[-1], {'line' : LOTP}, \
#             option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx3_inner = Line([psiMaxInner__imidLine.p + imidLine__LITP.p][0])
#         sptrx3_outer = Line([psiMaxOuter__omidLine.p + omidLine__LOTP.p][0])

#         # Computing angle between line tangent to sptrx1 in NW direction and unit vector < 1, 0 >
#         v = np.array([xpt1NW__sptrx1imidLine.p[1].x - self.xpt1[0], xpt1NW__sptrx1imidLine.p[1].y - self.xpt1[1]])
#         tilt = np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) + np.pi/6
#         xpt1W__sptrx2Lower = self.eq.draw_line(xpt1_dict['W'], {'line' : (sptrx2imidLine__LITP, tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx2__LowerPsiMaxInner = self.eq.draw_line(xpt1W__sptrx2Lower.p[-1], {'line' : (imidLine__LITP, tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

#         # Computing angle between line tangent to sptrx1 in NE direction and unit vector < 1, 0 >
#         v = np.array([xpt1NE__sptrx1omidLine.p[1].x - self.xpt1[0], xpt1NE__sptrx1omidLine.p[1].y - self.xpt1[1]])
#         tilt = np.arccos(np.dot(v, np.array([1, 0])) / np.linalg.norm(v)) - np.pi/6
#         xpt1E__sptrx2Lower = self.eq.draw_line(xpt1_dict['E'], {'line': (sptrx2omidLine__LOTP, tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx2__LowerPsiMaxOuter = self.eq.draw_line(xpt1E__sptrx2Lower.p[-1], {'line' : (omidLine__LOTP, tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

#         sptrx1__core_top = self.eq.draw_line(xpt2__sptrx1.p[-1], {'line' : (core, topLine_tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx1__core_inner = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : (core, inner_tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx1__core_outer = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : (core, outer_tilt)}, \
#             option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
#         sptrx1__sptrx2_inner = self.eq.draw_line(xpt1NW__sptrx1imidLine.p[-1], {'line' : (sptrx2_inner, inner_tilt)}, \
#             option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
#         sptrx1__sptrx2_outer = self.eq.draw_line(xpt1NE__sptrx1omidLine.p[-1], {'line' : (sptrx2_outer, outer_tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)
#         sptrx2_inner__sptrx3_inner = self.eq.draw_line(xpt2NE__sptrx2imidLine.p[-1], {'line' : (sptrx3_inner, inner_tilt)}, \
#             option = 'z_const', direction = 'ccw', show_plot = visual, text = verbose)
#         sptrx2_outer__sptrx3_outer = self.eq.draw_line(xpt2NW__sptrx2omidLine.p[-1], {'line' : (sptrx3_outer, outer_tilt)}, \
#             option = 'z_const', direction = 'cw', show_plot = visual, text = verbose)

#         imidLine__LowerPsiMaxInner, LowerPsiMaxInner__LITP = imidLine__LITP.split(sptrx2__LowerPsiMaxInner.p[-1], add_split_point = True)
#         sptrx2imidLine__sptrx2Lower, sptrx2Lower__LITP = sptrx2imidLine__LITP.split(xpt1W__sptrx2Lower.p[-1], add_split_point = True)

#         # ============== Patch A1 ==============
#         location = 'W'
#         A1_N = LowerPsiMaxInner__LITP.reverse_copy()
#         A1_S = sptrx2Lower__LITP
#         A1_E = sptrx2__LowerPsiMaxInner.reverse_copy()
#         # =====================================================================================
#         # Trimming the target_plate to conform to the patch boundary.
#         # -------------------------------------------------------------------------------------
#         # Recall LITP has a clockwise orientation.
#         #
#         # The inner 'split' trims all Point objects BEFORE the point of intersection of LITP
#         # and A1_S. Call this new Line object Line_A.
#         #
#         # The outer 'split' trims all Point objects AFTER the point of intersection of Line_A
#         # and A1_N. This new Line object is the plate facing boundary of the Patch.
#         # =====================================================================================
#         A1_W = (LITP.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
#         A1 = DNL_Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'A1', platePatch = True, plateLocation = location)

#         # ============== Patch B1 ==============
#         location = 'W'
#         B1_N = A1_S.reverse_copy()
#         B1_E = xpt1W__sptrx2Lower.reverse_copy()
#         B1_S = xpt1__LITP
#         B1_W = (LITP.split(B1_S.p[-1])[1]).split(B1_N.p[0], add_split_point = True)[0]
#         B1 = DNL_Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'B1', platePatch = True, plateLocation = location)

#         # ============== Patch C1 ==============
#         location = 'W'
#         C1_N = B1_S.reverse_copy()
#         C1_E = xpt1__psiMinPF1
#         C1_S = psiMinPF1__LITP
#         C1_W = (LITP.split(C1_S.p[-1])[1]).split(C1_N.p[0], add_split_point = True)[0]
#         C1 = DNL_Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'C1', platePatch = True, plateLocation = location)

#         # ============== Patch A2 ==============
#         A2_N = imidLine__LowerPsiMaxInner.reverse_copy()
#         A2_E = sptrx2_inner__sptrx3_inner.reverse_copy()
#         A2_S = sptrx2imidLine__sptrx2Lower
#         A2_W = A1_E.reverse_copy()
#         A2 = DNL_Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'A2')

#         # ============== Patch B2 ==============
#         B2_N = A2_S.reverse_copy()
#         B2_E = sptrx1__sptrx2_inner.reverse_copy()
#         B2_S = xpt1NW__sptrx1imidLine.reverse_copy()
#         B2_W = B1_E.reverse_copy()
#         B2 = DNL_Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'B2')

#         # ============== Patch C2 ==============
#         C2_N = B2_S.reverse_copy()
#         C2_E = sptrx1__core_inner
#         C2_S = psiMinCore__imidLineCore.reverse_copy()
#         C2_W = xpt1N__psiMinCore.reverse_copy()
#         C2 = DNL_Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'C2')

#         # ============== Patch A3 ==============
#         A3_N = psiMaxInner__imidLine.reverse_copy()
#         A3_E = xpt2__psiMaxInner.reverse_copy()
#         A3_S = xpt2NE__sptrx2imidLine
#         A3_W = A2_E.reverse_copy()
#         A3 = DNL_Patch([A3_N, A3_E, A3_S, A3_W], patchName = 'A3')

#         # ============== Patch B3 ==============
#         B3_N = A3_S.reverse_copy()
#         B3_E = xpt2__sptrx1
#         B3_S = sptrx1imidLine__topLine.reverse_copy()
#         B3_W = B2_E.reverse_copy()
#         B3 = DNL_Patch([B3_N, B3_E, B3_S, B3_W], patchName = 'B3')

#         # ============== Patch C3 ==============
#         C3_N = B3_S.reverse_copy()
#         C3_E = sptrx1__core_top
#         C3_S = imidLineCore__topLine.reverse_copy()
#         C3_W = C2_E.reverse_copy()
#         C3 = DNL_Patch([C3_N, C3_E, C3_S, C3_W], patchName = 'C3')

#         # ============== Patch A4 ==============
#         location = 'E'
#         A4_N = psiMaxInner__UITP
#         A4_S = xpt2__UITP.reverse_copy()
#         A4_E = (UITP.split(A4_N.p[-1])[1]).split(A4_S.p[0], add_split_point = True)[0]
#         A4_W = A3_E.reverse()
#         A4 = DNL_Patch([A4_N, A4_E, A4_S, A4_W], patchName = 'A4', plateLocation = location)

#         # ============== Patch B4 ==============
#         location = 'E'
#         B4_N = A4_S.reverse_copy()
#         B4_S = psiMinPF2_A__UITP.reverse_copy()
#         B4_E = (UITP.split(B4_N.p[-1]))[1].split(B4_S.p[0], add_split_point = True)[0]
#         B4_W = xpt2__psiMinPF2_A.reverse_copy()
#         B4 = DNL_Patch([B4_N, B4_E, B4_S, B4_W], patchName = 'B4', plateLocation = location)
#         # ============== Patch C4 ==============
#         location = 'E'
#         C4_N = B4_S.reverse_copy()
#         C4_S = psiMinPF2_B__UITP.reverse_copy()
#         C4_E = (UITP.split(C4_N.p[-1])[1]).split(C4_S.p[0], add_split_point = True)[0]
#         C4_W = xpt2__psiMinPF2_B.reverse_copy()
#         C4 = DNL_Patch([C4_N, C4_E, C4_S, C4_W], patchName = 'C4', plateLocation = location)
#         # ============== Patch A5 ==============
#         location = 'W'
#         A5_N = psiMaxOuter__UOTP.reverse_copy()
#         A5_S = xpt2__UOTP
#         A5_E = xpt2__psiMaxOuter.reverse_copy()
#         A5_W = (UOTP.split(A5_S.p[-1])[1]).split(A5_N.p[0], add_split_point = True)[0]
#         A5 = DNL_Patch([A5_N, A5_E, A5_S, A5_W], patchName = 'A5', plateLocation = location)

#         # ============== Patch B5 ==============
#         location = 'W'
#         B5_N = A5_S.reverse_copy()
#         B5_S = psiMinPF2_A__UOTP
#         B5_E = xpt2__psiMinPF2_A
#         B5_W = (UOTP.split(B5_S.p[-1])[1]).split(B5_N.p[0], add_split_point = True)[0]
#         B5 = DNL_Patch([B5_N, B5_E, B5_S, B5_W], patchName = 'B5', plateLocation = location)

#         # ============== Patch C5 ==============
#         location = 'W'
#         C5_N = B5_S.reverse_copy()
#         C5_S = psiMinPF2_B__UOTP
#         C5_E = xpt2__psiMinPF2_B
#         C5_W = (UOTP.split(C5_S.p[-1])[1]).split(C5_N.p[0], add_split_point = True)[0]
#         C5 = DNL_Patch([C5_N, C5_E, C5_S, C5_W], patchName = 'C5', plateLocation = location)

#         # ============== Patch A6 ==============
#         A6_N = psiMaxOuter__omidLine
#         A6_S = xpt2NW__sptrx2omidLine.reverse_copy()
#         A6_E = sptrx2_outer__sptrx3_outer.reverse_copy()
#         A6_W = A5_E.reverse_copy()
#         A6 = DNL_Patch([A6_N, A6_E, A6_S, A6_W], patchName = 'A6')

#         # ============== Patch B6 ==============
#         B6_N = A6_S.reverse_copy()
#         B6_S = sptrx1omidLine__topLine
#         B6_E = sptrx1__sptrx2_outer.reverse_copy()
#         B6_W = B3_E.reverse_copy()
#         B6 = DNL_Patch([B6_N, B6_E, B6_S, B6_W], patchName = 'B6')

#         # ============== Patch C6 ==============
#         C6_N = B6_S.reverse_copy()
#         C6_S = omidLineCore__topLine
#         C6_E = sptrx1__core_outer
#         C6_W = C3_E.reverse_copy()
#         C6 = DNL_Patch([C6_N, C6_E, C6_S, C6_W], patchName = 'C6')

#         # ============== Patch A7 ==============
#         omidLine__LowerPsiMaxOuter, LowerPsiMaxOuter__LOTP = omidLine__LOTP.split(sptrx2__LowerPsiMaxOuter.p[-1], add_split_point = True)
#         sptrx2imidLine__sptrx2Lower, sptrx2Lower__LOTP = sptrx2omidLine__LOTP.split(xpt1E__sptrx2Lower.p[-1], add_split_point = True)
#         A7_N = omidLine__LowerPsiMaxOuter
#         A7_S = sptrx2imidLine__sptrx2Lower.reverse_copy()
#         A7_E = sptrx2__LowerPsiMaxOuter.reverse_copy()
#         A7_W = A6_E.reverse_copy()
#         A7 = DNL_Patch([A7_N, A7_E, A7_S, A7_W], patchName = 'A7')

#         # ============== Patch B7 ==============
#         B7_N = A7_S.reverse_copy()
#         B7_S = xpt1NE__sptrx1omidLine
#         B7_E = xpt1E__sptrx2Lower.reverse_copy()
#         B7_W = B6_E.reverse_copy()
#         B7 = DNL_Patch([B7_N, B7_E, B7_S, B7_W], patchName = 'B7')

#         # ============== Patch C7 ==============
#         C7_N = B7_S.reverse_copy()
#         C7_S = psiMinCore__omidLineCore
#         C7_E = C2_W.reverse_copy()
#         C7_W = C6_E.reverse_copy()
#         C7 = DNL_Patch([C7_N, C7_E, C7_S, C7_W], patchName = 'C7')

#         # ============== Patch A8 ==============
#         location = 'E'
#         A8_N = LowerPsiMaxOuter__LOTP
#         A8_S = sptrx2Lower__LOTP.reverse_copy()
#         A8_E = (LOTP.split(A8_N.p[-1])[1]).split(A8_S.p[0], add_split_point = True)[0]
#         A8_W = A7_E.reverse_copy()
#         A8 = DNL_Patch([A8_N, A8_E, A8_S, A8_W], patchName = 'A8', plateLocation = location)

#         # ============== Patch B8 ==============
#         location = 'E'
#         B8_N = A8_S.reverse_copy()
#         B8_S = xpt1__LOTP.reverse_copy()
#         B8_E = (LOTP.split(B8_N.p[-1])[1]).split(B8_S.p[0], add_split_point = True)[0]
#         B8_W = B7_E.reverse_copy()
#         B8 = DNL_Patch([B8_N, B8_E, B8_S, B8_W], patchName = 'B8', plateLocation = location)

#         # ============== Patch C8 ==============
#         location = 'E'
#         C8_N = B8_S.reverse_copy()
#         C8_S = psiMinPF1__LOTP.reverse_copy()
#         C8_E = (LOTP.split(C8_N.p[-1])[1]).split(C8_S.p[0], add_split_point = True)[0]
#         C8_W = C1_E.reverse_copy()
#         C8 = DNL_Patch([C8_N, C8_E, C8_S, C8_W], patchName = 'C8', plateLocation = location)

#         self.patches = [A1, B1, C1, A2, B2, C2, A3, B3, C3, A4, B4, C4, A5, B5, C5, A6, B6, C6, A7, B7, C7, A8, B8, C8]
#         patch_lookup = {}
#         for patch in self.patches:
#             patch_lookup[patch.patchName] = patch
#             patch.plot_border()
#             patch.fill()
#         self.patch_lookup = patch_lookup

#         p = self.patch_lookup
#         self.patch_matrix = [[[None],  [None],  [None],  [None],  [None], [None], [None], [None], [None], [None], [None], [None]], \
#                             [[None], p['A1'], p['A2'], p['A3'], p['A4'], [None], [None], p['A5'], p['A6'], p['A7'], p['A8'], [None]],  \
#                             [[None], p['B1'], p['B2'], p['B3'], p['B4'], [None], [None], p['B5'], p['B6'], p['B7'], p['B8'], [None]],  \
#                             [[None], p['C1'], p['C2'], p['C3'], p['C4'], [None], [None], p['C5'], p['C6'], p['C7'], p['C8'], [None]], \
#                             [[None],  [None],  [None],  [None],  [None], [None], [None], [None], [None], [None], [None], [None]]  \
#                             ]

#         m = self.patch_matrix
#         self.PRIMARY_SOL = m[2][1:4] + m[2][8:11]
#         self.SECONDARY_SOL = m[1][1:5] + m[1][7:11]
#         self.CORE = m[3][2:4] + m[3][8:10]
#         self.PRIMARY_PF = [m[3][1], m[3][-2]]
#         self.SECONDARY_PF = [m[3][4], m[2][4], m[3][7], m[2][7]]


