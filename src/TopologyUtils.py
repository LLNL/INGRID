import numpy as np
import matplotlib
import pathlib
import inspect
import yaml as yml
try:
    matplotlib.use("TkAgg")
except:
    pass
import matplotlib.pyplot as plt

class TopologyUtils():
    def __init__(self, Ingrid_obj, config):
        self.parent = Ingrid_obj
        self.config = config
        self.settings = Ingrid_obj.settings
        self.plate_data = Ingrid_obj.plate_data
        self.PatchTagMap = self.parent.GetPatchTagMap(config=config)
        self.eq = Ingrid_obj.eq
        self.efit_psi = Ingrid_obj.efit_psi
        self.psi_norm = Ingrid_obj.psi_norm
        self.eq = Ingrid_obj.eq
        self.CurrentListPatch={}
        self.Verbose=False
        self.plate_data = Ingrid_obj.plate_data
        self.CorrectDistortion={}

    def RefreshSettings(self):
        self.settings = self.parent.settings

    def patch_diagram(self, fig=None, ax=None):
        """
        Generates the patch diagram for a given configuration.
        @author: watkins35, garcia299
        """

        # colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
        #           'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
        #           'seagreen', 'firebrick', 'saddlebrown', 'dodgerblue',
        #           'violet', 'magenta', 'olivedrab', 'darkorange', 'mediumslateblue',
        #           'palevioletred', 'yellow', 'lightpink', 'plum', 'lawngreen',
        #           'tan', 'hotpink', 'lightgray', 'darkcyan', 'navy']
        colors = {'1' : 'blue', '2' : 'blue', '3' : 'blue'}
        alpha = {'3' : 1.0, '2' : 0.45, '1' : 0.15}

        f=fig if fig else plt.figure('INGRID: Patch Map', figsize=(6, 10))
        a=ax if ax else f.subplots(1,1)
        a.set_xlim([self.efit_psi.rmin, self.efit_psi.rmax])
        a.set_ylim([self.efit_psi.zmin, self.efit_psi.zmax])
        a.set_aspect('equal', adjustable='box')

        a.set_xlabel('R')
        a.set_ylabel('Z')
        f.suptitle(f'{self.config} Patch Diagram')

        for i, patch in enumerate(self.patches.values()):
            patch.plot_border(color='midnightblue', ax=a)
            patch.fill(colors[patch.get_tag()[-1]], ax=a, alpha=alpha[patch.get_tag()[-1]])
            patch.color=colors[patch.get_tag()[-1]]
        # ax.legend()
        f.show()


    def grid_diagram(self,ax=None):
        """
        Create Grid matplotlib figure for an SNL object.
        @author: watkins35, garcia299
        """
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
        plt.title(f'{self.config} Subgrid')
        plt.show()

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

    def get_config(self):
        """
        Returns a string indicating whether the
        """
        return self.config

    def concat_grid(self):
        """
        Concatenate all local grids on individual patches into a single
        array with branch cuts
        Parameters:
        ----------
            config : str
                Type of SNL grid to concat.
        """
        patch_matrix = self.patch_matrix

        for patch in self.patches.values():
            patch.npol = len(patch.cell_grid[0]) + 1
            patch.nrad = len(patch.cell_grid) + 1

        # Total number of poloidal indices in all subgrids.
        np_total1 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][1:5]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total1 = int(np.sum([patch[1].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        # Total number of poloidal indices in all subgrids.
        np_total2 = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:11]])) + 2

        # Total number of radial indices in all subgrids.
        nr_total2 = int(np.sum([patch[7].nrad - 1 for patch in patch_matrix[1:4]])) + 2

        rm1 = np.zeros((np_total1, nr_total1, 5), order = 'F')
        zm1  = np.zeros((np_total1, nr_total1, 5), order = 'F')
        rm2 = np.zeros((np_total2, nr_total2, 5), order = 'F')
        zm2  = np.zeros((np_total2, nr_total2, 5), order = 'F')

        ixcell = 0
        jycell = 0

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')
        for ixp in range(1, 5):

            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

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
                            rm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm1[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Iterate over all the patches in our DNL configuration (we exclude guard cells denoted by '[None]')

        ixcell = 0
        jycell = 0

        for ixp in range(7, 11):

            nr_sum = 0
            for jyp in range(1, 4):
                # Point to the current patch we are operating on.
                local_patch = patch_matrix[jyp][ixp]

                if local_patch == [None]:
                    continue

                nr_sum += local_patch.nrad - 1

                # Access the grid that is contained within this local_patch.
                # ixl - number of poloidal cells in the patch.
                for ixl in range(len(local_patch.cell_grid[0])):
                    # jyl - number of radial cells in the patch
                    for jyl in range(len(local_patch.cell_grid)):

                        ixcell = int(np.sum([patch.npol - 1 for patch in patch_matrix[1][7:ixp+1]])) \
                                - len(local_patch.cell_grid[0]) + ixl + 1

                        jycell = nr_sum - (local_patch.nrad - 1) + jyl + 1

                        ind = 0
                        for coor in ['CENTER', 'SW', 'SE', 'NW', 'NE']:
                            rm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].x
                            zm2[ixcell][jycell][ind] = local_patch.cell_grid[jyl][ixl].vertices[coor].y
                            ind += 1
                            if Verbose: print('Populated RM/ZM entry ({}, {}) by accessing cell ({}, {}) from patch "{}"'.format(ixcell, jycell, jyl, ixl, local_patch.patchName))

        # Flip indices into gridue format.
        for i in range(len(rm1)):
            rm1[i] = rm1[i][::-1]
        for i in range(len(zm1)):
            zm1[i] = zm1[i][::-1]
        for i in range(len(rm2)):
            rm2[i] = rm2[i][::-1]
        for i in range(len(zm2)):
            zm2[i] = zm2[i][::-1]

        # Add guard cells to the concatenated grid.
        ixrb1 = len(rm1) - 2
        ixlb1 = 0
        ixrb2 = len(rm2) - 2
        ixlb2 = 0

        rm1 = self.add_guardc(rm1, ixlb1, ixrb1)
        zm1 = self.add_guardc(zm1, ixlb1, ixrb1)
        rm2 = self.add_guardc(rm2, ixlb2, ixrb2)
        zm2 = self.add_guardc(zm2, ixlb2, ixrb2)

        self.rm = np.concatenate((rm1, rm2))
        self.zm = np.concatenate((zm1, zm2))

        try:
            debug = self.settings['DEBUG']['visual']['gridue']
        except:
            debug = False

        if debug:
            self.animate_grid()

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
        if Patch in self.PatchGroup['SOL']:
            try:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_sol']
            except:
                _radial_f = self.settings['grid_params']['grid_generation']['radial_f_primary_sol']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('SOL radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.PatchGroup['CORE']:
            _radial_f = self.settings['grid_params']['grid_generation']['radial_f_core']
            valid_function=self.CheckFunction(_radial_f,Enforce)
            if Verbose: print('Core radial transformation: "{}"'.format(_radial_f))
        elif Patch in self.PatchGroup['PF']:
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

        if Patch in self.PatchGroup['SOL']:
            _poloidal_f = self.settings['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('SOL poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.PatchGroup['CORE']:
            _poloidal_f = self.settings['grid_params']['grid_generation']['poloidal_f']
            valid_function=self.CheckFunction(_poloidal_f,Enforce)
            if Verbose: print('Core poloidal transformation: "{}"'.format(_poloidal_f))
        elif Patch in self.PatchGroup['PF']:
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

        if Patch.platePatch:
            if Patch.plateLocation == 'W':
                np_cells = self.settings['target_plates']['plate_W1']['np_local']
            elif Patch.plateLocation == 'E':
                np_cells = self.settings['target_plates']['plate_E1']['np_local']

        return (nr_cells,np_cells)


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
                (nr_cells,np_cells)=self.settings['grid_params']['grid_generation']['nr_global'], self.settings['grid_params']['grid_generation']['np_global']#self.GetNpoints(patch, Enforce=Enforce)
                (_radial_f,_poloidal_f)= lambda x: x, lambda x: x # self.GetFunctions(patch,ExtraSettings=ExtraSettings,Enforce=Enforce)
                print('>>> Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(name, np_cells, nr_cells, inspect.getsource(_poloidal_f), inspect.getsource(_radial_f)))
                patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f,verbose = verbose, visual = visual,ShowVertices=ShowVertices,OptionTrace=OptionTrace)
                self.AdjustPatch(patch)
                patch.plot_subgrid()
                self.CurrentListPatch[name] = patch

    # def construct_grid(self, np_cells = 1, nr_cells = 1,Verbose=False,ShowVertices=False,RestartScratch=False,OptionTrace='theta',ExtraSettings={},ListPatches='all', Enforce=True):

        # Straighten up East and West segments of our patches,
        # Plot borders and fill patches.
        
        # self.RefreshSettings()

        # if Verbose: print('Construct Grid')
        # try:
        #     visual = self.settings['DEBUG']['visual']['subgrid']
        # except:
        #     visual = False
        # try:
        #     verbose = self.settings['DEBUG']['verbose']['grid_generation']
        # except:
        #     verbose = False

        # verbose=Verbose or verbose
        
            
        # print('>>> Patches:', [k for k in self.patches.keys()])
        # if RestartScratch:
        #     self.CurrentListPatch={}
    
        # for name, patch in self.patches.items():
            
        #     if self.CorrectDistortion.get(name) is not None:
        #        patch.CorrectDistortion=self.CorrectDistortion.get(name)
        #     elif self.CorrectDistortion.get('all') is not None:
        #         patch.CorrectDistortion=self.CorrectDistortion.get('all')
        #     else:
        #         patch.CorrectDistortion={'Active':False}
        #     if (ListPatches=='all' and patch not in self.CurrentListPatch) or (ListPatches!='all' and name in ListPatches):
        #         self.SetPatchBoundaryPoints(patch)
        #         (nr_cells,np_cells)=self.GetNpoints(patch, Enforce=Enforce)
        #         (_radial_f,_poloidal_f)=self.GetFunctions(patch,ExtraSettings=ExtraSettings,Enforce=Enforce)
        #         print('>>> Making subgrid in patch:{} with np={},nr={},fp={},fr={}'.format(name, np_cells, nr_cells, inspect.getsource(_poloidal_f), inspect.getsource(_radial_f)))
        #         patch.make_subgrid(self, np_cells, nr_cells, _poloidal_f=_poloidal_f,_radial_f=_radial_f,verbose = verbose, visual = visual,ShowVertices=ShowVertices,OptionTrace=OptionTrace)
        #         self.AdjustPatch(patch)
        #         patch.plot_subgrid()
        #         self.CurrentListPatch[name] = patch


        # if all(['cell_grid' in patch.__dict__ for patch in self.patches.values()]):
        #     self.concat_grid()
        #     self.set_gridue()

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

    def CheckPatches(self,verbose=False):
        for name, patch in self.patches.items():
            if patch.platePatch:
                print(' # Checking patch: ', name)
                patch.CheckPatch(self)

    def GroupPatches(self):
        # p = self.patches
        # self.PatchGroup = {'SOL' : [], 
        # 'CORE' : (p['ICB'], p['ICT'], p['OCT'], p['OCB']), 
        # 'PF' : (p['IPF'], p['OPF'])}
        pass

    def SetupPatchMatrix(self):
        # p = self.patches
        # self.patch_matrix = [[[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]], \
        #                 [[None], p['IDL'], p['ISB'], p['IST'], p['OST'], p['OSB'], p['ODL'], [None]], \
        #                 [[None], p['IPF'], p['ICB'], p['ICT'], p['OCT'], p['OCB'], p['OPF'], [None]], \
        #                 [[None],   [None],   [None],   [None],   [None],   [None],   [None], [None]]  \
        #                 ]
        pass
