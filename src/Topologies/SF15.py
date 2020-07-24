"""
SF15.py

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


class SF15():
    def __init__(self, Ingrid_obj, config):

        self.parent = Ingrid_obj
        self.config = config
        self.yaml = Ingrid_obj.yaml
        self.plate_data = Ingrid_obj.plate_data

        self.parent.order_target_plates()
        self.PatchTagMap = self.parent.GetPatchTagMap(config='SF15')

        self.eq = Ingrid_obj.eq
        self.efit_psi = Ingrid_obj.efit_psi
        self.psi_norm = Ingrid_obj.psi_norm
        self.eq = Ingrid_obj.eq
        self.CurrentListPatch={}
        self.Verbose=False
        self.plate_data = Ingrid_obj.plate_data
        self.CorrectDistortion={}

    def patch_diagram(self):
        """ Generates the patch diagram for a given configuration. """

        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown', 'c',
                  'm', 'dodgerblue', 'darkorchid', 'crimson',
                  'darkorange', 'lightgreen', 'lightseagreen', 'indigo',
                  'mediumvioletred', 'mistyrose', 'darkolivegreen', 'rebeccapurple']

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
        plt.title('INGRID SNL Subgrid')
        plt.show()

    def get_config(self):
        """
        Returns a string indicating whether the
        """
        return self.config

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
            debug = self.yaml['DEBUG']['visual']['gridue']
        except:
            debug = False

        if debug:
            self.animate_grid()

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

    def construct_patches(self):
        """
        """
        def reorder_limiter(new_start, limiter):
            start_index, = find_split_index(new_start, limiter)
            return limiter
        def limiter_split(start, end, limiter):
            start_index, sls = find_split_index(start, limiter)
            end_index, sls = find_split_index(end, limiter)
            if end_index <= start_index:
                limiter.p = limiter.p[start_index:] + limiter.p[:start_index]
            return limiter

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


        xpt1 = self.eq.NSEW_lookup['xpt1']['coor']
        xpt2 = self.eq.NSEW_lookup['xpt2']['coor']

        magx = np.array([self.yaml['grid_params']['rmagx'] + self.yaml['grid_params']['patch_generation']['rmagx_shift'], \
            self.yaml['grid_params']['zmagx'] + self.yaml['grid_params']['patch_generation']['zmagx_shift']])

        psi_max_west = self.yaml['grid_params']['psi_max_west']
        psi_max_east = self.yaml['grid_params']['psi_max_east']
        psi_min_core = self.yaml['grid_params']['psi_min_core']
        psi_min_pf = self.yaml['grid_params']['psi_min_pf']
        psi_min_pf_2 = self.yaml['grid_params']['psi_pf2']

        if self.yaml['limiter']['use_limiter']:
            WestPlate1 = self.parent.limiter_data.copy()
            WestPlate2 = self.parent.limiter_data.copy()

            EastPlate1 = self.parent.limiter_data.copy()
            EastPlate2 = self.parent.limiter_data.copy()

        else:
            WestPlate1 = Line([Point(i) for i in self.plate_W1])
            WestPlate2 = Line([Point(i) for i in self.plate_W2])

            EastPlate1 = Line([Point(i) for i in self.plate_E1])
            EastPlate2 = Line([Point(i) for i in self.plate_E2])

        # Generate Horizontal Mid-Plane lines
        LHS_Point = Point(magx[0] - 1e6 * np.cos(inner_tilt), magx[1] - 1e6 * np.sin(inner_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(inner_tilt), magx[1] + 1e6 * np.sin(inner_tilt))
        west_midLine = Line([LHS_Point, RHS_Point])
        # inner_midLine.plot()

        LHS_Point = Point(magx[0] - 1e6 * np.cos(outer_tilt), magx[1] - 1e6 * np.sin(outer_tilt))
        RHS_Point = Point(magx[0] + 1e6 * np.cos(outer_tilt), magx[1] + 1e6 * np.sin(outer_tilt))
        east_midLine = Line([LHS_Point, RHS_Point])
        # outer_midLine.plot()

        # Generate Vertical Mid-Plane line
        Lower_Point = Point(magx[0], magx[1] - 1e6)
        Upper_Point = Point(magx[0], magx[1] + 1e6)
        topLine = Line([Lower_Point, Upper_Point])
        # topLine.plot()


        # Tracing primary-separatrix: core-boundary

        # H1_E / B1_W 
        xpt1N__psiMinCore = self.eq.draw_line(xpt1['N'], {'psi' : psi_min_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        H1_E = xpt1N__psiMinCore
        B1_W = H1_E.reverse_copy()

        # B1_N / B2_S
        xpt1NW__west_midLine = self.eq.draw_line(xpt1['NW'], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_N = xpt1NW__west_midLine
        B2_S = B1_N.reverse_copy()

        # H2_S / H1_N
        xpt1NE__east_midLine = self.eq.draw_line(xpt1['NE'], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)

        # C1_N / C2_S
        west_midLine__topLine = self.eq.draw_line(xpt1NW__west_midLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C1_N = west_midLine__topLine
        C2_S = C1_N.reverse_copy()

        # D1_N / D2_S
        east_midLine__topLine = self.eq.draw_line(xpt1NE__east_midLine.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_N = east_midLine__topLine.reverse_copy()
        D2_S = D1_N.reverse_copy()

        # Tracing core: psi-min-core region (rho = 1)
        
        # / B1_S
        psiMinCore__west_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : west_midLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        B1_S = psiMinCore__west_midLine_core.reverse_copy()

        # H1_E1
        psiMinCore__east_midLine_core = self.eq.draw_line(xpt1N__psiMinCore.p[-1], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        # / C1_S
        west_midLine_core__topLine = self.eq.draw_line(psiMinCore__west_midLine_core.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        C1_S = west_midLine_core__topLine.reverse_copy()

        # D1_S
        east_midLine_core__topLine = self.eq.draw_line(psiMinCore__east_midLine_core.p[-1], {'line' : topLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D1_S = east_midLine_core__topLine

        # A1_E / I1_W
        xpt1__psiMinPF1 = self.eq.draw_line(xpt1['S'], {'psi' : psi_min_pf}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        A1_E = xpt1__psiMinPF1
        I1_W = A1_E.reverse_copy()

        # A1_S
        psiMinPF1__WestPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A1_S = psiMinPF1__WestPlate1

        # / I1_S
        psiMinPF1__EastPlate1 = self.eq.draw_line(xpt1__psiMinPF1.p[-1], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        I1_S = psiMinPF1__EastPlate1.reverse_copy()

        # A2_S / A1_N
        xpt1__WestPlate1 = self.eq.draw_line(xpt1['SW'], {'line' : WestPlate1}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        A2_S = xpt1__WestPlate1
        A1_N = A2_S.reverse_copy()

        # I1_N / I2_S
        xpt1__EastPlate1 = self.eq.draw_line(xpt1['SE'], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        I1_N = xpt1__EastPlate1
        I2_S = I1_N.reverse_copy()

        # Drawing the portion of the separatrix found in the double-null configuration

        # E2_E / H2_W
        xpt2N__outer_core = self.eq.draw_line(xpt2['N'], {'line' : xpt1NE__east_midLine}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        E2_E = xpt2N__outer_core
        H2_W = E2_E.reverse_copy()

        # E1_E / H1_W
        xpt2N__outer_core__inner_core = self.eq.draw_line(xpt2N__outer_core.p[-1], \
            {'line' : psiMinCore__east_midLine_core}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)
        E1_E = xpt2N__outer_core__inner_core
        H1_W = E1_E.reverse_copy()
        
        E1_N, H1_N = xpt1NE__east_midLine.reverse_copy().split(xpt2N__outer_core.p[-1], add_split_point=True)
        E2_S = E1_N.reverse_copy()
        H2_S = H1_N.reverse_copy()


        H1_S, E1_S = psiMinCore__east_midLine_core.split(xpt2N__outer_core__inner_core.p[-1], add_split_point=True)
        
        # H2_N__I2_N
        xpt2NW__EastPlate1 = self.eq.draw_line(xpt2['NW'], {'line' : EastPlate1}, \
            option = 'theta', direction = 'cw', show_plot = visual, text = verbose)
        # E3_S / E2_N
        xpt2NE__east_midLine = self.eq.draw_line(xpt2['NE'], {'line' : east_midLine}, \
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        E3_S = xpt2NE__east_midLine
        E2_N = E3_S.reverse_copy()

        # D3_S / D2_N
        xpt2__east_midLine__topLine = self.eq.draw_line(xpt2NE__east_midLine.p[-1], {'line' : topLine}, 
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        D3_S = xpt2__east_midLine__topLine
        D2_N = D3_S.reverse_copy()

        # C3_S / C2_N
        xpt2__topLine__west_midLine = self.eq.draw_line(xpt2__east_midLine__topLine.p[-1], {'line' : west_midLine},
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)
        C3_S = xpt2__topLine__west_midLine
        C2_N = C3_S.reverse_copy()

        xpt2__west_midLine__WestPlate1 = self.eq.draw_line(xpt2__topLine__west_midLine.p[-1], {'line' : WestPlate1},
            option = 'theta', direction = 'ccw', show_plot = visual, text = verbose)


        B2_W = self.eq.draw_line(xpt1['W'], {'line' : xpt2__west_midLine__WestPlate1},
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        A2_E = B2_W.reverse_copy()

        I2_W = self.eq.draw_line(xpt1['E'], {'line' : xpt2NW__EastPlate1},
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        H2_E = I2_W.reverse_copy()

        # / A3_S, B3_S
        A2_N, B2_N = xpt2__west_midLine__WestPlate1.reverse_copy().split(B2_W.p[-1], add_split_point=True)
        A3_S = A2_N.reverse_copy()
        B3_S = B2_N.reverse_copy()

        # / H3_S, I3_S
        H2_N, I2_N = xpt2NW__EastPlate1.split(I2_W.p[-1], add_split_point=True)
        H3_S = H2_N.reverse_copy()
        I3_S = I2_N.reverse_copy()


        xpt2__psiMinPF2 = self.eq.draw_line(xpt2['S'], {'psi' : psi_min_pf_2}, \
            option = 'rho', direction = 'cw', show_plot = visual, text = verbose)

        F1_W, F2_W = xpt2__psiMinPF2.reverse_copy().split(xpt2__psiMinPF2.p[len(xpt2__psiMinPF2.p)//2], add_split_point = True)
        G1_E = F1_W.reverse_copy()
        G2_E = F2_W.reverse_copy()

        F1_S = self.eq.draw_line(F1_W.p[0], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose).reverse_copy()
        G1_S = self.eq.draw_line(F1_W.p[0], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose)
        
        F1_N = self.eq.draw_line(F2_W.p[0], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        F2_S = F1_N.reverse_copy()

        G1_N = self.eq.draw_line(F2_W.p[0], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        G2_S = G1_N.reverse_copy()

        F2_N = self.eq.draw_line(xpt2['SE'], {'line' : EastPlate2}, option = 'theta', direction = 'cw',
            show_plot = visual, text = verbose)
        F3_S = F2_N.reverse_copy()

        G2_N = self.eq.draw_line(xpt2['SW'], {'line' : WestPlate2}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        G3_S = G2_N.reverse_copy()


        H3_W = self.eq.draw_line(xpt2['W'], {'psi' : psi_max_east}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        G3_E = H3_W.reverse_copy()

        H3_W__EastPlate1 = self.eq.draw_line(H3_W.p[-1], {'line' : EastPlate1}, option = 'theta',
            direction = 'cw', show_plot = visual, text = verbose)

        I3_W = self.eq.draw_line(I2_W.p[-1], {'line' : H3_W__EastPlate1}, 
            option = 'rho', direction = 'ccw', show_plot = visual, text = verbose)
        H3_E = I3_W.reverse_copy()
        
        H3_N, I3_N = H3_W__EastPlate1.split(I3_W.p[-1], add_split_point=True)

        G3_N = self.eq.draw_line(H3_W.p[-1], {'line' : WestPlate2}, option = 'theta', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()

        F3_W = self.eq.draw_line(xpt2['E'], {'psi' : psi_max_west}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        E3_E = F3_W.reverse_copy()

        F3_N = self.eq.draw_line(F3_W.p[-1], {'line' : EastPlate2}, option = 'theta', direction = 'cw', 
            show_plot = visual, text = verbose)
        E3_N = self.eq.draw_line(F3_W.p[-1], {'line' : east_midLine}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()
        D3_N = self.eq.draw_line(E3_N.p[0], {'line' : topLine}, option = 'theta', direction = 'ccw', 
            show_plot = visual, text = verbose).reverse_copy()
        C3_N = self.eq.draw_line(D3_N.p[0], {'line' : west_midLine}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose).reverse_copy()

        west_midLine__WestPlate1 = self.eq.draw_line(C3_N.p[0], {'line' : WestPlate1}, option = 'theta', direction = 'ccw',
            show_plot = visual, text = verbose)
        B3_W = self.eq.draw_line(B2_W.p[-1], {'line' : west_midLine__WestPlate1}, option = 'rho', direction = 'ccw',
            show_plot = visual, text = verbose)
        A3_E = B3_W.reverse_copy()

        A3_N, B3_N = west_midLine__WestPlate1.reverse_copy().split(B3_W.p[-1], add_split_point=True)

        B1_E = self.eq.draw_line(B1_N.p[-1], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C1_W = B1_E.reverse_copy()

        B2_E = self.eq.draw_line(B2_N.p[-1], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C2_W = B2_E.reverse_copy()

        B3_E = self.eq.draw_line(B3_N.p[-1], {'psi_horizontal' : Point(xpt2['center']).psi(self)}, option = 'z_const', 
            direction = 'cw', show_plot = visual, text = verbose)
        C3_W = B3_E.reverse_copy()

        C1_E = self.eq.draw_line(C1_N.p[-1], {'psi_vertical' : psi_min_core}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D1_W = C1_E.reverse_copy()

        C2_E = self.eq.draw_line(C2_N.p[-1], {'psi_vertical' : 1.0}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D2_W = C2_E.reverse_copy()

        C3_E = self.eq.draw_line(C3_N.p[-1], {'psi_vertical' : Point(xpt2['center']).psi(self)}, option = 'r_const', 
            direction = 'ccw', show_plot = visual, text = verbose)
        D3_W = C3_E.reverse_copy()

        E1_W = self.eq.draw_line(E1_N.p[0], {'psi_horizontal' : psi_min_core}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D1_E = E1_W.reverse_copy()

        E2_W = self.eq.draw_line(E2_N.p[0], {'psi_horizontal' : 1.0}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D2_E = E2_W.reverse_copy()

        E3_W = self.eq.draw_line(E3_N.p[0], {'psi_horizontal' : Point(xpt2['center']).psi(self)}, option = 'z_const', 
            direction = 'ccw', show_plot = visual, text = verbose).reverse_copy()
        D3_E = E3_W.reverse_copy()

        A1_W = (WestPlate1.split(A1_S.p[-1])[1]).split(A1_N.p[0], add_split_point = True)[0]
        A2_W = (WestPlate1.split(A2_S.p[-1])[1]).split(A2_N.p[0], add_split_point = True)[0]
        A3_W = (WestPlate1.split(A3_S.p[-1])[1]).split(A3_N.p[0], add_split_point = True)[0]

        G1_W = (WestPlate2.split(G1_S.p[-1])[1]).split(G1_N.p[0], add_split_point = True)[0]
        G2_W = (WestPlate2.split(G2_S.p[-1])[1]).split(G2_N.p[0], add_split_point = True)[0]
        G3_W = (WestPlate2.split(G3_S.p[-1])[1]).split(G3_N.p[0], add_split_point = True)[0]

        I1_E = (EastPlate1.split(I1_N.p[-1])[1]).split(I1_S.p[0], add_split_point = True)[0]
        I2_E = (EastPlate1.split(I2_N.p[-1])[1]).split(I2_S.p[0], add_split_point = True)[0]
        I3_E = (EastPlate1.split(I3_N.p[-1])[1]).split(I3_S.p[0], add_split_point = True)[0]

        F1_E = (EastPlate2.split(F1_N.p[-1])[1]).split(F1_S.p[0], add_split_point = True)[0]
        F2_E = (EastPlate2.split(F2_N.p[-1])[1]).split(F2_S.p[0], add_split_point = True)[0]
        F3_E = (EastPlate2.split(F3_N.p[-1])[1]).split(F3_S.p[0], add_split_point = True)[0]

        # ============== Patch A1 ==============
        A1 = Patch([A1_N, A1_E, A1_S, A1_W], patchName = 'A1', platePatch = True, plateLocation = 'W')
        # ============== Patch A2 ==============
        A2 = Patch([A2_N, A2_E, A2_S, A2_W], patchName = 'A2', platePatch = True, plateLocation = 'W')
        # ============== Patch A3 ==============
        A3 = Patch([A3_N, A3_E, A3_S, A3_W], patchName = 'A3', platePatch = True, plateLocation = 'W')


        # ============== Patch B1 ==============
        B1 = Patch([B1_N, B1_E, B1_S, B1_W], patchName = 'B1')
        # ============== Patch B2 ==============
        B2 = Patch([B2_N, B2_E, B2_S, B2_W], patchName = 'B2')
        # ============== Patch B3 ==============
        B3 = Patch([B3_N, B3_E, B3_S, B3_W], patchName = 'B3')

        # ============== Patch C1 ==============
        C1 = Patch([C1_N, C1_E, C1_S, C1_W], patchName = 'C1')
        # ============== Patch C2 ==============
        C2 = Patch([C2_N, C2_E, C2_S, C2_W], patchName = 'C2')
        # ============== Patch C3 ==============
        C3 = Patch([C3_N, C3_E, C3_S, C3_W], patchName = 'C3')

        # ============== Patch D1 ==============
        D1 = Patch([D1_N, D1_E, D1_S, D1_W], patchName = 'D1')
        # ============== Patch D2 ==============
        D2 = Patch([D2_N, D2_E, D2_S, D2_W], patchName = 'D2')
        # ============== Patch D3 ==============
        D3 = Patch([D3_N, D3_E, D3_S, D3_W], patchName = 'D3')

        # ============== Patch E1 ==============
        E1 = Patch([E1_N, E1_E, E1_S, E1_W], patchName = 'E1')
        # ============== Patch E2 ==============
        E2 = Patch([E2_N, E2_E, E2_S, E2_W], patchName = 'E2')
        # ============== Patch E3 ==============
        E3 = Patch([E3_N, E3_E, E3_S, E3_W], patchName = 'E3')

        # ============== Patch F1 ==============
        F1 = Patch([F1_N, F1_E, F1_S, F1_W], patchName = 'F1')
        # ============== Patch F2 ==============
        F2 = Patch([F2_N, F2_E, F2_S, F2_W], patchName = 'F2')
        # ============== Patch F3 ==============
        F3 = Patch([F3_N, F3_E, F3_S, F3_W], patchName = 'F3')

        # ============== Patch G1 ==============
        G1 = Patch([G1_N, G1_E, G1_S, G1_W], patchName = 'G1', platePatch = True, plateLocation = 'W')
        # ============== Patch G2 ==============
        G2 = Patch([G2_N, G2_E, G2_S, G2_W], patchName = 'G2', platePatch = True, plateLocation = 'W')
        # ============== Patch G3 ==============
        G3 = Patch([G3_N, G3_E, G3_S, G3_W], patchName = 'G3', platePatch = True, plateLocation = 'W')

        # ============== Patch H1 ==============
        H1 = Patch([H1_N, H1_E, H1_S, H1_W], patchName = 'H1')
        # ============== Patch H2 ==============
        H2 = Patch([H2_N, H2_E, H2_S, H2_W], patchName = 'H2')
        # ============== Patch H3 ==============
        H3 = Patch([H3_N, H3_E, H3_S, H3_W], patchName = 'H3')

        # ============== Patch I1 ==============
        I1 = Patch([I1_N, I1_E, I1_S, I1_W], patchName = 'I1', platePatch = True, plateLocation = 'E')
        # ============== Patch I2 ==============
        I2 = Patch([I2_N, I2_E, I2_S, I2_W], patchName = 'I2', platePatch = True, plateLocation = 'E')
        # ============== Patch I3 ==============
        I3 = Patch([I3_N, I3_E, I3_S, I3_W], patchName = 'I3', platePatch = True, plateLocation = 'E')

        patches = [A3, A2, A1, B3, B2, B1, C3, C2, C1, D3, D2, D1, E3, E2, E1,
                   F3, F2, F1, G3, G2, G1, H3, H2, H1, I3, I2, I1]

        self.patches = {}
        for patch in patches:
            patch.parent = self
            patch.NameTagMap = self.PatchTagMap
            self.patches[patch.patchName] = patch

    def patch_diagram(self):
        """
        Generates the patch diagram for a given configuration.
        @author: watkins35, garcia299
        """

        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown', 'dodgerblue',
                  'violet', 'magenta', 'olivedrab', 'darkorange', 'mediumslateblue',
                  'palevioletred', 'yellow', 'lightpink', 'plum', 'lawngreen',
                  'tan', 'hotpink', 'lightgray', 'darkcyan', 'navy']


        self.FigPatch=plt.figure('INGRID: Patch Map', figsize=(6, 10))
        self.FigPatch.clf()
        ax=self.FigPatch.subplots(1,1)
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        ax.set_aspect('equal', adjustable='box')

        ax.set_xlabel('R')
        ax.set_ylabel('Z')
        self.FigPatch.suptitle('DNL Patch Diagram: ' + self.config)

        for i, patch in enumerate(self.patches.values()):
            patch.plot_border('green')
            patch.fill(colors[i])
            patch.color=colors[i]
        plt.show()

    def CheckPatches(self,verbose=False):
        for name, patch in self.patches.items():
            if patch.platePatch:
                print(' # Checking patch: ', name)
                patch.CheckPatch(self)