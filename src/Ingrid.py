#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:17:21 2019

@author: watkins35
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import f90nml


class Ingrid:
    """ An interactive grid generator for edge plasmas in a tokamak
    Accepts a dictionary generated from a namelist file that contains
    the paramters.
    
    Parameters
    ----------
    nml : dict
        params dictionary object with shape
        files :
            geqdsk, itp, otp
        grid params :
            psi_max, psi_min_core, psi_min_pf, Rmagx, Zmagx, Rxpt, Zxpt
    
    """

    def __init__(self, nml):
        self.files = nml['files']
        self.grid_params = nml['grid_params']
        # TODO: include this in the file for parameters
        self.grid_params['step_ratio'] = 0.02

        print('Welcome to Ingrid!\n')
        
    def setup(self):
        """ Add the magnetic axis and the x-point """
        self.add_magx(self.grid_params['Rmagx'], self.grid_params['Zmagx'])
        self.add_xpt1(self.grid_params['Rxpt'], self.grid_params['Zxpt'])

    # TODO update this method to be able to change parameters.
    def set_param(self, key=None, value=None):
        """ Sets the parameters that we want, without
        adding new and unwanted options.
        User can pass in new values directly in,
        or can be prompted.
        
        Parameters
        ----------
        key : str, optional
            keyword for the dicitonary object. Specify which value we
            wish to change. Will be prompted for it if is None.
        value : varies, optional
            value associated with a key. Must match the type of the
            value currently tied to the key.
        
        
        """
        def check_key(self, key=None, first_time=True):
            while True:
                try:
                    if not first_time:
                        # defaults to a string
                        key = raw_input('Enter Key: ')
                    if key not in self._grid_params:
                        raise ValueError
                    print('Key: "{}" is in Params'.format(key))
                    break

                except ValueError:
                    print('That key is not in Grid Parameters.')
                    first_time = False
            return key

        def check_value(self, key, value=None, first_time=True):
            while True:
                try:
                    if not first_time or value is None:
                        # Need to preserve the type of the value entered
                        value = input('Enter Value: ')
                    if not isinstance(self._grid_params[key],
                                      type(value)):
                        raise TypeError
                    self._grid_params[key] = value
                    print('New value recorded as "{}".'.format(value))
                    break

                except TypeError:
                    print('That value is the wrong type.')
                    first_time = False

        # in case the user didn't enter the key
        if key is None:
            key = check_key(first_time=False)
            check_value(key, first_time=False)

        # check if the key the user entered is in grid params
        else:
            key = check_key(key)
            check_value(key, value)

    # TODO: update this method.
    def get_param(self, key=None):
        """ Returns the value associated with a key. """
        first_time = True
        while True:
            try:
                if not first_time or key is None:
                    key = raw_input('Enter Key: ')
                if key not in self._grid_params:
                    raise ValueError
                return self._grid_params[key]

            except ValueError:
                print('Key not a grid parameter.')
                first_time = False

    def print_params(self):
        """ Lists the parameters that the user should know about. """
        for key, value in self._grid_params.items():
            print(key, value)


    def OMFIT_read_psi(self):
        """
        Python class to read the psi data in from an ascii file.
        Saves the boundary information and generated efit_data instance
        """
        
        from OMFITgeqdsk import OMFITgeqdsk
        from Interpol.Setup_Grid_Data import Efit_Data

        g = OMFITgeqdsk(self.files['geqdsk'])

        nxefit = g['NW']
        nyefit = g['NH']
        rdim = g['RDIM']
        zdim = g['ZDIM']
        zmid = g['ZMID']
        rgrid1 = g['RLEFT']

        psi = g['PSIRZ'].T

        # TODO: possibly use the limiters to determine where the strke plates
        # are located. rlim and zlim contain lists of the r and z coordinates
        # for the limiter surrounding the plasma
        # rlim = g['RLIM']  # limiter - similar to strike plates
        # zlim = g['ZLIM']


        # calc array for r and z
        rmin = rgrid1
        rmax = rmin + rdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        # reproduce efit grid
        self.efit_psi = Efit_Data(rmin, rmax, nxefit,
                                  zmin, zmax, nyefit,
                                  name='Efit Data')
        self.efit_psi.set_v(psi)


    def import_psi_data(self):
        """ Same as OMFIT_read_psi, but calls a fortran function to
        read in the psi data. Saves and creates the same instance.
        """
        import Efit.efit as efit
        from Interpol.Setup_Grid_Data import Efit_Data
        efit.readg(self.files['geqdsk'])  # defines the variables
        # extract dimensions and boundaries
        nxefit, nyefit, nbdry, nlim = efit.get_nxy()

        fold, fpol, pres, qpsi, \
            rbdry, zbdry, xlim, ylim, xdim, zdim, \
            rcentr, rgrid1, zmid, Rmagx, Zmagx, \
            simagx, sibdry, bcentr = efit.get_psi(nxefit, nyefit,
                                                  nbdry, nlim)
        # calc array for r and z
        rmin = rgrid1
        rmax = rmin + xdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        # reproduce efit grid
        self.efit_psi = Efit_Data(rmin, rmax, nxefit,
                                  zmin, zmax, nyefit,
                                  name='Efit Data')
        self.efit_psi.set_v(fold)  # this is psi


    def read_target_plate(self):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        The lines to define target plates end up with the shape ((x,y),(x,y)).
        These files read can contain more complicated lines than this.
        """

        self.itp = []
        with open(self.files['itp']) as f:
            for line in f:
                point = line.strip()
                if point.startswith('#'):
                    # move along to the next iteration
                    # this simulates a comment
                    continue
                x = float(point.split(',')[0])
                y = float(point.split(',')[1])
                self.itp.append((x, y))

        self.otp = []
        with open(self.files['otp']) as f:
            for line in f:
                point = line.strip()
                if point.startswith('#'):
                    # move along to the next iteration
                    continue
                x = float(point.split(',')[0])
                y = float(point.split(',')[1])
                self.otp.append((x, y))
        print('Using inner target plate', self.itp)
        print('Using outer target plate', self.otp)

    def plot_target_plate(self):
        """ Plots the inner and outer target plates on the current figure """
        
        itp = np.array(self.itp)
        otp = np.array(self.otp)

        plt.plot(itp[:, 0], itp[:, 1], label='itp')
        plt.plot(otp[:, 0], otp[:, 1], label='otp')
        plt.draw()

    def calc_efit_derivs(self):
        """ Calculates the partial derivatives using finite differences.
        Wrapper for the member function in the efit class.
        """
        self.efit_psi.Calculate_PDeriv()

    def plot_efit_data(self):
        """ generates the plot that we will be able to manipulate
        using the root finder """
        self.efit_psi.plot_data()

    def find_roots(self):
        """ displays a plot, and has the user click on an approximate
        zero point. Uses a root finder to adjust to the more exact point.
        Right click to disable.
        """
        # All we need to do is pass in the crude grid with its derivatives.
        # The interpolation for derivatives is handled by our newton functon
        # make plot - click on point - refine null point
        # the roots are saved at self.root_finder.roots
        from Root_Finder import RootFinder
        self.root_finder = RootFinder(self.efit_psi)

    def save_root(self):
        """ Returns the value of the root finder. This will either be the
        root, or the click, depending on the root finding settings.
        """
        return self.root_finder.final_root

    def toggle_root_finder(self):
        """ Activates or deactivates the root finder ability. Enables
        the user to save the location where they last clicked.
        """
        self.root_finder.toggle_root_finding()

    # TODO: remove this function and replace with read in value
    def add_magx(self, r=None, z=None):
        """ Adds the magnetic axis using the coordinates of last point
        that the user has clicked.
        """
        if r is None and z is None:
            self.magx = self.root_finder.final_root
        else:
            self.magx = (r, z)
        print('Added {} as the magnetic axis.'.format(self.magx))

    # TODO: remove this function and replace with read in value
    def add_xpt1(self, r=None, z=None):
        """ Adds the X-point using the coordinates of last point
        that the user has clicked.
        """
        if r is None and z is None:
            self.xpt1 = self.root_finder.final_root
        else:
            self.xpt1 = (r, z)
        print('Added {} as the primary x point.'.format(self.xpt1))

    def calc_psinorm(self):
        """ Uses magx and xpt1 to normalize the psi data. Furthur calculations
        will use this information """
        from Interpol.Setup_Grid_Data import Efit_Data
        # use the same bounds as the efit data
        self.psi_norm = Efit_Data(self.efit_psi.rmin, self.efit_psi.rmax,
                                  self.efit_psi.nr, self.efit_psi.zmin,
                                  self.efit_psi.zmax, self.efit_psi.nz,
                                  name='psi norm')
        psi = self.efit_psi.v
        psi_magx = self.efit_psi.get_psi(self.magx[0], self.magx[1])
        psi_xpt1 = self.efit_psi.get_psi(self.xpt1[0], self.xpt1[1])
        psinorm = (psi - np.full_like(psi, psi_magx))/(psi_xpt1 - psi_magx)
        self.psi_norm.set_v(psinorm)
        self.psi_norm.Calculate_PDeriv()
        self.psi_norm.plot_data()



    def draw_polodal_lines(self):
        """ trace contour lines anywhere you click
        saves the most recent line that was draw.
        To keep track of lines, call the add_line function"""
        from line_tracing import LineTracing
        self.new_line = LineTracing(self.psi_norm, self.grid_params,
                                    option='theta')

    def draw_radial_lines(self):
        """ trace the orthogonal lines to the stream function
        saves most recent line"""
        from line_tracing import LineTracing
        self.new_line = LineTracing(self.psi_norm, self.grid_params,
                                    option='rho', direction='ccw')

    def compute_eq_psi(self):
        """ Initializes the line tracing class for the construction
        of the grid. Necessary to be separate so we can wiat until the
        user clicks on the root.
        """
        from line_tracing import LineTracing
        self.eq = LineTracing(self.psi_norm, self.grid_params,
                              option='xpt_circ')
        self.eq.calc_equal_psi_points(self.xpt1[0], self.xpt1[1])
        self.eq.disconnect()

    def construct_SNL_patches(self):
        """ more general format to construct the grid for SNL using patches
        for theta direction, cw is positive
        for rho direction, ccw is positive, which is away from the magx

        lines are drawn in the form:
        self.eq.draw_line(, , option='', direction='')

        Patch Labeling Key:
            I: Inner
            O: Outer
            DL: Divertor Leg
            PF: Private Flux
            T: Top
            B: Bottom
            S: Scrape Off Layer
            C: Core
        """
        plt.show()
        from geometry import Point, Patch, Line

        xpt = self.eq.eq_psi
        psi_max = self.grid_params['psi_max']
        psi_min_core = self.grid_params['psi_min_core']
        psi_min_pf = self.grid_params['psi_min_pf']


        # IDL ===================================================
        E = self.eq.draw_line(xpt['W'], {'psi': psi_max}, option='rho', direction='ccw')
        N = self.eq.draw_line(E.p[-1], {'line': self.itp}, option='theta', direction='ccw').reverse()
        S = self.eq.draw_line(xpt['SW'], {'line': self.itp}, option='theta', direction='ccw')
        E = Line([N.p[-1], S.p[0]]) # straighten it up
        W = Line([S.p[-1], N.p[0]])
        IDL = Patch([N, E, S, W])
        # Lines are now saved inside of the patch

        # IPF ===================================================
        N = IDL.S.reverse()
        E = self.eq.draw_line(xpt['S'], {'psi': psi_min_pf}, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(E.p[-1], {'line': self.itp}, option='theta', direction='ccw')
        W = Line([S.p[-1], N.p[0]])
        IPF = Patch([N, E, S, W])

        # OPF ===================================================
        N = self.eq.draw_line(xpt['SE'], {'line': self.otp}, option='theta', direction='cw')
        W = IPF.E.reverse()
        S = self.eq.draw_line(W.p[0], {'line': self.otp}, option='theta', direction='cw').reverse()
        E = Line([N.p[-1], S.p[0]])
        OPF = Patch([N, E, S, W])

        # ODL ===================================================
        W = self.eq.draw_line(xpt['E'], {'psi': psi_max}, option='rho', direction='ccw')
        W = Line([W.p[0], W.p[-1]])
        N = self.eq.draw_line(W.p[-1], {'line': self.otp}, option='theta', direction='cw')
        S = OPF.N.reverse()
        E = Line([N.p[-1], S.p[0]])
        ODL = Patch([N, E, S, W])

        # need the mid and top points of the separatrix
        sep = self.eq.draw_line(Point(xpt['NW']), {'point': Point(xpt['NE'])}, option='theta', direction='cw')
        top_index = np.argmax(sep.yval)
        midpoint = np.median(sep.yval)
        mid_index, = np.where(np.abs(sep.yval-midpoint) < 1e-3)
        dp = {}  # defining points
        dp['top'] = Point(sep.xval[top_index], sep.yval[top_index])
        dp['impt'] = Point(sep.xval[mid_index[0]], sep.yval[mid_index[0]])
        dp['ompt'] = Point(sep.xval[mid_index[1]], sep.yval[mid_index[1]])

        # ISB ===================================================
        E = self.eq.draw_line(dp['impt'], {'psi': psi_max}, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(IDL.E.p[0], {'point': E.p[0]}, option='theta', direction='cw')
        S = Line(sep.p[:mid_index[0]+1]).reverse()
        W = IDL.E.reverse()
        ISB = Patch([N, E, S, W])

        # ICB ===================================================
        N = Line(sep.p[:mid_index[0]+1])
        E = self.eq.draw_line(dp['impt'], {'psi': psi_min_core}, option='rho', direction='cw')
        W = self.eq.draw_line(xpt['N'], {'psi': psi_min_core}, option='rho', direction='cw').reverse()
        W = Line([W.p[0], W.p[-1]])
        S = self.eq.draw_line(W.p[0], {'point': E.p[-1]}, option='theta', direction='cw').reverse()
        ICB = Patch([N, E, S, W])

        # IST ===================================================
        E = self.eq.draw_line(dp['top'], {'psi': psi_max}, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(ISB.N.p[-1], {'point': E.p[0]}, option='theta', direction='cw')
        S = Line(sep.p[mid_index[0]:top_index+1]).reverse()
        W = ISB.W.reverse()
        IST = Patch([N, E, S, W])

        # ICT ===================================================
        E = self.eq.draw_line(dp['top'], {'psi': psi_min_core}, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(ICB.S.p[0], {'point': E.p[-1]}, option='theta', direction='cw').reverse()
        N = IST.S.reverse()
        W = ICB.E.reverse()
        ICT = Patch([N, E, S, W])

        # OST ===================================================
        E = self.eq.draw_line(dp['ompt'], {'psi': psi_max}, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(IST.N.p[-1], {'point': E.p[0]}, option='theta', direction='cw')
        S = Line(sep.p[top_index: mid_index[1]+1]).reverse()
        W = IST.W.reverse()
        OST = Patch([N, E, S, W])

        # OCT ===================================================
        E = self.eq.draw_line(dp['ompt'], {'psi': psi_min_core}, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(ICT.E.p[-1], {'point': E.p[-1]}, option='theta', direction='cw').reverse()
        N = Line(sep.p[top_index: mid_index[1]+1])
        W = ICT.E.reverse()
        OCT = Patch([N, E, S, W])

        # OCB ===================================================
        W = OCT.E.reverse()
        N = Line(sep.p[mid_index[1]:])
        S = self.eq.draw_line(W.p[0], {'point': ICB.W.p[0]}, option='theta', direction='cw').reverse()
        E = ICB.W.reverse()
        OCB = Patch([N, E, S, W])

        # OSB ===================================================
        W = OST.E.reverse()
        N = self.eq.draw_line(W.p[-1], {'point': ODL.W.p[-1]}, option='theta', direction='cw')
        S = OCB.N.reverse()
        E = ODL.W.reverse()
        OSB = Patch([N, E, S, W])

        self.patches = [IDL, IPF, OPF, ODL, ISB, ICB, IST, ICT, OST, OCT, OCB, OSB]
        for patch in self.patches:
            patch.plot_border()
            patch.fill()

    def patch_diagram(self):
        """ Generates the patch diagram for a given configuration. """
        
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']

        plt.figure('patches', figsize=(6, 10))
        for i in range(len(self.patches)):
            self.patches[i].plot_border('green')
            self.patches[i].fill(colors[i])
        plt.xlim(self.efit_psi.rmin, self.efit_psi.rmax)
        plt.ylim(self.efit_psi.zmin, self.efit_psi.zmax)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('Example Patch')
        plt.show()

    def refine_patches(self):
        """ break each patch into smaller grids based of psi"""
        # TODO: use scipy.optimize.curve_fit to generate a polynomial
        # fit for the two curved sections of each patch,
        # then use the length to break into even subsections
        # For the horizontal division use the psi levels to define subsections
        print('Refining patches')
        
        # self.patches[0].refine(self)  # test with a single patch
        
        for patch in self.patches:
            patch.refine(self)
    
    def export(self):
        """ Saves the grid as an ascii file """
        # TODO export the points the patches contain, but don't overlap
        # any points
        pass

    def test_interpol(self, option=2, nfine=100, ncrude=10, tag='v'):
        """ Provides a demonstration and test of the bicubic interpolation
        methods used. Wrapper for the real code from another module.
        
        Parameters
        ----------
        option : int, optional
            Specify which of a set of function to use as the test.
            Accepts 1, 11, 12, 13, 2, 3, 4, 5, 51, 52
            See Test_Functions.get_f for more detail
        nfine : int, optional
            Density of the fine grid we will interpolate onto.
        ncrude : int, optional
            Density of the crude grid the sample data will be generated
            onto.
        tag : str, optional
            Specify if it is wanted to test the derivative interpolation
            methods. Accepts 'v', 'vr', 'vz', 'vrz'.
        
        """
        from Interpol.Test_Interpol import test_interpol
        test_interpol(option, nfine, ncrude, tag)



def prep_input():
    """ Interactive mode for Ingrid. Gets the files used, and opens
    a plot of the Efit data. Prompts the user for magx, xpt, and psi levels.
    Saves the data in a namelist file."""
    
    def paws():
        # Helps the code to wait for the user to select points
        raw_input("Press the <ENTER> key to continue...")
    
    nml = {'files': {}, 'grid_params': {}}

    # get files from the user for data, inner and outer strike plates
    # This is an acceptable form - also takes full path names
    #    Geqdsk_file = “../data/SNL/neqdsk”
    #    Inner_plate_file  = “../data/SNL/itp1”
    #    Outer_plate_file  = “../data/SNL/itp1”
    nml['files']['geqdsk'] = raw_input("Enter the geqsdk filename: ").strip()  # remove any whitespace
    nml['files']['itp'] = raw_input("Enter the inner strike plate filename: ").strip()
    nml['files']['otp'] = raw_input("Enter the outer strike plate filename: ").strip()

    # now we can read the data and fine tune our parameters
    grid = Ingrid(nml)
    grid.OMFIT_read_psi()
    grid.read_target_plate()
    grid.calc_efit_derivs()
    grid.plot_efit_data()
    grid.plot_target_plate()
    grid.find_roots()

    # find the null-point (magnetic axis)
    print("Click on the magnetic axis") ##-need blocking here!
    paws()
    magx = grid.save_root()
    nml['grid_params']['Rmagx'] = magx[0]
    nml['grid_params']['Zmagx'] = magx[1]

    # find the x-point
    print("Click on the x-point")
    paws()
    xpt = grid.save_root()
    nml['grid_params']['Rxpt'] = xpt[0]
    nml['grid_params']['Zxpt'] = xpt[1]

    # enter or click on the location of the psi values
    psi_magx = grid.efit_psi.get_psi(magx[0], magx[1])
    psi_xpt1 = grid.efit_psi.get_psi(xpt[0], xpt[1])

    # get max psi level
    max_psi = raw_input("Enter max psi level, or press <ENTER> to click on the location: ")
    if max_psi == '':
        grid.toggle_root_finder()
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(max_psi)
    nml['grid_params']['psi_max'] = psi

    # get min psi value near the core
    core_psi = raw_input("Enter min psi level near the core plasma, or press <ENTER> to click on the location: ")
    if core_psi == '':
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(core_psi)
    nml['grid_params']['psi_min_core'] = psi

    # get min psi value near pf region
    pf_psi = raw_input("Enter min psi level near the private flux region, or press <ENTER> to click on the location: ")
    if pf_psi == '':
        paws()
        x, y = grid.save_root()
        psi_efit = grid.efit_psi.get_psi(x, y)
        psi = (psi_efit - np.full_like(psi_efit, psi_magx))/(psi_xpt1 - psi_magx)
    else:
        psi = float(pf_psi)
    nml['grid_params']['psi_min_pf'] = psi

    outFile = raw_input("Enter output filename [default is 'grid_params.nml']: ").strip()
    if outFile == '':
        outFile = 'grid_params.nml'

    f90nml.write(nml, outFile, force=True)  # force tag is to overwrite the previous file
    print("Saved paramters to '{}'.".format(outFile))
    plt.close('Efit Data') # finish with the raw Psi data


def run():
    """ Reads a namelist file containing the parameters for Ingrid, and
    runs the grid generator for a single null configuration.
    """
    nml_file = raw_input("Enter params filename [params.nml]: ")
    if nml_file == '':
        # default is an example case for snl
        nml_file = 'params.nml'
    
    nml = f90nml.read(nml_file)

    grid = Ingrid(nml)
    grid.setup()
    grid.OMFIT_read_psi()
    grid.read_target_plate()
    grid.calc_psinorm()
    grid.plot_target_plate()
    grid.compute_eq_psi()

    #-construct the patch-map for SNL
    from time import time
    start = time()
    grid.construct_SNL_patches()
    end = time()
    
    # TODO: finish writing the refine patches method.
#    grid.refine_patches()
    
    
    grid.patch_diagram()

    print("Time for grid: {} seconds.".format(end-start))


if __name__ == '__main__':
    prep_input()
