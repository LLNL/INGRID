#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:17:21 2019

@author: watkins35
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt


class Ingrid:
    """ An interactive grid generator """
    def __init__(self,
                 option='bicubic',
                 epsilon=1e-9,
                 rmagx=0.0,
                 rmagy=0.0,
                 zmagx=0.0,
                 zmagy=0.0,
                 manual=True,
                 step_ratio=0.02,
                 psi_max=1.1,
                 psi_min=0.9):

        print('Welcome to Ingrid!\n')
        self._grid_params = {'option': option,
                             'epsilon': epsilon,
                             'rmagx': rmagx,
                             'rmagy': rmagy,
                             'zmagx': zmagx,
                             'zmagy': zmagy,
                             'manual': manual,
                             'step_ratio': step_ratio,
                             'psi_max': psi_max,
                             'psi_min': psi_min}

        # parameters: what is this epsilon for?
        # which names for setting psi max and psi min

    def set_param(self, key=None, value=None):
        """ Sets the parameters that we want, without
        adding new and unwanted options.
        User can pass in new values directly in,
        or can be prompted.
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

    def get_param(self, key=None):
        """ Returns a the value associated with a key. """
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

    def import_psi_data(self, plot=False):
        """ Read psi data using fortran code """
        import Efit.efit as efit
        from Interpol.Setup_Grid_Data import Efit_Data
        efit.readg()  # defines the variables
        # extract dimensions and boundaries
        nxefit, nyefit, nbdry, nlim = efit.get_nxy()

        fold, fpol, pres, qpsi, \
            rbdry, zbdry, xlim, ylim, xdim, zdim, \
            rcentr, rgrid1, zmid, rmagx, zmagx, \
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

        if plot:
            plt.figure('efit data')
            plt.subplot(1, 2, 1)
            plt.plot(rbdry, zbdry, label='plasma boundary')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.legend()

            plt.subplot(1, 2, 2)
            plt.contour(self.efit_psi.x, self.efit_psi.y, self.efit_psi.v, 30)
            plt.title('Original Cross-Section')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlim(rmin, rmax)
            plt.ylim(zmin, zmax)
            plt.show()

    def read_target_plate(self):
        """ Reads the coordinates for a line defining the inner
        and outer target plates.
        """
        self.itp = []
        with open('../data/SNL/itp1.txt') as f:
            for line in f:
                point = line.strip()
                if point.startswith('#'):
                    # move along to the next iteration
                    continue
                x = float(point.split(',')[0])
                y = float(point.split(',')[1])
                self.itp.append((x, y))

        self.otp = []
        with open('../data/SNL/otp1.txt') as f:
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

    def add_magx(self, r=None, z=None):
        """ Adds the magnetic axis using the coordinates of last point
        that the user has clicked.
        """
        if r is None and z is None:
            self.magx = self.root_finder.final_root
        else:
            self.magx = (r, z)
        print('Added {} as the magnetic axis.'.format(self.magx))

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
        self.new_line = LineTracing(self.psi_norm, self._grid_params,
                                    option='theta')

    def draw_radial_lines(self):
        """ trace the orthogonal lines to the stream function
        saves most recent line"""
        from line_tracing import LineTracing
        self.new_line = LineTracing(self.psi_norm, self._grid_params,
                                    option='rho', direction='ccw')

    def compute_eq_psi(self):
        """ Initializes the line tracing class for the construction
        of the grid. Necessary to be separate so we can wiat until the
        user clicks on the root.
        """
        from line_tracing import LineTracing
        self.eq = LineTracing(self.psi_norm, self._grid_params,
                              option='xpt_circ')
        self.eq.calc_equal_psi_points(self.xpt1[0], self.xpt1[1])
        self.eq.disconnect()

    def construct_SNL_patches(self, movie=False):
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
#        dp = {}  # defining points
        psi_max = self._grid_params['psi_max']
        psi_min = self._grid_params['psi_min']

        # IDL ===================================================
#        idl_e = self.eq.draw_line(xpt['W'], psi_max, option='rho',
#                                  direction='ccw', show_plot=False)
##        idl_e = Line([idl_e_temp[-1], idl_e_temp[0]])
##        idl_e.plot('red')
##        dp['ixpt'] = idl_e.p[0]
#        idl_n = self.eq.draw_line(dp['ixpt'], self.itp, option='theta',
#                                  direction='ccw', show_plot=False)
##        idl_nn = Line(idl_n)
##        idl_n.plot('red')
#
#
#        idl_s = self.eq.draw_line(xpt['SW'], self.itp, option='theta',
#                                  direction='ccw', show_plot=False)
##        idl_ss = Line(idl_s)
##        idl_s.plot('red')
#
#
#
#        idl_w = Line([idl_n.p[-1], idl_s.p[-1]])
##        idl_w.plot('red')
#        IDL = Patch(idl_n, idl_s[::-1])
#        IDL.fill()

        # +++ #### ========
        # trying a new way of defining
        E = self.eq.draw_line(xpt['W'], psi_max, option='rho', direction='ccw')
        N = self.eq.draw_line(E.p[-1], self.itp, option='theta', direction='ccw').reverse()
        S = self.eq.draw_line(xpt['SW'], self.itp, option='theta', direction='ccw')
        E = Line([N.p[-1], S.p[0]]) # straighten it up
        W = Line([S.p[-1], N.p[0]])
        IDL = Patch([N, E, S, W])
        if movie:
            IDL.plot_border()
            IDL.fill()
            plt.draw()
        # Lines are now saved inside of the patch


        # IPF ===================================================
#        ipf_n = idl_s
#        ipf_e_temp = self.eq.draw_line(xpt['S'], psi_min, option='rho',
#                                  direction='cw', show_plot=False)
#        ipf_e = Line([ipf_e_temp[0], ipf_e_temp[-1]])
#        ipf_e.plot('red')
#        ipf_s = self.eq.draw_line(ipf_e.p[-1], self.itp, option='theta',
#                                  direction='ccw')
#        ipf_w = Line([ipf_n[-1], ipf_s[-1]])
#        ipf_w.plot('red')
#        IPF = Patch(ipf_n, ipf_s[::-1])
#        IPF.fill()
        N = IDL.S.reverse()
        E = self.eq.draw_line(xpt['S'], psi_min, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(E.p[-1], self.itp, option='theta', direction='ccw')
        W = Line([S.p[-1], N.p[0]])
        IPF = Patch([N, E, S, W])
        if movie:
            IPF.plot_border()
            IPF.fill()
            plt.draw()
        
#        # OPF ===================================================
#        opf_w = ipf_e
#        opf_s = self.eq.draw_line(opf_w[-1], self.otp, option='theta',
#                                  direction='cw')
#        opf_n = self.eq.draw_line(xpt['SE'], self.otp, option='theta',
#                                  direction='cw')
#        opf_e = Line([opf_n[-1], opf_s[-1]])
#        opf_e.plot('red')
#        OPF = Patch(opf_s, opf_n[::-1])
#        OPF.fill()
        N = self.eq.draw_line(xpt['SE'], self.otp, option='theta', direction='cw')
        W = IPF.E.reverse()
        S = self.eq.draw_line(W.p[0], self.otp, option='theta', direction='cw').reverse()
        E = Line([N.p[-1], S.p[0]])
        OPF = Patch([N, E, S, W])
        if movie:
            OPF.plot_border('orange')
            OPF.fill()
            plt.draw()

#        # ODL ===================================================
#        odl_w = self.eq.draw_line(xpt['E'], psi_max, option='rho',
#                                  direction='ccw')
#        dp['oxpt'] = odl_w[-1]
#        odl_n = self.eq.draw_line(dp['oxpt'], self.otp, option='theta',
#                                  direction='cw')
#        odl_s = opf_n
#        odl_e = Line([odl_n[-1], odl_s[-1]])
#        odl_e.plot('red')
#        ODL = Patch(odl_n, odl_s[::-1])
#        ODL.fill()
        W = self.eq.draw_line(xpt['E'], psi_max, option='rho', direction='ccw')
        W = Line([W.p[0], W.p[-1]])
        N = self.eq.draw_line(W.p[-1], self.otp, option='theta', direction='cw')
        S = OPF.N.reverse()
        E = Line([N.p[-1], S.p[0]])
        ODL = Patch([N, E, S, W])
        if movie:
            ODL.plot_border('blue')
            ODL.fill()
            plt.draw()
        

        # need the mid and top points of the separatrix
        sep = self.eq.draw_line(xpt['NW'], xpt['NE'], option='theta', direction='cw')
#        sep_x = np.array([p.x for p in sep])
#        sep_y = np.array([p.y for p in sep])
        sep_x = sep.xval
        sep_y = sep.yval
        top_index = np.argmax(sep_y)
        midpoint = np.median(sep_y)
        mid_index, = np.where(np.abs(sep_y-midpoint) < 1e-3)
        dp = {}
        dp['top'] = Point(sep_x[top_index], sep_y[top_index])
        dp['impt'] = Point(sep_x[mid_index[0]], sep_y[mid_index[0]])
        dp['ompt'] = Point(sep_x[mid_index[1]], sep_y[mid_index[1]])

#        # ISB ===================================================
#        isb_e = self.eq.draw_line(dp['impt'], psi_max, option='rho',
#                                  direction='ccw')
#        isb_n = self.eq.draw_line(dp['ixpt'], isb_e[-1], option='theta',
#                                  direction='cw')
#        isb_s = sep[:mid_index[0]+1]
#        ISB = Patch(isb_n, isb_s[::-1])
#        ISB.fill()
        E = self.eq.draw_line(dp['impt'], psi_max, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(IDL.E.p[0], E.p[0], option='theta', direction='cw')
        S = Line(sep.p[:mid_index[0]+1]).reverse()
        W = IDL.E.reverse()
        ISB = Patch([N, E, S, W])
        

#        # ICB ===================================================
#        icb_e = self.eq.draw_line(dp['impt'], psi_min, option='rho',
#                                  direction='cw')
#        icb_w = self.eq.draw_line(xpt['N'], psi_min, option='rho',
#                                  direction='cw')
#        icb_s = self.eq.draw_line(icb_w[-1], icb_e[-1], option='theta',
#                                  direction='cw')
#        icb_n = isb_s
#        ICB = Patch(icb_s[::-1], icb_n)
#        ICB.fill()
        N = Line(sep.p[:mid_index[0]+1])
        E = self.eq.draw_line(dp['impt'], psi_min, option='rho', direction='cw')
        W = self.eq.draw_line(xpt['N'], psi_min, option='rho', direction='cw').reverse()
        W = Line([W.p[0], W.p[-1]])
        S = self.eq.draw_line(W.p[0], E.p[-1], option='theta', direction='cw').reverse()
        ICB = Patch([N, E, S, W])
        

#        # IST ===================================================
#        ist_e = self.eq.draw_line(dp['top'], psi_max, option='rho',
#                                  direction='ccw')
#        ist_n = self.eq.draw_line(isb_n[-1], ist_e[-1], option='theta',
#                                  direction='cw')
#        ist_s = sep[mid_index[0]:top_index+1]
#        IST = Patch(ist_n, ist_s[::-1])
#        IST.fill()
        E = self.eq.draw_line(dp['top'], psi_max, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(ISB.N.p[-1], E.p[0], option='theta', direction='cw')
        S = Line(sep.p[mid_index[0]:top_index+1]).reverse()
        W = ISB.W.reverse()
        IST = Patch([N, E, S, W])

#        # ICT ===================================================
#        ict_e = self.eq.draw_line(dp['top'], psi_min, option='rho',
#                                  direction='cw')
#        ict_s = self.eq.draw_line(icb_s[-1], ict_e[-1], option='theta',
#                                  direction='cw')
#        ict_n = ist_s
#        ICT = Patch(ict_n, ict_s[::-1])
#        ICT.fill()
        E = self.eq.draw_line(dp['top'], psi_min, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(ICB.S.p[0], E.p[-1], option='theta', direction='cw').reverse()
        N = IST.S.reverse()
        W = ICB.E.reverse()
        ICT = Patch([N, E, S, W])
        
#        # OST ===================================================
#        ost_e = self.eq.draw_line(dp['ompt'], psi_max, option='rho',
#                                  direction='ccw')
#        ost_n = self.eq.draw_line(ist_n[-1], ost_e[-1], option='theta',
#                                  direction='cw')
#        ost_s = sep[top_index: mid_index[1]+1]
#        OST = Patch(ost_n, ost_s[::-1])
#        OST.fill()
        E = self.eq.draw_line(dp['ompt'], psi_max, option='rho', direction='ccw').reverse()
        E = Line([E.p[0], E.p[-1]])
        N = self.eq.draw_line(IST.N.p[-1], E.p[0], option='theta', direction='cw')
        S = Line(sep.p[top_index: mid_index[1]+1]).reverse()
        W = IST.W.reverse()
        OST = Patch([N, E, S, W])

#        # OCT ===================================================
#        oct_e = self.eq.draw_line(dp['ompt'], psi_min, option='rho',
#                                  direction='cw')
#        oct_s = self.eq.draw_line(ict_e[-1], oct_e[-1], option='theta',
#                                  direction='cw')
#        oct_n = ost_s
#        OCT = Patch(oct_n, oct_s[::-1])
#        OCT.fill()
        E = self.eq.draw_line(dp['ompt'], psi_min, option='rho', direction='cw')
        E = Line([E.p[0], E.p[-1]])
        S = self.eq.draw_line(ICT.E.p[-1], E.p[-1], option='theta', direction='cw').reverse()
        N = Line(sep.p[top_index: mid_index[1]+1])
        W = ICT.E.reverse()
        OCT = Patch([N, E, S, W])

#        # OCB ===================================================
#        ocb_w = oct_e
#        ocb_s = self.eq.draw_line(ocb_w[-1], icb_w[-1], option='theta',
#                                  direction='cw')
#        ocb_n = sep[mid_index[1]:]
#        OCB = Patch(ocb_s, ocb_n[::-1])
#        OCB.fill()
        W = OCT.E.reverse()
        N = Line(sep.p[mid_index[1]:])
        S = self.eq.draw_line(W.p[0], ICB.W.p[0], option='theta', direction='cw').reverse()
        E = ICB.W.reverse()
        OCB = Patch([N, E, S, W])

#        # OSB ===================================================
#        osb_n = self.eq.draw_line(ost_e[-1], dp['oxpt'], option='theta',
#                                  direction='cw')
#        osb_s = ocb_n
#        OSB = Patch(osb_n, osb_s[::-1])
#        OSB.fill()
        W = OST.E.reverse()
        N = self.eq.draw_line(W.p[-1], ODL.W.p[-1], option='theta', direction='cw')
        S = OCB.N.reverse()
        E = ODL.W.reverse()
        OSB = Patch([N, E, S, W])

        self.patches = [IDL, IPF, OPF, ODL, ISB, ICB,
                        IST, ICT, OST, OCT, OCB, OSB]
        if not movie:
            for patch in self.patches:
                patch.plot_border()
                patch.fill()


    def patch_diagram(self):
        colors = ['salmon', 'skyblue', 'mediumpurple', 'mediumaquamarine',
                  'sienna', 'orchid', 'lightblue', 'gold', 'steelblue',
                  'seagreen', 'firebrick', 'saddlebrown']

        plt.figure('patches')
        for i in range(len(self.patches)):
            self.patches[i].plot_borders('red')
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
        pass

    def test_interpol(self, option=2, nfine=100, ncrude=10, tag='v'):
        """ Provides a demonstration and test of the bicubic interpolation
        methods used. Wrapper for the real code from another module.
        """
        from Interpol.Test_Interpol import test_interpol
        test_interpol(option, nfine, ncrude, tag)


if __name__ == '__main__':
    plt.close('all')
    grid = Ingrid()
    grid.import_psi_data()
    grid.calc_efit_derivs()
    grid.plot_efit_data()
    grid.find_roots()

    # for the snl neqdsk
    grid.read_target_plate()
    grid.add_magx(0.6818276108184476, -0.0036146834671545633)
    grid.add_xpt1(0.563043117220232, -0.39498650311360006)
    grid.calc_psinorm()
    plt.close('Efit Data')
    grid.plot_target_plate()
    grid.compute_eq_psi()
