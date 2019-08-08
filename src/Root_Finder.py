#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:42 2019

@author: watkins35

    Goal of this module: we want to create
    a plot of the psi function,
    then click on the plot at some point and
    calculate where the nearest zero is.

"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root


class RootFinder:
    """ Finds the root closest to a point the user clicks on. 
    Saves the value as self.final_root
    If active is set to false the find_roots method must be passed
    an initial guess.
    """
    def __init__(self, grid, active=True):
        self.grid = grid
        self.root_finding = True
        if active:
            self.cid = grid.ax.figure.canvas.mpl_connect('button_press_event', self)
            print("Entering Root Finder. "
                  + "Right click to disable.")
        else:
            print("Root Finder on standby.")

    def func(self, xy):
        # combine the deriv functions to solve the system
        x, y = xy
        F = np.zeros(2)
        F[0] = self.grid.get_psi(x, y, tag='vr')
        F[1] = self.grid.get_psi(x, y, tag='vz')
        return F

    def __call__(self, event):
        """ Executes when the user clicks on the plot """
        if event.button == 3:
            print('Root Finder Disabled.')
            event.canvas.mpl_disconnect(self.cid)
            return
        if event.inaxes != self.grid.ax.axes:
            # safe gaurd against clicking outside the figure
            return
        
        x, y = event.xdata, event.ydata
        if self.root_finding: 
            self.find_root(x, y)
        else:
            print("You chose ({0:.5f}, {1:.5f}). ".format(x, y))
            self.final_root = (x, y)
    
    def toggle_root_finding(self):
        """ Activates or deactivates the root finding capacity.
        Leaves the click ability active so the user can save clicks."""
        if self.root_finding:
            self.root_finding = False
        else:
            self.root_finding = True
        
    def disconnect(self):
        """ turn of the click functionality fo the root finder """
        self.grid.ax.figure.canvas.mpl_disconnect(self.cid)
    
    def find_root(self, x, y):
        plt.plot(x, y, 'x')
        plt.draw()
        sol = root(self.func, [x, y])
        r, z = sol.x[0], sol.x[1]

        if (not self.grid.rmin < r < self.grid.rmax or not
                self.grid.zmin < z < self.grid.zmax):
            print("You clicked too far from the true zero point")
        else:
            print("You chose ({0:.5f}, {1:.5f}). ".format(x, y) +
                  "The zero point is ({0:.5f},{1:.5f})".format(r, z))
            plt.plot(r, z, '1')  # the '1' determines the shape of the marker
            plt.draw()
            self.final_root = (r, z)   


if __name__ == "__main__":
#    from Interpol.Setup_Grid_Data import Efit_Data
    # need a driver to test the findRoots and newtons method
    # this is the global scope. Things defined here are defined everywhere
    # the code here demonstrates how to call the runNewton method
    plt.close('all')  # keep the screen clean
    data = True

    if data:
        from Read_Psi_Data import read_psi_data
        g = read_psi_data()
    else:
        from Interpol.Test_Functions import get_f
        from Interpol.Setup_Grid_Data import Efit_Data
        
        # set up the crude grid
        g = Efit_Data(-3, 3, 50, -3, 3, 50)
        # gaussian function
        g.set_v(get_f(g['r'], g['z'], option=5, x0=1, y0=-1))
    g.Calculate_PDeriv()  # this step must be done

    g.plot_data()

#    run_root_finder(fig, g)
    root_finder = RootFinder(g, active=False)
