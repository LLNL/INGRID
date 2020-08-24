#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:42 2019

@author: Joey Watkins
"""

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

import GUI.SimpleGUI as sg

class RootFinder:
    """
    Finds the root closest to a point the user clicks on. 
    If active is set to false the find_roots method must be passed
    an initial guess.
    
    Parameters
    ----------
    grid : Setup_Grid_Data.Efit_Data
        Uses the grid definitions defined in the Efit_Data class
        to find the roots.
    
    """
    def __init__(self, grid, active=True, mode = 'root_finder', controller=None):
        self.grid = grid
        self.root_finding = True if mode == 'root_finder' else False
        self.psi_finding = True if mode == 'psi_finder' else False
        self.controller = controller

        self.curr_x = 0.0
        self.curr_y = 0.0
        self.final_root = ()
        self.root_trace={'default' : {'point' : None, 'separatrix' : None},
            'xpt1' : {'point' : None, 'separatrix' : None},
            'xpt2' : {'point' : None, 'separatrix' : None},
            'magx' : {'point' : None, 'separatrix' : None}}

        if active:
            # self.cid = grid.ax.figure.canvas.mpl_connect('button_press_event', self)
            # print("Entering Root Finder. "
            #       + "Right click to disable.")
            pass
        else:
            print("Root Finder on standby.")

    def func(self, xy):
        """ Combines the first partial derivatives to solve the system for
        maximum, minimum, and saddle locations.     
        
        Parameters
        ----------
        xy : array-like
            Contains x and y. Ex: xy = (x0, y0).
            
        Returns
        -------
        F : array
            Vector function to be used in find root.
        """
        # combine the deriv functions to solve the system
        x, y = xy
        F = np.zeros(2)
        F[0] = self.grid.get_psi(x, y, tag='vr')
        F[1] = self.grid.get_psi(x, y, tag='vz')
        return F

    def __call__(self, event):
        """ Executes when the user clicks on the plot. """
        if event.button == 3:
            print('Root Finder Disabled.')
            event.canvas.mpl_disconnect(self.cid)
            return
        if event.inaxes != self.grid.ax.axes:
            # safe gaurd against clicking outside the figure
            return
        
        self.curr_x, self.curr_y = event.xdata, event.ydata
        if self.root_finding: 
            self.find_root(self.curr_x, self.curr_y)
        elif self.psi_finding:
            self.find_psi(self.curr_x, self.curr_y, self.topofeature)
        else:
            print("You chose ({0:.5f}, {1:.5f}). ".format(self.curr_x, self.curr_y))
            self.final_root = (self.curr_x, self.curr_y)
    
    def toggle_root_finding(self):
        """ Activates or deactivates the root finding capacity.
        Leaves the click ability active so the user can save clicks.
        """
        if self.root_finding:
            self.root_finding = False
        else:
            self.root_finding = True
        
    def disconnect(self):
        """ Turn off the click functionality for the root finder. """
        self.grid.ax.figure.canvas.mpl_disconnect(self.cid)
        print('RF Disabled')
        self.topofeature=None
   
    def connect(self, topofeature=None):
        """ Turn on the click functionality for the root finder. """
        self.topofeature = topofeature
        self.cid = self.grid.ax.figure.canvas.mpl_connect('button_press_event', self)
        print('RF Enabled')
    
    def find_root(self, x, y, topofeature='default'):
        """        
        Accepts an initial guess for the root, then uses scipy.optimize.root
        to refine the zero point. Saves the root as a tuple
        called self.final_root.
        
        Parameters
        ---------
        x : float
            x or r value for the guess
        y : float
            y or z value for the guess

        """
        # plt.cla()
        # self.controller.IngridSession.plot_efit_data()
        # self.controller.IngridSession.plot_target_plates()
        #plt.plot(x, y, '*', color = 'blue')
        sol = root(self.func, [x, y])
        r, z = sol.x[0], sol.x[1]

        if (not self.grid.rmin < r < self.grid.rmax or not
                self.grid.zmin < z < self.grid.zmax):
            print("You clicked too far from the true zero point")
        else:
            print("You chose ({0:.5f}, {1:.5f}). ".format(x, y) +
                  "The zero point is ({0:.5f},{1:.5f})".format(r, z))
        
        try:
            self.root_trace[topofeature]['separatrix'].collections[0].remove()
            self.grid.ax.lines[-1].remove()
        except:
            pass
        self.root_trace[topofeature]['point'] = plt.plot(r, z, '.', color = 'blue')
        self.root_trace[topofeature]['separatrix'] = plt.contour(self.grid.r, self.grid.z, self.grid.v, [self.grid.get_psi(r,z)], colors = 'purple')
        plt.draw()
        self.final_root = (r, z)

        self.controller.frames[sg.ParamPicker].curr_click = self.final_root
        # self.controller.curr_root = self.final_root
        self.controller.frames[sg.ParamPicker].update_root_finder()

    def find_psi(self, x, y, topofeature=None):
        self.controller.frames[sg.ParamPicker].curr_click = (x, y)
        self.controller.frames[sg.ParamPicker].update_root_finder()

        Dic={'psi_max':'blue',
            'psi_core':'magenta',
            'psi_max_east':'blue',
            'psi_max_west':'blue',
            'psi_pf_1':'green',
            'psi_pf_2':'yellow'}

        if topofeature in [k for k in Dic.keys()]:
            self.grid.PlotLevel(self.grid.get_psi(x, y), color=Dic[topofeature], label=topofeature)
        else:
            plt.contour(self.grid.r, self.grid.z, self.grid.v, [self.grid.get_psi(x,y)], colors = 'red', label = 'psi_line')

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
