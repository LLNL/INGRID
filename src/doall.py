###-in Python2 run it by doing execfile("doall.py")-###
#######################################################

from Ingrid import Ingrid
import matplotlib.pyplot as plt
import sys

if __name__ == '__main__':
        
    #raw_input = input # for python 3
    
    def paws():
        programPause = raw_input("Press the <ENTER> key to continue...")
    
    
    if sys.argv[1:]:
        # parse the stuff
        # filename, magx, xpt
        # g129883.05000 1.6,0 1.6,-.5
        gfile = sys.argv[1]
        interactive = False
        magx = tuple(map(float,sys.argv[2].split(',')))
        xpt = tuple(map(float,sys.argv[3].split(',')))
    
        grid = Ingrid(gfile)
        grid.add_magx(*magx)
        grid.add_xpt1(*xpt)
    
    else:
        interactive = True
        gfile = raw_input("Enter filename: ")
        #gfile = 'neqdsk'
    
        grid = Ingrid(gfile)
    
    #-read data
    #grid.import_psi_data()
    
    grid.OMFIT_read_psi()
    grid.read_target_plate()
    
    #initial processing of the input data
    grid.calc_efit_derivs()
    grid.plot_efit_data()
    
    
    #-now find the two roots: null-point and X-point
    grid.find_roots()
    
    if interactive:
        #-first find the null-point (magnetic axis)
        print("Click on the magnetic axis") ##-need blocking here!
        paws()
        grid.add_magx()
    
        #-next, find the primary X-point
        print("Click on the X-point") ##-need blocking here!
        paws()
        grid.add_xpt1()
    
    plt.close('Efit Data') #-finish with the raw Psi data
    
    
    
    #-calculate the normalized flux and plot normalized Psi data
    grid.calc_psinorm()
    
    # for the snl neqdsk
    grid.plot_target_plate()
    grid.compute_eq_psi()
    
    #-here we need to block until clicked
    #print("Click on the X-point")
    #paws()
    
    #-construct the patch-map for SNL
    from time import time
    start = time()
    grid.construct_SNL_patches()
    end = time()
    grid.patch_diagram()
    
    print("Time for grid: {} seconds.".format(end-start))
    
