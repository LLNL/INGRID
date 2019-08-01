###-in Python2 run it by doing execfile("doall.py")-###
#######################################################

from Ingrid import Ingrid
import matplotlib.pyplot as plt

def paws():
    programPause = raw_input("Press the <ENTER> key to continue...")




grid = Ingrid()

#-read data
grid.import_psi_data()
grid.read_target_plate()

#initial processing of the input data
grid.calc_efit_derivs()
grid.plot_efit_data()


#-now find the two roots: null-point and X-point
grid.find_roots()

#-first find the null-point (magnetic axis)
print("Click on the null-point") ##-need blocking here!
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
print("Click on the X-point")
paws()

print("Click on the null-point") ##-need blocking here!
paws()


#-construct the patch-map for SNL
grid.construct_SNL_patches()
