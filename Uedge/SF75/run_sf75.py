from uedge import *
from uedge.hdf5 import *
import uedge_mvu.plot as mp
import uedge_mvu.utils as mu


bbb.gengrid=0 #-don't generate grid
bbb.mhdgeo=1 #-toroidal MHD equilibrium
com.geometry="snowflake75"
com.nxpt=2 #-how many X-points in the domain

#-non-orthogonal grid settings
isnonog=1
methg=66

bbb.gallot("Xpoint_indices",0)
grd.readgrid("gridue",com.runid)

com.nx=com.nxm
com.ny=com.nym

mp.plotmesh()
mu.paws()


#-set flat initial profiles
bbb.allocate()
bbb.tes=10*bbb.ev
bbb.tis=10*bbb.ev
bbb.nis=1e20
bbb.ngs=1e16
bbb.ups=0

#-set up physics model
bbb.isteon=1
bbb.istion=1
bbb.isnion=1
bbb.isupon=1
bbb.isupgon=0
bbb.isngon=0

#-do a short initial run first, to check if everything is set up ok
bbb.restart=1
bbb.isbcwdt=1
bbb.dtreal=1e-12
bbb.ftol=1e-6
bbb.exmain()

#-run to steady state
#bbb.rundt()

#-or restore a solution
#hdf5_restore("tmp.h5")
#bbb.exmain()
