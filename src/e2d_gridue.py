from Ingrid import Ingrid,LSN
from Interpol.Setup_Grid_Data import *
import e2dgrid as e
from  matplotlib.pyplot import figure,ion
import numpy as np
from copy import copy
ion()

i=Ingrid()
i.yaml['eqdsk']='/Users/holma2/Dropbox (Aalto)/UEDGE/INGRID/INGRID/seq#1/eqdsk'
i.OMFIT_read_psi()
i.psi_norm=0
i.eq=0
i.plate_data=0
l=LSN(i)
l.set_gridue_manual(*e.grid('/Users/holma2/Dropbox (Aalto)/UEDGE/E2D-EIR_integration/exporte2dgrid/aholm_mar1719_11'))
i.write_gridue(l.gridue_params,fname = 'gridue_e2d')

''' Test interpolation '''
rmin=i.OMFIT_psi['RLEFT']
rmax=rmin+i.OMFIT_psi['RDIM']
zmin=i.OMFIT_psi['ZMID']-0.5*i.OMFIT_psi['ZDIM']
zmax=i.OMFIT_psi['ZDIM']+zmin
nr=i.OMFIT_psi['NW']
nz=i.OMFIT_psi['NH']

tag='v'
grid0=i.efit_psi
grid0.Calculate_PDeriv(unit_spacing=False)
nlev = 30
lev = (grid0.get_v(tag).min() + (grid0.get_v(tag).max()
        - grid0.get_v(tag).min()) * np.arange(nlev) / (nlev-1))
xpt=6.05e-2
lev[np.argmin(abs(lev-xpt))]=xpt # Go close to sep


for n in [50,64,65,66,100]:
    f=figure("Interpol Demo EQDSK-{} res={}x{}.".format(tag,n,n),figsize=(6,8))
    ax=f.add_subplot(111)
    CS=ax.contour(grid0.r, grid0.z, grid0.get_v(tag), lev, colors='black',linewidths=1.5,linestyles='solid')


    ''' Create finer grid for interpolation '''
    grid1=Efit_Data(nr=n,nz=n,rmin=rmin,rmax=rmax,zmin=zmin,zmax=zmax)
    for i in range(grid1.nr):
        for j in range(grid1.nz):
            grid1.set_v(grid0.get_psi(grid1.r[i,j],grid1.z[i,j],tag=tag),(i,j),tag)


    ax.contour(grid1.r, grid1.z, grid1.get_v(tag), lev, colors='r',linewidths=0.5,linestyles='solid')

    ax.set_title('n_orig=65x65 vs n_interpol={}x{}'.format(n,n))
    f.show()
    f.savefig('interpolcomp_{}x{}.png'.format(n,n))

