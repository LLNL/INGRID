from Ingrid import Ingrid
import matplotlib.pyplot as plt

fpath = '/Users/torvaltz/Desktop/SPARC_XPTD/V2_FREEGS_geqdsk_LSNX'
nlevs = 75

EfitPlot = Ingrid(EqFile=fpath)
EfitPlot.yaml['grid_params']['nlevs'] = nlevs
EfitPlot.OMFIT_read_psi()
EfitPlot.plot_efit_data()
plt.ioff()
plt.show()