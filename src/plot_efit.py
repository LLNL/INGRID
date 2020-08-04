from Ingrid import Ingrid
import matplotlib.pyplot as plt

fpath = '../data/SNL/DIII-D/neqdsk'
nlevs = 35

EfitPlot = Ingrid(EqFile=fpath)
EfitPlot.settings['grid_params']['nlevs'] = nlevs
EfitPlot.OMFIT_read_psi()
EfitPlot.plot_efit_data()
plt.ioff()
plt.show()