from Ingrid import Ingrid
import matplotlib.pyplot as plt

fpath = '../data/SF95/neqdsk'
nlevs = 200

EfitPlot = Ingrid(EqFile=fpath)
EfitPlot.yaml['grid_params']['nlevs'] = nlevs
EfitPlot.OMFIT_read_psi()
EfitPlot.plot_efit_data()
plt.ioff()
plt.show()