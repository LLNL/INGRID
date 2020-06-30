from Ingrid import Ingrid
import matplotlib.pyplot as plt

fpath = '../data/SF45/neqdsk'#../data/SF95/neqdsk'
nlevs = 50

EfitPlot = Ingrid(EqFile=fpath)
EfitPlot.yaml['grid_params']['nlevs'] = nlevs
EfitPlot.OMFIT_read_psi()
EfitPlot.plot_efit_data()
import pdb
pdb.set_trace()