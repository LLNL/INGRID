"""
CL_startup

    Description:
        A script demonstrating basic command line workflow with Ingrid when utilizing the YAML
        parameter file as a driver.

"""
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import sys
import pathlib
import numpy as np
sys.path.append('../ingrid/')

from ingrid import Ingrid  # noqa: E402


IG = Ingrid()

IG.LoadEFIT('../data/SNL/DIII-D/neqdsk')

geo_W = '../data/SNL/DIII-D/d3d_itp.txt'
geo_e = '../data/SNL/DIII-D/d3d_otp.txt'

IG.PlotPsiUNorm()
import pdb
pdb.set_trace()
IG.SetGeometry({'limiter': 'default'}, zshift=-1)
IG.SetGeometry({'W1': {'file': geo_W, 'zshift': 1.6}, 'E1': {'file': geo_e, 'zshift': 1.6}})
IG.PlotTargetPlates()
IG.SetGeometry({'E1': {'x': np.linspace(1.8, 2.4), 'y': np.linspace(0.5, 1.0)}})
IG.PlotTargetPlates()
IG.PlotLimiter()
IG.SetGeometry({'W1': {'file': geo_W, 'zshift': 1}, 'E1': {'file': geo_e, 'zshift': 1}})
IG.PlotLimiter()
IG.PlotTargetPlates()
IG.SetGeometry({'limiter': 'eq'}, zshift=1)
IG.PlotLimiter()
IG.SetGeometry({'limiter': 'eq'})
IG.SetGeometry({'W1': {'file': geo_W, 'zshift': 1.6}, 'E1': {'file': geo_e, 'zshift': 1.6}})
IG.PlotStrikeGeometry()

IG = Ingrid(InputFile='../Parameter Files/SNL/DIIID_SNL.yml')
IG.StartSetup()
IG.ShowSetup()