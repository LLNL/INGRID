from Ingrid import Ingrid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from geometry import Point, Line, Patch, segment_intersect
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from matplotlib.patches import Polygon
import pathlib

# def SF15_test():
#     SF15_case = "../Parameter Files/SF15/SF15.yml"
#     fname = SF15_case
#     SF15 = Ingrid(InputFile=fname)
#     try:
#         import pdb
#         pdb.set_trace()
#         SF15.Setup(topology='DNL', use_efit_bounds=True)
#     except:
#         pass
# def SF45_test():
#     SF45_case = "../Parameter Files/SF45/SF45.yml"
#     fname = SF45_case
#     SF45 = Ingrid(InputFile=fname)
#     coor=[(1.0, -1), (2.4, -1), (2.4, 1.0), (1.0, 1.0), (1.0, -1)]
#     try:
#         SF45.Setup(topology='DNL', limiter_coordinates=coor)
#     except:
#         pass

def SF75_test():
    SF75_case = "../Parameter Files/SF75/SF75.yml"
    fname = SF75_case
    SF75 = Ingrid(InputFile=fname)
    try:
        SF75.Setup(topology='DNL', use_efit_bounds=True)
    except:
        return SF75.eq.config

def SF105_test():
    SF105_case = "../Parameter Files/SF105/SF105.yml"
    fname = SF105_case
    SF105 = Ingrid(InputFile=fname)
    try:
        SF105.Setup(topology='DNL', use_efit_bounds=True)
    except:
        return SF105.eq.config

def ADX_test():
    ADX_case = "../Parameter Files/SF45/ADX.yml"
    fname = ADX_case
    ADX = Ingrid(InputFile=fname)
    try:
        ADX.Setup(topology='DNL', user_efit_bounds=True)
    except:
        return ADX.eq.config

def SPARC_test():
    SPARC_case = "../Parameter Files/SF45/SPARC_SF.yml"
    fname = SPARC_case
    SPARC = Ingrid(InputFile=fname)
    try:
        SPARC.Setup(topology='DNL')
    except:
        return SPARC.eq.config

sf75=SF75_test()
sf105=SF105_test()
adx=ADX_test()
sparc=SPARC_test()

print('\n-----------------------------------------------------')
print(f"# SF75_test() classified '{sf75}' configuration...")
print(f"# SF105_test() classified '{sf105}' configuration...")
print(f"# ADX_test() classified '{adx}' configuration...")
print(f"# SPARC_test() classified '{sparc}' configuration...")
print('-----------------------------------------------------\n')




