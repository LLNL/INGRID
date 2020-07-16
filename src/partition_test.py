from Ingrid import Ingrid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
from geometry import Point, Line, Patch, segment_intersect
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from matplotlib.patches import Polygon
import pathlib

def SF15_test():
    SF15_case = "../Parameter Files/SF15/SF15.yml"
    fname = SF15_case
    SF15 = Ingrid(InputFile=fname)
    try:
        SF15.Setup()
    except:
        return SF15.eq.config
def SF45_test():
    SF45_case = "../Parameter Files/SF45/SF45.yml"
    fname = SF45_case
    SF45 = Ingrid(InputFile=fname)
    try:
        SF45.Setup()
    except:
        return SF45.eq.config

def SF75_test():
    SF75_case = "../Parameter Files/SF75/SF75.yml"
    fname = SF75_case
    SF75 = Ingrid(InputFile=fname)
    try:
        SF75.Setup()
    except:
        return SF75.eq.config

def SF105_test():
    SF105_case = "../Parameter Files/SF105/SF105.yml"
    fname = SF105_case
    SF105 = Ingrid(InputFile=fname)
    try:
        SF105.Setup()
    except:
        return SF105.eq.config

def ADX_test():
    ADX_case = "../Parameter Files/SF45/ADX.yml"
    fname = ADX_case
    ADX = Ingrid(InputFile=fname)
    try:
        ADX.Setup()
    except:
        return ADX.eq.config

def SPARC_test():
    SPARC_case = "../Parameter Files/SF45/SPARC_SF.yml"
    fname = SPARC_case
    SPARC = Ingrid(InputFile=fname)
    try:
        SPARC.Setup()
    except:
        return SPARC.eq.config

sf15=SF15_test()
sf45=SF45_test()
sf75=SF75_test()
sf105=SF105_test()
adx=ADX_test()
sparc=SPARC_test()

print('\n-----------------------------------------------------')
print(f"# SF15_test() classified '{sf15}' configuration...")
print(f"# SF45_test() classified '{sf45}' configuration...")
print(f"# SF75_test() classified '{sf75}' configuration...")
print(f"# SF105_test() classified '{sf105}' configuration...")
print(f"# ADX_test() classified '{adx}' configuration...")
print(f"# SPARC_test() classified '{sparc}' configuration...")
print('-----------------------------------------------------\n')




