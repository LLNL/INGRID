"""
GUI_startup

    Description:
        Starts Ingrid in GUI mode.
"""
import sys
sys.path.append('../src/')
from Ingrid import Ingrid
import pathlib

if __name__ == '__main__':
    GridDemo = Ingrid()
    GridDemo.StartGUI()