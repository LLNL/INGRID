"""
CL_startup

    Description:
        A script demonstrating basic command line workflow with Ingrid when utilizing the YAML
        parameter file as a driver.

"""

import sys
sys.path.append('../src/')
from Ingrid import Ingrid
import pathlib

# Path to parameter file cases
LSN_case = "../Parameter Files/SNL/LSN_YAML_EXAMPLE.yml"
CMOD_case = "../Parameter Files/SNL/cmod_param.yml"
DIIID_case = "../Parameter Files/DIIID_SNL.yml"
SAS_case = "../Parameter Files/SNL/SAS1_modif.yml"
SPARC_case = "../Parameter Files/DNL/SPARC_DNL.yml"
# < Your path to Ingrid formatted parameter file here >

# Select a case.
case = SPARC_case


if __name__ == '__main__':

    # Construction of Ingrid object with parameter file parsing.
    GridDemo = Ingrid(InputFile=case)

    # Read EFIT data, target plate data, and plot data.
    GridDemo.Setup()

    # Begin patch construction with parameter file psi values.
    GridDemo.ConstructPatches()

    # Plot constructed patches.
    GridDemo.ShowSetup()

    # Begin patch refinement and actively plot grid.
    GridDemo.CreateSubgrid()

    # Export gridue file.
    fname = 'gridue'
    GridDemo.export(fname=fname)