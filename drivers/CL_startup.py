"""
CL_startup

    Description:
        A script demonstrating basic command line workflow with Ingrid when utilizing the YAML
        parameter file as a driver.

"""

import sys
sys.path.append('../ingrid/')
from ingrid import Ingrid
import pathlib
import numpy as np

# Path to parameter file cases
LSN_case = "../Parameter Files/SNL/LSN_YAML_EXAMPLE.yml"
CMOD_case = "../Parameter Files/SNL/cmod_param.yml"
DIIID_case = "../Parameter Files/SNL/DIIID_SNL.yml"
SAS_case = "../Parameter Files/SNL/SAS1_modif.yml"
SPARC_case = "../Parameter Files/DNL/SPARC_DNL.yml"

SF45_maxim = '../Parameter Files/NEQDSK_maxim/neqdsk_6.yml'
# < Your path to Ingrid formatted parameter file here >

# Select a case.
case = DIIID_case


if __name__ == '__main__':

    # Construction of Ingrid object with parameter file parsing.
    GridDemo = Ingrid(InputFile=case)
    

    # Read EFIT data, target plate data, and plot data.
    GridDemo.OMFIT_read_psi()
    GridDemo.calc_efit_derivs()
    GridDemo.AutoRefineMagAxis()
    GridDemo.AutoRefineXPoint()

    GeoSettings = {
    
        'W1' : {
            'file' : GridDemo.settings['target_plates']['plate_W1']['file'], 
            'zshift' : GridDemo.settings['target_plates']['plate_W1']['zshift']
            },

        'E1' : {
            'file' : GridDemo.settings['target_plates']['plate_E1']['file'],
            'zshift' : GridDemo.settings['target_plates']['plate_E1']['zshift']
        }
    }

    GridDemo.SetGeometry(GeoSettings)
    GridDemo.SetLimiter(GridDemo.settings['limiter']['file'])
    GridDemo.SetMagReference()
    GridDemo.CalcPsiNorm()
    GridDemo.AnalyzeTopology()

    # Begin patch construction with parameter file psi values.
    GridDemo.ConstructPatches()
    GridDemo.ShowSetup()

    # Begin patch refinement and actively plot grid.

    #CellCorrection={'ThetaMin':60,'ThetaMax':120,'Resolution':1000,'Active':True}
    #GridDemo.CurrentTopology.distortion_correction={'all' : CellCorrection}
    GridDemo.CreateSubgrid(NewFig=True)
    # Export gridue file.
    fname = 'gridue'
    GridDemo.ExportGridue(fname=fname)
