import sys
sys.path.append('../src/')
from Ingrid import Ingrid
import pathlib

SPARC_case = "../Parameter Files/SF105/SF105.yml"
# < Your path to Ingrid formatted parameter file here >

# Select a case.
case = SPARC_case

GridDemo = Ingrid(InputFile=case)

# Read EFIT data, target plate data, and plot data.s
GridDemo.Setup(topology='DNL', use_efit_bounds=True)

# Begin patch construction with parameter file psi values.
GridDemo.ConstructPatches()

# Plot constructed patches.
GridDemo.ShowSetup()

# Begin patch refinement and actively plot grid.
GridDemo.CreateSubgrid()

# Export gridue file.
fname = 'gridue'
GridDemo.export(fname=fname)