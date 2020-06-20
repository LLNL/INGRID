from Ingrid import Ingrid
import pathlib

USN_case = "../Parameter Files/USN_YAML_EXAMPLE.yml"
DIIID_case = "../Parameter Files/DIIID_SNL.yml"
SAS_case = "/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS1_modif.yml"

fname = DIIID_case

GridDemo = Ingrid(InputFile=fname)
GridDemo.Setup()
GridDemo.ConstructPatches()
GridDemo.ShowSetup()
GridDemo.CreateSubgrid()
GridDemo.export()