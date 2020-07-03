from Ingrid import Ingrid
import pathlib

USN_case = "../Parameter Files/USN_YAML_EXAMPLE.yml"
DIIID_case = "../Parameter Files/DIIID_SNL.yml"
SAS_case = "/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS1_modif.yml"
SPARC_case = "../Parameter Files/SPARC_DNL.yml"

fname = SAS_case

GridDemo = Ingrid(InputFile=fname)
GridDemo.yaml['target_plates']['plate_W1']['file'] = "/Users/torvaltz/Desktop/JeromeGridGenerator/GridGenerator/D3DSAS/SAS_odt.txt"
GridDemo.yaml['target_plates']['plate_E1']['file'] = "/Users/torvaltz/Desktop/INGRID/data/SNL/USN/itp4.txt"
import pdb
pdb.set_trace()
GridDemo.Setup()
GridDemo.ConstructPatches()
GridDemo.ShowSetup()
import pdb
pdb.set_trace()
GridDemo.CreateSubgrid(Enforce=False)
GridDemo.export()