from Ingrid import Ingrid
import pathlib

GridDemo = Ingrid(InputFile="../Parameter Files/USN_YAML_EXAMPLE.yml")
GridDemo.Setup()
GridDemo.ConstructPatches()
GridDemo.CreateSubgrid(NewFig=False)
GridDemo.export()