from pathlib import Path
from INGRID.ingrid import Ingrid


gridue_dict = Ingrid.ImportGridue(str(Path(__file__).parent / "gridue"))
Ingrid.PlotGridue(gridue_dict)