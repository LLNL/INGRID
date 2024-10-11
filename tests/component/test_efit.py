from pathlib import Path
from typing import List, Callable
import numpy as np
from freeqdsk import geqdsk
from ingrid.interface import Ingrid
from ingrid.utils.helper import resolve_path

def test_load_geqdsk(data_dir: Path) -> None:
    """
    Test loading of eqdsk files via the Ingrid class
 
    Parameters
    ----------
    data_dir : pathlib.Path
        Pytest fixture for the data dir
    """
    eqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True)

    #
    # Load data using FreeQDSK
    #
    with open(file=eqdsk_path, mode='r') as f:
        baseline: geqdsk.GEQDSKFile = geqdsk.read(f)

    #
    # Load data using Ingrid class
    #
    session = Ingrid()
    session.LoadGEQDSK(geqdsk_path=eqdsk_path)

    #
    # Check keys and values of loaded SNL eqdsk are consistent.
    # The geqdsk dataclass has fields that need to be accessed via __dataclass_fields__
    #
    baseline_keys: List[str] = sorted(baseline.__dataclass_fields__.keys())
    test_keys: List[str] = sorted(session.geqdsk_data.__dataclass_fields__.keys())
    explicit_comparisons: List[str] = ['comment']  # Keys to use "==" operator with
    for baseline_k, test_k in zip(baseline_keys, test_keys):
        assert baseline_k == test_k, "Key mismatch in loaded data."
        k: str = baseline_k
        #
        # Use a lambda to select the appropriate comparison operator
        #
        compare: Callable = lambda x, y: x == y if k in explicit_comparisons else np.allclose
        failure_message: str = "Numerics mismatch in data."
        assert compare(getattr(baseline, k), getattr(session.geqdsk_data, k)), failure_message
