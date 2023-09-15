import numpy as np
import pytest
import utils
from INGRID.ingrid import Ingrid
from freeqdsk import geqdsk

def test_load_geqdsk(data_dir):
    """
    Test loading of eqdsk files via the Ingrid class
 
    Parameters
    ----------
    data_dir : pathlib.Path
        Pytest fixture for the data dir
    """
    eqdsk_path = data_dir / 'SNL' / 'DIII-D' / 'neqdsk'
    eqdsk_path = utils.resolve_path(eqdsk_path, as_str=True)

    #
    # Load data using FreeQDSK
    #
    with open(eqdsk_path, 'r') as f:
        baseline = geqdsk.read(f)

    #
    # Load data using Ingrid class
    #
    session = Ingrid()
    session.LoadGEQDSK(geqdsk_path=eqdsk_path)

    #
    # Check keys and values of loaded SNL eqdsk are consistent
    #
    for baseline_k, test_k in zip(sorted(baseline), sorted(session.geqdsk_data)):
        assert baseline_k == test_k, "Key mismatch in loaded data."
        k = baseline_k
        assert np.allclose(baseline[k], session.geqdsk_data[k]), "Numerics mismatch in data."

