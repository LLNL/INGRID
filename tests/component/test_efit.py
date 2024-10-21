from pathlib import Path
from typing import List, Callable
import numpy as np
from freeqdsk import geqdsk
from ingrid.utils.helper import resolve_path
from ingrid.core.interpol import GEQDSKInterpolator

def test_loading_geqdsk(data_dir: Path) -> None:
    """
    Test loading of eqdsk files via the GEQDSKInterpolator class
 
    Parameters
    ----------
    data_dir : pathlib.Path
        Pytest fixture for the data dir
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True)

    #
    # Load data using FreeQDSK
    #
    with open(file=geqdsk_path, mode='r') as f:
        baseline: geqdsk.GEQDSKFile = geqdsk.read(f)

    #
    # Load data using GEQDSKInterpolator class
    #
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)

    #
    # Check keys and values of loaded SNL eqdsk are consistent.
    # The geqdsk dataclass has fields that need to be accessed via __dataclass_fields__
    #
    baseline_keys: List[str] = sorted(baseline.__dataclass_fields__.keys())
    test_keys: List[str] = sorted(geqdsk_interpolator.geqdsk_data.__dataclass_fields__.keys())
    explicit_comparisons: List[str] = ['comment']  # Keys to use "==" operator with
    for baseline_k, test_k in zip(baseline_keys, test_keys):
        assert baseline_k == test_k, "Key mismatch in loaded data."
        k: str = baseline_k
        #
        # Use a lambda to select the appropriate comparison operator
        #
        compare: Callable = lambda x, y: x == y if k in explicit_comparisons else np.allclose
        failure_message: str = "Numerics mismatch in data."
        assert compare(getattr(baseline, k), getattr(geqdsk_interpolator.geqdsk_data, k)), failure_message

def test_geqdsk_interpolator_initialization_from_path(data_dir: Path) -> None:
    """
    Test the __call__ method of the GEQDSKInterpolator class

    Parameters
    ----------
    data_dir : pathlib.Path
        Pytest fixture for the data dir
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True) 
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    assert geqdsk_interpolator is not None
    assert geqdsk_interpolator.geqdsk_data is not None

def test_geqdsk_interpolator_refine_with_root_finder(data_dir: Path) -> None:
    """
    Test the refine_with_root_finder method of the GEQDSKInterpolator class
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True) 
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    #
    # Baseline magx coordinates
    #
    baseline_rmagx: float = 1.7578560423053675
    baseline_zmagx: float = -0.029247875129327684
    #
    # Baseline xpt coordinates
    #
    baseline_rxpt: float = 1.300088209826728
    baseline_zxpt: float = -1.133074411967326
    #
    # Define the guess coordinates
    #
    guess_rmagx: float = 1.75
    guess_zmagx: float = -0.02
    guess_rxpt: float = 1.30
    guess_zxpt: float = -1.13
    #
    # Refine the coordinates
    #
    refined_rmagx, refined_zmagx = geqdsk_interpolator.refine_with_root_finder(guess_rmagx, guess_zmagx)
    refined_rxpt, refined_zxpt = geqdsk_interpolator.refine_with_root_finder(guess_rxpt, guess_zxpt)
    assert np.isclose(refined_rmagx, baseline_rmagx)
    assert np.isclose(refined_zmagx, baseline_zmagx)
    assert np.isclose(refined_rxpt, baseline_rxpt)
    assert np.isclose(refined_zxpt, baseline_zxpt)

def test_geqdsk_interpolator_generates_normalized_interpolator(data_dir: Path) -> None:
    """
    Test the generate_normalized_interpolator method of the GEQDSKInterpolator class
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True) 
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    #
    # Baseline magx coordinates
    #
    rmagx: float = 1.7578560423053675
    zmagx: float = -0.029247875129327684
    #
    # Baseline xpt coordinates
    #
    rxpt: float = 1.300088209826728
    zxpt: float = -1.133074411967326
    #
    # Get psi at the magnetic axis and x-point
    #
    magx_psi: float = geqdsk_interpolator(rmagx, zmagx)
    xpt_psi: float = geqdsk_interpolator(rxpt, zxpt)
    #
    # Generate a normalized interpolator
    #
    normalized_geqdsk_interpolator = geqdsk_interpolator.generate_normalized_interpolator(psi_min=magx_psi, psi_max=xpt_psi)
    #
    # Check that the normalized interpolator returns the correct values
    #
    magx_psi: float = normalized_geqdsk_interpolator(rmagx, zmagx)
    assert np.isclose(magx_psi, 0.0)
    #
    # X-point separatrix psi should be one
    #
    xpt_psi: float = normalized_geqdsk_interpolator(rxpt, zxpt)
    assert np.isclose(xpt_psi, 1.0)

def test_geqdsk_interpolator_gradients_at_magx_and_xpt(data_dir: Path) -> None:
    """
    Test the normalized gradients of the GEQDSKInterpolator class
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True) 
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    #
    # Baseline magx coordinates
    #
    rmagx: float = 1.7578560423053675
    zmagx: float = -0.029247875129327684
    #
    # Baseline xpt coordinates
    #
    rxpt: float = 1.300088209826728
    zxpt: float = -1.133074411967326
    #
    # Get psi at the magnetic axis and x-point
    #
    magx_psi: float = geqdsk_interpolator(rmagx, zmagx)
    xpt_psi: float = geqdsk_interpolator(rxpt, zxpt)
    #
    # Generate a normalized interpolator
    #
    normalized_geqdsk_interpolator = geqdsk_interpolator.generate_normalized_interpolator(psi_min=magx_psi, psi_max=xpt_psi)
    #
    # Check that the normalized interpolator returns the correct gradients at the magnetic axis and x-point
    #
    magx_psi_vz: float = normalized_geqdsk_interpolator(rmagx, zmagx, tag='vz')
    magx_psi_vr: float = normalized_geqdsk_interpolator(rmagx, zmagx, tag='vr')
    assert np.isclose(magx_psi_vz, 0.0)
    assert np.isclose(magx_psi_vr, 0.0)
    #
    # Check that the normalized interpolator returns the correct gradients at the x-point
    #
    xpt_psi_vz: float = normalized_geqdsk_interpolator(rxpt, zxpt, tag='vz')
    xpt_psi_vr: float = normalized_geqdsk_interpolator(rxpt, zxpt, tag='vr')
    assert np.isclose(xpt_psi_vz, 0.0)
    assert np.isclose(xpt_psi_vr, 0.0)

    #
    # Check that the original interpolator returns the correct gradients at the magnetic axis and x-point
    #
    magx_psi_vz: float = geqdsk_interpolator(rmagx, zmagx, tag='vz')
    magx_psi_vr: float = geqdsk_interpolator(rmagx, zmagx, tag='vr')
    assert np.isclose(magx_psi_vz, 0.0)
    assert np.isclose(magx_psi_vr, 0.0)

    #
    # Check that the original interpolator returns the correct gradients at the x-point
    #
    rxpt_psi_vz: float = geqdsk_interpolator(rxpt, zxpt, tag='vz')
    rxpt_psi_vr: float = geqdsk_interpolator(rxpt, zxpt, tag='vr')
    assert np.isclose(rxpt_psi_vz, 0.0)
    assert np.isclose(rxpt_psi_vr, 0.0)
