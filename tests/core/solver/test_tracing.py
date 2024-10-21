import numpy as np
from pathlib import Path
from ingrid.core.geometry import Point
from ingrid.core.interpol import GEQDSKInterpolator
from ingrid.utils.helper import resolve_path
from ingrid.core.solver.experimental_trace import TracerSettings, TracerOptions, TerminationCriteria, IntegratorSettings, LineTracer

def test_tracing_SNL_core_plasma_separatrix_with_point_convergence(data_dir: Path) -> None:
    """
    Test tracing the core plasma separatrix in SNL with point convergence.
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True)
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)

    #
    # Xpt coordinates
    #
    rxpt: float = 1.300088209826728
    zxpt: float = -1.133074411967326

    #
    # Baseline magx coordinates
    #
    rmagx: float = 1.7578560423053675
    zmagx: float = -0.029247875129327684

    magx_psi: float = geqdsk_interpolator(rmagx, zmagx)
    xpt_psi: float = geqdsk_interpolator(rxpt, zxpt)
    normalized_geqdsk_interpolator = geqdsk_interpolator.generate_normalized_interpolator(psi_min=magx_psi, psi_max=xpt_psi)

    primary_xpoint: Point = Point([rxpt, zxpt])
    tracer_settings = TracerSettings(
        option=TracerOptions.THETA,
        direction='cw',
        termination_criteria=TerminationCriteria.POINT_CONVERGENCE,
        termination_target=primary_xpoint,
        args=(normalized_geqdsk_interpolator,),
        convergence_atol=1e-3,
        convergence_rtol=1e-3
    )

    integrator_settings = IntegratorSettings(
        atol=1e-12,
        rtol=1e-13,
        dt=0.001,
        first_step=1e-5
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings
    )

    traced_points: np.ndarray = tracer.trace(
        start_point=primary_xpoint + 1e-5,  # Start a bit away from the x-point to avoid singularity during integration
        geqdsk_interpolator=normalized_geqdsk_interpolator,
        tracer_settings=tracer_settings
    )

    #
    # Check that the traced points are on the separatrix (psi = 1.0)
    # NOTE: The bivariate spline interpolation only supports vectorized inputs when they are in
    #       non-decreasing order. Therefore, we pass points in individually since the separatrix
    #       is not monotonic in r or z.
    #
    psi_values: np.ndarray = np.array([
        normalized_geqdsk_interpolator(r, z)
        for r, z in zip(traced_points[:, 0], traced_points[:, 1])
    ])
    assert np.allclose(psi_values,  1.0)
