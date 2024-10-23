import pytest
import numpy as np
from pathlib import Path
from ingrid.core.geometry import Point, Line
from ingrid.core.interpol import GEQDSKInterpolator
from ingrid.utils.helper import resolve_path
from ingrid.core.solver.experimental_trace import TracerSettings, TracerOptions, TerminationCriteria, IntegratorSettings, LineTracer
from ingrid.core.solver.convergence import ConvergenceSettings
from ingrid.core.solver.exceptions import BoundaryCrossedError

@pytest.fixture
def xpt() -> Point:
    #
    # Xpt coordinates
    #
    rxpt: float = 1.300088209826728
    zxpt: float = -1.133074411967326
    return Point([rxpt, zxpt])

@pytest.fixture
def magx() -> Point:
    #
    # Magx coordinates
    #
    rmagx: float = 1.7578560423053675
    zmagx: float = -0.029247875129327684
    return Point([rmagx, zmagx])

@pytest.fixture
def normalized_geqdsk_interpolator(data_dir: Path, xpt: Point, magx: Point) -> GEQDSKInterpolator:
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True)
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    magx_psi: float = geqdsk_interpolator(*magx)
    xpt_psi: float = geqdsk_interpolator(*xpt)
    normalized_geqdsk_interpolator = geqdsk_interpolator.generate_normalized_interpolator(psi_min=magx_psi, psi_max=xpt_psi)
    return normalized_geqdsk_interpolator

def test_tracing_SNL_core_plasma_separatrix_with_point_convergence(normalized_geqdsk_interpolator: GEQDSKInterpolator, xpt: Point, magx: Point) -> None:
    """
    Test tracing the core plasma separatrix in SNL with point convergence.
    """
    primary_xpoint = xpt

    tracer_settings = TracerSettings(
        option=TracerOptions.THETA,
        direction='cw',
        termination_criteria=TerminationCriteria.POINT_CONVERGENCE,
        terminal_event=primary_xpoint,
        geqdsk_interpolator=normalized_geqdsk_interpolator
    )

    integrator_settings = IntegratorSettings(
        atol=1e-10,
        rtol=1e-10,
        dt=0.01,
        first_step=1e-5
    )

    convergence_settings = ConvergenceSettings(
        atol=1e-2,
        rtol=1e-2,
        norm=1,
        escape_epsilon=1e-5
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings,
        convergence_settings=convergence_settings
    )

    traced_points: np.ndarray = tracer.trace(
        start_point=primary_xpoint + integrator_settings.first_step,  # Start a bit away from the x-point to avoid singularity during integration
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

    #
    # Ensure that we have traced a non-trivial line (i.e. no early termination)
    #
    assert tracer.last_solver_state.count > 100

def test_tracing_SNL_leg_with_line_intersection(normalized_geqdsk_interpolator: GEQDSKInterpolator, xpt: Point, magx: Point) -> None:
    """
    Test tracing a leg in SNL with point convergence.
    """
    SE_target_plate: Line = Line(np.loadtxt(SE_target_plate_path, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift

    tracer_settings = TracerSettings(
        option=TracerOptions.THETA,
        direction='cw',
        termination_criteria=TerminationCriteria.LINE_INTERSECTION,
        terminal_event=SE_target_plate,
        geqdsk_interpolator=normalized_geqdsk_interpolator
    )

    integrator_settings = IntegratorSettings(
        atol=1e-10,
        rtol=1e-10,
        dt=0.01,
        first_step=1e-5
    )

    convergence_settings = ConvergenceSettings(
        atol=1e-8,
        rtol=1e-8,
        norm=1
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings,
        convergence_settings=convergence_settings
    )

    traced_points: np.ndarray = tracer.trace(
        start_point=xpt - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
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
    assert tracer.last_solver_state.count > 20

def test_tracing_SNL_leg_with_line_group_intersection(normalized_geqdsk_interpolator: GEQDSKInterpolator, xpt: Point, magx: Point) -> None:
    """
    Test tracing a leg in SNL with line group convergence.
    """
    SW_target_plate: Line = Line(np.loadtxt(SW_target_plate_path, delimiter=',') - np.array([0.5, 1.6]))  # Applying an rz-shift
    SE_target_plate: Line = Line(np.loadtxt(SE_target_plate_path, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift

    tracer_settings = TracerSettings(
        option=TracerOptions.THETA,
        direction='cw',
        termination_criteria=TerminationCriteria.LINE_GROUP_INTERSECTION,
        terminal_event=[SW_target_plate, SE_target_plate],
        geqdsk_interpolator=normalized_geqdsk_interpolator
    )

    convergence_settings = ConvergenceSettings(
        atol=1e-8,
        rtol=1e-8,
        norm=1
    )

    integrator_settings = IntegratorSettings(
        atol=1e-10,
        rtol=1e-10,
        dt=0.01,
        first_step=1e-5
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings,
        convergence_settings=convergence_settings
    )

    traced_points: np.ndarray = tracer.trace(
        start_point=xpt - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
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
    assert tracer.last_solver_state.count > 20
    #
    # The tracer should intersect the SE target plate. This corresponds to index 1 of
    # the terminal_event list
    #
    assert tracer.last_convergence_checker.cache['event_index'] == 1

def test_tracing_SNL_leg_intersects_boundaries(normalized_geqdsk_interpolator: GEQDSKInterpolator, xpt: Point, magx: Point) -> None:
    """
    Test tracing a leg in without target plates in the path will result in
    going out of bounds
    """
    SW_target_plate: Line = Line(np.loadtxt(SW_target_plate_path, delimiter=','))  # Apply no shifts
    SE_target_plate: Line = Line(np.loadtxt(SE_target_plate_path, delimiter=','))  # Apply no shifts

    tracer_settings = TracerSettings(
        option=TracerOptions.THETA,
        direction='cw',
        termination_criteria=TerminationCriteria.LINE_GROUP_INTERSECTION,
        terminal_event=[SW_target_plate, SE_target_plate],
        geqdsk_interpolator=normalized_geqdsk_interpolator
    )

    integrator_settings = IntegratorSettings(
        atol=1e-10,
        rtol=1e-10,
        dt=0.01,
        first_step=1e-5
    )

    convergence_settings = ConvergenceSettings(
        atol=1e-8,
        rtol=1e-8,
        norm=1
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings,
        convergence_settings=convergence_settings
    )

    with pytest.raises(BoundaryCrossedError):
        traced_points: np.ndarray = tracer.trace(
            start_point=xpt - integrator_settings.first_step,
            tracer_settings=tracer_settings
        )

    with pytest.raises(BoundaryCrossedError):
        tracer_settings.direction = 'ccw'  # Trace in the other direction as well
        traced_points: np.ndarray = tracer.trace(
            start_point=xpt - integrator_settings.first_step,
            tracer_settings=tracer_settings
        )

def test_tracing_SNL_psi_convergence(normalized_geqdsk_interpolator: GEQDSKInterpolator, xpt: Point, magx: Point) -> None:
    """
    Test tracing to psi value defining the SOL boundary
    """
    SE_target_plate: Line = Line(np.loadtxt(SE_target_plate_path, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift
    psi_1: float = 1.1

    tracer_settings = TracerSettings(
        option=TracerOptions.RHO,
        direction='ccw',
        termination_criteria=TerminationCriteria.PSI_VALUE,
        terminal_event=psi_1,
        geqdsk_interpolator=normalized_geqdsk_interpolator
    )

    integrator_settings = IntegratorSettings(
        atol=1e-10,
        rtol=1e-10,
        dt=0.01,
        first_step=1e-5
    )

    convergence_settings = ConvergenceSettings(
        atol=1e-8,
        rtol=1e-8,
        norm=1
    )

    tracer: LineTracer = LineTracer(
        tracer_settings=tracer_settings, 
        integrator_settings=integrator_settings,
        convergence_settings=convergence_settings
    )

    segment_to_sol: np.ndarray = tracer.trace(
        start_point=xpt - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
        tracer_settings=tracer_settings
    )

    #
    # First let's ensure the line we traced in the RHO direction terminated at the SOL boundary psi_1
    #
    assert np.isclose(normalized_geqdsk_interpolator(segment_to_sol[-1, 0], segment_to_sol[-1, 1]), psi_1)

    #
    # Update the tracer settings to trace in the poloidal THETA direction, terminating at the SE target plate
    #
    tracer_settings.option = TracerOptions.THETA
    tracer_settings.termination_criteria = TerminationCriteria.LINE_INTERSECTION
    tracer_settings.terminal_event = SE_target_plate
    tracer_settings.direction = 'cw'

    sol_boundary: np.ndarray = tracer.trace(
        start_point=segment_to_sol[-1],  # Start a bit away from the x-point in the approx SW direction
        tracer_settings=tracer_settings
    )

    #
    # Check that the traced points are all on SOL boundary we defined as terminal in the first tracing
    #
    psi_values: np.ndarray = np.array([
        normalized_geqdsk_interpolator(r, z)
        for r, z in zip(sol_boundary[:, 0], sol_boundary[:, 1])
    ])

    assert np.allclose(psi_values,  psi_1)
    assert tracer.last_solver_state.count > 20
