from __future__ import annotations
import pytest
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from ingrid.core.geometry import Point, Line
from ingrid.core.interpol import GEQDSKInterpolator
from ingrid.utils.helper import resolve_path
from ingrid.core.solver.experimental_trace import TracerSettings, TracerOptions, TerminationCriteria, IntegratorSettings, LineTracer
from ingrid.core.solver.convergence import ConvergenceSettings
from ingrid.core.solver.exceptions import BoundaryCrossedError

@dataclass
class SNLTracingResources:
    """
    Resources needed for Single Null (SNL) magnetic field line tracing tests.

    Attributes
    ----------
    primary_xpoint : Point
        Coordinates of the primary X-point.
    magnetic_axis : Point
        Coordinates of the magnetic axis.
    normalized_geqdsk_interpolator : GEQDSKInterpolator
        Interpolator for normalized magnetic flux values.
    plate_w1 : Path
        Path to the southwest target plate geometry file.
    plate_e1 : Path
        Path to the southeast target plate geometry file.
    """
    primary_xpoint: Point
    magnetic_axis: Point
    normalized_geqdsk_interpolator: GEQDSKInterpolator
    plate_w1: Path
    plate_e1: Path

@pytest.fixture
def primary_xpoint() -> Point:
    """
    Fixture providing the coordinates of the primary X-point.

    Returns
    -------
    Point
        X-point coordinates [R,Z].
    """
    #
    # Xpt coordinates
    #
    rxpt: float = 1.300088209826728
    zxpt: float = -1.133074411967326
    return Point([rxpt, zxpt])

@pytest.fixture
def magnetic_axis() -> Point:
    """
    Fixture providing the coordinates of the magnetic axis.

    Returns
    -------
    Point
        Magnetic axis coordinates [R,Z].
    """
    #
    # Magx coordinates
    #
    rmagx: float = 1.7578560423053675
    zmagx: float = -0.029247875129327684
    return Point([rmagx, zmagx])

@pytest.fixture
def normalized_geqdsk_interpolator(data_dir: Path, primary_xpoint: Point, magnetic_axis: Point) -> GEQDSKInterpolator:
    """
    Fixture providing a normalized GEQDSK interpolator.

    Parameters
    ----------
    data_dir : Path
        Base directory containing test data.
    primary_xpoint : Point
        Primary X-point coordinates.
    magnetic_axis : Point
        Magnetic axis coordinates.

    Returns
    -------
    GEQDSKInterpolator
        Interpolator normalized between magnetic axis and X-point psi values.
    """
    geqdsk_path: str = resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'neqdsk', as_str=True)
    geqdsk_interpolator = GEQDSKInterpolator(geqdsk_data=geqdsk_path)
    magx_psi: float = geqdsk_interpolator(*magnetic_axis)
    xpt_psi: float = geqdsk_interpolator(*primary_xpoint)
    normalized_geqdsk_interpolator = geqdsk_interpolator.generate_normalized_interpolator(psi_min=magx_psi, psi_max=xpt_psi)
    return normalized_geqdsk_interpolator

@pytest.fixture
def plate_w1(data_dir: Path) -> Path:
    """
    Fixture providing path to southwest target plate geometry file.

    Parameters
    ----------
    data_dir : Path
        Base directory containing test data.

    Returns
    -------
    Path
        Path to southwest target plate file.
    """
    return resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'd3d_itp.txt', as_str=True)

@pytest.fixture
def plate_e1(data_dir: Path) -> Path:
    """
    Fixture providing path to southeast target plate geometry file.

    Parameters
    ----------
    data_dir : Path
        Base directory containing test data.

    Returns
    -------
    Path
        Path to southeast target plate file.
    """
    return resolve_path(path=data_dir / 'SNL' / 'DIII-D' / 'd3d_otp.txt', as_str=True)

@pytest.fixture
def snl_tracing_resources(
        primary_xpoint: Point,
        magnetic_axis: Point,
        normalized_geqdsk_interpolator: GEQDSKInterpolator, 
        plate_w1: Path, 
        plate_e1: Path
    ) -> SNLTracingResources:
    """
    Fixture combining all resources needed for SNL tracing tests.

    Parameters
    ----------
    primary_xpoint : Point
        Primary X-point coordinates.
    magnetic_axis : Point
        Magnetic axis coordinates.
    normalized_geqdsk_interpolator : GEQDSKInterpolator
        Normalized flux interpolator.
    plate_w1 : Path
        Path to southwest target plate file.
    plate_e1 : Path
        Path to southeast target plate file.

    Returns
    -------
    SNLTracingResources
        Combined resources for SNL tracing tests.
    """
    return SNLTracingResources(
        primary_xpoint=primary_xpoint,
        magnetic_axis=magnetic_axis,
        normalized_geqdsk_interpolator=normalized_geqdsk_interpolator,
        plate_w1=plate_w1,
        plate_e1=plate_e1
    )

class TestSNLTracing:

    def test_tracing_SNL_core_plasma_separatrix_with_point_convergence(self, snl_tracing_resources: SNLTracingResources) -> None:
        """
        Test tracing the core plasma separatrix in SNL with point convergence.

        Parameters
        ----------
        snl_tracing_resources : SNLTracingResources
            Resources for SNL tracing.

        Returns
        -------
        None
        """

        tracer_settings = TracerSettings(
            option=TracerOptions.THETA,
            direction='cw',
            termination_criteria=TerminationCriteria.POINT_CONVERGENCE,
            terminal_event=snl_tracing_resources.primary_xpoint,
            geqdsk_interpolator=snl_tracing_resources.normalized_geqdsk_interpolator
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
            start_point=snl_tracing_resources.primary_xpoint + integrator_settings.first_step,  # Start a bit away from the x-point to avoid singularity during integration
            tracer_settings=tracer_settings
        )

        #
        # Check that the traced points are on the separatrix (psi = 1.0)
        # NOTE: The bivariate spline interpolation only supports vectorized inputs when they are in
        #       non-decreasing order. Therefore, we pass points in individually since the separatrix
        #       is not monotonic in r or z.
        #
        psi_values: np.ndarray = np.array([
            snl_tracing_resources.normalized_geqdsk_interpolator(r, z)
            for r, z in zip(traced_points[:, 0], traced_points[:, 1])
        ])

        assert np.allclose(psi_values,  1.0)

        #
        # Ensure that we have traced a non-trivial line (i.e. no early termination)
        #
        assert tracer.last_solver_state.count > 100

    def test_tracing_SNL_leg_with_line_intersection(self, snl_tracing_resources: 'SNLTracingResources') -> None:
        """
        Test tracing a leg in SNL with point convergence.

        Parameters
        ----------
        snl_tracing_resources : SNLTracingResources
            Resources for SNL tracing.

        Returns
        -------
        None
        """
        SE_target_plate: Line = Line(np.loadtxt(snl_tracing_resources.plate_e1, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift

        tracer_settings = TracerSettings(
            option=TracerOptions.THETA,
            direction='cw',
            termination_criteria=TerminationCriteria.LINE_INTERSECTION,
            terminal_event=SE_target_plate,
            geqdsk_interpolator=snl_tracing_resources.normalized_geqdsk_interpolator
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
            start_point=snl_tracing_resources.primary_xpoint - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
            tracer_settings=tracer_settings
        )

        #
        # Check that the traced points are on the separatrix (psi = 1.0)
        # NOTE: The bivariate spline interpolation only supports vectorized inputs when they are in
        #       non-decreasing order. Therefore, we pass points in individually since the separatrix
        #       is not monotonic in r or z.
        #
        psi_values: np.ndarray = np.array([
            snl_tracing_resources.normalized_geqdsk_interpolator(r, z)
            for r, z in zip(traced_points[:, 0], traced_points[:, 1])
        ])

        assert np.allclose(psi_values,  1.0)
        assert tracer.last_solver_state.count > 20

    def test_tracing_SNL_leg_with_line_group_intersection(self, snl_tracing_resources: 'SNLTracingResources') -> None:
        """
        Test tracing a leg in SNL with line group convergence.

        Parameters
        ----------
        snl_tracing_resources : SNLTracingResources
            Resources for SNL tracing.

        Returns
        -------
        None
        """
        plate_w1: Line = Line(np.loadtxt(snl_tracing_resources.plate_w1, delimiter=',') - np.array([0.5, 1.6]))  # Applying an rz-shift
        plate_e1: Line = Line(np.loadtxt(snl_tracing_resources.plate_e1, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift

        tracer_settings = TracerSettings(
            option=TracerOptions.THETA,
            direction='cw',
            termination_criteria=TerminationCriteria.LINE_GROUP_INTERSECTION,
            terminal_event=[plate_w1, plate_e1],
            geqdsk_interpolator=snl_tracing_resources.normalized_geqdsk_interpolator
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
            start_point=snl_tracing_resources.primary_xpoint - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
            tracer_settings=tracer_settings
        )

        #
        # Check that the traced points are on the separatrix (psi = 1.0)
        # NOTE: The bivariate spline interpolation only supports vectorized inputs when they are in
        #       non-decreasing order. Therefore, we pass points in individually since the separatrix
        #       is not monotonic in r or z.
        #
        psi_values: np.ndarray = np.array([
            snl_tracing_resources.normalized_geqdsk_interpolator(r, z)
            for r, z in zip(traced_points[:, 0], traced_points[:, 1])
        ])

        assert np.allclose(psi_values,  1.0)
        assert tracer.last_solver_state.count > 20
        #
        # The tracer should intersect the SE target plate. This corresponds to index 1 of
        # the terminal_event list
        #
        assert tracer.last_convergence_checker.cache['event_index'] == 1

    def test_tracing_SNL_leg_intersects_boundaries(self, snl_tracing_resources: 'SNLTracingResources') -> None:
        """
        Test tracing a leg in without target plates in the path will result in
        going out of bounds

        Parameters
        ----------
        snl_tracing_resources : SNLTracingResources
            Resources for SNL tracing.

        Returns
        -------
        None
        """
        plate_w1: Line = Line(np.loadtxt(snl_tracing_resources.plate_w1, delimiter=','))  # Apply no shifts
        plate_e1: Line = Line(np.loadtxt(snl_tracing_resources.plate_e1, delimiter=','))  # Apply no shifts

        tracer_settings = TracerSettings(
            option=TracerOptions.THETA,
            direction='cw',
            termination_criteria=TerminationCriteria.LINE_GROUP_INTERSECTION,
                terminal_event=[plate_w1, plate_e1],
            geqdsk_interpolator=snl_tracing_resources.normalized_geqdsk_interpolator
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
                start_point=snl_tracing_resources.primary_xpoint - integrator_settings.first_step,
                tracer_settings=tracer_settings
            )

        with pytest.raises(BoundaryCrossedError):
            tracer_settings.direction = 'ccw'  # Trace in the other direction as well
            traced_points: np.ndarray = tracer.trace(
                start_point=snl_tracing_resources.primary_xpoint - integrator_settings.first_step,
                tracer_settings=tracer_settings
            )

    def test_tracing_SNL_psi_convergence(self, snl_tracing_resources: 'SNLTracingResources') -> None:
        """
        Test tracing to psi value defining the SOL boundary

        Parameters
        ----------
        snl_tracing_resources : SNLTracingResources
            Resources for SNL tracing.

        Returns
        -------
        None
        """
        SE_target_plate: Line = Line(np.loadtxt(snl_tracing_resources.plate_e1, delimiter=',') - np.array([0.0, 1.6]))  # Applying a z-shift
        psi_1: float = 1.1

        tracer_settings = TracerSettings(
            option=TracerOptions.RHO,
            direction='ccw',
            termination_criteria=TerminationCriteria.PSI_VALUE,
            terminal_event=psi_1,
            geqdsk_interpolator=snl_tracing_resources.normalized_geqdsk_interpolator
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
            start_point=snl_tracing_resources.primary_xpoint - integrator_settings.first_step,  # Start a bit away from the x-point in the approx SW direction
            tracer_settings=tracer_settings
        )

        #
        # First let's ensure the line we traced in the RHO direction terminated at the SOL boundary psi_1
        #
        assert np.isclose(snl_tracing_resources.normalized_geqdsk_interpolator(segment_to_sol[-1, 0], segment_to_sol[-1, 1]), psi_1)

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
            snl_tracing_resources.normalized_geqdsk_interpolator(r, z)
            for r, z in zip(sol_boundary[:, 0], sol_boundary[:, 1])
        ])

        assert np.allclose(psi_values,  psi_1)
        assert tracer.last_solver_state.count > 20
