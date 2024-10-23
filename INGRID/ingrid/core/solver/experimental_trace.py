import numpy as np
from typing import Optional, Union, List, Callable
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from time import time
from dataclasses import dataclass
import enum
from ingrid.core.geometry import Point, Line
from ingrid.core.interpol import GEQDSKInterpolator
from ingrid.core.solver.convergence import ConvergenceCriteria, PointConvergence, LineIntersection, PsiConvergence, LineGroupIntersection, ConvergenceSettings
from ingrid.core.solver.state import SolverState
from ingrid.core.solver.exceptions import BoundaryCrossedError

class TerminationCriteria(enum.Enum):
    """
    Enum for different termination criteria used in the tracing process.
    """
    LINE_INTERSECTION = enum.auto()
    POINT_CONVERGENCE = enum.auto()
    LINE_GROUP_INTERSECTION = enum.auto()
    PSI_VALUE = enum.auto()

class TracerOptions(enum.Enum):
    """
    Enum for different tracer options.
    """
    THETA = enum.auto()
    RHO = enum.auto()
    R_CONST = enum.auto()
    Z_CONST = enum.auto()
    RZ_CONST = enum.auto()

@dataclass
class TracerSettings:
    """
    Settings for the tracer, including termination criteria and target.
    
    Attributes
    ----------
    termination_criteria : TerminationCriteria
        The criteria for terminating the trace.
    terminal_event : Union[Point, Line, List[Line], float]
        The target for termination.
    option : str, optional
        The option for tracing, default is 'theta'.
    direction : str, optional
        The direction for tracing, default is 'cw' (clockwise).
    """
    geqdsk_interpolator: GEQDSKInterpolator
    termination_criteria: TerminationCriteria
    terminal_event: Union[Point, Line, List[Line], float]
    option: TracerOptions = TracerOptions.THETA
    direction: str = 'cw'

@dataclass
class IntegratorSettings:
    """
    Settings for the integrator used in the tracing process.
    
    Attributes
    ----------
    first_step : float
        The first step size for the integrator.
    step_ratio : float
        The ratio of step size.
    dt : float
        The time step size.
    max_step : float
        The maximum step size.
    dynamic_step_threshold : Optional[float]
        The threshold for dynamic step adjustment.
    atol : float
        Absolute tolerance for the integrator.
    rtol : float
        Relative tolerance for the integrator.
    method : str, optional
        The method used by the integrator, default is 'LSODA'.
    """
    first_step: float = 5e-5
    step_ratio: float = 0.02
    dt: float = 0.01
    max_step: float = 0.02
    dynamic_step_threshold: Optional[float] = None
    atol: float = 1e-10
    rtol: float = 1e-10
    method: str = 'LSODA'


class LineTracer:
    tracer_settings: TracerSettings
    integrator_settings: IntegratorSettings

    def __init__(self, tracer_settings: TracerSettings, integrator_settings: IntegratorSettings, convergence_settings: ConvergenceSettings) -> None:
        """
        Initialize the LineTracer with given tracer and integrator settings.

        Parameters
        ----------
        tracer_settings : TracerSettings
            Settings for the tracer.
        integrator_settings : IntegratorSettings
            Settings for the integrator.
        convergence_settings : ConvergenceSettings
            Settings for the convergence checker.
        """
        self.tracer_settings = tracer_settings
        self.integrator_settings = integrator_settings
        self.convergence_settings = convergence_settings

    def _differential_theta(self, t: float, xy: np.ndarray, geqdsk_interpolator: GEQDSKInterpolator) -> np.ndarray:
        """
        Coupled set of differential equations to trace the poloidal lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].
        efit_data : EfitData
            EFIT data for magnetic field calculations.

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        R: float
        Z: float
        R, Z = xy
        B_R: float = (1 / R) * geqdsk_interpolator.get_psi(R, Z, tag='vz')
        B_Z: float = -(1 / R) * geqdsk_interpolator.get_psi(R, Z, tag='vr')
        B: float = np.sqrt(B_R**2 + B_Z**2)
        dR: float = B_R / B
        dZ: float = B_Z / B
        
        coefficient: float = 1.0 if self.tracer_settings.direction == 'cw' else -1.0
        result: np.ndarray = np.array([dR, dZ]) * coefficient
        return result

    def _differential_rho(self, t: float, xy: np.ndarray, geqdsk_interpolator: GEQDSKInterpolator) -> np.ndarray:
        """
        Coupled set of differential equations to trace the radial lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].
        efit_data : EfitData
            EFIT data for magnetic field calculations.

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        R: float
        Z: float
        R, Z = xy
        B_R: float = (1 / R) * geqdsk_interpolator.get_psi(R, Z, tag='vz')
        B_Z: float = -(1 / R) * geqdsk_interpolator.get_psi(R, Z, tag='vr')
        B: float = np.sqrt(B_R**2 + B_Z**2)
        dR: float = B_Z / B
        dZ: float = -B_R / B
        
        coefficient: float = 1.0 if self.tracer_settings.direction == 'cw' else -1.0
        result: np.ndarray = np.array([dR, dZ]) * coefficient
        return result

    def _differential_r_const(self, t: float, xy: np.ndarray) -> np.ndarray:
        """
        Coupled set of differential equations to trace vertical lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        B_R: float = 0
        B_Z: float = 1
        B: float = np.sqrt(B_R**2 + B_Z**2)
        dR: float = B_R / B
        dZ: float = B_Z / B
        
        coefficient: float = 1.0 if self.tracer_settings.direction == 'cw' else -1.0
        result: np.ndarray = np.array([dR, dZ]) * coefficient
        return result

    def _differential_rz_const(self, t: float, xy: np.ndarray, tilt_angle: Optional[float] = None) -> np.ndarray:
        """
        Coupled set of differential equations to trace horizontal lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].
        tilt_angle : Optional[float], optional
            Tilt angle for the line. Defaults to None.

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        if tilt_angle is None:
            tilt_angle = self.tracer_settings.tilt_angle
        B_R: float = np.cos(tilt_angle)
        B_Z: float = np.sin(tilt_angle)
        B: float = np.sqrt(B_R**2 + B_Z**2)
        dR: float = B_R / B
        dZ: float = B_Z / B

        coefficient: float = 1.0 if self.tracer_settings.direction == 'cw' else -1.0
        result: np.ndarray = np.array([dR, dZ]) * coefficient
        return result
    
    def _differential_r_const(self, t: float, xy: np.ndarray) -> np.ndarray:
        """
        Coupled set of differential equations to trace vertical lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        return self._differential_rz_const(t=t, xy=xy, tilt_angle=np.pi / 2)
    
    def _differential_z_const(self, t: float, xy: np.ndarray) -> np.ndarray:
        """
        Coupled set of differential equations to trace horizontal lines.

        Parameters
        ----------
        t : float
            Time parameter.
        xy : np.ndarray
            Array containing the coordinates [R, Z].

        Returns
        -------
        np.ndarray
            Array containing the derivatives [dR, dZ].
        """
        return self._differential_rz_const(t=t, xy=xy, tilt_angle=0)

    def _get_function(self, tracer_settings: Optional[TracerSettings] = None) -> None:
        """
        Get the function to be used for tracing based on the tracer settings.

        Parameters
        ----------
        tracer_settings : Optional[TracerSettings], optional
            Settings for the tracer. Defaults to None.
        """
        if tracer_settings is None:
            tracer_settings = self.tracer_settings

        registered_functions: dict[str, Callable[..., np.ndarray]] = {
            TracerOptions.THETA: self._differential_theta,
            TracerOptions.R_CONST: self._differential_r_const,
            TracerOptions.Z_CONST: self._differential_z_const,
            TracerOptions.RHO: self._differential_rho,
            TracerOptions.RZ_CONST: self._differential_rz_const
        }
        try:
            to_utilize: Callable[..., np.ndarray] = registered_functions[tracer_settings.option]
        except KeyError:
            raise ValueError(f"Unsupported function requested: {tracer_settings.option = }")
        return to_utilize

    def configure_convergence_checker(self, convergence_settings: ConvergenceSettings, tracer_settings: Optional[TracerSettings] = None) -> ConvergenceCriteria:
        """
        Configure the convergence checker based on the tracer settings.

        Parameters
        ----------
        tracer_settings : Optional[TracerSettings], optional
            Settings for the tracer. Defaults to None.

        Returns
        -------
        ConvergenceCriteria
            The configured convergence checker.
        """
        if tracer_settings is None:
            tracer_settings = self.tracer_settings

        if tracer_settings.termination_criteria == TerminationCriteria.POINT_CONVERGENCE:
            message: str = '# Starting search for Point convergence...'
            if not isinstance(tracer_settings.terminal_event, (Point, np.ndarray)):
                raise ValueError(f"Invalid termination target type for {TerminationCriteria.POINT_CONVERGENCE.name} criteria. Received {tracer_settings.terminal_event = }")
            convergence_checker: ConvergenceCriteria = PointConvergence(tracer_settings=tracer_settings, convergence_settings=convergence_settings)
        elif tracer_settings.termination_criteria == TerminationCriteria.LINE_INTERSECTION:
            message = '# Starting search for Line intersection...'
            if not isinstance(tracer_settings.terminal_event, (Line, np.ndarray)):
                raise ValueError(f"Invalid termination target type for {TerminationCriteria.LINE_INTERSECTION.name} criteria. Received {tracer_settings.terminal_event = }")
            convergence_checker = LineIntersection(tracer_settings=tracer_settings, convergence_settings=convergence_settings)
        elif tracer_settings.termination_criteria == TerminationCriteria.LINE_GROUP_INTERSECTION:
            message = '# Starting search for Line group intersection...'
            if not isinstance(tracer_settings.terminal_event, list) or not all(isinstance(line, (Line, np.ndarray)) for line in tracer_settings.terminal_event):
                raise ValueError(f"Invalid termination target type for {TerminationCriteria.LINE_GROUP_INTERSECTION.name} criteria. Received {tracer_settings.terminal_event = }")
            convergence_checker = LineGroupIntersection(tracer_settings=tracer_settings, convergence_settings=convergence_settings)
        elif tracer_settings.termination_criteria == TerminationCriteria.PSI_VALUE:
            message = '# Starting search for Psi value convergence...'
            convergence_checker = PsiConvergence(tracer_settings=tracer_settings, convergence_settings=convergence_settings)
        else:
            raise ValueError(f"Unsupported termination criteria: {tracer_settings.termination_criteria.name = }. Required one of {TerminationCriteria.__members__.keys()}")

        print(message)
        return convergence_checker

    def trace(self, start_point: Point, tracer_settings: Optional[TracerSettings] = None, convergence_settings: Optional[ConvergenceSettings] = None) -> Line:
        """
        Trace a line starting from the given point using the provided EFIT data and tracer settings.

        Parameters
        ----------
        start_point : Point
            The starting point for the trace.
        tracer_settings : Optional[TracerSettings], optional
            Settings for the tracer. Defaults to None.
        convergence_settings : Optional[ConvergenceSettings], optional
            Settings for the convergence criteria. Defaults to None

        Returns
        -------
        Line
            The traced line.
        """
        if tracer_settings is None:
            tracer_settings = self.tracer_settings
        if convergence_settings is None:
            convergence_settings = self.convergence_settings

        convergence_checker: ConvergenceCriteria = self.configure_convergence_checker(tracer_settings=tracer_settings, convergence_settings=convergence_settings)
        max_range: float = max(tracer_settings.geqdsk_interpolator.rmax - tracer_settings.geqdsk_interpolator.rmin, tracer_settings.geqdsk_interpolator.zmax - tracer_settings.geqdsk_interpolator.zmin)
        self.integrator_settings.max_step = self.integrator_settings.step_ratio * max_range
        count: int = 0
        start: float = time()

        #
        # Use minimum of self.dt or dynamic_step value if provided.
        #
        dt: float = self.integrator_settings.dt
        if self.integrator_settings.dynamic_step_threshold is not None:
            dt = np.amin([self.integrator_settings.dt, self.integrator_settings.dynamic_step_threshold])
            if dt < self.integrator_settings.dt:
                print('Using dynamic_step value!\n')
        #
        # Initialize the solver state
        #
        solver_state: SolverState = SolverState(
            y_current=start_point,
            y_prev=start_point,
            function=self._get_function(tracer_settings=tracer_settings),
            tspan=np.array([0.0, dt])
        )

        bounding_box: np.ndarray = np.array([
            [tracer_settings.geqdsk_interpolator.rmin, tracer_settings.geqdsk_interpolator.zmin],
            [tracer_settings.geqdsk_interpolator.rmax, tracer_settings.geqdsk_interpolator.zmax]
            ])

        def in_bounds(solver_state: SolverState) -> bool:
            """
            Check if the current solver state is within the bounding box.

            Parameters
            ----------
            solver_state : SolverState
                The current state of the solver.

            Returns
            -------
            bool
                True if the solver state is within bounds, False otherwise.
            """
            test_point: np.ndarray = solver_state.y_current
            if (test_point < bounding_box)[0, :].any() or (test_point > bounding_box)[1, :].any():
                #
                # Determine if this is the correct behavior and if we should return False instead
                #
                raise BoundaryCrossedError(f'# Line intersected the boundary during tracing with TerminationCriteria.{tracer_settings.termination_criteria.name}')
            return True

        print(f"Beginning tracing with {tracer_settings.termination_criteria.name} criteria. Solver state and tracer settings:")
        print(f"{solver_state = }")
        print(f"{tracer_settings = }")

        generated_curve: List[np.ndarray] = []

        while in_bounds(solver_state=solver_state) and not convergence_checker.is_converged(solver_state=solver_state):
            first_iteration = False
            #
            # Solve the system of differential equations
            #
            sol = solve_ivp(
                            fun=solver_state.function, 
                            t_span=solver_state.tspan,
                            #t_eval=solver_state.tspan,  # Only store the solution at the endpoints
                            y0=solver_state.y_current,
                            method=self.integrator_settings.method,
                            first_step=self.integrator_settings.first_step, 
                            max_step=self.integrator_settings.max_step,
                            rtol=self.integrator_settings.rtol,
                            atol=self.integrator_settings.atol,
                            args=(tracer_settings.geqdsk_interpolator,)
                        )
            #
            # Update the solver state.
            #
            solver_state.y_prev = sol.y[:, 0]
            solver_state.y_current = sol.y[:, -1]
            solver_state.tspan += self.integrator_settings.dt
            # print(f'{solver_state = }')

            #
            # Store the current coordinates.
            #
            generated_curve.append(solver_state.y_current)
            solver_state.count += 1
        solver_state.running_time = time() - start
        print(f'Drew for {solver_state.running_time} seconds')
        print(f'Tracer took {solver_state.count} step(s)')

        print("Finalizing curve...")

        if tracer_settings.termination_criteria == TerminationCriteria.POINT_CONVERGENCE:
            #
            # Remove the last point since it may not be the exact termination point
            #
            _: np.ndarray = generated_curve.pop()
            #
            # Append the exact termination point
            #
            generated_curve.append(convergence_checker.terminal_event)
        elif tracer_settings.termination_criteria in [TerminationCriteria.LINE_INTERSECTION, TerminationCriteria.LINE_GROUP_INTERSECTION]:
            #
            # Remove the last point since it may not be the exact termination point
            #
            _: np.ndarray = generated_curve.pop()
            #
            # Trim the last segment if it intersects with the termination line for exact intersection
            #
            intersection_point: np.ndarray = convergence_checker.cache['latest_truncation'].points[-1]
            generated_curve.append(intersection_point)
        elif tracer_settings.termination_criteria == TerminationCriteria.PSI_VALUE:
            #
            # Remove the last point since it may not be the exact termination point
            #
            _: np.ndarray = generated_curve.pop()
            #
            # Compute the exact termination point via root finding
            #
            x1: float
            y1: float
            x2: float
            y2: float
            x1, y1 = solver_state.y_prev[0], solver_state.y_prev[1]
            x2, y2 = solver_state.y_current[0], solver_state.y_current[1]

            #
            # Refine the root finding to find the exact termination point
            #
            # NOTE: We need to be careful when computing slopes as we may encounter
            #       horizontal or vertical lines which cause division by zero.
            #

            #
            # Let's start by assuming the end point is the last point in the curve
            #
            x_end: float = x2 
            y_end: float = y2

            #
            # Define generic root finder
            #
            def find_root(f, bracket):
                return root_scalar(f, bracket=bracket, rtol=convergence_settings.rtol, method='bisect')

            #
            # Vertical line case
            #
            if np.isclose(x2, x1, atol=convergence_settings.atol):
                def refine_psi_root(y: float) -> float:
                    """
                    Refine the root finding for vertical lines.

                    Parameters
                    ----------
                    y : float
                        y-coordinate.

                    Returns
                    -------
                    float
                        Difference between terminal event and psi value.
                    """
                    return convergence_checker.terminal_event - tracer_settings.geqdsk_interpolator.get_psi(x2, y).squeeze()
                sol = find_root(refine_psi_root, bracket=[y1, y2])
                y_end = sol.root
            #
            # Horizontal line case
            #
            elif np.isclose(y2, y1, atol=convergence_settings.atol):
                def refine_psi_root(x: float) -> float:
                    """
                    Refine the root finding for horizontal lines.

                    Parameters
                    ----------
                    x : float
                        x-coordinate.

                    Returns
                    -------
                    float
                        Difference between terminal event and psi value.
                    """
                    return convergence_checker.terminal_event - tracer_settings.geqdsk_interpolator.get_psi(x, y2).squeeze()
                sol = find_root(refine_psi_root, bracket=[x1, x2])
                x_end = sol.root
            #
            # Every other case is a line with a non-degenerate slope
            #
            else:
                def refine_psi_root(x: float) -> float:
                    """
                    Refine the root finding for lines with non-degenerate slopes.

                    Parameters
                    ----------
                    x : float
                        x-coordinate.

                    Returns
                    -------
                    float
                        Difference between terminal event and psi value.
                    """
                    #
                    # Compute the corresponding y value for the given x value
                    #
                    y: float = (y2 - y1) / (x2 - x1) * (x - x1) + y1
                    #
                    # Compute the psi value at the given x, y coordinates
                    #
                    return convergence_checker.terminal_event - tracer_settings.geqdsk_interpolator.get_psi(x, y).squeeze()
                sol = find_root(refine_psi_root, bracket=[x1, x2])
                x_end = sol.root
                y_end = (y2 - y1) / (x2 - x1) * (x_end - x1) + y1
            generated_curve.append(np.array([x_end, y_end]))

        self.last_solver_state = solver_state
        self.last_convergence_checker = convergence_checker
        return Line(generated_curve)
