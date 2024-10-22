from abc import ABC, abstractmethod
from typing import Union, Optional, Callable, List, Dict, Any
from dataclasses import dataclass, replace
import numpy as np
from ingrid.core.geometry import Point, Line
from ingrid.core.geometry.utils import truncate_last_segment
from ingrid.core.solver.state import SolverState

@dataclass
class ConvergenceSettings:
    """
    Settings for convergence tolerances and behaviors.
    
    Attributes
    ----------
    norm: int, float
        The norm to pass to the linalg norm
    atol: float
        Absolute tolerance for identifying convergence
    rtol: float
        Relative tolerance for identifying convergence 
    """
    norm: Union[int, float] = 1
    atol: float = 1e-6
    rtol: float = 1e-6
    escape_epsilon: float = 1e-6

class ConvergenceCriteria(ABC):
    convergence_checker: ConvergenceSettings
    cache: Dict[str, Any]

    @abstractmethod
    def is_converged(self, solver_state: SolverState, **kwargs) -> bool:
        ...

class PointConvergence(ConvergenceCriteria):
    terminal_event: Union[Point, np.ndarray]

    def __init__(self, tracer_settings: 'TracerSettings', convergence_settings: ConvergenceSettings) -> None:
        self.terminal_event = tracer_settings.terminal_event
        self.convergence_settings = convergence_settings
        self.escaped = False
    
    def is_converged(self, solver_state: SolverState) -> bool:
        distance: np.floating[np.Any] = np.linalg.norm(x=self.terminal_event - solver_state.y_current)
        # print(f'{distance = }')
        if self.escaped:
            return distance <= self.convergence_settings.atol
        else:
            if distance >= self.convergence_settings.escape_epsilon:
                self.escaped = True
            return False

class LineIntersection(ConvergenceCriteria):
    terminal_event: Union[Line, np.ndarray]

    def __init__(self, tracer_settings: 'TracerSettings', convergence_settings: ConvergenceSettings) -> None:
        self.terminal_event = tracer_settings.terminal_event
        self.convergence_settings = convergence_settings
        self.cache = {}
    
    def is_converged(self, solver_state: SolverState) -> bool:
        test_segment: Line = Line([solver_state.y_prev, solver_state.y_current])
        result: Union[Line, None] = truncate_last_segment(last_segment=test_segment, other=self.terminal_event)
        self.cache['latest_truncation'] = result
        return result is not None
    
class LineGroupIntersection(LineIntersection):
    terminal_event: List[Union[Line, np.ndarray]]

    def __init__(self, tracer_settings: 'TracerSettings', convergence_settings: ConvergenceSettings) -> None:
        self.terminal_event = tracer_settings.terminal_event
        self.tracer_settings = tracer_settings
        self.convergence_settings = convergence_settings
        self.cache = {}
    
    def is_converged(self, solver_state: SolverState) -> bool:
        for event_index, line in enumerate(self.terminal_event):
            tracer_settings: 'TracerSettings' = replace(self.tracer_settings)  # Make a fresh copy
            tracer_settings.terminal_event = line
            convergence_checker = LineIntersection(
                tracer_settings=tracer_settings,
                convergence_settings=self.convergence_settings
            )
            if convergence_checker.is_converged(solver_state=solver_state):
                self.cache['latest_truncation'] = convergence_checker.cache['latest_truncation']
                self.cache['event_index'] = event_index
                return True
        return False

class PsiConvergence(ConvergenceCriteria):
    terminal_event: Union[float, np.ndarray]

    def __init__(self, tracer_settings: 'TracerSettings', convergence_settings: ConvergenceSettings) -> None:
        self.terminal_event = tracer_settings.terminal_event
        self.geqdsk_interpolator = tracer_settings.geqdsk_interpolator
        self.convergence_settings = convergence_settings
        self.cache = {}

    def is_converged(self, solver_state: SolverState) -> bool:
        x1, y1 = solver_state.y_prev[0], solver_state.y_prev[1]
        x2, y2 = solver_state.y_current[0], solver_state.y_current[1]

        psi1 = self.geqdsk_interpolator(x1, y1)
        psi2 = self.geqdsk_interpolator(x2, y2)

        test_psi = (psi1 - self.terminal_event) * (psi2 - self.terminal_event)
        return test_psi < 0
