from abc import ABC, abstractmethod
from typing import Union, Optional, Callable, List, Dict, Any
import numpy as np
from ingrid.core.geometry import Point, Line
from ingrid.core.geometry.utils import truncate_last_segment
from ingrid.core.solver.state import SolverState

class ConvergenceCriteria(ABC):
    atol: Optional[float]
    rtol: Optional[float]
    norm: Optional[Union[float, str]]
    cache: Dict[str, Any]

    @abstractmethod
    def is_converged(self, solver_state: SolverState, **kwargs) -> bool:
        ...

class PointConvergence(ConvergenceCriteria):
    terminal_event: Union[Point, np.ndarray]

    def __init__(self, terminal_event: Union[Point, np.ndarray], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_event = terminal_event
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
    
    def is_converged(self, solver_state: SolverState) -> bool:
        distance: np.floating[np.Any] = np.linalg.norm(x=self.terminal_event - solver_state.y_current)
        print(f'{distance = }')
        return distance <= self.atol

class LineIntersection(ConvergenceCriteria):
    terminal_event: Union[Line, np.ndarray]

    def __init__(self, terminal_event: Union[Line, np.ndarray], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_event = terminal_event
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
        self.cache = {}
    
    def is_converged(self, solver_state: SolverState) -> bool:
        test_segment: Line = Line([solver_state.y_prev, solver_state.y_current])
        result: Union[Line, None] = truncate_last_segment(last_segment=test_segment, other=self.terminal_event)
        self.cache['latest_truncation'] = result
        return result is not None
    
class LineGroupIntersection(LineIntersection):
    terminal_event: List[Union[Line, np.ndarray]]

    def __init__(self, terminal_event: List[Union[Line, np.ndarray]], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_event = terminal_event
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
        self.cache = {}
    
    def is_converged(self, solver_state: SolverState) -> bool:
        for line in self.terminal_event:
            convergence_checker = LineIntersection(terminal_event=line, atol=self.atol, rtol=self.rtol, norm=self.norm)
            if convergence_checker.is_converged(solver_state=solver_state):
                self.cache['latest_truncation'] = convergence_checker.cache['latest_truncation']
                return True
        return False

class PsiConvergence(ConvergenceCriteria):
    terminal_event: Union[float, np.ndarray]

    def __init__(self, terminal_event: Union[float, np.ndarray], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_event = terminal_event
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
        self.cache = {}

    def is_converged(self, solver_state: SolverState, scalar_func: Callable[..., float]) -> bool:
        x1, y1 = solver_state.y_prev[0], solver_state.y_prev[1]
        x2, y2 = solver_state.y_current[0], solver_state.y_current[1]

        psi1 = scalar_func(x1, y1)
        psi2 = scalar_func(x2, y2)

        return (psi1 - self.terminal_event) * (psi2 - self.terminal_event) < 0
