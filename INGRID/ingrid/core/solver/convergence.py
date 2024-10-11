from abc import ABC, abstractmethod
from typing import Union, Optional
import numpy as np
from ingrid.core.geometry import Point, Line
from ingrid.core.solver.state import SolverState

class ConvergenceCriteria(ABC):
    atol: Optional[float]
    rtol: Optional[float]
    norm: Optional[Union[float, str]]

    @abstractmethod
    def is_converged(self, solver_state: SolverState) -> bool:
        ...

class PointConvergence(ConvergenceCriteria):
    terminal_point: Union[Point, np.ndarray]

    def __init__(self, terminal_point: Union[Point, np.ndarray], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_point = terminal_point
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
    
    def is_converged(self, solver_state: SolverState) -> bool:
        distance: np.floating[np.Any] = np.linalg.norm(x=self.terminal_point - solver_state.y)
        return np.isclose(distance, 0.0, atol=self.atol, rtol=self.rtol)

    def update_terminal_point(self, new_terminal_point: Union[Point, np.ndarray]) -> None:
        self.terminal_point = new_terminal_point

class LineIntersection(ConvergenceCriteria):
    terminal_line: Union[Line, np.ndarray]

    def __init__(self, terminal_line: Union[Line, np.ndarray], atol=1.0e-9, rtol=1.0e-9, norm=2.) -> None:
        self.terminal_line = terminal_line
        self.atol = atol
        self.rtol = rtol
        self.norm = norm
    
    def is_converged(self, solver_state: SolverState) -> bool:
        test_segment = Line([solver_state.y_prev, solver_state.y_current])
        result = test_segment.compute_truncation(other=self.terminal_line)
        self.latest_truncation = result
        return result is not None
    
    def get_latest_truncation(self) -> Union[Line, None]:
        return self.latest_truncation

    def update_terminal_line(self, new_terminal_line: Union[Line, np.ndarray]) -> None:
        self.terminal_line = new_terminal_line
