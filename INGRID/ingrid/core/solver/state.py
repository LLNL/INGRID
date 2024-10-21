from dataclasses import dataclass, field
from typing import Union, Tuple, Optional
import numpy as np
from ingrid.core.geometry import Point

@dataclass
class SolverState:
    y_current: Union[Point, np.ndarray]
    y_prev: Union[Point, np.ndarray]
    function: callable
    tspan: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0]))
    count: int = field(default=0)
    running_time: float = field(default=0.0)