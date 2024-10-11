from dataclasses import dataclass
from typing import Union, Tuple, Optional
import numpy as np
from ingrid.core.geometry import Point

@dataclass()
class SolverState:
    y_current: Union[Point, np.ndarray]
    y_prev: Union[Point, np.ndarray]
    dt: float
    function: callable
    tspan: Tuple[float, float]
    first_step: Optional[float] = 1e-5
    step_ratio: Optional[float] = 1e-3
    max_step: Optional[float] = 1.0
    atol: Optional[float] = 1e-12
    rtol: Optional[float] = 1e-13
    count: int = 0
    running_time: float = 0.0