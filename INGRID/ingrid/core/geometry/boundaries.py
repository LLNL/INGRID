from typing import Literal, Union, Optional
import numpy as np
from ingrid.core.geometry.line import Line
from ingrid.core.geometry.exceptions import MisalignedBoundaryError

class ClosedBoundary:
    north: Union[Line, np.ndarray]
    south: Union[Line, np.ndarray]
    east: Union[Line, np.ndarray]
    west: Union[Line, np.ndarray]
    atol: float
    rtol: float

    def __init__(self,
                 *,
                 north: Union[Line, np.ndarray],
                 south: Union[Line, np.ndarray],
                 east: Union[Line, np.ndarray],
                 west: Union[Line, np.ndarray],
                 atol: Optional[float] = 1e-9,
                 rtol: Optional[float] = 1e-9) -> None:

        self.north = north
        self.south = south
        self.east = east
        self.west = west
        self.atol = atol
        self.rtol = rtol

        self._validate_aligned_boundary(atol=atol, rtol=rtol)
        self._remove_duplicate_points()

    def _validate_aligned_boundary(self, atol: Optional[float] = None, rtol: Optional[float] = None) -> None:
        atol = self.atol if atol is None else atol
        rtol = self.rtol if rtol is None else rtol

        aligned_corners: list[bool] = [
            # Check north-east corner alignment
            np.allclose(self.north[-1], self.east[0], rtol=self.rtol, atol=self.atol),
            # Check south-east corner alignment
            np.allclose(self.east[-1], self.south[0], rtol=self.rtol, atol=self.atol),
            # Check south-west corner alignment
            np.allclose(self.south[-1], self.west[0], rtol=self.rtol, atol=self.atol),
            # Check north-west corner alignment
            np.allclose(self.west[-1], self.north[0], rtol=self.rtol, atol=self.atol)
        ]
        if not all(aligned_corners):
            corners: tuple[str, str, str, str] = ('NE', 'SE', 'SW', 'NW')
            misaligned_corners: str = ', '.join(corners[i] for i in range(4) if not aligned_corners[i])
            raise MisalignedBoundaryError(f"Misaligned BoundaryLine corners detected: {misaligned_corners}")
        
    def _remove_duplicate_points(self) -> None:
        def unique(values: Union[Line, np.ndarray]) -> Union[Line, np.ndarray]:
            _, idx = np.unique(values, return_index=True, axis=0)
            return values[np.sort(a=idx)]
        self.north = unique(values=self.north)
        self.south = unique(values=self.south)
        self.east = unique(values=self.east)
        self.west = unique(values=self.west)
