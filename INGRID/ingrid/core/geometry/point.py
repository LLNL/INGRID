import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin
from dataclasses import dataclass
from typing import Union, Iterable, Any

@dataclass(frozen=True)
class Point(NDArrayOperatorsMixin):
    """
    Define a Point in N-dimensional space that supports NumPy operations.

    Parameters
    ----------
    coordinates : np.ndarray
        An array representing the coordinates of the point.
    """
    coordinates: np.ndarray

    def __init__(self, *coordinates: Union[float, Iterable[float]]):
        object.__setattr__(self, 'coordinates', self._process_coordinates(coordinates))

    def _process_coordinates(self, coordinates: Any) -> np.ndarray:
        if len(coordinates) == 1 and isinstance(coordinates[0], Iterable):
            coords = np.array(coordinates[0], dtype=float)
        else:
            coords = np.array(coordinates, dtype=float)
        return coords

    def __repr__(self) -> str:
        return f'Point({self.coordinates.tolist()})'
    
    def __getitem__(self, item) -> np.ndarray:
        return self.coordinates[item]
    
    def __len__(self) -> int:
        return len(self.coordinates)

    def __array__(self, dtype=None) -> np.ndarray:
        return np.array(self.coordinates, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Convert inputs to arrays or scalars
        arrays = []
        for x in inputs:
            if isinstance(x, Point):
                arrays.append(x.coordinates)
            else:
                arrays.append(x)
        # Perform the operation
        result = getattr(ufunc, method)(*arrays, **kwargs)
        # Handle output
        if isinstance(result, np.ndarray) and result.shape == self.coordinates.shape:
            return Point(result)
        else:
            return result

    def distance_to(self, other, ord=2.0, axis=0) -> np.floating[Any]:
        return np.linalg.norm(self - other, ord=ord, axis=axis)
