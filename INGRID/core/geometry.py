import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin
from dataclasses import dataclass
from typing import ListUnion, Iterable, Any

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


class Line(NDArrayOperatorsMixin):
    """
    Define a Line composed of multiple Point instances, supporting NumPy operations.

    Parameters
    ----------
    points : Iterable[Union[Point, Iterable[float]]]
        An iterable containing Point instances or coordinate iterables.
    """

    def __init__(self, points: Iterable[Union[Point, Iterable[float]]]):
        # Process the points and store them as a list of Point instances
        self.points = []
        for p in points:
            if isinstance(p, Point):
                self.points.append(p)
            else:
                self.points.append(Point(p))
        # Create an array representation for efficient computations
        self.array = np.array([p.coordinates for p in self.points])

    def __repr__(self) -> str:
        return f'Line({self.points})'

    def __getitem__(self, key):
        return self.points[key]

    def __len__(self):
        return len(self.points)

    def __iter__(self):
        return iter(self.points)

    def __array__(self, dtype=None) -> np.ndarray:
        return np.array(self.array, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        # Convert Line instances to arrays
        arrays = []
        for x in inputs:
            if isinstance(x, Line):
                arrays.append(x.array)
            else:
                arrays.append(x)
        # Perform the ufunc operation
        result = getattr(ufunc, method)(*arrays, **kwargs)
        # Handle the output
        if isinstance(result, np.ndarray) and result.shape == self.array.shape:
            # Create new Points from the result coordinates
            new_points = [Point(coords) for coords in result]
            return Line(new_points)
        else:
            return result

    def calculate_length(self) -> float:
        """
        Calculate the total arclength of the Line.

        Returns
        -------
        float
            The sum of the lengths of all segments in the Line.
        """
        diffs = np.diff(self.array, axis=0)
        segment_lengths = np.linalg.norm(diffs, axis=1)
        return np.sum(segment_lengths)

    def calculate_cumulative_length(self) -> np.ndarray:
        """
        Calculate the cumulative arclength at each point in the Line.

        Returns
        -------
        np.ndarray
            An array containing the cumulative lengths.
        """
        diffs = np.diff(self.array, axis=0)
        segment_lengths = np.linalg.norm(diffs, axis=1)
        cumulative_lengths = np.cumsum(segment_lengths)
        return cumulative_lengths

    def at(self, parametric_t: float) -> Point:
        """
        Compute the Point at a given parametric value along the Line.

        Parameters
        ----------
        parametric_t : float
            A value between 0.0 and 1.0 representing the position along the Line.

        Returns
        -------
        Point
            The Point at the specified parametric position.
        """
        parametric_t = np.clip(parametric_t, 0.0, 1.0)
        total_length = self.calculate_length()
        target_length = parametric_t * total_length
        cumulative_lengths = self.calculate_cumulative_length()
        # Find the segment where the target length falls
        idx = np.searchsorted(cumulative_lengths, target_length)
        if idx == 0:
            frac = target_length / cumulative_lengths[0]
            p1 = self.points[0]
            p2 = self.points[1]
        elif idx >= len(cumulative_lengths):
            return self.points[-1]
        else:
            segment_length = cumulative_lengths[idx] - cumulative_lengths[idx - 1]
            frac = (target_length - cumulative_lengths[idx - 1]) / segment_length
            p1 = self.points[idx]
            p2 = self.points[idx + 1]
        # Linear interpolation between points
        coords = (1 - frac) * p1.coordinates + frac * p2.coordinates
        return Point(coords)

    def __add__(self, other):
        """
        Concatenate two Line instances.

        Parameters
        ----------
        other : Line
            The other Line to concatenate.

        Returns
        -------
        Line
            A new Line containing points from both Lines.
        """
        if isinstance(other, Line):
            new_points = self.points + other.points
            return Line(new_points)
        else:
            return NotImplemented

    def __radd__(self, other):
        """
        Support for sum() with Line instances.

        Parameters
        ----------
        other : int or Line
            The initial value or another Line.

        Returns
        -------
        Line
            The resulting Line after addition.
        """
        if other == 0:
            return self
        elif isinstance(other, Line):
            return self.__add__(other)
        else:
            return NotImplemented



# import numpy as np
# import numpy.lib.mixins
# from collections.abc import Iterable
# from numbers import Number

# class Point(np.lib.mixins.NDArrayOperatorsMixin):
#     """
#     Define a Point.
#     Can be used to later define Line objects.

#     Inherits from `numpy.lib.mixins.NDArrayOperatorsMixin` to allow
#     numpy operations directly on a Point object

#     Parameters
#     ----------
#     *coordinates : array-like, positional
#         Accepts either an array like container of numerical values or
#         numerical arguments as positional arguments
    
#     """

#     def __init__(self, *coordinates, **kargs) -> None:

#         if len(coordinates) == 0:
#             # Check kargs
#             coordinates = kargs.get('coordinates', None)
#             if coordinates is None:
#                 error_message  = 'No positional or keyword arguments for '
#                 error_message += '"coordinates" were provided'
#                 raise ValueError(error_message)
#             coordinates = [coordinates]
#         if len(coordinates) == 1:
#             error_message  = 'The provided argument must be an iterable '
#             error_message += 'defining Point coordinates.\n'
#             error_message += f'Recieved type = {type(coordinates[0])}'
#             assert isinstance(coordinates[0], Iterable), error_message
#             self.coordinates = coordinates[0]
#         else:
#             self.coordinates = [coordinate for coordinate in coordinates]

#         self.array = self.__array__(dtype=kargs.get('dtype', None))
#         self.shape = self.array.shape

#     def __getitem__(self, slice, **kargs):
#         return self.__array__(kargs.get('dtype', None))[slice]
        
#     def __repr__(self) -> str:
#         return f'{self.__class__.__name__}{(*self.coordinates, )}'

#     def __str__(self) -> str:
#         return self.__repr__()
    
#     def __array__(self, dtype=None):
#         return np.array(self.coordinates, dtype=dtype)

#     def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
#         if method == '__call__':
#             scalars = []
#             for input in inputs:
#                 if isinstance(input, Number):
#                     scalars.append(input)
#                 elif isinstance(input, self.__class__):
#                     scalars.append(input.array)
#                 elif isinstance(input, np.ndarray):
#                     scalars.append(input)
#                 else:
#                     raise NotImplementedError
#             return self.__class__(ufunc(*scalars, **kwargs))
#         else:
#             return NotImplemented

#     def distance_to(self, point):
#         """
#         Compute the distance in unitless coordinates to another Point

#         Parameters:
#         -----------
#         point: Point, np.ndarray
#             Another Point object or numpy array of same shape
        
#         Returns:
#         --------
#             Euclidean distance to the other point
#         """
#         return np.linalg.norm(self - point)


# class Line:
#     """
#     Define a Line.
#     Can be used to later define Patch objects.

#     Parameters
#     ----------
#     points : array-like
#         An array like container containing either numerical values or Point
#         objects
    
#     """

#     def __init__(self, points, **kargs) -> None:
#         error_message  = 'The provided argument must be an iterable '
#         error_message += 'containing coordinates/Points defining a curve.\n'
#         error_message += f'Recieved type = {type(points[0])}'
#         assert isinstance(points, Iterable), error_message
        
#         if all(isinstance(point, Point) for point in points):
#             points      = points
#             array       = np.vstack([point.array for point in points])
#         elif all(isinstance(point, np.ndarray) for point in points):
#             points      = [Point(point) for point in points]
#             array       = points
#         else:
#             raise ValueError('Non-homogenous initializer list was provided.')

#         self.points = points
#         self.array  = array
#         self.cumulative_length = self.calculate_cumulative_length()
#         self.normalized_cumulative_length = \
#             self.cumulative_length / self.calculate_length()


#     def __getitem__(self, val, **kargs):
#         """
#         Define slicing/indexing capabilities
#         """
#         if not isinstance(val, (slice, Iterable, int)):
#             raise NotImplementedError
#         return self.points[val]

#     def __iter__(self):
#         """
#         Dispatch to the list of points
#         """
#         return self.points.__iter__()

#     def __next__(self):
#         """
#         Dispatch to the list of points
#         """
#         return self.points.__next__()
        
#     def __repr__(self) -> str:
#         _repr  = f'{self.__class__.__name__}['
#         _repr += '\n'.join('\t' + repr(p) for p in self.points)
#         _repr += '\n' + len(f'{self.__class__.__name__}') * ' ' + ']'
#         return _repr

#     def __str__(self) -> str:
#         return self.__repr__()
    
#     def __array__(self, dtype=None):
#         return np.array(self.array, dtype=dtype)

#     def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
#         if method == '__call__':
#             scalars = []
#             for input in inputs:
#                 if isinstance(input, Number):
#                     scalars.append(input)
#                 elif isinstance(input, self.__class__):
#                     scalars.append(input.array)
#                 elif isinstance(input, np.ndarray):
#                     scalars.append(input)
#                 else:
#                     raise NotImplementedError
#             return self.__class__(ufunc(*scalars, **kwargs))
#         else:
#             return NotImplementedError

#     def __add__(self, other):
#         """
#         Concatenate two line objects
#         """
#         if not isinstance(other, self.__class__):
#             error_message  = 'Line addition can only be performed between '
#             error_message += 'Line types. Recieved type = "{type(other)}"'
#             raise ValueError(error_message)
        
#         new_points = self.points + other.points
#         return Line(points=new_points)

#     def calculate_length(self):
#         """
#         Calculate the arclength of the Line object

#         Returns:
#         --------
#             Float computing the sum of all curve segments
#         """
#         return np.linalg.norm(np.diff(self.array, axis=0), axis=1).sum()

#     def calculate_cumulative_length(self):
#         """
#         Calculate the cumulative arclength of the Line object

#         Returns:
#         --------
#             Numpy array containing cumulative lengths
#         """
#         return np.linalg.norm(np.diff(self.array, axis=0), axis=1).cumsum()
    
#     def at(self, parametric_t):
#         """
#         Interface to compute_parametric_point
#         """
#         return self.compute_parametric_point(parametric_t=parametric_t)

#     def compute_parametric_point(self, parametric_t):
#         """
#         Infer the Point which would reside along a parametric representation
#         of the Line object.

#         Given a parametric value p in [0., 1.], compute the coordinates via
#         linear interpolation along line segments

#         Parameters:
#         -----------
#         parametric_t: float
#             A numeric representing the input to the parametric curve equation
        
#         Returns:
#         --------
#             A Point object
#         """

#         parametric_t = np.clip(parametric_t, a_min=0., a_max=1.)

#         if np.allclose(parametric_t, 1.):
#             result = self.points[-1]
#         elif np.allclose(parametric_t, 0.):
#             result = self.points[0]
#         else:
#             #
#             # Find indices where parametric request is less than normalized
#             # cumsum of line segments. 
#             #
#             # This allows us to find the Point in the Line *after* where our 
#             # parametric request should lie
#             #
#             where_below = \
#                 np.where(parametric_t < self.normalized_cumulative_length)[0]
#             endpoint_index = where_below[0]
#             endpoint = self.points[endpoint_index]
            
#             endpoint_parametric_t = \
#                 self.normalized_cumulative_length[endpoint_index]

#             parametric_t_diff = endpoint_parametric_t - parametric_t
            
#             #
#             # Compute new coordinates via linear interpolation
#             #
#             result = endpoint * (1 - parametric_t_diff)

#         return result

#     def nearest(self, parametric_t):
#         """
#         Find the Point closest to a parametric t value

#         Parameters:
#         -----------
#         parametric_t: float
#             A numeric representing the input to the parametric curve
#         """
         

# class Patch:
#     """
#     Define a Patch region 
#     """