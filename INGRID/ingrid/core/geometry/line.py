from typing import Union, Iterable, Any
import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin
from ingrid.core.geometry.point import Point

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
        return self.array[key]

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return self.array.__iter__()

    def __array__(self, dtype=None) -> np.ndarray:
        return np.array(self.array, dtype=dtype)
    
    def __getattribute__(self, name: str) -> Any:
        try:
            return super().__getattribute__(name)
        except AttributeError:
            return getattr(self.array, name)

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
    
    def truncate_at_intersection(self, other: Union['Line', np.ndarray]) -> Union['Line', None]:
        """
        Compute the truncation of this line with another line.
        
        Given a curve C parameterized with s, the truncation is defined as
        all points in C up to the intersection at parameterization value s'

        Parameters:
        ----------
        other: Another line

        Returns:
        --------
            A Line object or None if no intersection
        """
        #
        # TODO: REFACTOR ME!!!
        #
        for j, line1 in enumerate(iterable=zip(reversed(self.array[:-1]), reversed(self.array[1:]))):
            print(j)
            (xa, ya), (xb, yb) = line1

            for i in range(len(other) - 1):
                (xc, yc), (xd, yd) = other[i], other[i + 1]

                M = np.array([[xb - xa, -xd + xc], [yb - ya, -yd + yc]])
                r = np.array([xc - xa, yc - ya])
                try:
                    sol = np.linalg.solve(M, r)
                except np.linalg.LinAlgError:
                    continue

                if (sol[0] <= 1) and (sol[1] <= 1) \
                        and (sol[0] >= 0) and (sol[1] >= 0):
                    #return Line([Point((xc, yc)), Point((xc + sol[1] * (xd - xc), yc + sol[1] * (yd - yc)))])
                    print(xc, yc)
                    return Line(self.points[:i] + [Point((xc + sol[1] * (xd - xc), yc + sol[1] * (yd - yc)))])
        return None

    def __add__(self, other: Union['Line', np.ndarray]) -> 'Line':
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
