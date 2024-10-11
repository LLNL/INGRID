from typing import Union, List, Tuple
import numpy as np
from ingrid.core.geometry.line import Line
from ingrid.core.geometry.line import Point

def compute_bounding_box_from_points(points: Union[np.ndarray, List[Point]]) -> Tuple[float, float, float, float]:
    """
    Compute the bounding box from a set of points.

    This will compute the minimum and maximum x and y coordinates.

    Parameters:
    -----------
    points: List of points or numpy array of shape (N, 2) with N > 1
        This should be a list of points or a numpy array of shape (N, 2) with N > 1

    Returns:
    --------
        Tuple of floats (x_min, x_max, y_min, y_max)
    """
    if len(points) <= 1:
        raise ValueError(f"The input 'points' must contain more than one element. Found {len(points)} elements.")
    
    x_coords: np.ndarray = np.array([p[0] for p in points])
    y_coords: np.ndarray = np.array([p[1] for p in points])
    return np.min(x_coords), np.max(x_coords), np.min(y_coords), np.max(y_coords)

def truncate_at_intersection(to_truncate: Union['Line', np.ndarray], other: Union['Line', np.ndarray]) -> Union['Line', None]:
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
    for j, line1 in enumerate(iterable=zip(reversed(to_truncate.array[:-1]), reversed(to_truncate.array[1:]))):
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
                return Line(to_truncate.points[:i] + [Point((xc + sol[1] * (xd - xc), yc + sol[1] * (yd - yc)))])
    return None

def truncate_last_segment(last_segment: Union['Line', np.ndarray], other: Union['Line', np.ndarray]) -> Union['Line', None]:
    """
    Compute the truncation of this line segment with another line.
    
    Given a curve C parameterized with s, the truncation is defined as
    all points in C up to the intersection at parameterization value s'

    Parameters:
    ----------
    other: Another line

    Returns:
    --------
        A Line object or None if no intersection
    """
    if len(last_segment) > 2:
        last_segment = Line(last_segment[-2:])

    (xa, ya), (xb, yb) = last_segment

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
            return Line([(xa, ya), Point((xc + sol[1] * (xd - xc), yc + sol[1] * (yd - yc)))])
    return None
