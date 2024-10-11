import numpy as np
from typing import List, Optional
import ingrid.core.geometry.utils as geo_utils
from ingrid.core.geometry.line import Line, Point

np.random.seed(seed=0)

def test_line_truncation_perpendicular() -> None:
    """
    Test the truncation of a line segment with a perpendicular line segment.
    """
    num_segments: int = 10

    #
    # Horizontal line segment at y = 0 from x = [0, 1]
    #
    line_A: Line = Line(np.column_stack((np.linspace(0, 1, num=num_segments+1), np.zeros(num_segments + 1))))
    
    #
    # Vertical line segment at x = 0.5 from y = [-1, 1]
    #
    line_B: Line = Line(np.column_stack((0.5 * np.ones(num_segments + 1), np.linspace(-1, 1, num_segments + 1))))
    
    #
    # Should be line segment at y = 0 from x = [0, 0.5] with half the points of line_A
    #
    line_C: Optional[Line] = geo_utils.truncate_at_intersection(to_truncate=line_A, other=line_B)
    assert line_C is not None
    subset: np.ndarray = line_A[:line_A.shape[0] // 2]
    assert line_C.shape == subset.shape

def test_line_truncation_last_segment() -> None:
    """
    Test the truncation of a line segment at its last segment.
    """
    num_segments: int = 100
    
    #
    # Horizontal line segment at y = 0 from x = [0, 1]
    #
    line_A: Line = Line(np.column_stack((np.linspace(0, 1, num=num_segments+1), np.zeros(num_segments + 1))))
    
    #
    # Get the last segment  
    #
    last_segment: np.ndarray = line_A[-2:]
    
    #
    # Get the middle x value of the last segment
    #
    mid_x_val: float = last_segment[0, 0] + 0.5 * (last_segment[1, 0] - last_segment[0, 0])
    line_B: Line = Line(np.column_stack((np.full(shape=num_segments + 1, fill_value=mid_x_val), np.linspace(-1, 1, num_segments + 1))))
    
    #
    # Should be line segment at y = 0 from x = [0, mid_x_val]
    #
    line_C: Optional[Line] = geo_utils.truncate_last_segment(last_segment=line_A, other=line_B)
    
    #
    # Check that the truncation is not None
    #
    assert line_C is not None
    
    #
    # Check that the last point of the truncation is the middle x value of the last segment
    #
    assert line_C[-1, 0] == mid_x_val

def test_bounding_box_with_two_points() -> None:
    """
    Test the bounding box computation with two points.
    """
    points: List[Point] = [Point(0, 0), Point(1, 1)]
    
    #
    # Should be x_min = 0, x_max = 1, y_min = 0, y_max = 1 since its the unit square
    #
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    x_min, x_max, y_min, y_max = geo_utils.compute_bounding_box_from_points(points=points)
    assert x_min == 0 and x_max == 1 and y_min == 0 and y_max == 1