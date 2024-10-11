import numpy as np
from ingrid.core.geometry.point import Point
from ingrid.core.geometry.line import Line

np.random.seed(0)

def test_initialization_with_numpy_array() -> None:
    """
    Test the ability to initialize a Line via a numpy array
    """

    initializer: np.ndarray = np.random.rand(100, 2)
    L: Line = Line(points=initializer)
    assert np.array_equal(L.array, initializer)

def test_initialization_with_point_list() -> None:
    """
    Test the ability to initialize a Line via a list of
    Point objects
    """

    initializer: np.ndarray = np.random.rand(100, 2)
    points: list[Point] = [Point(p) for p in initializer]
    L: Line = Line(points=points)

    assert np.array_equal(L.array, initializer)
    assert all(lp == p for lp, p in zip(L.points, points))

def test_iterator_over_line_object() -> None:
    """
    Test our ability to iterate over a Line object
    """

    initializer: np.ndarray = np.random.rand(100, 2)
    points: list[Point] = [Point(p) for p in initializer]
    L: Line = Line(points=points)
    assert all(lp == p for lp, p in zip(L, points))

def test_line_addition() -> None:
    """
    Test ability to add Line objects and instantiate a new
    copy consisting of two concatenated curves
    """

    initializer: np.ndarray = np.random.rand(100, 2)
    slice_1: slice = slice(0, 50)
    slice_2: slice = slice(50, 100)
    L1: Line = Line(points=initializer[slice_1, :])
    L2: Line = Line(points=initializer[slice_2, :])
    L3: Line = L1 + L2
    assert np.array_equal(L3.array, np.vstack([L1.array, L2.array]))
    assert L3.points == (L1.points + L2.points)
    assert L3.points is not (L1.points + L2.points)

def test_line_calculate_length_unit_circle() -> None:
    """
    Test the ability to compute the arclength of an arbitrary Line
    object.

    In particular, computing the arglength of the unit-circle
    """
    baseline_length: float = 2 * np.pi
    period: np.ndarray = np.linspace(start=0., stop=2 * np.pi, num=1000)
    xy: np.ndarray = np.column_stack([np.cos(period), np.sin(period)])
    L: Line = Line(points=xy)
    assert np.allclose(L.calculate_length(), baseline_length)

def test_line_calculate_length_sin_curve() -> None:
    """
    Test the ability to compute the arclength of an arbitrary Line
    object.

    In particular, computing the arclength of sin(x)
    x in [-pi, 2 * pi]
    """
    baseline_length: float = 11.460593367083
    x: np.ndarray = np.linspace(start=-np.pi, stop=2 * np.pi, num=1000)
    y: np.ndarray = np.sin(x)
    xy: np.ndarray = np.column_stack([x, y])
    L: Line = Line(points=xy)
    assert np.allclose(L.calculate_length(), baseline_length)


def test_line_calculate_cumulative_length_unit_circle() -> None:
    """
    Test the ability to compute the cumulative length record of an arbitrary 
    Line object.

    In particular, obtain the cumsum of along a segmented unit-circle
    """
    #
    # Configure a curve representing the unit-circle
    #
    period: np.ndarray = np.linspace(start=0, stop=2 * np.pi, num=10000)
    xy: np.ndarray = np.column_stack([np.cos(period), np.sin(period)])
    L: Line = Line(points=xy)

    #
    # Compute the differential line segment length
    #
    baseline_differential_length: float = np.diff(period)[0]

    #
    # Generate the cumulative length via ds * arange trick
    #
    # NOTE: We start at 1 to omit a cumsum value of 0.0
    #
    segments: np.ndarray = np.arange(start=1, stop=xy.shape[0])
    baseline_cumulative_length: np.ndarray = baseline_differential_length * segments

    result: np.ndarray = L.calculate_cumulative_length()
    assert np.allclose(result, baseline_cumulative_length)

def test_line_point_at_parametric() -> None:
    """
    Test the ability to recover point coordinates by treating
    the Line object as a parametric curve

    In particular, test that we can recover points along a simple
    line eminating from the origin
    """

    x: np.ndarray = np.linspace(start=0., stop=2., num=100)
    y: np.ndarray = x
    xy: np.ndarray = np.column_stack([x, y])
    L: Line = Line(points=xy)

    parametric_t_vals: list[float] = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, -1.]
    for pt in parametric_t_vals:
        if pt > 1.:
            target: Point = Point(xy[-1])
        elif pt < 0.:
            target: Point = Point(xy[0])
        else:
            target: Point = Point(np.array([2. * pt, 2. * pt]))
        result: Point = L.at(parametric_t=pt)
        assert np.allclose(result, target, rtol=0.001)
