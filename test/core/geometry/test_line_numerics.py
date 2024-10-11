import numpy as np
from ingrid.core.geometry.point import Point
from ingrid.core.geometry.line import Line

np.random.seed(0)

def test_initialization_000():
    """
    Test the ability to initialize a Line via a numpy array
    """

    initializer = np.random.rand(100, 2)
    L           = Line(points=initializer)
    assert np.array_equal(L.array, initializer)

def test_initialization_001():
    """
    Test the ability to initialize a Line via a list of
    Point objects
    """

    initializer = np.random.rand(100, 2)
    points      = [Point(p) for p in initializer]
    L           = Line(points=points)

    assert np.array_equal(L.array, initializer)
    assert all(lp == p for lp, p in zip(L.points, points))

def test_iterator_000():
    """
    Test our ability to iterate over a Line object
    """

    initializer = np.random.rand(100, 2)
    points      = [Point(p) for p in initializer]
    L           = Line(points=points)
    assert all(lp == p for lp, p in zip(L, points))

def test_line_add_000():
    """
    Test ability to add Line objects and instantiate a new
    copy consisting of two concatenated curves
    """

    initializer = np.random.rand(100, 2)
    slice_1     = slice(0, 50)
    slice_2     = slice(50, 100)
    L1          = Line(points=initializer[slice_1, :])
    L2          = Line(points=initializer[slice_2, :])
    L3          = L1 + L2
    assert np.array_equal(L3.array, np.vstack([L1.array, L2.array]))
    assert L3.points == (L1.points + L2.points)
    assert L3.points is not (L1.points + L2.points)

def test_line_calculate_length_000():
    """
    Test the ability to compute the arclength of an arbitrary Line
    object.

    In particular, computing the arglength of the unit-circle
    """
    baseline_length = 2 * np.pi
    period          = np.linspace(start=0., stop=2 * np.pi, num=1000)
    xy              = np.column_stack([np.cos(period), np.sin(period)])
    L               = Line(points=xy)
    assert np.allclose(L.calculate_length(), baseline_length)

def test_line_calculate_length_001():
    """
    Test the ability to compute the arclength of an arbitrary Line
    object.

    In particular, computing the arclength of sin(x)
    x in [-pi, 2 * pi]
    """
    baseline_length = 11.460593367083
    x               = np.linspace(start=-np.pi, stop=2 * np.pi, num=1000)
    y               = np.sin(x)
    xy              = np.column_stack([x, y])
    L               = Line(points=xy)
    assert np.allclose(L.calculate_length(), baseline_length)


def test_line_calculate_cumulative_length_000():
    """
    Test the ability to compute the cumulative length record of an arbitrary 
    Line object.

    In particular, obtain the cumsum of along a segmented unit-circle
    """
    #
    # Configure a curve representing the unit-circle
    #
    period         = np.linspace(start=0, stop=2 * np.pi, num=10000)
    xy             = np.column_stack([np.cos(period), np.sin(period)])
    L              = Line(points=xy)

    #
    # Compute the differential line segment length
    #
    baseline_differential_length = np.diff(period)[0]

    #
    # Generate the cumulative length via ds * arange trick
    #
    # NOTE: We start at 1 to omit a cumsum value of 0.0
    #
    segments = np.arange(start=1, stop=xy.shape[0])
    baseline_cumulative_length = baseline_differential_length * segments

    result = L.calculate_cumulative_length()
    assert np.allclose(result, baseline_cumulative_length)

def test_line_point_at_parametric_000():
    """
    Test the ability to recover point coordinates by treating
    the Line object as a parametric curve

    In particular, test that we can recover points along a simple
    line eminating from the origin
    """

    x  = np.linspace(start=0., stop=2., num=10000)
    y  = x
    xy = np.column_stack([x, y])
    L  = Line(points=xy)

    parametric_t_vals = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, -1.]
    for pt in parametric_t_vals:
        if pt > 1.:
            target = Point(xy[-1])
        elif pt < 0.:
            target = Point(xy[0])
        else:
            target = Point(np.array([2. * pt, 2. * pt]))
        result = L.at(parametric_t=pt)
        assert np.allclose(result, target, rtol=0.001)

def test_line_truncation_perpendicular():
    num_segments = 10000
    # Horizontal line segment at y = 0 from x = [0, 1]
    line_A: Line = Line(np.column_stack([np.linspace(0, 1, num=num_segments+1), np.zeros(num_segments + 1)]))
    # Vertical line segment at x = 0.5 from y = [-1, 1] 
    line_B: Line = Line(np.column_stack([0.5 * np.ones(num_segments + 1), np.linspace(-1, 1, num_segments + 1)]))
    # Should be line segment at y = 0 from x = [0, 0.5] with half the points of line_A
    line_C: Line = line_A.truncate_at_intersection(other=line_B)
    subset = line_A[:line_A.shape[0] // 2]
    assert line_C.shape == subset.shape
    print(np.linalg.norm(line_C - subset))
