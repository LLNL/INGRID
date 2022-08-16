import numpy as np
from ingrid.core.geometry import Point

np.random.seed(0)

def test_add_000():
    """
    Test adding Point objects: Point + Point
    """

    initializers = np.random.rand(2, 2)
    A, B         = Point(initializers[0]), Point(initializers[1])
    result       = A + B
    baseline     = initializers.sum(axis=0)

    assert np.array_equal(result, baseline)
    assert isinstance(result, Point)

def test_add_001():
    """
    Test adding Point objects: Point + np.array
    """

    initializers = np.random.rand(2, 2)
    A            = Point(initializers[0])
    result       = A + initializers[1]
    baseline     = initializers.sum(axis=0)

    assert np.array_equal(result, baseline)
    assert isinstance(result, Point)

def test_add_002():
    """
    Test adding Point objects: Point + scalars

    In particular, test addition of int and float arguments

    """

    A = Point(np.random.rand(2))

    assert np.array_equal(A, A + 0)
    assert isinstance(A + 0, Point)

    assert np.array_equal(A, A + 0.)
    assert isinstance(A + 0., Point)

def test_scale_000():
    """
    Test scaling a Point object by scalers

    In particular, test scaling with int and float arguments

    """

    initializer    = np.random.rand(2).astype(float)
    A              = 2 * Point(initializer)
    int_baseline   = 2 * initializer

    initializer    = np.random.rand(2).astype(int)
    B              = 2. * Point(initializer)
    float_baseline = 2. * initializer

    assert np.array_equal(A, int_baseline)
    assert np.array_equal(B, float_baseline)

def test_slicing_000():
    """
    Test the ability to slice into a Point object like a numpy array
    """

    initializer = np.random.rand(100).astype(int)
    indices     = np.random.randint(0, 100, size=50)
    indices.sort()
    indices     = np.unique(indices)
    A           = Point(initializer)
    assert np.array_equal(A[indices], initializer[indices])

def test_initialization_000():
    """
    Test the ability to initialize a Point via an array-like container
    or via positional-arguments/unpacking
    """

    initializer = np.random.rand(100)

    A = Point(initializer)
    B = Point(*initializer)

    assert np.array_equal(A, B)

def test_initialization_001():
    """
    Test the ability to initialize a Point via an array-like container
    and recover data via the 'coordinates' attribute
    """

    initializer = np.random.rand(100)
    A           = Point(initializer)

    assert np.array_equal(A, A.coordinates)
