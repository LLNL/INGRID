from typing import Any
import numpy as np
from INGRID.core.geometry import Point

np.random.seed(0)

def test_add_000():
    """
    Test adding Point objects: Point + Point
    """

    initializer: np.ndarray = np.random.rand(2, 2)
    A: Point = Point(initializer[0])
    B: Point = Point(initializer[1])
    result: Point = A + B
    baseline: Any = initializer.sum(axis=0)

    assert np.array_equal(result, baseline)
    assert isinstance(result, Point)

def test_add_001():
    """
    Test adding Point objects: Point + np.array
    """

    initializers: np.ndarray = np.random.rand(2, 2)
    A: Point = Point(initializers[0])
    result: Point = A + initializers[1]
    baseline: Any = initializers.sum(axis=0)

    assert np.array_equal(result, baseline)
    assert isinstance(result, Point)

def test_add_002():
    """
    Test adding Point objects: Point + scalars

    In particular, test addition of int and float arguments

    """

    A: Point = Point(np.random.rand(2))

    assert np.array_equal(A, A + 0)
    assert isinstance(A + 0, Point)

    assert np.array_equal(A, A + 0.)
    assert isinstance(A + 0., Point)

def test_scale_000():
    """
    Test scaling a Point object by scalers

    In particular, test scaling with int and float arguments

    """

    initializer: np.ndarray = np.random.rand(2).astype(float)
    A: Point = 2 * Point(initializer)
    int_baseline: np.ndarray = 2 * initializer

    initializer: np.ndarray = np.random.rand(2).astype(int)
    B: Point = 2. * Point(initializer)
    float_baseline: np.ndarray = 2. * initializer

    assert np.array_equal(A, int_baseline)
    assert np.array_equal(B, float_baseline)

def test_slicing_000():
    """
    Test the ability to slice into a Point object like a numpy array
    """

    initializer: np.ndarray = np.random.rand(100).astype(int)
    indices: np.ndarray = np.random.randint(0, 100, size=50)
    indices.sort()
    indices: np.ndarray = np.unique(indices)
    A: Point = Point(initializer)
    assert np.array_equal(A[indices], initializer[indices])

def test_initialization_000():
    """
    Test the ability to initialize a Point via an array-like container
    """

    initializer: np.ndarray= np.random.rand(100)
    A: Point = Point(initializer)

def test_initialization_001():
    """
    Test the ability to initialize a Point via an array-like container
    and recover data via the 'coordinates' attribute
    """

    initializer: np.ndarray= np.random.rand(100)
    A: Point = Point(initializer)

    assert np.array_equal(A, A.coordinates)

def test_initialization_002():
    """
    Test the ability to initialize a Point via an array-like container
    and a keyword argument
    """

    initializer: np.ndarray = np.random.rand(100)
    A: Point = Point(initializer)

def test_distance_to_000():
    """
    Test the ability to compute the distance between two Point objects
    """

    initializer: np.ndarray = np.random.rand(2, 2)
    A: Point = Point(initializer[0])
    B: Point = Point(initializer[1])
    baseline_distance: np.floating[Any] = np.linalg.norm(np.diff(initializer, axis=0))
    result: np.floating[Any] = A.distance_to(B)

    print(result, baseline_distance)
    assert np.isclose(result, baseline_distance)
    assert isinstance(result, float)

def test_distance_to_001():
    """
    Test the ability to compute the distance between a Point object and
    a numpy array
    """

    initializer: np.ndarray = np.random.rand(2, 2)
    A: Point = Point(initializer[0])
    B: np.float = initializer[1]
    baseline_distance: np.floating[Any] = np.linalg.norm(np.diff(initializer, axis=0))
    result: np.floating[Any] = A.distance_to(B)

    print(result, baseline_distance)
    assert np.isclose(result, baseline_distance)
    assert isinstance(result, float)
