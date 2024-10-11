from ingrid.core.geometry.line import Line
from ingrid.core.geometry.boundaries import ClosedBoundary
from ingrid.core.geometry.exceptions import MisalignedBoundaryError
import numpy as np
import pytest

def test_ClosedBoundary_with_closed_circle() -> None:
    """
    Test ClosedBoundary accept a closed circle
    """
    north_t: np.ndarray = np.linspace(start=0, stop=0.5 * np.pi)
    east_t: np.ndarray = np.linspace(start=np.pi * 0.5, stop=np.pi)
    south_t: np.ndarray = np.linspace(start=np.pi, stop=1.5 * np.pi)
    west_t: np.ndarray = np.linspace(start=1.5 * np.pi, stop=2.0 * np.pi)

    north: Line = Line(np.column_stack((np.cos(north_t), np.sin(north_t))))
    east: Line = Line(np.column_stack((np.cos(east_t), np.sin(east_t))))
    south: Line = Line(np.column_stack((np.cos(south_t), np.sin(south_t))))
    west: Line = Line(np.column_stack((np.cos(west_t), np.sin(west_t))))

    _: ClosedBoundary = ClosedBoundary(north=north, east=east, south=south, west=west)

def test_ClosedBoundary_with_gapped_circle() -> None:
    """
    Test ClosedBoundary rejects a gapped circle
    """
    north_t: np.ndarray = np.linspace(start=0, stop=0.5 * np.pi)[:-2]
    east_t: np.ndarray = np.linspace(start=np.pi * 0.5, stop=np.pi)
    south_t: np.ndarray = np.linspace(start=np.pi, stop=1.5 * np.pi)
    west_t: np.ndarray = np.linspace(start=1.5 * np.pi, stop=2.0 * np.pi)

    north: Line = Line(np.column_stack((np.cos(north_t), np.sin(north_t))))
    east: Line = Line(np.column_stack((np.cos(east_t), np.sin(east_t))))
    south: Line = Line(np.column_stack((np.cos(south_t), np.sin(south_t))))
    west: Line = Line(np.column_stack((np.cos(west_t), np.sin(west_t))))

    with pytest.raises(expected_exception=MisalignedBoundaryError):
        _: ClosedBoundary = ClosedBoundary(north=north, east=east, south=south, west=west)

def test_ClosedBoundary_has_unique_points() -> None:
    """
    Test ClosedBoundary applies point uniqueness
    """
    north_t: np.ndarray = np.linspace(start=0, stop=0.5 * np.pi)
    east_t: np.ndarray = np.linspace(start=np.pi * 0.5, stop=np.pi)
    south_t: np.ndarray = np.linspace(start=np.pi, stop=1.5 * np.pi)
    west_t: np.ndarray = np.linspace(start=1.5 * np.pi, stop=2.0 * np.pi)

    north: Line = Line(np.column_stack((np.cos(north_t), np.sin(north_t))))
    east: Line = Line(np.column_stack((np.cos(east_t), np.sin(east_t))))
    south: Line = Line(np.column_stack((np.cos(south_t), np.sin(south_t))))
    west: Line = Line(np.column_stack((np.cos(west_t), np.sin(west_t))))

    # Create a line doubling up on all points
    double_north: np.ndarray = np.array([val for point in zip(north, north) for val in point])
    result: ClosedBoundary = ClosedBoundary(north=double_north, east=east, south=south, west=west)
    assert result.north.shape[0] == double_north.shape[0] / 2
    assert np.allclose(result.north, north)