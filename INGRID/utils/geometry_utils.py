import numpy as np
import ingrid.core.geometry as geometry

def find_split_index(line: geometry.Line, 
                     target_point: geometry.Point,
                     eps: float = 1e-9):
    """
    Attempt to find the index of a Point in a given Line
    object which is closest to an input target Point

    Parameters:
    -----------
    line: geometry.Line
        Line object to search for a split index
    target_point: geometry.Point
        Point object use for searching for an index
    eps: float
        Epsilon value to use for colinearity calculation

    Returns:
    --------
        The index of the nearest Point object, or None
        if not found
    """

    found_split_index = False

    for index, point in enumerate(line[:-2]):

        if np.allclose(point, target_point):
            found_split_index = True
            break

        vector_to_next_point   = line[index + 1] - point
        vector_to_target_point = target_point - point

        cross_vector = np.cross(vector_to_next_point, vector_to_target_point)

        colinear_criteria = [
            #
            # Vectors are colinear up to epsilon
            #
            np.linalg.norm(cross_vector) < eps,
            #
            # Angle between vectors are constrained to quadrants I & IV
            #
            np.dot(vector_to_next_point, vector_to_target_point) > 0,
            #
            # Vector to next point has greater magnitude than magnitude
            # of vector to target point
            #
            np.linalg.norm(vector_to_next_point) \
                > np.linalg.norm(vector_to_target_point)
        ]
        
        if all(colinear_criteria):
            found_split_index = True
            break

    index = index if found_split_index else None
    return index

def split_line_at_index(line: geometry.Line, 
                        index: int):

    n_points = line.array.shape[0]

    if index == n_points - 1:
        segment_A = line
        segment_B = None
    elif index == n_points - 2:
        segment_A = geometry.Line(points=line[:-2])
        segment_B = geometry.Line(points=line[-2:-1])
    else:
        segment_A_coordinates = line[:index + 1]
        segment_B_coordinates = line[index + 2:]
        segment_A = geometry.Line(points=segment_A_coordinates)
        segment_B = geometry.Line(points=segment_B_coordinates)

    return segment_A, segment_B

def split_line_at_point(line: geometry.Line, target_point: geometry.Point):
    """
    Split a Line object at a

    Parameters:
    -----------
    line: geometry.Line
        Line to operate on
    parametric_t: float
        Parametric t value to split a Line at

    Returns:
    --------
        A tuple of Line objects whose union defines an equivalent curve.
        The generated parametric_t point will be included in the first segment
    """ 

    index = find_split_index(line=line, target_point=target_point)

    if index is None:
        msg  = f"Attempted to split Line object {line} "
        msg += f"with target Point {target_point}. Could not find an index "
        msg += "that satisfies the colinearity criteria"
        raise ValueError(msg)

    segment_A, segment_B = split_line_at_index(line=line, index=index)

    return segment_A, segment_B

def split_line_at_parametric_t(line: geometry.Line, parametric_t: float):
    """
    Split a Line object at a parametric_t value.

    Parameters:
    -----------
    line: geometry.Line
        Line to operate on
    parametric_t: float
        Parametric t value to split a Line at

    Returns:
    --------
        A tuple of Line objects whose union defines an equivalent curve.
        The generated parametric_t point will be included in the first segment
    """
    generated_point      = line.at(parametric_t=parametric_t)
    segment_A, segment_B = split_line_at_point(
                                line  = line, 
                                point = generated_point
                            )
    
    return segment_A, segment_B
