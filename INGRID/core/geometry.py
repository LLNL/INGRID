import numpy as np
import numpy.lib.mixins
from collections.abc import Iterable
from numbers import Number

class Point(np.lib.mixins.NDArrayOperatorsMixin):
    """
    Define a Point.
    Can be used to later define Line objects.

    Inherits from `numpy.lib.mixins.NDArrayOperatorsMixin` to allow
    numpy operations directly on a Point object

    Parameters
    ----------
    *coordinates : array-like, positional
        Accepts either an array like container of numerical values or
        numerical arguments as positional arguments
    
    """

    def __init__(self, *coordinates, **kargs) -> None:

        if len(coordinates) == 1:
            error_message  = 'The provided argument must be an iterable '
            error_message += 'defining Point coordinates.\n'
            error_message += f'Recieved type = {type(coordinates[0])}'
            assert isinstance(coordinates[0], Iterable), error_message
            self.coordinates = coordinates[0]
        else:
            self.coordinates = [coordinate for coordinate in coordinates]

        self.coordinates = self.__array__(dtype=kargs.get('dtype', None))

    def __getitem__(self, slice, **kargs):
        return self.__array__(kargs.get('dtype', None))[slice]
        
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}{(*self.coordinates, )}'

    def __str__(self) -> str:
        return self.__repr__()
    
    def __array__(self, dtype=None):
        return np.array(self.coordinates, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            scalars = []
            for input in inputs:
                if isinstance(input, Number):
                    scalars.append(input)
                elif isinstance(input, self.__class__):
                    scalars.append(input.coordinates)
                elif isinstance(input, np.ndarray):
                    scalars.append(input)
                else:
                    raise NotImplementedError
            return self.__class__(ufunc(*scalars, **kwargs))
        else:
            return NotImplemented
