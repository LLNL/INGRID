from ast import Num
import numpy as np
import numpy.lib.mixins
from collections.abc import Iterable
from numbers import Number

class Point(np.lib.mixins.NDArrayOperatorsMixin):

    def __init__(self, *coordinates, **kargs) -> None:

        if len(coordinates) == 1:
            error_message  = 'The provided argument must be an iterable '
            error_message += 'defining Point coordinates.\n'
            error_message += f'Recieved type = {type(coordinates[0])}'
            assert isinstance(coordinates[0], Iterable), error_message
            self.coordinates = coordinates[0]
        else:
            self.coordinates = [coordinate for coordinate in coordinates]
        
    def __repr__(self) -> str:
        return self.__repr__()

    def __str__(self) -> str:
        return f'{self.__class__.__name__}{(*self.coordinates, )}'
    
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

if __name__ == '__main__':

    r = np.random.rand(10, 2)

    for xy in r:
        print(Point(xy))
        result =  2 * (Point(xy) + [2, 45]) + (xy + [1 ,2 ])
        print(result, type(result))
        print(result + result)
