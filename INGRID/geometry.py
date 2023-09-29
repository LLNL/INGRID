import numpy as np
import numpy.lib.mixins
from collections.abc import Iterable
from numbers import Number

class Point(np.lib.mixins.NDArrayOperatorsMixin):
    """
    Define a Point in the RZ plane.
    Can be used to later define Line objects.

    Inherits from `numpy.lib.mixins.NDArrayOperatorsMixin` to allow
    numpy operations directly on a Point object

    Parameters
    ----------
    rz : array-like
        Accepts either an array like container of numerical values or
        numerical arguments 
    
    """
    def __init__(self, rz, **kargs) -> None:
        self.dtype   = kargs.get('dtype', float)
        self.data    = rz
        self.data    = self.__array__(dtype=self.dtype)
        self.shape   = self.data.shape
        shape_check  = (2,)
        format_error = f'Point data must result in shape {shape_check}'
        assert self.shape == shape_check, format_error

    def __getitem__(self, slice, **kargs):
        return self.data[slice]
        
    def __repr__(self) -> str:
        return f'{self.__class__.__name__}{(*self.data, )}'

    def __str__(self) -> str:
        return self.__repr__()
    
    def __array__(self, dtype=None):
        return np.array(self.data, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            scalars = []
            for input in inputs:
                if isinstance(input, Number):
                    scalars.append(input)
                elif isinstance(input, self.__class__):
                    scalars.append(input.data)
                elif isinstance(input, np.ndarray):
                    scalars.append(input)
                else:
                    raise NotImplementedError
            return self.__class__(ufunc(*scalars, **kwargs))
        else:
            return NotImplemented

    def distance_to(self, point):
        """
        Compute the distance in unitless coordinates to another Point

        Parameters:
        -----------
        point: Point, np.ndarray
            Another Point object or numpy array of same shape
        
        Returns:
        --------
            Euclidean distance to the other point
        """
        return np.linalg.norm(self - point)
    
    def plot(self, ax: 'matplotlib.axes.Axes', **kargs) -> None:
        """
        Plot the Point.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The Axes instance to plot to.
        """
        ax.plot(*self.data, **kargs)

    def psi(self, grid: 'EfitData', tag: str = 'v') -> float:
        """
        Get the psi value of this Point from an EfitData instance.

        Parameters
        ----------
        grid : EfitData
            The grid upon which the value of psi is to be calculated on.
        tag : str, optional
            Char to specify the type of psi derivative.
            Defaults 'v' (no derivative).

        Returns
        -------
            The psi value at the Point.
        """
        return grid.PsiNorm.get_psi(self.x, self.y, tag)

class Line:
    """
    Define a Line.
    Can be used to later define Patch objects.

    Parameters
    ----------
    points : array-like
        An array like container containing either numerical values of shape (n, 2)
        An array like container containing Point objects
    """

    def __init__(self, points, **kargs) -> None:

        # Case container of Point objects
        if all(isinstance(point, Point) for point in points):
            data        = np.array(points, dtype=float)
        # Case numpy array input
        elif isinstance(points, np.ndarray):
            data        = points
        # Case list. Convert to numpy array and assume format is valid
        elif isinstance(points, list):
            data = np.array(points, dtype=float)
        else:
            raise ValueError('Invalid "points" argument type')
        #
        # Enforce shape (n, 2)
        #
        format_error  = f'Invalid format of points provided. Final shape = {data.shape}. \n'
        format_error += f'Points argument = {points}'
        assert data.shape[-1] == 2, format_error

        # If we didnt set self.points, do it now with validated data

        self.points = [Point(point) for point in data]
        
        # Finalize setup and compute quantities
        self.data                   = data
        self.arclen                 = self.calculate_length()
        self.cumlen_data            = self.calculate_cumulative_length()
        self.normalized_cumlen_data = self.cumlen_data / self.arclen

    def __getitem__(self, val, **kargs):
        """
        Define slicing/indexing capabilities
        """
        if not isinstance(val, (slice, Iterable, int)):
            raise NotImplementedError
        return self.points[val]

    def __iter__(self):
        """
        Dispatch to the list of points
        """
        return self.points.__iter__()

    def __next__(self):
        """
        Dispatch to the list of points
        """
        return self.points.__next__()
    
    def __len__(self):
        return len(self.data)
        
    def __repr__(self) -> str:
        _repr  = f'{self.__class__.__name__} (\n'
        _repr += '\n'.join('    ' + repr(p) for p in self.points)
        _repr += '\n)'
        return _repr

    def __str__(self) -> str:
        return self.__repr__()
    
    def __array__(self, dtype=float):
        return np.array(self.data, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method == '__call__':
            scalars = []
            for input in inputs:
                if isinstance(input, Number):
                    scalars.append(input)
                elif isinstance(input, self.__class__):
                    scalars.append(input.data)
                elif isinstance(input, np.ndarray):
                    scalars.append(input)
                else:
                    raise NotImplementedError
            return self.__class__(ufunc(*scalars, **kwargs))
        else:
            return NotImplementedError

    def __add__(self, other):
        """
        Concatenate two line objects
        """
        if not isinstance(other, self.__class__):
            error_message  = 'Line addition can only be performed between '
            error_message += 'Line types. Recieved type = "{type(other)}"'
            raise ValueError(error_message)
        
        # Concatenate values along first dimension
        new_points = np.vstack([self.points, other.points])
        return Line(points=new_points)

    def calculate_length(self):
        """
        Calculate the arclength of the Line object

        Returns:
        --------
            Float computing the sum of all curve segments
        """
        return np.linalg.norm(np.diff(self.data, axis=0), axis=1).sum()

    def calculate_cumulative_length(self):
        """
        Calculate the cumulative arclength of the Line object

        Returns:
        --------
            Numpy array containing cumulative lengths
        """
        return np.linalg.norm(np.diff(self.data, axis=0), axis=1).cumsum()
    
    def at(self, parametric_t):
        """
        Interface to compute_parametric_point
        """
        return self.compute_parametric_point(parametric_t=parametric_t)

    def compute_parametric_point(self, parametric_t):
        """
        Infer the Point which would reside along a parametric representation
        of the Line object.

        Given a parametric value p in [0., 1.], compute the coordinates via
        linear interpolation along line segments

        Parameters:
        -----------
        parametric_t: float
            A numeric representing the input to the parametric curve equation
        
        Returns:
        --------
            A Point object
        """

        parametric_t = np.clip(parametric_t, a_min=0., a_max=1.)

        if np.allclose(parametric_t, 1.):
            result = self.points[-1]
        elif np.allclose(parametric_t, 0.):
            result = self.points[0]
        else:
            #
            # Find indices where parametric request is less than normalized
            # cumsum of line segments. 
            #
            # This allows us to find the Point in the Line *after* where our 
            # parametric request should lie
            #
            where_below = \
                np.where(parametric_t < self.normalized_cumlen_data)[0]
            endpoint_index = where_below[0]
            endpoint = self.points[endpoint_index]
            
            endpoint_parametric_t = \
                self.normalized_cumlen_data[endpoint_index]

            parametric_t_diff = endpoint_parametric_t - parametric_t
            
            #
            # Compute new coordinates via linear interpolation
            #
            result = endpoint * (1 - parametric_t_diff)

        return result
    
    def refined_copy(self, num_points_between: int):
        """
        Refine the Line object with linear interpolation given a requested
        number of points

        Parameters
        ----------
        num_points_between : int
            The number of points place between each line segment of the
            refined copy
        
        Returns
        -------
            A Line object
        """
        refined_array = []
        for i in range(self.data.shape[0] - 1):
            refined_array.append(self.data[i])
            #
            # Perform the linear interpolation over each dimension
            #
            interpolated_values = np.column_stack(
                [
                    np.linspace(self.data[i][dim], self.data[i + 1][dim], num_points_between + 2)[1:-1]
                    for dim in range(self.data.shape[-1])
                ]
            )
            refined_array.extend(interpolated_values)
        refined_array.append(self.data[-1])
        return Line(points=np.array(refined_array))

    def plot(self, ax, **kargs):
        """
        Plot the Line

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The Axes instance to plot to.
        """
        ax.plot(*self.data.T, **kargs)

class Patch:
    """
    Define a Patch representing a portion of the tokamak domain.

    Each Patch can be refined into a subgrid that can then create a global grid.

    Parameters
    ----------
    lines : array-like
        The four Line objects defining the boundary of this Patch (N, E, S, W).

    patch_name : 

    """

    def __init__(self, lines: 'array-like', patch_name: str = '', PatchTagMap: dict = None,
                 plate_patch: bool = False, plate_location: str = None, color: str = 'blue'):
        self.lines = lines
        self.N = lines[0]
        self.E = lines[1]
        self.S = lines[2]
        self.W = lines[3]
        self.BoundaryPoints = {}

        # This is the border for the fill function
        # It need to only include N and S lines
        self.RemoveDuplicatePoints()
        self.p = list(self.N.p) + list(self.E.p) + list(self.S.p) + list(self.W.p)
        self.PatchTagMap = PatchTagMap
        self.plate_patch = plate_patch
        self.plate_location = plate_location
        self.patch_name = patch_name
        self.color = color
        self.Verbose = True

        # TODO: Populate these attributes when generating subgrid.
        #       Currently this is being done in the concat_grid method.
        self.npol = 0
        self.nrad = 0

        self.PatchLabelDoc = {
            'I': 'Inner',
            'O': ' Outer',
            'D': 'Divertor',
            'L': 'Leg',
            'P': 'Private',
            'F': 'Flux',
            'T': 'Top',
            'B': 'Bottom',
            'S': 'SOL',
            'C': 'Core'
        }

        self.ax = None
        self.cell_grid = None

    def plot_border(self, color='red', ax=None):
        """
        Draw solid borders around the patch.

        Parameters
        ----------
        color : str, optional
            Defaults to red.
        """

        _ax = ax
        for line in self.lines:
            line.plot(color, ax=_ax, linewidth=0.5)

    def fill(self, color='lightsalmon', ax=None, alpha=1.0):
        """
        Shades in the patch with a given color

        Parameters
        ----------
        color : str, optional
            Defaults to a light salmon.
        """
        x = np.array([p.x for p in self.p])
        y = np.array([p.y for p in self.p])
        arr = np.column_stack((x, y))
        patch = Polygon(arr, fill=True, closed=True, color=color, label=self.get_tag(), alpha=alpha)
        _ax = plt.gca() if ax is None else ax
        _ax.add_patch(patch)
        _ax.plot()

    def plot_subgrid(self, fig=None, ax=None, color='blue'):
        if fig is None:
            fig = plt.gcf()
        if ax is None:
            ax = plt.gca()
        plt.figure(fig.number)
        for row in self.cell_grid:
            for cell in row:
                cell.plot_border(color=color, ax=ax)
        plt.pause(0.1)

    def RemoveDuplicatePoints(self):
        for line in [self.N, self.E, self.S, self.W]:
            line.RemoveDuplicatePoints()
        self.p = list(self.N.p) + list(self.E.p) + list(self.S.p) + list(self.W.p)

    def adjust_corner(self, point, corner):
        if corner == 'NW':
            self.cell_grid[0][0].vertices[corner] = point
            self.cell_grid[0][0].vertices[corner] = point
        elif corner == 'NE':
            self.cell_grid[0][-1].vertices[corner] = point
            self.cell_grid[0][-1].vertices[corner] = point
        elif corner == 'SE':
            self.cell_grid[-1][-1].vertices[corner] = point
            self.cell_grid[-1][-1].vertices[corner] = point
        elif corner == 'SW':
            self.cell_grid[-1][0].vertices[corner] = point
            self.cell_grid[-1][0].vertices[corner] = point

    def AdjustBorder(self, face, patch):

        if face == 'E':
            for i in range(len(self.cell_grid)):
                self.cell_grid[i][-1].vertices['NE'] = patch.cell_grid[i][0].vertices['NW']
                self.cell_grid[i][-1].vertices['SE'] = patch.cell_grid[i][0].vertices['SW']
        elif face == 'W':
            for i in range(len(self.cell_grid)):
                self.cell_grid[i][0].vertices['NW'] = patch.cell_grid[i][-1].vertices['NE']
                self.cell_grid[i][0].vertices['SW'] = patch.cell_grid[i][-1].vertices['SE']
        elif face == 'N':
            for j in range(len(self.cell_grid[0])):
                self.cell_grid[0][j].vertices['NW'] = patch.cell_grid[-1][j].vertices['SW']
                self.cell_grid[0][j].vertices['NE'] = patch.cell_grid[-1][j].vertices['SE']
        elif face == 'S':
            for j in range(len(self.cell_grid[0])):
                self.cell_grid[-1][j].vertices['SW'] = patch.cell_grid[0][j].vertices['NW']
                self.cell_grid[-1][j].vertices['SE'] = patch.cell_grid[0][j].vertices['NE']
        else:
            raise ValueError(f"# Invalid face '{face}' provided for adjusting.")

    def get_tag(self):
        return self.PatchTagMap[self.patch_name]

    def get_settings(self):
        settings = {}
        settings['patch_name'] = self.patch_name
        settings['plate_patch'] = self.plate_patch
        settings['plate_location'] = self.plate_location
        settings['PatchTagMap'] = self.PatchTagMap
        return settings

    def cell_grid_as_np(self):
        if self.cell_grid is None:
            cg_np = np.array([])
        else:
            cg_np = []
            for row in self.cell_grid:
                r = []
                for cell in row:
                    r.append(cell.as_np())
                cg_np.append(r)
            cg_np = np.array(cg_np)
        return cg_np

    def as_np(self):

        patch_data = []
        for line in [self.N, self.E, self.S, self.W]:
            R, Z = line.as_np()
            patch_data.append(R)
            patch_data.append(Z)
        patch_data = np.array(patch_data, dtype=object)
        cell_data = self.cell_grid_as_np()
        patch_settings = self.get_settings()
        payload = np.array([patch_data, cell_data, patch_settings], dtype=object)
        return payload

    def make_subgrid(self, grid, np_cells=2, nr_cells=2, _poloidal_f=lambda x: x, _radial_f=lambda x: x, verbose=False, visual=False, ShowVertices=False):
        """
        Generate a refined grid within a patch.
        This 'refined-grid' within a Patch is a collection
        of num x num Cell objects

        Parameters
        ----------
        grid : Ingrid
                To be used for obtaining Efit data and all
                other class information.
        num  : int, optional
                Number to be used to generate num x num
                cells within our Patch.
        """

        def psi_parameterize(grid, r, z):
            """
            Helper function to be used to generate a
            list of values to parameterize a spline
            in Psi. Returns a list to be used for splprep only
            """
            vmax = grid.PsiNorm.get_psi(r[-1], z[-1])
            vmin = grid.PsiNorm.get_psi(r[0], z[0])

            vparameterization = np.empty(0)
            for i in range(len(r)):
                vcurr = grid.PsiNorm.get_psi(r[i], z[i])
                vparameterization = np.append(vparameterization, abs((vcurr - vmin) / (vmax - vmin)))

            return vparameterization

        def psi_test(leg, grid):
            """
            Determine whether a leg's orientation is increasing in psi
            by checking the endpoints of a given leg.

            Returns True if final Point entry is greater in psi than the
            first Point entry. Returns False otherwise.
            """
            p1 = leg.p[0]
            p2 = leg.p[-1]

            return True if grid.PsiNorm.get_psi(p2.x, p2.y) > grid.PsiNorm.get_psi(p1.x, p1.y) else False

        def set_dynamic_step(ratio=0.01):
            """
            Compute a dt value for line_tracing that occurs during subgrid generation of
            this Patch object. The value is taken to be a ratio of the average length of
            a Patch in the radial direction.
            """
            d = 0
            for i in range(nr_lines):
                d += np.sqrt((self.E_vertices[i].x - self.W_vertices[i].x)**2 + (self.E_vertices[i].y - self.W_vertices[i].y)**2)
            dynamic_step = d / nr_lines
            # print('Dynamic-Step value came out to be: {}\n'.format(dynamic_step * ratio))
            return dynamic_step * ratio

        # Allocate space for collection of cell objects.
        # Arbitrary 2D container for now.
        cell_grid = []

        if verbose:
            print('Constructing grid for patch "{}" with dimensions (np, nr) = ({}, {})'.format(self.patch_name, np_cells, nr_cells))
            print(np_cells)
            print(nr_cells)
        np_lines = np_cells + 1
        nr_lines = nr_cells + 1
        if verbose: print(' # Create B-Splines along the North and South boundaries.')
        # Create B-Splines along the North and South boundaries.
        N_vals = self.N.fluff()

        self.N_spl, uN = splprep([N_vals[0], N_vals[1]], s=0)
        # Reverse the orientation of the South line to line up with the North.

        S_vals = self.S.reverse_copy().fluff()
        self.S_spl, uS = splprep([S_vals[0], S_vals[1]], s=0)
        if verbose: print(' # Create B-Splines along West boundaries.')
        # Create B-Splines along the East and West boundaries.
        # Parameterize EW splines in Psi
        try:
            #Cannot fluff with too many points
            n = 500 if len(self.W.p) < 50 else 100
           # W_vals = self.W.reverse_copy().fluff(num = n)
            W_vals = self.W.reverse_copy().fluff(n, verbose=verbose)
            Psi = psi_parameterize(grid, W_vals[0], W_vals[1])
            self.W_spl, uW = splprep([W_vals[0], W_vals[1]], u=Psi, s=10)
        except Exception as e:
            exc_type, exc_obj, tb = sys.exc_info()
            f = tb.tb_frame
            lineno = tb.tb_lineno
            filename = f.f_code.co_filename
            linecache.checkcache(filename)
            line = linecache.getline(filename, lineno, f.f_globals)
            print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

        if verbose: print(' # Create B-Splines along the East boundaries.')
        try:
            n = 500 if len(self.E.p) < 50 else 100
            E_vals = self.E.fluff(num=n)
            self.E_spl, uE = splprep([E_vals[0], E_vals[1]], u=psi_parameterize(grid, E_vals[0], E_vals[1]), s=10)
        except Exception as e:
            print(' Number of points on the boundary:', len(self.E.p))
            plt.plot(E_vals[0], E_vals[1], '.', color='black')
            print(repr(e))
        if verbose: print(' #check plate_patch')

        if self.plate_patch:

            def f(u, *args):
                _y = splev(u, args[0])
                return grid.PsiNorm.get_psi(args[1], args[2]) - grid.PsiNorm.get_psi(_y[0], _y[1])

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return array[idx]

            if self.plate_location == 'W':
                _u = uW
                U_vals = W_vals
                U_spl = self.W_spl
                plate_north = [N_vals[0][0], N_vals[1][0]]
                plate_south = [S_vals[0][0], S_vals[1][0]]
            elif self.plate_location == 'E':
                _u = uE
                U_vals = E_vals
                U_spl = self.E_spl
                plate_north = [N_vals[0][-1], N_vals[1][-1]]
                plate_south = [S_vals[0][-1], S_vals[1][-1]]

            if verbose:
                print('=' * 80 + '\n')
                print('Parameterization of Target Plate in PSI:')
                print(_u)
                print('=' * 80 + '\n')

            lookup = {}
            for i in range(len(_u)):
                lookup[_u[i]] = i
            try:
                plate_north_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args=(U_spl, plate_north[0], plate_north[1])))]
            except ValueError:
                try:
                    plate_north_index = lookup[find_nearest(_u, fsolve(f, 0, args=(U_spl, plate_north[0], plate_north[1])))]
                except:
                    if verbose: print('NorthIndex: ERROR IN PARAMETERIZATION IN PSI')
            try:
                plate_south_index = lookup[find_nearest(_u, brentq(f, _u[0], _u[-1], args=(U_spl, plate_south[0], plate_south[1])))]
            except ValueError:
                try:
                    plate_south_index = lookup[find_nearest(_u, fsolve(f, 0, args=(U_spl, plate_south[0], plate_south[1])))]
                except:
                    if verbose: print('SouthIndex: ERROR IN PARAMETERIZATION IN PSI')
            if plate_south_index > plate_north_index:
                plate_north_index, plate_south_index = plate_south_index, plate_north_index

            U_vals = [U_vals[0][plate_south_index:plate_north_index + 1], U_vals[1][plate_south_index:plate_north_index + 1]]
            U_spl, _u = splprep([U_vals[0], U_vals[1]], u=psi_parameterize(grid, U_vals[0], U_vals[1]), s=0)

            if self.plate_location == 'W':
                W_vals = U_vals
                self.W_spl = U_spl
                uW = _u
            elif self.plate_location == 'E':
                E_vals = U_vals
                self.E_spl = U_spl
                uE = _u

        # Generate our sub-grid anchor points along the North
        # and South boundaries of our patch.
        self.N_vertices = []
        self.S_vertices = []
        self.E_vertices = []
        self.W_vertices = []

        if verbose: print('# Generate our sub-grid anchor points along the North and South boundaries of our patch.')
        # and South boundaries of our patch')

        if self.BoundaryPoints.get('N') is None:
            for i in range(np_lines):
                _n = splev(_poloidal_f(i / (np_lines - 1)), self.N_spl)
                self.N_vertices.append(Point((float(_n[0]), float(_n[1]))))
        else:
            if verbose: print('Find boundary points at face "N" for {}:{}'.format(self.patch_name, self.BoundaryPoints.get('N')))
            self.N_vertices = self.BoundaryPoints.get('N')

        if self.BoundaryPoints.get('S') is None:
            for i in range(np_lines):
                _s = splev(_poloidal_f(i / (np_lines - 1)), self.S_spl)
                self.S_vertices.append(Point((float(_s[0]), float(_s[1]))))
        else:
            self.S_vertices = self.BoundaryPoints.get('S')

        u = [_radial_f(i / (nr_lines - 1)) for i in range(nr_lines)]

        Npts = 1000
        xy = splev(np.linspace(0, 1, Npts), self.E_spl)

        if self.BoundaryPoints.get('W') is None:
            for i in range(nr_lines):
                _w = splev(u[i], self.W_spl)
                self.W_vertices.append(Point((float(_w[0]), float(_w[1]))))
        else:
            self.W_vertices = self.BoundaryPoints.get('W')

        if self.BoundaryPoints.get('E') is None:
            for i in range(nr_lines):
                _e = splev(u[i], self.E_spl)
                self.E_vertices.append(Point((float(_e[0]), float(_e[1]))))
        else:
            self.E_vertices = self.BoundaryPoints.get('E')

        DetailedVertices = False
        if not isinstance(ShowVertices, bool):
            if ShowVertices == 'detailed':
                ShowVertices = True
                DetailedVertices = True
            else:
                raise ValueError('Unknow type for ShowVertices: must be True,False or "detailed"', ShowVertices)

        if visual or ShowVertices:
            if not DetailedVertices:
                markers = ['o'] * 4
                ms = 8
                color = 'black'

            else:
                color = self.color
                ms = 6
                markers = ['o', 'X', 's', 'D']
            for vertices, mark in zip([self.W_vertices, self.E_vertices, self.N_vertices, self.S_vertices], markers):
                for p in vertices:
                    plt.plot(p.x, p.y, '.', color=color, markersize=ms, marker=mark, markeredgecolor='black')
        # Radial lines of Psi surfaces. Ordered with increasing magnitude, starting with
        # the South boundary of the current Patch, and ending with the North boundary of
        # this current Patch. These will serve as delimiters when constructing cells.
        radial_lines = [self.N]
        radial_vertices = [self.N_vertices]
        dynamic_step = set_dynamic_step()
        if verbose: print('# Interpolate radial lines between North and South patch boundaries..')
        # Interpolate radial lines between North and South patch boundaries.
        self.radial_spl = []
        temp_vertices = []
        temp_vertices.append(self.N.p[-1])
        for i in range(len(self.W_vertices) - 2):
            #TODO: parallelize tracing of radial lines (draw_line function must be "externalized" in the scope of the script)
            radial_lines.append(grid.LineTracer.draw_line(self.W_vertices[i + 1], {'line': self.E}, option='theta',
                direction='cw', show_plot=visual, dynamic_step=dynamic_step, text=verbose))
            temp_vertices.append(radial_lines[-1].p[-1])
            radial_vals = radial_lines[i + 1].fluff(1000)
            Radial_spl, uR = splprep([radial_vals[0], radial_vals[1]], s=0)
            self.radial_spl.append(Radial_spl)
            vertex_list = []
            for j in range(np_lines):
                u = _poloidal_f(j / (np_lines - 1))
                _r = splev(u, self.radial_spl[i])
                Pt = Point((float(_r[0]), float(_r[1])))
                if self.distortion_correction['active'] and j > 0 and j < np_lines - 1:
                    Res = self.distortion_correction['resolution']
                    ThetaMin = self.distortion_correction['theta_min']
                    ThetaMax = self.distortion_correction['theta_max']
                    umin = _poloidal_f((j - 1) / (np_lines - 1))
                    umax = _poloidal_f((j + 1) / (np_lines - 1))
                    Pt1 = radial_vertices[i][j]
                    Pt2 = radial_vertices[i][j - 1]
                    Tag = '---- Correcting points: {},{}'.format(i, j)
                    Pt = CorrectDistortion(u, Pt, Pt1, Pt2, self.radial_spl[i], ThetaMin, ThetaMax, umin, umax, Res, visual, Tag, verbose)

                vertex_list.append(Pt)
                if visual:
                    for p in vertex_list:
                        plt.plot(p.x, p.y, '.', color='black', markersize=8)

            radial_vertices.append(vertex_list)
        radial_lines.append(self.S)
        temp_vertices.append(self.S.p[0])
        self.E_vertices = temp_vertices

        #Correct point on south boundary
        if self.distortion_correction['active']:
            for j in range(1, np_lines - 1):
                u = _poloidal_f(j / (np_lines - 1))
                Pt = self.S_vertices[j]
                Res = self.distortion_correction['resolution']
                ThetaMin = self.distortion_correction['theta_min']
                ThetaMax = self.distortion_correction['theta_max']
                umin = _poloidal_f((j - 1) / (np_lines - 1))
                umax = _poloidal_f((j + 1) / (np_lines - 1))
                Pt1 = radial_vertices[-1][j]
                Pt2 = radial_vertices[-1][j - 1]
                Tag = '---- Correcting south boundary points:{}'.format(j)
                self.S_vertices[j] = CorrectDistortion(u, Pt, Pt1, Pt2, self.S_spl, ThetaMin, ThetaMax, umin, umax, Res, visual, Tag, verbose)
        radial_vertices.append(self.S_vertices)
        if verbose: print('# Create Cells: South boundary -> North Boundary')
        # Create Cells: South boundary -> North Boundary
        for i in range(len(radial_lines)):
            if radial_lines[i] is self.S:
                break
            cell_grid.append([])
            for j in range(len(radial_vertices[i]) - 1):
                NW = radial_vertices[i][j]
                NE = radial_vertices[i][j + 1]
                SW = radial_vertices[i + 1][j]
                SE = radial_vertices[i + 1][j + 1]
                cell_grid[i].append(Cell([Line([NW, NE]), Line([SW, SE]), Line([SE, NE]), Line([SW, NW])]))

        self.cell_grid = cell_grid

    def CheckPatch(self, grid, verbose=False):
        if verbose: print(' # Checking if patch boundaries can be interpolated wiht splines')

        def psi_parameterize(grid, r, z):
            """
            Helper function to be used to generate a
            list of values to parameterize a spline
            in Psi. Returns a list to be used for splprep only
            """
            vmax = grid.PsiNorm.get_psi(r[-1], z[-1])
            vmin = grid.PsiNorm.get_psi(r[0], z[0])

            vparameterization = np.empty(0)
            for i in range(len(r)):
                vcurr = grid.PsiNorm.get_psi(r[i], z[i])
                vparameterization = np.append(vparameterization, abs((vcurr - vmin) / (vmax - vmin)))

            return vparameterization

        def IsMonotonic(psi, x, y, Str):
            if not (non_increasing(psi) or non_decreasing(psi)):
                print(Str + ' is not monotonic')
                d = which_non_increasing(psi)
                u = which_increasing(psi)
                plt.plot(x[u], y[u], 'o', 'r')
                plt.plot(x[d], y[d], 'o', 'b')
                raise ValueError(Str + ' is not monotonic')
            else:
               print(Str + ' is monotonic')

       # N_vals = self.N.fluff()
        #S_vals = self.S.reverse_copy().fluff()
        n = 20 if len(self.W.p) > 500 else 100
        W_vals = self.W.reverse_copy().fluff(n)
        n = 20 if len(self.E.p) > 500 else 100
        E_vals = self.E.fluff(num=n)
        if verbose: print(' ## Getting Psi values along the boundaries')
        #PsiN=psi_parameterize(grid, N_vals[0], N_vals[1])
        #PsiS=psi_parameterize(grid, S_vals[0], S_vals[1])
        PsiW = psi_parameterize(grid, W_vals[0], W_vals[1])
        PsiE = psi_parameterize(grid, E_vals[0], E_vals[1])
        if verbose: print(' ## Checking monoticity of Psi values along the boundaries')
        IsMonotonic(PsiW, W_vals[0], W_vals[1], 'PsiW')
        IsMonotonic(PsiE, E_vals[0], E_vals[1], 'PsiE')


def strictly_decreasing(L: 'array-like') -> bool:
    """
    Determine if strictly decreasing.

    Parameters
    ----------
    L : array-like
        Values to test.

    Returns
    -------
        True if strictly decreasing and False otherwise
    """
    return all(x > y for x, y in zip(L, L[1:]))


def non_increasing(L: 'array-like') -> bool:
    """
    Determine if non-increasing.

    Parameters
    ----------
    L : array-like
        Values to test.

    Returns
    -------
        True if non-increasing and False otherwise
    """
    return all(x >= y for x, y in zip(L, L[1:]))

def which_non_increasing(L: 'array-like') -> list:
    """
    Determine non-increasing values.

    Parameters
    ----------
    L : array-like
        Values to test.

    Returns
    -------
        A list of 2-tuples containing index and non-increasing element.
    """
    return [i for i, (x, y) in enumerate(zip(L, L[1:])) if x > y]

def which_increasing(L: 'array-like') -> list:
    """
    Determine increasing values.

    Parameters
    ----------
    L : array-like
        Values to test.

    Returns
    -------
        A list of 2-tuples containing index and increasing element
    """
    return [i for i, (x, y) in enumerate(zip(L, L[1:])) if x < y]

def non_decreasing(L: 'array-like') -> bool:
    """
    Determine if non-decreasing.

    Parameters
    ----------
    L : array-like
        Values to test.

    Returns
    -------
        True if non-decreasing and False otherwise
    """
    return all(x <= y for x, y in zip(L, L[1:]))

def find_split_index(split_point: Point, line: Line) -> tuple:
    """
    Determine which index a Point would best split a Line.

    This method is useful for modifying Line objects during line tracing
    (searching for intersection of two line objects and trimming excess)

    Parameters
    ----------
    split_point : Point
        The candidate Point to find the split index with respect to.

    line : Line
        The Line to search for a split index within.

    Returns
    -------
        A 2-tuple containing the split-index and a boolean flag indicating
            whether the split_point was contained within the Line

    Notes
    -----
    Should no appropriate split index be found, the method will return a
    None value in place of an integer index.

    The second entry of the tuple return value would be a value of True if
    the `split_point` parameter was used to define the `line` parameter.
    """
    same_line_split = False
    for i in range(len(line) - 1):
        # Split point is exactly on a Line object's point. Occurs often
        # when splitting a Line object with itself.

        if split_point[-1] == line[i][-1] and split_point[0] == line[i][0]:
            same_line_split = True
            return i, same_line_split
        # Create two vectors.
        end_u   = line[i + 1] - line[i]
        split_v = split_point - line.p[i]

        if is_between(end_u, split_v):
            # store index corresponding to the start of the segment containing the split_point.
            return i, same_line_split
        else:
            continue
    return None, same_line_split

def is_between(end_u: 'array-like', split_v: 'array-like') -> bool:
    eps = 1e-9
    # check cross product vector norm against eps.
    if np.linalg.norm(np.cross(end_u, split_v)) < eps:
        # Colinear up to eps.
        # Ensure dot product is positive and vector v lies in distance of u norm.
        if (np.dot(end_u, split_v) > 0) and (np.linalg.norm(end_u) > np.linalg.norm(split_v)):
            return True
    else:
        return False

def rotmatrix(theta: float) -> 'numpy.ndarray':
    """
    Construct a rotation matrix

    Parameters
    ----------
    theta : float
        Angle in radians.

    Returns
    -------
        An ndarray with shape (2, 2) representing a 2D rotation matrix.
    """
    rot = np.zeros((2, 2))
    rot[0, 0] = np.cos(theta)
    rot[1, 1] = np.cos(theta)
    rot[0, 1] = -np.sin(theta)
    rot[1, 0] = np.sin(theta)
    return rot

def rotate(vec, theta, origin):
    return np.matmul(rotmatrix(theta), vec - origin) + origin

def unit_vector(v):
    """
    Returns unit vector
    """
    return v / np.linalg.norm(v)

def angle_between(u, v, origin, relative=False):
    """
    Compute angle in radians between vectors u and v
    """
    u_norm = unit_vector(u - origin)
    v_norm = unit_vector(v - origin)
    angle = np.arccos(np.clip(np.dot(u_norm, v_norm), -1, 1))
    return orientation_between(u, v, origin) * angle if relative else angle

def orientation_between(u, v, origin):
    """
    Compute angle in radians between vectors u and v
    """
    u_norm = unit_vector(u - origin)
    v_norm = unit_vector(v - origin)
    return np.sign(np.arctan2(u_norm[0] * v_norm[1] - u_norm[1] * v_norm[0],
                   u_norm[0] * v_norm[0] + u_norm[0] * v_norm[1]))

def reorder_limiter(new_start, limiter):
    start_index, = find_split_index(new_start, limiter)
    return limiter

def limiter_split(start, end, limiter):
    start_index, sls = find_split_index(start, limiter)
    end_index, sls = find_split_index(end, limiter)
    if end_index <= start_index:
        limiter = Line(limiter.points[start_index:] + limiter.points[:start_index + 1])
    return limiter

def trim_geometry(geoline, start, end):
    try:
        trim = (geoline.split(start)[1]).split(end, add_split_point=True)[0]
    except:
        trim = limiter_split(start, end, geoline)
    return trim

def CorrectDistortion(u, Pt, Pt1, Pt2, spl, ThetaMin, ThetaMax, umin, umax, Resolution, visual, Tag, MinTol=1.02, MaxTol=0.98, Verbose=False):
    MaxIter = Resolution * 10
    dumax = (umax - u) / Resolution
    dumin = (u - umin) / Resolution
    Theta = Line([Pt1, Pt]).GetAngle(Line([Pt1, Pt2]))
    if Verbose: print(f'Current theta value: {Theta}')
    if Theta < ThetaMin or Theta > ThetaMax:
        if Verbose: print('{}: u={};Theta={};ThetaMin={};ThetaMax={}'.format(Tag, u, Theta, ThetaMin, ThetaMax))
        if visual:
                plt.plot(Pt.x, Pt.y, '.', color='red', markersize=8, marker='o')
                plt.show()
                plt.draw()
        icount = 0
        color = 'purple'
        while Theta < ThetaMin or Theta > ThetaMax:
            icount += 1
            if Theta < ThetaMin:
               u = u + dumax
            elif Theta > ThetaMax:
               u = u - dumin
            if (u > umax * MaxTol or u < umin * MinTol) or icount >= MaxIter:
               if Verbose: print('[{}]>>>> umax={} umin={};u:{};Theta={}'.format(icount, umax, umin, u, Theta))
               color = 'gray'
               break
            _r = splev(u, spl)
            Pt = Point((_r[0], _r[1]))
            Theta = Line([Pt1, Pt]).GetAngle(Line([Pt1, Pt2]))
        if Verbose: print('[{}]>>>> u:{};Theta={}'.format(icount, u, Theta))
        if visual:
            plt.plot(Pt.x, Pt.y, '.', color=color, markersize=8, marker='s')
    return Pt

def calc_mid_point(v1, v2):
    """
    Calculates the bisection of two vectors of equal length,
    and returns the point on the circle at that angle.


    Parameters
    ----------
    v1 : geometry.Vector
        v1 must be furthest right in a counter clockwise direction.
    v2 : geometry.Vector
        Vector on the left.

    Returns
    -------
    tuple
        The point at the bisection of two vectors.

    """
    # bisection
    theta = np.arccos(np.dot(v1.arr(), v2.arr()) / (v1.mag() * v2.mag())) / 2.

    # check with quadrant the vectors are in and
    # compute the angles appropriately
    if v1.quadrant == (1, 1):
        # NE quadrant
        angle = np.arccos(v1.xnorm / v1.mag())
    elif v1.quadrant == (-1, 1):
        # NW quadrant
        angle = np.pi - np.arcsin(v1.ynorm / v1.mag())
    elif v1.quadrant == (-1, -1):
        # SW quadrant
        angle = np.pi + np.arctan(v1.ynorm / v1.xnorm)
    elif v1.quadrant == (1, -1):
        # SE quadrant
        angle = - np.arccos(v1.xnorm / v1.mag())
    else:
        print("Something went wrong")

    x = v1.xorigin + v1.mag() * np.cos(theta + angle)
    y = v1.yorigin + v1.mag() * np.sin(theta + angle)
    return x, y

def test2points(p1, p2, line):
    """
    Check if two points are on opposite sides of a given line.

    Parameters
    ----------
    p1 : tuple
        First point, (x, y)
    p2 : tuple
        Second point, (x, y)
    line : array-like
        The line is comprised of two points ((x, y), (x, y)).
    Returns
    -------
    tuple
        Returns two numbers, if the signs are different
        the points are on opposite sides of the line.

    """
    for ind in range(len(line) - 1):
        (x1, y1), (x2, y2) = line[ind], line[ind + 1]
        x, y = p1
        a, b = p2
        d1 = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
        d2 = (a - x1) * (y2 - y1) - (b - y1) * (x2 - x1)

        if (np.sign(d1) != np.sign(d2)):
            return True
    return False

def intersect(line1, line2, verbose=False):
    """ Finds the intersection of two line segments


    Parameters
    ----------
    line1 : array-like
    line2 : array-like
        Both lines of the form A = ((x, y), (x, y)).

    Returns
    -------
    tuple
        Coordinates of the intersection.

    """
    # TODO: Imporove this method so it can handle verticle lines.
    # there is a division by zero error in some DIII-D data caused by this.
    def line(x, line):
        """ Point slope form. """
        (x1, y1), (x2, y2) = line
        if x2 - x1 == 0:
            x1 += 1e-4

        return (y2 - y1) / (x2 - x1) * (x - x1) + y1

    def f(xy):
        # parse the line, access any point
        # normalized to help solve for zero
        x, y = xy
        return np.array([y - line(x, line1),
                         y - line(x, test_line)])

    for ind in range(len(line2) - 1):
        # use the mean for initial guess
        test_line = [line2[ind], line2[ind + 1]]
        (a, b), (c, d) = line1
        (i, j), (p, q) = test_line
        guess = (np.mean([a, c, i, p]), np.mean([b, d, j, q]))
        sol, infoDict, ier, mesg = fsolve(f, guess, full_output=True)
        if verbose:
            print('{}: {}'.format(ind, mesg))
        if ier == 1:
            break
    return sol[0], sol[1]

def segment_intersect(line1, line2, verbose=False):
    """ Finds the intersection of two FINITE line segments.
    Parameters
    ----------
    line1 : array-like
    line2 : array-like
        Both lines of the form line1 = (P1, P2), line2 = (P3, P4)
    Returns
    -------
    bool, tuple
        True/False of whether the segments intersect
        Coordinates of the intersection
    """
    (xa, ya), (xb, yb) = line1

    for i in range(len(line2) - 1):
        (xc, yc), (xd, yd) = line2[i], line2[i + 1]

        M = np.array([[xb - xa, -xd + xc], [yb - ya, -yd + yc]])
        r = np.array([xc - xa, yc - ya])
        try:
            sol = np.linalg.solve(M, r)
        except np.linalg.LinAlgError:
            continue

        if (sol[0] <= 1) and (sol[1] <= 1) \
                and (sol[0] >= 0) and (sol[1] >= 0):
            return True, [(xc, yc), (xc + sol[1] * (xd - xc), yc + sol[1] * (yd - yc))]
    return False, [(np.nan, np.nan), (np.nan, np.nan)]

def UnfoldLabel(Dic: dict, Name: str) -> str:
    """
    Unfold Patch label (e.g. "C1" -> "Inner Core Top")

    Parameters
    ----------
    Dic : dict
        Dictionnary containing description of acronym characters
    Name : str
        patch label

    Returns
    -------
    str
        Unfolded patch label.

    """
    Output = []
    for s in Name:
        if Dic.get(s) is not None:
            Output.append(Dic.get(s) + ' ')
        else:
            Output.append(s)
    if len(Output) > 0:
        return ''.join(Output)
    else:
        return ''
