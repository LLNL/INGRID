import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from scipy.interpolate import RectBivariateSpline as rbs
from typing import Optional, Tuple, Dict, Union
from pathlib import Path
from freeqdsk import geqdsk
from scipy.optimize._optimize import OptimizeResult
from scipy.optimize._optimize import OptimizeResult

class EfitData:
    """
    Structure to store the rectangular grid of psi data. It uses
    cylindrical coordinates, where R and Z are similar to the cartesian
    x and y. The phi component goes away due to the symmetry of a
    tokamak.

    Parameters
    ----------
    rmin: float, optional
        Left boundary of the grid.
    rmax: float, optional
        Right boundary of the grid.
    nr: int, optional
        Number of grid points in the R direction.
    zmin: float, optional
        Bottom boundary for the grid.
    zmax: float, optional
        Top boundary for the grid.
    nz: int, optional
        Number of grid points in the Z direction.
    rcenter: float, optional
        R coordinate of the magnetic axis.
    bcenter: float, optional
        Magnetic field at the magnetic axis.
    rlimiter: Optional[float], optional
        R coordinate of the limiter.
    zlimiter: Optional[float], optional
        Z coordinate of the limiter.
    rmagx: float, optional
        R coordinate of the magnetic axis.
    zmagx: float, optional
        Z coordinate of the magnetic axis.
    name: str, optional
        Title of the figure the data will be plotted on.
    parent: Optional[object], optional
        Parent object reference.
    """
    
    geqdsk_data: geqdsk.GEQDSKFile

    def __init__(self,
                 geqdsk_data: Union[str, Path, geqdsk.GEQDSKFile],
                 name: str = 'unnamed', 
                 parent: Optional[object] = None,
                 psi_min: Optional[float] = None,
                 psi_max: Optional[float] = None):

        #
        # Read the GEQDSK file if it is a path
        #
        if isinstance(geqdsk_data, (str, Path)):
            with open(geqdsk_data, 'r') as f:
                geqdsk_data = geqdsk.read(f)

        self.geqdsk_data = geqdsk_data

        #
        # Extract quantities needed to initialize EfitData class
        #
        nr       = geqdsk_data['nx']
        nz       = geqdsk_data['ny']
        rdim     = geqdsk_data['rdim']
        zdim     = geqdsk_data['zdim']
        zmid     = geqdsk_data['zmid']
        rleft    = geqdsk_data['rleft']
        psi      = geqdsk_data['psi']

        #
        # Derived values
        #
        rmin = rleft
        rmax = rmin + rdim
        zmin = (zmid - 0.5 * zdim)
        zmax = zmin + zdim

        #
        # Create the grid
        #
        r, dr = np.linspace(rmin, rmax, nr, retstep=True)
        z, dz = np.linspace(zmin, zmax, nz, retstep=True)
        rgrid, zgrid = np.meshgrid(r, z, indexing='ij')
        value = np.zeros([nr, nz])

        self.rmin: float = rmin
        self.rmax: float = rmax
        self.zmin: float = zmin
        self.zmax: float = zmax
        self.r: np.ndarray = rgrid
        self.z: np.ndarray = zgrid
        self.v: np.ndarray = value
        self.dr: float = dr
        self.dz: float = dz
        self.name: str = name
        self.parent: Optional[object] = parent
        self.psi_levels: Dict[str, plt.contour] = {}

        #
        # Initialize the bivariate spline
        #
        if psi_min is not None and psi_max is not None:
            psi = (psi - psi_min) / (psi_max - psi_min)
        self.init_bivariate_spline(rgrid, zgrid, psi)

    def init_bivariate_spline(self, r: np.ndarray, 
                              z: np.ndarray, 
                              v: np.ndarray) -> None:
        """
        Initialize scipy.interpolate.RectBivariateSpline
        object for Bicubic interpolation.

        Sets class member v to crude EFIT grid.

        Parameters
        ----------
        r: np.ndarray
            1-D array of r coordinates in strictly ascending order.
        z: np.ndarray
            1-D array of z coordinates in strictly ascending order.
        v: np.ndarray
            2-D array of EFIT data with shape (r.shape, z.shape).
        """
        self.v = v  # Crude EFIT grid.
        self.rbs = rbs(r, z, v)  # RectBivariateSpline object.

    def Gradient(self, xy: Tuple[float, float]) -> np.ndarray:
        """
        Combines the first partial derivatives to solve the system for
        maximum, minimum, and saddle locations.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates. Ex: xy = (x0, y0).

        Returns
        -------
        np.ndarray
            Vector function to be used in find root.
        """
        x, y = xy
        F = np.zeros(2)
        F[0] = self.get_psi(x, y, tag='vr')
        F[1] = self.get_psi(x, y, tag='vz')
        return F

    def Hessian(self, xy: Tuple[float, float]) -> np.ndarray:
        """
        Compute the Hessian at a point.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates. Ex: xy = (x0, y0).

        Returns
        -------
        np.ndarray
            Numpy array of shape (2, 2) representing the Hessian at xy.
        """
        x, y = xy
        H = np.zeros((2, 2))
        H[0, 0] = self.get_psi(x, y, 'vrr')
        H[1, 1] = self.get_psi(x, y, 'vzz')
        H[0, 1] = self.get_psi(x, y, 'vrz')
        H[1, 0] = self.get_psi(x, y, 'vrz')
        return H

    def PsiFunction(self, xy: Tuple[float, float]) -> float:
        """
        Compute the psi value at a given point.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates.

        Returns
        -------
        float
            Psi value at the given point.
        """
        x, y = xy
        return self.get_psi(x, y)

    def get_psi(self, r0: float, z0: float, tag: str = 'v') -> float:
        """
        Find grid cell encompassing (r0,z0) and use Bicubic Interpolation
        to calculate the exact value at the point. Useful for
        finding information between grid points.

        Parameters
        ----------
        r0: float
            R coordinate of the point of interest.
        z0: float
            Z coordinate of the point of interest.
        tag: str, optional
            Type of derivative we want: v, vr, vz, vrz.
            If nothing is provided, it assumes no derivative (v).

        Returns
        -------
        float
            Value of psi or its derivative at the specified coordinate.
        """
        lookup = {'v': (0, 0), 'vr': (1, 0), 'vrr': (2, 0),
                  'vz': (0, 1), 'vzz': (0, 2), 'vrz': (1, 1)}

        dx, dy = lookup[tag]
        return self.rbs(r0, z0, dx, dy)[0]
    
    def __call__(self, r0: float, z0: float, tag: str = 'v') -> float:
        """
        Call the bivariate spline to get the value of psi at a point.
        """
        return self.get_psi(r0=r0, z0=z0, tag=tag)

    def plot_levels(self, level: float = 1.0, color: str = 'red') -> None:
        """
        Quickly visualize a particular line of constant psi.

        This function is useful for a quick view but cannot store
        points of intersection and cannot be generalized. For a segment
        of psi, use the draw_lines method in the line tracing class.

        Parameters
        ----------
        level: float, optional
            Value of psi you wish to see.
        color: str, optional
            Color of the line.
        """
        #
        # Draw contour line on top of existing figure
        #
        level = float(level)
        self.ax.contour(self.r, self.z, self.v, level, colors=color)

    def PlotLevel(self, level: float = 1.0, color: str = 'red', label: str = '', 
                  linestyles: str = 'solid', refined: bool = True, 
                  refine_factor: int = 10) -> None:
        """
        Plot a psi level and provide it a label.

        This function is useful for management of psi boundaries
        such as 'psi_pf', 'psi_core', etc and ensuring the contour will
        be properly replotted (no duplicate of same label).

        Parameters
        ----------
        level: float, optional
            Psi level to plot. Default is 1.0 (separatrix of normalized psi).
        color: str, optional
            Color to pass to matplotlib contour function.
        label: str, optional
            Label to associate with the psi level.
        linestyles: str, optional
            Line style to pass to matplotlib contour function.
        refined: bool, optional
            Plot level with high-resolution cubic spline representation.
        refine_factor: int, optional
            Refinement factor to be passed to SciPy zoom method.
        """
        data = self.v
        rgrid = self.r
        zgrid = self.z

        if refined:
            data = zoom(input=self.v, zoom=refine_factor)
            rgrid, zgrid = np.meshgrid(np.linspace(self.rmin, self.rmax, data.shape[0]),
                                       np.linspace(self.zmin, self.zmax, data.shape[1]),
                                       indexing='ij')
        try:
            self.psi_levels[label].collections[0].remove()
            self.psi_levels[label] = plt.contour(rgrid, zgrid, data, [float(level)], 
                                                 colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)
        except KeyError:
            self.psi_levels[label] = plt.contour(rgrid, zgrid, data, [float(level)], 
                                                 colors=color, label=label, linestyles=linestyles)
            self.psi_levels[label].collections[0].set_label(label)

    def plot_data(self, nlevs: int = 30, interactive: bool = False, 
                  fig: Optional[plt.Figure] = None, ax: Optional[plt.Axes] = None, 
                  view_mode: str = 'filled', refined: bool = True, 
                  refine_factor: int = 10) -> None:
        """
        Plot the EFIT data.

        Visualizes eqdsk file with either contour lines or filled contours.

        Parameters
        ----------
        nlevs: int, optional
            Number of levels to be plotted.
        interactive: bool, optional
            Set matplotlib interactive mode on or off.
        fig: Optional[plt.Figure], optional
            Matplotlib figure handle.
        ax: Optional[plt.Axes], optional
            Matplotlib axes handle.
        view_mode: str, optional
            Represent EFIT data with standard contour lines or filled contour lines.
            'filled' enables filled contours, 'lines' omits filling of contours.
        refined: bool, optional
            Plot level with high-resolution cubic spline representation.
        refine_factor: int, optional
            Refinement factor to be passed to SciPy zoom method.
        """
        lev = self.v.min() + (self.v.max() - self.v.min()) * np.arange(nlevs) / (nlevs - 1)
        self.fig = fig if fig is not None else plt.figure('INGRID: ' + self.name, figsize=(8, 10))
        self.fig.subplots_adjust(bottom=0.075)
        self.ax = ax if ax is not None else self.fig.add_subplot(111)

        data = self.v
        rgrid = self.r
        zgrid = self.z

        if refined:
            data = zoom(input=self.v, zoom=refine_factor)
            rgrid, zgrid = np.meshgrid(np.linspace(self.rmin, self.rmax, data.shape[0]),
                                       np.linspace(self.zmin, self.zmax, data.shape[1]),
                                       indexing='ij')
        if view_mode == 'lines':
            self.ax.contour(rgrid, zgrid, data, lev, cmap='gist_gray')
        elif view_mode == 'filled':
            self.ax.contourf(rgrid, zgrid, data, lev, cmap='gist_gray')
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlabel('R')
        self.ax.set_ylabel('Z')
        self.ax.set_xlim(self.rmin, self.rmax)
        self.ax.set_ylim(self.zmin, self.zmax)
        if interactive:
            plt.ion()
        self.fig.show()

    def clear_plot(self) -> None:
        """
        Clear the current plot if it exists.
        """
        if plt.get_fignums():
            plt.clf()

class GEQDSKInterpolator:
    """
    Interpolator for GEQDSK data.
    """
    geqdsk_data: geqdsk.GEQDSKFile

    def __init__(self,
                 geqdsk_data: Union[str, Path, geqdsk.GEQDSKFile],
                 name: str = 'unnamed', 
                 parent: Optional[object] = None,
                 psi_min: Optional[float] = None,
                 psi_max: Optional[float] = None):

        #
        # Read the GEQDSK file if it is a path
        #
        if isinstance(geqdsk_data, (str, Path)):
            with open(geqdsk_data, 'r') as f:
                geqdsk_data = geqdsk.read(f)

        self.geqdsk_data = geqdsk_data

        #
        # Extract quantities needed to initialize EfitData class
        #]
        rdim     = geqdsk_data['rdim']
        zdim     = geqdsk_data['zdim']
        zmid     = geqdsk_data['zmid']
        rleft    = geqdsk_data['rleft']
        psi      = geqdsk_data['psi']

        #
        # Derived values
        #
        self.rmin = rleft
        self.rmax = self.rmin + rdim
        self.zmin = (zmid - 0.5 * zdim)
        self.zmax = self.zmin + zdim

        #
        # Create the grid
        #
        self.nx = geqdsk_data['nx']
        self.ny = geqdsk_data['ny']
        self.xrange = np.linspace(self.rmin, self.rmax, self.nx)
        self.yrange = np.linspace(self.zmin, self.zmax, self.ny)

        self.name: str = name
        self.parent: Optional[object] = parent

        #
        # Initialize the bivariate spline
        #
        if psi_min is not None and psi_max is not None:
            psi = (psi - psi_min) / (psi_max - psi_min)
        self.init_bivariate_spline(self.xrange, self.yrange, psi)


    def init_bivariate_spline(self, r: np.ndarray, 
                              z: np.ndarray, 
                              v: np.ndarray) -> None:
        """
        Initialize scipy.interpolate.RectBivariateSpline
        object for Bicubic interpolation.

        Sets class member v to crude EFIT grid.

        Parameters
        ----------
        r: np.ndarray
            1-D array of r coordinates in strictly ascending order.
        z: np.ndarray
            1-D array of z coordinates in strictly ascending order.
        v: np.ndarray
            2-D array of EFIT data with shape (r.shape, z.shape).
        """
        self.rbs = rbs(r, z, v)  # RectBivariateSpline object.

    def generate_normalized_interpolator(self, psi_min: float, psi_max: float) -> 'GEQDSKInterpolator':
        """
        Generate a normalized interpolator.
        """
        return GEQDSKInterpolator(geqdsk_data=self.geqdsk_data,
                                  name=self.name,
                                  parent=self.parent,
                                  psi_min=psi_min,
                                  psi_max=psi_max)

    def Gradient(self, xy: Tuple[float, float]) -> np.ndarray:
        """
        Combines the first partial derivatives to solve the system for
        maximum, minimum, and saddle locations.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates. Ex: xy = (x0, y0).

        Returns
        -------
        np.ndarray
            Vector function to be used in find root.
        """
        x, y = xy
        F = np.zeros(shape=(2, 1))
        F[0] = self.get_psi(x, y, tag='vr')
        F[1] = self.get_psi(x, y, tag='vz')
        return F.flatten()

    def Hessian(self, xy: Tuple[float, float]) -> np.ndarray:
        """
        Compute the Hessian at a point.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates. Ex: xy = (x0, y0).

        Returns
        -------
        np.ndarray
            Numpy array of shape (2, 2) representing the Hessian at xy.
        """
        x, y = xy
        H = np.zeros((2, 2))
        H[0, 0] = self.get_psi(x, y, 'vrr')
        H[1, 1] = self.get_psi(x, y, 'vzz')
        H[0, 1] = self.get_psi(x, y, 'vrz')
        H[1, 0] = self.get_psi(x, y, 'vrz')
        return H

    def PsiFunction(self, xy: Tuple[float, float]) -> float:
        """
        Compute the psi value at a given point.

        Parameters
        ----------
        xy: Tuple[float, float]
            Contains x and y coordinates.

        Returns
        -------
        float
            Psi value at the given point.
        """
        x, y = xy
        return self.get_psi(x, y)

    def get_psi(self, r0: float, z0: float, tag: str = 'v') -> float:
        """
        Find grid cell encompassing (r0,z0) and use Bicubic Interpolation
        to calculate the exact value at the point. Useful for
        finding information between grid points.

        Parameters
        ----------
        r0: float
            R coordinate of the point of interest.
        z0: float
            Z coordinate of the point of interest.
        tag: str, optional
            Type of derivative we want: v, vr, vz, vrz.
            If nothing is provided, it assumes no derivative (v).

        Returns
        -------
        float
            Value of psi or its derivative at the specified coordinate.
        """
        lookup = {'v': (0, 0), 'vr': (1, 0), 'vrr': (2, 0),
                  'vz': (0, 1), 'vzz': (0, 2), 'vrz': (1, 1)}

        dx, dy = lookup[tag] 

        value = self.rbs(r0, z0, dx, dy)[0]
        return value
    
    def __call__(self, r0: float, z0: float, tag: str = 'v') -> float:
        """
        Call the bivariate spline to get the value of psi at a point.
        """
        return self.get_psi(r0=r0, z0=z0, tag=tag)

    def refine_with_root_finder(self, r0: float, z0: float) -> Tuple[float, float]:
        """
        Refine the coordinates to the root of the gradient of psi.
        """
        from scipy.optimize import root
        sol: OptimizeResult = root(self.Gradient, [r0, z0])
        print(f'Found root at: {sol.x[0]}, {sol.x[1]}')
        print(f'Delta (refined - input): {sol.x[0] - r0}, {sol.x[1] - z0}')
        return sol.x[0], sol.x[1]
