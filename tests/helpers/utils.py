from pathlib import Path

def resolve_path(path: str, as_str: bool = True):
    """
    Helper function for resolving paths.

    Parameters
    ----------
    path : str, pathlib.Path
        Path to fully resolve
    as_str : bool, optional
        Bool indicating to return resolved path as a str
        or PosixPath obj
    
    Returns
    -------
        Resolved path
    """
    path = Path(path).resolve()
    return str(path) if as_str else path

