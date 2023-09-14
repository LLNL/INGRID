from freeqdsk import geqdsk
from pathlib import Path

class EQLoader:
    """
    Loader for G-EQDSK files.

    This class will leverage the FreeQDSK G-EQDSK loader

    """
    def __init__(self, fpath: str) -> None:
        self.fpath = fpath

    def validate_eqdsk(self, fpath: str):
        """
        Helper function for checking provided fpath
        exists

        Parameters
        ----------
        fpath : str
            The fpath to validate exists
        
        Returns
        -------
            Bool indicating file existence
        """
        raise FileNotFoundError
    
    def load(self, fpath: str = None):
        """
        Read the G-EQDSK data in the fpath provided.

        If fpath is None, use the fpath provided at
        initialization. Otherwise, we will clobber the
        stored fpath with the new fpath and reload the internal
        data.

        Parameters
        ----------
        fpath : optional
            Filepath to load. If None, use initialization time fpath
        
        Returns
        -------
            A dict of loaded G-EQDSK data
        """

        if fpath is None:
            fpath = self.eqdsk_path
