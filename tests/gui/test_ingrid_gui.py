import pytest
import time
from INGRID.ingrid import Ingrid
from INGRID.exceptions import TkInitializationSuccess

def test_ingrid_gui_initialization():
    """
    A simple tk intialization test that the IngridGUI class was
    initialized and was able to run the mainloop and destroy methods.

    This test expects the TkInitializationSuccess exception to be raised
    """
    session = Ingrid()
    with pytest.raises(TkInitializationSuccess):
        session.StartGUI(test_initialization=True)

