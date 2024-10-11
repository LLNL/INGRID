import os
import sys
import pytest
from pathlib import Path
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))

def pytest_addoption(parser):
    """
    Allow for adding options to pytest command.
    """
    parser.addoption(
        "--headless",
        action="store_true",
        help="Create a virtual display when running on headless servers"
    )

def pytest_configure(config):
    """
    Configure settings before running tests
    """
    #
    # We might want conditional virtual displays. Configure this when
    # appropriate using the "--headless" flag above.
    # For now, let's stay headless always
    #
    os.system('Xvfb :1 -screen 0 1600x1200x16 &')
    os.environ['DISPLAY']=':1.0'

@pytest.fixture
def data_dir():
    """
    Fixture for root of data dir.
    """
    return Path(__file__).parent.parent / 'resources' / 'data'
