from pathlib import Path
from ingrid.utils.helper import resolve_path

def test_utility_resolve_path() -> None:
    """
    Test the resolve_path function.
    """
    to_resolve: Path = Path("~")
    resolved: Path = resolve_path(to_resolve)
    assert resolved.exists()
