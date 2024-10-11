from pathlib import Path
from typing import Union

def resolve_path(path: Union[Path, str], as_str: bool = False) -> Union[Path, str]:
    """
    Resolve a path to an absolute path, taking into account the "~" that can occur.
    """
    expanded_path = Path(path).expanduser()
    resolved_path = expanded_path.resolve()
    return resolved_path if not as_str else str(resolved_path)
