import os
from pathlib import Path


def local_paths(value, root):
    """Resolve file references without dereferencing virtualenv interpreter symlinks."""
    if isinstance(value, list):
        return [local_paths(item, root) for item in value]
    if not isinstance(value, dict):
        return value
    output = {}
    for key, item in value.items():
        is_path = key.endswith(("_path", "_dir")) or key == "target_structure"
        if isinstance(item, str) and item and (is_path or key == "executable" and "/" in item):
            path = Path(item).expanduser()
            output[key] = os.path.abspath(root / path) if not path.is_absolute() else str(path)
        elif key == "protenix_msa_paths" and isinstance(item, dict):
            output[key] = {
                chain: os.path.abspath(root / Path(path).expanduser())
                for chain, path in item.items()
            }
        elif key == "required_files" and isinstance(item, list):
            output[key] = [os.path.abspath(root / Path(path).expanduser()) for path in item]
        else:
            output[key] = local_paths(item, root)
    return output
