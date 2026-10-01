"""Portable configuration paths for source checkouts and installed packages."""

import json
import os
from pathlib import Path
from string import Template

PACKAGE = Path(__file__).resolve().parent


def workspace():
    return Path(os.environ.get('OMICRAFT_WORKSPACE', Path.cwd())).resolve()


def config_path(name, env=None):
    if env and os.environ.get(env):
        return Path(os.environ[env]).expanduser()
    local = PACKAGE / 'configs' / f'{name}.local.json'
    return local if local.is_file() else PACKAGE / 'configs' / f'{name}.example.json'


def load_json(path):
    """Expand configured path placeholders without interpolating JSON syntax."""
    values = {**os.environ, 'OMICRAFT_WORKSPACE': str(workspace()),
              'OMICRAFT_PACKAGE': str(PACKAGE)}

    def expand(value):
        if isinstance(value, dict):
            return {key: expand(item) for key, item in value.items()}
        if isinstance(value, list):
            return [expand(item) for item in value]
        return Template(value).safe_substitute(values) if isinstance(value, str) else value

    return expand(json.loads(Path(path).read_text(encoding='utf-8')))
