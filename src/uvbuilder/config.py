'''Paths of the external programs. Set once in the user config file
(~/.config/uvbuilder/config.toml, written by `uvbuilder config`); an environment
variable takes precedence, e.g. to try another MG5 for one run.'''

import json
import os
import tomllib
from pathlib import Path
from typing import Optional, Tuple

VARIABLES = {'feynrules': 'UVB_FEYNRULES',          # FeynRules folder
             'wolframscript': 'UVB_WOLFRAMSCRIPT',  # wolframscript executable
             'mg5': 'UVB_MG5'}                      # mg5_aMC executable


def config_file() -> Path:
    base = os.environ.get('XDG_CONFIG_HOME') or Path.home() / '.config'
    return Path(base) / 'uvbuilder' / 'config.toml'


def _read() -> dict:
    path = config_file()
    if not path.is_file():
        return {}
    try:
        with open(path, 'rb') as f:
            return tomllib.load(f)
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f'{path}: {error}') from None


def lookup(name: str) -> Tuple[Optional[str], str]:
    '''(path or None, where it comes from) of an external program (see VARIABLES).'''
    variable = VARIABLES[name]
    if os.environ.get(variable):
        return os.environ[variable], f'environment variable {variable}'
    value = _read().get(name)
    if value:
        return value, str(config_file())
    return None, 'not set'


def tool_path(name: str) -> Path:
    '''Path of an external program; raises if it is not set or does not exist.'''
    value, source = lookup(name)
    if value is None:
        raise RuntimeError(f'{name} is not configured: run `uvbuilder config {name} <path>`')
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f'{name} = {value} ({source}): no such file or folder')
    return path


def set_path(name: str, value) -> Path:
    '''Write the absolute path of an external program into the config file.'''
    if name not in VARIABLES:
        raise ValueError(f'unknown program {name!r}: one of {", ".join(VARIABLES)}')
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f'{path}: no such file or folder')
    values = _read()
    values[name] = str(path)
    file = config_file()
    file.parent.mkdir(parents=True, exist_ok=True)
    file.write_text(''.join(f'{key} = {json.dumps(val)}\n' for key, val in values.items()))
    return file
