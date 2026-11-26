"""The examples are the reference: building each examples/<folder>/<name>.dict gives exactly
the .fr and .md files stored in examples/ (the date of the .fr aside)."""

import re
from pathlib import Path

import pytest

from uvbuilder.cli import main

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'
PROTOMODELS = sorted(EXAMPLES.glob('*/*.dict'))
AUTHORS = ['--author', 'Andre Lessa', '--author', 'Lucas Heck']


def _undated(text):
    return re.sub(r'Date    -> "[^"]*"', 'Date    -> ""', text)


@pytest.mark.parametrize('protomodel', PROTOMODELS, ids=lambda p: p.stem)
def test_example_is_reproduced(protomodel, tmp_path):
    assert main(['build', str(protomodel), '-o', str(tmp_path)] + AUTHORS) == 0
    built = sorted(p.name for p in tmp_path.iterdir())
    stored = sorted(p.name for p in protomodel.parent.iterdir() if p.is_dir())
    assert built == stored
    for name in built:
        for suffix in ('.fr', '.md'):
            new = (tmp_path / name / f'{name}{suffix}').read_text(encoding='utf-8')
            old = (protomodel.parent / name / f'{name}{suffix}').read_text(encoding='utf-8')
            assert _undated(new) == _undated(old), f'{name}{suffix}'
