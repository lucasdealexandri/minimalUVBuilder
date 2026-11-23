"""Regression test: the generated .fr of reference SDDM and stop models must match
the stored files in tests/golden (validated with FeynRules and MG5 check gauge).

After an intended change of the output, regenerate and review with git diff:
    UPDATE_GOLDEN=1 pytest tests/test_golden.py
"""

import difflib
import os
import re
import warnings
from pathlib import Path

import pytest

from tests import test_feynrules
from tests.test_conjugation import _sddm_fr

GOLDEN = Path(__file__).resolve().parent / 'golden'
DATE = re.compile(r'^\s*Date\s*->.*$', re.M)

SDDM = {
    'sddm_doublet_partner': {'PhiI1Y0', 'PhiI2Y3'},
    'sddm_real_triplet': {'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'},
    'sddm_fermion_majorana': {'ChiI1Y0Lw', 'ChiI2Y3Lw', 'ChiI2Ym3Lw'},
    'sddm_fermion_dirac': {'ChiI1Y0D1F1Lw', 'ChiI1Y0D1F2Lw', 'ChiI1Y0Dm1F1Lw',
                           'ChiI1Y0Dm1F2Lw', 'ChiI2Y3D1Lw', 'ChiI2Ym3Dm1Lw'},
}
STOP = {
    'stop_weylstop_scalargluino': {'ChiC3I1Y4B1Lw', 'ChiC3bI1Ym4Bm1Lw', 'PhiC8I1Y0', 'PhiI1Y0'},
    'stop_diracstop_scalargluino': {'ChiC3I1Y4B1D1Lw', 'ChiC3bI1Ym4Bm1Dm1Lw', 'PhiC8I1Y0D1', 'PhiI1Y0D1'},
    'stop_scalarstop_majorana': {'ChiC8I1Y0Lw', 'ChiI1Y0Lw', 'PhiC3I1Y4B1'},
    'stop_scalarstop_dirac': {'ChiC8I1Y0D1Lw', 'ChiC8I1Y0Dm1Lw', 'ChiI1Y0D1Lw', 'ChiI1Y0Dm1Lw',
                              'PhiC3I1Y4B1D1'},
}


def _stop_model(names):
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.tools.reps import GlobalU1s
    builder = ModelBuilder(modelData=test_feynrules.TestRenames.STOP,
                           possible_dm_reps={'globalU1s': [GlobalU1s(B=0), GlobalU1s(B=0, D=1)]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for pm in builder.generateParticleModels():
            for mm in PartitionsBuilder(pm).get_multipletModels():
                if mm.is_valid and {m.fr_symbol() for m in mm.non_conjugated_multiplets} == names:
                    return FeynRulesModel(mm, pm, mm.possible_mixings[0], 'stop.fr')
    raise AssertionError(f'no stop model with multiplets {names}')


def _stop_fr(names):
    return _stop_model(names).to_feynrules()


@pytest.mark.parametrize('name', [*SDDM, *STOP])
def test_golden_fr_outputs(name):
    fr = _sddm_fr(SDDM[name]) if name in SDDM else _stop_fr(STOP[name])
    path = GOLDEN / f'{name}.fr'
    if os.environ.get('UPDATE_GOLDEN'):
        GOLDEN.mkdir(exist_ok=True)
        path.write_text(fr)
    expected, actual = (DATE.sub('', text).splitlines() for text in (path.read_text(), fr))
    if actual != expected:
        diff = list(difflib.unified_diff(expected, actual, str(path), 'generated', lineterm='', n=1))
        pytest.fail('generated .fr differs from the golden file:\n' + '\n'.join(diff[:40]))
