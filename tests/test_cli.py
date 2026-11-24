"""Command line interface (cli.py)."""

import textwrap

from uvbuilder.cli import main
from uvbuilder.models.feynrulesmodel import FeynRulesModel
from uvbuilder.paramCard import ParamCardGenerator
from uvbuilder.models.slhaBuilder import write_card
from tests.test_slhaBuilder import MAJORANA, SDDM_POINT, _sddm_model, _template

PROTOMODEL = '''
    from uvbuilder.tools.reps import DIRAC_SPINOR, GlobalU1s
    protomodel = {"masses": {1000022: 100., 1000023: 300., 1000024: 200., 1000025: 400.},
                  "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0},
                             1000024: {(1000022, 24): 1.0}, 1000025: {(1000023, 25): 1.0}},
                  "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}
    possible_dm_reps = {'spin': [DIRAC_SPINOR], 'globalU1s': [GlobalU1s()]}
'''


def test_build_writes_the_fr_and_pickle_of_each_model(tmp_path, capsys):
    source = tmp_path / 'sddm.py'
    source.write_text(textwrap.dedent(PROTOMODEL))
    assert main(['build', str(source), '--out', str(tmp_path / 'out')]) == 0
    assert sorted(p.name for p in (tmp_path / 'out').iterdir()) == ['sddm.fr', 'sddm.pkl']
    assert 'M$ModelName = "sddm";' in (tmp_path / 'out' / 'sddm.fr').read_text()
    assert '[ChiI1Y0Lw, ChiI2Y3Lw, ChiI2Ym3Lw]' in capsys.readouterr().out


def test_fill_from_the_pickle_of_build(tmp_path):
    source = tmp_path / 'sddm.py'
    source.write_text(textwrap.dedent(PROTOMODEL))
    main(['build', str(source), '--out', str(tmp_path)])
    generator = ParamCardGenerator(FeynRulesModel.load(tmp_path / 'sddm.pkl'))
    template = _template(generator.fr_model, tmp_path / 'template.dat')
    reference = write_card(_sddm_model(MAJORANA), template, tmp_path / 'reference.dat',
                           SDDM_POINT).read_text()
    card = generator.fill(template, SDDM_POINT, tmp_path / 'card.dat')     # inputs set
    assert main(['fill', str(tmp_path / 'sddm.pkl'), str(card)]) == 0
    assert card.read_text() == reference


def test_build_reads_a_dict_file(tmp_path):
    (tmp_path / 'sddm-2.dict').write_text(repr(
        {'masses': {1000022: 100., 1000023: 300., 1000024: 200., 1000025: 400.},
         'decays': {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                    1000025: {(1000023, 25): 1.0}},
         'xsecs[fb]': {((1000022, 1000023), 13.0): 1.0}, 'walkerid': 'uphi'}))
    assert main(['build', str(tmp_path / 'sddm-2.dict'), '--out', str(tmp_path / 'out')]) == 0
    assert (tmp_path / 'out' / 'sddm_2_1.fr').is_file() and (tmp_path / 'out' / 'sddm_2_1.pkl').is_file()


def test_build_needs_a_python_file_with_the_protomodel(tmp_path, capsys):
    (tmp_path / 'model.slha').write_text('BLOCK MASS\n')
    (tmp_path / 'empty.py').write_text('x = 1\n')
    assert main(['build', str(tmp_path / 'model.slha')]) == 1
    assert main(['build', str(tmp_path / 'empty.py')]) == 1
    assert 'does not define the dictionary `protomodel`' in capsys.readouterr().err
    (tmp_path / 'code.dict').write_text("__import__('os').system('echo hi')")
    assert main(['build', str(tmp_path / 'code.dict')]) == 1
    assert 'is not a dictionary literal' in capsys.readouterr().err


def test_fill_rejects_a_pickle_of_another_object(tmp_path, capsys):
    import pickle
    (tmp_path / 'other.pkl').write_bytes(pickle.dumps({'masses': {}}))
    assert main(['fill', str(tmp_path / 'other.pkl'), str(tmp_path / 'card.dat')]) == 1
    assert 'holds a dict, not a FeynRulesModel' in capsys.readouterr().err
