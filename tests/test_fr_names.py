"""Tests for tools/fr_names.py: the static FeynRules model-file parser and the
reserved-name registry.

Real files: tests/data/SM.fr and the package SMMod.fr (uvbuilder/data). Synthetic snippets cover the
scanner (comments, strings, statement splitting) and the value forms FeynRules
allows (Mass as number / symbol / {sym, value} / {family, {m1, v1}, ...}).
"""

from fractions import Fraction
from pathlib import Path

import pytest

from uvbuilder.tools.fr_names import (FRModelFile, ReservedNames, is_valid_symbol, language_lists,
                            split_statements, strip_comments, symbols_in)

SM = Path(__file__).resolve().parent / 'data' / 'SM.fr'
from uvbuilder.tools.fr_names import DEFAULT_SM_FILE as SMMOD


@pytest.fixture(scope='module')
def smmod():
    return FRModelFile.parse(SMMOD)


@pytest.fixture(scope='module')
def reserved():
    return ReservedNames.from_fr(SMMOD)


# --- scanner -------------------------------------------------------------------------------

def test_strip_comments_nested_keeps_lines_and_strings():
    text = 'a = 1; (* x (* nested *) y *)\nb = "(* not a comment *)";'
    clean = strip_comments(text)
    assert len(clean) == len(text) and clean.count('\n') == 1
    assert 'nested' not in clean and '(* not a comment *)' in clean


def test_split_statements_semicolons_and_newlines():
    text = 'x = f[a,\n  b];\ny := g[c]\nz = 1 +\n  2\nw = {1, 2}'
    stmts = [s for _, s in split_statements(text)]
    assert stmts == ['x = f[a,\n  b]', 'y := g[c]', 'z = 1 +\n  2', 'w = {1, 2}']


def test_symbols_in_ignores_strings_patterns_and_named_characters():
    assert symbols_in('{MZ, 91.1876}') == ['MZ']
    assert symbols_in('{Ml, {Me,5.11*^-4}, {MMU,0.10566}}') == ['Ml', 'Me', 'MMU']
    assert symbols_in('Subscript[\\[Alpha], s] "e-" sp_ Index[Colour]') == ['Subscript', 's', 'Index', 'Colour']


@pytest.mark.parametrize('name, valid', [('Chi1', True), ('C+', False), ('1000006Lw', False),
                                         ('x_1', False), ('M$x', False), ('', False)])
def test_is_valid_symbol(name, valid):
    assert is_valid_symbol(name) is valid


def test_mass_forms():
    text = '''M$ClassesDescription = {
      S[1] == {ClassName -> a, SelfConjugate -> True, Mass -> 0, Width -> 0},
      S[2] == {ClassName -> b, SelfConjugate -> True, Mass -> Mb, Width -> Wb},
      S[3] == {ClassName -> c, SelfConjugate -> False, Mass -> {Mc, Internal}, Width -> {Wc, 1.}},
      F[4] == {ClassName -> d, ClassMembers -> {d1, d2}, Mass -> {Md, {Md1, 1}, {Md2, 2}}, Width -> {0, {Wd2, 3}}}
    };'''
    m = FRModelFile.from_text(text)
    assert [m.classes[n].mass_symbols for n in 'abcd'] == [(), ('Mb',), ('Mc',), ('Md', 'Md1', 'Md2')]
    assert m.classes['d'].width_symbols == ('Wd2',)
    assert m.classes['c'].conjugate_names == ('cbar',)
    assert m.classes['d'].conjugate_names == ('dbar', 'd1bar', 'd2bar')   # SelfConjugate missing -> False
    assert m.classes['c'].line == 4


# --- SMMod.fr ------------------------------------------------------------------------------

def test_files_parse_completely():
    for path in (SM, SMMOD):
        m = FRModelFile.parse(path)
        assert m.skipped == [], (path.name, m.skipped)
        assert set(m.gauge_groups) == {'U1Y', 'SU2L', 'SU3C'}
        assert 'LSM' in m.definitions


def test_smmod_weyl_classes(smmod):
    LLw, QLw, eRw, uRw, dRw = (smmod.classes[n] for n in ('LLw', 'QLw', 'eRw', 'uRw', 'dRw'))
    assert (LLw.kind, LLw.chirality, LLw.indices) == ('W', 'Left', ('SU2D', 'Generation'))
    assert QLw.indices == ('SU2D', 'Generation', 'Colour')
    assert QLw.quantum_numbers == {'Y': Fraction(1, 6), 'BaryonNumber': Fraction(1, 3)}
    assert (eRw.chirality, eRw.quantum_numbers['Y']) == ('Right', -1)
    assert uRw.indices == ('Generation', 'Colour') and dRw.quantum_numbers['Y'] == Fraction(-1, 3)
    assert eRw.conjugate_names == ('eRwbar',)


def test_smmod_dirac_classes(smmod):
    dq = smmod.classes['dq']
    assert dq.weyl_components == ('dLw', 'dRw')
    assert dq.members == ('d', 's', 'b') and dq.pdg == (1, 3, 5)
    assert dq.mass_symbols == ('Md', 'MD', 'MS', 'MB')
    assert dq.quantum_numbers == {'Q': Fraction(-1, 3), 'BaryonNumber': Fraction(1, 3)}
    assert smmod.classes['l'].particle_names == ('e-', 'mu-', 'ta-')


def test_find_classes_by_quantum_numbers(smmod):
    '''SM multiplets are matched to SMMod.fr classes by quantum numbers.'''
    assert [c.name for c in smmod.find_classes(kind='W', Y=Fraction(1, 6))] == ['QLw']
    assert [c.name for c in smmod.find_classes(kind='W', chirality='Right', Y=-1)] == ['eRw']
    assert {c.name for c in smmod.find_classes(kind='W', chirality='Right', BaryonNumber=Fraction(1, 3))} \
        == {'uRw', 'dRw'}
    assert [c.name for c in smmod.find_classes(kind='S', unphysical=True, Y=Fraction(1, 2))] == ['Phi']
    assert [c.name for c in smmod.find_classes(kind='S', Y=Fraction(-1, 2))] == ['Phitilde']


def test_smmod_indices_are_sm_definitions(smmod):
    '''The SM index declarations a generated file must not override.'''
    su2d, gen = smmod.indices['SU2D'], smmod.indices['Generation']
    assert (su2d.unfold, su2d.size, su2d.style) == (True, 2, 'k')
    assert (gen.unfold, gen.size, gen.style) == (False, 3, 'f')
    assert smmod.indices['Colour'].range == 'NoUnfold[Range[3]]'


def test_parameters(smmod):
    assert smmod.parameters['aEWM1'].external and smmod.parameters['aEWM1'].block == 'SMINPUTS'
    assert smmod.parameters['gs'].aliases == ('G',)
    assert smmod.parameters['yl'].aliases == ('ye', 'ym', 'ytau')
    assert smmod.parameters['CKM'].indices == ('Generation', 'Generation')


def test_gauge_group_symbols(smmod):
    su2 = smmod.gauge_groups['SU2L'].symbols
    assert su2[:5] == ('gw', 'Wi', 'Eps', 'Ta', 'SU2D') and su2[-1] == 'FSU2L'
    assert {f'TSU2d{d}' for d in range(3, 11)} | {f'SU2d{d}' for d in range(3, 11)} <= set(su2)
    assert set(smmod.gauge_groups['SU3C'].symbols) == {'gs', 'G', 'f', 'T', 'Colour', 'T6', 'Sextet', 'dSUN'}


# --- ReservedNames ---------------------------------------------------------------------------

def test_ms_incident_is_caught(reserved):
    '''A generated parameter MS (or ms) would redefine the strange-quark mass.'''
    (origin,) = reserved.clashes('MS')
    assert (origin.kind, origin.where.split(':')[0]) == ('mass', 'SMMod.fr')
    assert [o.symbol for o in reserved.clashes('ms')] == ['MS']     # MG5 is case-insensitive
    assert not reserved.is_free('Ms')


def test_case_rules(reserved):
    # FR-only names (Weyl classes) are exact: a differently-cased name is free ...
    assert reserved.clashes('eRw') and not reserved.clashes('ERw')
    # ... unless the new symbol itself goes to the UFO
    assert [o.symbol for o in reserved.clashes('ERw', ufo_bound=True)] == ['eRw']
    # language lists are always exact, even for UFO-bound symbols
    assert [o.where for o in reserved.clashes('lambda')] == ['python keywords']
    assert all(o.kind != 'language' for o in reserved.clashes('Lambda', ufo_bound=True))


def test_automatic_names_are_reserved(reserved):
    for name in ('uRwbar', 'ebar', 'Phitilde', 'LSM', 'LSMWeyl', 'Generation', 'ye', 'SMINPUTS', 'FSU2L'):
        assert reserved.clashes(name), name
    assert reserved.clashes('ta-')     # particle-name strings


def test_origin_points_to_line(reserved):
    (origin,) = [o for o in reserved.origins('eRw') if o.kind == 'class']
    line = int(origin.where.split(':')[1])
    lines = SMMOD.read_text().splitlines()                  # line is 1-based: the 'W[n] == {' line
    assert lines[line - 1].strip().startswith('W[12] ==') and 'eRw' in lines[line]


def test_language_lists_fallback_and_python(tmp_path):
    lists = language_lists(tmp_path)            # empty dir -> fallbacks
    assert 'mathematica (fallback)' in lists and 'feynrules (fallback)' in lists
    assert 'lambda' in lists['python keywords'] and 'abs' in lists['python builtins']
    (tmp_path / 'mathematica_14.1.txt').write_text('# Mathematica 14.1\nSystem`Pi\nFoo\n')
    lists = language_lists(tmp_path)
    assert lists['mathematica_14.1.txt'] == ['Pi', 'Foo'] and 'mathematica (fallback)' not in lists


# --- renaming ------------------------------------------------------------------------------

from uvbuilder.tools.fr_names import (Renamer, fix_symbol, particle_symbols, parameter_symbols,
                            reserved_names, warn_labels)


@pytest.mark.parametrize('name, fixed', [('C+', 'CP'), ('C-', 'CM'), ('1000022', 'X1000022'),
                                         ('st~', 'st'), ('N1', 'N1'), ('', 'X')])
def test_fix_symbol(name, fixed):
    assert fix_symbol(name) == fixed


def test_renamer_suffixes_and_generated_uniqueness():
    r = Renamer(reserved_names(SMMOD))
    assert r.claim('N1', particle_symbols, 'test') == 'N1'           # free: unchanged
    assert r.claim('B', particle_symbols, 'test') == 'B1'            # SM class B, mass MB
    assert r.claim('C+', particle_symbols, 'test') == 'CP'           # invalid: fixed
    assert r.claim('MN1', parameter_symbols, 'test') == 'MN11'       # taken by N1's mass
    assert r.claim('b1', particle_symbols, 'test') == 'b11'          # MB1 taken (case-insensitive); suffix appended


def test_renamer_does_not_touch_the_shared_registry():
    shared = reserved_names(SMMOD)
    Renamer(shared).claim('Zp', particle_symbols, 'test')
    assert not shared.clashes('Zp')


def test_warn_labels_only_warns():
    with pytest.warns(UserWarning, match=r"'B'.*will be renamed|will be renamed.*'B'"):
        bad = warn_labels(['N1', 'B'], particle_symbols, 'Particle', sm_file=SMMOD)
    assert [b.split(':')[0] for b in bad] == ["'B'"]


def test_smmod_declares_su3_indices(smmod):
    assert {'Colour', 'Sextet', 'Gluon'} <= set(smmod.indices)
    assert smmod.indices['Sextet'].size == 6 and smmod.indices['Sextet'].unfold is False
    assert 'T6' in smmod.gauge_groups['SU3C'].symbols


def test_smmod_declares_bsm_su2_indices(smmod):
    for d in range(3, 11):
        idx = smmod.indices[f'SU2d{d}']
        assert (idx.size, idx.unfold) == (d, True), d
    assert 'SU2d11' not in smmod.indices
