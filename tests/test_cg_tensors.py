"""Tests for the CG layer: tools/cg_tensors.py.

Covers:
  * the FeynRules index-type mapping (SU2D/SU2W/SU2d{dim}, Colour/Gluon/...),
  * the SU(2) singlet tensor (su2_singlet_coeff) against hand CG values,
  * SU2CGTensor naming/ordering and its FeynRules External definition,
  * SU3CGTensor head resolution from the color reps (K3/K6/T/... and the
    composite J heads), the trivial delta/trace contraction, and the escape
    warning.
"""

import warnings

import sympy as sp

import pytest

from uvbuilder.tools.reps import U1Rep, SU2Rep, SU3Rep
from uvbuilder.tools.cg_tensors import (
    SU2CGTensor, SU3CGTensor,
    cg_block, feyn_index_type, su2_singlet_coeff)


def test_index_type_map():
    assert feyn_index_type(SU2Rep(1,)) == 'SU2D'
    assert feyn_index_type(SU2Rep(2,)) == 'SU2d3'      # never SU2W (Cartesian SM adjoint)
    assert feyn_index_type(SU2Rep(3,)) == 'SU2d4'
    assert feyn_index_type(SU3Rep(1, 0)) == 'Colour'
    assert feyn_index_type(SU3Rep(0, 1)) == 'SU3d3bar'
    assert feyn_index_type(SU3Rep(2, 0)) == 'Sextet'      # declared in SMMod.fr
    assert feyn_index_type(SU3Rep(0, 2)) == 'SU3d6bar'
    assert feyn_index_type(SU3Rep(1, 1)) == 'Gluon'


def test_su2_index_type_limit():
    assert feyn_index_type(SU2Rep(9,)) == 'SU2d10'      # largest declared in SMMod.fr
    with pytest.raises(ValueError):
        feyn_index_type(SU2Rep(10,))                     # dimension 11


def test_su2_singlet_coeff_two_doublets():
    t = su2_singlet_coeff((SU2Rep(1,), SU2Rep(1,)))
    assert t[(U1Rep(1, 2), U1Rep(-1, 2))] == sp.sqrt(2) / 2
    assert t[(U1Rep(-1, 2), U1Rep(1, 2))] == -sp.sqrt(2) / 2
    assert len(t) == 2


def test_su2_singlet_coeff_sum_zero_only():
    t = su2_singlet_coeff((SU2Rep(1,), SU2Rep(1,), SU2Rep(2,)))
    for mtuple in t:
        total = sum(sp.Rational(m.numerator, m.denominator) for m in mtuple)
        assert total == 0
    assert t  # nonzero singlet tensor exists for 2 x 2 x 3


def test_su2_cg_naming_sorted_dims():
    # isospins 4 5 2 2 2 3 (SU2Rep dynkin = dim - 1)
    t = SU2CGTensor([SU2Rep(4,), SU2Rep(5,), SU2Rep(1,), SU2Rep(1,),
                     SU2Rep(1,), SU2Rep(2,)])
    assert t.dims == (2, 2, 2, 3, 5, 6)
    assert t.name == 'SU2CG222356'
    assert t.index_types == ('SU2D', 'SU2D', 'SU2D', 'SU2d3', 'SU2d5', 'SU2d6')


def test_su2_cg_drops_trivial_singlets():
    # isospin dim-1 entries are redundant: never part of the tensor/name.
    t = SU2CGTensor([SU2Rep(0,), SU2Rep(0,)])
    assert t.dims == ()
    assert t.name == ''
    assert t.render(()) == ''
    assert t.fr_definition() == ''
    assert t.numeric == {}
    mixed = SU2CGTensor([SU2Rep(0,), SU2Rep(0,), SU2Rep(1,), SU2Rep(1,)])
    assert mixed.dims == (2, 2)
    assert mixed.name == 'SU2CG22'
    assert 'SU2d1' not in mixed.fr_definition()


def test_su2_cg_numeric_and_def():
    t = SU2CGTensor([SU2Rep(1,), SU2Rep(1,)])
    assert t.name == 'SU2CG22'
    assert cg_block(SU2Rep(1,), SU2Rep(1,)) is not None
    # FR definition: rigid Internal parameter, matrix-style Definitions
    # (explicit non-zero entries, then a numeric-only default of 0)
    d = t.fr_definition()
    assert d.startswith('  SU2CG22 == {')
    assert 'Index[SU2D]' in d
    assert 'ParameterType    -> Internal' in d
    assert 'ComplexParameter -> False' in d
    assert 'SU2CG22[1, 2] -> 1' in d
    assert 'SU2CG22[2, 1] -> -1' in d
    assert 'SU2CG22[i1_?NumericQ, i2_?NumericQ] :> 0' in d
    assert d.index('SU2CG22[2, 1] ->') < d.index('_?NumericQ')   # explicit entries first


def test_su3_head_default_labels():
    assert SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 0), SU3Rep(1, 0)]).render(
        ('ii', 'jj', 'kk')) == 'K3[ii, jj, kk]'
    assert SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 1), SU3Rep(2, 0)]).render(
        ('m', 'ii', 'jj')) == 'K6bar[m, ii, jj]'
    assert SU3CGTensor([SU3Rep(0, 1), SU3Rep(1, 0), SU3Rep(1, 1)]).render(
        ('a', 'ii', 'jj')) == 'T[a, ii, jj]'


def test_su3_head_names_and_dims():
    assert SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 0), SU3Rep(1, 0)]).name == 'K3'
    assert SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 1), SU3Rep(2, 0)]).name == 'K6bar'
    assert SU3CGTensor([SU3Rep(0, 2), SU3Rep(1, 0), SU3Rep(1, 0)]).name == 'K6'
    assert SU3CGTensor([SU3Rep(0, 1), SU3Rep(1, 0), SU3Rep(1, 1)]).name == 'T'
    assert SU3CGTensor([SU3Rep(0, 2), SU3Rep(1, 1), SU3Rep(2, 0)]).name == 'T6'
    assert SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 1), SU3Rep(2, 0)]).name == 'J'
    assert SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 2), SU3Rep(1, 1)]).name == 'Jbar'
    assert SU3CGTensor([SU3Rep(1, 1), SU3Rep(1, 1)]).name == ''


def test_su3_composite_head_render_with_labels():
    # Dummy indices are optional: render auto-fills them from default_dummy
    # when omitted, or uses explicit labels when supplied (externals, dummies).
    j = SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 1), SU3Rep(2, 0)])
    assert j.dummy_count == 3
    explicit = '-I Sqrt[2] K3[i1, x, y] T[i2, x, z] K6bar[i3, z, y]'
    assert j.render(('i1', 'i2', 'i3', 'x', 'y', 'z')) == explicit
    assert j.render(('i1', 'i2', 'i3')) == \
        '-I Sqrt[2] K3[i1, xx, yy] T[i2, xx, zz] K6bar[i3, zz, yy]'
    jb = SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 2), SU3Rep(1, 1)])
    assert jb.render(('i1', 'i2', 'i3', 'x', 'y', 'z')) == \
        'I Sqrt[2] K6[i1, x, y] T[i2, y, z] K3bar[i3, z, x]'


def test_su3_render_too_few_indices_raises():
    import pytest
    j = SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 1), SU3Rep(2, 0)])
    with pytest.raises(ValueError, match='too few indices'):
        j.render(('i1', 'i2'))
    k3 = SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 0), SU3Rep(1, 0)])
    with pytest.raises(ValueError, match='too few indices'):
        k3.render(())
    with pytest.raises(ValueError, match='expected'):
        j.render(('i1', 'i2', 'i3', 'x'))


def test_su3_builtins_need_no_definition():
    k3 = SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 0), SU3Rep(1, 0)])
    assert k3.fr_symbol() == 'K3'
    assert k3.fr_definition() == ''
    k6bar = SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 1), SU3Rep(2, 0)])
    assert k6bar.fr_definition() == ''


def test_su3_trivial_conjugate_pairs():
    t = SU3CGTensor([SU3Rep(1, 0), SU3Rep(0, 1)])
    assert t.is_trivial
    assert not t.is_head
    assert t.name == ''
    assert t.render(()) == ''
    assert t.trivial_pairs == ((0, 1),)
    # multi-pair: each conjugate pair shares one label
    t2 = SU3CGTensor([SU3Rep(0, 1), SU3Rep(0, 1), SU3Rep(1, 0), SU3Rep(1, 0)])
    assert t2.is_trivial
    assert t2.trivial_pairs == ((0, 2), (1, 3))


def test_su3_no_color_singlets_is_trivial():
    t = SU3CGTensor([])
    assert t.is_trivial
    assert t.trivial_pairs == ()
    assert t.render(()) == ''
    # color singlets are dropped entirely
    assert SU3CGTensor([SU3Rep(0, 0), SU3Rep(0, 0)]).is_trivial


def test_su3_escape_warns():
    # 3 x 3 x 8 has no color singlet and no matching head: escape.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        t = SU3CGTensor([SU3Rep(1, 0), SU3Rep(1, 0), SU3Rep(1, 1)])
        assert t.is_notimplemented
        assert any('no color contraction' in str(x.message) for x in w)
    assert t.name == ''
    assert t.render(()) == ''


@pytest.mark.parametrize('dims', [(2, 2), (2, 2, 3), (2, 2, 3, 3), (3, 3), (2, 2, 2, 2)])
def test_fr_definition_indexed_like_multiplet_components(dims):
    # the written tensor, at the component positions used in the class definitions
    # (getSU2Index), equals the numeric singlet; with odd sum of isospins a reversed
    # index order would flip its sign
    import re
    from types import SimpleNamespace
    from uvbuilder.tools.fr_utils import getSU2Index
    from uvbuilder.tools.reps import SU2Rep
    from uvbuilder.tools.cg_tensors import su2_singlet_paths
    isospins = tuple(SU2Rep(d - 1) for d in dims)
    for path in (su2_singlet_paths(isospins) if len(dims) > 3 else [None]):
        t = SU2CGTensor(list(isospins), path=path)
        written = {tuple(map(int, k.split(', '))): v for k, v in
                   re.findall(rf'{t.name}\[([\d, ]+)\] -> ([^,\n}}]+)', t.fr_definition())}
        for ms, value in t.numeric.items():
            position = tuple(getSU2Index(SimpleNamespace(isospin=SU2Rep(d - 1), i3=m))
                             for m, d in zip(ms, t.dims))
            assert sp.sympify(written[position].replace('^', '**')) == value
