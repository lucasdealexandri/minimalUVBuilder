"""Formatting helpers of the FeynRules writer (tools/fr_utils.py)."""

import pytest

from uvbuilder.tools.fr_utils import (SECTION_KEYS, formatClassBlock, formatParameterBlock, formatVariable,
                            fr_quantum_numbers, merge_chunks, parameter_block, su2_index_type,
                            toggle_dagger)
from uvbuilder.tools.reps import GlobalU1s, SU2Rep, SU3Rep, U1Rep


def test_parameter_block_options_in_fixed_order_and_aligned():
    text = formatParameterBlock('y', {'Value': 1, 'ComplexParameter': False,
                                      'BlockName': 'BSMYUKAWA', 'ParameterType': 'External'})
    lines = text.splitlines()
    assert lines[0] == '  y == {' and lines[-1] == '  }'
    assert [line.split('->')[0].strip() for line in lines[1:-1]] == [
        'ParameterType', 'BlockName', 'Value', 'ComplexParameter']
    assert len({line.index('->') for line in lines[1:-1]}) == 1


def test_parameter_block_aligns_multi_line_values():
    text = formatParameterBlock('U', {'ParameterType': 'External', 'ComplexParameter': False,
                                      'Value': '{U[1, 1] -> 1,\n U[1, 2] -> 0}'})
    first, second = [line for line in text.splitlines() if 'U[1' in line]
    assert first.index('{') == second.index(' U[1, 2]')


@pytest.mark.parametrize('missing', ['ParameterType', 'ComplexParameter'])
def test_parameter_block_requires_type_and_complexity(missing):
    options = {'ParameterType': 'External', 'ComplexParameter': False}
    del options[missing]
    with pytest.raises(KeyError, match=missing):
        formatParameterBlock('y', options)


def test_parameter_entry_of_a_real_coupling():
    p = parameter_block('M', 'Mass', value=1, complex_parameter=False, block_name='BSMMASS')
    assert (p.name, p.sector, p.kind, p.sort_key) == ('M', 'Mass', 'External', ('M', 0))
    assert p.base == 'M'
    assert dict(p.options) == {'ParameterType': 'External', 'BlockName': 'BSMMASS',
                               'Value': 1, 'ComplexParameter': False}
    assert p.render(OrderBlock=3).count('OrderBlock') == 1


def test_parameter_entries_of_a_complex_coupling_stay_together():
    parts = [parameter_block('y', 'Yukawa', value=0, complex_parameter=False, part='I'),
             parameter_block('y', 'Yukawa', parameter_type='Internal', value='Ry + I*Iy'),
             parameter_block('y', 'Yukawa', value=1, complex_parameter=False, part='R')]
    assert [p.name for p in sorted(parts, key=lambda p: p.sort_key)] == ['y', 'Ry', 'Iy']
    assert [p.base for p in parts] == ['y', 'y', 'y']


def test_merge_chunks_keeps_order_and_drops_repeats():
    merged = merge_chunks({'Parameters': ['a', 'b', 'a']}, {}, {'Parameters': ['b', 'c']})
    assert merged['Parameters'] == ['a', 'b', 'c']
    assert set(SECTION_KEYS) <= set(merged)
    assert all(merged[key] == [] for key in SECTION_KEYS if key != 'Parameters')


def test_su2_index_type():
    assert [su2_index_type(d) for d in (2, 3, 10)] == ['SU2D', 'SU2d3', 'SU2d10']
    with pytest.raises(ValueError):
        su2_index_type(11)


def test_quantum_numbers_omit_vanishing_global_charges():
    assert fr_quantum_numbers({'Q': -1}, GlobalU1s(D=1, B=0)) == '{Q -> -1, DarkCharge -> 1}'
    assert fr_quantum_numbers({'Y': U1Rep(1, 2)}) == '{Y -> 1/2}'


def test_toggle_dagger_is_an_involution():
    assert toggle_dagger('X') == 'Xbar' and toggle_dagger('Xbar') == 'X'


def test_index_types_of_gauge_representations():
    colors = [SU3Rep(1, 0), SU3Rep(1, 1), SU3Rep(2, 0), SU3Rep(0, 1), SU3Rep(0, 2)]
    assert [formatVariable(c) for c in colors] == ['Colour', 'Gluon', 'Sextet', 'SU3d3bar', 'SU3d6bar']
    assert formatVariable(SU2Rep(2)) == 'SU2d3'
    assert formatVariable([SU3Rep(1, 0), 'x']) == '{ Colour, x }'
    with pytest.raises(ValueError):
        formatVariable(SU3Rep(3, 0))


def test_class_block_requires_name_and_skips_unset_options():
    text = formatClassBlock('S[1]', {'ClassName': 'X', 'SelfConjugate': 'False', 'Mass': None})
    assert 'ClassName      -> X' in text and 'Mass' not in text
    with pytest.raises(KeyError, match='ClassName'):
        formatClassBlock('S[1]', {'SelfConjugate': 'False'})
