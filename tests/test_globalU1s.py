"""Tests for GlobalU1s (tools/reps.py): the B, L and dark global U(1) charges.

Covers tracked-zero vs untracked, projection, arithmetic (and its refusal to
mix different tracked sets), conjugation, is_zero, the total order and repr.
"""

import pytest

from uvbuilder.tools.reps import GlobalU1s, U1Rep


def test_untracked_vs_tracked_zero():
    assert GlobalU1s() != GlobalU1s(L=0)
    assert GlobalU1s().tracked == frozenset()
    assert GlobalU1s(L=0).tracked == frozenset({'L'})
    assert GlobalU1s(L=0)['L'] == 0 and GlobalU1s()['L'] is None


def test_values_become_u1reps():
    q = GlobalU1s(B=U1Rep(1, 3), L=1)
    assert isinstance(q.B, U1Rep) and q.B == U1Rep(1, 3)
    assert isinstance(q.L, U1Rep) and q.L == 1
    assert GlobalU1s(B=0.5).B == U1Rep(1, 2)


def test_hash_and_equality():
    assert GlobalU1s(L=1) == GlobalU1s(L=U1Rep(1))
    assert len({GlobalU1s(L=1), GlobalU1s(L=U1Rep(1)), GlobalU1s(L=0), GlobalU1s()}) == 3


def test_project():
    sm_quark = GlobalU1s(B=U1Rep(1, 3))
    assert sm_quark.project({'L'}) == GlobalU1s(L=0)   # B dropped, L tracked and zero
    assert sm_quark.project(set()) == GlobalU1s()
    # a charge the state does not carry becomes a tracked zero
    assert sm_quark.project({'B', 'L'}) == GlobalU1s(B=U1Rep(1, 3), L=0)
    assert GlobalU1s().project({'L'}) == GlobalU1s(L=0)
    assert GlobalU1s(B=U1Rep(1, 3), L=0).project({'B'}) == GlobalU1s(B=U1Rep(1, 3))


def test_arithmetic_and_conjugation():
    a, b = GlobalU1s(L=1, D=1), GlobalU1s(L=-1, D=0)
    assert a + b == GlobalU1s(L=0, D=1)
    assert a - b == GlobalU1s(L=2, D=1)
    assert -a == GlobalU1s(L=-1, D=-1) == a.conjRep
    assert GlobalU1s() + GlobalU1s() == GlobalU1s()


def test_mixing_tracked_sets_raises():
    with pytest.raises(ValueError, match='different tracked sets'):
        GlobalU1s(L=1) + GlobalU1s()
    with pytest.raises(ValueError):
        GlobalU1s(B=U1Rep(1, 3)) - GlobalU1s(L=0)


def test_is_zero():
    assert GlobalU1s().is_zero                 # nothing tracked
    assert GlobalU1s(L=0, D=0).is_zero
    assert not GlobalU1s(L=0, D=1).is_zero


def test_total_order():
    reps = [GlobalU1s(), GlobalU1s(D=1), GlobalU1s(L=1), GlobalU1s(L=0), GlobalU1s(B=U1Rep(-1, 3))]
    ordered = sorted(reps)
    assert sorted(reversed(reps)) == ordered          # deterministic
    for x, y in zip(ordered, ordered[1:]):
        assert x < y and y > x and x <= y and not (y < x)


def test_repr():
    assert repr(GlobalU1s()) == '()'
    assert repr(GlobalU1s(B=U1Rep(1, 3), L=0)) == '(B=1/3, L=0)'


# -- operators: SM charges are projected onto the charges the BSM fields track --

def _bsm(label, spin, iso, y, q):
    from uvbuilder.states.multiplet import Multiplet
    from uvbuilder.tools.reps import BrokenRep, SU2Rep, SU3Rep
    rep = BrokenRep(spin=spin, colordim=SU3Rep(0, 0), isospin=SU2Rep(iso,),
                    hypercharge=y, globalU1s=q)
    return Multiplet(rep, label=label, fillMultiplet=True)


def _lepton_portal_is_singlet(q_scalar, q_weyl):
    '''S (scalar singlet) x Du (Weyl doublet, Y=1/2) x lL1 (SM lepton doublet, L=1).'''
    from uvbuilder.tools.operators import Operator
    from uvbuilder.tools.reps import SCALAR, WEYL_LEFT
    from uvbuilder.states.SMfields import lL1
    S = _bsm('S', SCALAR, 0, U1Rep(0), q_scalar)
    Du = _bsm('Du', WEYL_LEFT, 1, U1Rep(1, 2), q_weyl)
    return any(rep.is_singlet for rep in Operator(S, Du, lL1).possible_reps)


def test_operator_ignores_untracked_lepton_number():
    # L not conserved: the lepton's L = 1 is projected away
    assert _lepton_portal_is_singlet(GlobalU1s(), GlobalU1s())


def test_operator_enforces_conserved_lepton_number():
    assert not _lepton_portal_is_singlet(GlobalU1s(L=0), GlobalU1s(L=0))   # L: 0 + 0 + 1
    assert _lepton_portal_is_singlet(GlobalU1s(L=0), GlobalU1s(L=-1))      # L: 0 - 1 + 1


# -- FeynRules: conserved non-zero charges become QuantumNumbers --

def test_fr_quantum_numbers_helper():
    from uvbuilder.tools.fr_utils import fr_quantum_numbers
    assert fr_quantum_numbers({'Q': U1Rep(-1)}, GlobalU1s(L=1)) == '{Q -> -1, LeptonNumber -> 1}'
    assert fr_quantum_numbers({'Y': U1Rep(2, 3)}, GlobalU1s(B=U1Rep(1, 3), D=0)) \
        == '{Y -> 2/3, BaryonNumber -> 1/3}'                       # zero charges omitted
    assert fr_quantum_numbers({'Q': U1Rep(0)}, GlobalU1s(D=-1)) == '{Q -> 0, DarkCharge -> -1}'
    assert fr_quantum_numbers({'Q': U1Rep(1)}, GlobalU1s()) == '{Q -> 1}'


def test_stop_model_emits_baryon_number():
    import warnings
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    stop = {"masses": {1000022: 180.603, 1000006: 1164.264, 1000021: 1500.0},
            "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}, 1000021: {(1000006, -6): 1.0}},
            "labels": {1000022: 'B', 1000006: 'st', 1000021: 'Go'}}
    builder = ModelBuilder(modelData=stop, possible_dm_reps={'globalU1s': [GlobalU1s(B=0)]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        valid = [(pm, mm) for pm in builder.generateParticleModels()
                 for mm in PartitionsBuilder(pm).get_multipletModels() if mm.is_valid]
    assert valid
    pm, mm = valid[0]
    fr = FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr').to_feynrules()
    assert 'QuantumNumbers -> {Q -> 2/3, BaryonNumber -> 1/3}' in fr     # physical stop class
    assert 'BaryonNumber -> -1/3' in fr                                   # its conjugate gauge partner
    assert 'LeptonNumber' not in fr and 'DarkCharge' not in fr      # not conserved here
