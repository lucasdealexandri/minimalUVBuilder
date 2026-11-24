"""SU(2) conjugation convention (tools/reps.py: su2_conj_phase, su2_reality_phase) and its
use: multiplet phases, operators, tilde classes and real multiplets.

The generators below use the same formula as SU2Gen in uvbuilder/data/SMMod.fr (charge
basis, index 1 = highest i3, Condon-Shortley phases), so these tests tie the Python phases
to the matrices FeynRules uses. References in the convention block of tools/reps.py
(Dreiner, Haber, Martin, Phys. Rept. 494 (2010) 1, App. E.1 and eq. (4.3.12)).
"""

import re
import warnings

import pytest
import sympy as sp

from uvbuilder.tools.reps import (BrokenRep, SCALAR, WEYL_LEFT, SU2Rep, SU3Rep, U1Rep,
                        su2_conj_phase, su2_reality_phase)
from uvbuilder.states.multiplet import Multiplet
from uvbuilder.states.SMfields import SM_MULTIPLETS, HiggsDoublet
from uvbuilder.tools.operators import Operator
from uvbuilder.models.feynrulesmodel import IndexTracker

DIMS = range(2, 11)


def _ms(d):
    j = sp.Rational(d - 1, 2)
    return [j - k for k in range(d)]                  # index 1 = highest m


def generators(d):
    '''(J1, J2, J3) exactly as SU2Gen in SMMod.fr.'''
    j, ms = sp.Rational(d - 1, 2), _ms(d)
    Jp = sp.zeros(d)
    for col, m2 in enumerate(ms):
        for row, m1 in enumerate(ms):
            if m1 == m2 + 1:
                Jp[row, col] = sp.sqrt(j * (j + 1) - m2 * (m2 + 1))
    Jm = Jp.T
    return [(Jp + Jm) / 2, (Jp - Jm) / (2 * sp.I), sp.diag(*ms)]


def phase_matrix(d, phase):
    '''Matrix M with (M psi*)_{-m} = phase(m) psi*_m.'''
    ms, M = _ms(d), sp.zeros(d)
    for col, m in enumerate(ms):
        M[ms.index(-m), col] = phase(SU2Rep(d - 1), U1Rep(m.p, m.q))
    return M


@pytest.mark.parametrize('d', DIMS)
def test_generators_are_su2(d):
    J1, J2, J3 = generators(d)
    assert sp.simplify(J1 * J2 - J2 * J1 - sp.I * J3) == sp.zeros(d)
    assert all(sp.simplify(J - J.H) == sp.zeros(d) for J in (J1, J2, J3))   # hermitian


@pytest.mark.parametrize('d', DIMS)
def test_conjugation_matrix(d):
    '''psi^c = C psi* transforms with the same generators: C J* C^-1 = -J
    [Dreiner-Haber-Martin (E.1.4)]; C symmetric for integer j, antisymmetric for
    half-integer j [(E.1.5), (E.1.6)].'''
    C = phase_matrix(d, su2_conj_phase)
    for J in generators(d):
        assert sp.simplify(C * J.conjugate() * C.inv() + J) == sp.zeros(d)
    assert C.T == (C if d % 2 else -C)


def test_doublet_conjugation_is_i_tau2():
    '''j = 1/2: C = i tau^2 = epsilon (Phitilde = i tau^2 Phi*).'''
    assert phase_matrix(2, su2_conj_phase) == sp.Matrix([[0, 1], [-1, 0]])


@pytest.mark.parametrize('d', [3, 5, 7, 9])
def test_reality_condition(d):
    '''psi_{-m} = (-1)^m psi_m*: gauge covariant (R J* R^-1 = -J), an involution
    (R R* = 1) and the neutral component is real (R_00 = +1).'''
    R = phase_matrix(d, su2_reality_phase)
    for J in generators(d):
        assert sp.simplify(R * J.conjugate() * R.inv() + J) == sp.zeros(d)
    assert R * R.conjugate() == sp.eye(d)
    assert R[d // 2, d // 2] == 1
    assert su2_reality_phase(SU2Rep(2,), U1Rep(1)) == -1      # triplet: T^- = -(T^+)*


def test_reality_needs_integer_isospin():
    with pytest.raises(ValueError):
        su2_reality_phase(SU2Rep(1,), U1Rep(1, 2))


def test_labels_carry_no_sign_and_phases_live_on_the_multiplet():
    for m in SM_MULTIPLETS:
        for mm in (m, m.conjugateMultiplet):
            assert not any(f.label.startswith('-') for f in mm.fields), mm.label
    Hc = HiggsDoublet.conjugateMultiplet
    assert Hc.phase(U1Rep(1, 2)) == 1 and Hc.phase(U1Rep(-1, 2)) == -1   # H~ = (Phi0*, -Phi-)
    assert HiggsDoublet.phase(U1Rep(1, 2)) == 1
    assert Hc.conjugateMultiplet is HiggsDoublet


def test_operator_phases_give_phi_dagger_phi():
    '''Phi~^T eps Phi-type contraction of Phi with its conjugate is Phi^dagger Phi: both
    components with the same sign (the CG antisymmetry is undone by the i tau^2 phases).'''
    op = Operator(HiggsDoublet, HiggsDoublet.conjugateMultiplet)
    coeffs = {frozenset(f.label for f in t.fields): t.coeff for t in op.expand_singlets()}
    assert len(coeffs) == 2
    a, b = coeffs.values()
    assert sp.simplify(a - b) == 0 and a != 0


def _scalar(label, dim, y):
    return Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0), isospin=SU2Rep(dim - 1,),
                               hypercharge=U1Rep(y)), label=label, fillMultiplet=True)


def test_tilde_class_definitions():
    from uvbuilder.models.multipletModel import tilde_class
    X = _scalar('X', 2, U1Rep(1, 2))
    block = tilde_class(X, 'Xtilde', IndexTracker())
    assert 'Xtilde[1] :> Xbar[2]' in block and 'Xtilde[2] :> -Xbar[1]' in block
    assert 'Y -> -1/2' in block
    T = _scalar('T', 3, U1Rep(1))
    block = tilde_class(T, 'Ttilde', IndexTracker())
    # (T^c)_{-m} = (-1)^(1+m) T_m*: m = +1 -> +, 0 -> -, -1 -> +
    assert 'Ttilde[3] :> Tbar[1]' in block and 'Ttilde[2] :> -Tbar[2]' in block \
        and 'Ttilde[1] :> Tbar[3]' in block


def test_operator_uses_tilde_symbol_for_conjugate_partner():
    X = _scalar('X', 2, U1Rep(1, 2))
    Xc = Multiplet(X.rep.conjRep, label='Xbar', fillMultiplet=True)   # independent partner multiplet
    t = IndexTracker()
    t.conjugate_symbols = {Xc: 'Xtilde'}
    term = Operator(Xc, HiggsDoublet).fr_term(t)
    assert 'Xtilde[' in term and 'Xbar[' not in term


def test_real_triplet_is_complex_class_with_phased_components():
    '''SDDM builds a real scalar triplet (Y = 0): complex FeynRules class, invariant kinetic
    term, and component -1 defined as minus the conjugate of component +1.'''
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    sddm = {"masses": {1000022: 89.88887, 1000023: 126.046269, 1000024: 155.494055, 1000025: 451.886204},
            "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                       1000025: {(1000023, 25): 1.0}},
            "labels": {1000022: 'N1', 1000023: 'N2', 1000024: 'C1', 1000025: 'N3'}}
    found = 0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for pm in ModelBuilder(modelData=sddm).generateParticleModels():
            for mm in PartitionsBuilder(pm).get_multipletModels():
                reals = [m for m in mm.non_conjugated_multiplets if m.is_real_nontrivial]
                if not reals or not mm.is_valid:
                    continue
                fr = FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr').to_feynrules()
                for m in reals:
                    s = m.fr_symbol()
                    block = re.search(rf'ClassName\s*->\s*{s},.*?\n\s*\}}', fr, re.S).group()
                    assert 'SelfConjugate  -> False' in block
                    assert 'QuantumNumbers -> {Y -> 0}' in block
                    # real singlets are self-conjugate: no QuantumNumbers (FeynRules LoadModel::QN)
                    for single in mm.non_conjugated_multiplets:
                        if single.selfConjugate and single.isospin.dim == 1:
                            sblock = re.search(rf'ClassName\s*->\s*{single.fr_symbol()},.*?\n\s*\}}',
                                               fr, re.S).group()
                            assert 'SelfConjugate  -> True' in sblock and 'QuantumNumbers' not in sblock
                    defs = re.search(r'Definitions\s*->\s*\{(.*)\}', block).group(1)
                    plus = re.search(rf'{s}\[1\] :> ([^,]+?)(?:, {s}\[|$)', defs).group(1)
                    minus = re.search(rf'{s}\[3\] :> ([^,]+?)(?:, {s}\[|$)', defs).group(1)
                    assert minus.startswith('-'), defs
                    assert f'1/2 DC[{s}bar[' in fr
                    # its neutral component mixes through a real orthogonal matrix
                    for mix in set(re.findall(r'UMixR\w+', defs)):
                        ext = re.search(rf'\n\s*{mix} == \{{.*?\n\s*\}}', fr, re.S).group()
                        assert 'ComplexParameter -> False' in ext.replace('  ', ''), ext
                        assert f'I{mix} ==' not in fr and f'{mix}[i_, j_]' not in fr
                    found += 1
    assert found, 'no valid SDDM model with a real scalar triplet'


def _sddm_doublet_partner_fr():
    return _sddm_fr({'PhiI1Y0', 'PhiI2Y3'})


def _sddm_fr(names):
    return _sddm_model(names).to_feynrules()


def _sddm_model(names):
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    sddm = {"masses": {1000022: 89.88887, 1000023: 126.046269, 1000024: 155.494055, 1000025: 451.886204},
            "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                       1000025: {(1000023, 25): 1.0}},
            "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        from uvbuilder.tools.reps import GlobalU1s
        # with and without a dark charge (Majorana and Dirac DM)
        dm = {'globalU1s': [GlobalU1s(), GlobalU1s(D=1)]}
        for pm in ModelBuilder(modelData=sddm, possible_dm_reps=dm).generateParticleModels():
            for mm in PartitionsBuilder(pm).get_multipletModels():
                if mm.is_valid and {m.fr_symbol() for m in mm.non_conjugated_multiplets} == names:
                    return FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr')
    raise AssertionError(f'no SDDM model with multiplets {names}')


def test_complex_component_in_real_block_uses_two_rows():
    '''X0 = (phi_R + i phi_I)/sqrt(2): with any orthogonal U the neutral kinetic terms
    1/2 (dS)^2 + |dX0|^2 are canonical, 1/2 sum_j (dZ_j)^2.'''
    import numpy as np
    fr = _sddm_doublet_partner_fr()
    assert 'PhiI1Y0 :> Module[{j}, UMixRZ1Z2Z3[1, j] PhiQ0C1[j]]' in fr
    x0 = 'PhiI2Y3[2] :> Module[{j}, (UMixRZ1Z2Z3[2, j] + I*UMixRZ1Z2Z3[3, j])/Sqrt[2] PhiQ0C1[j]]'
    assert x0 in fr
    U, _ = np.linalg.qr(np.arange(1, 10, dtype=float).reshape(3, 3) + np.eye(3))
    s, x = U[0], (U[1] + 1j * U[2]) / np.sqrt(2)
    K = 0.5 * np.outer(s, s) + np.real(np.outer(x.conj(), x))
    assert np.allclose(K, 0.5 * np.eye(3))


# -- mass terms, two-index invariants, hermitian operators --

def _lint_lines(fr):
    '''Mass and interaction terms (LBSMMass, LBSMInt), one per line, without comments.'''
    lines = []
    for piece in ('LBSMMass', 'LBSMInt'):
        body = re.search(rf'{piece} := ExpandIndices\[\n(.*?)\n\s*, FlavorExpand', fr, re.S)
        if body:
            lines += [l.strip() for l in body.group(1).splitlines() if not l.strip().startswith('(*')]
    return lines


def test_scalar_mass_terms_have_no_cg_and_no_hc():
    lines = _lint_lines(_sddm_doublet_partner_fr())
    assert '- (1/2)*MPhiI1Y0^2 PhiI1Y0 PhiI1Y0' in lines
    assert '- MPhiI2Y3^2 PhiI2Y3bar[i1] PhiI2Y3[i1]' in lines
    lines = _lint_lines(_sddm_fr({'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'}))
    assert '- (1/2)*MPhiI3Y0^2 PhiI3Y0bar[i1] PhiI3Y0[i1]' in lines


def test_hermitian_operators_without_hc_and_conjugates_emitted_once():
    fr = _sddm_doublet_partner_fr()
    (quartic,) = [l for l in _lint_lines(fr) if l.startswith('- lambdaPhiI1Y0PhiI1Y0 ')]
    assert 'HC[' not in quartic                                   # S S H^dagger H
    assert 'muPhiI1Y0PhiI2Y3bar' not in fr and 'lambdaPhiI2Y3barPhiI2Y3bar' not in fr
    (mu,) = [l for l in _lint_lines(fr) if l.startswith('- muPhiI1Y0PhiI2Y3 ')]
    assert 'HC[' in mu
    fr = _sddm_fr({'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'})
    assert not any('HC[' in l for l in _lint_lines(fr))          # real fields only


def test_two_index_invariants_have_unit_entries():
    op = Operator(HiggsDoublet, HiggsDoublet.conjugateMultiplet)
    assert set(op.singlets.values()) == {1}                       # H^dagger H
    assert op.is_hermitian
    Du = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(1,),
                             hypercharge=U1Rep(1, 2)), label='Du', fillMultiplet=True)
    Dd = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(1,),
                             hypercharge=U1Rep(-1, 2)), label='Dd', fillMultiplet=True)
    dirac = Operator(Du, Dd)
    assert {abs(c) for c in dirac.singlets.values()} == {1}       # charged and neutral mass |M|
    assert not dirac.is_hermitian



# -- several invariants for the same fields --

def test_su2_singlet_paths_count_invariants():
    from uvbuilder.tools.cg_tensors import su2_singlet_paths
    d2, d3 = SU2Rep(1,), SU2Rep(2,)
    assert len(su2_singlet_paths((d2, d2, d2, d2))) == 2
    assert len(su2_singlet_paths((d2, d2))) == 1
    assert len(su2_singlet_paths((d2, d2, d3))) == 1
    assert len(su2_singlet_paths((d2, d2, d3, d3))) == 2


def test_inert_doublet_masses_are_independent():
    '''Two invariants for X X^dagger H H^dagger plus (X H~)^2 and M^2 X^dagger X give
    independent H+, H0, A0 masses and a free Higgs coupling (Jacobian of rank 3).'''
    from uvbuilder.tools.cg_tensors import SU2CGTensor, su2_singlet_paths
    X = _scalar('X', 2, U1Rep(1, 2))
    Xc, H, Hc = X.conjugateMultiplet, HiggsDoublet, HiggsDoublet.conjugateMultiplet
    v, R, Im, Xp, Xq = sp.symbols('v R Im Xp Xq', real=True)
    vals = {X.components[U1Rep(1, 2)]: Xp, X.components[U1Rep(-1, 2)]: (R + sp.I * Im) / sp.sqrt(2),
            Xc.components[U1Rep(-1, 2)]: Xq, Xc.components[U1Rep(1, 2)]: (R - sp.I * Im) / sp.sqrt(2),
            H.components[U1Rep(1, 2)]: 0, H.components[U1Rep(-1, 2)]: v / sp.sqrt(2),
            Hc.components[U1Rep(-1, 2)]: 0, Hc.components[U1Rep(1, 2)]: v / sp.sqrt(2)}

    def value(op):
        return sp.expand(sum(c * sp.Mul(*[vals[f] for f in k]) for k, c in op.singlets.items()))

    paths = su2_singlet_paths(SU2CGTensor([X.isospin] * 4).isospins)
    ops = [Operator(Xc, X, Hc, H, path=p) for p in paths]
    assert all(op.is_hermitian for op in ops)
    l5 = [Operator(X, X, Hc, Hc, path=p) for p in paths]
    l5 = [op for op in l5 if op.singlets]
    assert len(l5) == 1                                  # identical fields: one invariant
    M2, a, b, c = sp.symbols('M2 a b c', real=True)
    o5 = value(l5[0])
    V = M2 * (Xp * Xq + (R**2 + Im**2) / 2) + a * value(ops[0]) + b * value(ops[1]) \
        + c * (o5 + sp.conjugate(o5))
    masses = [sp.diff(V, Xp, Xq), sp.diff(V, R, 2), sp.diff(V, Im, 2)]
    jac = sp.Matrix([[sp.diff(m, x) for x in (a, b, c)] for m in masses])
    assert jac.rank() == 3


def test_several_invariants_and_complex_couplings_in_output():
    fr = _sddm_doublet_partner_fr()
    lines = _lint_lines(fr)
    for k in (1, 2):
        (line,) = [l for l in lines if l.startswith(f'- lambdaPhiI2Y3PhiI2Y3barP{k} ')]
        assert f'SU2CG2222P{k}[' in line and 'HC[' not in line          # hermitian, real
    (l5,) = [l for l in lines if l.startswith('- lambdaPhiI2Y3PhiI2Y3 ')]
    assert 'HC[' in l5                                                   # single invariant, no suffix
    assert 'muPhiI1Y0PhiI2Y3 == {' in fr and 'RmuPhiI1Y0PhiI2Y3 + I*ImuPhiI1Y0PhiI2Y3' in fr
    fr = _sddm_fr({'ChiI1Y0Lw', 'ChiI2Y3Lw', 'ChiI2Ym3Lw'})
    assert 'P1' not in fr                                                # one invariant everywhere
    assert re.search(r'MChiI2Y3ChiI2Ym3 == \{[^}]*ComplexParameter -> False', fr.replace('  ', ''))


def test_conjugate_component():
    from tests._sddm import DOUBLET_PARTNER, REAL_TRIPLET, multiplet_model
    triplet = next(m for m in multiplet_model(REAL_TRIPLET).multiplets if m.label == 'PhiI3Y0')
    assert triplet.conjugate_component(U1Rep(1)) == (triplet[U1Rep(-1)], -1)     # T+* = -T-
    assert triplet.conjugate_component(U1Rep(0)) == (triplet[U1Rep(0)], 1)       # T0 real
    doublet = next(m for m in multiplet_model(DOUBLET_PARTNER).multiplets if m.label == 'PhiI2Y3')
    upper = doublet[U1Rep(1, 2)]
    assert doublet.conjugate_component(U1Rep(1, 2)) == (upper.conjugateField, 1)
