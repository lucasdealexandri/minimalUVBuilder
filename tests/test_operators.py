"""Operators of the multiplet models (tools/operators.py)."""

import sympy as sp

from tests._sddm import DOUBLET_PARTNER, MAJORANA, REAL_TRIPLET, multiplet_model


def _by_description(names):
    return {op._description(): op for op in multiplet_model(names).operators}


def test_sectors_and_descriptions_of_the_fermion_model():
    operators = _by_description(MAJORANA)
    assert {text: op.sector for text, op in operators.items()} == {
        'Majorana mass of ChiI1Y0': 'Mass',
        'Dirac mass of ChiI2Ym3 ChiI2Y3': 'Mass',
        'Yukawa coupling ChiI1Y0 ChiI2Y3 Phibar': 'Yukawa',
        'Yukawa coupling ChiI1Y0 ChiI2Ym3 Phi': 'Yukawa'}


def test_sectors_and_descriptions_of_the_scalar_models():
    operators = _by_description(DOUBLET_PARTNER)
    assert operators['Mass of PhiI1Y0'].sector == operators['Mass of PhiI2Y3'].sector == 'Mass'
    assert operators['Trilinear coupling PhiI1Y0 PhiI2Y3 Phibar'].sector == 'Trilinear'
    # the invariant is named only when the same fields have several
    assert operators['Quartic coupling PhiI2Y3 PhiI2Y3 Phibar Phibar'].sector == 'Quartic'
    assert {'Quartic coupling PhiI2Y3bar PhiI2Y3 Phibar Phi (SU(2) invariant P1)',
            'Quartic coupling PhiI2Y3bar PhiI2Y3 Phibar Phi (SU(2) invariant P2)'} <= set(operators)
    assert 'Mass of PhiI3Y0' in _by_description(REAL_TRIPLET)


def test_only_the_neutral_higgs_component_survives_at_the_vev():
    yukawa = _by_description(MAJORANA)['Yukawa coupling ChiI1Y0 ChiI2Ym3 Phi']
    assert len(yukawa.expand_singlets()) == 2
    (term,) = yukawa.expand_singlets(at_vev=True)
    assert [f.label for f in term.fields] == ['ChiI1Y0q0', 'ChiI2Ym3q0']


def test_mass_terms_of_the_singlet_doublet_model():
    '''-L = sum of the terms below + h.c., i.e. the neutral mass matrix in the basis
    (singlet, Y = -1/2, Y = +1/2) is [[M_s, y_d v/sqrt2, -y_u v/sqrt2], [., 0, -mu], [., -mu, 0]]
    and the charged state has Dirac mass +mu (Multiplet.mass_partner).'''
    v, Ms, mu, yu, yd = sp.symbols('v MChiI1Y0 MChiI2Y3ChiI2Ym3 yChiI1Y0ChiI2Y3 yChiI1Y0ChiI2Ym3')
    terms = {tuple(f.label for f in term.fields): term.coeff
             for op in multiplet_model(MAJORANA).operators
             for term in op.expand_singlets(at_vev=True)}
    expected = {('ChiI1Y0q0', 'ChiI1Y0q0'): Ms / 2,
                ('ChiI1Y0q0', 'ChiI2Ym3q0'): yd * v / sp.sqrt(2),
                ('ChiI1Y0q0', 'ChiI2Y3q0'): -yu * v / sp.sqrt(2),
                ('ChiI2Ym3q0', 'ChiI2Y3q0'): -mu,
                ('ChiI2Ym3qm3', 'ChiI2Y3q3'): mu}
    assert terms.keys() == expected.keys()
    assert all(sp.simplify(terms[key] - expected[key]) == 0 for key in expected)


def test_lagrangian_operators_keep_one_of_each_conjugate_pair():
    mm = multiplet_model(DOUBLET_PARTNER)
    key = lambda op, conj: tuple(sorted((m.isSM, m.rep.conjRep if conj else m.rep, m.index)
                                        for m in op.multiplets))
    ops = mm.lagrangian_operators
    keys = {key(op, False) for op in ops}
    assert all(key(op, True) == key(op, False) or key(op, True) not in keys for op in ops)
    assert all(any(op is o for o in mm.operators) for op in ops)
    assert len(ops) < len(mm.operators)          # complex scalars: O and O^dagger both generated


def _terms(terms):
    return {tuple(f.label for f in t.fields): t.coeff for t in terms}


def test_real_triplet_mass_term_is_positive_definite():
    '''sum_i T^*_i T_i with T_- = -(T_+)^*: -M^2 T_- T_+ + 1/2 M^2 T_0 T_0 = M^2 |T_+|^2 + 1/2 M^2 T_0^2,
    the term fr_mass_term writes (the CG singlet has the opposite sign).'''
    M = sp.Symbol('MPhiI3Y0')
    op = _by_description(REAL_TRIPLET)['Mass of PhiI3Y0']
    assert _terms(op.expand_singlets()) == {('PhiI3Y0qm3', 'PhiI3Y0q3'): -M**2,
                                            ('PhiI3Y0q0', 'PhiI3Y0q0'): M**2 / 2}
    assert op.is_hermitian and op.lagrangian_terms() == op.expand_singlets()


def test_hermitian_conjugate_of_a_trilinear_at_the_vev():
    v, mu = sp.symbols('v muPhiI1Y0PhiI2Y3')
    op = _by_description(DOUBLET_PARTNER)['Trilinear coupling PhiI1Y0 PhiI2Y3 Phibar']
    assert not op.is_hermitian
    assert _terms(op.expand_singlets(at_vev=True)) == {
        ('PhiI1Y0q0', 'PhiI2Y3q0'): -sp.sqrt(2) * mu * v / 2}
    assert _terms(op.hermitian_conjugate(at_vev=True)) == {
        ('PhiI1Y0q0', 'PhiI2Y3q0bar'): -sp.sqrt(2) * sp.conjugate(mu) * v / 2}
    assert op.lagrangian_terms(at_vev=True) == (op.expand_singlets(at_vev=True)
                                                + op.hermitian_conjugate(at_vev=True))
    assert op._conjugated(op._conjugated(op.singlets)) == op.singlets


def test_lagrangian_terms_add_the_hc_exactly_where_the_fr_writes_it():
    from uvbuilder.models.feynrulesmodel import IndexTracker
    for names in (DOUBLET_PARTNER, REAL_TRIPLET, MAJORANA):
        for op in multiplet_model(names).lagrangian_operators:
            (interaction,) = next(v for k, v in op.to_feynrules(IndexTracker()).items()
                                  if k.endswith('Terms'))
            with_hc = len(op.lagrangian_terms()) > len(op.expand_singlets())
            assert with_hc == ('HC[' in interaction), op._description()


def test_majorana_mass_of_a_fermion_triplet_is_positive_for_every_component():
    '''The CG singlet of a real triplet is 2 X+ X- - X0 X0; the mass term pairs each component
    with its partner (Multiplet.mass_partner), so a positive mass gives +M to the neutral
    component, which is its own partner.'''
    from uvbuilder.models.feynrulesmodel import IndexTracker
    from uvbuilder.states.multiplet import Multiplet
    from uvbuilder.tools.operators import Operator
    from uvbuilder.tools.reps import WEYL_LEFT, BrokenRep, SU2Rep, SU3Rep, U1Rep
    W = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(2,),
                            hypercharge=U1Rep(0,)), label='W', fillMultiplet=True)
    M = sp.Symbol('M')
    op = Operator(W, W, parameter=M / 2)
    plus, zero, minus = (W.components[U1Rep(q)] for q in (1, 0, -1))
    assert {tuple(sorted(t.fields, key=str)): t.coeff for t in op.expand_singlets()} == {
        tuple(sorted((minus, plus), key=str)): -M, (zero, zero): M / 2}
    interaction = op.to_feynrules(IndexTracker())['MassTerms'][0]
    assert interaction.startswith('+ (1/2)*M SU2CG33') and '+ HC[' in interaction
