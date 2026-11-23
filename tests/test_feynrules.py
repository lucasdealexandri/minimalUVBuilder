"""Tests for the FeynRules class declarations of particles and their
Weyl/scalar components (tools/components.py, states/particle.py).

Covers:
  * the <particle label>Lw / <particle label>Rw component label convention,
  * the ``bar`` conjugation marker (eLw.dagger == eLwbar, involution),
  * WeylComponent.to_feynrules (W[%i] class declaration),
  * Particle.to_feynrules for Dirac, Weyl and scalar particles.
"""

import re

import pytest

from uvbuilder.tools.reps import (UnbrokenRep, BrokenRep, DIRAC_SPINOR, WEYL_LEFT,
                        WEYL_RIGHT, SCALAR, SU3Rep, SU2Rep, U1Rep)
from uvbuilder.states.particle import Particle
from uvbuilder.states.physicsUnits import GeV
from uvbuilder.states.multiplet import Multiplet
from uvbuilder.tools.components import WeylComponent, ScalarComponent
from uvbuilder.tools.operators import Operator
from uvbuilder.models.feynrulesmodel import IndexTracker


def _dirac(label, charge=U1Rep(-1,), pdg=11, mass=None, width=None):
    return Particle(
        label=label, pdg=pdg, isSM=True,
        rep=UnbrokenRep(spin=DIRAC_SPINOR, colordim=SU3Rep(0, 0),
                        eCharge=charge),
        mass=mass, totalwidth=width)


def _weyl(label, left=True, charge=U1Rep(0,), pdg=12):
    return Particle(
        label=label, pdg=pdg, isSM=True,
        rep=UnbrokenRep(spin=WEYL_LEFT if left else WEYL_RIGHT,
                        colordim=SU3Rep(0, 0), eCharge=charge))


def _scalar(label, pdg=25, charge=U1Rep(0,)):
    return Particle(
        label=label, pdg=pdg, isSM=True,
        rep=UnbrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0),
                        eCharge=charge))


class TestComponentLabels:

    def test_dirac_components_lw_rw(self):
        e = _dirac('e-')
        comps = e.components
        assert [c.label for c in comps] == ['e-Lw', 'e-Rw']
        assert [c.chirality for c in comps] == ['LH', 'RH']

    def test_weyl_component_suffix_by_chirality(self):
        assert _weyl('n', left=True).components[0].label == 'nLw'
        assert _weyl('nr', left=False).components[0].label == 'nrRw'

    def test_dagger_appends_bar(self):
        e = _dirac('e-')
        lw = e.components[0]
        assert lw.dagger.label == 'e-Lwbar'
        assert lw.dagger.dagger.label == 'e-Lw'  # involution

    def test_scalar_dagger_appends_bar(self):
        s = _scalar('h', charge=U1Rep(1,))
        assert s.components[0].dagger.label == 'hbar'
        assert s.components[0].dagger.dagger.label == 'h'


class TestWeylComponentFr:

    def test_weyl_block_left(self):
        t = IndexTracker()
        e = _dirac('e-')
        block = e.components[0].to_feynrules(t)['Classes'][0]
        assert 'W[101] == {' in block   # class numbers start at 101 (SM.fr uses the low ones)
        assert 'ClassName      -> e-Lw' in block
        assert 'Chirality      -> Left' in block
        assert 'Unphysical     -> True' in block
        assert 'QuantumNumbers -> {Q -> -1}' in block

    def test_weyl_block_right(self):
        t = IndexTracker()
        block = _weyl('nr', left=False).components[0].to_feynrules(t)['Classes'][0]
        assert 'W[101] == {' in block   # class numbers start at 101 (SM.fr uses the low ones)
        assert 'Chirality      -> Right' in block

    def test_weyl_block_with_family(self):
        t = IndexTracker()
        family = t.family('FGenQ0C1', 2)
        e = _dirac('e-')
        block = e.components[0].to_feynrules(t, family)['Classes'][0]
        assert 'Indices        -> {Index[FGenQ0C1]}' in block
        assert 'FlavorIndex    -> FGenQ0C1' in block


class TestParticleFr:

    def test_dirac_emits_weyl_components_and_f_class(self):
        t = IndexTracker()
        e = _dirac('e-', mass=0.000511 * GeV, width=0.0 * GeV)
        classes = e.to_feynrules(t)['Classes']
        joined = '\n'.join(classes)
        assert 'ClassName      -> e-Lw' in joined
        assert 'ClassName      -> e-Rw' in joined
        assert 'WeylComponents -> {e-Lw, e-Rw}' in joined
        assert 'PDG            -> {11}' in joined
        assert 'Mass           -> {Me-, 0.000511}' in joined

    def test_dirac_omits_unset_mass_and_width(self):
        t = IndexTracker()
        joined = '\n'.join(_dirac('e-').to_feynrules(t)['Classes'])
        assert 'Mass' not in joined
        assert 'Width' not in joined

    def test_weyl_emits_single_w_class(self):
        t = IndexTracker()
        joined = '\n'.join(_weyl('n').to_feynrules(t)['Classes'])
        assert 'ClassName      -> nLw' in joined
        assert 'F[' not in joined

    def test_scalar_emits_nothing(self):
        t = IndexTracker()
        assert _scalar('h').to_feynrules(t) == {'Classes': []}


class TestIndexTracker:

    def test_new_symbol_unique_per_style(self):
        t = IndexTracker()
        assert [t.new_symbol('sp') for _ in range(3)] == ['sp1', 'sp2', 'sp3']
        assert t.new_symbol('i') == 'i1'
        assert t.new_symbol('c') == 'c1'

    def test_reset_symbols_per_term(self):
        t = IndexTracker()
        assert [t.new_symbol('i') for _ in range(2)] == ['i1', 'i2']
        t.reset_symbols()
        assert t.new_symbol('i') == 'i1'

    def test_index_consistency_check(self):
        t = IndexTracker()
        t.index('SU2D', 2, SU2Rep)
        assert t.index('SU2D', 2, SU2Rep).range == 2
        with pytest.raises(ValueError):
            t.index('SU2D', 3, SU2Rep)

    def test_gauge_and_lorentz(self):
        t = IndexTracker()
        m = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(1, 0),
                                isospin=SU2Rep(1,), hypercharge=U1Rep(1, 2)),
                      label='Q', fillMultiplet=True)
        assert [i.name for i in t.gauge(m.rep)] == ['Colour', 'SU2D']
        assert t.lorentz(m.rep.spin).name == 'W'

    def test_indices_block(self):
        t = IndexTracker()
        t.index('Colour', 3, SU3Rep)
        block = t.indices_block()
        assert 'IndexRange[Index[Colour]] = Range[3];' in block
        assert 'IndexStyle[Colour,c];' in block


class TestMultipletConjugate:

    def test_conjugate_uses_bar(self):
        m = Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0),
                                isospin=SU2Rep(0,), hypercharge=U1Rep(1,)),
                      label='X', fillMultiplet=True)
        assert m.conjugateMultiplet.label == 'Xbar'
        assert m.conjugateMultiplet.conjugateMultiplet is m


class TestMultipletKinetic:

    def _multiplet(self, spin, iso, y, label):
        return Multiplet(BrokenRep(spin=spin, colordim=SU3Rep(0, 0),
                                   isospin=iso, hypercharge=y),
                         label=label, fillMultiplet=True)

    def test_weyl_kinetic_si_form(self):
        t = IndexTracker()
        m = self._multiplet(WEYL_LEFT, SU2Rep(0,), U1Rep(0,), 'S')
        term = m.generate_fr_kinetic(t)
        assert 'sibar[mu, sp1, sp2]' in term
        assert 'SLwbar[sp1].DC[SLw[sp2], mu]' in term
        assert 'Index[mu]' not in term

    def test_scalar_kinetic_dc_form(self):
        t = IndexTracker()
        m = self._multiplet(SCALAR, SU2Rep(0,), U1Rep(1,), 'H')
        assert m.generate_fr_kinetic(t) == '+ DC[Hbar, mu].DC[H, mu]'

    def test_kinetic_gauge_labels(self):
        t = IndexTracker()
        m = self._multiplet(WEYL_LEFT, SU2Rep(1,), U1Rep(1, 2), 'Du')
        term = m.generate_fr_kinetic(t)
        assert 'DuLwbar[sp1, i1]' in term
        assert 'DuLw[sp2, i1]' in term


class TestOperatorTracker:

    def _multiplet(self, label, iso, y, index=1):
        return Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0),
                                   isospin=iso, hypercharge=y),
                         label=label, fillMultiplet=True, index=index)

    def test_fr_term_uses_tracker_labels(self):
        t = IndexTracker()
        A = self._multiplet('A', SU2Rep(1,), U1Rep(0,), 1)
        B = self._multiplet('B', SU2Rep(1,), U1Rep(0,), 2)
        term = Operator(A, B).fr_term(t)
        assert term == '1 SU2CG22[i1, i2] A[i1] B[i2]'

    def test_tracker_labels_do_not_collide_across_operators(self):
        t = IndexTracker()
        A = self._multiplet('A', SU2Rep(1,), U1Rep(0,), 1)
        B = self._multiplet('B', SU2Rep(1,), U1Rep(0,), 2)
        C = self._multiplet('C', SU2Rep(1,), U1Rep(0,), 3)
        assert Operator(A, B).fr_term(t) == '1 SU2CG22[i1, i2] A[i1] B[i2]'
        assert Operator(A, C).fr_term(t) == '1 SU2CG22[i1, i2] A[i1] C[i2]'

    def test_fr_definition_alias(self):
        t = IndexTracker()
        A = self._multiplet('A', SU2Rep(1,), U1Rep(0,), 1)
        B = self._multiplet('B', SU2Rep(1,), U1Rep(0,), 2)
        op = Operator(A, B)
        definition = op.fr_definition(t)
        assert 'SU2CG22[i' in definition
        assert 'A[i' in definition and 'B[i' in definition


class TestMixingStub:

    def test_generate_mixing_matrices_stub(self):
        from uvbuilder.models.multipletModel import MultipletModel
        m = Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0),
                                isospin=SU2Rep(0,), hypercharge=U1Rep(0,)),
                      label='X', fillMultiplet=True)
        mm = MultipletModel((m,))
        assert mm.generate_mixing_matrices(IndexTracker()) == []


def _sddm_model():
    from uvbuilder.states.SMparticles import nue, d, u, higgs
    from uvbuilder.states.particle import Particle
    from uvbuilder.states.physicsUnits import GeV
    from uvbuilder.vertices.particleVertex import ParticleVertex
    from uvbuilder.models.particleModel import ParticleModel
    from uvbuilder.models.multipletModel import MultipletModel

    charged_rep = UnbrokenRep(spin=DIRAC_SPINOR, colordim=SU3Rep(0, 0),
                              eCharge=U1Rep(1,))
    neutral_rep = UnbrokenRep(spin=DIRAC_SPINOR, colordim=SU3Rep(0, 0),
                              eCharge=U1Rep(0,))
    C = Particle(label='CPlus', pdg=1000024, mass=155.494055 * GeV,
                 rep=charged_rep, isSM=False)
    N1 = Particle(label='N1', pdg=1000022, mass=89.88887 * GeV,
                  rep=neutral_rep, isSM=False)
    N2 = Particle(label='N2', pdg=1000023, mass=126.046269 * GeV,
                  rep=neutral_rep, isSM=False)
    N3 = Particle(label='N3', pdg=1000025, mass=451.886204 * GeV,
                  rep=neutral_rep, isSM=False)
    vertices = [ParticleVertex(N2, [N1, nue, nue]),
                ParticleVertex(C, [N1, u, d]),
                ParticleVertex(N3, [N2, higgs])]
    pm = ParticleModel(vertices)

    Du = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0),
                             isospin=SU2Rep(1,), hypercharge=U1Rep(1, 2)),
                   label='Du', labels={U1Rep(-1, 2): 'Du0', U1Rep(1, 2): 'DuP'},
                   fillMultiplet=True)
    Dd = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0),
                             isospin=SU2Rep(1,), hypercharge=U1Rep(-1, 2)),
                   label='Dd', labels={U1Rep(-1, 2): 'DdM', U1Rep(1, 2): 'Dd0'},
                   fillMultiplet=True)
    S = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0),
                            isospin=SU2Rep(0,), hypercharge=U1Rep(0,)),
                  label='S0', labels={U1Rep(0): 'S0'}, fillMultiplet=True)
    return MultipletModel((Du, Dd, S), pm), pm


def _sddm_frm(**kwargs):
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    mm, pm = _sddm_model()
    mixings = mm.possible_mixings[0]
    return FeynRulesModel(mm, pm, mixings, 'test.fr', **kwargs)


class TestPhysicalWeylComponents:

    def test_multi_member_class_components(self):
        classes = _sddm_frm().get_classes_definitions()
        # Majorana (self-conjugate) family: a single left-handed Weyl
        # component, SelfConjugate True, no right-handed partner.
        assert 'ClassName      -> ChiQ0C1Lw' in classes
        assert 'ClassName      -> ChiQ0C1Rw' not in classes
        assert 'Chirality      -> Left' in classes
        assert 'WeylComponents -> ChiQ0C1Lw' in classes
        assert 'SelfConjugate  -> True' in classes
        assert 'ClassMembers   -> {N1,N2,N3}' in classes
        assert 'FlavorIndex    -> FGenQ0C1' in classes

    def test_single_particle_components(self):
        classes = _sddm_frm().get_classes_definitions()
        assert 'ClassName      -> CPlusLw' in classes
        assert 'ClassName      -> CPlusRw' in classes
        assert 'WeylComponents -> {CPlusLw, CPlusRw}' in classes


class TestGaugeDefinitions:
    """Gauge multiplet components rotated into the physical classes."""

    def test_implicit_family_sum(self):
        classes = _sddm_frm().get_classes_definitions()
        # 1-based mixing row, implicit sum over the family index j
        assert 'S0Lw[s_] :> Module[{j}, UMixMN1N2N3[3, j] ChiQ0C1Lw[s, j]]' in classes
        assert 'DdLw[s_, 1] :> Module[{j}, UMixMN1N2N3[1, j] ChiQ0C1Lw[s, j]]' in classes
        assert 'DuLw[s_, 2] :> Module[{j}, UMixMN1N2N3[2, j] ChiQ0C1Lw[s, j]]' in classes

    def test_single_member_direct_map(self):
        classes = _sddm_frm().get_classes_definitions()
        # one-member class: no Module, no family index; RH component enters as Rwbar,
        # with the phase that makes the Dirac mass positive
        assert 'DuLw[s_, 1] :> CPlusLw[s]' in classes
        assert 'DdLw[s_, 2] :> PCPlus CPlusRwbar[s]' in classes

    def test_higgs_tilde_class(self):
        # SM.fr has no H~ class: the generated file declares it ...
        from pathlib import Path
        classes = _sddm_frm(sm_file=Path(__file__).resolve().parent / 'data' / 'SM.fr').get_classes_definitions()
        assert 'ClassName      -> Phitilde' in classes
        assert 'Phitilde[1] :> Phibar[2]' in classes
        assert 'Phitilde[2] :> -Phibar[1]' in classes

    def test_higgs_tilde_not_redeclared_with_smmod(self):
        # ... SMMod.fr (the default) defines it, so it must not be declared again
        frm = _sddm_frm()
        assert 'ClassName      -> Phitilde' not in frm.get_classes_definitions()
        assert any('Phitilde[' in t for key, terms in frm._chunks.items() if key.endswith('Terms')
                   for t in terms)


class TestCGDefinitions:

    def test_su2_cg_tensor_defined(self):
        params = _sddm_frm().get_parameters()
        assert 'SU2CG22 == {' in params
        assert 'ParameterType    -> Internal' in params
        assert 'Indices          -> {Index[SU2D], Index[SU2D]}' in params
        assert 'SU2CG22[1, 2] -> 1' in params
        assert 'SU2CG22[2, 1] -> -1' in params

    def test_cg_definitions_deduped(self):
        params = _sddm_frm().get_parameters()
        assert params.count('SU2CG22 == {') == 1


class TestParameters:

    def test_operator_parameters_declared(self):
        params = _sddm_frm().get_parameters()
        for name in ('MDuDd', 'yDuS0', 'yDdS0', 'MS0'):
            assert f'{name} == {{' in params
        assert 'ParameterType    -> External' in params
        assert 'Value            -> 1' in params

    def test_interaction_orders(self):
        frm = _sddm_frm()
        params = frm.get_parameters()
        blocks = dict(re.findall(r'(\w+) == \{(.*?)\n\s*\}', params, re.S))
        # Yukawas (3 fields, complex couplings): {QED, 1}; bilinears (mass terms only): no order
        assert 'InteractionOrder -> {QED, 1}' in blocks['yDuS0']
        assert 'InteractionOrder -> {QED, 1}' in blocks['yDdS0']
        assert 'RyDuS0' in blocks and 'IyDuS0' in blocks
        assert 'InteractionOrder' not in blocks['MDuDd']
        assert 'InteractionOrder' not in blocks['MS0']

    def test_mixing_matrix_blocks(self):
        frm = _sddm_frm()
        params = frm.get_parameters()
        assert 'RUMixMN1N2N3 == {' in params and 'BlockName        -> UMIXMN1N2N3' in params
        assert 'IUMixMN1N2N3 == {' in params and 'BlockName        -> IMUMIXMN1N2N3' in params
        assert 'UMixMN1N2N3[i_, j_] :> RUMixMN1N2N3[i, j] + I*IUMixMN1N2N3[i, j]' in params
        assert 'Unitary          -> True' in params

    def test_parameter_sections_in_order_and_comma_separated(self):
        params = _sddm_frm().get_parameters()
        headers = ['(* Mass terms: external parameters *)', '(* Yukawa couplings: external parameters *)',
                   '(* Yukawa couplings: internal parameters *)', '(* Mixing matrices: external parameters *)',
                   '(* Mixing matrices: internal parameters *)',
                   '(* SU(2) Clebsch-Gordan tensors: internal parameters *)']
        positions = [params.index(h) for h in headers]
        assert positions == sorted(positions)
        entries = params.count(' == {')
        assert params.count('\n  },') == entries - 1 and params.rstrip().endswith('}\n\n};')
        assert '\n,' not in params

    def test_complex_coupling_parts_share_order_block(self):
        from tests.test_conjugation import _sddm_doublet_partner_fr
        fr = _sddm_doublet_partner_fr()
        blocks = dict(re.findall(r'(\w+) == \{(.*?)\n\s*\}', fr, re.S))
        real, imag = blocks['RlambdaPhiI2Y3PhiI2Y3'], blocks['IlambdaPhiI2Y3PhiI2Y3']
        assert 'BlockName        -> BSMQUARTIC' in real and 'BlockName        -> IMBSMQUARTIC' in imag
        order = re.search(r'OrderBlock\s+-> (\d+)', real).group(1)
        assert f'OrderBlock       -> {order}' in imag
        assert '(SU(2) invariant P1)' in blocks['lambdaPhiI2Y3PhiI2Y3barP1']
        assert 'invariant' not in blocks['lambdaPhiI2Y3PhiI2Y3']           # single invariant left


class TestInteractionTerms:

    def test_terms_are_dot_products_with_hc(self):
        frm = _sddm_frm()
        terms = frm._chunks['MassTerms'] + frm._chunks['YukawaTerms']
        joined = ''.join(terms)
        bilinear = [t for t in terms if 'MDuDd' in t][0]
        assert '(DdLw[sp' in bilinear and '.DuLw[sp' in bilinear
        assert ' - HC[' in bilinear
        assert 'ncc[' not in joined
        yukawa = [t for t in terms if 'yDuS0' in t][0]
        assert '(S0Lw[sp' in yukawa and '.DuLw[sp' in yukawa
        assert 'Phitilde[' in yukawa        # conjugate Higgs as H~, never Phibar
        assert ' - HC[' in yukawa
        lint = frm.get_lagrangian()
        assert '- MDuDd SU2CG22[i1, i2] (DdLw[sp1, i1].DuLw[sp1, i2]) - HC[' \
               'MDuDd SU2CG22[i1, i2] (DdLw[sp1, i1].DuLw[sp1, i2])]' in lint

    def test_flavor_expand_includes_su2w(self):
        frm = _sddm_frm()
        expand = re.findall(r'FlavorExpand -> \{([^}]*)\}', frm.get_lagrangian())
        assert len(expand) == 3 and all('SU2W' in e for e in expand)

    def test_lint_block_and_lbsm(self):
        fr = _sddm_frm().to_feynrules()
        for piece in ('LBSMKin', 'LBSMMass', 'LBSMInt'):
            assert f'{piece} := ' in fr
        assert 'ExpandIndices[' in fr
        assert 'LBSM := LBSMKin + LBSMMass + LBSMInt;' in fr



class TestValidMixing:

    def test_sddm_has_a_valid_mixing(self):
        mm, _ = _sddm_model()
        assert len(mm.possible_mixings) == 1
        assert mm.valid_mixing(mm.possible_mixings[0])

    def test_z2_even_decay_rejected(self):
        # X -> SM SM cannot come from operators built on BSM pairs (implicit Z2)
        from uvbuilder.states.SMparticles import higgs
        from uvbuilder.vertices.particleVertex import ParticleVertex
        mm, pm = _sddm_model()
        mixing = mm.possible_mixings[0]          # valid before the extra decay
        n3 = pm.pdgDict[1000025]
        pm.vertices.append(ParticleVertex(n3, [higgs, higgs]))
        assert not mm.valid_mixing(mixing)


class TestMixChunks:

    def test_only_matrices_with_emitted_rows_are_declared(self):
        # SDDM with a dark charge: complex-scalar variants have conjugate partner
        # multiplets that are never declared; their mixing matrices must not be
        # emitted (regression: `and` instead of `&` emitted every matrix)
        import warnings
        from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
        from uvbuilder.partition.partitions_builder import PartitionsBuilder
        from uvbuilder.models.feynrulesmodel import FeynRulesModel
        sddm = {"masses": {1000022: 89.88887, 1000023: 126.046269, 1000024: 155.494055, 1000025: 451.886204},
                "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                           1000025: {(1000023, 25): 1.0}},
                "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}
        checked = skipped = 0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for pm in ModelBuilder(modelData=sddm).generateParticleModels():
                for mm in PartitionsBuilder(pm).get_multipletModels():
                    if not mm.is_valid:
                        continue
                    mixing = mm.possible_mixings[0]
                    fr = FeynRulesModel(mm, pm, mixing, 'test.fr').to_feynrules()
                    emitted = {f for m in mm.non_conjugated_multiplets for f in m.fields}
                    for mat in mixing.get_mixing_matrices():
                        if mat.dim == 1:
                            continue
                        # complex: R<sym> real part; real scalars: <sym> itself (real orthogonal)
                        declared = (f'R{mat.standard_symbol} == {{' in fr
                                    or f' {mat.standard_symbol} == {{' in fr)
                        assert declared == bool(emitted & set(mat.f_components)), mat.standard_symbol
                        checked += 1
                        skipped += not declared
        assert checked and skipped      # both cases actually exercised


@pytest.fixture(scope='module')
def sm_map():
    from uvbuilder.tools.fr_names import SMClassMap, load_model
    return SMClassMap(load_model())


class TestSMTranslation:
    """SM multiplets are written as the classes of the SM model file (SMMod.fr),
    matched by quantum numbers; SM indices are not re-declared."""

    def test_every_sm_multiplet_has_a_class(self, sm_map):
        from uvbuilder.states.SMfields import SM_MULTIPLETS
        names = {m.label: sm_map.symbol(m) for m in SM_MULTIPLETS}
        assert {names[l] for l in ('lL1', 'lL2', 'lL3')} == {'LLw'}
        assert {names[l] for l in ('QL1', 'QL2', 'QL3')} == {'QLw'}
        assert {names[l] for l in ('eRbar', 'muRbar', 'tauRbar')} == {'eRwbar'}
        assert {names[l] for l in ('uRbar', 'cRbar', 'tRbar')} == {'uRwbar'}
        assert {names[l] for l in ('dRbar', 'sRbar', 'bRbar')} == {'dRwbar'}
        assert names['Phi'] == 'Phi'
        higgs_conj = [m for m in SM_MULTIPLETS if m.label != 'Phi' and m.rep.spin == SCALAR]
        assert [sm_map.symbol(m) for m in higgs_conj] == ['Phitilde']

    def test_render_index_order_and_generation(self, sm_map):
        from uvbuilder.states.SMfields import SM_MULTIPLETS
        by_label = {m.label: m for m in SM_MULTIPLETS}
        assert sm_map.render(by_label['QL3'], spin='sp1', colour='c1', isospin='i1') == 'QLw[sp1, i1, 3, c1]'
        assert sm_map.render(by_label['lL2'], spin='sp1', isospin='i2') == 'LLw[sp1, i2, 2]'
        assert sm_map.render(by_label['tRbar'], spin='sp2', colour='c1') == 'uRwbar[sp2, 3, c1]'
        assert sm_map.render(by_label['tRbar'], spin='sp2', colour='c1', bar=True) == 'uRw[sp2, 3, c1]'
        with pytest.raises(ValueError):
            sm_map.render(by_label['QL1'], spin='sp1', isospin='i1')     # colour index missing

    def test_no_match_is_an_error(self, sm_map):
        weird = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(1,),
                                    hypercharge=U1Rep(3, 2)), label='X', fillMultiplet=True, isSM=True)
        with pytest.raises(LookupError):
            sm_map.symbol(weird)

    def test_operator_emits_sm_classes(self, sm_map):
        from uvbuilder.states.SMfields import SM_MULTIPLETS, HiggsDoublet
        by_label = {m.label: m for m in SM_MULTIPLETS}
        # lepton portal of an inert doublet: eta N L (operators always hold two BSM multiplets)
        eta = Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0), isospin=SU2Rep(1,),
                                  hypercharge=U1Rep(1, 2)), label='Eta', fillMultiplet=True)
        N = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(0,),
                                hypercharge=U1Rep(0,)), label='N', fillMultiplet=True)
        term = Operator(eta, N, by_label['lL2']).fr_term(IndexTracker(sm=sm_map))
        assert re.search(r'\(NLw\[sp1\]\.LLw\[sp1, i\d, 2\]\)', term), term
        assert 'lL2Lw' not in term
        # coloured: a BSM triplet with Y = 2/3 pairs with tRbar (3bar, Y = -2/3)
        X = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(1, 0), isospin=SU2Rep(0,),
                                hypercharge=U1Rep(2, 3)), label='X', fillMultiplet=True)
        term = Operator(X, by_label['tRbar']).fr_term(IndexTracker(sm=sm_map))
        assert term == '1 (XLw[sp1, c1].uRwbar[sp1, 3, c1])'
        # without an SM map (no SM file) the old internal names are kept
        assert 'tRbarLw' in Operator(X, by_label['tRbar']).fr_term(IndexTracker())

    def test_sm_indices_not_redeclared(self):
        block = _sddm_frm().get_indices()
        assert 'Index[SU2D]' not in block and 'IndexStyle[SU2D' not in block
        assert 'IndexRange[Index[FGenQ0C1]] = Range[3];' in block      # BSM family index still declared


def test_lorentz_index_not_declared():
    block = _sddm_frm().get_indices()
    assert 'Index[mu]' not in block and 'IndexStyle[mu' not in block


def test_output_parameter_renames():
    import sympy as sp
    t = IndexTracker()
    A = Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0), isospin=SU2Rep(0,),
                            hypercharge=U1Rep(0,)), label='A', fillMultiplet=True)
    op = Operator(A, A, parameter=sp.Symbol('mAA') / 2)
    t.parameter_names = {'mAA': 'mAA1'}
    term = op.fr_term(t)
    assert 'mAA1' in term and 'mAA ' not in term and op.parameter == sp.Symbol('mAA') / 2
    assert 'mAA1 == {' in op.to_feynrules(t)['Parameters'][0].render()


class TestRenames:
    """Invalid/reserved names are fixed when the FeynRules file is written (and only
    there): STOP labels its neutralino 'B', which is SM.fr's hypercharge boson and
    whose mass MB would overwrite the b-quark mass."""

    STOP = {"masses": {1000022: 180.603, 1000006: 1164.264, 1000021: 1500.0},
            "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}, 1000021: {(1000006, -6): 1.0}},
            "labels": {1000022: 'B', 1000006: 'st', 1000021: 'Go'}}

    @pytest.fixture(scope='class')
    @classmethod
    def stop(cls):
        import warnings
        from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
        from uvbuilder.partition.partitions_builder import PartitionsBuilder
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for pm in ModelBuilder(modelData=cls.STOP).generateParticleModels():
                for mm in PartitionsBuilder(pm).get_multipletModels():
                    if mm.is_valid:
                        return mm, pm
        pytest.fail('no valid STOP multiplet model')

    def test_early_stage_only_warns(self):
        from uvbuilder.models.finalStateModel import FinalStateModel
        with pytest.warns(UserWarning, match='Particle labels will be renamed'):
            fsm = FinalStateModel.fromDict(self.STOP)
        assert 'B' in {p.label for fs in fsm.pdgDict.values() for p in fs.getParticles()}

    def test_renamed_in_output_and_restored(self, stop):
        from uvbuilder.models.feynrulesmodel import FeynRulesModel
        mm, pm = stop
        with pytest.warns(UserWarning, match='particle B -> B1'):
            frm = FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr')
        fr = frm.to_feynrules()
        assert frm.renames['particle'] == {'B': 'B1'}
        assert 'ClassName      -> B,' not in fr and '{MB,' not in fr
        assert '{MB1, 180.603}' in fr
        assert pm.pdgDict[1000022].label == 'B'                    # the model keeps its label

    def test_strict_raises(self, stop):
        from uvbuilder.models.feynrulesmodel import FeynRulesModel
        mm, pm = stop
        with pytest.raises(ValueError, match='particle B -> B1'):
            FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr', strict=True)
        assert pm.pdgDict[1000022].label == 'B'


def test_bsm_triplet_uses_spherical_su2d3_declared_in_smmod():
    from uvbuilder.tools.fr_names import SMClassMap, load_model
    t = IndexTracker(sm=SMClassMap(load_model()))
    T = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 0), isospin=SU2Rep(2,),
                            hypercharge=U1Rep(0,)), label='T', fillMultiplet=True)
    assert [i.name for i in t.gauge(T.rep)] == ['SU2d3']
    assert 'TLwbar[sp1, i1].DC[TLw[sp2, i1], mu]' in T.generate_fr_kinetic(t)
    assert 'SU2d3' not in t.indices_block() and 'SU2W' not in t.indices_block()


# -- naming: default gauge labels and physical-class tags --

def test_default_gauge_labels():
    from uvbuilder.models.multipletModel import fr_component_label
    from uvbuilder.tools.reps import SCALAR, GlobalU1s
    doublet = BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0), isospin=SU2Rep(1,),
                        hypercharge=U1Rep(1, 2))
    assert fr_component_label(doublet) == 'PhiI2Y3'
    assert fr_component_label(doublet, U1Rep(1, 2)) == 'PhiI2Y3q6'
    assert fr_component_label(doublet, family_index=2) == 'PhiI2Y3F2'
    quark = BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(1, 0), isospin=SU2Rep(0,),
                      hypercharge=U1Rep(2, 3), globalU1s=GlobalU1s(B=U1Rep(1, 3)))
    assert fr_component_label(quark) == 'ChiC3I1Y4B1'
    assert fr_component_label(quark.conjRep) == 'ChiC3bI1Ym4Bm1'
    T = Multiplet(quark, label=fr_component_label(quark), fillMultiplet=True)
    assert T.fr_symbol() == 'ChiC3I1Y4B1Lw'


@pytest.mark.parametrize('charge, colour, glob, tag', [
    (U1Rep(2, 3), SU3Rep(1, 0), None, 'Q2C3'),
    (U1Rep(-1, 3), SU3Rep(0, 1), None, 'Qm1C3b'),
    (U1Rep(1), SU3Rep(0, 0), None, 'Q3C1'),
    (U1Rep(0), SU3Rep(0, 0), {'D': U1Rep(1)}, 'Q0C1D1'),
    (U1Rep(2, 3), SU3Rep(1, 0), {'B': U1Rep(1, 3)}, 'Q2C3B1'),
])
def test_physical_class_tag(charge, colour, glob, tag):
    from types import SimpleNamespace
    from uvbuilder.models.particleModel import ParticleModel
    from uvbuilder.tools.reps import SCALAR, GlobalU1s
    rep = UnbrokenRep(spin=SCALAR, colordim=colour, eCharge=charge,
                      globalU1s=GlobalU1s(**(glob or {})))
    assert ParticleModel._gauge_name([SimpleNamespace(rep=rep)]) == tag


def test_physical_class_tag_rejects_non_third_charges():
    from types import SimpleNamespace
    from uvbuilder.models.particleModel import ParticleModel
    from uvbuilder.tools.reps import SCALAR
    rep = UnbrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0), eCharge=U1Rep(1, 6))
    with pytest.raises(ValueError):
        ParticleModel._gauge_name([SimpleNamespace(rep=rep)])


def test_unsupported_colour_and_lone_complex_scalar_raise():
    from uvbuilder.tools.cg_tensors import feyn_index_type
    from uvbuilder.tools.fr_utils import formatVariable
    from uvbuilder.models.multipletModel import MultipletModel
    from uvbuilder.tools.reps import SCALAR
    for f in (feyn_index_type, formatVariable):
        with pytest.raises(ValueError):
            f(SU3Rep(3, 0))                                  # colour 10
    X = Multiplet(BrokenRep(spin=SCALAR, colordim=SU3Rep(1, 0), isospin=SU2Rep(0,),
                            hypercharge=U1Rep(2, 3)), label='X', fillMultiplet=True)
    with pytest.raises(ValueError, match='without its conjugate'):
        MultipletModel([X])


# -- LH Weyl fields in 3bar/6bar: declared through their RH conjugate in 3/6 --

def _stop_dirac_model():
    import warnings
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.tools.reps import GlobalU1s
    stop = {"masses": {1000022: 180.603, 1000006: 1164.264, 1000021: 1500.0},
            "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}, 1000021: {(1000006, -6): 1.0}},
            "labels": {1000022: 'B', 1000006: 'st', 1000021: 'Go'}}
    builder = ModelBuilder(modelData=stop, possible_dm_reps={'globalU1s': [GlobalU1s(B=0)]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for pm in builder.generateParticleModels():
            for mm in PartitionsBuilder(pm).get_multipletModels():
                if mm.is_valid and any(m.is_colour_conjugated for m in mm.non_conjugated_multiplets):
                    return FeynRulesModel(mm, pm, mm.possible_mixings[0], 'test.fr')
    raise AssertionError('no stop model with a colour-antitriplet Weyl field')


def test_antitriplet_weyl_declared_through_rh_conjugate():
    import warnings
    fm = _stop_dirac_model()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fr = fm.to_feynrules()
    assert 'SU3d3bar' not in fr and 'DC[ChiC3bI1Ym4Bm1Lw' not in fr
    rh = re.search(r'ClassName\s*->\s*ChiC3I1Y4B1Rw,.*?\n\s*\}', fr, re.S).group()
    assert 'Chirality      -> Right' in rh and 'Indices        -> {Index[Colour]}' in rh
    assert 'QuantumNumbers -> {Y -> 2/3, BaryonNumber -> 1/3}' in rh
    assert 'ChiC3I1Y4B1Rw[s_, c_] :> Conjugate[Pst] stRw[s, c]' in rh
    assert 'ChiC3bI1Ym4Bm1Lw[s_, c_] :> ChiC3I1Y4B1Rwbar[s, c]' in fr
    assert ('I sibar[mu, sp1, sp2] (ChiC3I1Y4B1Rw[sp1, c1].DC[ChiC3I1Y4B1Rwbar[sp2, c1], mu])'
            in fr)


def test_sm_su2_indices_expanded_without_bsm_doublets():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fr = _stop_dirac_model().to_feynrules()
    assert 'Phitilde[i1] Phi[i2]' in fr                            # Higgs portal, SU(2) singlets only
    expand = re.findall(r'FlavorExpand -> \{([^}]*)\}', fr)
    assert len(expand) == 3 and all({'SU2D', 'SU2W'} <= {i.strip() for i in e.split(',')} for e in expand)


def test_rh_conjugate_name_is_checked():
    fm = _stop_dirac_model()
    captured = []

    class Recorder:
        def check(self, derived, where):
            captured.extend(derived)
            return []

    fm._check_fixed_names(Recorder())
    assert ('ChiC3I1Y4B1Rw', False) in captured and ('ChiC3I1Y4B1Rwbar', False) in captured


def test_antitriplet_doublet_helper_phases():
    X = Multiplet(BrokenRep(spin=WEYL_LEFT, colordim=SU3Rep(0, 1), isospin=SU2Rep(1,),
                            hypercharge=U1Rep(1, 6)), label='X', fillMultiplet=True)
    assert X.is_colour_conjugated
    # X_m = (-1)^(j+m) (X^c_{-m})^*
    assert X.conjugate_helper_definitions('XRw') == \
        '{XLw[s_, c_, 1] :> -XRwbar[s, c, 2], XLw[s_, c_, 2] :> XRwbar[s, c, 1]}'


def test_weyl_kinetic_terms_use_sibar():
    R = Multiplet(BrokenRep(spin=WEYL_RIGHT, colordim=SU3Rep(0, 0), isospin=SU2Rep(0,),
                            hypercharge=U1Rep(0,)), label='R', fillMultiplet=True)
    assert not R.is_colour_conjugated
    assert R.generate_fr_kinetic(IndexTracker()) == \
        ' + I sibar[mu, sp1, sp2] (RRw[sp1].DC[RRwbar[sp2], mu])'


@pytest.mark.slow
def test_fr_output_independent_of_hash_seed():
    import os, subprocess, sys
    from pathlib import Path
    names = {'ChiI1Y0D1F1Lw', 'ChiI1Y0D1F2Lw', 'ChiI1Y0Dm1F1Lw',
             'ChiI1Y0Dm1F2Lw', 'ChiI2Y3D1Lw', 'ChiI2Ym3Dm1Lw'}             # Dirac SDDM: U and V matrices
    code = f'from tests.test_conjugation import _sddm_fr; print(_sddm_fr({names!r}))'
    root = Path(__file__).resolve().parents[1]
    outs = {subprocess.run([sys.executable, '-c', code], cwd=root, capture_output=True, text=True,
                           check=True, env={**os.environ, 'PYTHONHASHSEED': seed}).stdout
            for seed in ('0', '1')}
    assert len(outs) == 1 and 'VMixDZ1Z2Z3' in outs.pop()


def test_parameters_carry_block_and_order_block():
    # the OrderBlock numbers rendered in the .fr are those of parameters()
    from tests.test_conjugation import _sddm_model
    frm = _sddm_model({'ChiI1Y0Lw', 'ChiI2Y3Lw', 'ChiI2Ym3Lw'})
    lagrangian = [(p.name, p.block, p.order_block, p.part) for p in frm.parameters()
                  if p.kind == 'External' and p.sector in ('Mass', 'Yukawa')]
    assert lagrangian == [
        ('MChiI1Y0', 'BSMMASS', 1, ''), ('MChiI2Y3ChiI2Ym3', 'BSMMASS', 2, ''),
        ('RyChiI1Y0ChiI2Y3', 'BSMYUKAWA', 1, 'R'), ('IyChiI1Y0ChiI2Y3', 'IMBSMYUKAWA', 1, 'I'),
        ('RyChiI1Y0ChiI2Ym3', 'BSMYUKAWA', 2, 'R'), ('IyChiI1Y0ChiI2Ym3', 'IMBSMYUKAWA', 2, 'I')]
