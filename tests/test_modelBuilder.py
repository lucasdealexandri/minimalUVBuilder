#!/usr/bin/env python3

import pytest

from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
from uvbuilder.models.feynrulesmodel import FeynRulesModel
from uvbuilder.models.finalStateModel import FinalStateModel
from uvbuilder.models.particleModel import ParticleModel
from uvbuilder.states.finalState import FinalState
from uvbuilder.states.particle import Particle
from uvbuilder.states.physicsUnits import GeV
from uvbuilder.tools.reps import DIRAC_SPINOR, GlobalU1s
from uvbuilder.tools.reps import SCALAR, SU3Rep, U1Rep, UnbrokenRep
from uvbuilder.vertices.particleVertex import ParticleVertex


def _final_state_model_with_dm():
    model = FinalStateModel(vertices=[])
    dm = Particle(pdg=1000022, label='DM')
    dm.mass = 100*GeV   # every protomodel particle has a mass
    model.pdgDict[1000022] = FinalState(label='DM', particleList=[dm])
    return model


def test_dm_particles_use_complete_unbroken_reps():
    builder = ModelBuilder()
    particles = builder.getDMParticles(_final_state_model_with_dm())

    assert len(particles) == 2       # 2 spins, no global charge
    assert all(isinstance(p.rep, UnbrokenRep) for p in particles)
    assert {p.rep.spin for p in particles} == {SCALAR, builder.possible_dm_reps['spin'][1]}
    assert {p.rep.globalU1s for p in particles} == {GlobalU1s()}
    assert all(p.selfConjugate == (p.rep.conjRep == p.rep) for p in particles)


# -- global U(1)s: which charges are conserved is decided by the DM candidate --

_STOP = {"masses": {1000022: 180.603, 1000006: 1164.264, 1000021: 1500.0},
         "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}, 1000021: {(1000006, -6): 1.0}},
         "labels": {1000022: 'B', 1000006: 'st', 1000021: 'Go'}}


def _bsm_charges(dm_options):
    from uvbuilder.tools.reps import GlobalU1s  # noqa: F401  (documents the type in use)
    builder = ModelBuilder(modelData=_STOP, possible_dm_reps={'globalU1s': dm_options})
    return [{p.label: p.globalU1s for p in pm.bsm_particles}
            for pm in builder.generateParticleModels()]


def test_nothing_conserved_fills_nothing():
    from uvbuilder.tools.reps import GlobalU1s
    models = _bsm_charges([GlobalU1s()])
    assert models
    assert all(q == GlobalU1s() for charges in models for q in charges.values())


def test_baryon_number_conserved_is_propagated():
    from uvbuilder.tools.reps import GlobalU1s
    models = _bsm_charges([GlobalU1s(B=0)])        # B set (to 0) on the DM: B conserved
    assert models
    for charges in models:
        assert charges['B'] == GlobalU1s(B=0)
        assert charges['st'] == GlobalU1s(B=U1Rep(1, 3))   # st -> t B
        assert charges['Go'] == GlobalU1s(B=0)             # Go -> st t~


def test_untracked_charges_do_not_leak_from_sm():
    from uvbuilder.tools.reps import GlobalU1s
    # only L conserved: the top's B = 1/3 must not reach the stop
    for charges in _bsm_charges([GlobalU1s(L=0)]):
        assert all(q == GlobalU1s(L=0) for q in charges.values())


def test_neutral_decay_to_charged_bsm_particle():
    # Z2 -> W1+ W- together with W1+ -> Z1 W+: both W signs must be kept
    data = {"masses": {1000022: 100., 1000023: 300., 1000024: 200.},
            "decays": {1000022: {},
                       1000023: {(1000022, 23): 0.9, (1000024, -24): 0.1},
                       1000024: {(1000022, 24): 1.0}},
            "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1'}}
    models = [{p.label: p.eCharge for p in pm.bsm_particles}
              for pm in ModelBuilder(modelData=data).generateParticleModels()]
    assert models
    assert all(charges == {'Z1': 0, 'Z2': 0, 'W1': 1} for charges in models)


# -- multiplet and FeynRules models straight from the builder --

def _sddm_builder():
    data = {"masses": {1000022: 100., 1000023: 300., 1000024: 200., 1000025: 400.},
            "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                       1000025: {(1000023, 25): 1.0}},
            "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}
    return ModelBuilder(modelData=data)


def _multiplet_symbols(multipletModel):
    return sorted(m.fr_symbol() for m in multipletModel.non_conjugated_multiplets)


def test_generate_multiplet_models_matches_partitions_builder():
    import warnings
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    builder = _sddm_builder()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        expected = [mm for pm in builder.generateParticleModels()
                    for mm in PartitionsBuilder(pm).get_multipletModels()]
        valid = [_multiplet_symbols(mm) for mm in expected if mm.is_valid]
        generated = list(builder.generateMultipletModels())
        unfiltered = list(builder.generateMultipletModels(validOnly=False))
    assert valid and len(valid) < len(expected)
    assert [_multiplet_symbols(mm) for mm in generated] == valid
    assert [_multiplet_symbols(mm) for mm in unfiltered] == [_multiplet_symbols(mm) for mm in expected]
    assert all(mm.particle_model is not None for mm in generated)


def test_labeler_is_applied_to_valid_models_only():
    from uvbuilder.models.multipletModel import default_component_label
    def labeler(rep, i3=None, family_index=None):
        # charged singlets appear only in the invalid partitions
        if rep.isospin.dim == 1 and rep.hypercharge != U1Rep(0):
            raise ValueError(f'no label for {rep}')
        return 'S' + default_component_label(rep, i3, family_index)
    builder = _sddm_builder()
    valid = list(builder.generateMultipletModels(labeler=labeler))
    unfiltered = list(builder.generateMultipletModels(labeler=labeler, validOnly=False))
    assert valid and all(s.startswith('S') for mm in valid for s in _multiplet_symbols(mm))
    assert len(unfiltered) > len(valid)
    assert [_multiplet_symbols(mm) for mm in unfiltered if mm.is_valid] == \
           [_multiplet_symbols(mm) for mm in valid]


def test_decays_need_a_direct_coupling():
    # colored MSSM: b1 -> t1 W needs the W coupling of a squark doublet; SU(2)-singlet
    # squarks are joined only through the gluino Yukawas (loop level)
    T1, T2, B1, GO, DM = 1000006, 2000006, 1000005, 1000021, 1000022
    data = {'masses': {DM: 100., T1: 540., T2: 790., B1: 750., GO: 1400.},
            'decays': {DM: {}, T1: {(DM, 6): 1.0},
                       T2: {(T1, 25): 0.4, (T1, 23): 0.3, (DM, 6): 0.3},
                       B1: {(DM, 5): 0.3, (T1, -24): 0.7},
                       GO: {(T1, -6): 0.5, (T2, -6): 0.3, (B1, -5): 0.2}},
            'labels': {DM: 'X1', T1: 'T1', T2: 'T2', B1: 'S1', GO: 'G1'}}
    builder = ModelBuilder(modelData=data, possible_dm_reps={'globalU1s': [GlobalU1s(B=0)]})
    scalars = [sorted((m.rep.isospin.dim, m.rep.hypercharge) for m in mm.non_conjugated_multiplets
                      if m.rep.spin == SCALAR and m.rep.hypercharge != U1Rep(0))
               for mm in builder.generateMultipletModels()]
    assert [(1, U1Rep(2, 3)), (2, U1Rep(1, 6))] in scalars
    assert all(any(dim == 2 for dim, _ in s) for s in scalars if s)


def test_mixing_columns_follow_the_mass_order():
    data = {"masses": {1000022: 100., 1000023: 400., 1000024: 200., 1000025: 300.},
            "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                       1000025: {(1000022, 25): 1.0}},
            "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}
    multipletModel = next(ModelBuilder(modelData=data).generateMultipletModels())
    columns = {m.standard_symbol: [p.pdg for p in m.particles]
               for m in multipletModel.possible_mixings[0].get_mixing_matrices()}
    assert columns['UMixRZ1Z3Z2'] == [1000022, 1000025, 1000023]


def test_feynrules_model_survives_pickling(tmp_path):
    import pickle, subprocess, sys
    model = next(_sddm_builder().generateFeynRulesModels())
    path = tmp_path / 'model.pkl'
    path.write_bytes(pickle.dumps(model))
    # a fresh interpreter, so no FeynRules index is cached
    loaded = subprocess.run([sys.executable, '-c',
                             f'import pickle; print(pickle.load(open({str(path)!r}, "rb")).to_feynrules(), end="")'],
                            capture_output=True, text=True, check=True).stdout
    assert loaded == model.to_feynrules()


def test_feynrules_model_saved_next_to_its_fr_and_loaded(tmp_path):
    import pickle
    model = next(_sddm_builder().generateFeynRulesModels(outputFile=str(tmp_path / 'sddm.fr')))
    path = model.save2pickle()
    assert path == model.output_file.removesuffix('.fr') + '.pkl'
    assert FeynRulesModel.load(path).to_feynrules() == model.to_feynrules()
    (tmp_path / 'other.pkl').write_bytes(pickle.dumps({'masses': {}}))
    with pytest.raises(TypeError, match='holds a dict, not a FeynRulesModel'):
        FeynRulesModel.load(tmp_path / 'other.pkl')


def test_generate_feynrules_models():
    import warnings
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    builder = _sddm_builder()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        multipletModels = list(builder.generateMultipletModels())
        models = list(builder.generateFeynRulesModels(outputFile='out/sddm.fr'))
        first = multipletModels[0]
        direct = FeynRulesModel(first, first.particle_model, first.possible_mixings[0], 'direct.fr')
        assert models[0].to_feynrules() == direct.to_feynrules()
    assert len(models) == sum(len(mm.possible_mixings) for mm in multipletModels)
    assert [m.output_file for m in models] == [f'out/sddm_{n}.fr' for n in range(1, len(models) + 1)]


def test_feynrules_model_is_saved_to_its_output_file(tmp_path):
    import warnings
    target = tmp_path / 'models' / 'sddm.fr'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = next(_sddm_builder().generateFeynRulesModels(outputFile=str(target)))
        model.save2file()
        written = tmp_path / 'models' / 'sddm_1.fr'
        assert written.read_text() == model.to_feynrules()


# -- operators without an implemented color contraction are dropped --

_M368 = {"masses": {1000022: 100.0, 1000003: 400.0, 1000021: 800.0, 9000006: 1200.0},
         "decays": {1000022: {}, 1000003: {(1000022, 6): 1.0}, 1000021: {(1000003, -6): 1.0},
                    9000006: {(1000021, -6): 1.0}},
         "labels": {1000022: 'X1', 1000003: 'S3', 1000021: 'O8', 9000006: 'P6'}}


def test_unimplemented_color_contraction_drops_operator_and_model():
    '''P6 -> O8 tbar allows P6 in 3bar, 6 or 15. No color contraction exists for (15, 8, 3bar),
    so the 15 models lose their decay operator and are invalid; the 6 model keeps its Jbar operator.'''
    import warnings
    from uvbuilder.tools.reps import DIRAC_SPINOR
    builder = ModelBuilder(modelData=_M368,
                           possible_dm_reps={'spin': [DIRAC_SPINOR, SCALAR], 'eCharge': U1Rep(0),
                                             'colordim': SU3Rep(0, 0), 'globalU1s': [GlobalU1s(B=0)]})
    with pytest.warns(UserWarning, match='not implemented; operator dropped'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            warnings.filterwarnings('always', message='.*not implemented; operator dropped')
            models = list(builder.generateMultipletModels())
    contents = [{m.fr_symbol() for m in mm.non_conjugated_multiplets} for mm in models]
    assert not any('C15' in name for names in contents for name in names)
    target = {'ChiC8I1Y0Lw', 'ChiI1Y0Lw', 'PhiC3I1Y4B1', 'PhiC6I1Ym4Bm1'}
    assert target in contents
    m368 = models[contents.index(target)]
    def has(op, prefix):    # the sextet may enter as PhiC6I1Ym4Bm1 or as its conjugate
        return any(m.fr_symbol().startswith(prefix) for m in op.multiplets)
    sextet_octet_ops = [op for op in m368.operators if has(op, 'PhiC6I1Ym4Bm1') and has(op, 'ChiC8I1Y0')]
    assert len(sextet_octet_ops) == 3                     # u, c, t
    assert len({op.su3_tensor.name for op in sextet_octet_ops} & {'J', 'Jbar'}) == 1   # 3-8-6 head


def test_custom_name_author_and_version():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        models = list(_sddm_builder().generateFeynRulesModels(name='SDDM', author='Me',
                                                              version='2.0'))
    assert len(models) > 1
    assert [m.name for m in models] == [f'SDDM_{n}' for n in range(1, len(models) + 1)]
    assert [m.output_file for m in models] == [f'SDDM_{n}.fr' for n in range(1, len(models) + 1)]
    fr = models[0].to_feynrules()
    assert 'M$ModelName = "SDDM_1";' in fr
    assert 'Authors -> {"Me"}' in fr and 'Version -> "2.0"' in fr


def test_single_model_is_unnumbered():
    import warnings
    from uvbuilder.tools.reps import DIRAC_SPINOR
    from tests.test_particle_mixing import _PDG_NOT_MASS_ORDERED
    builder = ModelBuilder(modelData=_PDG_NOT_MASS_ORDERED,
                           possible_dm_reps={'spin': [DIRAC_SPINOR], 'globalU1s': [GlobalU1s()]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        (model,) = builder.generateFeynRulesModels(name='Single')
    assert (model.name, model.output_file) == ('Single', 'Single.fr')


def test_invalid_model_name_is_rejected():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(ValueError, match='invalid model name'):
            next(_sddm_builder().generateFeynRulesModels(name='3-6-8'))


# -- electric charges from the decay vertices --

def _charges(decays, masses, globalU1s=GlobalU1s(), labels=None):
    data = {'masses': masses, 'decays': decays, **({'labels': labels} if labels else {})}
    builder = ModelBuilder(modelData=data, possible_dm_reps={'spin': [DIRAC_SPINOR],
                                                             'globalU1s': [globalU1s]})
    return [{p.pdg: p.eCharge for p in pm.bsm_particles} for pm in builder.generateParticleModels()]


_EWINO_MASSES = {1000022: 100., 1000024: 200., 1000023: 300.}


@pytest.mark.parametrize('z2_decays', [
    {(-1000024, 24): 1.0},                           # conjugated BSM daughter
    {(-1000024, 24): 0.5, (1000024, -24): 0.5},      # both charge states
])
def test_charge_of_a_neutral_mother_of_a_conjugated_daughter(z2_decays):
    charges = _charges({1000022: {}, 1000024: {(1000022, 24): 1.0}, 1000023: z2_decays},
                       _EWINO_MASSES)
    assert charges == [{1000022: U1Rep(0), 1000024: U1Rep(1), 1000023: U1Rep(0)}]


def test_charge_of_a_doubly_charged_mother():
    charges = _charges({1000022: {}, 1000024: {(1000022, 24): 1.0}, 1000023: {(1000024, 24): 1.0}},
                       _EWINO_MASSES)
    assert charges == [{1000022: U1Rep(0), 1000024: U1Rep(1), 1000023: U1Rep(2)}]


def test_charges_with_standard_pdgs_not_ordered_by_mass():
    # gluino -> stop2 tbar with stop2 = 2000006 (Z2 vertex visited after a copy of the model)
    charges = _charges({1000022: {}, 2000006: {(1000022, 6): 1.0}, 1000021: {(2000006, -6): 1.0}},
                       {1000022: 100., 2000006: 800., 1000021: 1400.}, globalU1s=GlobalU1s(B=0))
    assert charges and all(c == {1000022: U1Rep(0), 2000006: U1Rep(2, 3), 1000021: U1Rep(0)}
                           for c in charges)


def test_charges_of_a_two_step_chain_against_pdg_order():
    charges = _charges({1000022: {}, 1000025: {(1000022, 23): 1.0}, 1000023: {(1000025, 25): 1.0}},
                       {1000022: 100., 1000025: 200., 1000023: 400.})
    assert charges == [{1000022: U1Rep(0), 1000025: U1Rep(0), 1000023: U1Rep(0)}]


def test_charges_of_a_three_step_chain_against_pdg_order():
    charges = _charges({1000022: {}, 1000035: {(1000022, 25): 1.0},
                        1000025: {(1000035, 25): 1.0}, 1000023: {(1000025, 25): 1.0}},
                       {1000022: 100., 1000035: 200., 1000025: 300., 1000023: 400.})
    assert charges == [{pdg: U1Rep(0) for pdg in (1000022, 1000035, 1000025, 1000023)}]


def test_conjugated_daughters_give_the_same_multiplet_models():
    def models(z1, z3):
        data = {'masses': {1000022: 109.3, 1000023: 289.8, 1000024: 200., 1000025: 430.6},
                'decays': {1000022: {}, 1000024: {(1000022, 24): 1.0},
                           1000023: {(z1, 23): 0.8, (1000022, 25): 0.2}, 1000025: z3}}
        builder = ModelBuilder(modelData=data, possible_dm_reps={'spin': [DIRAC_SPINOR],
                                                                 'globalU1s': [GlobalU1s()]})
        return sorted(str(m) for m in builder.generateMultipletModels())
    reference = models(1000022, {(1000024, -24): 1.0})
    assert len(reference) == 1
    assert models(-1000022, {(1000024, -24): 1.0}) == reference                         # neutral
    assert models(1000022, {(-1000024, 24): 0.5, (1000024, -24): 0.5}) == reference    # charged


def test_states_without_labels_are_named_by_their_pdg():
    import warnings
    data = {'masses': {1000022: 100., 1000024: 200.},
            'decays': {1000022: {}, 1000024: {(1000022, 24): 1.0}}, 'labels': {1000024: 'W1'}}
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        model = FinalStateModel.fromDict(data)
    labels = {p.pdg: p.label for fs in model.pdgDict.values() for p in fs.getParticles()
              if not p.isSM and p.pdg > 0}
    assert labels == {1000022: 'X1000022', 1000024: 'W1'}


# -- the SM particle of a decay must be emitted by a tree-level coupling --

def _valid(protomodel, **options):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return {frozenset(m.fr_symbol() for m in mm.non_conjugated_multiplets)
                for mm in ModelBuilder(modelData=protomodel, **options).generateMultipletModels()}


# inert doublet: H0 (DM), A0, H+
_IDM = {'masses': {1000022: 100., 1000023: 200., 1000024: 250.},
        'decays': {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0}}}


def test_inert_doublet_and_no_singlet_triplet():
    # A0 -> H0 Z: the singlet and the neutral triplet component have T3 = 0, no Z coupling
    assert _valid(_IDM) == {frozenset({'PhiI2Y3'}),
                                        frozenset({'ChiI2Y3Lw', 'ChiI2Ym3Lw'})}


def test_no_model_emits_a_photon():
    # (a gluon is read as a jet, which a colored state emits as a quark)
    data = {'masses': {1000022: 100., 1000023: 300.},
            'decays': {1000022: {}, 1000023: {(1000022, 22): 1.0}}}
    assert _valid(data) == set()


def test_real_triplet_emits_h_but_not_z():
    from tests._sddm import PROTOMODEL, PROTOMODEL_H, REAL_TRIPLET
    triplet = frozenset(f'{name}' for name in REAL_TRIPLET)
    assert triplet not in _valid(PROTOMODEL, possible_dm_reps={'spin': [SCALAR]})
    assert triplet in _valid(PROTOMODEL_H, possible_dm_reps={'spin': [SCALAR]})


def test_three_body_decays_are_rejected_by_name():
    data = {'masses': {1000022: 90., 1000023: 126., 1000024: 155.},
            'decays': {1000022: {}, 1000023: {(1000022, 12, 12): 1.0},
                       1000024: {(1000022, 2, 1): 1.0}}}
    with pytest.raises(ValueError, match='only two-body decays.*: '
                                         '1000023 -> 1000022 12 12, 1000024 -> 1000022 2 1$'):
        FinalStateModel.fromDict(data)


def test_sm_decays_are_skipped():
    data = {'masses': {1000022: 100., 1000024: 200.},
            'decays': {1000022: {}, 1000024: {(1000022, 24): 1.0}, 6: {(5, 24): 1.0}}}
    assert len(list(ModelBuilder(modelData=data).generateMultipletModels())) > 0


def test_protomodel_labels_and_masses_do_not_change_sm_particles():
    from uvbuilder.states.SMparticles import SMparticlesDict
    w = SMparticlesDict[24]
    label, mass = w.label, w.mass
    _charges({1000022: {}, 1000024: {(1000022, 24): 1.0}}, {1000022: 100., 1000024: 200., 24: 1.0},
             labels={24: 'Wboson'})
    assert (w.label, w.mass) == (label, mass)
