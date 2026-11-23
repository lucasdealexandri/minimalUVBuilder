"""Particle <-> field mixing and the mixing matrices (partition/particle_mixing.py)."""

import warnings

import pytest

from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
from uvbuilder.models.feynrulesmodel import FeynRulesModel, IndexTracker
from uvbuilder.models.multipletModel import default_component_label
from uvbuilder.partition.particle_mixing import MixingMatrix, ParticleMixing
from uvbuilder.tools.reps import DIRAC_SPINOR, GlobalU1s, U1Rep
from tests._sddm import DIRAC, DOUBLET_PARTNER, MAJORANA, REAL_TRIPLET, multiplet_model


def _matrices(names, labeler=default_component_label):
    mixing = multiplet_model(names, labeler).possible_mixings[0]
    return {m.standard_symbol: m for m in mixing.get_mixing_matrices()}


def test_mixing_values_are_sets_and_inverse_collects_particles():
    mixing = ParticleMixing()
    mixing['p'] = 'f'
    mixing['q'] = ['f', 'g']
    assert mixing['p'] == {'f'} and mixing['q'] == {'f', 'g'}
    assert dict(mixing.inverse) == {'f': {'p', 'q'}, 'g': {'q'}}


def test_majorana_states_mix_through_one_matrix():
    matrices = _matrices(MAJORANA)
    assert list(matrices) == ['UMixDW1', 'UMixMZ1Z2Z3', 'VMixDW1']       # sorted by name
    neutral = matrices['UMixMZ1Z2Z3']
    assert (neutral.type, neutral.dim, neutral.is_right_block) == ('MajoranaSpinor', 3, False)
    assert [p.label for p in neutral.particles] == ['Z1', 'Z2', 'Z3']
    assert [f.label for f in neutral.f_components] == ['ChiI1Y0q0', 'ChiI2Y3q0', 'ChiI2Ym3q0']
    assert [neutral.get_index(f) for f in neutral.f_components] == [0, 1, 2]
    assert neutral.get_index(matrices['UMixDW1'].f_components[0]) is None


def test_dirac_states_have_left_and_right_matrices():
    matrices = _matrices(DIRAC)
    left, right = matrices['UMixDZ1Z2Z3'], matrices['VMixDZ1Z2Z3']
    assert left.type == right.type == 'DiracSpinor'
    assert (left.is_right_block, right.is_right_block) == (False, True)
    assert [f.label for f in left.f_components] == ['ChiI1Y0D1q0F1', 'ChiI1Y0D1q0F2', 'ChiI2Y3D1q0']
    assert [f.label for f in right.f_components] == ['ChiI1Y0Dm1q0F1', 'ChiI1Y0Dm1q0F2', 'ChiI2Ym3Dm1q0']


def test_real_scalars_mix_through_one_real_matrix():
    neutral = _matrices(REAL_TRIPLET)['UMixRZ1Z2Z3']
    assert neutral.type == 'RealScalar'
    (parameter,) = neutral.to_feynrules(IndexTracker())['Parameters']
    options = dict(parameter.options)
    assert (parameter.name, parameter.kind) == ('UMixRZ1Z2Z3', 'External')
    assert options['ComplexParameter'] is False and options['Unitary'] is True


def test_complex_scalar_component_has_two_rows_in_a_real_block():
    neutral = _matrices(DOUBLET_PARTNER)['UMixRZ1Z2Z3']
    assert [f.label for f in neutral.f_components] == ['PhiI1Y0q0', 'PhiI2Y3q0', 'PhiI2Y3q0bar']


def test_complex_matrix_is_declared_by_real_and_imaginary_parts():
    real, imaginary, full = _matrices(MAJORANA)['UMixMZ1Z2Z3'].to_feynrules(IndexTracker())['Parameters']
    assert [(p.name, p.kind) for p in (real, imaginary, full)] == [
        ('RUMixMZ1Z2Z3', 'External'), ('IUMixMZ1Z2Z3', 'External'), ('UMixMZ1Z2Z3', 'Internal')]
    assert dict(real.options)['BlockName'] == 'UMIXMZ1Z2Z3'
    assert dict(imaginary.options)['BlockName'] == 'IMUMIXMZ1Z2Z3'
    # identity by default
    value = dict(real.options)['Value']
    assert all(f'RUMixMZ1Z2Z3[{i}, {j}] -> {int(i == j)}' in value
               for i in (1, 2, 3) for j in (1, 2, 3))
    assert '-> 1' not in dict(imaginary.options)['Value']
    assert dict(full.options)['Value'] == '{UMixMZ1Z2Z3[i_, j_] :> RUMixMZ1Z2Z3[i, j] + I*IUMixMZ1Z2Z3[i, j]}'


def test_single_state_needs_no_matrix():
    assert _matrices(MAJORANA)['UMixDW1'].to_feynrules(IndexTracker()) == {}


def test_matrix_must_be_square():
    neutral = _matrices(MAJORANA)['UMixMZ1Z2Z3']
    with pytest.raises(ValueError, match='Different number'):
        MixingMatrix(neutral.p_components[:2], neutral.f_components, neutral.particles)


def test_matrix_rows_follow_the_field_labels():
    '''The rows are sorted by field label: renaming the multiplets reorders them.'''
    def labeler(rep, i3=None, family_index=None):
        names = {(1, U1Rep(0)): 'B0', (2, U1Rep(1, 2)): 'Hu', (2, U1Rep(-1, 2)): 'Hd'}
        name = names.get((rep.isospin.dim, rep.hypercharge))
        if name is None or family_index is not None:        # multiplets of the other partitions
            return default_component_label(rep, i3, family_index)
        return name if i3 in (None, U1Rep(0)) else name + ('Up' if i3 == U1Rep(1, 2) else 'Down')

    hypercharges = lambda matrices: [f.rep.hypercharge for f in matrices['UMixMZ1Z2Z3'].f_components]
    assert hypercharges(_matrices(MAJORANA)) == [U1Rep(0), U1Rep(1, 2), U1Rep(-1, 2)]
    assert hypercharges(_matrices(MAJORANA, labeler)) == [U1Rep(0), U1Rep(-1, 2), U1Rep(1, 2)]


_PDG_NOT_MASS_ORDERED = {
    "masses": {1000022: 100., 1000023: 400., 1000024: 200., 1000025: 300.},
    "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
               1000025: {(1000022, 23): 1.0}},
    "labels": {1000022: 'Z1', 1000023: 'Z3', 1000024: 'W1', 1000025: 'Z2'}}


def test_family_index_is_the_mass_rank():
    # Z3 (1000023) is heavier than Z2 (1000025): the physical class and the
    # mixing columns are ordered by mass, not by PDG
    builder = ModelBuilder(modelData=_PDG_NOT_MASS_ORDERED,
                           possible_dm_reps={'spin': [DIRAC_SPINOR], 'globalU1s': [GlobalU1s()]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        (mm,) = builder.generateMultipletModels()
    mixing = mm.possible_mixings[0]
    neutral = next(m for m in mixing.get_mixing_matrices() if m.dim == 3)
    assert [p.pdg for p in neutral.particles] == [1000022, 1000025, 1000023]
    assert neutral.standard_symbol == 'UMixMZ1Z2Z3'
    fr = FeynRulesModel(mm, mm.particle_model, mixing, 'x.fr').to_feynrules()
    assert 'ClassMembers   -> {Z1,Z2,Z3}' in fr
    assert 'PDG            -> {1000022, 1000025, 1000023}' in fr
