"""Partitions of the particles into gauge multiplets (partition/partitions_builder.py)
and the labels of the resulting multiplet models (models/multipletModel.py)."""

import warnings

import pytest

from uvbuilder.models.multipletModel import default_component_label, multiplet_model_from_partition
from uvbuilder.partition.partitions_builder import PartitionsBuilder
from uvbuilder.tools.reps import U1Rep
from tests._sddm import builder


@pytest.fixture(scope='module')
def partitions():
    '''Particle model and partitions of the fermion singlet-doublet protomodel, by content.'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        (particle_model,) = builder().generateParticleModels()
        found = PartitionsBuilder(particle_model).get_partitions()
    content = lambda partition: tuple(sorted((rep.isospin.dim, str(rep.hypercharge)) for rep in partition))
    return particle_model, {content(partition): partition for partition in found}


SINGLETS = ((1, '-1'), (1, '0'), (1, '0'), (1, '0'), (1, '1'))
TRIPLET = ((1, '0'), (1, '0'), (3, '0'))
DOUBLETS = ((1, '0'), (2, '-1/2'), (2, '1/2'))


def test_partitions_of_the_singlet_doublet_protomodel(partitions):
    assert set(partitions[1]) == {SINGLETS, TRIPLET, DOUBLETS}


def test_partition_needs_the_conjugate_of_every_multiplet(partitions):
    particle_model, found = partitions
    check = PartitionsBuilder(particle_model).valid_partition
    singlet, down, up = sorted(found[DOUBLETS], key=lambda rep: rep.hypercharge + rep.isospin.dim)
    assert check([singlet]) and check([up, down]) and check([singlet, up, down])
    assert not check([up]) and not check([up, up, down])


def _labelled(partitions, content, labeler):
    particle_model, found = partitions
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return multiplet_model_from_partition(found[content], labeler=labeler,
                                              particle_model=particle_model)


def test_automatic_labels_are_accepted(partitions):
    for content in (SINGLETS, TRIPLET, DOUBLETS):
        assert _labelled(partitions, content, default_component_label) is not None


def test_label_of_a_standard_model_field_is_rejected(partitions):
    def labeler(rep, i3=None, family_index=None):
        if rep.isospin.dim == 1 and i3 is not None:
            return 'Phi'
        return default_component_label(rep, i3, family_index)
    with pytest.raises(ValueError, match='SM field labels'):
        _labelled(partitions, DOUBLETS, labeler)


def test_label_of_a_particle_is_rejected(partitions):
    def labeler(rep, i3=None, family_index=None):
        if rep.isospin.dim == 1 and i3 is not None:
            return 'Z1'
        return default_component_label(rep, i3, family_index)
    with pytest.raises(ValueError, match='particle labels'):
        _labelled(partitions, DOUBLETS, labeler)


def test_labeler_must_use_the_family_index(partitions):
    labeler = lambda rep, i3=None, family_index=None: default_component_label(rep, i3, None)
    with pytest.raises(ValueError, match='family_index'):
        _labelled(partitions, TRIPLET, labeler)


def test_two_multiplets_cannot_share_a_component_label(partitions):
    def labeler(rep, i3=None, family_index=None):
        if rep.isospin.dim == 2 and i3 is not None and rep.hypercharge + i3 == U1Rep(0):
            return 'Neutral'
        return default_component_label(rep, i3, family_index)
    with pytest.raises(ValueError, match='Duplicate'):
        _labelled(partitions, DOUBLETS, labeler)
