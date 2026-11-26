#!/usr/bin/env python3

def test_SM():
    from uvbuilder.states.SMparticles import em

    assert em.isSM
    assert em.antiParticle.isSM

    # Flavor
    assert em.antiParticle.pdg < 0
    assert isinstance(em.antiParticle.pdgList, list)
    assert em.antiParticle.pdgList[0] < 0


def test_state_dict_stores_the_particle_and_finds_the_antiparticle():
    '''A StateDict keeps the particle (positive PDG, its label) whatever is stored, and
    returns the antiparticle for a negative PDG.'''
    import pytest
    from uvbuilder.states.SMparticles import em
    from uvbuilder.states.stateDict import StateDict
    states = StateDict({-11: em.antiParticle})
    states['e+'] = em.antiParticle
    assert list(dict.keys(states)) == [11, 'e-'] and dict.__getitem__(states, 11) is em
    assert states[-11] is em.antiParticle and states.get(-11) is em.antiParticle
    with pytest.raises(ValueError, match='PDGs'):
        states[1.5] = em
