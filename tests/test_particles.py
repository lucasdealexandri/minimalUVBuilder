#!/usr/bin/env python3

def test_SM():
    from uvbuilder.states.SMparticles import em

    assert em.isSM
    assert em.chargeConjugate().isSM

    # Flavour
    assert em.chargeConjugate().pdg < 0
    assert isinstance(em.chargeConjugate().pdgList, list)
    assert em.chargeConjugate().pdgList[0] < 0
