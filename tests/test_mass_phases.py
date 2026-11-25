"""Phases of fermion mass eigenstates alone in their block (MultipletModel.mass_phases):
P = m*/|m| on the right-handed partner of a Dirac state, P = sqrt(m*/|m|) on a
Majorana state, with m the mass coefficient in the Lagrangian."""

import re

from tests._reference import SDDM, STOP, _sddm_fr, _stop_fr


def _phases(fr):
    return dict(re.findall(r'  (P\w+) == \{.*?Value\s+-> ([^\n]+),', fr, re.S))


def test_dirac_phase_is_trivial_for_a_positive_charged_mass():
    # + M SU2CG22[i1, i2] (Hd[i1].Hu[i2]): the charged pair has Dirac mass +M (mass_partner)
    fr = _sddm_fr(SDDM['sddm_fermion_majorana'])
    assert _phases(fr) == {'PW1': 'MChiI2Y3ChiI2Ym3/Abs[MChiI2Y3ChiI2Ym3]'}
    assert 'ChiI2Ym3Lw[s_, 2] :> PW1 W1Rwbar[s]' in fr
    assert 'ChiI2Y3Lw[s_, 1] :> W1Lw[s]' in fr                         # LH side unchanged
    assert fr.index('(* Mass-eigenstate phases: internal parameters *)') < fr.index('SU2CG22 == {')


def test_majorana_phase_is_square_root():
    fr = _stop_fr(STOP['stop_scalarstop_majorana'])
    assert _phases(fr) == {'PB1': 'Sqrt[MChiI1Y0/Abs[MChiI1Y0]]',
                           'PGo': 'Sqrt[MChiC8I1Y0/Abs[MChiC8I1Y0]]'}
    assert 'ChiI1Y0Lw[s_] :> PB1 B1Lw[s]' in fr and 'ChiC8I1Y0Lw[s_, c_] :> PGo GoLw[s, c]' in fr


def test_phase_conjugated_in_rh_conjugate_class():
    fr = _stop_fr(STOP['stop_weylstop_scalargluino'])
    assert 'ChiC3I1Y4B1Rw[s_, c_] :> Conjugate[Pst] stRw[s, c]' in fr


def test_no_phase_for_scalars_or_mixing_blocks():
    assert _phases(_sddm_fr(SDDM['sddm_doublet_partner'])) == {}
    assert 'PZ' not in ''.join(_phases(_sddm_fr(SDDM['sddm_fermion_majorana'])))


def test_pure_wino_neutral_phase_is_trivial_for_a_positive_mass():
    # Majorana triplet: + 1/2 M (X0 X0 - 2 X+ X-) (Multiplet.mass_partner): the neutral state
    # has the mass +M, the charged pair the sign SU(2) fixes, -M (Dirac phase -1)
    import warnings
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.tools.reps import DIRAC_SPINOR
    wino = {'masses': {1000022: 500., 1000024: 500.2},
            'decays': {1000022: {}, 1000024: {(1000022, 24): 1.0}},
            'labels': {1000022: 'N1', 1000024: 'C1'}}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mm = next(mm for mm in ModelBuilder(modelData=wino, possible_dm_reps={'spin': [DIRAC_SPINOR]})
                  .generateMultipletModels()
                  if {m.fr_symbol() for m in mm.non_conjugated_multiplets} == {'ChiI3Y0Lw'})
        fr = FeynRulesModel(mm, mm.particle_model, mm.possible_mixings[0], 'wino.fr').to_feynrules()
    assert _phases(fr) == {'PN1': 'Sqrt[MChiI3Y0/Abs[MChiI3Y0]]', 'PC1': '-MChiI3Y0/Abs[MChiI3Y0]'}
