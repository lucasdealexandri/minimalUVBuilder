"""Phases of fermion mass eigenstates alone in their block (MultipletModel.mass_phases):
P = m*/|m| on the right-handed partner of a Dirac state, P = sqrt(m*/|m|) on a
Majorana state, with m the mass coefficient in the Lagrangian."""

import re

from tests.test_golden import SDDM, STOP, _sddm_fr, _stop_fr


def _phases(fr):
    return dict(re.findall(r'  (P\w+) == \{.*?Value\s+-> ([^\n]+),', fr, re.S))


def test_dirac_phase_follows_su2_sign_of_the_charged_mass():
    # -M SU2CG22[i1, i2] (Hd[i1].Hu[i2]): charged part +M (Hd-.Hu+), Dirac mass -M
    fr = _sddm_fr(SDDM['sddm_fermion_majorana'])
    assert _phases(fr) == {'PW1': '-Conjugate[MChiI2Y3ChiI2Ym3]/Abs[MChiI2Y3ChiI2Ym3]'}
    assert 'ChiI2Ym3Lw[s_, 2] :> PW1 W1Rwbar[s]' in fr
    assert 'ChiI2Y3Lw[s_, 1] :> W1Lw[s]' in fr                         # LH side unchanged
    assert fr.index('(* Mass-eigenstate phases: internal parameters *)') < fr.index('SU2CG22 == {')


def test_majorana_phase_is_square_root():
    fr = _stop_fr(STOP['stop_scalarstop_majorana'])
    assert _phases(fr) == {'PB1': '(Conjugate[MChiI1Y0]/Abs[MChiI1Y0])^(1/2)',
                           'PGo': '(Conjugate[MChiC8I1Y0]/Abs[MChiC8I1Y0])^(1/2)'}
    assert 'ChiI1Y0Lw[s_] :> PB1 B1Lw[s]' in fr and 'ChiC8I1Y0Lw[s_, c_] :> PGo GoLw[s, c]' in fr


def test_phase_conjugated_in_rh_conjugate_class():
    fr = _stop_fr(STOP['stop_weylstop_scalargluino'])
    assert 'ChiC3I1Y4B1Rw[s_, c_] :> Conjugate[Pst] stRw[s, c]' in fr


def test_no_phase_for_scalars_or_mixing_blocks():
    assert _phases(_sddm_fr(SDDM['sddm_doublet_partner'])) == {}
    assert 'PZ' not in ''.join(_phases(_sddm_fr(SDDM['sddm_fermion_majorana'])))
