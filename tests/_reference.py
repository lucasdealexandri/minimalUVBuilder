"""Reference SDDM and stop multiplet models shared by the tests."""

import warnings

from tests import test_feynrules
from tests.test_conjugation import _sddm_fr

SDDM = {
    'sddm_doublet_partner': {'PhiI1Y0', 'PhiI2Y3'},
    'sddm_real_triplet': {'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'},
    'sddm_fermion_majorana': {'ChiI1Y0Lw', 'ChiI2Y3Lw', 'ChiI2Ym3Lw'},
    'sddm_fermion_dirac': {'ChiI1Y0D1F1Lw', 'ChiI1Y0D1F2Lw', 'ChiI1Y0Dm1F1Lw',
                           'ChiI1Y0Dm1F2Lw', 'ChiI2Y3D1Lw', 'ChiI2Ym3Dm1Lw'},
}
STOP = {
    'stop_weylstop_scalargluino': {'ChiC3I1Y4B1Lw', 'ChiC3bI1Ym4Bm1Lw', 'PhiC8I1Y0', 'PhiI1Y0'},
    'stop_diracstop_scalargluino': {'ChiC3I1Y4B1D1Lw', 'ChiC3bI1Ym4Bm1Dm1Lw', 'PhiC8I1Y0D1', 'PhiI1Y0D1'},
    'stop_scalarstop_majorana': {'ChiC8I1Y0Lw', 'ChiI1Y0Lw', 'PhiC3I1Y4B1'},
    'stop_scalarstop_dirac': {'ChiC8I1Y0D1Lw', 'ChiC8I1Y0Dm1Lw', 'ChiI1Y0D1Lw', 'ChiI1Y0Dm1Lw',
                              'PhiC3I1Y4B1D1'},
}


def _stop_model(names):
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.partition.partitions_builder import PartitionsBuilder
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.tools.reps import GlobalU1s
    builder = ModelBuilder(modelData=test_feynrules.TestRenames.STOP,
                           possible_dm_reps={'globalU1s': [GlobalU1s(B=0), GlobalU1s(B=0, D=1)]})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for pm in builder.generateParticleModels():
            for mm in PartitionsBuilder(pm).get_multipletModels():
                if mm.is_valid and {m.fr_symbol() for m in mm.non_conjugated_multiplets} == names:
                    return FeynRulesModel(mm, pm, mm.possible_mixings[0], 'stop.fr')
    raise AssertionError(f'no stop model with multiplets {names}')


def _stop_fr(names):
    return _stop_model(names).to_feynrules()
