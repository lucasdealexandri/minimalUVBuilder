"""Singlet-doublet protomodel shared by the tests."""

import warnings
from functools import lru_cache

from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
from uvbuilder.models.multipletModel import default_component_label
from uvbuilder.tools.reps import DIRAC_SPINOR, SCALAR, GlobalU1s

PROTOMODEL = {"masses": {1000022: 100., 1000023: 300., 1000024: 200., 1000025: 400.},
              "decays": {1000022: {}, 1000023: {(1000022, 23): 1.0}, 1000024: {(1000022, 24): 1.0},
                         1000025: {(1000023, 25): 1.0}},
              "labels": {1000022: 'Z1', 1000023: 'Z2', 1000024: 'W1', 1000025: 'Z3'}}

MAJORANA = frozenset({'ChiI1Y0', 'ChiI2Y3', 'ChiI2Ym3'})
DIRAC = frozenset({'ChiI1Y0D1F1', 'ChiI1Y0D1F2', 'ChiI1Y0Dm1F1', 'ChiI1Y0Dm1F2',
                   'ChiI2Y3D1', 'ChiI2Ym3Dm1'})
REAL_TRIPLET = frozenset({'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'})
DOUBLET_PARTNER = frozenset({'PhiI1Y0', 'PhiI2Y3'})


def builder(spin=DIRAC_SPINOR, charges=GlobalU1s()):
    return ModelBuilder(modelData=PROTOMODEL,
                        possible_dm_reps={'spin': [spin], 'globalU1s': [charges]})


@lru_cache(maxsize=None)
def multiplet_model(names, labeler=default_component_label):
    """The valid multiplet model whose (automatic) multiplet labels are `names`."""
    spin = SCALAR if names in (REAL_TRIPLET, DOUBLET_PARTNER) else DIRAC_SPINOR
    charges = GlobalU1s(D=1) if names == DIRAC else GlobalU1s()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        automatic = list(builder(spin, charges).generateMultipletModels())
        position = [frozenset(m.label for m in mm.non_conjugated_multiplets)
                    for mm in automatic].index(names)
        if labeler is default_component_label:
            return automatic[position]
        return list(builder(spin, charges).generateMultipletModels(labeler=labeler))[position]
