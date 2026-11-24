'''MG5 parameter cards of a FeynRules model built by the UV builder.'''

from pathlib import Path
from typing import Dict, List, Optional

from uvbuilder.models.feynrulesmodel import FeynRulesModel
from uvbuilder.models.slhaBuilder import input_parameters, numeric_spectrum, write_card


class ParamCardGenerator:
    '''Fills the param_card of the model's UFO from its Lagrangian inputs: masses and
    mixing matrices of the BSM states (numeric_spectrum) and DECAY Auto (write_card).
    The model is read with FeynRulesModel.load from the pickle that `uvbuilder build`
    writes next to its .fr.'''

    def __init__(self, fr_model: FeynRulesModel):
        if not isinstance(fr_model, FeynRulesModel):
            raise TypeError(f'expected a FeynRulesModel, got {type(fr_model).__name__}')
        self.fr_model = fr_model
        self._spectrum = numeric_spectrum(fr_model)

    @property
    def inputs(self) -> List[str]:
        '''Names of the input parameters, as in the param_card.'''
        return [p.name for p in input_parameters(self.fr_model)]

    def spectrum(self, values: Dict[str, float]):
        '''(masses, mixings) for the input values (slhaBuilder.numeric_spectrum).'''
        return self._spectrum(values)

    def fill(self, card, values: Optional[Dict[str, float]] = None, output=None) -> Path:
        '''Fill card from its input blocks, or from values (then written into them);
        written to output, by default over card (slhaBuilder.write_card).'''
        return write_card(self.fr_model, card, output, values)
