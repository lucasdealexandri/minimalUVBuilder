"""ParamCardGenerator (paramCard.py)."""

import pytest

from uvbuilder.models.feynrulesmodel import FeynRulesModel
from uvbuilder.models.slhaBuilder import write_card
from uvbuilder.paramCard import ParamCardGenerator
from tests.test_slhaBuilder import MAJORANA, SDDM_POINT, _sddm_model, _template


def test_generator_fills_like_write_card_after_a_round_trip(tmp_path):
    frm = _sddm_model(MAJORANA)
    reference = write_card(frm, _template(frm, tmp_path / 'template.dat'),
                           tmp_path / 'reference.dat', SDDM_POINT).read_text()
    generator = ParamCardGenerator(FeynRulesModel.load(frm.save2pickle(str(tmp_path / 'model.pkl'))))
    assert set(generator.inputs) == set(SDDM_POINT)
    card = generator.fill(tmp_path / 'template.dat', SDDM_POINT, tmp_path / 'card.dat')
    assert card.read_text() == reference


def test_generator_needs_a_feynrules_model():
    with pytest.raises(TypeError, match='expected a FeynRulesModel'):
        ParamCardGenerator({'masses': {}})
