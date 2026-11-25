'''Components of a multiplet in the FeynRules output: scalar and Weyl components
(ScalarComponent, WeylComponent) and their class definitions.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
'''
from __future__ import annotations
from typing import Union, Sequence, List, Optional, Dict
from .reps import Scalar, WeylSpinor, WEYL_LEFT, WEYL_RIGHT, U1Rep
from .reps import UnbrokenRep, BrokenRep
from uvbuilder.tools.fr_utils import (formatClassBlock, fr_quantum_numbers, LEFT_SUFFIX, RIGHT_SUFFIX,
                            DAGGER_SUFFIX, toggle_dagger)


def _fr_charge(charge: Optional[U1Rep]) -> str:
    '''Render a U(1) charge as Mathematica text, e.g. U1Rep(2,3) -> "2/3".'''
    if charge is None:
        return '0'
    if charge.denominator == 1:
        return str(charge.numerator)
    return f'{charge.numerator}/{charge.denominator}'


class ScalarComponent:

    def __init__(self, label: str, rep: Union[UnbrokenRep, BrokenRep]):
        assert isinstance(rep.spin, Scalar), f'spin should be {Scalar} instead of {type(rep.spin).__name__}.'
        self.rep: Union[UnbrokenRep, BrokenRep] = rep
        self.label: str = label
        self.spin = rep.spin

    @property
    def conjugate_label(self) -> str:

        # For selfConjugate Internal Reps, do not change the label
        if self.rep.is_realInternalRep:
            return self.label

        # If the label carries the conjugation marker, remove it
        if self.label.endswith(DAGGER_SUFFIX):
            # Assures conjugate_label operator is an involution.
            return toggle_dagger(self.label)

        charge_flip = str.maketrans({'+': '-', '-': '+'})
        label = self.label.translate(charge_flip)

        # If there is a sign in the label, flip the sign
        if label != self.label:
            return label

        # If none of the above, simply add the conjugation marker at the end
        return toggle_dagger(self.label)

    @property
    def dagger(self) -> ScalarComponent:

        return ScalarComponent(self.conjugate_label, self.rep.conjRep)

    def scalar_product(
            self, 
            other: Union[Sequence[ScalarComponent], ScalarComponent]) -> List[ScalarComponent]:

        if isinstance(other, Sequence):
            final_list = []
            for component in other:
                final_list += self * component
            return final_list

        elif isinstance(other, ScalarComponent):
            reps = self.rep * other.rep
            return [
                ScalarComponent(f'{self.label} {other.label}', rep) for rep in reps
            ]

    def __mul__(self, other: Union[Sequence[ScalarComponent], ScalarComponent]):

        return self.scalar_product(other)

    def __repr__(self) -> str:
        return f'{self.label}'

    def __eq__(self, other: ScalarComponent) -> bool:
        return self.label == other.label and self.rep == other.rep

    def __lt__(self, other: ScalarComponent) -> bool:
        return self.label.__lt__(other.label)

    def __le__(self, other: ScalarComponent) -> bool:
        return self.label.__le__(other.label)

    def __gt__(self, other: ScalarComponent) -> bool:
        return self.label.__gt__(other.label)

    def __ge__(self, other: ScalarComponent) -> bool:
        return self.label.__ge__(other.label)

    def __hash__(self):
        return hash((self.label, self.rep))


class WeylComponent:
    '''
    Class for Weyl spinors. Left Weyls should represent chi, eta or xi,
    while Right Weyls chi^dagger, eta^dagger or xi^dagger (rendered with a
    trailing ``bar`` in FeynRules, e.g. chibar).
    '''
    def __init__(self, label: str, rep: Union[UnbrokenRep, BrokenRep]):
        assert isinstance(rep.spin, WeylSpinor), f'spin should be WeylSpinor instead of {type(rep.spin).__name__}.'
        self.rep: Union[UnbrokenRep, BrokenRep] = rep
        self.label: str = label
        self.chirality: str = "LH" if self.rep.spin == WEYL_LEFT else "RH"
        self.spin = rep.spin

    @property
    def conjugate_label(self) -> str:
        # Involution: X <-> Xbar
        return toggle_dagger(self.label)

    @property
    def dagger(self) -> WeylComponent:
        return WeylComponent(self.conjugate_label, self.rep.conjRep)

    def scalar_product(self, other: WeylComponent):
        '''
        Intepreted as chi1_a chi2^a for two left handed Weyls,
        or chi1^adot chi2_adot for two right handed Weyls. 
        '''
        if self.rep.spin != other.rep.spin:
            raise ValueError(f'Both Weyls should either be LH or RH.')
        reps = self.rep * other.rep
        label = f'{self.label} {other.label}'
        return [ScalarComponent(label, rep) for rep in reps]

    def describe(self) -> str:
        return f'{self.chirality} Weyl {self.label} of rep {self.rep}'

    def to_feynrules(self, tracker, family=None) -> Dict[str, List[str]]:
        '''FeynRules class declaration for this Weyl component, as an FRChunk::

            W[3] == {
                ClassName -> uqLw,
                Unphysical -> True,
                Chirality -> Left,
                SelfConjugate -> False,
                Indices -> {Index[Colour], Index[Gen]},
                FlavorIndex -> Gen,
                QuantumNumbers -> {Q -> 2/3}
            }

        The index ``tracker`` provides the class number and the gauge
        indices; ``family`` (an optional FRIndex) adds the generation index.
        '''
        header = f"W[{tracker.class_number('W')}]"
        chirality = 'Left' if self.chirality == 'LH' else 'Right'

        blockDict = {
            'ClassName': self.label,
            'Unphysical': 'True',
            'Chirality': chirality,
            'SelfConjugate': 'False',
        }
        indices = [str(index) for index in tracker.gauge(self.rep)]
        if family is not None:
            indices.append(str(family))
        if indices:
            blockDict['Indices'] = f"{{{','.join(indices)}}}"
        if family is not None:
            blockDict['FlavorIndex'] = family.name
        charge = getattr(self.rep, 'eCharge', None)
        blockDict['QuantumNumbers'] = fr_quantum_numbers(
            {'Q': charge if charge is not None else U1Rep(0)}, getattr(self.rep, 'globalU1s', None))

        return {'Classes': [formatClassBlock(header, blockDict)]}

    @property
    def left_hand(self) -> WeylComponent:
        '''Returns LH itself or RHbar'''

        if self.rep.spin == WEYL_LEFT: return self
        elif self.rep.spin == WEYL_RIGHT:
            if not self.rep.is_realInternalRep:
                return self.dagger
            else: # Majorana Case
                if self.label.endswith(RIGHT_SUFFIX):
                    return WeylComponent(self.label.removesuffix(RIGHT_SUFFIX) + LEFT_SUFFIX, self.dagger.rep)
                else: raise NameError(f"Label should end with {RIGHT_SUFFIX}, otherwise there may be unintended consequences.")

        else: raise TypeError(f"Wrong spin: {self.rep.spin}. Should be LH Weyl or RH Weyl.")

    def __repr__(self):
        return f'{self.label}'

    def __eq__(self, other: WeylComponent) -> bool:
        return self.label == other.label and self.rep == other.rep

    def __lt__(self, other: WeylComponent) -> bool:
        return self.label.__lt__(other.label)

    def __le__(self, other: WeylComponent) -> bool:
        return self.label.__le__(other.label)

    def __gt__(self, other: WeylComponent) -> bool:
        return self.label.__gt__(other.label)

    def __ge__(self, other: WeylComponent) -> bool:
        return self.label.__ge__(other.label)

    def __hash__(self):
        return hash((self.label, self.rep))