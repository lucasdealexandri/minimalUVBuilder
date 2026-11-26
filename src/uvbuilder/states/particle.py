"""
.. module:: particle
   :synopsis: Defines the Particle class, which represents a mass eigenstate

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from __future__ import annotations

from typing import Union, Dict, List, Optional, Self
from collections import namedtuple
from dataclasses import replace

from unum import Unum

from uvbuilder.tools.reps import (DIRAC_SPINOR, LEFT_SPINOR, RIGHT_SPINOR, WEYL_LEFT,
                        WEYL_RIGHT, LorentzRep, U1Rep, SU3Rep, UnbrokenRep,
                        DiracSpinor, WeylSpinor, Scalar, AllReps)
from uvbuilder.tools.reps import GlobalU1s
from uvbuilder.tools.components import ScalarComponent, WeylComponent, _fr_charge
from uvbuilder.tools.fr_utils import (formatClassBlock, LEFT_SUFFIX, RIGHT_SUFFIX,
                            chirality_suffix)
from uvbuilder.states.field import Field
from uvbuilder.states.physicsUnits import GeV

fTerm = namedtuple('fTerm',['coeff','field'])


class Particle(object):
    """
    An instance of this class represents a single particle (mass eigenstate).

    The properties are: label, pdg, mass, spin, electric charge, colordim, width.
    All the quantum numbers are stored in a single UnbrokenRep object (self.rep),
    with the flat attributes (spin, colordim, eCharge, globalU1s) acting as
    read-through properties. The particle can be defined progressively: a bare
    Particle (with an empty UnbrokenRep) can be created and its quantum numbers
    filled in one at a time.
    """
    __slots__ = ["_rep", "_pdg", "_label", "_isSM", "_antiParticle",
                 "_mass", "_totalwidth", "_fieldExpansion"]

    # This distinguishes particles necessary only to fill
    # out multiplets from those present in the decays.
    QUANTUM_NUMBERS = list(UnbrokenRep.__dataclass_fields__)

    def __init__(self, pdg: int, rep: Optional[UnbrokenRep] = None, label: Optional[str] = None,
                 mass: Optional[Unum] = None, totalwidth: Optional[Unum] = None,
                 isSM: bool = False, fields: Optional[List[Field]] = None):
        """
        Creates a particle.

        :param pdg: number in pdg (required)
        :param rep: UnbrokenRep holding all the quantum numbers (defaults to an empty rep)
        :param label: str, e.g. 'e-'
        :param mass: mass of the particle
        :param totalwidth: total width
        :param isSM: whether the particle belongs to the SM
        :param fields: list of gauge-eigenstate fields making up the particle

        The self-conjugacy of the particle is derived from the rep: a particle
        is its own anti-particle iff the rep equals its conjugate
        (rep.conjRep == rep), i.e. all charges (including the additional UQNs)
        are real.
        """

        if pdg is None:
            raise ValueError("Particles/anti-particles can only be created with a PDG number")
        if rep is None:
            rep = UnbrokenRep()
        if not isinstance(rep, UnbrokenRep):
            raise ValueError(f"rep should be an UnbrokenRep object and not {type(rep)}")

        self._rep = rep
        self._pdg = pdg
        self._label = label if label is not None else str(pdg)
        self._isSM = isSM
        self._mass = mass
        self._totalwidth = totalwidth
        self._fieldExpansion = None

        # Define the antiParticle. It is always created (even for a currently
        # self-conjugate particle) so that progressive filling of the quantum
        # numbers - which may turn a real rep into a complex one - can later
        # expose it through the antiParticle property.
        anti = Particle.__new__(Particle)
        anti._rep = rep.conjRep
        anti._pdg = -pdg
        anti._label = self._conjugate_label_string(self._label, pdg)
        anti._isSM = isSM
        anti._mass = mass
        anti._totalwidth = totalwidth
        anti._fieldExpansion = None
        anti._antiParticle = self
        self._antiParticle = anti

        if fields is not None:
            self.fields = fields

    def _set_qnumber(self, qnumber: str, value) -> None:
        """
        Set the quantum number `qnumber` of the particle's rep to `value`.
        The conjugate value is mirrored to the anti-particle. Raises an
        AttributeError if the quantum number has already been set to a
        different value (used to discard inconsistent possibilities).
        """
        if value is None:
            return
        current = getattr(self._rep, qnumber)
        if current is not None and current != value:
            raise AttributeError(f"{qnumber} for {self} has already been set. Trying to set it to value {value}")
        elif current is not None:
            return
        new_rep = replace(self._rep, **{qnumber: value})
        self._rep = new_rep
        if not self.selfConjugate and self._antiParticle is not None:
            self._antiParticle._rep = new_rep.conjRep

    @property
    def rep(self) -> UnbrokenRep:
        return self._rep

    @rep.setter
    def rep(self, value: UnbrokenRep):
        if not isinstance(value, UnbrokenRep):
            raise ValueError(f"rep should be an UnbrokenRep object and not {type(value)}")
        self._rep = value
        if not self.selfConjugate and self._antiParticle is not None:
            self._antiParticle._rep = value.conjRep

    @property
    def spin(self) -> Optional[LorentzRep]:
        return self._rep.spin

    @spin.setter
    def spin(self, value: LorentzRep):
        if not isinstance(value, LorentzRep):
            raise ValueError(f"Spin should be of type LorentzRep and not {type(value)} for {self}")
        self._set_qnumber('spin', value)

    @property
    def colordim(self) -> Optional[SU3Rep]:
        return self._rep.colordim

    @colordim.setter
    def colordim(self, value: SU3Rep):
        if not isinstance(value, SU3Rep):
            raise ValueError(f"colordim should be of type SU3Rep and not {type(value)} for {self}")
        self._set_qnumber('colordim', value)

    @property
    def eCharge(self) -> Optional[U1Rep]:
        return self._rep.eCharge

    @eCharge.setter
    def eCharge(self, value: U1Rep):
        if not isinstance(value, U1Rep):
            raise ValueError(f"eCharge should be of type U1Rep and not {type(value)} for {self}")
        self._set_qnumber('eCharge', value)

    @property
    def globalU1s(self) -> GlobalU1s:
        return self._rep.globalU1s

    @globalU1s.setter
    def globalU1s(self, value: GlobalU1s):
        if not isinstance(value, GlobalU1s):
            raise ValueError(f"globalU1s should be an GlobalU1s object and not {type(value)}")
        self._set_qnumber('globalU1s', value)

    @property
    def antiParticle(self) -> Union[None, "Particle"]:
        """
        Returns the anti-particle. For a self-conjugate particle this is the
        particle itself; otherwise it is the automatically-created conjugate.
        """
        if self.selfConjugate:
            return self
        return self._antiParticle

    @antiParticle.setter
    def antiParticle(self, value):
        """
        Make sure the antiParticle attribute is protected
        """
        raise AttributeError("The anti-particle of a particle is automatically defined when the particle is created "
                             "and it should not be changed directly.")

    @property
    def label(self) -> Union[str, None]:
        return self._label

    @label.setter
    def label(self, value: str):
        """
        Labels are allowed to be modified after particle creation, but they should be consistent with the anti-particle's label.
        """
        if not isinstance(value, str):
            raise ValueError(f"label should be of type str and not {type(value)}")

        self._label = value
        if not self.selfConjugate and self.antiParticle is not None:
            self.antiParticle._label = self.conjugate_label()

    @property
    def mass(self) -> Union[Unum, None]:
        return self._mass

    @mass.setter
    def mass(self, value: Unum):
        if not isinstance(value, Unum):
            raise ValueError(f"mass should be of type Unum and not {type(value)}")
        elif self.mass is not None and self.mass != value:
            raise AttributeError(f"Mass for {self} has already been set. Trying to set it to value {value}")
        elif self.mass is None:
            self._mass = value
            if not self.selfConjugate and self.antiParticle is not None:
                self.antiParticle._mass = value

    @property
    def totalwidth(self) -> Union[Unum, None]:
        return self._totalwidth

    @totalwidth.setter
    def totalwidth(self, value: Unum):
        if not isinstance(value, Unum):
            raise ValueError(f"totalwidth should be of type Unum and not {type(value)}")
        elif self.totalwidth is not None and self.totalwidth != value:
            raise AttributeError(f"Width for {self} has already been set. Trying to set it to value {value}")
        elif self.totalwidth is None:
            self._totalwidth = value
            if not self.selfConjugate and self.antiParticle is not None:
                self.antiParticle._totalwidth = value

    @property
    def pdg(self) -> Union[None, int]:
        return self._pdg

    @property
    def pdgList(self) -> List[int]:
        return [self._pdg]

    @property
    def isSM(self) -> Union[bool, None]:
        return self._isSM

    @isSM.setter
    def isSM(self, value: bool):
        if not isinstance(value, bool):
            raise ValueError(f"isSM should be of type bool and not {type(value)}")
        elif self.isSM is not None and self.isSM != value:
            raise AttributeError(f"isSM for {self} has already been set ({self.isSM}). Trying to set it to value {value}")
        elif self.isSM is None:
            self._isSM = value
            if not self.selfConjugate and self.antiParticle is not None:
                self.antiParticle._isSM = value

    @property
    def selfConjugate(self) -> bool:
        """
        Whether the particle is its own anti-particle. Derived from the rep:
        a particle is self-conjugate iff the rep equals its conjugate
        (rep.conjRep == rep), i.e. all charges (including the additional UQNs)
        are real. While a charge (eCharge, colordim or globalU1s) is unset, the
        particle is not self-conjugate: its anti-particle stays a distinct state
        until the charges are known.
        """
        if any(getattr(self._rep, q) is None for q in ('eCharge', 'colordim', 'globalU1s')):
            return False
        return self._rep.conjRep == self._rep

    @property
    def components(self) -> List[Union[WeylComponent, ScalarComponent]]:
        """
        Returns the Weyl/scalar components of the particle based on its spin:
          - Dirac/Weyl spinors: one WeylComponent per non-empty chirality slot,
            with the same gauge charges as the particle. Labels follow the
            convention ``<particle label>Lw`` / ``<particle label>Rw`` for the
            left- and right-handed Weyls (e.g. e- -> [e-Lw, e-Rw]).
          - Scalars: a complex scalar gives (phi, phi^dagger), a real (self-conjugate)
            scalar gives a single component.
          - Vectors: not implemented yet, returns an empty list.
        Fresh objects are generated on each access.
        """
        rep = self._rep
        if self.selfConjugate:
            assert rep.conjRep == rep, (f"Self-conjugate particle {self} has a non-real representation {rep}.")

        spin = rep.spin
        if isinstance(spin, DiracSpinor):
            comps = []
            if len(spin.left):
                comps.append(WeylComponent(f'{self.label}{LEFT_SUFFIX}', replace(rep, spin=WEYL_LEFT)))
            if len(spin.right):
                comps.append(WeylComponent(f'{self.label}{RIGHT_SUFFIX}', replace(rep, spin=WEYL_RIGHT)))
            return comps
        elif isinstance(spin, WeylSpinor):
            return [WeylComponent(f'{self.label}{chirality_suffix(rep.spin)}', rep)]
        elif isinstance(spin, Scalar):
            component = ScalarComponent(self.label, rep)
            if self.selfConjugate:
                return [component]
            return [component, component.dagger]
        return []

    @property
    def fieldExpansion(self) -> Union[None, List[fTerm]]:
        return self._fieldExpansion

    @fieldExpansion.setter
    def fieldExpansion(self, value: List[fTerm]):
        if not isinstance(value, list) or not all(isinstance(ft, fTerm) for ft in value):
            raise ValueError(f"fieldExpansion should be a list of fieldTerms and not {type(value)}")
        self.setFieldExpansion(value)

    @property
    def fields(self) -> List[Field]:
        """
        Get the list of fields in the particle's field expansion.
        """
        if self.fieldExpansion is None:
            return []
        return [fTerm.field for fTerm in self.fieldExpansion]

    @fields.setter
    def fields(self, value: List[Field]):
        """
        Convenience method to allow for setting the field expansion from a list of fields (instead of a list of fieldTerms).
        The coefficient for each field will be set to 1.0.
        """
        if not isinstance(value, list) or not all(isinstance(f, Field) for f in value):
            raise ValueError(f"fields should be a list of Field instances and not {type(value)}")

        fTerms = [fTerm(coeff=1.0, field=f) for f in value]
        self.setFieldExpansion(fTerms)

    @property
    def representation(self) -> Dict[str, AllReps]:
        """
        Returns a dictionary with the particle representation (quantum numbers)
        """

        repDict = {}
        for qnumber in Particle.QUANTUM_NUMBERS:
            repDict[qnumber] = getattr(self, qnumber)

        return repDict

    @staticmethod
    def _conjugate_label_string(label: str, pdg: int) -> str:
        """
        Conjugates a label string (ignoring self-conjugacy of the particle).
        """

        if label == str(pdg):
            return str(-pdg)
        elif label.endswith("+"):
            return label[:-1] + "-"
        elif label.endswith("-"):
            return label[:-1] + "+"
        elif label.endswith("bar"):
            return label[:-3]
        elif label.endswith("~"):
            return label[:-1]
        else:
            return label + "~"

    def conjugate_label(self) -> str:
        """
        Returns the label of the conjugate particle.
        """

        if self.selfConjugate:
            return self.label
        return self._conjugate_label_string(self.label, self.pdg)

    def __str__(self) -> str:
        if hasattr(self, "label"):
            return str(self.label)
        elif hasattr(self, "pdg"):
            return str(self.pdg)
        else:
            return ""

    def __repr__(self) -> str:
        return self.__str__()

    def describe(self) -> str:
        """
        Returns a string describing the particle and its quantum numbers.
        """

        description = f"{self} (PDG: {self.pdg})"
        qnumbers = []
        for qnumber in Particle.QUANTUM_NUMBERS:
            if getattr(self, qnumber) is not None:
                qnumbers.append(f"{qnumber}={getattr(self, qnumber)}")
        if qnumbers:
            description += " : " + ", ".join(qnumbers)

        return description

    def to_feynrules(self, tracker, family=None) -> Dict[str, List[str]]:
        """FeynRules class declarations for this particle, as an FRChunk.

        For Dirac particles, emits the two ``W[...]`` Weyl component classes
        (via ``WeylComponent.to_feynrules``) and the ``F[...]`` Dirac class
        with ``WeylComponents``. For Weyl particles, emits the single
        ``W[...]`` class. Scalars emit no component declaration.

        ``family`` (an optional FRIndex) adds a generation index to the
        component classes and the F class. The index ``tracker`` provides
        class numbers and gauge indices.
        """

        from uvbuilder.tools.fr_utils import merge_chunks

        spin = self.rep.spin
        if isinstance(spin, Scalar):
            return {'Classes': []}

        comps = self.components
        if self.selfConjugate: comps = [self.components[0]]
        chunk = merge_chunks(*[
            comp.to_feynrules(tracker, family) for comp in comps])
        if isinstance(spin, WeylSpinor):
            return chunk

        # Dirac: the F class on top of the two W component classes
        header = f"F[{tracker.class_number('F')}]"
        weyl_members = ', '.join(comp.label for comp in comps)
        blockDict = {
            'ClassName': self.label,
            'WeylComponents': '{' + weyl_members + '}',
            'SelfConjugate': 'True' if self.selfConjugate else 'False',
            'QuantumNumbers': f"{{Q -> {_fr_charge(self.rep.eCharge)}}}",
        }
        indices = [str(index) for index in tracker.gauge(self.rep)]
        if family is not None:
            indices.append(str(family))
        if indices:
            blockDict['Indices'] = f"{{{','.join(indices)}}}"
        if family is not None:
            blockDict['FlavorIndex'] = family.name
        if self.mass is not None:
            blockDict['Mass'] = f'{{M{self.label}, {self.mass.asNumber(GeV)}}}'
        if self.totalwidth is not None:
            blockDict['Width'] = f'{{W{self.label}, {self.totalwidth.asNumber(GeV)}}}'
        blockDict['PDG'] = f'{{{self.pdg}}}'

        chunk.setdefault('Classes', []).append(
            formatClassBlock(header, blockDict))
        return chunk

    def __eq__(self, other) -> bool:
        """
        Compare particles according to their PDG and quantum numbers (the rep
        fully encodes self-conjugacy through rep.conjRep == rep).
        """

        if type(self) != type(other):
            return False
        elif not hasattr(other, 'pdg'):
            return False
        elif self.pdg != other.pdg:
            return False
        elif self.representation != other.representation:
            return False

        return True

    def __ne__(self, other) -> bool:
        return not (self == other)

    def __lt__(self, other) -> bool:
        if isinstance(other.pdg, list):
            other_pdg = max(other.pdg)
        else:
            other_pdg = other.pdg

        return other_pdg > self.pdg

    def __gt__(self, other) -> bool:
        if isinstance(other.pdg, list):
            other_pdg = max(other.pdg)
        else:
            other_pdg = other.pdg

        return self.pdg > other_pdg

    def __hash__(self) -> int:
        return id(self)

    def __len__(self) -> int:
        """
        Defines the "length" of a particle based on how many of its quantum numbers have been defined.
        """

        return len(self.hasQNumbers())

    def hasQNumbers(self) -> List[str]:
        """
        Returns a list of the qnumbers which have been defined for the particle
        """

        good_qnumbers = []
        for qnumber in Particle.QUANTUM_NUMBERS:
            if getattr(self, qnumber) is not None:
                good_qnumbers.append(qnumber)
        return good_qnumbers

    def copy(self) -> "Particle":
        """
        Create a copy of the particle.
        """

        fieldExpansion = self._fieldExpansion
        if fieldExpansion is None:
            fieldExpansion_copy = None
        else:
            # A copy of the field expansion should be used to guarantee that the new particle
            # points to distinct field objects.
            fieldExpansion_copy = [fTerm(coeff=ft.coeff, field=ft.field.copy())
                                   for ft in fieldExpansion]
        self_copy = Particle(pdg=self._pdg, rep=self._rep, label=self._label,
                             mass=self._mass, totalwidth=self._totalwidth,
                             isSM=self._isSM)
        self_copy._fieldExpansion = fieldExpansion_copy

        return self_copy

    def setFieldExpansion(self, fieldTerms: List[fTerm] | List[Field]):
        """
        Set the particle's field expansion to the specified list of field terms.
        If a list of fields is provided, the coefficient for each field will be set to 1.0.
        Checks if the fields match the particle quantum numbers.
        """

        if self.fieldExpansion is None:
            self._fieldExpansion = []
        self._addToFieldExpansion(fieldTerms)

    def _addToFieldExpansion(self, fieldTerm: Union[fTerm, Field, List[fTerm], List[Field]],
                             addToAntiParticle=True):
        """
        Add a field term to the particle's field expansion. If the representations
        of the field(s) to be added already appear in the expansion, add to the corresponding expansion coefficient(s).
        If the representations do not appear in the expansion, append the new field term(s) to the expansion.
        If addToAntiParticle is True, the same field term will be added to the anti-particle's field expansion,
        but with the conjugate field and the same coefficient.
        """

        if isinstance(fieldTerm, List):
            for f in fieldTerm:
                if isinstance(f, Field):
                    f = fTerm(coeff=1.0, field=f)
                self._addToFieldExpansion(f, addToAntiParticle=addToAntiParticle)

        elif isinstance(fieldTerm, Field):
            fieldTerm = fTerm(coeff=1.0, field=fieldTerm)
            self._addToFieldExpansion(fieldTerm, addToAntiParticle=addToAntiParticle)
        elif not isinstance(fieldTerm, fTerm):
            raise TypeError(f"Expected a fieldTerm instance, got {type(fTerm)}")
        else:
            if self._fieldExpansion is None:
                self._fieldExpansion = []
            if not self.qnumbersMatch(fieldTerm.field):
                raise ValueError(f"Quantum numbers of the field {fieldTerm.field.label} do not match those of the particle {self.label}.")
            f_index = self.containsFieldRepAt(fieldTerm.field)
            # If field is not found, append it to list, otherwise add the coefficients
            if f_index is None:
                self.fieldExpansion.append(fieldTerm)
            else:
                self.fieldExpansion[f_index] = fTerm(coeff=self.fieldExpansion[f_index].coeff + fieldTerm.coeff,
                                                    field=self.fieldExpansion[f_index].field)
            if addToAntiParticle and not self.selfConjugate and self.antiParticle is not None:
                conjugate_field = fieldTerm.field.conjugateField
                self.antiParticle._addToFieldExpansion(fTerm(coeff=fieldTerm.coeff,
                                                            field=conjugate_field),
                                                            addToAntiParticle=False)

    def qnumbersMatch(self, field: Field | Self) -> bool:
        """
        Check if the quantum numbers of the particle match those of the field/particle.
        """
        for qnumber in Particle.QUANTUM_NUMBERS:
            if hasattr(field, qnumber) and hasattr(self, qnumber):
                if getattr(field, qnumber) is None or getattr(self, qnumber) is None:
                    continue
                elif qnumber == 'spin':
                    # For spin, we need to check if the the spin representation is contained
                    # in the particle's spin representation (since the particle can be a sum of left/right spinor).
                    if self.spin == DIRAC_SPINOR:
                        if field.spin not in [LEFT_SPINOR, RIGHT_SPINOR, DIRAC_SPINOR]:
                            return False
                    elif self.spin != field.spin:
                        return False
                elif getattr(field, qnumber) != getattr(self, qnumber):
                    return False
        return True

    def fieldExpansionContains(self, field: Field | List[Field]) -> bool:
        """
        Check if the particle's field expansion contains the representation of the specified field.
        """

        if isinstance(field, List):
            return any(self.fieldExpansionContains(f) for f in field)
        else:
            return (self.containsFieldRepAt(field) is not None)

    def containsFieldRepAt(self, field: Field) -> Union[int, None]:
        """
        Returns the field position in the particle field expansion (using the field representation).
        If the field representation is not found, returns None.
        """
        if self.fieldExpansion is None:
            return None
        for i, f in enumerate(self.fieldExpansion):
            if f.field.representation == field.representation:
                return i
        return None

def family_key(particle) -> tuple:
    '''Order of the members of a physical class (its family index): by mass, then PDG.'''
    if particle.mass is None:
        raise ValueError(f'{particle} has no mass: the family index is its mass rank')
    return (particle.mass.asNumber(GeV), particle.pdg)
