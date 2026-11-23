"""
.. module:: Field
   :synopsis: Defines the Field class, which represents a gauge eigenstate

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import copy
from typing import Union, Dict, List, Tuple, Optional
from itertools import product
from uvbuilder.tools.reps import DIRAC_SPINOR, SCALAR, LEFT_SPINOR, RIGHT_SPINOR, WEYL_LEFT, WEYL_RIGHT, VECTOR, LorentzRep, U1Rep, SU3Rep, SU2Rep, AllReps, BrokenRep
from uvbuilder.tools.reps import GlobalU1s

class Field(object):
    """
    An instance of this class represents a single gauge-eigenstate field.

    All the quantum numbers are stored in a single BrokenRep object (self.rep),
    with the flat attributes (spin, colordim, isospin, hypercharge, i3, eCharge,
    globalU1s) acting as read-only read-through properties. The field is
    rigid: it is always created with a complete, well-defined BrokenRep and its
    quantum numbers cannot be modified afterwards.

    The self-conjugacy of the field is derived from the rep: a field is its own
    conjugate iff the rep equals its conjugate (rep.conjRep == rep), i.e. all
    charges are real. The conjugateField is defined once at construction and can
    never be changed.
    """
    __slots__ = ["_rep", "_label", "_isSM", "_conjugateField", "_index"]

    QUANTUM_NUMBERS = ["spin", "colordim", "isospin", "hypercharge", "i3", "globalU1s"]


    def __init__(self, rep: BrokenRep, label: str = "Field", isSM: bool = False, index: int = 1):
        """
        Creates a field.

        :param rep: complete BrokenRep holding all the quantum numbers (required)
        :param label: str, e.g. 'eL'
        :param isSM: whether the field belongs to the SM
        :param index: family index
        """
        if not isinstance(rep, BrokenRep):
            raise ValueError(f"rep should be a BrokenRep object and not {type(rep)}")
        self._rep = rep
        self._label = label
        self._isSM = isSM
        self._index = index

        # Define the conjugate field. The field is rigid, so self-conjugacy is
        # known at construction and never changes: a self-conjugate field is
        # its own conjugate, otherwise a distinct conjugate is created once
        # (recursion-safe via __new__).
        if self.selfConjugate:
            self._conjugateField = self
        else:
            conj = Field.__new__(Field)
            conj._rep = rep.conjRep
            conj._label = self.conjugate_label()
            conj._isSM = isSM
            conj._index = index
            conj._conjugateField = self
            self._conjugateField = conj

    @property
    def rep(self) -> BrokenRep:
        return self._rep

    @property
    def spin(self) -> Optional[LorentzRep]:
        return self._rep.spin

    @property
    def colordim(self) -> Optional[SU3Rep]:
        return self._rep.colordim

    @property
    def isospin(self) -> SU2Rep:
        return self._rep.isospin

    @property
    def hypercharge(self) -> Optional[U1Rep]:
        return self._rep.hypercharge

    @property
    def i3(self) -> Optional[U1Rep]:
        """
        Physical isospin projection i3 (e.g. -1/2 for the electron doublet component).
        """
        return self._rep.i3

    @property
    def eCharge(self) -> Optional[U1Rep]:
        return self._rep.eCharge

    @property
    def globalU1s(self) -> GlobalU1s:
        return self._rep.globalU1s

    @property
    def conjugateField(self) -> Union["Field", None]:
        return self if self.selfConjugate else self._conjugateField

    @conjugateField.setter
    def conjugateField(self, value):
        """
        Make sure the conjugateField attribute is protected
        """
        raise AttributeError("The conjugate field of a field is automatically defined when the field is created "
                             "and it should not be changed directly.")

    @property
    def selfConjugate(self) -> bool:
        """
        Whether the field is its own conjugate. Derived from the rep:
        a field is self-conjugate iff the rep equals its conjugate
        (rep.conjRep == rep), i.e. all charges are real.
        """
        return self._rep.conjRep == self._rep

    @property
    def label(self) -> str:
        return self._label

    @label.setter
    def label(self, value: str):
        """
        Labels are allowed to be modified after field creation, but they should be consistent with the conjugate field's label.
        """
        if not isinstance(value, str):
            raise ValueError(f"label should be of type str and not {type(value)}")

        self._label = value
        if not self.selfConjugate:
            self._conjugateField._label = self.conjugate_label()

    @property
    def isSM(self) -> bool:
        return self._isSM

    @isSM.setter
    def isSM(self, value: bool):
        if not isinstance(value, bool):
            raise ValueError(f"isSM should be of type bool and not {type(value)}")
        self._isSM = value
        if not self.selfConjugate:
            self._conjugateField._isSM = value

    @property
    def index(self) -> int:
        return self._index

    @staticmethod
    def _conjugate_label_string(label: str) -> str:
        """
        Conjugates a label string (ignoring self-conjugacy of the field).
        """

        if label.endswith("+"):
            return label[:-1] + "-"
        elif label.endswith("-"):
            return label[:-1] + "+"
        elif label.endswith("bar"):
            return label[:-3]
        elif label.endswith("~") or label.endswith("*"):
            return label[:-1]
        else:
            return label + "bar"
        

    def conjugate_label(self) -> str:
        """
        Returns the label of the conjugate field (a name only: no phase). Phases of
        conjugated components live on the multiplet (Multiplet.phase), from
        tools.reps.su2_conj_phase.
        """
        if self.selfConjugate:
            return self.label
        return self._conjugate_label_string(self.label)

    def __str__(self) -> str:
        return str(self.label)

    def describe(self) -> str:
        """
        Returns a string describing the field and its quantum numbers.
        """

        description = f"{self}"
        qnumbers = []
        for qnumber in Field.QUANTUM_NUMBERS:
            if hasattr(self, qnumber) and getattr(self, qnumber) is not None:
                qnumbers.append(f"{qnumber}={getattr(self, qnumber)}")
        if qnumbers:
            description += " : " + ", ".join(qnumbers)

        return description

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other) -> bool:
        """
        Compare fields according to their representation and then label.
        """

        # Handle other = None cases.
        if not isinstance(other, Field):
            return False
        return self.__cmp__(other) == 0

    def __ne__(self, other) -> bool:
        return not self.__eq__(other)

    def __cmp__(self, other) -> int:
        """
        Compare fields according to their labels and then representation.
        """

        if self.label != other.label:
            return -1 if self.label < other.label else 1
        elif self.rep != other.rep:
            return -1 if self.dim < other.dim else 1

        return 0

    def __lt__(self, other):

        if not isinstance(other, Field):
            return NotImplemented
        return self.__cmp__(other) == -1

    def __gt__(self, other):
        if not isinstance(other, Field):
            return NotImplemented
        return self.__cmp__(other) == 1

    def __len__(self) -> int:
        """
        Defines the "length" of a field based on how many of its quantum numbers have been defined.
        """

        return len(self.hasQNumbers())

    def __hash__(self) -> int:
        return hash((self.label, self.rep))

    @property
    def dim(self) -> Tuple[int]:
        """
        Defines the dimension of a field as tuple with the dimensions of its representations.
        """

        dim = tuple([])
        for qnumber in Field.QUANTUM_NUMBERS:
            if qnumber == 'i3':
                continue
            if hasattr(self, qnumber) and getattr(self, qnumber) is not None:
                value = getattr(self, qnumber)
                if isinstance(value, (U1Rep, SU2Rep, SU3Rep, LorentzRep)):
                    dim += (value.dim,)
        return dim

    def hasQNumbers(self) -> List[str]:
        """
        Returns a list of the qnumbers which have been defined for the particle
        """

        good_qnumbers = []
        for qnumber in Field.QUANTUM_NUMBERS:
            if hasattr(self, qnumber) and getattr(self, qnumber) is not None:
                good_qnumbers.append(qnumber)
        return good_qnumbers

    @property
    def representation(self) -> Dict[str, AllReps]:
        """
        Returns a dictionary with the field representation (quantum numbers)
        """

        repDict = {}
        for qnumber in Field.QUANTUM_NUMBERS:
            repDict[qnumber] = getattr(self, qnumber)

        return repDict

    def copy(self):
        """
        Create a copy of the field (with a fresh conjugate pair).
        """

        return Field(rep=self._rep, label=self._label, isSM=self._isSM)

    def representationProduct(self, other: "Field") -> List[Dict[str, AllReps]]:
        """
        Returns the representation of the product of self and other.
        """

        repDict = {}
        for qnumber in Field.QUANTUM_NUMBERS:
            value1 = getattr(self, qnumber)
            value2 = getattr(other, qnumber)
            if value1 is not None and value2 is not None:
                if isinstance(value1, (SU2Rep, SU3Rep, LorentzRep)):
                    value_prod = value1 * value2
                    # Make sure the product is always a list (even if it contains only one representation)
                    if isinstance(value_prod, AllReps):
                        value_prod = [value_prod]
                    repDict[qnumber] = list(value_prod)
                elif qnumber == 'i3' or isinstance(value1, U1Rep):
                    # Handle the particular case where the product of quantum numbers is just the
                    # sum of the charges (e.g. for U(1) representations)
                    repDict[qnumber] = [value1 + value2]
                elif isinstance(value1, GlobalU1s) and isinstance(value2, GlobalU1s):
                    # This is the case when dealing with globalU1s
                    repDict[qnumber] = [value1 + value2]
                else:
                    raise TypeError(f"Cannot compute representation product for qnumber {qnumber} with values {value1} and {value2}")
            elif value1 is not None:
                repDict[qnumber] = [value1]
            elif value2 is not None:
                repDict[qnumber] = [value2]

        if not repDict:
            return [repDict]

        qnumbers = list(repDict.keys())
        repProducts = product(*(repDict[qnumber] for qnumber in qnumbers))

        return [dict(zip(qnumbers, repProduct)) for repProduct in repProducts]

    def __deepcopy__(self, memo):
        """
        Natively supports copy.deepcopy.
        Bypasses __init__ to ensure no dummy fields are spawned and
        preserves the exact internal state.
        """
        if id(self) in memo:
            return memo[id(self)]

        # Instantiate a completely blank Field (bypasses __init__ side-effects)
        new_field = Field.__new__(Field)

        # Register in memo to safely handle Majorana cyclic references
        # (where conjugateField points to self)
        memo[id(self)] = new_field

        # Copy the internal state exactly as it is, nothing stripped.
        for attr in self.__slots__:
            if hasattr(self, attr):
                val = getattr(self, attr)
                setattr(new_field, attr, copy.deepcopy(val, memo))

        return new_field
