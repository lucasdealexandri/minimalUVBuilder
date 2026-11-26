"""
.. module:: Multiplet
   :synopsis: Defines the Multiplet class, which represents a SU(2)_L multiplet

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import copy
from typing import Union, Dict, Sequence, List, Optional, Tuple
from uvbuilder.tools.reps import (LorentzRep, U1Rep, SU3Rep, SU2Rep, BrokenRep,
                        WeylSpinor, Scalar, Vector, WEYL_LEFT)
from uvbuilder.tools.fr_utils import (formatVariable, getGaugeDefinitions, getSU2Index,
                            DAGGER_SUFFIX, chirality_suffix, toggle_dagger)
from uvbuilder.states.field import Field
from uvbuilder.tools.reps import GlobalU1s, su2_conj_phase, su2_reality_phase
from uvbuilder.partition.particle_mixing import ParticleMixing


class Multiplet(object):
    """
    An instance of this class represents a single SU(2)_L gauge multiplet.

    All the quantum numbers are stored in a single BrokenRep object describing
    the whole multiplet (i3=None), with the flat attributes (spin, colordim,
    isospin, hypercharge, globalU1s) acting as read-only read-through
    properties. The multiplet is rigid: it is always created with a complete,
    well-defined BrokenRep and neither its quantum numbers nor its components
    can be modified afterwards.

    The self-conjugacy of the multiplet is derived from its rep: a multiplet
    is self-conjugate iff the internal rep equals its conjugate
    (rep == rep.conjInternalRep), i.e. it transforms under a real
    representation of the gauge group (spin is not conjugated). This is a
    property of the multiplet as a whole: individual components need not be
    self-conjugate (e.g. the W triplet is self-conjugate, but its W+/W-
    charge eigenstates swap under conjugation). The conjugateMultiplet is
    defined once at construction and can never be changed.
    """
    __slots__ = ["_rep", "_components", "_label", "_labels", "_isSM",
                 "_conjugateMultiplet", "index", "_phases"]

    _EXO_COUNTER = 0
    _EXOTIC_LABEL = "M"

    def __init__(
            self,
            rep: BrokenRep,
            label: Optional[str] = None,
            fields: Sequence[Field] = (),
            labels: Optional[Dict[U1Rep, str]] = None,
            fillMultiplet: bool = False,
            isSM: bool = False,
            index: int = 1
    ):
        """
        Creates a multiplet.

        :param rep: complete whole-multiplet BrokenRep (i3=None) holding all the quantum numbers (required)
        :param label: str to name the whole multiplet
        :param fields: component fields to be placed in the multiplet
        :param labels: dict mapping i3 to a component label (used for the fields created by fillMultiplet)
        :param fillMultiplet: whether to fill the empty components with new, exotic fields
        :param isSM: whether the multiplet belongs to the SM
        :param index: family index
        """
        if not isinstance(rep, BrokenRep):
            raise ValueError(f"rep should be a BrokenRep object and not {type(rep)}")
        if rep.i3 is not None:
            raise ValueError("Multiplet requires a whole-multiplet BrokenRep (i3=None)")

        self._rep = rep
        self._isSM = isSM
        self.index = index
        self._label = label
        self._labels: Dict[U1Rep, str] = {U1Rep(i3): label for i3, label in (labels or {}).items()}
        self._components: Dict[U1Rep, Union[Field, None]] = {
            i3: None for i3 in rep.isospin.i3_range(reverse=True)
        }
        # Phase with which each component field enters the multiplet: 1 here; the
        # conjugate multiplet carries the SU(2) conjugation phases (see phase()).
        self._phases: Dict[U1Rep, int] = {i3: 1 for i3 in self._components}

        for field in fields:
            self._add_field(field)

        if fillMultiplet:
            self._fill()

        # Define the conjugate multiplet. The multiplet is rigid, so
        # self-conjugacy is known at construction and never changes: a
        # self-conjugate multiplet is its own conjugate, otherwise a distinct
        # conjugate is created once (recursion-safe via __new__).
        if self.selfConjugate:
            self._conjugateMultiplet = self
        else:
            conj = Multiplet.__new__(Multiplet)
            conj._rep = rep.conjRep
            conj._isSM = isSM
            conj.index = index
            conj._label = None if label is None else toggle_dagger(label)
            conj._labels = {U1Rep(-i3): field.conjugateField.label
                            for i3, field in self._components.items() if field is not None}
            conj._components = {U1Rep(-i3): (field.conjugateField if field is not None else None)
                                for i3, field in self._components.items()}
            # (psi^c)_{-m} = (-1)^(j+m) (psi_m)^*  (tools.reps.su2_conj_phase)
            conj._phases = {U1Rep(-i3): su2_conj_phase(rep.isospin, i3)
                            for i3 in self._components}
            conj._conjugateMultiplet = self
            self._conjugateMultiplet = conj

    def __repr__(self) -> str:
        return f"Multiplet({self._components})"

    @property
    def rep(self) -> BrokenRep:
        return self._rep

    @property
    def label(self) -> str:
        return self._label

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
    def globalU1s(self) -> GlobalU1s:
        return self._rep.globalU1s

    @property
    def selfConjugate(self) -> bool:
        """
        Whether the multiplet is its own conjugate. Derived from the rep:
        a multiplet is self-conjugate iff rep == rep.conjInternalRep, i.e.
        it transforms under a real representation of the internal gauge group
        (spin is not conjugated). This is a property of the multiplet as a
        whole: individual components need not be self-conjugate.
        """
        return self._rep == self._rep.conjInternalRep

    @property
    def isSM(self) -> bool:
        return self._isSM

    @isSM.setter
    def isSM(self, value: bool):
        if not isinstance(value, bool):
            raise ValueError(f"isSM should be of type bool and not {type(value)}")
        self._isSM = value
        if not self.selfConjugate:
            self._conjugateMultiplet._isSM = value

    @property
    def is_real_nontrivial(self) -> bool:
        """
        Self-conjugate scalar with isospin dim > 1 (e.g. a Y = 0 triplet). Its
        reality condition in the charge basis is psi_{-m} = (-1)^m (psi_m)^*, so it
        is declared as a complex class with that relation in its definitions.
        """
        return self.selfConjugate and isinstance(self.rep.spin, Scalar) and self.isospin.dim > 1

    @property
    def is_color_conjugated(self) -> bool:
        '''LH Weyl multiplet in 3bar or 6bar. FeynRules only has generators for 3 and 6,
        so it is declared through its RH conjugate (see to_feynrules).'''
        return self.rep.spin == WEYL_LEFT and self.colordim.lower > self.colordim.upper

    def fr_symbol(self) -> str:
        '''FeynRules/Mathematica symbol for this multiplet.

        Follows the component naming convention: Weyl multiplets carry a
        ``Lw``/``Rw`` chirality suffix (e.g. ``SLw``, ``DuLw``), scalars and
        vectors keep the bare label.
        '''
        label = self.label
        if isinstance(self.rep.spin, WeylSpinor):
            return f'{label}{chirality_suffix(self.rep.spin)}'
        return label

    def generate_fr_kinetic(self, tracker, symbol: Optional[str] = None) -> str:
        '''Raw FeynRules kinetic term for this multiplet.

        The caller (FeynRulesModel) wraps the returned term inside
        ``Block``/``ExpandIndices``. The index ``tracker`` is duck-typed: it
        must provide ``index(name, range, index_type)`` (canonical index
        object), ``gauge(rep)`` (list of gauge FRIndex objects) and
        ``new_symbol(style)`` (unique contraction label). Weyl fermions get
        two unique spin indices (``sp1``, ``sp2``) and the ``si``/``sibar``
        projector; scalars get the ``DC`` form. Vectors/Dirac return '' for
        now. Contraction labels are reset per term (one call = one term).
        '''
        tracker.reset_symbols()

        spin = self.rep.spin
        symbol = symbol or self.fr_symbol()
        bar = symbol + DAGGER_SUFFIX   # FeynRules conjugate: always append
        factor = '1/2 ' if self.selfConjugate else ''
        # Real multiplet in the charge basis: 1/2 DC[Xbar_i].DC[X_i], not 1/2 DC[X_i].DC[X_i]
        use_bar = not self.selfConjugate or self.is_real_nontrivial

        gauge = [tracker.new_symbol(index.style) for index in tracker.gauge(self.rep)]
        gauge_str = ', '.join(gauge) if gauge else ''

        if isinstance(spin, WeylSpinor):
            mu = tracker.index('mu', 4, Vector).name
            sp1 = tracker.new_symbol('sp')
            sp2 = tracker.new_symbol('sp')
            lh = self.rep.spin == WEYL_LEFT

            # i chibar sigmabar^mu D_mu chi, with chi = X (LH) or chi = Xbar (RH X)
            if len(gauge) == 1: gauge_str = ', ' + gauge_str
            term = f'I sibar[{mu}, {sp1}, {sp2}] '
            term += f'({bar if lh else symbol}[{sp1}{gauge_str}].DC[{symbol if lh else bar}[{sp2}{gauge_str}], {mu}])'

            return f' + {term}'

        if isinstance(spin, Scalar):
            mu = tracker.index('mu', 4, Vector).name
            if gauge_str: gauge_str = f'[{gauge_str}]'
            return f'+ {factor}DC[{bar if use_bar else symbol}{gauge_str}, {mu}].DC[{symbol}{gauge_str}, {mu}]'

        return ''

    def to_feynrules(self, tracker, mixings: ParticleMixing,
                     rh_symbol: Optional[str] = None) -> Dict[str, List[str]]:
        '''FeynRules chunk for this multiplet: the W/S class block and its
        kinetic term. The index ``tracker`` is duck-typed (class numbers,
        gauge indices). Color-conjugated multiplets need ``rh_symbol``, the
        name of their RH conjugate class.'''

        from uvbuilder.tools.fr_utils import formatClassBlock, fr_quantum_numbers

        if self.is_color_conjugated:
            return self._color_conjugated_feynrules(tracker, mixings, rh_symbol)

        spin_label = formatVariable(self.rep.spin)
        header = f"{spin_label}[{tracker.class_number(spin_label)}]"
        self_conjugate = (self.selfConjugate and not isinstance(self.rep.spin, WeylSpinor)
                          and not self.is_real_nontrivial)
        blockDict = {
            'ClassName': self.fr_symbol(),
            'Unphysical': 'True',
            'SelfConjugate': 'True' if self_conjugate else 'False',
        }
        # Self-conjugate classes carry no quantum numbers (as in ParticleModel)
        if not self_conjugate:
            blockDict['QuantumNumbers'] = fr_quantum_numbers({'Y': self.rep.hypercharge},
                                                             self.rep.globalU1s)
        if isinstance(self.rep.spin, WeylSpinor):
            blockDict['Chirality'] = ('Left' if self.rep.spin == WEYL_LEFT
                                      else 'Right')
        indices = [str(index) for index in tracker.gauge(self.rep)]
        if indices:
            blockDict['Indices'] = f"{{{','.join(indices)}}}"
        if self.rep.isospin.dim > 1:
            blockDict['FlavorIndex'] = tracker.gauge(self.rep)[-1].name

        blockDict['Definitions'] = getGaugeDefinitions(self, mixings,
                                                       phases=getattr(tracker, 'mass_phases', None))

        kinetic = self.generate_fr_kinetic(tracker)
        return {
            'Classes': [formatClassBlock(header, blockDict)],
            'KineticTerms': [kinetic] if kinetic else [],
        }

    def _color_conjugated_feynrules(self, tracker, mixings, rh_symbol) -> Dict[str, List[str]]:
        '''Two classes: the RH conjugate rh_symbol = X^c in 3/6, which carries the
        kinetic term and the physical definitions, and X itself as an unphysical
        class defined through rh_symbol-bar, used by the operators.'''
        from uvbuilder.tools.fr_utils import formatClassBlock, fr_quantum_numbers
        conj = self.conjugateMultiplet
        gauge = tracker.gauge(conj.rep)
        common = {'Unphysical': 'True', 'SelfConjugate': 'False',
                  'Indices': '{' + ','.join(str(i) for i in gauge) + '}'}
        if self.isospin.dim > 1:
            common['FlavorIndex'] = gauge[-1].name
        rh = {'ClassName': rh_symbol, 'Chirality': 'Right', **common,
              'QuantumNumbers': fr_quantum_numbers({'Y': conj.rep.hypercharge}, conj.rep.globalU1s),
              'Definitions': getGaugeDefinitions(self, mixings, symbol=rh_symbol, conjugate=True,
                                                 phases=getattr(tracker, 'mass_phases', None))}
        lh = {'ClassName': self.fr_symbol(), 'Chirality': 'Left', **common,
              'QuantumNumbers': fr_quantum_numbers({'Y': self.rep.hypercharge}, self.rep.globalU1s),
              'Definitions': self.conjugate_helper_definitions(rh_symbol)}
        return {'Classes': [formatClassBlock(f"W[{tracker.class_number('W')}]", rh),
                            formatClassBlock(f"W[{tracker.class_number('W')}]", lh)],
                'KineticTerms': [conj.generate_fr_kinetic(tracker, symbol=rh_symbol)]}

    def conjugate_helper_definitions(self, rh_symbol: str) -> str:
        '''X_m = (-1)^(j+m) (X^c_{-m})^*, with X^c = rh_symbol, e.g. for a doublet
        {X[s_, c_, 1] :> -XRwbar[s, c, 2], X[s_, c_, 2] :> XRwbar[s, c, 1]}.'''
        conj = self.conjugateMultiplet
        shared = ['s'] + (['c'] if self.colordim.dim > 1 else [])
        defs = []
        for f in self.fields:
            lhs, rhs = [f'{i}_' for i in shared], list(shared)
            if self.isospin.dim > 1:
                lhs.append(str(getSU2Index(f)))
                rhs.append(str(getSU2Index(f.conjugateField)))
            sign = '-' if conj.phase(f.conjugateField.i3) == -1 else ''
            defs.append(f'{self.fr_symbol()}[{", ".join(lhs)}] :> '
                        f'{sign}{rh_symbol}{DAGGER_SUFFIX}[{", ".join(rhs)}]')
        return '{' + ', '.join(defs) + '}'

    @property
    def components(self) -> Dict[U1Rep, Union[Field, None]]:
        """
        Read-only view of the multiplet components (a copy, so it cannot be
        used to modify the multiplet).
        """
        return dict(self._components)

    @property
    def labels(self) -> Dict[U1Rep, str]:
        """Read-only view of the component labels (a copy)."""
        return dict(self._labels)

    def phase(self, i3) -> int:
        """
        Phase with which the component field at ``i3`` enters this multiplet:
        component(i3) = phase(i3) * self[i3]. It is 1 for every multiplet built from
        fields, and (-1)^(j+m) for the conjugate multiplet, whose component -m is
        (-1)^(j+m) times the conjugate of the original component m
        (tools.reps.su2_conj_phase; e.g. H~ = i tau^2 H*: H~_2 = -H^-).
        Conjugating twice returns the original multiplet object (phases 1): the
        (-1)^(2j) of (psi^c)^c = W W* psi is not tracked (never needed).
        """
        return self._phases[U1Rep(i3)]

    def conjugate_component(self, i3) -> Tuple[Field, int]:
        '''(field, phase) with (component i3)^* = phase * field: for real multiplets with
        isospin dim > 1 the component -i3 (tools.reps.su2_reality_phase), otherwise the
        conjugate field with phase 1.'''
        if self.is_real_nontrivial:
            return self.components[U1Rep(-i3)], su2_reality_phase(self.isospin, i3)
        return self.components[U1Rep(i3)].conjugateField, 1

    def mass_partner(self, other: 'Multiplet', i3) -> Tuple[Field, int]:
        '''(field, phase) paired with component i3 = m in the mass term of self and other,
        the one convention of the mass terms (docs/conventions.md, Mass terms):
            -L = M sum_m phase self[m] field + h.c.   (1/2 M if other is self).
        Scalars: the conjugate component (conjugate_component), -L = M^2 X* X.
        Fermions: the component -m of other with phase c_m = s (-1)^(j+m) (su2_conj_phase),
        which makes the sum an SU(2) singlet; s = (-1)^(j+m0) gives c_m0 = +1 to the pair
        m0: the neutral component for a Majorana mass (other is self), the component of
        self with the largest |Q| (the larger m on a tie) for a Dirac mass. The other pairs
        have c_m = (-1)^(m-m0), e.g. -M for the charged pair of a Majorana triplet. A
        positive M gives positive masses; the negative entries, and negative inputs, are
        absorbed by the phases of the mass eigenstates (MultipletModel.mass_phases, and the
        Takagi or singular value decomposition of the mixing blocks).'''
        if isinstance(self.rep.spin, Scalar):
            return self.conjugate_component(i3)
        if other is self:
            lead = U1Rep(0)
        else:
            lead = max(self.components, key=lambda m: (abs(self.components[m].eCharge), m))
        sign = su2_conj_phase(self.isospin, lead)
        return other.components[U1Rep(-i3)], sign * su2_conj_phase(self.isospin, i3)

    @property
    def fields(self) -> List[Field]:
        return [f for f in self._components.values() if f is not None]

    @property
    def conjugateMultiplet(self) -> "Multiplet":
        return self if self.selfConjugate else self._conjugateMultiplet

    @conjugateMultiplet.setter
    def conjugateMultiplet(self, value):
        """
        Make sure the conjugateMultiplet attribute is protected
        """
        raise AttributeError("The conjugate multiplet of a multiplet is automatically defined when the multiplet is created "
                             "and it should not be changed directly.")

    def __getitem__(self, i3) -> Union[Field, None]:
        return self._components[i3]

    def __setitem__(self, i3, field):
        raise AttributeError("The components of a multiplet are fixed when the multiplet is created "
                             "and they should not be changed directly.")

    def __iter__(self):
        return iter(self._components)

    def __len__(self) -> int:
        return len(self._components)

    def keys(self) -> List[U1Rep]:
        return list(self._components)

    def values(self) -> List[Union[Field, None]]:
        return [self._components[i3] for i3 in self._components]

    def items(self):
        return zip(self.keys(), self.values())

    def is_filled(self) -> bool:
        '''True if every component has a field.'''
        return all(comp is not None for comp in self._components.values())

    def _add_field(self, field: Field) -> None:
        '''Place field at its component, checking its quantum numbers.'''
        assert field.isSM == self.isSM
        assert self.spin == field.spin
        assert self.colordim == field.colordim
        assert self.isospin == field.isospin
        assert self.hypercharge == field.hypercharge
        assert self.globalU1s == field.globalU1s
        i3 = field.i3
        assert i3 is not None
        assert (i3 in self._components) and (self._components[i3] is None)
        self._components[i3] = field
        if i3 not in self._labels:
            self._labels[i3] = field.label
        return None

    def _fill(self) -> None:
        """Fills the empty components of the multiplet in place with new, exotic fields.
        If the multiplet is already full, the function does nothing.
        """
        for i3, comp in self._components.items():
            if comp is not None:
                continue
            if i3 in self._labels:
                newLabel = self._labels[i3]
            else:
                Multiplet._EXO_COUNTER += 1
                newLabel = Multiplet._EXOTIC_LABEL + str(Multiplet._EXO_COUNTER)

            newField = Field(rep=BrokenRep.from_qnumbers(spin=self.spin,
                                                         colordim=self.colordim,
                                                         isospin=self.isospin,
                                                         hypercharge=self.hypercharge,
                                                         i3=i3,
                                                         globalU1s=self.globalU1s),
                             label=newLabel,
                             index=self.index,
                             isSM=self.isSM)

            self._add_field(newField)

    def __deepcopy__(self, memo):
        if id(self) in memo:
            return memo[id(self)]

        new_multiplet = Multiplet.__new__(Multiplet)
        memo[id(self)] = new_multiplet

        for attr in self.__slots__:
            if hasattr(self, attr):
                val = getattr(self, attr)
                setattr(new_multiplet, attr, copy.deepcopy(val, memo))

        return new_multiplet
