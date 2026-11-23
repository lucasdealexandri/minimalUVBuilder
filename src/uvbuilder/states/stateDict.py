"""
.. module:: stateDict
   :synopsis: Defines a dictionary of particles/final sates/fields which stores particles and not their conjugates.
              However, it also allows to retrieve the conjugate particle/final state/field by accessing 
              the antiParticle attribute of the particle/final state or the conjugateField attribute of the field.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from typing import Union, Any
from uvbuilder.states.particle import Particle
from uvbuilder.states.finalState import FinalState
from uvbuilder.states.field import Field
States = Union[Particle, FinalState, Field]


class StateDict(dict):
    """
    An instance of this class represents a dictionary for storing particles/final states/fields.

    The properties are: label, pdg, mass, spin, electric charge, colordim, width
    """

    def __init__(self, other=(), /, **kwargs):
        """
        Initialize the dictionary while enforcing canonical storage rules.
        """

        dict.__init__(self)
        self.update(other, **kwargs)

    

    @classmethod
    def _get_conjugate(cls, state: States) -> Union[None,States]:
        """
        Return the conjugate object if available.
        """

        if hasattr(state, "antiParticle"):
            return getattr(state, "antiParticle")
        if hasattr(state, "antiFinalState"):
            return getattr(state, "antiFinalState")
        if hasattr(state, "conjugateField"):
            return getattr(state, "conjugateField")
        return None

    @classmethod
    def _canonical_state(cls, state: States) -> States:
        """
        Return a canonical representative (particle/field, not anti-particle/conjugate) for storage.
        Particles are defined as the ones with pdg > 0.
        FinalStates are defined as the ones with more positive than negative pdgs.
        Fields are defined as the ones with label <= conjugate label.
        """

        if state.selfConjugate:
            return state
        else:
            conjugate = cls._get_conjugate(state)
            if conjugate is None:
                raise ValueError(f"Trying to create stateDict, but conjugate for {state} is not defined.")

        if isinstance(state, Particle) and state.pdg < 0:
            return conjugate
        elif isinstance(state, FinalState):
            pdgList = state.getPDGs()
            if sum([pdg < 0 for pdg in pdgList]) > sum([pdg > 0 for pdg in pdgList]):
                return conjugate
        elif isinstance(state, Field):
            if state.label > conjugate.label:
                return conjugate

        return state

    def __setitem__(self, key: int | str, state: States) -> None:
        """
        Store only the canonical state. 
        For integer keys, keep only the positive key, for string keys, keep only the label of the canonical state.
        """

        if not isinstance(key,(str,int)):
            raise ValueError(f"Cannot store keys which are not PDGs (int) or labels (str).")    

        canonical_state = self._canonical_state(state)
        if isinstance(key, int):
            canonical_key = abs(key)
        elif isinstance(key, str) and isinstance(getattr(canonical_state, "label", None), str):
            canonical_key = canonical_state.label
        else:
            canonical_key = key
        dict.__setitem__(self, canonical_key, canonical_state)

    def __getitem__(self, key: int | str) -> States:
        """
        Retrieve canonical states directly. For negative integer keys, return the conjugate.
        """

        if dict.__contains__(self, key):
            return dict.__getitem__(self, key)

        if isinstance(key, int) and key < 0:
            state = dict.__getitem__(self, abs(key))
            conjugate = self._get_conjugate(state)
            if conjugate is not None:
                return conjugate
            raise KeyError(key)

        elif isinstance(key, str):
            for state in dict.values(self):
                conjugate = self._get_conjugate(state)
                if conjugate is None:
                    continue
                if getattr(conjugate, "label", None) == key:
                    return conjugate

        raise KeyError(key)

    def get(self, key: str | int, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key: Any) -> bool:
        if dict.__contains__(self, key):
            return True
        if isinstance(key, int) and key < 0:
            return dict.__contains__(self, abs(key))
        if isinstance(key, str):
            for state in dict.values(self):
                conjugate = self._get_conjugate(state)
                if conjugate is None:
                    continue
                if getattr(conjugate, "label", None) == key:
                    return True
        return False

    def setdefault(self, key: str | int, default: Any = None) -> Any:
        if key in self:
            return self[key]
        self[key] = default
        return self[key]

    def update(self, other=(), /, **kwargs) -> None:
        if hasattr(other, "items"):
            iterable = other.items()
        else:
            iterable = other

        for key, value in iterable:
            self[key] = value

        for key, value in kwargs.items():
            self[key] = value

    def copy(self) -> "StateDict":
        """
        Creates a copy of the StateDict, where all states are also copied.
        """
        newDict = StateDict()
        for key, state in self.items():
            newDict[key] = state.copy()
        return newDict
