"""
.. module:: reps
   :synopsis: Defines the representations: U(1), global U(1), SU(2), SU(3) and Lorentz, and those of states and multiplets.

.. moduleauthor:: John Gargalionis <johngargalionis@gmail.com>
.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""
from __future__ import annotations

import numpy as np
from collections import namedtuple
from dataclasses import dataclass, replace
from fractions import Fraction
from typing import List, Tuple, Sequence
from typing import Union, Set, Dict, Optional
from functools import lru_cache
from itertools import product

# U(1), SU(2), SU(3) reps

class U1Rep(Fraction):
    '''
    Generic U1 representation
    '''

    def __new__(cls,numerator=0, denominator=None):
        frac = Fraction.__new__(cls,numerator,denominator)
        frac = frac.limit_denominator(100)
        self = super(U1Rep, cls).__new__(cls,frac)
        return self

    def __repr__(self):
        return str(self)

    def __neg__(self):
        return U1Rep(-self.numerator, self.denominator)
    
    def __sub__(self, other: U1Rep) -> U1Rep:
        r = Fraction.__sub__(self,other)
        r = U1Rep(r._numerator,r._denominator)
        return r
    
    def __mul__(self, other: U1Rep | int | float) -> U1Rep:
        r = Fraction.__mul__(self,other)
        r = U1Rep(r._numerator,r._denominator)
        return r
    
    def __rmul__(self, other: U1Rep | int | float) -> U1Rep:
        r = Fraction.__rmul__(self,other)
        r = U1Rep(r._numerator,r._denominator)
        return r
    
    def __add__(self, other: U1Rep) -> U1Rep:
        r = Fraction.__add__(self,other)
        r = U1Rep(r._numerator,r._denominator)
        return r
    
    def __radd__(self, other: U1Rep) -> U1Rep:        
        r = Fraction.__radd__(self,other)
        r = U1Rep(r._numerator,r._denominator)
        return r

    @property
    def conjRep(self):
        return -self
    
    @property
    def dim(self) -> int:
        '''
        Dimension of U1Rep (=1)
        '''
        return 1



@dataclass(frozen=True)
class GlobalU1s:
    '''
    Global U(1) charges: baryon number B, lepton number L and a dark charge D.

    A charge set to None is not tracked in the model; a charge set to 0 is
    tracked (conserved) and zero. Which charges are tracked is decided per
    model by the DM candidate (see ParticleModel.conserved); SM states carry
    all their charges and are projected onto the tracked ones where they meet
    BSM states.

    Example:
        GlobalU1s(B=U1Rep(1, 3))      # quark-like, only B tracked
        GlobalU1s(L=0)                # L conserved, this state has L = 0
    '''
    B: Optional[U1Rep] = None
    L: Optional[U1Rep] = None
    D: Optional[U1Rep] = None

    NAMES = ('B', 'L', 'D')

    def __post_init__(self):
        for name in self.NAMES:
            value = getattr(self, name)
            if value is not None and not isinstance(value, U1Rep):
                value = Fraction(value)
                object.__setattr__(self, name, U1Rep(value.numerator, value.denominator))

    @property
    def tracked(self) -> frozenset:
        '''Names of the tracked charges.'''
        return frozenset(n for n in self.NAMES if getattr(self, n) is not None)

    def project(self, names) -> GlobalU1s:
        '''Track exactly the charges in `names`: charges outside `names` become
        untracked, charges in `names` this state does not carry are set to 0
        (a state is neutral under a charge it does not carry, e.g. the Z under L).'''
        return GlobalU1s(**{n: getattr(self, n) or U1Rep(0) for n in self.NAMES if n in names})

    def items(self):
        '''(name, charge) for the tracked charges, in B, L, D order.'''
        return [(n, getattr(self, n)) for n in self.NAMES if getattr(self, n) is not None]

    def __getitem__(self, name: str) -> Optional[U1Rep]:
        return getattr(self, name)

    @property
    def is_zero(self) -> bool:
        '''True if every tracked charge is zero (vacuously true if none is tracked).'''
        return all(v == 0 for _, v in self.items())

    def _combine(self, other: GlobalU1s, sign: int) -> GlobalU1s:
        '''Sum (sign 1) or difference (sign -1) of the charges; both must track the same
        charges.'''
        if not isinstance(other, GlobalU1s):
            return NotImplemented
        if self.tracked != other.tracked:
            raise ValueError(
                f'Cannot combine global charges with different tracked sets: {self} and {other}. '
                f'Project both onto the conserved charges first.')
        return GlobalU1s(**{n: getattr(self, n) + sign * getattr(other, n) for n in self.tracked})

    def __add__(self, other: GlobalU1s) -> GlobalU1s:
        return self._combine(other, +1)

    def __sub__(self, other: GlobalU1s) -> GlobalU1s:
        return self._combine(other, -1)

    def __neg__(self) -> GlobalU1s:
        return GlobalU1s(**{n: -v for n, v in self.items()})

    @property
    def conjRep(self) -> GlobalU1s:
        return -self

    def _key(self) -> tuple:
        # total order: untracked sorts after tracked, then by value
        return tuple((getattr(self, n) is None, getattr(self, n) or 0) for n in self.NAMES)

    def __lt__(self, other: GlobalU1s) -> bool:
        if not isinstance(other, GlobalU1s): return NotImplemented
        return self._key() < other._key()

    def __le__(self, other: GlobalU1s) -> bool:
        if not isinstance(other, GlobalU1s): return NotImplemented
        return self._key() <= other._key()

    def __gt__(self, other: GlobalU1s) -> bool:
        if not isinstance(other, GlobalU1s): return NotImplemented
        return self._key() > other._key()

    def __ge__(self, other: GlobalU1s) -> bool:
        if not isinstance(other, GlobalU1s): return NotImplemented
        return self._key() >= other._key()

    def __repr__(self) -> str:
        return '(' + ', '.join(f'{n}={v}' for n, v in self.items()) + ')'


class SU2Rep(namedtuple("SU2", "dynkin")):
    '''SU(2) representation, by its Dynkin index (dimension dynkin + 1).'''

    def __repr__(self):
        return f"{self.dim}"

    @property
    def dim(self) -> int:
        '''
        Compute the SU2 dimension based on the
        dynkin index.
        '''
        d = self.dynkin+1
        
        return d
    
    def i3_range(self, reverse: bool = False) -> Tuple:
        '''
        Physical i3 eigenvalues of the rep, e.g. (-1/2, +1/2) for a doublet.
        '''
        r = tuple(U1Rep(i3, 2) for i3 in range(-self.dynkin, self.dynkin+1, 2))
        if reverse:
            r = r[::-1]
        return r

    @property
    def conjRep(self):
        return self
        
    def __mul__(self, other: SU2Rep) ->List[SU2Rep]:
        return su2_prod(self,other)

    def __rmul__(self, other: SU2Rep) ->List[SU2Rep]:
        return self.__mul__(other)
    

class SU3Rep(namedtuple("SU3", "upper lower")):
    '''SU(3) representation, by its Dynkin labels (upper, lower): 3 is (1, 0), 3bar (0, 1).'''

    def __repr__(self):        
        if self.upper == 0 and self.lower != 0:
            return f"{self.dim}*"
        else:
            return f"{self.dim}"
    
    @property
    def dim(self) -> int:
        '''
        Compute the SU3 dimension based on the
        number of upper and lower indices.
        '''
        m = self.upper
        n = self.lower
        d = int(0.5*(m+1)*(n+1)*(m+n+2))
        return d

    # maybe remove this one...
    def inv(self):
        '''
        Return the conjugate representation (upper <-> lower)
        '''

        return SU3Rep(upper=self.lower,lower=self.upper)    

    @property
    def conjRep(self):
        return SU3Rep(upper=self.lower,lower=self.upper)
    
    def __mul__(self, other: SU3Rep) ->List[SU3Rep]:
        return su3_prod(self,other)

    def __rmul__(self, other: SU3Rep) ->List[SU3Rep]:
        return self.__mul__(other)
        

# Umbrella Class for Unbroken and Broken Reps. These describe reps under the full group.
@dataclass(frozen=True, kw_only=True)
class Rep:
    '''Shared gauge data of both descriptions'''
    spin: Optional[LorentzRep] = None
    colordim: Optional[SU3Rep] = None
    globalU1s: GlobalU1s = GlobalU1s()

    @property
    def conjRep(self) -> Rep:
        '''Conjugate of the full rep: every factor conjugated by its own conjRep'''
        return replace(self, **{
            name: value.conjRep
            for name, value in vars(self).items()
            if value is not None and hasattr(value, 'conjRep')
        })

    @property
    def conjInternalRep(self) -> Rep:
        '''Conjugates every rep except for spin. Useful for Spinors'''
        return replace(self, **{
            name: value.conjRep
            for name, value in vars(self).items()
            if value is not None and hasattr(value, 'conjRep') and name != 'spin'
        })

    @property
    def is_realInternalRep(self) -> bool:
        '''Checks if a given rep is real (equal to its conjugate). Does not include Spin'''
        return self == self.conjInternalRep

    @property
    def is_singlet(self) -> bool:
        '''Checks if given Rep is scalar (singlet under all groups)'''
        scalar_rep: Dict[str, AllReps] = {
            'spin': SCALAR,
            'colordim': SU3Rep(0,0),
            'eCharge': U1Rep(0,),
            'isospin': SU2Rep(0,),
            'hypercharge': U1Rep(0,),
            'i3': U1Rep(0,)}

        # global U(1)s: a singlet has every *tracked* charge equal to zero
        if self.globalU1s is not None and not self.globalU1s.is_zero:
            return False
        return all(
            getattr(self, name) == scalar
            for name, scalar in scalar_rep.items()
            if name in type(self).__dataclass_fields__ and
            getattr(self, name) is not None
        )

    def scalar_product(self, other: Rep | List[Rep]) -> List[Rep]:
        '''Spin-scalar product: contracts spin to SCALAR, multiplies the gauge reps.'''
        return self._product(other, scalar=True)

    def tensor_product(self, other: Rep | List[Rep]) -> List[Rep]:
        '''Full tensor product: every factor (spin included) multiplies.'''
        return self._product(other, scalar=False)

    def _product(self, other, *, scalar: bool) -> List[Rep]:
        '''The representations in the products with other (one representation or several).'''
        others = other if isinstance(other, Sequence) else [other]
        for o in others:
            if not (isinstance(o, type(self)) and isinstance(self, type(o))):
                raise TypeError(f'Can only multiply {type(self).__name__} with '
                                f'{type(self).__name__}; got {type(o).__name__}.')
        results = []
        for o in others:
            factors = self._combine_factors(o, scalar=scalar)
            if factors is None:
                continue
            names = list(factors)
            results += [
                type(self)(**dict(zip(names, combo)))
                for combo in product(*(factors[name] for name in names))
            ]
        return results

    def _combine_factors(self, other: Rep, *, scalar: bool) -> Optional[Dict[str, List]]:
        '''One combined value-list per quantum number of `type(self)`.
        Returns None when a quantum number is undefined, i.e. the pair
        produces no terms. Charges add; irreps multiply into a list; the
        caller builds the Cartesian product over the per-factor lists.'''
        factors: Dict[str, List] = {}
        for name in type(self).__dataclass_fields__:
            va = getattr(self, name)
            vb = getattr(other, name)
            if name == 'spin':
                if scalar:
                    if va != vb:
                        raise ValueError('Can only scalar multiply reps with same spin.')
                    factors[name] = [SCALAR]
                elif va is None or vb is None:
                    return None
                else:
                    factors[name] = list(va * vb)
            elif name == 'i3':
                assert (va is None) == (vb is None), \
                    'Can only multiply two field reps or two multiplet reps.'
                factors[name] = [None] if va is None else [va + vb]
            elif va is None or vb is None:
                return None
            elif isinstance(va, (U1Rep, GlobalU1s)):
                factors[name] = [va + vb]
            else:
                factors[name] = list(va * vb)
        return factors

    @property
    def canonical_order(self) -> tuple:
        return (self.spin, self.colordim, self.globalU1s)

    def __mul__(self, other: Rep | List[Rep]) -> List[Rep]:
        return self.tensor_product(other)

    def __rmul__(self, other):
        return self.__mul__(other)

# Rep under Lorentz x SU(3)_c x U(1)_EM X prod_X U(1)_X
@dataclass(frozen=True, kw_only=True)
class UnbrokenRep(Rep):
    '''A single component: definite electric charge, multiplet structure unknown.'''
    eCharge: Optional[U1Rep] = None

    # TODO: maybe keep track of energy_dim on this new Scalar? 
    def __mul__(self, other: UnbrokenRep) -> List[UnbrokenRep]:
        return self.scalar_product(other)

    def can_connect(self, other: UnbrokenRep) -> bool:
        '''True if other can be a neighboring component of the same multiplet: same
        spin, color and global charges, electric charge differing by one.'''
        return (self.spin == other.spin and
                self.colordim == other.colordim and
                self.globalU1s == other.globalU1s and
                abs(self.eCharge - other.eCharge) == U1Rep(1,))

    @classmethod
    def from_particle(cls, particle) -> UnbrokenRep:
        '''The representation of a particle.'''
        return cls(
            spin=particle.spin,
            colordim=particle.colordim,
            eCharge=particle.eCharge,
            globalU1s=particle.globalU1s or GlobalU1s(),
        )

    @property
    def canonical_order(self):
        return super().canonical_order + (self.eCharge,)

    def __lt__(self, other: UnbrokenRep):
        if not isinstance(other, UnbrokenRep): return NotImplemented
        return self.canonical_order.__lt__(other.canonical_order)

    def __le__(self, other: UnbrokenRep):
        if not isinstance(other, UnbrokenRep): return NotImplemented
        return self.canonical_order.__le__(other.canonical_order)

    def __gt__(self, other: UnbrokenRep):
        if not isinstance(other, UnbrokenRep): return NotImplemented
        return self.canonical_order.__gt__(other.canonical_order)

    def __ge__(self, other: UnbrokenRep):
        if not isinstance(other, UnbrokenRep): return NotImplemented
        return self.canonical_order.__ge__(other.canonical_order)

    def __eq__(self, other: UnbrokenRep):
        if not isinstance(other, UnbrokenRep): return NotImplemented
        return self.canonical_order.__eq__(other.canonical_order)

    def __repr__(self) -> str:
        adUQNs = f', {self.globalU1s}' if self.globalU1s and self.globalU1s.tracked else ''
        color = f', c: {self.colordim}' if self.colordim != SU3Rep(0,0) else ''
        return f'({self.spin}, Q: {self.eCharge}{color}{adUQNs})'


# Rep under Lorentz x SU(3)_c x SU(2)_L x U(1)_Y X prod_X U(1)_X
@dataclass(frozen=True, kw_only=True)
class BrokenRep(Rep):
    '''An SU(2)_L x U(1)_Y multiplet; i3 set -> one component of it.'''
    isospin: SU2Rep
    hypercharge: U1Rep
    i3: Optional[U1Rep] = None

    def __post_init__(self) -> None:
        # Normalize i3 to a physical U1Rep value (e.g. -1/2 for a doublet component).
        if self.i3 is not None and not isinstance(self.i3, U1Rep):
            object.__setattr__(self, 'i3', U1Rep(self.i3))

    @property
    def is_multiplet_rep(self) -> bool:
        return self.i3 is None

    @property
    def is_field_rep(self) -> bool:
        return self.i3 is not None

    @property
    def multiplet_rep(self) -> BrokenRep:
        return replace(self, i3=None)

    def same_multiplet(self, other: BrokenRep) -> bool:
        '''True if both are components of the same multiplet.'''

        if self.multiplet_rep == other.multiplet_rep:
            return True

        return False


    @property
    def eCharge(self) -> Optional[U1Rep]:
        return None if self.i3 is None else self.i3 + self.hypercharge

    @property
    def components(self) -> List[BrokenRep]:
        assert self.is_multiplet_rep, 'Field reps are themselves "components", they do not contain components.'
        return [replace(self, i3=i3) for i3 in self.isospin.i3_range(reverse=True)]

    @property
    def unbroken(self) -> UnbrokenRep:
        '''The unbroken representation of this component.'''
        if self.i3 is None:
            raise ValueError('A whole multiplet (i3=None) has no definite electric charge.')
        return UnbrokenRep(spin=self.spin, colordim=self.colordim,
                           eCharge=self.eCharge,
                           globalU1s=self.globalU1s)

    @classmethod
    def from_block(cls, block: List[UnbrokenRep]) -> BrokenRep:
        '''The multiplet whose components are the chain block: isospin j = (len - 1)/2,
        hypercharge Q_max - j.'''
        anchor = block[0]
        max_charge = max(rep.eCharge for rep in block)
        isospin = SU2Rep(len(block) - 1,)
        hypercharge = max_charge - U1Rep(isospin.dynkin, 2)
        return cls(spin=anchor.spin, colordim=anchor.colordim,
                   globalU1s=anchor.globalU1s,
                   isospin=isospin, hypercharge=hypercharge)

    @classmethod
    def from_qnumbers(cls, *, spin: Optional[LorentzRep] = None,
                      colordim: Optional[SU3Rep] = None,
                      isospin: SU2Rep = SU2Rep(0,),
                      hypercharge: Optional[U1Rep] = None,
                      i3: Optional[Union[int, U1Rep]] = None,
                      eCharge: Optional[U1Rep] = None,
                      globalU1s: Optional[GlobalU1s] = None) -> BrokenRep:
        '''
        Build a complete BrokenRep from field quantum numbers.

        i3 is the physical value, given as a U1Rep (U1Rep(-1,2) = -1/2) or as
        an int/float/Fraction of the same (e.g. 0, -1/2, 1). The missing one of
        {eCharge, i3, hypercharge} is derived from the other two via
        eCharge = i3 + hypercharge (Gell-Mann-Nishijima).
        '''
        if globalU1s is None:
            globalU1s = GlobalU1s()

        i3_rep = None if i3 is None else U1Rep(i3)

        known = sum(x is not None for x in (eCharge, i3_rep, hypercharge))
        if known == 3:
            if eCharge != i3_rep + hypercharge:
                raise AttributeError(
                    f"Inconsistent electroweak quantum numbers: eCharge={eCharge}, "
                    f"hypercharge={hypercharge}, i3={i3_rep}. Expected eCharge={i3_rep + hypercharge}")
        elif known == 2:
            if hypercharge is None:
                hypercharge = eCharge - i3_rep
            elif i3_rep is None:
                i3_rep = eCharge - hypercharge
                if (i3_rep * 2).denominator != 1:
                    raise AttributeError(
                        f"Cannot derive half-integer i3 from eCharge={eCharge} and hypercharge={hypercharge}")
            else:
                eCharge = i3_rep + hypercharge
        else:
            raise AttributeError(
                f"Not enough electroweak quantum numbers to define a complete BrokenRep "
                f"(need at least two of eCharge, i3, hypercharge); got eCharge={eCharge}, "
                f"i3={i3_rep}, hypercharge={hypercharge}")

        return cls(spin=spin, colordim=colordim, globalU1s=globalU1s,
                   isospin=isospin, hypercharge=hypercharge, i3=i3_rep)

    @property
    def canonical_order(self):
        order = super().canonical_order + (self.isospin, self.hypercharge)
        if self.i3 != None: order = order + (self.eCharge,)
        return order

    def __lt__(self, other: BrokenRep):
        if not isinstance(other, BrokenRep): return NotImplemented
        return self.canonical_order.__lt__(other.canonical_order)

    def __le__(self, other: BrokenRep):
        if not isinstance(other, BrokenRep): return NotImplemented
        return self.canonical_order.__le__(other.canonical_order)

    def __gt__(self, other: BrokenRep):
        if not isinstance(other, BrokenRep): return NotImplemented
        return self.canonical_order.__gt__(other.canonical_order)

    def __ge__(self, other: BrokenRep):
        if not isinstance(other, BrokenRep): return NotImplemented
        return self.canonical_order.__ge__(other.canonical_order)

    def __eq__(self, other: BrokenRep):
        if not isinstance(other, BrokenRep): return NotImplemented
        return self.canonical_order.__eq__(other.canonical_order)

    def __repr__(self) -> str:
        adUQNs = f', {self.globalU1s}' if self.globalU1s and self.globalU1s.tracked else ''
        i3 = f', i3: {self.i3}' if self.i3 is not None else ''
        color = f', c: {self.colordim}' if self.colordim != SU3Rep(0,0) else ''
        return f'({self.spin}, I: {self.isospin}, Y: {self.hypercharge}{i3}{color}{adUQNs})'


# Aux functions for su3 multiplication

# Multiplication is implemented using the tensor method for SU3. Input tensors
# to be multiplied are assumed to be traceless and symmetric in upper and lower
# indices (but not across indices). First, traces on indices across the tensors
# are taken using deltas, then all possible contractions with the epsilon
# symbol are made for both the raised and lowered case.

# The function su3_multiplication_rule_deltas returns a list of intermediate
# tuple objects containing the index structure of the multiplied indices for
# the two tensors: (tensor1_raised, tensor1_lowered, tensor2_raised,
# tensor2_lowered) and the number of traces needed to get to that index
# structure. These traces are of two types: (1) a delta with a lowered index
# that is raised on the first tensor and a raised index that is lowered on the
# second tensor, and (2) a delta with a lowered index that is raised on the
# second tensor and a raised index that is lowered on the first tensor. The
# numbers of these deltas required is returned as a regular tuple of integers,
# given in the order presented in the previous sentence.

# The function su3_multiplication_rule_epsilons returns a list of EpsTensorMul
# tuples, that contain the final representations from the multiplication along
# with the deltas and epsilons required to attain them. As with the deltas, the
# epsilons are represented as a tuple of integers. The first entry corresponds
# to the number of epsilons with lowered indices required to construct the
# representation, while the second is the number of epsilons with raised
# indices. From this object, the explicit tensor corresponding to the
# representation can be constructed.

DeltaTensorMul = namedtuple("DeltaTensorMul", ["tensor", "deltas"])
EpsTensorMul = namedtuple("EpsTensorMul", ["tensor", "deltas", "epsilons"])


def su3_multiplication_rule_deltas(indices):
    '''Returns a list of tuples with the structure [(tensor, (deltas)), ...] where
    the tensors are represented as a numpy array of indices.

    '''
    #         first kind of delta,      second kind
    deltas = [np.array([-1, 0, 0, -1]), np.array([0, -1, -1, 0])]
    index_array = np.array(indices)
    max_deltas = min(indices[0] + indices[2], indices[1] + indices[3]) + 1

    delta_coeffs = []
    for a in range(max_deltas):
        for b in range(max_deltas):
            if np.all(index_array + a * deltas[0] + b * deltas[1] >= 0):
                delta_coeffs.append((a, b))

    tensor_indices = [
        index_array + a * deltas[0] + b * deltas[1] for a, b in delta_coeffs
    ]
    return [
        DeltaTensorMul(tensor=t, deltas=d) for t, d in zip(tensor_indices, delta_coeffs)
    ]


def su3_multiplication_rule_epsilons(delta_tensor_mul):
    tensor, deltas = delta_tensor_mul
    ups = tensor[0] + tensor[2]
    downs = tensor[1] + tensor[3]
    index_array = np.array([ups, downs])
    #      first kind of eps, second kind
    eps = [np.array([-2, 1]), np.array([1, -2])]
    eps1_coeffs, eps2_coeffs = [], []
    max_eps_down = min(tensor[0], tensor[2]) + 1
    max_eps_up = min(tensor[1], tensor[3]) + 1

    for a in range(max_eps_down):
        if np.all(index_array + a * eps[0] >= 0):
            eps1_coeffs.append(a)

    for a in range(max_eps_up):
        if np.all(index_array + a * eps[1] >= 0):
            eps2_coeffs.append(a)

    # the number of indices you would have after multiplying the different epsilons
    tensor_indices_eps1 = [(index_array + a * eps[0]).tolist() for a in eps1_coeffs]
    tensor_indices_eps2 = [(index_array + a * eps[1]).tolist() for a in eps2_coeffs]

    result = []
    for tensor, eps in zip(tensor_indices_eps1, eps1_coeffs):
        result.append(EpsTensorMul(tensor=tuple(tensor), deltas=deltas, epsilons=(eps, 0)))

    for tensor, eps in zip(tensor_indices_eps2, eps2_coeffs):
        if eps == 0:
            # don't count trivial case twice
            continue
        result.append(EpsTensorMul(tensor=tuple(tensor), deltas=deltas, epsilons=(0, eps)))

    return result


def su3_multiplication_rule(tensor1_idxs, tensor2_idxs):
    '''Implementation of procedure followed in equation X in Zee V.2.11.'''
    indices = (*tensor1_idxs, *tensor2_idxs)
    intermediate_result = su3_multiplication_rule_deltas(indices)
    return [c for r in intermediate_result for c in su3_multiplication_rule_epsilons(r)]


@lru_cache
def su3_prod(
    a: SU3Rep, b: SU3Rep) -> List[SU3Rep]:
    '''In principle the function `su3_multiplication_rule` returns enough
    information to construct the Lagrangian term. We don't need all of that
    information for now.

    '''
    return [
        SU3Rep(*(x.tensor)) for x in su3_multiplication_rule(a, b)
    ]


# ----------------------------------------------------------------------
# SU(2) conjugation convention (single source of truth)
# ----------------------------------------------------------------------
#
# Basis of every SU(2) multiplet: the charge (i3) eigenstates |j, m>, index 1 = highest m,
# with Condon-Shortley phases (J+ and J- have real positive matrix elements). The
# generators in this basis are the TSU2d<d> of SMMod.fr (and PauliSigma/2 for d = 2).
#
# Conjugation matrix. For a self-conjugate representation there is a unitary W with
# D(U)* = W D(U) W^-1, i.e. (i T^a)* = W (i T^a) W^-1; W is symmetric (W W* = 1) for real
# and antisymmetric (W W* = -1) for pseudo-real representations, and for the SU(2) doublet
# C^{ab} = (i tau^2)^{ab} = epsilon^{ab} [Dreiner, Haber, Martin, Phys. Rept. 494 (2010) 1,
# arXiv:0812.1594, App. E.1, eqs. (E.1.3)-(E.1.9)]. For general j the conjugation of the
# Wigner D matrices, D^j_{m'm}(R)* = (-1)^{m'-m} D^j_{-m',-m}(R) [Condon & Shortley, The
# Theory of Atomic Spectra (1935); Edmonds, Angular Momentum in Quantum Mechanics (1957),
# ch. 4; Varshalovich, Moskalev, Khersonskii, Quantum Theory of Angular Momentum (1988),
# sec. 4.4], gives W_{m,-m} proportional to (-1)^{j-m}; W is fixed only up to an overall phase.
# We use
#
#     (psi^c)_{-m} = (-1)^(j+m) (psi_m)^*          [su2_conj_phase]
#
# which for j = 1/2 is exactly i tau^2 (our H~ = i tau^2 H*, the Phitilde class of SMMod.fr),
# and is symmetric for integer j and antisymmetric for half-integer j, as (E.1.5)/(E.1.6) require.
#
# Real multiplets (integer j, Y = 0, no conserved charge) obey psi = eta W psi^*; eta is the
# overall-phase freedom of W and is fixed, eta = (-1)^j, so that the neutral component is real:
#
#     psi_{-m} = (-1)^m (psi_m)^*   e.g. T^- = -(T^+)^* for a real triplet   [su2_reality_phase]
#
# Gauge couplings: g W^a T^a = (g/sqrt 2)(W+ T+ + W- T-) + ..., T+- = T1 +- i T2
# [Dreiner, Haber, Martin, eq. (4.3.12)], i.e. W+ raises m by one unit.


def _integer(q: U1Rep, what: str) -> int:
    '''q as an int; ValueError naming what if it is not an integer.'''
    q = U1Rep(q)
    if q.denominator != 1:
        raise ValueError(f'{what} must be an integer, got {q}')
    return q.numerator


def su2_conj_phase(isospin: SU2Rep, i3) -> int:
    '''Phase (-1)^(j+m) with which (psi_m)^* enters the component -m of the conjugate
    multiplet: (psi^c)_{-m} = su2_conj_phase(isospin, m) (psi_m)^*. See the convention
    block above (Dreiner-Haber-Martin App. E.1; Condon-Shortley phases).'''
    n = _integer(U1Rep(isospin.dynkin, 2) + U1Rep(i3), 'j + m')
    return -1 if n % 2 else 1


def su2_reality_phase(isospin: SU2Rep, i3) -> int:
    '''For a real multiplet (integer j): psi_{-m} = su2_reality_phase(isospin, m) (psi_m)^*,
    i.e. (-1)^m, so that the neutral component is real (eta = (-1)^j).'''
    _integer(U1Rep(isospin.dynkin, 2), 'j (real multiplets need integer isospin)')
    n = _integer(U1Rep(i3), 'm')
    return -1 if n % 2 else 1


@lru_cache
def su2_prod(a: SU2Rep, b: SU2Rep) -> List[SU2Rep]:
    '''Function to calculate products of SU(2). The integers represent the dynkin
    indices of the SU(2) representations (add one for the dimension).

    '''
    decomposition = []
    difference = abs(a.dynkin - b.dynkin)
    while difference <= a.dynkin + b.dynkin:
        decomposition.append(SU2Rep(difference,))
        difference += 2
    return decomposition

@lru_cache
def weyl_spin_prod(spin1: LorentzRep, spin2: LorentzRep):
    '''List su2_prod but implements additional steps:
          (1) Removes exotic spins
          (2) Promotes Weyl-spinor reps to Dirac spinors

    '''
    # No inv needed for SU(2)
    tuple1 = tuple([s for s in spin1 if len(s) != 0])
    tuple2 = tuple([s for s in spin2 if len(s) != 0])
    
    product_sum: List[Tuple[SU2Rep, SU2Rep]] = su2su2_prod(tuple1,tuple2)

    # Mapping of allowed irreps to consistent representations
    allowed_irreps: Dict[Tuple[SU2Rep, SU2Rep], LorentzRep] = {
        (SU2Rep(1),SU2Rep(0)): WEYL_LEFT,
        (SU2Rep(0),SU2Rep(1)): WEYL_RIGHT,
        (SU2Rep(1),SU2Rep(1)): VECTOR,
        (SU2Rep(0),SU2Rep(0)): SCALAR,
    }

    # Don't allow spin 2, spin 3/2, or anything else exotic
    reduced_sum = [term for term in product_sum if term in allowed_irreps]

    # Promote every spinor to a Dirac fermion and remove multiplicity from the
    # list with set
    processed_sum: Set[LorentzRep] = set(allowed_irreps[term] for term in reduced_sum)

    # Return a list output that is sorted so as to not allow indeterminacy that might come from the set call above
    return sorted(processed_sum)


@lru_cache
def su2su2_prod(
        sum_a: Sequence[Tuple[SU2Rep, SU2Rep]],
        sum_b: Sequence[Tuple[SU2Rep, SU2Rep]]
) -> List[Tuple[SU2Rep, SU2Rep]]:
    """Function to calculate products of SU(2)_{+} \\otimes SU(2)_{-}, enumerating
    representations of the Lorentz group. The integers represent the dynkin
    indices of the SU(2) representations.

    Lists are used to represent terms in a direct sum. Use `itertools.product`
    to expand out these terms. This choice means the inputs cannot be cached
    though, so use additional, custom decorator to convert lists into tuples.

    E.g.
    ```
    su2su2_prod(((0,1),), ((0,1), (1,0)))
    # => [(0, 0), (0, 2), (1, 1)]
    ```

    """
    result = []
    for a,b in product(sum_a, sum_b):
        (a_1, a_2), (b_1, b_2) = a, b
        first_prod = (tup[0] for tup in a_1*b_1)
        second_prod = (tup[0] for tup in a_2*b_2)
        result += [(SU2Rep(x),SU2Rep(y)) for x,y in product(first_prod, second_prod)]

    return result

@lru_cache
def spin_prod(spin_1: LorentzRep, spin_2: LorentzRep) -> List[LorentzRep]:
    """
    Allowed products of two Lorentz reps. 
    Considers only products of LEFT_SPINOR, RIGHT_SPINOR, DIRAC_SPINOR, SCALAR and VECTOR.
    """

    
    allowedProduct = {
                        (SCALAR,SCALAR) : [SCALAR],                       
                        (SCALAR,VECTOR) : [SCALAR] 
                    }
    
    spinors = (LEFT_SPINOR, RIGHT_SPINOR, DIRAC_SPINOR)
    for s1,s2 in product(spinors,repeat=2):
        allowedProduct[(s1,s2)] = [SCALAR]
    for s in spinors:
        allowedProduct[(s,SCALAR)] = [DIRAC_SPINOR]
        allowedProduct[(s,VECTOR)] = [DIRAC_SPINOR]

                        
    allowedProduct = {tuple(sorted(k)) : v for k,v in list(allowedProduct.items())}

    spins = tuple(sorted([spin_1,spin_2]))
    return allowedProduct.get(spins,[])


def spin_prod_full(spin_1: LorentzRep, spin_2: LorentzRep) -> List[LorentzRep]:
    """
    List su2su2_prod but implements additional steps:
          (1) Removes exotic spins
          (2) Promotes Weyl-spinor reps to Dirac spinors
    It does not allow for SCALAR*VECTOR = SCALAR!

    """

    # No inv needed for SU(2)
    tuple1 = tuple([s for s in spin_1 if len(s) != 0])
    tuple2 = tuple([s for s in spin_2 if len(s) != 0])
    
    product_sum: List[Tuple[SU2Rep, SU2Rep]] = su2su2_prod(tuple1,tuple2)

    # Mapping of allowed irreps to consistent representations
    allowed_irreps: Dict[Tuple[SU2Rep, SU2Rep], LorentzRep] = {
        (SU2Rep(1),SU2Rep(0)): DIRAC_SPINOR,
        (SU2Rep(0),SU2Rep(1)): DIRAC_SPINOR,
        (SU2Rep(0),SU2Rep(0)): SCALAR,
    }

    # Don't allow spin 2, spin 3/2, or anything else exotic
    reduced_sum = [term for term in product_sum if term in allowed_irreps]

    # Promote every spinor to a Dirac fermion and remove multiplicity from the
    # list with set
    processed_sum: Set[LorentzRep] = set(allowed_irreps[term] for term in reduced_sum)

    # Return a list output that is sorted so as to not allow indeterminacy that might come from the set call above
    return sorted(processed_sum)

# Spin Classes 

def _swap(rep):
    '''Reverse a (j+, j-) pair; returns () for a missing Weyl component.'''
    return (rep[1], rep[0]) if rep else rep


def _spin_mul(self, other, prod):
    '''LorentzRep product against a single rep or a list of reps (flattened).'''
    if isinstance(other, (DiracSpinor, WeylSpinor, Scalar, Vector)):
        return list(prod(self, other))
    return [r for o in other for r in prod(self, o)]


class DiracSpinor(namedtuple("DiracSpinor", "left right")):
    '''Dirac spinor, (1/2, 0) + (0, 1/2).'''
    def __repr__(self):
        return "DiracSpinor"
        
    @property
    def dim(self) -> int:
        return 4

    @property
    def energy_dim(self) -> U1Rep:
        return U1Rep(3,2)

    @property
    def conjRep(self):
        return DiracSpinor(left=_swap(self.right), right=_swap(self.left))

    def __mul__(self, other: "LorentzRep") -> List["LorentzRep"]:
        return _spin_mul(self, other, spin_prod_full)
    
    def __rmul__(self,other: "LorentzRep") ->List["LorentzRep"]:
        return self.__mul__(other)
    
    
class WeylSpinor(namedtuple("WeylSpinor", "left right")):
    '''Weyl spinor: left-handed (1/2, 0) or right-handed (0, 1/2).'''

    def __repr__(self):
        chirality = "LH" if self.left else "RH"
        return f"{chirality} WeylSpinor"
        
    @property
    def dim(self) -> int:
        return 2
    
    @property
    def dynkin(self) -> int:
        return self.dim - 1

    @property
    def energy_dim(self) -> U1Rep:
        return U1Rep(3, 2)

    def scalar_product(self, other: WeylSpinor) -> Scalar:
        '''Lorentz scalar of two spinors of the same chirality.'''
        assert type(self) == type(other), 'Can only scalar multiply two spinors of the same type.'
        return Scalar

    def tensor_product(self, other: WeylSpinor) -> List[LorentzRep]:
        '''Tensor product of two spins, generates list of spin reps.'''
        return weyl_spin_prod(self, other)

    def __mul__(self, other: "LorentzRep") -> List["LorentzRep"]:
        '''Interpreted as the scalar product of 2 Weyl Spinors'''
        return _spin_mul(self, other, weyl_spin_prod)
    
    def __rmul__(self,other: "LorentzRep") ->List["LorentzRep"]:
        return self.__mul__(other)

    @property
    def conjRep(self):
        return WeylSpinor(left=_swap(self.right), right=_swap(self.left))
    

class Scalar(namedtuple("Scalar", "arg")):
    '''Scalar, (0, 0).'''
    def __repr__(self):
        return "Scalar"
    
    @property
    def dim(self) -> int:
        return 0

    @property
    def energy_dim(self) -> U1Rep:
        return U1Rep(1, 1)
    
    def __mul__(self, other: "LorentzRep") -> List["LorentzRep"]:
        return _spin_mul(self, other,
                         lambda a, b: weyl_spin_prod(b, a)
                         if isinstance(b, WeylSpinor) else spin_prod_full(a, b))
    
    def __rmul__(self,other: "LorentzRep") ->List["LorentzRep"]:
        return self.__mul__(other)

    @property
    def conjRep(self):
        return self
    

class Vector(namedtuple("Vector", "arg")):
    '''Vector, (1/2, 1/2).'''
    def __repr__(self):
        return "Vector"
    
    @property
    def dim(self) -> int:
        return 1

    @property
    def energy_dim(self) -> U1Rep:
        return U1Rep(1, 1)
    
    def __mul__(self, other: "LorentzRep") -> List["LorentzRep"]:
        return _spin_mul(self, other, spin_prod_full)
    
    def __rmul__(self,other: "LorentzRep") ->List["LorentzRep"]:
        return self.__mul__(other)

    @property
    def conjRep(self):
        return self

DIRAC_SPINOR = DiracSpinor(left=(SU2Rep(1),SU2Rep(0)), right=(SU2Rep(0),SU2Rep(1)))
SCALAR =       Scalar(arg=(SU2Rep(0),SU2Rep(0)))
VECTOR =       Vector(arg=(SU2Rep(1),SU2Rep(1)))
LEFT_SPINOR =  DiracSpinor(left=(SU2Rep(1),SU2Rep(0)), right=())
RIGHT_SPINOR = DiracSpinor(left=(), right=(SU2Rep(0),SU2Rep(1)))
WEYL_LEFT =    WeylSpinor(left=(SU2Rep(1),SU2Rep(0)), right=())
WEYL_RIGHT =   WeylSpinor(left=(), right=(SU2Rep(0),SU2Rep(1)))

# Define Lorentz reps
LorentzRep = Union[DiracSpinor, WeylSpinor, Scalar, Vector]
AllReps = Union[LorentzRep, SU3Rep, SU2Rep, U1Rep, int]
