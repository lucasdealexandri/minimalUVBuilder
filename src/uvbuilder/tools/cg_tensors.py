#!/usr/bin/env python3
'''SU(2) and SU(3) Clebsch-Gordan tensor objects.

The CG layer treats both groups symmetrically: a CGTensor knows its own
FeynRules symbol, the indices passed to it, and (SU(2) only) the numeric
singlet values that must be declared with a ``Definitions`` rule in
FeynRules. The SU(3) heads (K3, K6, T, ...) are FeynRules built-ins, so no
numeric data lives here (K3/K3bar alias the built-in ``Eps``).

Field/index assignment is NOT this module's job - the Operator knows how to
render a FeynRules term with its multiplets and distributes the indices.
'''

import warnings

import sympy as sp
from functools import lru_cache
from sympy.printing.mathematica import mathematica_code
from typing import Dict, List, Optional, Sequence, Tuple

from uvbuilder.tools.reps import U1Rep, SU2Rep, SU3Rep
from uvbuilder.tools.fr_utils import su2_index_type, su2_fr_index


# ----------------------------------------------------------------------
# SU(2) Clebsch-Gordan table (shared cache)
# ----------------------------------------------------------------------

def _cg_val(m: U1Rep) -> sp.Expr:
    '''Exact half-integer value of a U1Rep i3 charge, e.g. U1Rep(-1, 2) -> -1/2.'''
    return sp.Rational(m.numerator, m.denominator)


@lru_cache(maxsize=None)
def cg_block(j1: SU2Rep, j2: SU2Rep) -> Dict[Tuple[U1Rep, U1Rep, SU2Rep, U1Rep], sp.Expr]:
    '''All nonzero SU(2) Clebsch-Gordan coefficients for j1 x j2 -> J.

    Returns a dict {(m1, m2, J, M): exact sympy coefficient} with zeros omitted.
    Shared by every Operator with the same isospin pair.
    '''
    from sympy.physics.quantum.cg import CG
    j1h, j2h = sp.Rational(j1.dynkin, 2), sp.Rational(j2.dynkin, 2)
    block: Dict[Tuple[U1Rep, U1Rep, SU2Rep, U1Rep], sp.Expr] = {}
    for m1 in j1.i3_range():
        m1h = _cg_val(m1)
        for m2 in j2.i3_range():
            m2h = _cg_val(m2)
            for J in j1 * j2:
                Jh = sp.Rational(J.dynkin, 2)
                for M in J.i3_range():
                    Mh = _cg_val(M)
                    c = CG(j1h, m1h, j2h, m2h, Jh, Mh).doit()
                    if c != 0:
                        block[(m1, m2, J, M)] = sp.simplify(c)
    return block


@lru_cache(maxsize=None)
def su2_singlet_paths(isospins: Tuple[SU2Rep, ...]) -> List[Tuple[SU2Rep, ...]]:
    '''Intermediate spins (J12, J123, ...) of every successive-coupling path of the
    isospins that ends in a singlet. Each path gives one independent invariant.'''
    paths = [((), isospins[0])]
    for j in isospins[1:-1]:
        paths = [(path + (J,), J) for path, J0 in paths for J in J0 * j]
    return sorted(path for path, J in paths if SU2Rep(0,) in J * isospins[-1])


@lru_cache(maxsize=None)
def su2_singlet_coeff(isospins: Tuple[SU2Rep, ...],
                      path: Optional[Tuple[SU2Rep, ...]] = None) -> Dict[Tuple[U1Rep, ...], sp.Expr]:
    '''Rank-n SU(2) singlet CG tensor for the successive-coupling chain of the
    given isospins (order = slot order)::

        T[m1, ..., mn] = <j1 m1; ...; jn mn | 0 0>

    With a path (see su2_singlet_paths) only those intermediate spins are kept;
    without one, all paths are summed.

    Raw projection values (no field collapse / identical-multiplet rules; the
    Operator machinery applies those on its own). Nonzero only for m-tuples
    summing to zero. This is the numeric matrix declared in FeynRules via
    ``fr_definition``.
    '''
    j0 = isospins[0]
    state: Dict[Tuple[SU2Rep, U1Rep], Dict[Tuple[U1Rep, ...], sp.Expr]] = {
        (j0, m): {(m,): sp.Integer(1)} for m in j0.i3_range()}
    for step, j in enumerate(isospins[1:]):
        nxt: Dict[Tuple[SU2Rep, U1Rep], Dict[Tuple[U1Rep, ...], sp.Expr]] = {}
        for (J, M), monomials in state.items():
            block = cg_block(J, j)
            for m in j.i3_range():
                for Jout in J * j:
                    if path is not None and step < len(path) and Jout != path[step]:
                        continue
                    Mout = M + m
                    if Mout not in Jout.i3_range():
                        continue
                    cg = block.get((M, m, Jout, Mout), sp.Integer(0))
                    if cg == 0:
                        continue
                    bucket = nxt.setdefault((Jout, Mout), {})
                    for mt, c in monomials.items():
                        key = mt + (m,)
                        prev = bucket.get(key)
                        if prev is None:
                            bucket[key] = c * cg
                        else:
                            bucket[key] = sp.simplify(prev + c * cg)
        state = nxt
    out = state.get((SU2Rep(0,), U1Rep(0)), {})
    return {k: c for k, c in out.items() if sp.simplify(c) != 0}


# ----------------------------------------------------------------------
# FeynRules index types
# ----------------------------------------------------------------------

def feyn_index_type(rep) -> str:
    '''FeynRules index-type name for a rep: SU2Rep -> SU2D / SU2d{dim} (d <= 10);
    SU3Rep -> Colour / Gluon / Sextet / SU3d3bar / SU3d6bar (other reps: ValueError).'''
    if isinstance(rep, SU3Rep):
        if rep == SU3Rep(1, 1):
            return 'Gluon'
        if rep == SU3Rep(1, 0):
            return 'Colour'
        if rep == SU3Rep(2, 0):
            return 'Sextet'     # declared in SMMod.fr (representation T6)
        if (rep.upper, rep.lower) in ((0, 1), (0, 2)):
            return f'SU3d{rep.dim}bar'
        raise ValueError(f'No FeynRules index type for color {rep!r} '
                         f'(supported: 3, 3bar, 6, 6bar, 8)')
    if isinstance(rep, SU2Rep):
        return su2_index_type(rep.dim)
    raise TypeError(f'no FeynRules index type for {rep!r}')


# ----------------------------------------------------------------------
# CG tensor objects
# ----------------------------------------------------------------------

class CGTensor:
    '''A FeynRules-compatible CG object: knows only its own symbol 
    and the indices passed to it.'''

    def render(self, labels: Sequence[str]) -> str:
        return f'{self.name}[{", ".join(labels)}]'

    def fr_symbol(self) -> str:
        return self.name

    def fr_definition(self) -> str:
        return ''


class SU2CGTensor(CGTensor):
    '''SU(2) singlet CG tensor, ordered by isospin dimension ascending.

    The name follows the dim-sorted convention: isospins 4 5 2 2 2 3 give
    ``SU2CG222345``. Values are numeric (used for Python contractions) and
    declared in FeynRules through ``fr_definition``. Slot/index assignment
    to fields is done by the Operator/FeynRulesModel. Two-index tensors are the
    plain invariant (entries +-1, e.g. X^dagger X); larger ones are unit-normalized.
    '''

    def __init__(self, isospins: Sequence[SU2Rep], path: Optional[Tuple[SU2Rep, ...]] = None):
        kept = [j for j in isospins if j.dim > 1]   # isospin singlets are trivial
        self.isospins = tuple(sorted(kept, key=lambda j: j.dim))
        self.dims = tuple(j.dim for j in self.isospins)
        self.path = path
        self.name = 'SU2CG' + ''.join(map(str, self.dims)) if self.dims else ''
        if path is not None:   # one of several invariants: P<k>, k = position in su2_singlet_paths
            self.name += f'P{su2_singlet_paths(self.isospins).index(path) + 1}'
        self._numeric: Optional[Dict[Tuple[U1Rep, ...], sp.Expr]] = None

    def render(self, labels: Sequence[str]) -> str:
        if not self.dims:
            return ''
        return super().render(labels)

    @property
    def numeric(self) -> Dict[Tuple[U1Rep, ...], sp.Expr]:
        if self._numeric is None:
            self._numeric = (su2_singlet_coeff(self.isospins, self.path) if self.dims else {})
            if len(self.dims) == 2:
                norm = sp.sqrt(self.dims[0])
                self._numeric = {k: sp.simplify(norm * v) for k, v in self._numeric.items()}
        return self._numeric

    @property
    def index_types(self) -> Tuple[str, ...]:
        return tuple(feyn_index_type(SU2Rep(d - 1)) for d in self.dims)

    def fr_parameter(self):
        '''FeynRules Internal-parameter entry of the tensor (rigid CG
        coefficients, never set from outside): callable as
        ``SU2CG223[m1, m2, m3]`` with the SU(2) indices. Matrix-style
        Definitions: the non-zero entries explicitly, in index order, then a
        default of 0. The default only matches numeric indices (_?NumericQ),
        so the tensor stays symbolic (with its declared index types) until
        ExpandIndices has fixed the indices; zero entries then vanish from
        the Lagrangian before Feynman rules are computed.'''
        from uvbuilder.tools.fr_utils import FRParameter
        if not self.dims:
            return None
        labels = [f'i{k + 1}' for k in range(len(self.dims))]
        indices = ', '.join(f'Index[{t}]' for t in self.index_types)
        entries = []
        for mtuple, value in sorted(
                self.numeric.items(),
                key=lambda kv: tuple(su2_fr_index(d, m)
                                     for m, d in zip(kv[0], self.dims))):
            position = ', '.join(str(su2_fr_index(d, m))
                                 for m, d in zip(mtuple, self.dims))
            entries.append(f'{self.name}[{position}] -> {mathematica_code(value)}')
        entries.append(f'{self.name}[{", ".join(f"{l}_?NumericQ" for l in labels)}] :> 0')
        options = (('ParameterType', 'Internal'), ('Indices', f'{{{indices}}}'),
                   ('Definitions', '{' + ',\n '.join(entries) + '}'),
                   ('ComplexParameter', False))
        return FRParameter(self.name, 'CG', options, (self.name, 0))

    def fr_definition(self) -> str:
        '''Rendered entry of ``fr_parameter`` ('' for a trivial tensor).'''
        parameter = self.fr_parameter()
        return parameter.render() if parameter else ''


class SU3CGTensor(CGTensor):
    '''FeynRules SU(3) color contraction, resolved from the color reps.

    Mirrors ``SU2CGTensor``: takes the reps and derives everything itself,
    carrying no index names. Resolution is one of: a genuine head (K3, K3bar,
    K6, K6bar, T, T6, J, Jbar), a trivial delta/trace contraction (pure
    conjugate pairs, no symbol, one shared index per pair), or a notimplemented
    (color singlet but no known head / not pure pairs -> ``UserWarning``).

    ``external_slots`` are the field reps in head-argument order (render order
    and per-rep label iterators); ``dummy_count`` is the number of internal
    contracted indices of a composite head. ``render`` takes the external
    labels, optionally followed by the dummy labels (auto-filled from
    ``default_dummy`` when omitted); it raises ``ValueError`` when too few or
    an ambiguous number of labels is given.
    '''

    def __init__(self, reps: Sequence[SU3Rep]):
        self.reps = tuple(sorted(r for r in reps if r != SU3Rep(0, 0)))
        self.kind, self.data = _su3_case(self.reps)

    @property
    def is_head(self) -> bool:
        return self.kind == 'head'

    @property
    def is_trivial(self) -> bool:
        return self.kind == 'trivial'

    @property
    def is_notimplemented(self) -> bool:
        return self.kind == 'notimplemented'

    @property
    def name(self) -> str:
        return self.data.get('name', '') if self.is_head else ''

    @property
    def external_slots(self) -> Tuple[Tuple[SU3Rep, str], ...]:
        '''(rep, default label) per external field slot, in head-arg order.'''
        return self.data['external_slots'] if self.is_head else ()

    @property
    def dummy_count(self) -> int:
        return self.data['dummy_count'] if self.is_head else 0

    @property
    def default_dummy(self) -> Tuple[str, ...]:
        return self.data.get('default_dummy', ()) if self.is_head else ()

    @property
    def trivial_pairs(self) -> Tuple[Tuple[int, ...], ...]:
        '''Slot groups that share one label for a trivial contraction.'''
        return self.data.get('pairs', ()) if self.is_trivial else ()

    def render(self, labels: Sequence[str]) -> str:
        '''Render the head. ``labels`` holds the external labels, optionally
        followed by the dummy labels (auto-filled when omitted); too few or an
        ambiguous count raises ``ValueError``. Returns '' for trivial/notimplemented.'''
        if not self.is_head:
            return ''
        n, d = len(self.external_slots), self.dummy_count
        if len(labels) < n:
            raise ValueError(
                f'SU3CGTensor {self.name}: too few indices, expected at '
                f'least {n}, got {len(labels)}')
        if len(labels) not in (n, n + d):
            raise ValueError(
                f'SU3CGTensor {self.name}: expected {n} or {n + d} indices, '
                f'got {len(labels)}')
        external = labels[:n]
        dummy = labels[n:] if len(labels) == n + d else self.default_dummy[:d]
        data = self.data
        if data['factors'] is None:
            return f'{self.name}[{", ".join(external)}]'
        prefactor = data.get('prefactor', '')
        parts = [f'{symbol}[{", ".join(labels[p] if p < n else dummy[p - n] for p in pins)}]'
                 for symbol, pins in data['factors']]
        body = ' '.join(parts)
        return f'{prefactor} {body}'.strip() if prefactor else body

    def fr_definition(self) -> str:
        '''SU(3) heads are FeynRules built-ins: no definition.'''
        return ''


# ----------------------------------------------------------------------
# SU(3) contraction-head (FeynRules built-ins)
# ----------------------------------------------------------------------

SINGLET = SU3Rep(0,0)
TRIPLET = SU3Rep(1,0)
ANTITRIPLET = SU3Rep(0,1)
SEXTET = SU3Rep(2,0)
ANTISEXTET = SU3Rep(0,2)
OCTET = SU3Rep(1,1)

_HEAD_CASES = {
    (TRIPLET, TRIPLET, TRIPLET): dict(
        name='K3', external_slots=((TRIPLET, 'ii'), (TRIPLET, 'jj'), (TRIPLET, 'kk')),
        factors=None, dummy_count=0),
    (ANTITRIPLET, ANTITRIPLET, ANTITRIPLET): dict(
        name='K3bar', external_slots=((ANTITRIPLET, 'ii'), (ANTITRIPLET, 'jj'), (ANTITRIPLET, 'kk')),
        factors=None, dummy_count=0),
    # K6 multiplies an antisextet field, K6bar a sextet field (FeynRules UFO writer convention).
    (ANTITRIPLET, ANTITRIPLET, SEXTET): dict(
        name='K6bar', external_slots=((SEXTET, 'm'), (ANTITRIPLET, 'ii'), (ANTITRIPLET, 'jj')),
        factors=None, dummy_count=0),
    (ANTISEXTET, TRIPLET, TRIPLET): dict(
        name='K6', external_slots=((ANTISEXTET, 'm'), (TRIPLET, 'ii'), (TRIPLET, 'jj')),
        factors=None, dummy_count=0),
    (ANTITRIPLET, TRIPLET, OCTET): dict(
        name='T', external_slots=((OCTET, 'a'), (ANTITRIPLET, 'ii'), (TRIPLET, 'jj')),
        factors=None, dummy_count=0),
    (ANTISEXTET, OCTET, SEXTET): dict(
        name='T6', external_slots=((OCTET, 'a'), (ANTISEXTET, 'm'), (SEXTET, 'n')),
        factors=None, dummy_count=0),
    # J: K3[ii,xx,yy] T[a,xx,zz] K6bar[m,zz,yy]
    (TRIPLET, OCTET, SEXTET): dict(
        name='J', external_slots=((TRIPLET, 'ii'), (OCTET, 'a'), (SEXTET, 'm')),
        prefactor='-I Sqrt[2]',
        factors=(('K3', (0, 3, 4)), ('T', (1, 3, 5)), ('K6bar', (2, 5, 4))),
        dummy_count=3, default_dummy=('xx', 'yy', 'zz')),
    # Jbar: K6[m,yy,zz] T[a,zz,xx] K3bar[ii,xx,yy]
    (ANTITRIPLET, ANTISEXTET, OCTET): dict(
        name='Jbar', external_slots=((ANTISEXTET, 'm'), (OCTET, 'a'), (ANTITRIPLET, 'ii')),
        prefactor='I Sqrt[2]',
        factors=(('K6', (0, 3, 4)), ('T', (1, 4, 5)), ('K3bar', (2, 5, 3))),
        dummy_count=3, default_dummy=('yy', 'zz', 'xx')),
}


@lru_cache(maxsize=None)
def _su3_case(reps: Tuple[SU3Rep, ...]) -> Tuple[str, dict]:
    '''Resolve a sorted color-rep tuple into ('head'|'trivial'|'notimplemented', data).'''
    if not reps:
        return 'trivial', {'pairs': ()}
    if reps in _HEAD_CASES:
        return 'head', _HEAD_CASES[reps]
    counts = {r: reps.count(r) for r in set(reps)}
    if all(counts[r] == counts.get(r.conjRep, 0) for r in counts):
        pairs: List[Tuple[int, int]] = []
        used: Dict[int, bool] = {}
        for i, r in enumerate(reps):
            if used.get(i):
                continue
            j = next(k for k in range(i + 1, len(reps))
                     if not used.get(k) and reps[k] == r.conjRep)
            pairs.append((i, j))
            used[i] = used[j] = True
        return 'trivial', {'pairs': tuple(tuple(p) for p in pairs)}
    warnings.warn(
        f'SU3CGTensor: no color contraction known for {tuple(reps)}; '
        f'no term produced.')
    return 'notimplemented', {}
