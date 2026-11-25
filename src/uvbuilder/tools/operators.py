"""
.. module:: operators
   :synopsis: Defines the Operator class, a coupling of multiplets contracted to an SU(2) singlet.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""

import sympy as sp
from sympy.printing.mathematica import mathematica_code
from typing import NamedTuple, Tuple, Dict, List, Optional, Union
import copy
import dataclasses
import warnings

from uvbuilder.states.multiplet import Multiplet
from uvbuilder.states.field import Field
from uvbuilder.states.SMfields import HiggsDoublet
from uvbuilder.tools.reps import (U1Rep, BrokenRep,
                        WeylSpinor, DiracSpinor, Scalar, WEYL_RIGHT)
from uvbuilder.tools.cg_tensors import SU2CGTensor, SU3CGTensor
from uvbuilder.tools.fr_utils import DAGGER_SUFFIX, HIGGS_TILDE


# ----------------------------------------------------------------------
# Canonical ordering
# ----------------------------------------------------------------------

def canonical_sort(multiplet_or_field: Multiplet | Field):
    '''
    Sorting function to canonically sort multiplets/fields inside an operator.
    Follows BSM > spin > globalU1s > isospin > hypercharge > family_index > eCharge (fields).
    '''
    rep = multiplet_or_field.rep

    multiplet_order = (
        int(multiplet_or_field.isSM),
        rep.spin,
        rep.colordim,
        rep.globalU1s,
        rep.isospin,
        rep.hypercharge,
        multiplet_or_field.index)

    if rep.is_multiplet_rep:
        return multiplet_order
    field_order = multiplet_order + (rep.eCharge,)
    return field_order


# ----------------------------------------------------------------------
# Records returned by Operator.expand / expand_singlets
# ----------------------------------------------------------------------

class Term(NamedTuple):
    '''One term of an operator expansion (an SU(2) singlet monomial).

    coeff            - sympy coefficient (already includes the operator parameter)
    fields           - Field objects, canonically ordered (unique per monomial)

    Color contraction is resolved only at the multiplet level (Operator.fr_term);
    the algebraic expansion carries no SU(3) data.
    '''
    coeff: sp.Expr
    fields: Tuple[Field]

    def to_sympy_expr(self) -> sp.Expr:
        return self.coeff * sp.prod([sp.Symbol(f.label) for f in self.fields])

    @property
    def dim(self) -> int:
        vevs = int(self.coeff.as_coeff_exponent(sp.Symbol('v'))[1])
        return sum(f.spin.energy_dim for f in self.fields) + vevs

    def at_vev(self) -> Optional['Term']:
        '''Replaces the neutral components of the Higgs doublet(s) with v/sqrt(2)
        and the charged components with 0. Returns a new Term, or None if the
        monomial vanishes after the vev substitution.'''
        Phi0 = HiggsDoublet.components[U1Rep(-1, 2)]
        Phi0bar = HiggsDoublet.conjugateMultiplet.components[U1Rep(1, 2)]
        PhiP = HiggsDoublet.components[U1Rep(1, 2)]
        PhiM = HiggsDoublet.conjugateMultiplet.components[U1Rep(-1, 2)]

        if PhiP in self.fields or PhiM in self.fields:
            return None

        counts = sum(1 for f in self.fields if f in (Phi0, Phi0bar))
        vev = (sp.Symbol('v') / sp.sqrt(2)) ** counts

        fields = tuple(f for f in self.fields if f not in (Phi0, Phi0bar, PhiP, PhiM))
        return Term(coeff=vev * self.coeff, fields=fields)


# ----------------------------------------------------------------------
# FeynRules / Mathematica emission helpers
# ----------------------------------------------------------------------
#
# Sympy -> Mathematica conversion is delegated to sympy's
# ``mathematica_code`` (which also parenthesizes additive sub-factors
# correctly, so no extra Mul handling is needed here).


# ----------------------------------------------------------------------
# The main class
# ----------------------------------------------------------------------

class Operator:
    '''A single coupling operator: a collection of canonically ordered multiplets
    contracted to SU(2) singlets via their Clebsch-Gordan coefficients.

    The only class in this module carrying behavior; output is plain Term records.
    '''

    def __init__(self, *multiplets: Multiplet, parameter=1, path=None):
        self.multiplets = tuple(sorted(multiplets, key=canonical_sort))
        self.parameter = sp.sympify(parameter)
        self.path = path   # SU(2) coupling path, when the fields form several invariants
        self._singlets: Optional[Dict[Tuple[Field, ...], sp.Expr]] = None

    @property
    def dim(self) -> int:
        return sum(m.spin.energy_dim for m in self.multiplets)

    @property
    def conserved(self) -> frozenset:
        '''Global U(1)s tracked in this model: the ones carried by the BSM
        multiplets (all BSM multiplets of a model track the same charges).'''
        return frozenset().union(*(m.rep.globalU1s.tracked
                                   for m in self.multiplets if not m.isSM))

    def _projected_rep(self, m: Multiplet) -> BrokenRep:
        '''SM multiplets carry all their global charges; keep only the conserved ones.'''
        if not m.isSM:
            return m.rep
        return dataclasses.replace(m.rep, globalU1s=m.rep.globalU1s.project(self.conserved))

    @property
    def possible_reps(self) -> List[BrokenRep]:
        reps = [self._projected_rep(m) for m in self.multiplets]
        result = reps[0]
        for rep in reps[1:]:
            result = result * rep
        if isinstance(result, BrokenRep):
            result = [result]
        return list(result)

    @property
    def contains_singlet(self) -> bool:
        '''True iff the couplings can be contracted to the (0, 0) singlet with a
        surviving (non-canceling) monomial.'''
        any_singlet_rep = any(rep.is_singlet for rep in self.possible_reps)
        return any_singlet_rep and bool(self.singlets)

    @property
    def su2_tensor(self) -> SU2CGTensor:
        '''The SU(2) singlet CG tensor of this operator (dim-1 isospins dropped).'''
        return SU2CGTensor([m.isospin for m in self.multiplets], path=self.path)

    @property
    def su3_tensor(self) -> SU3CGTensor:
        '''The SU(3) color contraction of this operator (color singlets dropped).'''
        return SU3CGTensor([m.colordim for m in self.multiplets])

    @property
    def singlets(self) -> Dict[Tuple[Field, ...], sp.Expr]:
        '''Merged singlet monomials: {
        canonical-field-tuple: coeff}. Empty dict => no surviving singlet.

        Maps the SU(2) tensor's numeric values onto the component fields
        (tensor slot order = isospin dim ascending, ties by canonical order),
        times the phase with which each field enters its multiplet
        (Multiplet.phase: the SU(2) conjugation phases of conjugate multiplets),
        canonicalizes the field tuples and merges coefficients, so repeated
        identical multiplets add (or cancel) as in the successive-coupling
        contraction.'''
        if self._singlets is not None:
            return self._singlets
        if self.is_scalar_mass_term:
            self._singlets = self._mass_singlets()
            return self._singlets
        tensor = self.su2_tensor
        non_singlet = [m for m in self.multiplets if m.isospin.dim > 1]
        slots = sorted(range(len(non_singlet)),
                       key=lambda j: non_singlet[j].isospin.dim)
        slot_of = {j: k for k, j in enumerate(slots)}
        singletons = [m.components[U1Rep(0)]
                      for m in self.multiplets if m.isospin.dim == 1]
        singleton_phase = sp.Integer(1)
        for m in self.multiplets:
            if m.isospin.dim == 1:
                singleton_phase *= m.phase(U1Rep(0))

        merged: Dict[Tuple[Field, ...], sp.Expr] = {}
        if not tensor.dims:
            merged[tuple(sorted(singletons, key=canonical_sort))] = singleton_phase
        else:
            for i3_tuple, coeff in tensor.numeric.items():
                i3s = [i3_tuple[slot_of[j]] for j in range(len(non_singlet))]
                fields = [non_singlet[j].components[i3] for j, i3 in enumerate(i3s)]
                phase = singleton_phase
                for j, i3 in enumerate(i3s):
                    phase *= non_singlet[j].phase(i3)
                key = tuple(sorted(fields + singletons, key=canonical_sort))
                merged[key] = sp.simplify(merged.get(key, 0) + phase * coeff)
        merged = {k: c for k, c in merged.items() if c != 0}
        self._flipped = False
        if self.is_fermion_mass_term:
            # the convention of the mass terms (Multiplet.mass_partner); the .fr writes the
            # CG singlet, which is the same up to its sign
            partner = self._partner_singlets()
            ratios = {sp.simplify(partner.get(k, 0) / merged[k]) for k in merged}
            assert partner.keys() == merged.keys() and len(ratios) == 1 and abs(*ratios) == 1, \
                f'{self.pretty_print()}: mass term is not the CG singlet'
            self._flipped, merged = ratios == {-1}, partner
        self._singlets = merged
        return self._singlets

    def _partner_singlets(self) -> Dict[Tuple[Field, ...], sp.Expr]:
        '''sum_i3 phase self[i3] partner of a fermion mass term (Multiplet.mass_partner).'''
        first, second = self.multiplets
        merged: Dict[Tuple[Field, ...], sp.Expr] = {}
        for i3, f in first.components.items():
            g, phase = first.mass_partner(second, i3)
            phase *= first.phase(i3) * second.phase(-i3)
            key = tuple(sorted((f, g), key=canonical_sort))
            merged[key] = merged.get(key, 0) + phase
        return {k: c for k, c in merged.items() if c != 0}

    @property
    def is_fermion_mass_term(self) -> bool:
        '''Majorana or Dirac mass of BSM Weyl multiplets.'''
        return (len(self.multiplets) == 2 and not any(m.isSM for m in self.multiplets)
                and all(isinstance(m.rep.spin, WeylSpinor) for m in self.multiplets))

    def _mass_singlets(self) -> Dict[Tuple[Field, ...], sp.Expr]:
        '''sum_i X^*_i X_i of a scalar mass term, the term fr_mass_term writes (the CG
        singlet of a real multiplet of isospin dim > 1 has the opposite sign).'''
        declared = max(self.multiplets, key=lambda m: m.rep)
        merged: Dict[Tuple[Field, ...], sp.Expr] = {}
        for i3, f in declared.components.items():
            partner, phase = declared.conjugate_component(i3)
            key = tuple(sorted((partner, f), key=canonical_sort))
            merged[key] = merged.get(key, 0) + phase
        return {k: c for k, c in merged.items() if c != 0}

    @property
    def forms_bilinear(self) -> bool:
        multiplets = [
            m for m in self.multiplets
            if m != HiggsDoublet and m != HiggsDoublet.conjugateMultiplet]
        if len(multiplets) != 2:
            return False
        m1, m2 = multiplets
        family_off_diagonal = (
            m1.rep == m2.rep.conjInternalRep and 
            m1.index != m2.index and
            len(self.multiplets) == 2
        )
        not_empty = self.expand_singlets(at_vev=True) != []
        return not family_off_diagonal and not_empty

    @property
    def is_scalar_mass_term(self) -> bool:
        '''Bilinear of a scalar multiplet with its conjugate (or a real one with itself).
        A self-conjugate multiplet of half-integer isospin is pseudo-real: no mass term.'''
        if len(self.multiplets) != 2:
            return False
        m1, m2 = self.multiplets
        pseudo_real = m1.selfConjugate and m1.isospin.dim % 2 == 0
        return (isinstance(m1.rep.spin, Scalar) and not m1.isSM and not pseudo_real
                and m2.rep == m1.rep.conjRep)

    @property
    def is_hermitian(self) -> bool:
        '''O^dagger = O: conjugating every field maps the singlet expansion onto itself.'''
        if any(isinstance(m.rep.spin, (WeylSpinor, DiracSpinor)) for m in self.multiplets):
            return False
        conj = self._conjugated(self.singlets)
        return (conj.keys() == self.singlets.keys()
                and all(sp.simplify(conj[k] - self.singlets[k]) == 0 for k in conj))

    def _conjugated(self, monomials) -> Dict[Tuple[Field, ...], sp.Expr]:
        '''{fields: coeff} of the hermitian conjugate of a sum of monomials of the fields
        of this operator or of its conjugate: each field replaced by its conjugate
        (Multiplet.conjugate_component), coefficients conjugated with v real.'''
        conj_of = {f: mm.conjugate_component(i3)
                   for m in self.multiplets for mm in (m, m.conjugateMultiplet)
                   for i3, f in mm.components.items()}
        v = sp.Symbol('v')
        out: Dict[Tuple[Field, ...], sp.Expr] = {}
        for fields, c in monomials.items():
            phase = sp.prod([conj_of[f][1] for f in fields])
            key = tuple(sorted((conj_of[f][0] for f in fields), key=canonical_sort))
            out[key] = out.get(key, 0) + sp.conjugate(c).xreplace({sp.conjugate(v): v}) * phase
        return out

    def hermitian_conjugate(self, at_vev: bool = False) -> List[Term]:
        '''The terms of O^dagger, as expand_singlets gives those of O.'''
        terms = {t.fields: t.coeff for t in self.expand_singlets(at_vev)}
        return [Term(coeff=c, fields=f) for f, c in self._conjugated(terms).items() if c != 0]

    def lagrangian_terms(self, at_vev: bool = False) -> List[Term]:
        '''The terms of -L from this operator, as to_feynrules writes them: O, plus
        O^dagger unless O is hermitian.'''
        terms = self.expand_singlets(at_vev)
        return terms if self.is_hermitian else terms + self.hermitian_conjugate(at_vev)

    # -- SU(2) expansion --------------------------------------------------

    def expand_singlets(self, at_vev: bool = False) -> List[Term]:
        '''The (0, 0) singlet monomials as a flat list of canonical Terms.'''
        if not self.contains_singlet:
            return []
        terms = [Term(coeff=sp.simplify(self.parameter * coeff), fields=fields)
                 for fields, coeff in self.singlets.items()]
        if not at_vev:
            return terms
        out = []
        for term in terms:
            t = term.at_vev()
            if t is not None and t.coeff != 0:
                out.append(t)
        return out

    # -- FeynRules / Mathematica (multiplet level) ---------------

    def output_parameter(self, tracker=None) -> sp.Expr:
        '''The coupling with the output names of its symbols (parameter renames
        decided by FeynRulesModel travel on the tracker).'''
        names = getattr(tracker, 'parameter_names', None)
        if not names:
            return self.parameter
        return self.parameter.subs({sp.Symbol(old): sp.Symbol(new) for old, new in names.items()
                                    if sp.Symbol(old) in self.parameter.free_symbols})

    def fr_term(self, tracker=None) -> str:
        '''Multiplet-level FeynRules term::

            <param> <SU3CG[...]> <SU2CG[...]> <multiplet symbols with indices>

        The SU(3) contraction head (if any) is emitted first, then the SU(2)
        singlet CG tensor (if any non-singlet isospin), then the multiplets in
        canonical order. Each multiplet symbol carries ``[color, isospin]``
        index labels. The SU(2) tensor is ordered by isospin dimension
        ascending (ties by canonical field order), so e.g. canonical
        A(2), B(3), C(2) gives ``SU2CG223[m1, m2, m3] A[m1] B[m3] C[m2]``.
        Index labels use the CG-tensor defaults unless a ``tracker`` is given,
        in which case fresh unique symbols are drawn from
        ``tracker.new_symbol(style)`` so the caller (FeynRulesModel) owns
        index uniqueness. Adjacent spinor multiplets are rendered as Weyl dot
        products with a shared spin index, e.g.
        ``(DdLw[sp, i1].DuLw[sp, i2])``; left+right pairs conjugate the
        right-handed Weyl. The returned term is the bare product (the caller
        adds the hermitian conjugate, e.g. ``- term - HC[term]``).'''

        if not self.contains_singlet:
            warnings.warn(
                f'{self.pretty_print()}: No Singlet. No Lagrangian Term can be produced.')
            return ''

        if tracker is not None:
            tracker.reset_symbols()

        su3 = self.su3_tensor
        if su3.is_notimplemented:
            warnings.warn(
                f'{self.pretty_print()}: Unresolved color contraction; '
                f'no term produced.')
            return ''

        non_singlet = [m for m in self.multiplets if m.isospin.dim > 1]
        su2 = self.su2_tensor

        if tracker is not None:
            # Unique labels from the index tracker. Every color index,
            # including the dummy indices of a composite head, comes from the
            # tracker's 'c' pool (externals first, then dummies).
            # Field j gets iso_labels[j]; tensor slot k belongs to the field slots[k]
            # (slots sorted by isospin dim, canonical ties), as in singlets().
            iso_labels = [tracker.new_symbol('i') for _ in non_singlet]
            iso = iter(iso_labels)
            slots = sorted(range(len(non_singlet)),
                           key=lambda j: non_singlet[j].isospin.dim)
            su2_labels = [iso_labels[j] for j in slots]
            if su3.is_head:
                fresh = [tracker.new_symbol('c') for _ in su3.external_slots]
                dummy = [tracker.new_symbol('c') for _ in range(su3.dummy_count)]
            elif su3.is_trivial:
                fresh = [tracker.new_symbol('c') for _ in su3.trivial_pairs]
                dummy = []
        else:
            # Deterministic names for tests/debug: default isospin labels and
            # the head's per-arg default labels; dummies auto-fill in render.
            slots = sorted(range(len(non_singlet)),
                           key=lambda i: non_singlet[i].isospin.dim)   # dim asc, canonical ties
            iso = iter(f'm{k + 1}' for k in slots)
            su2_labels = [f'm{k + 1}' for k in range(len(su2.dims))]
            if su3.is_head:
                fresh = [lab for _, lab in su3.external_slots]
                dummy = []
            elif su3.is_trivial:
                fresh = [f'c{k + 1}' for k in range(len(su3.trivial_pairs))]
                dummy = []

        su3_render = ''
        color: Dict = {}
        if su3.is_head:
            su3_render = su3.render(fresh + dummy)
            for (rep, _), lab in zip(su3.external_slots, fresh):
                color.setdefault(rep, []).append(lab)
            color = {rep: iter(labs) for rep, labs in color.items()}
        elif su3.is_trivial:
            # Pure conjugate pairs contract through the identity/trace: each
            # pair shares one index (delta_{i j*} / adjoint trace a == a).
            for (i, j), lab in zip(su3.trivial_pairs, fresh):
                color.setdefault(su3.reps[i], []).append(lab)
                color.setdefault(su3.reps[j], []).append(lab)
            color = {rep: iter(labs) for rep, labs in color.items()}

        parts = [mathematica_code(self.output_parameter(tracker))]
        if su3_render:
            parts.append(su3_render)
        if su2.dims:
            parts.append(su2.render(su2_labels))

        items: List[str] = []
        spinor_run: List[Multiplet] = []

        def emit(m: Multiplet, sp: Union[str, None] = None,
                 bar: bool = False) -> str:
            iso_idx = next(iso, None) if m.isospin.dim > 1 else None
            idx = [next(color[m.colordim], None) if m.colordim in color else None,
                   iso_idx]
            symbol = m.fr_symbol() + (DAGGER_SUFFIX if bar else '')   # FeynRules conjugate: always append
            # SM multiplets are written as the classes of the SM model file the
            # tracker knows about (e.g. QLw[sp, i, 3, c], uRwbar[sp, 3, c]).
            sm_map = getattr(tracker, 'sm', None)
            if m.isSM and sm_map is not None:
                return sm_map.render(m, spin=sp, color=idx[0], isospin=iso_idx, bar=bar)
            # A conjugate-partner multiplet enters as its tilde class psi^c = C psi-bar
            tilde = (getattr(tracker, 'conjugate_symbols', None) or {}).get(m)
            if tilde is not None:
                symbol = tilde + (DAGGER_SUFFIX if bar else '')
            # FeynRules translation: our conjugate Higgs is H~ = i sigma2 H* (components
            # ordered by i3), declared as the unphysical class HIGGS_TILDE whose
            # Definitions map it onto FeynRules' Phibar (= H*); see
            # multipletModel.tilde_class. Keeps every index on a field.
            if m is HiggsDoublet.conjugateMultiplet:
                return f'{HIGGS_TILDE}[{iso_idx}]'
            if sp is not None:
                idx = [sp] + idx
            if any(idx):
                return f'{symbol}[{", ".join(filter(None, idx))}]'
            return symbol

        def new_sp() -> str:
            return tracker.new_symbol('sp') if tracker is not None else 'sp'

        def spinor_product(pair: List[Multiplet]) -> str:
            '''Weyl dot product of a spinor pair with a shared spin index,
            e.g. ``(DdLw[sp, i1].DuLw[sp, i2])``. For a left+right pair the
            right-handed Weyl is conjugated (``psiL[sp].psiRbar[sp]``).'''
            a, b = pair
            sp = new_sp()
            a_rh = a.rep.spin == WEYL_RIGHT
            b_rh = b.rep.spin == WEYL_RIGHT
            if a_rh == b_rh:
                return f'({emit(a, sp)}.{emit(b, sp)})'
            if b_rh:
                return f'({emit(a, sp)}.{emit(b, sp, bar=True)})'
            return f'({emit(a, sp, bar=True)}.{emit(b, sp)})'

        for m in self.multiplets:
            if isinstance(m.rep.spin, (WeylSpinor, DiracSpinor)):
                spinor_run.append(m)
                continue
            if spinor_run:
                items.append(' '.join(
                    spinor_product(spinor_run[i:i + 2])
                    for i in range(0, len(spinor_run), 2)))
                spinor_run.clear()
            items.append(emit(m))
        if spinor_run:
            items.append(' '.join(
                spinor_product(spinor_run[i:i + 2])
                for i in range(0, len(spinor_run), 2)))

        parts.extend(items)
        return ' '.join(parts)

    def fr_mass_term(self, tracker) -> str:
        '''Scalar mass term without CG tensor: M^2 Xbar[i] X[i] (complex),
        1/2 M^2 Xbar[i] X[i] (real, isospin dim > 1), 1/2 M^2 X[i] X[i] (real).'''
        tracker.reset_symbols()
        declared = max(self.multiplets, key=lambda m: m.rep)
        symbol = declared.fr_symbol()
        labels = [tracker.new_symbol(index.style) for index in tracker.gauge(declared.rep)]
        idx = f'[{", ".join(labels)}]' if labels else ''
        real = declared.selfConjugate and not declared.is_real_nontrivial
        left = symbol if real else symbol + DAGGER_SUFFIX
        return f'{mathematica_code(self.output_parameter(tracker))} {left}{idx} {symbol}{idx}'

    def fr_definition(self, tracker=None) -> str:
        '''FeynRules interaction definition for this operator (the term
        string; the model supplies the index tracker for unique labels).'''
        return self.fr_term(tracker)

    @property
    def sector(self) -> str:
        '''Mass (two fields), Trilinear (three scalars), Yukawa (two fermions and
        a scalar), Quartic (four scalars) or Other.'''
        n = len(self.multiplets)
        n_fermions = sum(not isinstance(m.rep.spin, Scalar) for m in self.multiplets)
        return {(2, 0): 'Mass', (2, 2): 'Mass', (3, 0): 'Trilinear',
                (3, 2): 'Yukawa', (4, 0): 'Quartic'}.get((n, n_fermions), 'Other')

    def _description(self) -> str:
        labels = ' '.join(m.label for m in self.multiplets)
        if self.sector != 'Mass':
            # label only when the field set keeps several invariants (parameter suffix P<k>)
            name = self.su2_tensor.name
            suffix = name[name.rindex('P'):] if self.path is not None else None
            several = suffix and any(str(s).endswith(suffix) for s in self.parameter.free_symbols)
            invariant = f' (SU(2) invariant {suffix})' if several else ''
            return f'{self.sector} coupling {labels}{invariant}'
        if self.is_scalar_mass_term:
            return f'Mass of {max(self.multiplets, key=lambda m: m.rep).label}'
        first, second = self.multiplets
        return (f'Majorana mass of {first.label}' if first is second
                else f'Dirac mass of {first.label} {second.label}')

    def to_feynrules(self, tracker) -> Dict[str, list]:
        '''FeynRules chunk for this operator: the Lagrangian term(s) with its
        hermitian conjugate, the coupling parameter(s) and the CG tensor.'''

        from uvbuilder.tools.fr_utils import parameter_block, SECTOR_BLOCK

        if self.is_scalar_mass_term:
            term = self.fr_mass_term(tracker)
            interaction = f'- {term}'
        elif self.singlets and self._flipped:
            # the CG singlet with the sign of the mass convention (Multiplet.mass_partner)
            term = self.fr_term(tracker)
            interaction = f'+ {term} + HC[{term}]'
        else:
            term = self.fr_term(tracker)
            interaction = f'- {term}' if self.is_hermitian else f'- {term} - HC[{term}]'
        sector, block = self.sector, SECTOR_BLOCK[self.sector]
        chunk: Dict[str, list] = {f'{sector}Terms': [interaction] if term else [],
                                  'Parameters': []}
        # SM convention: QED order = (number of fields) - 2, i.e. Yukawas and
        # scalar cubics {QED, 1}, quartics {QED, 2}; bilinears get none (they
        # only enter mass terms, never a vertex). MG5 rejects order-less couplings.
        n_fields = len(self.multiplets)
        order = f'{{QED, {n_fields - 2}}}' if n_fields > 2 else None
        # Hermitian operators and fermion masses have real couplings; the others are
        # complex, written as real External parts and an Internal complex parameter
        fermion_mass = n_fields == 2 and not isinstance(self.multiplets[0].rep.spin, Scalar)
        real = self.is_scalar_mass_term or fermion_mass or self.is_hermitian
        description = f'"{self._description()}"'
        for symbol in self.output_parameter(tracker).free_symbols:
            name = str(symbol)
            if real:
                chunk['Parameters'].append(parameter_block(
                    name, sector, value=1, description=description, complex_parameter=False,
                    interaction_order=order, block_name=block))
                continue
            for part, value, label in (('R', 1, 'real'), ('I', 0, 'imaginary')):
                chunk['Parameters'].append(parameter_block(
                    name, sector, value=value, complex_parameter=False,
                    description=f'"{self._description()} ({label} part)"',
                    block_name=block if part == 'R' else f'IM{block}', part=part))
            chunk['Parameters'].append(parameter_block(
                name, sector, parameter_type='Internal', value=f'R{name} + I*I{name}',
                description=description, complex_parameter=True, interaction_order=order))
        tensor = self.su2_tensor
        if tensor.dims and not self.is_scalar_mass_term:
            # Rigid CG coefficients: Internal parameter
            chunk['Parameters'].append(tensor.fr_parameter())
        return chunk

    # -- printing ----------------------------------------------------------

    def pretty_print(self) -> str:
        return f'{self.parameter} ' + ' '.join(m.label for m in self.multiplets)

    def __repr__(self):
        return self.pretty_print()

    def __eq__(self, other: "Operator") -> bool:
        return self.multiplets == other.multiplets and self.parameter == other.parameter
