"""Tests for the CG-tensor expansion in tools/operators.py.

Covers:
  * the cached CG block against hand-computed Clebsch-Gordan values,
  * the Operator.expand / expand_singlets contract-and-join behaviour,
  * equal-field monomial merging (and zero cancellation),
  * repeated identical multiplets (representative, no over-counting),
  * the 0-tolerant coupling (a vacuum-check that missing lines vanish).
"""

import re

import sympy as sp
import pytest

from uvbuilder.tools.reps import (BrokenRep, SU3Rep, SU2Rep, U1Rep, SCALAR,
                        WEYL_LEFT)
from uvbuilder.states.multiplet import Multiplet
from uvbuilder.states.SMfields import HiggsDoublet
from uvbuilder.tools.cg_tensors import cg_block
from uvbuilder.tools.operators import Operator, Term


def _doublet(label, hypercharge, i3_hi, i3_lo, y_hi, y_lo):
    return Multiplet(
        BrokenRep(spin=SCALAR, colordim=SU3Rep(0, 0),
                  isospin=SU2Rep(1,), hypercharge=hypercharge),
        label=label,
        labels={U1Rep(1, 2): i3_hi, U1Rep(-1, 2): i3_lo},
        fillMultiplet=True,
    )


# ----------------------------------------------------------------------
# CG block vs. hand values
# ----------------------------------------------------------------------

class TestCGBlock:
    def test_half_half_singlet(self):
        """1/2 x 1/2 -> 0 singlet must be the antisymmetric (+/-)1/sqrt(2) pair.
        <+,-|00> = +1/sqrt(2), <-,+|00> = -1/sqrt(2) (Condon-Shortley)."""
        block = cg_block(SU2Rep(1,), SU2Rep(1,))
        J0, M0 = SU2Rep(0,), U1Rep(0)
        assert block[(U1Rep(1, 2), U1Rep(-1, 2), J0, M0)] == sp.sqrt(2) / 2
        assert block[(U1Rep(-1, 2), U1Rep(1, 2), J0, M0)] == -sp.sqrt(2) / 2

    def test_half_half_triplet_edges(self):
        """1/2 x 1/2 -> J=1, M=+1 is the (++ ) component (coefficient 1)."""
        block = cg_block(SU2Rep(1,), SU2Rep(1,))
        J1, M1 = SU2Rep(2,), U1Rep(1)
        assert block[(U1Rep(1, 2), U1Rep(1, 2), J1, M1)] == sp.Integer(1)

    def test_forbidden_magnetic_line_absent(self):
        """M = m1 + m2 mismatch: no entry (never fabricated zero)."""
        block = cg_block(SU2Rep(1,), SU2Rep(1,))
        # m1 + M2 = 1, but J=0 requires M = 0 -> no such entry.
        assert (U1Rep(1, 2), U1Rep(1, 2), SU2Rep(0,), U1Rep(0)) not in block


# ----------------------------------------------------------------------
# Operator expansion
# ----------------------------------------------------------------------

Su = _doublet('Su', U1Rep(1, 2), 'Su+', 'Su0', None, None)
Sd = _doublet('Sd', U1Rep(-1, 2), 'Sd0', 'Sd-', None, None)


class TestExpand:
    def test_singlet_contains(self):
        op = Operator(Su, Sd, HiggsDoublet, HiggsDoublet.conjugateMultiplet)
        assert op.contains_singlet

    def test_singlets_from_su2_tensor(self):
        """The singlet monomials come from su2_tensor; repeated identical
        multiplets cancel to nothing (antisymmetric singlet of a symmetric
        bosonic operator)."""
        D = _doublet('D', U1Rep(0,), 'D+', 'D0', None, None)
        op = Operator(D, D)
        assert op.su2_tensor.name == 'SU2CG22'
        assert op.singlets == {}
        assert not op.contains_singlet

    def test_singlet_value_matches_hand(self):
        """Su0 Sd- Phi+ Phi0bar carries sqrt(3)/3 in the singlet (canonical fold)."""
        op = Operator(Su, Sd, HiggsDoublet, HiggsDoublet.conjugateMultiplet)
        terms = op.expand_singlets()
        by_label = {tuple(sorted(f.label for f in t.fields)): sp.simplify(t.coeff)
                    for t in terms}
        target = tuple(sorted(
            (Su.components[U1Rep(-1, 2)].label,
             Sd.components[U1Rep(-1, 2)].label,
             HiggsDoublet.components[U1Rep(1, 2)].label,
             HiggsDoublet.conjugateMultiplet.components[U1Rep(1, 2)].label)))
        coeff = by_label.get(target)
        assert coeff is not None
        assert sp.simplify(coeff - sp.sqrt(3) / 3) == 0

    def test_terms_are_canonical_terms(self):
        op = Operator(Su, Sd, HiggsDoublet, HiggsDoublet.conjugateMultiplet)
        for term in op.expand_singlets():
            assert isinstance(term, Term)
            assert isinstance(term.coeff, sp.Expr)

    def test_join_merges_and_cancels(self):
        """A repeated identical doublet couples antisymmetrically; for an
        identical commuting (bosonic) pair the antisymmetric singlet vanishes.
        The doublet must be real (Y=0) so that D (x) D conserves U(1)."""
        D = _doublet('D', U1Rep(0,), 'D+', 'D0', None, None)
        op = Operator(D, D)
        assert op.expand_singlets() == []


# ----------------------------------------------------------------------
# Merging helper (equal-field join)
# ----------------------------------------------------------------------

class TestJoining:
    # (test_strip_sign removed: conjugation phases are no longer carried as '-' labels,
    #  see Multiplet.phase and tests/test_conjugation.py)

    def test_parameter_applied(self):
        op = Operator(Su, Sd, HiggsDoublet, HiggsDoublet.conjugateMultiplet,
                      parameter=sp.Symbol('lambda'))
        for term in op.expand_singlets():
            assert 'lambda' in str(term.coeff)


# ----------------------------------------------------------------------
# Multiplet-level FeynRules term (Operator.fr_term)
# ----------------------------------------------------------------------

class TestFrTerm:
    def _multiplet(self, label, spin, col, iso, y, index=1):
        rep = BrokenRep(spin=spin, colordim=col, isospin=iso, hypercharge=y)
        labels = {i3: label for i3 in iso.i3_range(reverse=True)}
        return Multiplet(rep, label=label, labels=labels,
                         fillMultiplet=True, index=index)

    def test_su2_dim_sorted_indices(self):
        # Canonical order is A (dim2), C (dim2), B (dim3): dims sorted with
        # canonical tie-break give SU2CG223[m1, m2, m3] A[m1] C[m2] B[m3].
        A = self._multiplet('A', SCALAR, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 1)
        B = self._multiplet('B', SCALAR, SU3Rep(0, 0), SU2Rep(2,), U1Rep(0,), 2)
        C = self._multiplet('C', SCALAR, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 3)
        assert Operator(A, B, C).fr_term() == \
            '1 SU2CG223[m1, m2, m3] A[m1] C[m2] B[m3]'

    def test_colour_head_and_su2_none(self):
        T3 = self._multiplet('T3', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 1)
        A3 = self._multiplet('A3', SCALAR, SU3Rep(0, 1), SU2Rep(0,), U1Rep(0,), 2)
        O8 = self._multiplet('O8', SCALAR, SU3Rep(1, 1), SU2Rep(0,), U1Rep(0,), 3)
        assert Operator(T3, A3, O8, parameter=sp.Symbol('gs')).fr_term() == \
            'gs T[a, ii, jj] A3[ii] T3[jj] O8[a]'

    def test_weyl_pair_dot_product(self):
        chi = self._multiplet('chi', WEYL_LEFT, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 1)
        psi = self._multiplet('psibar', WEYL_LEFT, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 2)
        assert Operator(chi, psi).fr_term() == \
            '1 SU2CG22[m1, m2] (chiLw[sp, m1].psibarLw[sp, m2])'

    def test_pure_singlets_no_cg(self):
        P1 = self._multiplet('P1', SCALAR, SU3Rep(0, 0), SU2Rep(0,), U1Rep(0,), 1)
        P2 = self._multiplet('P2', SCALAR, SU3Rep(0, 0), SU2Rep(0,), U1Rep(0,), 2)
        assert Operator(P1, P2).fr_term() == '1 P1 P2'

    def test_no_singlet_warns_and_returns_empty(self):
        # 2 (x) 2 (x) 2 has no SU(2) singlet: fr_term() warns and yields ''. 
        D1 = self._multiplet('D1', SCALAR, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 1)
        D2 = self._multiplet('D2', SCALAR, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 2)
        D3 = self._multiplet('D3', SCALAR, SU3Rep(0, 0), SU2Rep(1,), U1Rep(0,), 3)
        with pytest.warns(UserWarning, match='No Singlet'):
            assert Operator(D1, D2, D3).fr_term() == ''

    def test_repeated_multiplet_distinct_indices(self):
        # Two identical triplets + Higgs pair: each occurrence gets its own index.
        # BSM (SI3Y0) sorts before SM (H, Hbar) in canonical order; the two SI3Y0
        # occurrences take the dim-3 slots (m3, m4), the Higgs pair m1, m2.
        T = self._multiplet('SI3Y0', SCALAR, SU3Rep(0, 0), SU2Rep(2,), U1Rep(0,), 1)
        op = Operator(T, T, HiggsDoublet, HiggsDoublet.conjugateMultiplet)
        assert op.fr_term() == \
            '1 SU2CG2233[m1, m2, m3, m4] SI3Y0[m3] SI3Y0[m4] Phitilde[m1] Phi[m2]'

    def test_repeated_colored_multiplet_distinct_indices(self):
        # Two identical colour triplets under K3: per-occurrence colour labels.
        P = self._multiplet('P', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 1)
        R = self._multiplet('R', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 2)
        op = Operator(P, P, R)
        assert op.fr_term() == '1 K3[ii, jj, kk] P[ii] P[jj] R[kk]'

    def test_tracker_no_head_pair_shares_index(self):
        # No SU(3) head (e.g. 8 x 8): a conjugate pair contracts through the
        # identity/trace, so one fresh index is shared by both fields.
        from uvbuilder.models.feynrulesmodel import IndexTracker
        O = self._multiplet('O', WEYL_LEFT, SU3Rep(1, 1), SU2Rep(0,), U1Rep(0,), 1)
        term = Operator(O, O).fr_term(IndexTracker())
        assert term == '1 (OLw[sp1, c1].OLw[sp1, c1])'

    def test_tracker_antitriplet_pair_shares_index(self):
        # 3 x 3* without a head contracts via delta: one shared index.
        from uvbuilder.models.feynrulesmodel import IndexTracker
        T = self._multiplet('T', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 1)
        Tb = self._multiplet('Tb', SCALAR, SU3Rep(0, 1), SU2Rep(0,), U1Rep(0,), 2)
        term = Operator(T, Tb).fr_term(IndexTracker())
        assert term == '1 Tb[c1] T[c1]'

    def test_tracker_head_gets_distinct_labels_per_slot(self):
        # A real SU(3) head (K3) still gets one distinct fresh label per slot.
        from uvbuilder.models.feynrulesmodel import IndexTracker
        P = self._multiplet('P', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 1)
        R = self._multiplet('R', SCALAR, SU3Rep(1, 0), SU2Rep(0,), U1Rep(0,), 2)
        op = Operator(P, P, R)
        assert op.fr_term(IndexTracker()) == '1 K3[c1, c2, c3] P[c1] P[c2] R[c3]'

    def test_tracker_su2_labels_follow_tensor_slots(self):
        # SU2CG slots are ordered by isospin dim (2, 2, 3); with a tracker the
        # triplet must carry the label of the dim-3 slot, whatever its position.
        from uvbuilder.models.feynrulesmodel import IndexTracker
        T = self._multiplet('T', SCALAR, SU3Rep(0, 0), SU2Rep(2,), U1Rep(0,))
        term = Operator(T, HiggsDoublet, HiggsDoublet.conjugateMultiplet).fr_term(IndexTracker())
        (slots,) = re.findall(r'SU2CG223\[(\w+), (\w+), (\w+)\]', term)
        assert re.search(rf'\bT\[{slots[2]}\]', term), term
        assert not re.search(rf'\bT\[({slots[0]}|{slots[1]})\]', term), term
