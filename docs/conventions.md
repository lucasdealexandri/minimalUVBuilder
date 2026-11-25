# uvbuilder conventions

The conventions of the models `uvbuilder build` writes and of the masses and mixing
matrices `uvbuilder fill` computes. Each model report (`<model>.md`) lists the model's
states, input parameters and mass matrices; this file explains them.

## Mass terms

- A positive mass parameter gives positive masses. Each mass term pairs every component
  with its partner in the SU(2) invariant (`Multiplet.mass_partner`): a Majorana mass
  gives $+M$ to the neutral component, a Dirac mass $+M$ to its most charged pair. The
  other pairs keep the relative sign SU(2) fixes, e.g. $-M$ for the neutral pair of two
  doublets, $-M$ for the charged pair of a Majorana triplet.
- These signs, and those of negative inputs, are absorbed into the mass eigenstates: the
  phase $P = \mathrm{sign}(M)$ (Dirac, on the right-handed partner) or
  $P = \sqrt{\mathrm{sign}(M)}$ (Majorana) of a state alone in its block, and the phases
  of the mixing matrices otherwise.

## Mass matrices

| block | matrix | $-\mathcal{L}$ | diagonalisation |
|---|---|---|---|
| Majorana fermion | complex symmetric | $\frac12 \psi^T M \psi + \mathrm{h.c.}$ | Takagi: $U^T M U = \mathrm{diag}(m)$, $U$ unitary |
| Dirac fermion | complex | $\psi^T M \chi + \mathrm{h.c.}$ ($\psi$: rows, $\chi$: columns) | SVD: $U^T M V = \mathrm{diag}(m)$, $U$, $V$ unitary |
| real scalar | real symmetric | $\frac12 \rho^T M^2 \rho$ | $O^T M^2 O = \mathrm{diag}(m^2)$, $O$ real orthogonal |
| complex scalar | Hermitian | $\phi^\dagger M^2 \phi$ (rows: $\phi^*$, columns: $\phi$) | $U^\dagger M^2 U = \mathrm{diag}(m^2)$, $U$ unitary |

- Gauge eigenstates $= U \times$ mass eigenstates; columns in ascending mass (the $j$-th
  lightest state is the family index $j$). Each column is rephased so that its entry of
  largest modulus is real and positive (real scalars and Majorana fermions: only a sign
  is free). Masses must be positive and non-degenerate within a block.
- A state alone in its block has mass $|M|$ (fermions) or $\sqrt{M^2}$ (scalars); for
  fermions the sign of $M$ goes into the phase $P$ above. Its param_card entry is only
  its MASS.
- The param_card block of a mixing matrix is its name in capitals (`UMIX...`), with the
  imaginary parts in `IMUMIX...`; Dirac blocks have two matrices, $U$ and $V$.
- A neutral complex scalar $\Phi$ mixing with real states is written in the real basis
  $\Phi = (\Phi_R + i\,\Phi_I)/\sqrt2$: its block has the rows $\Phi_R$ and $\Phi_I$.

## Parameters

- Higgs vev $v = 246.220569$ GeV, fixed by $G_F$ (SMINPUTS 2).
- Couplings to the Higgs carry the Clebsch-Gordan normalisation of their SU(2) invariant:
  factors such as $\sqrt2/2$ (doublet), $\sqrt3/6$ (neutral triplet), $\sqrt6/6$ (charged
  triplet) in the mass matrices.
- Complex couplings are $y = \mathrm{Re}\,y + i\,\mathrm{Im}\,y$, with `Ry...` and `Iy...`
  separate param_card entries (blocks `<NAME>` and `IM<NAME>`). The mass matrices are
  written in $y$; where $y$ enters together with its conjugate (real scalar blocks), in
  $\mathrm{Re}\,y$ and $\mathrm{Im}\,y$.
- In the reports a parameter `MChiI1Y0` is written $M_{\mathrm{ChiI1Y0}}$, `yChiI1Y0ChiI2Y3`
  is $y_{\mathrm{ChiI1Y0ChiI2Y3}}$, and the rows and columns are the gauge components of
  the FeynRules model.
