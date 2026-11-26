# uvbuilder conventions

The conventions of the models `uvbuilder build` writes and of the masses and mixing
matrices `uvbuilder fill` computes. Each model report (`<model>.md`) lists the model's
states, input parameters and mass matrices; this file explains them.

## Mass terms

A multiplet of isospin $j$ and hypercharge $Y$ has the components $\psi_m$,
$m = j, j-1, \dots, -j$, of electric charge $Q = m + Y$. Every mass term is an SU(2)
singlet; the formulas below are those of `Multiplet.mass_partner`.

### Scalars

```math
-\mathcal{L} \supset M^2 \sum_m \phi_m^* \phi_m \quad \text{(complex multiplet)}, \qquad
-\mathcal{L} \supset \frac{1}{2} M^2 \sum_m (-1)^m \phi_m \phi_{-m} \quad \text{(real multiplet)},
```

where a real multiplet ($Y = 0$, integer $j$) satisfies $\phi_m^* = (-1)^m \phi_{-m}$, so
that the second form is also $\frac12 M^2 \sum_m |\phi_m|^2$. A real singlet has
$-\mathcal{L} \supset \frac12 M^2 \phi^2$. The input is $M$, which enters squared.

### Fermions

For left-handed Weyl multiplets $\psi$ (isospin $j$, hypercharge $Y$) and $\chi$
(isospin $j$, hypercharge $-Y$), the Dirac mass term is

```math
-\mathcal{L} \supset M \sum_m c_m\, \psi_m \chi_{-m} + \mathrm{h.c.}, \qquad
c_m = s\,(-1)^{j+m},
```

and the Majorana mass term of a multiplet $\psi$ with $Y = 0$ and integer $j$ is

```math
-\mathcal{L} \supset \frac{1}{2} M \sum_m c_m\, \psi_m \psi_{-m} + \mathrm{h.c.}, \qquad
c_m = s\,(-1)^{j+m}.
```

The factors $(-1)^{j+m}$ make the sum an SU(2) singlet (Condon-Shortley phases; Dreiner,
Haber, Martin, Phys. Rept. 494 (2010) 1, App. E.1) and fix the relative signs of the
pairs. The overall sign $s = \pm 1$ is the convention: it is chosen so that
$c_{m_0} = +1$, i.e. $s = (-1)^{j+m_0}$, for

- a Majorana mass: $m_0 = 0$, the neutral component, paired with itself;
- a Dirac mass: $m_0$ the component of $\psi$ with the largest $|Q|$ (the larger $m$ if
  two have the same $|Q|$), i.e. the pair with the largest electric charge.

The other pairs then have $c_m = (-1)^{m - m_0}$. The result does not depend on which of
the two multiplets of a Dirac mass is called $\psi$. Examples:

- Two doublets, $\psi = (\psi^+, \psi^0)$ with $Y = 1/2$ and $\chi = (\chi^0, \chi^-)$
  with $Y = -1/2$: $m_0 = 1/2$, and
  $-\mathcal{L} \supset M(\psi^+\chi^- - \psi^0\chi^0) + \mathrm{h.c.}$: $+M$ for the charged
  pair, $-M$ for the neutral one (as $\mu$ in the SLHA chargino and neutralino matrices).
- A Majorana triplet $(\psi^+, \psi^0, \psi^-)$: $m_0 = 0$, and
  $-\mathcal{L} \supset \frac12 M \psi^0\psi^0 - M \psi^+\psi^- + \mathrm{h.c.}$: $+M$ for the
  neutral component, $-M$ for the charged pair.

### Positive masses

With these signs a positive $M$ gives positive masses. The negative entries ($-M$ above)
and negative inputs are absorbed into the mass eigenstates:

- a state alone in its block gets the phase $P = \mathrm{sign}(cM)$ (Dirac, on the
  right-handed partner) or $P = \sqrt{\mathrm{sign}(cM)}$ (Majorana), where $cM$ is
  its entry in the mass matrix. E.g. a pure Majorana triplet: $P = 1$ for the neutral
  state, $P = -1$ for the charged one;
- in a block of several states, by the phases of the mixing matrices (Takagi or
  singular value decomposition, below).

## Mass matrices

| block | matrix | $-\mathcal{L}$ | diagonalization |
|---|---|---|---|
| Majorana fermion | complex symmetric | $\frac12 \psi^T M \psi + \mathrm{h.c.}$ | Takagi: $U^T M U = \mathrm{diag}(m)$, $U$ unitary |
| Dirac fermion | complex | $\psi^T M \chi + \mathrm{h.c.}$ ($\psi$: rows, $\chi$: columns) | SVD: $U^T M V = \mathrm{diag}(m)$, $U$, $V$ unitary |
| real scalar | real symmetric | $\frac12 \rho^T M^2 \rho$ | $O^T M^2 O = \mathrm{diag}(m^2)$, $O$ real orthogonal |
| complex scalar | Hermitian | $\phi^\dagger M^2 \phi$ (rows: $\phi^*$, columns: $\phi$) | $U^\dagger M^2 U = \mathrm{diag}(m^2)$, $U$ unitary |

- Gauge eigenstates $= U \times$ mass eigenstates, with the columns in ascending mass:
  counting from 1, as the param_card entries (`i j`, row $i$, column $j$) and the
  FeynRules family indices do, column $j$ is the $j$-th lightest state. Each column is rephased so that its entry of
  largest modulus is real and positive (real scalars and Majorana fermions: only a sign
  is free). Masses must be positive and non-degenerate within a block.
- A state alone in its block has mass $|M|$ (fermions) or $\sqrt{M^2}$ (scalars); for
  fermions the sign of $M$ goes into the phase $P$ above. Its param_card entry is only
  its MASS.
- The param_card block of a mixing matrix is its name in capitals (`UMIX...`), with the
  imaginary parts in `IMUMIX...`; Dirac blocks have two matrices, $U$ and $V$.
- A neutral complex scalar $\Phi$ mixing with real states is written in the real basis
  $\Phi = (\Phi_R + i\Phi_I)/\sqrt2$: its block has the rows $\Phi_R$ and $\Phi_I$.

## Parameters

- Higgs vev $v = 246.220569$ GeV, fixed by $G_F$ (SMINPUTS 2).
- Couplings to the Higgs carry the Clebsch-Gordan normalization of their SU(2) invariant:
  factors such as $\sqrt2/2$ (doublet), $\sqrt3/6$ (neutral triplet), $\sqrt6/6$ (charged
  triplet) in the mass matrices.
- Complex couplings are $y = \operatorname{Re} y + i \operatorname{Im} y$, with `Ry...` and `Iy...`
  separate param_card entries (blocks `<NAME>` and `IM<NAME>`). The mass matrices are
  written in $y$; where $y$ enters together with its conjugate (real scalar blocks), in
  $\operatorname{Re} y$ and $\operatorname{Im} y$.
- In the reports a parameter `MChiI1Y0` is written $M_{\mathrm{ChiI1Y0}}$, `yChiI1Y0ChiI2Y3`
  is $y_{\mathrm{ChiI1Y0ChiI2Y3}}$, and the rows and columns are the gauge components of
  the FeynRules model.

## Names

### Field names

A multiplet name lists its quantum numbers, e.g. `ChiC3I1Y4Lw`, `PhiI2Y3`:

| part | meaning |
|---|---|
| `Phi`, `Chi` | scalar, fermion |
| `C3`, `C3b`, ... | color representation (`b`: conjugate); absent for a singlet |
| `I<n>` | SU(2) dimension |
| `Y<k>` | hypercharge $k/6$; `m` stands for a minus sign (`Ym3`: $-1/2$) |
| `B<k>`, `L<k>`, `D<k>` | global charges, when tracked (`-q`): baryon number $k/3$, lepton number and dark charge $k$ (`Bm1`: $B = -1/3$) |
| `F1`, `F2`, ... | copies of the same representation |
| `Lw`, `Rw` | left- and right-handed Weyl spinor |

A component adds its electric charge in units of 1/3, `q<k>` (`ChiI2Y3q3`: charge +1,
`ChiI2Ym3qm3`: charge -1), and a conjugate field ends in `bar`. SM fields keep the names
of `SMMod.fr`: `Phi` is the Higgs doublet, `tRbar` the conjugate right-handed top.

### Parameter names

The name of a Lagrangian parameter is a prefix followed by the names of the fields of its
term:

| prefix | term |
|---|---|
| `M` | mass: of a scalar (it enters squared, $M^2$) or of a Majorana fermion (`MChiI1Y0`), or Dirac mass of a fermion pair (`MChiI2Y3ChiI2Ym3`) |
| `mu` | cubic coupling of scalars with the Higgs doublet |
| `y` | Yukawa coupling of two fermions with the Higgs doublet |
| `z` | coupling of three fields without the Higgs doublet, e.g. two new fields and an SM fermion (`zPhiC3I1Y4ChiI1Y0tRbar`) |
| `lambda` | quartic coupling |

- A suffix `P1`, `P2`, ... distinguishes the independent SU(2) contractions of the
  same fields (e.g. `lambdaPhiI2Y3PhiI2Y3barP1`, `...P2`).
- Masses are real. The coupling of a term that is not its own Hermitian conjugate is
  complex: its real and imaginary parts are separate inputs, prefixed `R` and `I`
  (`RyChiI1Y0ChiI2Y3`, `IyChiI1Y0ChiI2Y3`).
- Internal parameters: the mixing matrices (below) and the phases `P<state>` of states
  alone in their block (see Mass terms).

### param_card blocks

The inputs are grouped by the fields of their term:

| block | terms |
|---|---|
| `BSMMASS` | mass terms (two fields) |
| `BSMTRILINEAR` | three scalars |
| `BSMYUKAWA` | two fermions and a scalar |
| `BSMQUARTIC` | four scalars |
| `BSMOTHER` | any other term |

The imaginary parts of complex couplings are in the same entry of `IM<block>`
(`IMBSMYUKAWA`, ...).

The mixing matrix of a block of states is named `UMix<type><states>`, with the type `M`
(Majorana), `D` (Dirac), `R` (real scalar) or `C` (complex scalar) and the labels of its
mass eigenstates, e.g. `UMixMZ1Z2Z3`; a Dirac block has a second matrix `VMix...`. Like
the complex parameters, a complex mixing matrix is split into its real and imaginary
parts, in the blocks `UMIX...` and `IMUMIX...` (as the FeynRules MSSM does with `NMIX`
and `IMNMIX`); a real scalar block (orthogonal matrix) has only `UMIX...`. These blocks
and the masses are written by `uvbuilder fill`.
