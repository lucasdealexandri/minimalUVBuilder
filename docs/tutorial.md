# Tutorial

## Goal

Starting from a dictionary of new particles, with their masses and decay channels, UV
Builder finds the minimal UV Lagrangians that produce every decay at tree level and
writes each as a FeynRules model file. For each model it then fills param_cards
self-consistently: from the Lagrangian parameters you choose, it computes the masses and
mixing matrices. Choosing parameters that reproduce the masses, branching ratios and
cross sections of the dictionary is up to you; the report written for each model lists
the parameters and the mass matrices to help with it.

The examples, from [`examples/`](../examples/):

1. a stop decaying to a top and the dark matter: the full workflow, from the dictionary
   to a filled param_card;
2. singlet-doublet dark matter: a mass matrix with mixing;
3. a dictionary whose valid model cannot reproduce its masses.

They assume `uvbuilder` is installed (see the [README](../README.md)), and FeynRules and MadGraph5 for
the UFO steps. The conventions are in [`docs/conventions.md`](conventions.md).

## 1. Stop -> top + DM

### Writing the protomodel

[`examples/stop/stop.dict`](../examples/stop/stop.dict):

```
{
    "masses": {1000022: 200.0, 1000006: 1000.0},          # GeV
    "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}}, # {mother: {(daughters): BR}}
    "labels": {1000022: 'X1', 1000006: 'T1'},              # default X<pdg>
    "name": 'stop',                                        # (optional)
    "author": ['Andre Lessa', 'Lucas Heck'],               # (optional) default UVBuilder
    "version": '1.0',                                      # (optional)
    "dm_spin": 'both',                                     # (optional) scalar, fermion or both
    "dm_charges": ['B=0'],                                 # (optional) e.g. ['B=0,D=1']
}
```

A state T1 of 1000 GeV decays to a top quark and the dark matter (DM) candidate X1 of
200 GeV; PDG 1000022 is always the DM candidate. The last two entries set the options
`-s` (`--dm-spin`) and `-q` (`--dm-charges`), which the command line overrides:

- `dm_spin: 'both'` (`-s both`, the default): the DM candidate may be a scalar or a
  fermion; `-s scalar` or `-s fermion` keeps one;
- `dm_charges: ['B=0']` (`-q B=0`): the DM has baryon number 0 and baryon number is
  tracked: every new multiplet is assigned the baryon number its couplings require, and
  only couplings that conserve it are written. The default, `-q none`, tracks no global
  charge.

### Building the models

```
cd examples/stop
uvbuilder build stop.dict
```

```
stop_1
  multiplet         spin           SU(3)  SU(2)  Y     charges
  ChiC3I1Y4B1Lw     LH WeylSpinor  3      1      2/3   B=1/3
  ChiC3bI1Ym4Bm1Lw  LH WeylSpinor  3*     1      -2/3  B=-1/3
  PhiI1Y0           Scalar         1      1      0     B=0
  state  kind           components
  T1     Dirac fermion  ChiC3I1Y4B1q2 | ChiC3bI1Ym4Bm1qm2
  X1     real scalar    PhiI1Y0q0

stop_2
  multiplet    spin           SU(3)  SU(2)  Y    charges
  ChiI1Y0Lw    LH WeylSpinor  1      1      0    B=0
  PhiC3I1Y4B1  Scalar         3      1      2/3  B=1/3
  state  kind              components
  T1     complex scalar    PhiC3I1Y4B1q2
  X1     Majorana fermion  ChiI1Y0q0

nothing written: use -o FOLDER to write each model's .fr, .pkl and .md
```

Each model lists its multiplets, with their quantum numbers, and the states of the
protomodel, with the gauge components they are made of (for a Dirac fermion, the two
Weyl spinors separated by `|`; states that mix share a line). The builder finds two
models in which T1 -> t X1 comes from a renormalizable coupling:

- `stop_1`: T1 is a Dirac fermion, a color triplet with hypercharge 2/3 and baryon
  number 1/3 (the pair of Weyl spinors `ChiC3I1Y4B1Lw`, `ChiC3bI1Ym4Bm1Lw`), and X1 a
  real scalar singlet (`PhiI1Y0`);
- `stop_2`: T1 is a complex scalar with the quantum numbers of a right-handed stop
  (`PhiC3I1Y4B1`), and X1 a Majorana fermion singlet (`ChiI1Y0Lw`).

The names encode the representations (see [`docs/conventions.md`](conventions.md#field-names)): `Phi` a scalar, `Chi`
a fermion, `Lw` a left-handed Weyl spinor, `C3` (`C3b`) a color triplet (antitriplet),
`I<n>` the SU(2) dimension, `Y<k>` the hypercharge in units of 1/6, `B<k>` the baryon
number in units of 1/3, `m` a minus sign.

Without baryon number, the same models are found, without the `B` tags, e.g.

```
uvbuilder build stop.dict -q none
```

```
...
stop_2
  multiplet  spin           SU(3)  SU(2)  Y    charges
  ChiI1Y0Lw  LH WeylSpinor  1      1      0
  PhiC3I1Y4  Scalar         3      1      2/3
  state  kind              components
  T1     complex scalar    PhiC3I1Y4q2
  X1     Majorana fermion  ChiI1Y0q0
```

Here the Lagrangians are the same: no coupling of these fields violates baryon number.
The rest of this section uses `stop.dict` as it is. Write the models:

```
uvbuilder build stop.dict -o .
```

```
stop_1/  stop_1.fr  stop_1.pkl  stop_1.md  [ChiC3I1Y4B1Lw, ChiC3bI1Ym4Bm1Lw, PhiI1Y0]
stop_2/  stop_2.fr  stop_2.pkl  stop_2.md  [ChiI1Y0Lw, PhiC3I1Y4B1]
```

Each folder holds the FeynRules model file (`.fr`), its object representation (`.pkl`,
read by `uvbuilder fill`) and the model report (`.md`). The rest of this section follows
`stop_2`.

### Reading the model report

[`stop_2/stop_2.md`](../examples/stop/stop_2/stop_2.md) lists the mass eigenstates, the input parameters with their
param_card entries, and the mass matrices. The inputs:

| name | card | type | enters |
| --- | --- | --- | --- |
| MChiI1Y0 | BSMMASS 1 | mass | MX1 |
| MPhiC3I1Y4B1 | BSMMASS 2 | M² | CT1 |
| RzPhiC3I1Y4B1ChiI1Y0cRbar | BSMYUKAWA 1 | Re | - |
| IzPhiC3I1Y4B1ChiI1Y0cRbar | IMBSMYUKAWA 1 | Im | - |
| RzPhiC3I1Y4B1ChiI1Y0tRbar | BSMYUKAWA 2 | Re | - |
| IzPhiC3I1Y4B1ChiI1Y0tRbar | IMBSMYUKAWA 2 | Im | - |
| RzPhiC3I1Y4B1ChiI1Y0uRbar | BSMYUKAWA 3 | Re | - |
| IzPhiC3I1Y4B1ChiI1Y0uRbar | IMBSMYUKAWA 3 | Im | - |
| lambdaPhiC3I1Y4B1PhiC3I1Y4B1bar | BSMQUARTIC 1 | coupling | CT1 |

- `MChiI1Y0`: the Majorana mass of X1.
- `MPhiC3I1Y4B1`: the mass parameter of the stop, which enters squared.
- `z...tRbar`, `z...cRbar`, `z...uRbar`: the couplings of the stop to X1 and the
  right-handed top, charm and up quarks. They are complex, so their real (`R`) and
  imaginary (`I`) parts are separate entries, in the blocks `BSMYUKAWA` and
  `IMBSMYUKAWA`. The decay to a top needs only the first.
- `lambda...`: the coupling of the stop to the Higgs, which shifts its mass.

The parameter names and blocks are described in [`docs/conventions.md`](conventions.md#parameter-names). Each state is
alone in its mass matrix, so its mass is a single entry of the report:

```math
m_{X1} = |M_{\mathrm{ChiI1Y0}}|, \qquad
m_{T1}^2 = M_{\mathrm{PhiC3I1Y4B1}}^2 + \frac{\lambda_{\mathrm{PhiC3I1Y4B1PhiC3I1Y4B1bar}}\, v^2}{2},
```

with the Higgs vev v = 246.22 GeV.

### Generating the UFO

[FeynRules](https://feynrules.irmp.ucl.ac.be/) writes the [UFO](https://arxiv.org/abs/1108.2040) of the model. In Mathematica, with the SM model file of
UV Builder ([`src/uvbuilder/data/SMMod.fr`](../src/uvbuilder/data/SMMod.fr)) copied to the FeynRules folder:

```
Unprotect[Commutator, NonCommutativeMultiply, MatrixSymbol];
FR$Parallel = False;
$FeynRulesPath = SetDirectory["/path/to/feynrules-current"];
<< FeynRules`
LoadModel["SMMod.fr", "/path/to/examples/stop/stop_2/stop_2.fr"];
WriteUFO[LSM + WeylToDirac[LBSM], Output -> "/path/to/examples/stop/stop_2/stop_2_UFO",
         AddDecays -> False];
```

`LBSM` is the Lagrangian of the new fields, written with Weyl spinors; `WeylToDirac`
turns it into Dirac and Majorana spinors.

### Generating the param_card

`WriteUFO` does not write the param_card. The UFO written by FeynRules uses Python 2
syntax, which recent MadGraph5 versions (Python 3) convert with the command
`convert model` ([MadGraph5 pull request #446](https://github.com/mg5amcnlo/mg5amcnlo/pull/446)). In MadGraph5 (`./bin/mg5_aMC`):

```
convert model /path/to/examples/stop/stop_2/stop_2_UFO
```

then write the param_card of the UFO:

```
cd /path/to/examples/stop/stop_2/stop_2_UFO
python write_param_card.py
```

### Choosing the inputs

To reproduce the protomodel masses, take X1 at 200 GeV and the couplings of the stop to
the top and to the Higgs at 0.5; the stop mass parameter then follows from the formula
above:

```math
M_{\mathrm{PhiC3I1Y4B1}} = \sqrt{1000^2 - 0.5\, v^2/2} = 992.392996 \ \mathrm{GeV}.
```

The couplings to the charm and up quarks are set to zero, so that T1 decays only to the
top. In [`stop_2_UFO/param_card.dat`](../examples/stop/stop_2/stop_2_UFO/param_card.dat):

```
Block BSMMASS
    1 2.00000000e+02 # MChiI1Y0
    2 9.92392996e+02 # MPhiC3I1Y4B1
Block BSMQUARTIC
    1 5.00000000e-01 # lambdaPhiC3I1Y4B1PhiC3I1Y4B1bar
Block BSMYUKAWA
    1 0.00000000e+00 # RzPhiC3I1Y4B1ChiI1Y0cRbar
    2 5.00000000e-01 # RzPhiC3I1Y4B1ChiI1Y0tRbar
    3 0.00000000e+00 # RzPhiC3I1Y4B1ChiI1Y0uRbar
Block IMBSMYUKAWA
    1 0.000000e+00 # IzPhiC3I1Y4B1ChiI1Y0cRbar
    2 0.000000e+00 # IzPhiC3I1Y4B1ChiI1Y0tRbar
    3 0.000000e+00 # IzPhiC3I1Y4B1ChiI1Y0uRbar
```

### Filling the param_card

```
uvbuilder fill stop_2/stop_2.pkl stop_2/stop_2_UFO/param_card.dat
```

computes the masses (and, for states that mix, the mixing matrices) from the inputs,
writes them into the card and sets the widths of the new states to `Auto`:

```
Block MASS
  1000006 1.0000000003342869e+03 # MT1
  1000022 2.0000000000000000e+02 # MX1
...
DECAY 1000006 Auto
DECAY 1000022 Auto
```

T1 is at 1000 GeV up to the rounding of the input. The card is ready for
[MadGraph5_aMC@NLO](https://launchpad.net/mg5amcnlo), which computes the widths, branching ratios and cross sections.

For `stop_1` the steps are the same. There X1 is a real scalar with

```math
m_{X1}^2 = M_{\mathrm{PhiI1Y0}}^2 + \lambda_{\mathrm{PhiI1Y0PhiI1Y0}}\, v^2,
```

so with λ = 0.1 the input `MPhiI1Y0` = 184.221451 gives 200 GeV; the Dirac mass
`MChiC3I1Y4B1ChiC3bI1Ym4Bm1` = 1000 gives T1, and the coupling to the top
`RzChiC3I1Y4B1PhiI1Y0tRbar` is 0.5, the others zero
([`stop_1_UFO/param_card.dat`](../examples/stop/stop_1/stop_1_UFO/param_card.dat)).

## 2. Singlet-doublet dark matter: mixing

Singlet-doublet dark matter ([Mahbubani, Senatore](https://arxiv.org/abs/hep-ph/0510064);
[Cohen et al.](https://arxiv.org/abs/1109.2604)) is a simple model in which the DM candidate is a
mixture of states from different multiplets, mixed by their couplings to the Higgs. It
is used here to look at mixing and at the consistent filling of the mixing matrices in
the param_card.

[`examples/sddm/sddm.dict`](../examples/sddm/sddm.dict) has three neutral states Z1, Z2, Z3 and a charged one W1, with
Z2 -> Z1 Z, W1 -> Z1 W and Z3 -> Z2 h. The builder finds a scalar and a fermion model:

```
cd examples/sddm
uvbuilder build sddm.dict -o .
```

```
sddm_1/  ...  [PhiI1Y0, PhiI2Y3]
sddm_2/  ...  [ChiI1Y0Lw, ChiI2Y3Lw, ChiI2Ym3Lw]
```

`sddm_2` is a fermion singlet and a pair of fermion doublets with opposite hypercharge,
like a bino and two higgsinos. Its report ([`sddm_2/sddm_2.md`](../examples/sddm/sddm_2/sddm_2.md)) gives the mass matrix of
the neutral states, a complex symmetric matrix in the basis of the neutral components:

```math
\begin{array}{c|ccc}
 & \mathrm{ChiI1Y0q0} & \mathrm{ChiI2Y3q0} & \mathrm{ChiI2Ym3q0} \\ \hline
\mathrm{ChiI1Y0q0} & M_{\mathrm{ChiI1Y0}} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Y3}}}{2} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Ym3}}}{2} \\
\mathrm{ChiI2Y3q0} & - \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Y3}}}{2} & 0 & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\mathrm{ChiI2Ym3q0} & \frac{\sqrt{2} v y_{\mathrm{ChiI1Y0ChiI2Ym3}}}{2} & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} & 0 \\
\end{array}
```

It is diagonalized by a Takagi factorization, U^T M U = diag(m) with U unitary and m
non-negative ([Dreiner, Haber, Martin](https://arxiv.org/abs/0812.1594);
[Hahn](https://arxiv.org/abs/physics/0607103)). The charged state W1 has the Dirac mass `MChiI2Y3ChiI2Ym3` alone.
With the inputs

| name | value |
| --- | --- |
| MChiI1Y0 | 250 |
| MChiI2Y3ChiI2Ym3 | 200 |
| RyChiI1Y0ChiI2Y3 | 0.2 |
| RyChiI1Y0ChiI2Ym3 | -1.5 |
| IyChiI1Y0ChiI2Y3, IyChiI1Y0ChiI2Ym3 | 0 |

`uvbuilder fill` (after the UFO and param_card steps of section 1) gives the masses of
the protomodel, Z1 = 109.281, Z2 = 289.840, W1 = 200 and Z3 = 430.559 GeV, and the mixing
matrix U, with the gauge eigenstates given by U times the mass eigenstates:

```
             Z1        Z2         Z3
ChiI1Y0q0  -0.4678   0.3912 i   0.7926
ChiI2Y3q0   0.8073   0.5541 i   0.2030
ChiI2Ym3q0 -0.3597   0.7348 i  -0.5750
```

As for the complex parameters, the mixing matrices are split into their real and
imaginary parts, written in two blocks of the param_card, here `UMIXMZ1Z2Z3` and
`IMUMIXMZ1Z2Z3`, following the [MSSM implementation of FeynRules](https://arxiv.org/abs/1102.4191) (e.g. its
neutralino mixing matrix in `NMIX` and `IMNMIX`). The column of Z2 is imaginary: the
real symmetric mass matrix has a negative eigenvalue there, and the factor i makes the
physical mass positive.

## 3. A valid model that cannot reproduce the masses

[`examples/idm/idm.dict`](../examples/idm/idm.dict) has the states of an inert doublet ([Deshpande, Ma](https://doi.org/10.1103/PhysRevD.18.2574);
[Barbieri, Hall, Rychkov](https://arxiv.org/abs/hep-ph/0603188)): H0 (100 GeV),
A0 (200 GeV) and Hp (250 GeV), with A0 -> H0 Z and Hp -> H0 W. Besides the inert doublet
itself (`idm_1`), the builder finds a pair of fermion doublets (`idm_2`,
`[ChiI2Y3Lw, ChiI2Ym3Lw]`). Its report ([`idm_2/idm_2.md`](../examples/idm/idm_2/idm_2.md)) shows that it cannot reproduce the protomodel: a
single parameter fixes every mass,

```math
\begin{array}{c|cc}
 & \mathrm{ChiI2Y3q0} & \mathrm{ChiI2Ym3q0} \\ \hline
\mathrm{ChiI2Y3q0} & 0 & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} \\
\mathrm{ChiI2Ym3q0} & - M_{\mathrm{ChiI2Y3ChiI2Ym3}} & 0 \\
\end{array}
```

so H0 and A0 are degenerate, and Hp has the same mass: two doublets with opposite
hypercharge have no renormalizable coupling to the Higgs that would split them. The
model is valid, since validity checks the couplings at tree level and not the masses,
but `fill` stops:

```
uvbuilder fill: UMixMH0A0: degenerate masses (100, 100 GeV)
```
