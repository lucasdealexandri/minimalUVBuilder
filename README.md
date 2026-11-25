# minimalUVBuilder

Minimal UV Lagrangians from LHC signatures. Given a *protomodel* (the masses of new
states and their decays, e.g. from the SModelS protomodel builder), `uvbuilder` finds
the minimal sets of SU(3) x SU(2) x U(1) multiplets whose renormalisable couplings
produce every decay at tree level, and writes each as a FeynRules (`.fr`) file, along
with its object representation (`.pkl`) and a model report (`.md`).

## Installation

Python >= 3.12:

```
git clone https://github.com/UVBuilder/minimalUVBuilder.git
cd minimalUVBuilder
pip install -e .            # pip install -e '.[test]' for the tests
```

FeynRules (with Mathematica 14) and MadGraph5_aMC@NLO are needed to turn a model into a
UFO and generate its param_card.

## Workflow

```
protomodel (.dict or .py)
   |  uvbuilder build -o FOLDER
   v
FOLDER/<model>/<model>.fr .pkl .md         one set per model
   |  FeynRules: WriteUFO; MadGraph5: param_card.dat
   v
<model>_UFO/param_card.dat                 set the input parameters
   |  uvbuilder fill <model>.pkl param_card.dat
   v
param_card with masses, mixing matrices, DECAY Auto  ->  MadGraph5
```

1. **Build.** `uvbuilder build PROTOMODEL` lists the models; with `-o FOLDER` it writes,
   for each, a folder `<model>` with the FeynRules file `.fr`, its object representation
   `.pkl` (the input of `fill`) and the model report `.md`.
2. **UFO.** Load the `.fr` in FeynRules together with the SM file `SMMod.fr` (shipped in
   `src/uvbuilder/data/`; copy it to the FeynRules folder or give its full path) and
   write the UFO:
   ```
   Unprotect[Commutator, NonCommutativeMultiply, MatrixSymbol];
   FR$Parallel = False;
   $FeynRulesPath = SetDirectory["/path/to/feynrules-current"];
   << FeynRules`
   LoadModel["SMMod.fr", "/path/to/<model>.fr"];
   WriteUFO[LSM + WeylToDirac[LBSM], Output -> "/path/to/<model>_UFO", AddDecays -> False];
   ```
   Then generate the UFO's `param_card.dat` with MadGraph5 (see `docs/tutorial.md`).
3. **Choose the inputs.** The model report lists the input parameters (param_card entry,
   type, the mass matrices each one enters) and the symbolic mass matrices. Set the
   inputs in the param_card.
4. **Fill.** `uvbuilder fill <model>.pkl <model>_UFO/param_card.dat` computes the masses
   and mixing matrices from the inputs and writes them into the card (over it, or to
   `-o FILE`), with the widths set to `DECAY Auto`; the branching ratios are computed by
   MadGraph5.

## Protomodels

A protomodel is a dictionary of masses, decays and labels, with optional settings of the
models. The stop example of the tutorial (`examples/stop/stop.dict`), with every option:

```
{
    "masses": {1000022: 200.0, 1000006: 1000.0},          # GeV
    "decays": {1000022: {}, 1000006: {(1000022, 6): 1.0}}, # {mother: {(daughters): BR}}
    "labels": {1000022: 'X1', 1000006: 'T1'},              # default X<pdg>
    "name": 'stop',                                        # (optional)
    "author": ['Andre Lessa', 'Lucas Heck'],               # (optional) default UVBuilder
    "version": '1.0',                                      # (optional)
    "dm_spin": 'both',                                     # (optional) scalar, fermion or both
    "dm_charges": ['none'],                                # (optional) e.g. ['B=0,D=1']
}
```

- PDG 1000022 is always the dark matter candidate.
- Decays are two-body; an off-shell decay is written with its virtual boson, e.g.
  `(1000022, 23)` instead of `(1000022, 12, -12)`. Decays to SM states only are skipped.
- The options take the values of the command-line options of the same name, which
  override them.

It is given as a `.dict` file, or as a Python file defining the dictionary `protomodel`
and optionally `possible_dm_reps`, `dm_pdg`, `labeler`, `name`, `author` (a name or a
list) and `version`.

## Command line

```
uvbuilder build PROTOMODEL [-o FOLDER] [--name NAME] [-s scalar|fermion|both]
                           [-q none|B=0,D=1 ...] [--author NAME ...] [--version V]
uvbuilder fill MODEL.pkl CARD [-o OUTPUT]
uvbuilder -v ...             debug messages, e.g. why candidate models were rejected
```

- `-o`, `--out`: folder of the models, one subfolder per model.
- `--name`: model name, default the protomodel file name (`model-1.dict` -> `model_1`);
  the models are numbered `<name>_1`, `<name>_2`, ...
- `-s`, `--dm-spin`: spin of the DM candidate, default both.
- `-q`, `--dm-charges`: global U(1) charges of the DM (`B`, `L`, `D`); default `none`, no
  charge tracked. Repeat the option for alternatives.
- `--author`, `--version`: `M$Information` of the `.fr` files, default `UVBuilder` and
  `1.0`; repeat `--author` for several authors.
- These options override those of the protomodel file.

Note: nothing is written without `--out`.

## Conventions

The conventions are listed in `docs/conventions.md`. They comprise the convention of the
mass terms, under which positive mass parameters give positive masses; the
diagonalisation of each type of mass matrix (Majorana and Dirac fermions, real and
complex scalars) and the ordering and phases of its mixing matrix, with the gauge
eigenstates given by U times the mass eigenstates and the columns in ascending mass; the
param_card blocks of the mixing matrices; the real basis of neutral complex scalars; the
treatment of complex couplings; and the notation of the model reports.

Validity is checked at tree level: a decay of one new state into another and an SM
particle needs a coupling between their fields that emits that particle. The masses are
not checked: a valid model may be unable to reproduce the protomodel's spectrum (e.g. a
degenerate one), which `fill` reports.

## Examples

`examples/` holds the protomodels and their models (`.fr`, `.pkl`, `.md`, UFOs); the
tutorial (`docs/tutorial.md`) goes through them.

| folder | protomodel |
|---|---|
| `stop` | stop -> top + DM |
| `sddm` | singlet-doublet dark matter |
| `idm` | inert doublet |
| `ewinos` | electroweakinos from an MSSM SLHA point (the EW-MSSM among the models) |

## Tests

```
pytest              # fast suite
pytest -m slow      # long-running tests
```

## Layout

```
src/uvbuilder/
  cli.py              command line interface
  modelBuilder/       protomodel -> candidate multiplet models (ModelBuilder)
  models/             final states, particles, multiplets, FeynRules output; param_card
                      builder
  partition/          partitions of states into multiplets, mixing blocks
  states/             SM and BSM fields, particles, multiplets
  vertices/           decay and particle vertices
  tools/              representations, operators, Clebsch-Gordan tensors, FeynRules names,
                      spectrum (diagonalisation)
  paramCard.py        ParamCardGenerator (uvbuilder fill)
  report.py           model report (.md)
  data/SMMod.fr       SM model file for FeynRules
docs/conventions.md   conventions
docs/tutorial.md      tutorial
examples/             protomodels and their models
scripts/dump_reserved_symbols.wl
                      dumps the Mathematica and FeynRules symbol lists into
                      src/uvbuilder/tools/reserved/ (names a model must not use)
```

## Authors

Andre Lessa and Lucas Heck.
