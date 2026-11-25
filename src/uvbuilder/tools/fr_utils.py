#!/usr/bin/env python3
"""
.. module:: fr_utils
   :synopsis: Helpers to render the FeynRules model file.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

from uvbuilder.tools.reps import (U1Rep, SU2Rep, SU3Rep, DiracSpinor, WeylSpinor,
                        Scalar, Vector, WEYL_LEFT)


# Naming conventions: the single source of truth
LEFT_SUFFIX, RIGHT_SUFFIX, DAGGER_SUFFIX = 'Lw', 'Rw', 'bar'
# Unphysical FeynRules class for H~ = i sigma2 H* (our conjugate Higgs multiplet)
HIGGS_TILDE = 'Phitilde'
# Largest SU(2) dimension whose representation SMMod.fr declares (SU2d3 ... SU2d10)
MAX_SU2_DIM = 10


def su2_index_type(d: int) -> str:
    '''FeynRules index type of a BSM SU(2) multiplet of dimension d: SU2D for
    doublets, SU2d<d> (declared in SMMod.fr, spherical basis) for 3 <= d <= 10.
    Triplets never use SU2W, SM's Cartesian adjoint index of Wi.'''
    if d == 2:
        return 'SU2D'
    if d > MAX_SU2_DIM:
        raise ValueError(f'SU(2) multiplets of dimension {d} > {MAX_SU2_DIM} have no '
                         f'representation in SMMod.fr (extend SU2L there and MAX_SU2_DIM)')
    return f'SU2d{d}'


def chirality_suffix(spin) -> str:
    '''Weyl chirality suffix: LEFT_SUFFIX for WEYL_LEFT, RIGHT_SUFFIX otherwise.'''
    return LEFT_SUFFIX if spin == WEYL_LEFT else RIGHT_SUFFIX


# FeynRules names of the conserved global U(1)s (LeptonNumber as in SM.fr)
GLOBAL_U1_QN = {'B': 'BaryonNumber', 'L': 'LeptonNumber', 'D': 'DarkCharge'}


def fr_quantum_numbers(charges: Dict[str, object], globalU1s=None) -> str:
    '''QuantumNumbers value, e.g. {Q -> -1, LeptonNumber -> 1}: the given
    charges followed by the non-zero tracked global U(1)s (a charge that is
    zero or not conserved is omitted; FeynRules treats a missing one as 0).'''
    entries = [f'{name} -> {formatVariable(value)}' for name, value in charges.items()]
    if globalU1s is not None:
        entries += [f'{GLOBAL_U1_QN[name]} -> {formatVariable(value)}'
                    for name, value in globalU1s.items() if value != 0]
    return '{' + ', '.join(entries) + '}'


def toggle_dagger(label: str) -> str:
    '''Label involution X <-> Xbar. (FeynRules' own conjugate symbol always
    appends DAGGER_SUFFIX instead, e.g. Phibar -> Phibarbar.)'''
    return (label.removesuffix(DAGGER_SUFFIX) if label.endswith(DAGGER_SUFFIX)
            else label + DAGGER_SUFFIX)


FRChunk = Dict[str, list]

# Parameter and Lagrangian sectors, in output order: (key, title, SLHA block of
# the External parameters; imaginary parts go to IM<block>)
SECTORS = (
    ('Mass', 'Mass terms', 'BSMMASS'),
    ('Trilinear', 'Trilinear couplings', 'BSMTRILINEAR'),
    ('Yukawa', 'Yukawa couplings', 'BSMYUKAWA'),
    ('Quartic', 'Quartic couplings', 'BSMQUARTIC'),
    ('Other', 'Other couplings', 'BSMOTHER'),
    ('Mixing', 'Mixing matrices', None),
    ('Phase', 'Mass-eigenstate phases', None),
    ('CG', 'SU(2) Clebsch-Gordan tensors', None),
)
SECTOR_BLOCK = {key: block for key, _, block in SECTORS}

SECTION_KEYS = (
    'Classes',
    'Parameters',
    'KineticTerms',
) + tuple(f'{key}Terms' for key, _, _ in SECTORS[:5])


def merge_chunks(*chunks: FRChunk) -> FRChunk:
    '''Merge chunks: concatenate lists per section, deduplicating identical
    rendered strings (so repeated parameter/CG blocks appear once).'''
    out: FRChunk = {key: [] for key in SECTION_KEYS}
    for chunk in chunks:
        if not chunk:
            continue
        for key, items in chunk.items():
            out.setdefault(key, []).extend(items)
    for key, items in out.items():
        seen = set()
        deduped = []
        for item in items:
            if item not in seen:
                seen.add(item)
                deduped.append(item)
        out[key] = deduped
    return out


def formatVariable(obj) -> str:
    '''Format a variable depending on the object type to a string to be used
    when writing the FeynRules file. If no dedicated method is found use the
    string representation for the object.'''
    if isinstance(obj, DiracSpinor):
        return "F"
    elif isinstance(obj, WeylSpinor):
        return "W"
    elif isinstance(obj, Scalar):
        return "S"
    elif isinstance(obj, Vector):
        return "V"
    elif isinstance(obj, SU3Rep):
        if obj.upper == 1 and obj.lower == 1:
            return "Gluon"  # Default from SM definitions
        elif obj.upper == 1 and obj.lower == 0:
            return "Colour"  # Default from SM definitions
        elif obj.upper == 2 and obj.lower == 0:
            return "Sextet"  # Declared in SMMod.fr (representation T6)
        elif (obj.upper, obj.lower) in ((0, 1), (0, 2)):
            return f"SU3d{obj.dim}bar"
        raise ValueError(f'No FeynRules index type for color {obj!r} '
                         f'(supported: 3, 3bar, 6, 6bar, 8)')
    elif isinstance(obj, SU2Rep):
        return su2_index_type(obj.dim)
    elif isinstance(obj, list):
        return f"{{ {', '.join([formatVariable(x) for x in obj])} }}"
    else:
        return str(obj)


def formatClassBlock(header: str,
                     blockDict: Dict[str, Union[None, str]]) -> str:
    '''Format data corresponding to a single FeynRules class block.'''
    headerPadding = 3
    fieldRequirements = {
        'ClassName': True,
        'ClassMembers': False,
        'QuantumNumbers': False,
        'SelfConjugate': True,
        'Unphysical': False,
        'Chirality': False,
        'WeylComponents': False,
        'Mass': False,
        'Width': False,
        'PDG': False,
        'Indices': False,
        'FlavorIndex': False,
        'Definitions': False,
    }
    largestStr = max([len(key) for key in fieldRequirements])
    headerStr = f"{'':<{headerPadding}}{header} == {{"
    lines = []
    for key, required in fieldRequirements.items():
        if key not in blockDict:
            if required:
                raise KeyError(
                    f"Field {key} is required for a Class block, but it was "
                    f"not found in:\n{blockDict}\n")
            continue
        value = blockDict[key]
        if value is None:
            continue
        value = formatVariable(value)
        lines.append(
            f"{' ':<{len(headerStr)+1}}{key:<{largestStr}} -> {value}")

    strBlock = ',\n'.join(lines)
    strBlock = f"{headerStr}\n{strBlock}\n{' ':<{len(headerStr)}}}}"
    return strBlock


PARAMETER_OPTIONS = ('ParameterType', 'BlockName', 'OrderBlock', 'Indices', 'Unitary',
                     'Value', 'Definitions', 'ComplexParameter', 'InteractionOrder',
                     'TeX', 'Description')
_OPTION_WIDTH = max(map(len, PARAMETER_OPTIONS))


def formatParameterBlock(header: str, blockDict: Dict[str, object], indent: str = 2*' ') -> str:
    '''One M$Parameters entry: options in PARAMETER_OPTIONS order, names padded
    to a common width; continuation lines of multi-line values are aligned
    under the first value character.'''
    for key in ('ParameterType', 'ComplexParameter'):
        if key not in blockDict:
            raise KeyError(f"Field {key} is required for a Parameter block, but it "
                           f"was not found in:\n{blockDict}\n")
    lines = []
    for key in PARAMETER_OPTIONS:
        value = blockDict.get(key)
        if value is None:
            continue
        prefix = f'{2*indent}{key:<{_OPTION_WIDTH}} -> '
        lines.append(prefix + formatVariable(value).replace('\n', '\n' + ' '*len(prefix)))
    return f'{indent}{header} == {{\n' + ',\n'.join(lines) + f'\n{indent}}}'


@dataclass(frozen=True)
class FRParameter:
    '''A parameter entry of M$Parameters, rendered by the model assembler
    (which groups by sector and numbers OrderBlock within each block).
    ``sort_key`` keeps the real and imaginary parts of a coupling adjacent.'''
    name: str
    sector: str
    options: Tuple[Tuple[str, object], ...]
    sort_key: Tuple[str, int] = ('', 0)

    @property
    def kind(self) -> str:
        return dict(self.options)['ParameterType']

    @property
    def block(self) -> Optional[str]:
        return dict(self.options).get('BlockName')

    @property
    def order_block(self) -> Optional[int]:
        return dict(self.options).get('OrderBlock')

    @property
    def value(self):
        return dict(self.options).get('Value')

    @property
    def description(self) -> str:
        return str(dict(self.options).get('Description', '')).strip('"')

    @property
    def base(self) -> str:
        '''Name of the parameter this entry belongs to: the complex parameter for its
        R / I parts, else its own name.'''
        return self.sort_key[0]

    @property
    def part(self) -> str:
        '''R or I for the real / imaginary External part of a complex coupling, else "".'''
        return ('', 'R', 'I')[self.sort_key[1]]

    def render(self, **extra) -> str:
        '''The entry of the parameter in M$Parameters.'''
        return formatParameterBlock(self.name, {**dict(self.options), **extra})


def parameter_block(name: str, sector: str, parameter_type: str = 'External',
                    value=None, description: Optional[str] = None,
                    complex_parameter: bool = True,
                    interaction_order: Optional[str] = None,
                    block_name: Optional[str] = None, part: str = '') -> FRParameter:
    '''Parameter entry; ``part`` is 'R' or 'I' for the real/imaginary External
    parts of a complex coupling (``name`` is then the complex parameter's name).'''
    options = {'ParameterType': parameter_type, 'BlockName': block_name, 'Value': value,
               'ComplexParameter': complex_parameter, 'InteractionOrder': interaction_order,
               'Description': description}
    return FRParameter(part + name, sector,
                       tuple((k, v) for k, v in options.items() if v is not None),
                       (name, {'': 0, 'R': 1, 'I': 2}[part]))


def su2_fr_index(dim: int, i3) -> int:
    '''FeynRules index (1 to dim) of the SU(2) component with projection i3:
    1 = highest i3, as in the SMMod.fr generators (j = 1/2: i3 = 1/2 -> 1,
    i3 = -1/2 -> 2). Used for multiplet components and CG tensor entries alike.'''
    return int(U1Rep(dim - 1, 2) + 1 - i3)


def getSU2Index(field) -> int:
    '''FeynRules index of a multiplet component (see su2_fr_index).'''
    return su2_fr_index(field.isospin.dim, field.i3)


def _indexed(symbol: str, indices: List[str]) -> str:
    '''symbol[indices], or symbol alone without indices.'''
    return f"{symbol}[{', '.join(indices)}]" if indices else symbol


def getGaugeDefinitions(multiplet, mixings, symbol=None, conjugate=False, phases=None) -> str:
    '''
    Rendered definitions of Gauge Eigenstates -> Mass Eigenstates, e.g.
        DdLw[s_, 1] :> Module[{j}, UMixM[1, j] ChiQ0C1Lw[s, j]]
        DuLw[s_, 1] :> CPlusLw[s]                    (single-member class)
    Row = 1-based position of the field in its mixing matrix; the implicit sum
    over j runs over the family index of the physical class, which equals the
    matrix column index (PDG order, one block per class: see
    ParticleMixing.get_cluster and MultipletModel.valid_mixing).

    Real multiplets with isospin dim > 1 (Multiplet.is_real_nontrivial): component
    -m is (-1)^m times the conjugate of the definition of component m, e.g.
        PhiI3Y0[3] :> -W1bar                         (real triplet: T^- = -(T^+)^*)

    With conjugate=True, the definitions of the conjugate multiplet X^c under
    `symbol`: component -m is (-1)^(j+m) times the conjugate of the definition
    of component m of X, e.g. XRw[s_, c_] :> stRw[s, c].

    A complex component in a block of real states (real orthogonal matrix) uses its
    own row as the real part and its conjugate's row as the imaginary part,
    X = (phi_R + i phi_I)/sqrt(2), e.g.
        PhiI2Y3[2] :> Module[{j}, (UMixR[2, j] + I*UMixR[3, j])/Sqrt[2] PhiQ0C1[j]]
    [Dreiner, Haber, Martin, Phys. Rept. 494 (2010) 1, eqs. (3.2.41)-(3.2.42) and,
    for the MSSM Higgs sector, (K.1.3)-(K.1.8)].

    A fermion alone in its block carries its phase from `phases` ({particle: symbol},
    see MultipletModel.mass_phases): on the field for a Majorana state, on the
    right-handed partner for a Dirac state, e.g.
        HdLw[s_, 2] :> PW1 W1Rwbar[s]
    '''
    from uvbuilder.models.particleModel import ParticleModel  # local: particleModel imports fr_utils
    from uvbuilder.states.particle import family_key          # local: particle imports fr_utils
    from uvbuilder.tools.reps import su2_reality_phase

    mixing_matrices = mixings.get_mixing_matrices()
    shared = ((['s'] if isinstance(multiplet.rep.spin, WeylSpinor) else [])
              + (['c'] if multiplet.colordim.dim > 1 else []))

    def definition(field):
        '''(name, matrix symbol or None, row, imaginary-part row or None) of the
        physical expansion of field.'''
        mixing_matrix = next(mm for mm in mixing_matrices if field in mm.f_components)
        row = mixing_matrix.get_index(field) + 1
        # the components of this multiplet are the Phi of their (Phi, PhiBar) pairs;
        # Phi = (Phi_R + i Phi_I)/sqrt(2) puts the row of Phi_I in its row of C
        _, C = mixing_matrix.real_basis(multiplet.fields)
        row_im = next((k + 1 for k in range(mixing_matrix.dim)
                       if k != row - 1 and C[row - 1, k] != 0), None)

        # Physical class this field mixes into (mass order = family order)
        members = sorted(mixings.inverse[field], key=family_key)
        first = members[0].components[0]           # LH Weyl (fermions) / phi (scalars)
        conj = field.rep.unbroken != first.rep     # field is RH-bar / phi^dagger
        if isinstance(first.spin, WeylSpinor):
            suffix = RIGHT_SUFFIX + DAGGER_SUFFIX if conj else LEFT_SUFFIX
        else:
            suffix = DAGGER_SUFFIX if conj else ''
        single = len(members) == 1
        name = (str(members[0]) if single
                else ParticleModel._container_name(members)) + suffix
        phase = None
        if single and phases and (members[0].selfConjugate or conj):
            phase = phases.get(members[0])
        return name, (None if single else mixing_matrix.standard_symbol), row, row_im, phase

    def render(name, matrix, row, row_im=None, phase=None, conjugate=False, sign=1):
        if conjugate:
            name = toggle_dagger(name)
        minus = '-' if sign == -1 else ''
        if matrix is None:
            factor = '' if phase is None else (f'Conjugate[{phase}] ' if conjugate else f'{phase} ')
            return minus + factor + _indexed(name, shared)
        if row_im is not None:
            i = '-' if conjugate else '+'
            element = f'({matrix}[{row}, j] {i} I*{matrix}[{row_im}, j])/Sqrt[2]'
        else:
            element = f'{matrix}[{row}, j]'
            if conjugate:
                element = f'Conjugate[{element}]'
        return f'{minus}Module[{{j}}, {element} {_indexed(name, shared + ["j"])}]'

    if conjugate:
        conj, defs = multiplet.conjugateMultiplet, []
        for field in conj.fields:
            lhs = [f'{i}_' for i in shared]
            if multiplet.isospin.dim > 1:
                lhs.append(str(getSU2Index(field)))
            rhs = render(*definition(field.conjugateField), conjugate=True,
                         sign=conj.phase(field.i3))
            defs.append(f'{_indexed(symbol, lhs)} :> {rhs}')
        return '{' + ', '.join(defs) + '}'

    real = getattr(multiplet, 'is_real_nontrivial', False)
    defs = []
    for field in multiplet.fields:
        lhs = [f'{i}_' for i in shared]
        if multiplet.isospin.dim > 1:
            lhs.append(str(getSU2Index(field)))
        if real and field.i3 < 0:
            m = -field.i3
            rhs = render(*definition(multiplet[m]), conjugate=True,
                         sign=su2_reality_phase(multiplet.isospin, m))
        else:
            rhs = render(*definition(field))
        defs.append(f'{_indexed(multiplet.fr_symbol(), lhs)} :> {rhs}')

    return '{' + ', '.join(defs) + '}'
