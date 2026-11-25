'''SLHA builder: the input parameters of a FeynRules model, the mass matrices of its
mixing blocks, their spectrum and the MG5 parameter card.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
'''

import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

from uvbuilder.tools.fr_utils import SECTORS, FRParameter
from uvbuilder.tools.spectrum import (diagonalize_complex_scalar, diagonalize_dirac,
                            diagonalize_real_scalar, takagi)
from uvbuilder.tools.reps import WeylSpinor

# Sectors whose External parameters are Lagrangian inputs (mixing matrices,
# phases and CG tensors are derived)
INPUT_SECTORS = {key for key, _, block in SECTORS if block is not None}


def input_parameters(fr_model) -> List[FRParameter]:
    '''External Lagrangian parameters of the model, in M$Parameters order
    (real and imaginary parts of a complex coupling as separate entries).'''
    return [p for p in fr_model.parameters()
            if p.kind == 'External' and p.sector in INPUT_SECTORS]


def _number(text):
    '''float of text (Fortran exponents d, D accepted), or None.'''
    try:
        return float(text.replace('d', 'e').replace('D', 'E'))
    except ValueError:
        return None


def _finite(x) -> bool:
    try:
        return math.isfinite(x)
    except TypeError:
        return False


def _symbol(field) -> sp.Symbol:
    return sp.Symbol(field.label)


def _mass_potential(fr_model) -> sp.Expr:
    '''V = sum of the two-field terms at the vev: for scalars all terms of -L
    (Operator.lagrangian_terms), for fermions the operators without their h.c.
    (-L = V + h.c.); components -m of real multiplets are written as conjugates of
    components m; parameters carry the output names.'''
    V = 0
    for op in fr_model.multiplet_model.lagrangian_operators:
        scalar = not any(isinstance(m.rep.spin, WeylSpinor) for m in op.multiplets)
        terms = op.lagrangian_terms(at_vev=True) if scalar else op.expand_singlets(at_vev=True)
        V += sum(t.to_sympy_expr() for t in terms if len(t.fields) == 2)
    real_components = {}
    for m in fr_model.multiplet_model.multiplets:
        if m.is_real_nontrivial:
            for i3, f in m.components.items():
                if i3 > 0:
                    partner, phase = m.conjugate_component(i3)
                    real_components[_symbol(partner)] = phase * _symbol(f.conjugateField)
    rename = {sp.Symbol(old): sp.Symbol(new)
              for old, new in fr_model.renames['parameter'].items()}
    return sp.expand(V.xreplace(real_components).xreplace(rename))


def _declared(fr_model) -> set:
    '''Fields of the non-conjugated multiplets.'''
    return {f for m in fr_model.multiplet_model.non_conjugated_multiplets for f in m.fields}


def mass_matrices(fr_model) -> Dict[str, Tuple[tuple, sp.Matrix]]:
    '''Mass matrix of each mixing block, {U block symbol: (blocks, matrix)}, with blocks
    (U,) or, for a Dirac pair, (U, V), from V = sum of the two-field terms at the vev
    (MixingMatrix.mass_matrix). Real blocks are written in their real basis
    (MixingMatrix.real_basis). Raises if V has a term outside the blocks.'''
    declared = _declared(fr_model)
    matrices = fr_model.mixings.get_mixing_matrices()
    rights = {m.particles: m for m in matrices if m.is_right_block}
    V = _mass_potential(fr_model)
    for m in matrices:
        V = m.to_real_basis(V, declared)
    V = sp.expand(V)
    result, quadratic = {}, 0
    for m in matrices:
        if m.is_right_block:
            continue
        right = rights.get(m.particles) if m.type == 'DiracSpinor' else None
        M, form = m.mass_matrix(V, declared, right)
        result[m.standard_symbol] = ((m,) if right is None else (m, right), M)
        quadratic += form
    rest = sp.simplify(sp.expand(V - quadratic))
    if rest != 0:
        raise ValueError(f'mass terms outside the mixing blocks: {rest}')
    return dict(sorted(result.items()))


# Higgs vev of SMMod.fr (vev = 2 MW sw / ee), numerical value from FeynRules
VEV = 246.220569073486


def numeric_mass_matrices(fr_model):
    '''evaluate(values) -> {U block symbol: numpy matrix}, with values {name: value} of
    exactly the model's inputs (input_parameters), e.g. from a param_card; raises
    ValueError listing missing and unknown names and values that are not finite numbers.
    A complex coupling y takes R<y> + i I<y>; v = VEV. Each matrix of mass_matrices is
    lambdified once (conjugate -> numpy.conj).'''
    matrices = mass_matrices(fr_model)
    entries = input_parameters(fr_model)
    complex_names = {p.base for p in entries if p.part}
    names = sorted({p.base for p in entries})
    symbols = [sp.Symbol(n) for n in names] + [sp.Symbol('v')]
    unknown = set().union(*(M.free_symbols for _, M in matrices.values())) - set(symbols)
    if unknown:
        raise ValueError(f'mass matrix symbols that are not inputs: {sorted(map(str, unknown))}')
    functions = {key: sp.lambdify(symbols, M, 'numpy') for key, (_, M) in matrices.items()}
    expected = {p.name for p in entries}

    def evaluate(values):
        missing = [p.name for p in entries if p.name not in values]
        unknown = sorted(set(values) - expected)
        bad = sorted(n for n in expected & set(values) if not _finite(values[n]))
        if missing or unknown or bad:
            raise ValueError('input values:\n  ' + '\n  '.join(
                [f'missing parameter {n}' for n in missing] +
                [f'unknown parameter {n}' for n in unknown] +
                [f'{n} = {values[n]!r} is not a finite number' for n in bad]))
        args = [complex(values[f'R{n}'], values[f'I{n}']) if n in complex_names else values[n]
                for n in names] + [VEV]
        return {key: np.array(f(*args), dtype=complex) for key, f in functions.items()}
    return evaluate


_DIAGONALIZE = {'RealScalar': diagonalize_real_scalar, 'ComplexScalar': diagonalize_complex_scalar,
                'MajoranaSpinor': takagi, 'DiracSpinor': diagonalize_dirac}


def numeric_spectrum(fr_model):
    '''compute(values) -> (masses, mixings) for the input values (as in evaluate):
    masses   {pdg: mass} of every BSM state of the mixing blocks;
    mixings  {matrix symbol: numpy matrix} of the blocks with more than one state (the
             mixing matrices of the FeynRules model): U, and V for Dirac pairs.
    Column j of a block is its j-th lightest state, the particle with family index j.
    Blocks with one state give only the mass; the sign of their Lagrangian mass is the
    phase P<particle> of the FeynRules model.'''
    evaluate = numeric_mass_matrices(fr_model)
    blocks = {key: b for key, (b, _) in mass_matrices(fr_model).items()}

    def compute(values):
        masses, mixings = {}, {}
        for key, matrix in evaluate(values).items():
            block = blocks[key]
            m, *rotations = _DIAGONALIZE[block[0].type](matrix, block=key)
            masses.update({p.pdg: float(mass) for p, mass in zip(block[0].particles, m)})
            if block[0].dim > 1:
                mixings.update({b.standard_symbol: r for b, r in zip(block, rotations)})
        return masses, mixings
    return compute


def _set_value(line, value) -> str:
    '''line with the last field before # replaced by value; layout and comment kept.'''
    code, sep, comment = line.partition('#')
    return re.sub(r'\S+(\s*)$', lambda m: value + m.group(1), code, count=1) + sep + comment


_ROWS = '# rows (gauge): '


def write_card(fr_model, card, path=None, values=None) -> Path:
    '''Fill a param_card of this model's UFO from its input blocks (BSMMASS, BSMYUKAWA,
    IMBSMYUKAWA, ...), or from values ({name: value}, then written into them):
      MASS entries of the BSM states, mixing blocks <NAME> and IM<NAME> (real blocks:
      <NAME>, each headed by a comment naming its rows and columns)  <- numeric_spectrum
      DECAY of the BSM states present in the card                     <- Auto (BR lines dropped)
    Written to path, by default over card. Everything else is kept; block names are
    compared without case. Raises if an entry is missing or its name (the comment MG5
    writes on each entry) is not the parameter it holds here (a card of another model),
    or if Gf (SMINPUTS 2) gives a vev other than VEV.'''
    card = Path(card)
    path = card if path is None else Path(path)
    blocks = {m.standard_symbol: m for m in fr_model.mixings.get_mixing_matrices()}
    pdgs = {q.pdg for m in blocks.values() for q in m.particles}
    inputs = {(p.block.upper(), (p.order_block,)): p.name for p in input_parameters(fr_model)}
    names = dict(inputs)
    label = {q.label: fr_model.renames['particle'].get(q.label, q.label)       # output names
             for m in blocks.values() for q in m.particles}
    names.update({('MASS', (q.pdg,)): f'M{label[q.label]}' for m in blocks.values() for q in m.particles})
    mixing, headers, declared = {}, {}, _declared(fr_model)
    for p in fr_model.parameters():
        if p.sector != 'Mixing' or p.kind != 'External':
            continue
        m = blocks[p.base]
        for i, j in np.ndindex(m.dim, m.dim):
            key = (p.block.upper(), (i + 1, j + 1))
            names[key], mixing[key] = f'{p.name}{i + 1}x{j + 1}', (p.base, p.part, i, j)
        rows = m.real_basis(declared)[0]
        headers[p.block.upper()] = (f'{_ROWS}{" ".join(map(str, rows))}; '
                                    f'columns (mass): {" ".join(label[q.label] for q in m.particles)}')

    lines = card.read_text().splitlines()
    found, wrong, block, gf = {}, [], None, None
    for n, line in enumerate(lines):
        fields = line.split('#', 1)[0].split()
        if fields and fields[0].upper() in ('BLOCK', 'DECAY'):
            block = fields[1].upper() if fields[0].upper() == 'BLOCK' else None
            continue
        try:
            key = (block, tuple(int(f) for f in fields[:-1]))
        except ValueError:
            continue
        if fields and key == ('SMINPUTS', (2,)):
            gf = _number(fields[-1])
        if fields and key in names:
            written = line.partition('#')[2].split()[:1]
            if written != [names[key]]:
                wrong.append(f'{block} {" ".join(map(str, key[1]))} is '
                             f'{written[0] if written else "unnamed"}, expected {names[key]}')
            found[key] = n
    missing = [f'{b} {" ".join(map(str, i))}' for b, i in names if (b, i) not in found]
    problems = ([f'entries missing: {", ".join(missing)}'] if missing else []) + wrong
    if problems:
        raise ValueError(f'{card}:\n  ' + '\n  '.join(problems))
    if gf is None or abs((math.sqrt(2) * gf) ** -0.5 - VEV) > 1e-6 * VEV:
        raise ValueError(f'{card}: Gf (SMINPUTS 2) = {gf} does not give the vev {VEV}')

    given = values is not None
    if not given:
        values = {name: _number(lines[found[key]].split('#', 1)[0].split()[-1])
                  for key, name in inputs.items()}
    masses, mixings = numeric_spectrum(fr_model)(values)
    numbers = {key: values[name] for key, name in inputs.items()} if given else {}
    numbers.update({('MASS', (pdg,)): mass for pdg, mass in masses.items()})
    for key, (base, part, i, j) in mixing.items():
        numbers[key] = (mixings[base].imag if part == 'I' else mixings[base].real)[i, j]
    replace = {found[key]: _set_value(lines[found[key]], f'{x:.16e}') for key, x in numbers.items()}

    out, block = [], None
    for n, line in enumerate(lines):
        fields = line.split('#', 1)[0].split()
        head = fields[0].upper() if fields else ''
        if line.startswith(_ROWS):
            continue
        if head in ('BLOCK', 'DECAY'):
            block = fields[1].upper() if head == 'BLOCK' else None
            if head == 'DECAY' and int(fields[1]) in pdgs:
                block, line = 'DECAY', _set_value(line, 'Auto')
            out.append(line)
            if block in headers:
                out.append(headers[block])
            continue
        if block == 'DECAY' and fields:
            continue
        out.append(replace.get(n, line))
    path.write_text('\n'.join(out) + '\n')
    return path
