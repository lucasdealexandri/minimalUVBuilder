'''Markdown report of a FeynRules model: its mass eigenstates, input parameters and the
mass matrices of its mixing blocks (LaTeX), to help choose the parameters of a
param_card. The conventions are in docs/conventions.md.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
'''

import re

import sympy as sp

from uvbuilder.models.slhaBuilder import VEV, _declared, input_parameters, mass_matrices

KIND = {'MajoranaSpinor': 'Majorana fermion', 'DiracSpinor': 'Dirac fermion',
        'RealScalar': 'real scalar', 'ComplexScalar': 'complex scalar'}
DIAGONALIZATION = {'MajoranaSpinor': 'Takagi', 'DiracSpinor': 'SVD',
                   'RealScalar': 'orthogonal', 'ComplexScalar': 'unitary'}


def _mathrm(name) -> str:
    '''name in upright LaTeX, underscores escaped.'''
    return r'\mathrm{%s}' % name.replace('_', r'\_')


def _parameter_latex(name) -> str:
    '''MChiI1Y0 -> M_{ChiI1Y0}, muPhi... -> \\mu_{Phi...}: the leading M or lowercase word
    is the symbol, the rest (field labels) its subscript.'''
    head, tail = re.match(r'(M|[a-z]+)?(.*)', name).groups()
    if not head:
        return _mathrm(name)
    head = sp.latex(sp.Symbol(head))
    return f'{head}_{{{_mathrm(tail)}}}' if tail else head


def _card_symbols(fr_model):
    '''({symbol: expression in the param_card entries}, {symbol: LaTeX}): y -> Ry + i Iy
    for complex couplings, every entry and v real.'''
    entries = {}
    for p in input_parameters(fr_model):
        entries.setdefault(p.base, []).append(p)
    subs = {sp.Symbol('v'): sp.Symbol('v', real=True)}
    names = {}
    for base, ps in entries.items():
        names[sp.Symbol(base)] = base_latex = _parameter_latex(base)
        card = {p.part: sp.Symbol(p.name, real=True) for p in ps}
        for part, s in card.items():
            names[s] = {'R': r'\mathrm{Re}\,', 'I': r'\mathrm{Im}\,', '': ''}[part] + base_latex
        subs[sp.Symbol(base)] = card['R'] + sp.I * card['I'] if 'R' in card else card['']
    return subs, names


def _matrix_latex(M, rows, columns, names) -> list:
    '''M as a LaTeX array, the rows and columns labeled by the gauge components.'''
    entry = lambda x: sp.latex(x, symbol_names=names)
    lines = [r'\begin{array}{c|' + 'c' * len(columns) + '}',
             ' & ' + ' & '.join(_mathrm(str(c)) for c in columns) + r' \\ \hline']
    lines += [_mathrm(str(r)) + ' & ' + ' & '.join(entry(M[i, j]) for j in range(M.shape[1])) + r' \\'
              for i, r in enumerate(rows)]
    # a math block, not $$: GitHub applies the Markdown escapes (\\ -> \) inside $$
    return ['```math'] + lines + [r'\end{array}', '```']


def _table(rows, header) -> list:
    '''Markdown table.'''
    line = lambda r: '| ' + ' | '.join(r) + ' |'
    return [line(header), line(['---'] * len(header))] + [line(r) for r in rows]


def _parameter_type(p) -> str:
    '''Type of an input in the report: mass, M² (a scalar mass, entering squared), Re or Im
    (part of a complex coupling), or coupling.'''
    if p.part:
        return 'Re' if p.part == 'R' else 'Im'
    if p.sector == 'Mass':
        return 'M²' if p.description.strip('"').startswith('Mass of') else 'mass'
    return 'coupling'


def model_report(fr_model) -> str:
    '''Markdown report of the model: mass eigenstates, input parameters (param_card entry,
    type, the mass matrices they enter), mass matrices in LaTeX, multiplets and notes.'''
    mm = fr_model.multiplet_model
    label = lambda q: fr_model.renames['particle'].get(q.label, q.label)
    # where a coupling enters with its conjugate (real scalar blocks), the matrix is written
    # in the param_card entries Ry, Iy; elsewhere with y = Ry + i Iy
    card, names = _card_symbols(fr_model)
    blocks = {key: (b, sp.simplify(M.xreplace(card)) if M.has(sp.conjugate) else M)
              for key, (b, M) in mass_matrices(fr_model).items()}
    declared = _declared(fr_model)
    short = lambda key: key[4:] if key.startswith('UMix') else key
    out = [f'# {fr_model.name}', '', 'Inputs for `uvbuilder fill`. Conventions: '
           '`docs/conventions.md` (uvbuilder repository).', '', '## Mass eigenstates', '']
    kind = {q: KIND[b.type] for (b, *_), _ in blocks.values() for q in b.particles}
    states = sorted(kind, key=lambda q: q.pdg)
    out += _table([[label(q), str(q.pdg), kind[q], str(q.eCharge), f'{q.mass}'] for q in states],
                  ['state', 'PDG', 'kind', 'charge', 'mass (protomodel)']) + ['']
    # which mass matrices each parameter enters
    entering = {}
    for key, (_, M) in blocks.items():
        for s in M.free_symbols:
            entering.setdefault(str(s), set()).add(short(key))
    rows = [[p.name, f'{p.block.upper()} {p.order_block}', _parameter_type(p),
             ', '.join(sorted(entering.get(p.name, set()) | entering.get(p.base, set()))) or '-']
            for p in input_parameters(fr_model)]
    out += ['## Input parameters', ''] + _table(rows, ['name', 'card', 'type', 'enters']) + ['']
    out += ['## Mass matrices', '']
    for key, (blocks_, M) in blocks.items():
        block = blocks_[0]
        rows_ = block.real_basis(declared)[0]
        if block.type == 'DiracSpinor':
            columns = [f.label for f in blocks_[1].f_components]
        elif block.type == 'ComplexScalar':
            rows_ = [f.conjugateField.label for f in block.f_components]
            columns = [f.label for f in block.f_components]
        else:
            columns = rows_
        how = DIAGONALIZATION[block.type] if block.dim > 1 else 'single state'
        states = ', '.join(label(q) for q in block.particles)
        out += [f'### {key}: {KIND[block.type]}, {how} → {states}', '']
        out += _matrix_latex(M, rows_, columns, names) + ['']
    out += ['## Multiplets', '']
    out += _table([[m.fr_symbol(), str(m.rep.spin), str(m.rep.isospin.dim), str(m.rep.hypercharge),
                    str(m.rep.colordim), str(m.rep.globalU1s)] for m in mm.non_conjugated_multiplets],
                  ['multiplet', 'spin', 'SU(2)', 'Y', 'color', 'global charges']) + ['']
    notes = [f'- v = {VEV} GeV (G_F).']
    if any(p.part for p in input_parameters(fr_model)):
        notes.append('- Complex couplings: y = Re y + i Im y, card entries `Ry...`, `Iy...`.')
    if any(b.type == 'RealScalar' for (b, *_), _ in blocks.values()):
        notes.append('- Real scalar blocks are in the real basis Phi = (Phi_R + i Phi_I)/sqrt(2).')
    return '\n'.join(out + ['## Notes', ''] + notes) + '\n'
