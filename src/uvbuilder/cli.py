"""
.. module:: cli
   :synopsis: Command line interface: uvbuilder build and uvbuilder fill.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""

import argparse
import ast
import logging
import re
import sys
from pathlib import Path

USAGE = '''
  uvbuilder build PROTOMODEL.dict [-o FOLDER] [-s SPIN] [-q CHARGES]
                                 [--name N] [--author A] [--version V]
                                  valid models of a protomodel (<model>/<model>.fr, .pkl, .md with --out)
  uvbuilder fill MODEL.pkl CARD   masses, mixing matrices and DECAY Auto from the inputs
'''


def _protomodel(path: Path) -> dict:
    '''ModelBuilder arguments of a .dict protomodel file: the dictionary of masses, decays
    and labels (see FinalStateModel.fromDict), as written by the SModelS protomodel
    builder, with the optional keys `name`, `author`, `version`, `dm_spin` and
    `dm_charges` (values as on the command line). Read as a literal, never executed.'''
    if path.suffix != '.dict':
        raise ValueError(f'{path}: the protomodel is a .dict file')
    try:
        data = ast.literal_eval(path.read_text())
    except (SyntaxError, ValueError) as error:
        raise ValueError(f'{path} is not a dictionary literal: {error}') from None
    if not isinstance(data, dict):
        raise ValueError(f'{path} is not a dictionary literal')
    return _dict_options(data)


SPINS = {'scalar': 'SCALAR', 'fermion': 'DIRAC_SPINOR'}


def _dict_options(data: dict) -> dict:
    '''ModelBuilder arguments of a .dict protomodel: its options taken out, the rest is
    the protomodel.'''
    from uvbuilder.tools import reps
    data = dict(data)
    options = {f'_{key}': data.pop(key) for key in ('name', 'author', 'version') if key in data}
    dm = {}
    if 'dm_spin' in data:
        spin = data.pop('dm_spin')
        if spin not in (*SPINS, 'both'):
            raise ValueError(f'dm_spin {spin!r}: scalar, fermion or both')
        dm['spin'] = [getattr(reps, SPINS[s]) for s in (SPINS if spin == 'both' else [spin])]
    if 'dm_charges' in data:
        charges = data.pop('dm_charges')
        try:
            dm['globalU1s'] = [_dm_charges(c) for c in ([charges] if isinstance(charges, str) else charges)]
        except argparse.ArgumentTypeError as error:
            raise ValueError(f'dm_charges {error}') from None
    if dm:
        options['possible_dm_reps'] = dm
    return dict(modelData=data, **options)


def _dm_charges(text: str):
    '''GlobalU1s of the DM from 'none' (no charge tracked) or e.g. 'B=0,D=1'.'''
    from uvbuilder.tools.reps import GlobalU1s
    if text.strip().lower() == 'none':
        return GlobalU1s()
    try:
        charges = dict(item.split('=') for item in text.replace(' ', '').split(','))
        return GlobalU1s(**charges)
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError(
            f'{text!r}: none, or charges of B, L, D as in B=0,D=1') from None


def _dm_reps(args, given: dict) -> dict:
    '''possible_dm_reps: spins and global charges of the command line, else those of
    the protomodel file (given), else both spins and no charge.'''
    from uvbuilder.tools import reps
    dm = dict(given)
    if args.dm_spin or 'spin' not in dm:
        spins = list(SPINS) if args.dm_spin in (None, 'both') else [args.dm_spin]
        dm['spin'] = [getattr(reps, SPINS[s]) for s in spins]
    if args.dm_charges or 'globalU1s' not in dm:
        dm['globalU1s'] = args.dm_charges or [_dm_charges('none')]
    return dm


def _columns(rows) -> list:
    '''rows as lines of left-aligned columns.'''
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    return ['  ' + '  '.join(c.ljust(w) for c, w in zip(r, widths)).rstrip() for r in rows]


def _describe(fr_model) -> str:
    '''The model's multiplets (spin, SU(3), SU(2), Y, global charges) and its states (kind
    and the gauge components of their mixing block; the two sides of a Dirac block are
    separated by |).'''
    from uvbuilder.models.slhaBuilder import _declared
    from uvbuilder.report import KIND
    multiplets = [['multiplet', 'spin', 'SU(3)', 'SU(2)', 'Y', 'charges']]
    for m in sorted(fr_model.multiplet_model.non_conjugated_multiplets, key=lambda m: m.fr_symbol()):
        r = m.rep
        multiplets.append([m.fr_symbol(), str(r.spin), str(r.colordim), str(r.isospin.dim),
                           str(r.hypercharge), str(r.globalU1s).strip('()')])
    label = lambda q: fr_model.renames['particle'].get(q.label, q.label)
    declared = _declared(fr_model)
    matrices = fr_model.mixings.get_mixing_matrices()
    rights = {m.particles: m for m in matrices if m.is_right_block}
    states = [['state', 'kind', 'components']]
    for m in sorted((m for m in matrices if not m.is_right_block), key=lambda m: m.standard_symbol):
        if m.type == 'ComplexScalar':
            components = [f.label for f in m.f_components]
        else:
            components = [str(f) for f in m.real_basis(declared)[0]]
        components = ', '.join(components)
        if m.type == 'DiracSpinor':
            components += ' | ' + ', '.join(f.label for f in rights[m.particles].f_components)
        states.append([' '.join(label(q) for q in m.particles), KIND[m.type], components])
    return '\n'.join([fr_model.name] + _columns(multiplets) + _columns(states))


def build(args):
    '''`uvbuilder build`: lists the valid models of the protomodel; with --out, writes each
    to <out>/<model>/ as .fr, .pkl and .md. Exit code 1 if there is no valid model.'''
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.report import model_report
    options = _protomodel(args.protomodel)
    options['possible_dm_reps'] = _dm_reps(args, options.get('possible_dm_reps', {}))
    # the default name is the file name made a valid model name: pmodel-2 -> pmodel_2
    stem = re.sub(r'\W', '_', args.protomodel.stem)
    name = args.name or options.pop('_name', None) or (stem if stem[:1].isalpha() else f'M{stem}')
    options.pop('_name', None)
    # author and version: the command line, else the protomodel file, else the defaults
    # of FeynRulesModel
    info = {key: getattr(args, key) or options.pop(f'_{key}') for key in ('author', 'version')
            if getattr(args, key) or f'_{key}' in options}
    options.pop('_author', None), options.pop('_version', None)
    out = Path(args.out or '.')
    models = ModelBuilder(**options).generateFeynRulesModels(
        outputFile=str(out / f'{name}.fr'), name=name, **info)
    count = 0
    for fr_model in models:
        multiplets = ', '.join(sorted(m.fr_symbol() for m in
                                      fr_model.multiplet_model.non_conjugated_multiplets))
        if args.out is None:
            print(_describe(fr_model) + '\n')
        else:
            # each model in its own folder: <out>/<model>/<model>.fr, .pkl, .md
            folder = out / fr_model.name
            fr_model.output_file = str(folder / Path(fr_model.output_file).name)
            fr_model.save2file()
            pickle = Path(fr_model.save2pickle())
            report = pickle.with_suffix('.md')
            report.write_text(model_report(fr_model))
            print(f'{folder}/  {Path(fr_model.output_file).name}  {pickle.name}  {report.name}  [{multiplets}]')
        count += 1
    if not count:
        print('no valid model (run with -v to see why the candidates were rejected)')
        return 1
    if args.out is None:
        print("nothing written: use -o FOLDER to write each model's .fr, .pkl and .md")
    return 0


def fill(args):
    '''`uvbuilder fill`: writes the masses, mixing matrices and DECAY Auto into the card,
    from its inputs and the pickled model.'''
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.paramCard import ParamCardGenerator
    print(ParamCardGenerator(FeynRulesModel.load(args.model)).fill(args.card, output=args.output))
    return 0


def main(argv=None) -> int:
    '''Entry point of the `uvbuilder` command. Returns the exit code; errors in the input
    are printed, not raised.'''
    parser = argparse.ArgumentParser(prog='uvbuilder', usage=USAGE,
                                     description='Minimal UV Lagrangians from protomodels.')
    parser.add_argument('-v', '--verbose', action='store_true', help='debug messages')
    commands = parser.add_subparsers(dest='command', required=True)

    p = commands.add_parser('build', help='FeynRules models of a protomodel')
    p.add_argument('protomodel', type=Path,
                   help='.dict file with the protomodel')
    p.add_argument('--name', help='model name (default: the protomodel file name)')
    p.add_argument('--author', action='append',
                   help='author of the models (default: UVBuilder); repeat for several')
    p.add_argument('--version', help='version of the models (default: 1.0)')
    p.add_argument('-o', '--out', help='folder for the models, one subfolder each with its .fr, .pkl '
                   'and .md (default: none written)')
    p.add_argument('-s', '--dm-spin', choices=['scalar', 'fermion', 'both'],
                   help='spin of the DM candidate (default: both)')
    p.add_argument('-q', '--dm-charges', type=_dm_charges, action='append', metavar='CHARGES',
                   help='global charges of the DM: none (default), or e.g. B=0,D=1 '
                        '(charges not given are not tracked); repeat for alternatives')
    p.set_defaults(run=build)

    p = commands.add_parser('fill', help='fill a param_card from its inputs')
    p.add_argument('model', type=Path, help='.pkl file written by build')
    p.add_argument('card', type=Path, help='param_card with the inputs set')
    p.add_argument('-o', '--output', default=None, help='default: over the card')
    p.set_defaults(run=fill)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format='%(levelname)s %(name)s: %(message)s')
    try:
        return args.run(args)
    except (FileNotFoundError, ValueError, TypeError) as error:
        print(f'uvbuilder {args.command}: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
