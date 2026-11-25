'''Command line interface: uvbuilder build | fill.'''

import argparse
import ast
import logging
import re
import runpy
import sys
from pathlib import Path

USAGE = '''
  uvbuilder build PROTOMODEL [--out FOLDER] [--dm-spin S] [--dm-charges C]
                            [--name N] [--author A] [--version V]
                                  valid models of a protomodel (.fr, .pkl and .md with --out)
  uvbuilder fill MODEL.pkl CARD   masses, mixing matrices and DECAY Auto from the inputs
'''


def _protomodel(path: Path) -> dict:
    '''ModelBuilder arguments of a protomodel file (a dictionary with masses, decays and
    labels, see FinalStateModel.fromDict): a .dict file holding the dictionary (as
    written by the SModelS protomodel builder; read as a literal, never executed), or a
    Python file defining `protomodel` and optionally `possible_dm_reps`, `dm_pdg`,
    `labeler`, `name`, `author` and `version`.'''
    if path.suffix == '.dict':
        try:
            data = ast.literal_eval(path.read_text())
        except (SyntaxError, ValueError) as error:
            raise ValueError(f'{path} is not a dictionary literal: {error}') from None
        if not isinstance(data, dict):
            raise ValueError(f'{path} is not a dictionary literal')
        return dict(modelData=data)
    if path.suffix != '.py':
        raise ValueError(f'{path}: the protomodel is a .dict file or a Python file')
    namespace = runpy.run_path(str(path))
    if not isinstance(namespace.get('protomodel'), dict):
        raise ValueError(f'{path} does not define the dictionary `protomodel`')
    options = dict(modelData=namespace['protomodel'])
    options.update({key: namespace[key] for key in ('possible_dm_reps', 'dm_pdg')
                    if key in namespace})
    for key in ('labeler', 'name', 'author', 'version'):
        if key in namespace:
            options[f'_{key}'] = namespace[key]
    return options


SPINS = {'scalar': 'SCALAR', 'fermion': 'DIRAC_SPINOR'}


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


def build(args):
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.models.multipletModel import default_component_label
    from uvbuilder.report import model_report
    options = _protomodel(args.protomodel)
    options['possible_dm_reps'] = _dm_reps(args, options.get('possible_dm_reps', {}))
    labeler = options.pop('_labeler', default_component_label)
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
        outputFile=str(out / f'{name}.fr'), name=name, labeler=labeler, **info)
    count = 0
    for fr_model in models:
        multiplets = ', '.join(sorted(m.fr_symbol() for m in
                                      fr_model.multiplet_model.non_conjugated_multiplets))
        if args.out is None:
            print(f'{fr_model.name}  [{multiplets}]')
        else:
            fr_model.save2file()
            pickle = Path(fr_model.save2pickle())
            report = pickle.with_suffix('.md')
            report.write_text(model_report(fr_model))
            print(f'{fr_model.output_file}  {pickle.name}  {report.name}  [{multiplets}]')
        count += 1
    if not count:
        print('no valid model (run with -v to see why the candidates were rejected)')
        return 1
    if args.out is None:
        print('nothing written: use --out FOLDER to write the .fr, .pkl and .md files')
    return 0


def fill(args):
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.paramCard import ParamCardGenerator
    print(ParamCardGenerator(FeynRulesModel.load(args.model)).fill(args.card, output=args.output))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog='uvbuilder', usage=USAGE,
                                     description='Minimal UV Lagrangians from protomodels.')
    parser.add_argument('-v', '--verbose', action='store_true', help='debug messages')
    commands = parser.add_subparsers(dest='command', required=True)

    p = commands.add_parser('build', help='FeynRules models of a protomodel')
    p.add_argument('protomodel', type=Path,
                   help='.dict file, or Python file defining the dictionary `protomodel`')
    p.add_argument('--name', help='model name (default: the protomodel file name)')
    p.add_argument('--author', action='append',
                   help='author of the models (default: UVBuilder); repeat for several')
    p.add_argument('--version', help='version of the models (default: 1.0)')
    p.add_argument('--out', help='folder for the .fr, .pkl and .md files (default: none written)')
    p.add_argument('--dm-spin', choices=['scalar', 'fermion', 'both'],
                   help='spin of the DM candidate (default: both)')
    p.add_argument('--dm-charges', type=_dm_charges, action='append', metavar='CHARGES',
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
