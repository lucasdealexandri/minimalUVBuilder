'''Command line interface: uvbuilder build | ufo | card | fill | check.'''

import argparse
import logging
import runpy
import sys
from pathlib import Path

USAGE = '''
  uvbuilder build PROTOMODEL      FeynRules models (.fr and .pkl) of a protomodel
  uvbuilder ufo MODEL.fr          UFO of a model (FeynRules)
  uvbuilder card MODEL_UFO        default param_card of a UFO (MG5)
  uvbuilder fill MODEL.pkl CARD   masses, mixing matrices and DECAY Auto from the inputs
  uvbuilder check MODEL_UFO       static checks, MG5 check gauge and widths of a UFO
  uvbuilder config [NAME [PATH]]  show or set the paths of FeynRules, wolframscript, MG5
'''


def _protomodel(path: Path) -> dict:
    '''ModelBuilder arguments of a protomodel file: an SLHA file, or a Python file
    defining `protomodel` (a dictionary with masses, decays and labels, or the path of
    an SLHA file) and optionally `possible_dm_reps`, `dm_pdg`, `labels`, `labeler` and
    `name`.'''
    if path.suffix != '.py':
        return dict(modelData=str(path))
    namespace = runpy.run_path(str(path))
    if 'protomodel' not in namespace:
        raise SystemExit(f'{path} does not define `protomodel`')
    data = namespace['protomodel']
    if isinstance(data, (str, Path)):
        data = str((path.parent / data).resolve())
    options = dict(modelData=data)
    options.update({key: namespace[key] for key in ('possible_dm_reps', 'dm_pdg', 'labels')
                    if key in namespace})
    for key in ('labeler', 'name'):
        if key in namespace:
            options[f'_{key}'] = namespace[key]
    return options


def build(args):
    from uvbuilder.modelBuilder.modelBuilder import ModelBuilder
    from uvbuilder.models.multipletModel import default_component_label
    options = _protomodel(args.protomodel)
    labeler = options.pop('_labeler', default_component_label)
    name = args.name or options.pop('_name', None) or args.protomodel.stem
    options.pop('_name', None)
    out = Path(args.out)
    models = ModelBuilder(**options).generateFeynRulesModels(
        outputFile=str(out / f'{name}.fr'), name=name, labeler=labeler)
    count = 0
    for fr_model in models:
        fr_model.save2file()
        pickle = Path(fr_model.save2pickle())
        multiplets = ', '.join(sorted(m.fr_symbol() for m in
                                      fr_model.multiplet_model.non_conjugated_multiplets))
        print(f'{fr_model.output_file}  {pickle.name}  [{multiplets}]')
        count += 1
    if not count:
        print('no valid model (run with -v to see why the candidates were rejected)')
        return 1
    return 0


def ufo(args):
    from uvbuilder.tools.ufo import create_ufo
    result = create_ufo(args.model, args.name or args.model.stem, out_dir=args.out,
                        fr_checks=not args.no_checks)
    return 0 if result['ok'] else 1


def card(args):
    from uvbuilder.tools.ufo import ufo_param_card
    print(ufo_param_card(args.ufo, args.output))
    return 0


def fill(args):
    from uvbuilder.models.feynrulesmodel import FeynRulesModel
    from uvbuilder.paramCard import ParamCardGenerator
    print(ParamCardGenerator(FeynRulesModel.load(args.model)).fill(args.card, output=args.output))
    return 0


def check(args):
    from uvbuilder.tools.ufo import test_ufo
    report = test_ufo(args.ufo, gauge=not args.no_gauge, widths=not args.no_widths)
    failed = report['static']['issues'] or any(
        r['status'] != 'Passed' for r in report['gauge'].values())
    return 1 if failed else 0


def config(args):
    from uvbuilder.config import VARIABLES, config_file, lookup, set_path
    if args.path is not None:
        set_path(args.name, args.path)
    print(f'config file: {config_file()}')
    for name in ([args.name] if args.name else VARIABLES):
        value, source = lookup(name)
        print(f'  {name:14s} {value or "-"}  ({source})')
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog='uvbuilder', usage=USAGE,
                                     description='Minimal UV Lagrangians from protomodels.')
    parser.add_argument('-v', '--verbose', action='store_true', help='debug messages')
    commands = parser.add_subparsers(dest='command', required=True)

    p = commands.add_parser('build', help='FeynRules models of a protomodel')
    p.add_argument('protomodel', type=Path, help='SLHA file, or Python file defining `protomodel`')
    p.add_argument('--name', help='model name (default: the protomodel file name)')
    p.add_argument('--out', default='.', help='output folder (default: current)')
    p.set_defaults(run=build)

    p = commands.add_parser('ufo', help='UFO of a model (FeynRules)')
    p.add_argument('model', type=Path, help='.fr file')
    p.add_argument('--name', help='UFO name, <name>_UFO (default: the .fr file name)')
    p.add_argument('--out', default=None, help='output folder (default: current)')
    p.add_argument('--no-checks', action='store_true', help='skip the FeynRules checks')
    p.set_defaults(run=ufo)

    p = commands.add_parser('card', help='default param_card of a UFO (MG5)')
    p.add_argument('ufo', type=Path, help='UFO folder')
    p.add_argument('-o', '--output', default='param_card.dat', help='default: param_card.dat')
    p.set_defaults(run=card)

    p = commands.add_parser('fill', help='fill a param_card from its inputs')
    p.add_argument('model', type=Path, help='.pkl file written by build')
    p.add_argument('card', type=Path, help='param_card with the inputs set')
    p.add_argument('-o', '--output', default=None, help='default: over the card')
    p.set_defaults(run=fill)

    p = commands.add_parser('check', help='checks of a UFO (MG5)')
    p.add_argument('ufo', type=Path, help='UFO folder')
    p.add_argument('--no-gauge', action='store_true', help='skip MG5 check gauge')
    p.add_argument('--no-widths', action='store_true', help='skip MG5 compute_widths')
    p.set_defaults(run=check)

    p = commands.add_parser('config', help='show or set the paths of the external programs')
    p.add_argument('name', nargs='?', choices=['feynrules', 'wolframscript', 'mg5'])
    p.add_argument('path', nargs='?', help='new path (written into the config file)')
    p.set_defaults(run=config)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format='%(levelname)s %(name)s: %(message)s')
    try:
        return args.run(args)
    except (RuntimeError, FileNotFoundError, FileExistsError, ValueError, TypeError) as error:
        print(f'uvbuilder {args.command}: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
