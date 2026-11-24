'''Command line interface: uvbuilder build | fill.'''

import argparse
import logging
import runpy
import sys
from pathlib import Path

USAGE = '''
  uvbuilder build PROTOMODEL      FeynRules models (.fr and .pkl) of a protomodel
  uvbuilder fill MODEL.pkl CARD   masses, mixing matrices and DECAY Auto from the inputs
'''


def _protomodel(path: Path) -> dict:
    '''ModelBuilder arguments of a protomodel file: a Python file defining `protomodel`
    (a dictionary with masses, decays and labels, see FinalStateModel.fromDict) and
    optionally `possible_dm_reps`, `dm_pdg`, `labeler` and `name`.'''
    if path.suffix != '.py':
        raise ValueError(f'{path}: the protomodel is a Python file defining `protomodel`')
    namespace = runpy.run_path(str(path))
    if not isinstance(namespace.get('protomodel'), dict):
        raise ValueError(f'{path} does not define the dictionary `protomodel`')
    options = dict(modelData=namespace['protomodel'])
    options.update({key: namespace[key] for key in ('possible_dm_reps', 'dm_pdg')
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
    p.add_argument('protomodel', type=Path, help='Python file defining the dictionary `protomodel`')
    p.add_argument('--name', help='model name (default: the protomodel file name)')
    p.add_argument('--out', default='.', help='output folder (default: current)')
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
