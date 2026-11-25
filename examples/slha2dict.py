#!/usr/bin/env python3
'''SLHA file (spectrum and decay tables, e.g. from the SModelS electroweakino scans) to a
protomodel .dict for `uvbuilder build`.

    python slha2dict.py point.slha [-o point.dict] [--max-mass 3000]

Kept: the BSM states (|PDG| > 1000000) lighter than --max-mass; the lightest neutral
one, 1000022, is the DM. Masses are taken as |m| (SLHA neutralino masses carry signs).
Decays of the kept states to kept states are converted to two-body decays:
  - two-body channels are kept as they are;
  - an off-shell three-body decay X -> Y f f' of colorless states, with f f' an SM
    fermion pair, becomes X -> Y V with the virtual boson V that carries the charge
    Q(X) - Q(Y) of the MSSM states (PDG convention: 1000024 and 1000037 are the
    positive charginos): a W+- for an SU(2) doublet pair (u d~, e- nu_e~, ...), a Z
    for a fermion and its antifermion (nu nu~, l+ l-, q q~; a virtual h or photon
    would give the same final states, Z is taken). The charge of the pair itself is
    not used: SOFTSUSY writes e.g. chi_20 -> chi_1- u d~ with the PDG 1000024;
  - any other channel is dropped (reported on stderr).
Branching ratios of channels that map to the same two-body decay are added.
'''

import argparse
import sys
from pathlib import Path

# electric charge (in units of e) of the SM fermions; antiparticles have the opposite sign
CHARGE = {1: -1/3, 3: -1/3, 5: -1/3, 2: 2/3, 4: 2/3, 6: 2/3,
          11: -1, 13: -1, 15: -1, 12: 0, 14: 0, 16: 0}
# SU(2) doublet partners (up-type, down-type) of each generation
DOUBLETS = {frozenset(p) for p in ((2, 1), (4, 3), (6, 5), (12, 11), (14, 13), (16, 15))}
# electric charge of the colorless MSSM states (particles; antiparticles opposite)
BSM_CHARGE = {**{p: 0 for p in (1000022, 1000023, 1000025, 1000035, 1000012, 1000014, 1000016)},
              **{p: 1 for p in (1000024, 1000037)},
              **{p: -1 for p in (1000011, 1000013, 1000015, 2000011, 2000013, 2000015)}}
LABELS = {1000022: 'N1', 1000023: 'N2', 1000025: 'N3', 1000035: 'N4',
          1000024: 'C1', 1000037: 'C2'}


def read_slha(path):
    '''({pdg: mass}, {pdg: [(br, daughters)]}) of the MASS block and DECAY tables.'''
    masses, decays, block, mother = {}, {}, None, None
    for line in Path(path).read_text().splitlines():
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        head = fields[0].upper()
        if head in ('BLOCK', 'XSECTION'):
            block, mother = (fields[1].upper() if head == 'BLOCK' else None), None
        elif head == 'DECAY':
            block, mother = 'DECAY', int(fields[1])
            decays[mother] = []
        elif block == 'MASS' and len(fields) >= 2:
            masses[int(fields[0])] = abs(float(fields[1]))
        elif block == 'DECAY' and len(fields) >= 3:
            n = int(fields[1])
            decays[mother].append((float(fields[0]), tuple(int(x) for x in fields[2:2 + n])))
    return masses, decays


def bsm_charge(pdg):
    return BSM_CHARGE[abs(pdg)] * (1 if pdg > 0 else -1)


def virtual_boson(mother, bsm, sm):
    '''PDG of the boson of an off-shell X -> Y f f' decay, or None.'''
    if abs(mother) not in BSM_CHARGE or abs(bsm) not in BSM_CHARGE or len(sm) != 2:
        return None
    if any(abs(p) not in CHARGE for p in sm):
        return None
    q = bsm_charge(mother) - bsm_charge(bsm)
    if q == 0 and sm[0] == -sm[1]:
        return 23
    if abs(q) == 1 and frozenset(abs(p) for p in sm) in DOUBLETS and sm[0] * sm[1] < 0:
        return 24 * q
    return None


def protomodel(path, max_mass=3000.0):
    masses, decays = read_slha(path)
    kept = {p for p, m in masses.items() if p > 1000000 and m < max_mass}
    out = {'masses': {p: masses[p] for p in sorted(kept)},
           'decays': {p: {} for p in sorted(kept)},
           'labels': {p: LABELS[p] for p in sorted(kept) if p in LABELS}}
    for mother in sorted(kept):
        for br, daughters in decays.get(mother, []):
            bsm = [d for d in daughters if abs(d) > 1000000]
            sm = [d for d in daughters if abs(d) <= 1000000]
            if len(bsm) != 1 or abs(bsm[0]) not in kept:
                channel = None
            elif len(daughters) == 2:
                channel = tuple(daughters)
            else:
                boson = virtual_boson(mother, bsm[0], sm)
                channel = (bsm[0], boson) if boson else None
            if channel is None:
                print(f'{Path(path).name}: dropped {mother} -> {daughters} (BR {br:.3g})',
                      file=sys.stderr)
                continue
            out['decays'][mother][channel] = out['decays'][mother].get(channel, 0.0) + br
    return out


def write_dict(data, path):
    lines = ['{', '    "masses": {']
    lines += [f'        {p}: {m!r},' for p, m in data['masses'].items()]
    lines += ['    },', '    "decays": {']
    for p, channels in data['decays'].items():
        lines += [f'        {p}: {{'] + [f'            {c}: {br!r},' for c, br in channels.items()]
        lines += ['        },']
    lines += ['    },', '    "labels": {']
    lines += [f'        {p}: {label!r},' for p, label in data['labels'].items()]
    lines += ['    }', '}']
    Path(path).write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('slha', type=Path)
    parser.add_argument('-o', '--output', type=Path, help='default: the SLHA name with .dict')
    parser.add_argument('--max-mass', type=float, default=3000.0,
                        help='heavier states are decoupled (default: 3000 GeV)')
    args = parser.parse_args()
    output = args.output or args.slha.with_suffix('.dict')
    write_dict(protomodel(args.slha, args.max_mass), output)
    print(output)


if __name__ == '__main__':
    main()
