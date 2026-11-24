'''UFO of a FeynRules model (FeynRules via wolframscript), its checks and default
param_card (MG5).'''

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from uvbuilder.config import tool_path
from uvbuilder.tools.fr_names import DEFAULT_SM_FILE

_WL_TEMPLATE = '''Unprotect[Commutator, NonCommutativeMultiply, MatrixSymbol];
FR$Parallel = False;
$FeynRulesPath = SetDirectory[{fr_path}];
<< FeynRules`;
LoadModel[{sm_fr}, {model_fr}];
If[{fr_checks},
  Print["@@BEGIN hermiticity"]; Print[ToString[CheckHermiticity[WeylToDirac[LBSM]], InputForm]]; Print["@@END"];
  Print["@@BEGIN kinetic"]; Print[ToString[CheckKineticTermNormalisation[WeylToDirac[LBSM]], InputForm]]; Print["@@END"];
];
SetDirectory[{out_dir}];
WriteUFO[LSM + WeylToDirac[LBSM], Output -> {ufo_name}, AddDecays -> False];
Print["@@WRITEUFO DONE"];
'''


def _wl_string(s):
    '''Wolfram string literal for a path or name.'''
    return json.dumps(str(s))


def create_ufo(fr_model, name, out_dir=None, fr_checks=True, timeout=3600):
    '''Write the UFO of a FeynRules model in a fresh Wolfram kernel.

    fr_model: path to a .fr file, or an object with to_feynrules(), which is
    written to <out_dir>/<name>.fr. The UFO goes to <out_dir>/<name>_UFO
    (default out_dir: the current folder) and is never written over an existing
    folder. The SM is the model's sm_file (for a .fr path, SMMod.fr of the package).
    Lagrangian LSM + WeylToDirac[LBSM]; no decays.py. With fr_checks, the
    outputs of CheckHermiticity and CheckKineticTermNormalisation are returned.
    '''
    out_dir = Path(out_dir or '.').expanduser().resolve()
    ufo = out_dir / f'{name}_UFO'
    if ufo.exists():
        raise FileExistsError(f'{ufo} already exists: remove it or choose another name.')
    out_dir.mkdir(parents=True, exist_ok=True)

    if hasattr(fr_model, 'to_feynrules'):
        fr_file = out_dir / f'{name}.fr'
        text = fr_model.to_feynrules()
        if fr_file.exists() and fr_file.read_text() != text:
            raise FileExistsError(f'{fr_file} exists with different content.')
        fr_file.write_text(text)
        sm_file = fr_model.sm_file
    else:
        fr_file = Path(fr_model).expanduser().resolve()
        if not fr_file.is_file():
            raise FileNotFoundError(fr_file)
        sm_file = DEFAULT_SM_FILE

    script = _WL_TEMPLATE.format(
        fr_path=_wl_string(tool_path('feynrules')),
        sm_fr=_wl_string(Path(sm_file).resolve()),
        model_fr=_wl_string(fr_file),
        fr_checks='True' if fr_checks else 'False',
        out_dir=_wl_string(out_dir),
        ufo_name=_wl_string(ufo.name))

    with tempfile.TemporaryDirectory() as tmp:
        wl = Path(tmp) / 'create_ufo.wl'
        wl.write_text(script)
        run = subprocess.run([str(tool_path('wolframscript')), '-file', str(wl)], cwd=tmp,
                             capture_output=True, text=True, timeout=timeout)
    log = run.stdout + run.stderr

    sections = dict(re.findall(r'@@BEGIN (\w+)\n(.*?)\n@@END', log, re.S))
    result = dict(
        ufo=str(ufo), fr_file=str(fr_file),
        ok=(ufo / 'vertices.py').is_file() and '@@WRITEUFO DONE' in log,
        hermiticity=sections.get('hermiticity'),
        kinetic=sections.get('kinetic'),
        messages=[line for line in log.splitlines() if re.match(r'\S+::\S+', line)],
        log=log)
    if ufo.is_dir():
        (ufo / 'create_ufo.log').write_text(log)

    print(f"{'OK' if result['ok'] else 'FAILED'}: {ufo}")
    for key in ('hermiticity', 'kinetic'):
        if result[key] is not None:
            print(f'  {key}: {result[key][:300]}')
    for message in result['messages']:
        print(f'  message: {message[:300]}')
    if not result['ok']:
        print(log[-3000:])
    return result


SM_PDGS = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25}

# Particle lists and vertices of a UFO, read in a separate interpreter so that
# its modules (particles, parameters, ...) do not enter this kernel.
_UFO_DUMP = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
import particles, vertices

def number(x):
    try:
        return float(x)
    except Exception:
        return None

P = [dict(pdg=p.pdg_code, name=p.name, antiname=p.antiname, spin=p.spin, color=p.color,
          charge=number(p.charge), ghost=getattr(p, 'GhostNumber', 0) != 0,
          goldstone=bool(getattr(p, 'goldstone', False)))
     for p in particles.all_particles]
V = [dict(name=v.name, pdgs=[p.pdg_code for p in v.particles], color=list(v.color))
     for v in vertices.all_vertices]
print(json.dumps(dict(particles=P, vertices=V)))
'''

# MG5 colour objects: required colour rep of each external (positive) index.
_COLOUR_SLOTS = {'T': (8, 3, -3), 'T6': (8, 6, -6), 'f': (8, 8, 8), 'd': (8, 8, 8),
                 'Epsilon': (3, 3, 3), 'EpsilonBar': (-3, -3, -3),
                 'K6': (6, -3, -3), 'K6Bar': (-6, 3, 3)}
_TRIALITY = {1: 0, 8: 0, 3: 1, -3: 2, 6: 2, -6: 1}
_ANSI = re.compile(r'\x1b\[[0-9;]*m')


def _ufo_info(ufo):
    run = subprocess.run([sys.executable, '-I', '-c', _UFO_DUMP, str(ufo)],
                         capture_output=True, text=True, timeout=300)
    if run.returncode != 0:
        raise RuntimeError(f'cannot import the UFO {ufo}:\n{run.stderr[-2000:]}')
    info = json.loads(run.stdout.strip().splitlines()[-1])
    info['by_pdg'] = {p['pdg']: p for p in info['particles']}
    return info


def _is_bsm(p):
    return abs(p['pdg']) not in SM_PDGS and not p['ghost'] and not p['goldstone']


def _conjugate(p):
    return p['pdg'] if p['name'] == p['antiname'] else -p['pdg']


def _colour_issues(colour, reps):
    issues = []
    for obj, args in re.findall(r'(\w+)\(([^)]*)\)', colour):
        if obj == 'Identity':
            i, j = (int(a) for a in args.split(','))
            if i > 0 and j > 0 and not (reps[i - 1] == -reps[j - 1] or reps[i - 1] == reps[j - 1] == 8):
                issues.append(f'{obj}({args}) on reps {reps[i - 1]}, {reps[j - 1]}')
            continue
        if obj not in _COLOUR_SLOTS:
            issues.append(f'unknown colour object {obj}({args})')
            continue
        for slot, a in zip(_COLOUR_SLOTS[obj], (int(a) for a in args.split(','))):
            if a > 0 and reps[a - 1] != slot:
                issues.append(f'{obj}({args}): index {a} has rep {reps[a - 1]}, expected {slot}')
    return issues


def _static_checks(ufo, info):
    P = info['by_pdg']
    keys = {tuple(sorted(v['pdgs'])) for v in info['vertices']}
    issues = []
    for v in info['vertices']:
        ps = [P[pdg] for pdg in v['pdgs']]
        label = f"{v['name']} [{', '.join(p['name'] for p in ps)}]"
        if any(p['charge'] is None for p in ps) or abs(sum(p['charge'] for p in ps)) > 1e-6:
            issues.append(f'{label}: electric charge not conserved')
        if sum(_is_bsm(p) for p in ps) % 2:
            issues.append(f'{label}: odd number of BSM fields (Z2)')
        reps = [p['color'] for p in ps]
        if sum(_TRIALITY.get(r, 0) for r in reps) % 3:
            issues.append(f'{label}: colour reps {reps} cannot form a singlet')
        for colour in v['color']:
            issues += [f'{label}: {i}' for i in _colour_issues(colour, reps)]
        if not any(p['ghost'] for p in ps) and tuple(sorted(_conjugate(p) for p in ps)) not in keys:
            issues.append(f'{label}: conjugate vertex missing')
    decays = Path(ufo) / 'decays.py'
    decays_py = None
    if decays.is_file():
        decays_py = 'unsummed colour indices' if 'Index(' in decays.read_text() else 'present'
    return dict(issues=issues, n_vertices=len(info['vertices']), decays_py=decays_py)


def _auto_processes(info):
    '''Pair production of every BSM state, and q g > B1 B2 for every q-B1-B2 vertex
    (one quark flavour per BSM pair).'''
    P = info['by_pdg']
    bsm = sorted({abs(p['pdg']) for p in info['particles'] if _is_bsm(p)})
    procs = []
    for a in bsm:
        anti = _conjugate(P[a])
        if P[a]['color'] != 1:
            procs.append(f'21 21 > {a} {anti}')
        procs.append(f'2 -2 > {a} {anti}')
    seen = set()
    for v in info['vertices']:
        ps = [P[pdg] for pdg in v['pdgs']]
        b = [p['pdg'] for p in ps if _is_bsm(p)]
        s = [p for p in ps if not _is_bsm(p)]
        if len(ps) != 3 or len(b) != 2 or s[0]['spin'] != 2 or s[0]['color'] == 1:
            continue
        key = tuple(sorted(b))
        if key not in seen:
            seen.add(key)
            procs.append(f"{_conjugate(s[0])} 21 > {b[0]} {b[1]}")
    return procs


def _mg5(commands, model, workdir, timeout, load=True):
    script = Path(workdir) / 'run.mg5'
    header = ['set automatic_html_opening False'] + ([f'import model {model}'] if load else [])
    script.write_text('\n'.join(header + list(commands)) + '\n')
    run = subprocess.run([str(tool_path('mg5')), str(script)], cwd=workdir, capture_output=True,
                         text=True, timeout=timeout)
    return _ANSI.sub('', run.stdout + run.stderr)


_GAUGE_ROW = re.compile(r'^(\S.*?\S)\s+((?:[-+]?\d\.\d+e[-+]\d+\s+)+)(Passed|Failed)\s*$', re.M)


def _parse_gauge(log, procs):
    chunks = re.split(r'^check gauge (.+)$', log, flags=re.M)
    found = dict(zip(chunks[1::2], chunks[2::2]))
    out = {}
    for proc in procs:
        chunk = found.get(proc)
        if chunk is None:
            continue
        row = _GAUGE_ROW.search(chunk)
        if row:
            out[proc] = dict(mg_process=row.group(1), status=row.group(3),
                             reldiff=float(row.group(2).split()[-1]), detail='')
        elif 'No diagrams' in chunk or 'No amplitudes' in chunk:
            out[proc] = dict(mg_process='', status='no diagrams', reldiff=None, detail='')
        else:
            error = [line for line in chunk.splitlines() if 'rror' in line]
            out[proc] = dict(mg_process='', status='error', reldiff=None,
                             detail=(' | '.join(error[:2]) or chunk.strip()[-300:])[:200])
    return out


def _read_decays(path, pdgs):
    out, current = {}, None
    for line in Path(path).read_text().splitlines():
        t = line.split('#')[0].split()
        if not t:
            continue
        if t[0].upper() == 'DECAY':
            current = int(t[1])
            out[current] = dict(width=float(t[2]), br={})
        elif t[0].upper() == 'BLOCK':
            current = None
        elif current is not None and len(t) >= 3:
            n = int(t[1])
            out[current]['br'][tuple(int(x) for x in t[2:2 + n])] = float(t[0])
    return {pdg: out[pdg] for pdg in pdgs if pdg in out}


def test_ufo(ufo, processes='auto', gauge=True, widths=True, timeout=1800, verbose=True):
    '''Checks of a UFO; failures are reported, not raised.

    static : per vertex, charge, Z2 (even number of BSM fields), colour reps
             (singlet triality and the slots of T, T6, K6, K6Bar, Epsilon, ...),
             presence of the conjugate vertex; state of decays.py.
    gauge  : MG5 check gauge on `processes` (PDG codes; 'auto': pair production
             of each BSM state and q g > B1 B2 for each q-B1-B2 vertex).
    widths : MG5 compute_widths of all BSM states, two-body, {pdg: {width, br}}.
    Everything runs on a temporary copy of the UFO without decays.py (converted to
    Python 3 by MG5; widths computed numerically), so the folder is left unchanged.
    '''
    ufo = Path(ufo).expanduser().resolve()
    with tempfile.TemporaryDirectory() as tmp:
        model = Path(tmp) / ufo.name
        shutil.copytree(ufo, model, ignore=shutil.ignore_patterns('*.pkl', '__pycache__', 'decays.py'))
        convert_log = _mg5([f'convert model {model}'], model, tmp, timeout, load=False)
        info = _ufo_info(model)
        report = dict(ufo=str(ufo), static=_static_checks(ufo, info), gauge={}, widths=None,
                      logs=dict(convert=convert_log))
        procs = _auto_processes(info) if processes == 'auto' else list(processes)
        bsm = sorted({abs(p['pdg']) for p in info['particles'] if _is_bsm(p)})

        if gauge:
            pending = list(procs)
            while pending:      # MG5 stops a script at the first error: rerun the rest
                log = _mg5([f'check gauge {p}' for p in pending], model, tmp, timeout)
                report['logs'].setdefault('gauge', []).append(log)
                done = _parse_gauge(log, pending)
                report['gauge'].update(done)
                if not done:
                    report['gauge'].update({p: dict(mg_process='', status='not run', reldiff=None,
                                                    detail='MG5 stopped') for p in pending})
                    break
                pending = [p for p in pending if p not in report['gauge']]
        if widths and bsm:
            card = Path(tmp) / 'widths.dat'
            log = _mg5([f"compute_widths {' '.join(map(str, bsm))} --body_decay=2 --output={card}"],
                       model, tmp, timeout)
            report['logs']['widths'] = log
            failed = 'Traceback' in log or 'with error' in log or not card.is_file()
            report['widths'] = dict(error=[line.strip()[:200] for line in log.splitlines()
                                           if 'Error' in line][-2:]) if failed else _read_decays(card, bsm)

    if verbose:
        print_ufo_report(report, info)
    return report


def print_ufo_report(report, info=None):
    names = {p['pdg']: p['name'] for p in info['particles']} if info else {}
    static = report['static']
    print(f"UFO: {report['ufo']}")
    print(f"static checks: {static['n_vertices']} vertices, {len(static['issues'])} issues; "
          f"decays.py: {static['decays_py'] or 'absent'}")
    for issue in static['issues']:
        print(f'   {issue}')
    if report['gauge']:
        print('check gauge:')
        for proc, r in report['gauge'].items():
            rd = f"{r['reldiff']:.1e}" if r['reldiff'] is not None else ''
            print(f"   {proc:32s} {r['mg_process']:24s} {rd:>8s}  {r['status']} {r['detail']}")
    if report['widths'] is not None:
        print('widths (GeV) and BRs:')
        if 'error' in report['widths']:
            print(f"   error: {report['widths']['error']}")
        else:
            for pdg, d in report['widths'].items():
                channels = ', '.join(f"{' '.join(names.get(x, str(x)) for x in ch)} {br:.3f}"
                                     for ch, br in sorted(d['br'].items(), key=lambda kv: -kv[1]))
                print(f"   {names.get(pdg, pdg)} ({pdg}): {d['width']:.6g}  [{channels}]")


_WRITE_CARD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import write_param_card
write_param_card.ParamCardWriter(sys.argv[2], generic=True)
'''


def ufo_param_card(ufo, output='param_card.dat', timeout=600):
    '''Default param_card of a UFO (its write_param_card.py), written to output.
    Runs on a temporary copy converted to Python 3 by MG5, so the folder is left
    unchanged.'''
    ufo = Path(ufo).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    with tempfile.TemporaryDirectory() as tmp:
        model = Path(tmp) / ufo.name
        shutil.copytree(ufo, model, ignore=shutil.ignore_patterns('*.pkl', '__pycache__', 'decays.py'))
        log = _mg5([f'convert model {model}'], model, tmp, timeout, load=False)
        run = subprocess.run([sys.executable, '-I', '-c', _WRITE_CARD, str(model), str(output)],
                             cwd=tmp, capture_output=True, text=True, timeout=timeout)
    if run.returncode != 0 or not output.is_file():
        raise RuntimeError(f'cannot write the param_card of {ufo}:\n{run.stderr[-2000:]}\n{log[-1000:]}')
    return output
