"""SLHA builder input template (models/slhaBuilder.py)."""

import pytest

import sympy as sp

import numpy as np

from uvbuilder.models.slhaBuilder import (VEV, input_parameters, mass_matrices, numeric_mass_matrices,
                                 numeric_spectrum, read_input, write_card, write_input_template)
from tests.test_conjugation import _sddm_model

MAJORANA = {'ChiI1Y0Lw', 'ChiI2Y3Lw', 'ChiI2Ym3Lw'}


def test_input_parameters_are_the_external_lagrangian_entries():
    names = [p.name for p in input_parameters(_sddm_model(MAJORANA))]
    assert names == ['MChiI1Y0', 'MChiI2Y3ChiI2Ym3', 'RyChiI1Y0ChiI2Y3', 'IyChiI1Y0ChiI2Y3',
                     'RyChiI1Y0ChiI2Ym3', 'IyChiI1Y0ChiI2Ym3']


def test_template_has_one_line_per_input(tmp_path):
    path = tmp_path / 'input.txt'
    write_input_template(_sddm_model(MAJORANA), path)
    lines = [l for l in path.read_text().splitlines() if not l.startswith('#')]
    assert [l.split()[:2] for l in lines] == [
        ['MChiI1Y0', '1.000000e+00'], ['MChiI2Y3ChiI2Ym3', '1.000000e+00'],
        ['RyChiI1Y0ChiI2Y3', '1.000000e+00'], ['IyChiI1Y0ChiI2Y3', '0.000000e+00'],
        ['RyChiI1Y0ChiI2Ym3', '1.000000e+00'], ['IyChiI1Y0ChiI2Ym3', '0.000000e+00']]
    assert lines[0].endswith('# [GeV] Majorana mass of ChiI1Y0 (BSMMASS 1)')
    assert lines[3].endswith('(imaginary part) (IMBSMYUKAWA 1)')


def test_template_is_not_overwritten_by_default(tmp_path):
    path = tmp_path / 'input.txt'
    model = _sddm_model(MAJORANA)
    write_input_template(model, path)
    with pytest.raises(FileExistsError):
        write_input_template(model, path)
    write_input_template(model, path, overwrite=True)


def test_template_reads_back_to_the_defaults(tmp_path):
    path = tmp_path / 'input.txt'
    model = _sddm_model(MAJORANA)
    write_input_template(model, path)
    assert read_input(path) == {
        'MChiI1Y0': 1.0, 'MChiI2Y3ChiI2Ym3': 1.0, 'RyChiI1Y0ChiI2Y3': 1.0,
        'IyChiI1Y0ChiI2Y3': 0.0, 'RyChiI1Y0ChiI2Ym3': 1.0, 'IyChiI1Y0ChiI2Ym3': 0.0}


def test_edited_values_are_read(tmp_path):
    path = tmp_path / 'input.txt'
    model = _sddm_model(MAJORANA)
    write_input_template(model, path)
    text = path.read_text().replace('MChiI1Y0            1.000000e+00', 'MChiI1Y0 2.5e2')
    text = text.replace('MChiI2Y3ChiI2Ym3    1.000000e+00', 'MChiI2Y3ChiI2Ym3 1.0d3')
    path.write_text(text)
    values = read_input(path)
    assert (values['MChiI1Y0'], values['MChiI2Y3ChiI2Ym3']) == (250.0, 1000.0)


def test_all_input_problems_are_reported_together(tmp_path):
    path = tmp_path / 'input.txt'
    path.write_text('MChiI1Y0 100   # ok\n'
                    'MChiI1Y0 200\n'
                    'RyChiI1Y0ChiI2Y3 abc\n'
                    'RyChiI1Y0ChiI2Ym3 nan\n'
                    'IyChiI1Y0ChiI2Y3 0 extra\n')
    with pytest.raises(ValueError) as error:
        read_input(path)
    message = str(error.value)
    for problem in ('MChiI1Y0 given twice',
                    "RyChiI1Y0ChiI2Y3 = 'abc' is not a finite number",
                    "RyChiI1Y0ChiI2Ym3 = 'nan' is not a finite number",
                    'expected "name value"'):
        assert problem in message


def test_evaluate_reports_missing_and_unknown_names():
    values = {p.name: 1.0 for p in input_parameters(_sddm_model(MAJORANA))}
    del values['MChiI1Y0']
    values.update(MChiI1Y00=1.0, RyChiI1Y0ChiI2Y3=float('nan'))
    with pytest.raises(ValueError) as error:
        numeric_mass_matrices(_sddm_model(MAJORANA))(values)
    message = str(error.value)
    for problem in ('missing parameter MChiI1Y0', 'unknown parameter MChiI1Y00',
                    'RyChiI1Y0ChiI2Y3 = nan is not a finite number'):
        assert problem in message


def _matrix(names, block):
    return mass_matrices(_sddm_model(names))[block][1]


def _same(actual, expected):
    return sp.simplify(actual - sp.Matrix(expected)) == sp.zeros(*actual.shape)


S = sp.Symbol
v = S('v')
r2, r3, r6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)


def test_majorana_block_and_dirac_pair():
    y, yp, M = S('yChiI1Y0ChiI2Y3'), S('yChiI1Y0ChiI2Ym3'), S('MChiI2Y3ChiI2Ym3')
    assert _same(_matrix(MAJORANA, 'UMixMZ1Z2Z3'), [
        [S('MChiI1Y0'), -v * y / r2, v * yp / r2],
        [-v * y / r2, 0, M],
        [v * yp / r2, M, 0]])
    blocks, matrix = mass_matrices(_sddm_model(MAJORANA))['UMixDW1']
    assert [b.standard_symbol for b in blocks] == ['UMixDW1', 'VMixDW1']
    assert _same(matrix, [[-M]])


def test_dirac_block_rows_left_columns_right():
    names = {'ChiI1Y0D1F1Lw', 'ChiI1Y0D1F2Lw', 'ChiI1Y0Dm1F1Lw', 'ChiI1Y0Dm1F2Lw',
             'ChiI2Y3D1Lw', 'ChiI2Ym3Dm1Lw'}
    assert _same(_matrix(names, 'UMixDZ1Z2Z3'), [
        [S('MChiI1Y0D1F1ChiI1Y0Dm1F1'), 0, v * S('yChiI1Y0D1F1ChiI2Ym3Dm1') / r2],
        [0, S('MChiI1Y0D1F2ChiI1Y0Dm1F2'), v * S('yChiI1Y0D1F2ChiI2Ym3Dm1') / r2],
        [-v * S('yChiI2Y3D1ChiI1Y0Dm1F1') / r2, -v * S('yChiI2Y3D1ChiI1Y0Dm1F2') / r2,
         S('MChiI2Y3D1ChiI2Ym3Dm1')]])


def test_scalar_blocks_with_hermitian_conjugate_and_real_parts():
    # real block rows: S, phi_R, phi_I with Phi0 = (phi_R + i phi_I)/sqrt(2)
    names = {'PhiI1Y0', 'PhiI2Y3'}
    mu, lam = S('muPhiI1Y0PhiI2Y3'), S('lambdaPhiI2Y3PhiI2Y3')
    p1, p2 = S('lambdaPhiI2Y3PhiI2Y3barP1'), S('lambdaPhiI2Y3PhiI2Y3barP2')
    m2 = S('MPhiI2Y3')**2 + p1 * v**2 / 4 - r3 * p2 * v**2 / 12
    c = sp.conjugate
    a, abar = r3 * lam * v**2 / 6, r3 * c(lam) * v**2 / 6      # v real
    assert _same(_matrix(names, 'UMixRZ1Z2Z3'), [
        [S('MPhiI1Y0')**2 + S('lambdaPhiI1Y0PhiI1Y0') * v**2,
         -v * (mu + c(mu)) / 2, -sp.I * v * (mu - c(mu)) / 2],
        [-v * (mu + c(mu)) / 2, m2 + a + abar, sp.I * (a - abar)],
        [-sp.I * v * (mu - c(mu)) / 2, sp.I * (a - abar), m2 - a - abar]])
    assert _same(_matrix(names, 'UMixCW1'),
                 [[S('MPhiI2Y3')**2 + p1 * v**2 / 4 + r3 * p2 * v**2 / 12]])


def test_real_triplet_components_are_degenerate():
    # .fr: - 1/2 M^2 T^*_i T_i, the only splitting would come from H^dagger T H (absent)
    names = {'PhiI1Y0F1', 'PhiI1Y0F2', 'PhiI3Y0'}
    m2 = S('MPhiI3Y0')**2 - r6 * S('lambdaPhiI3Y0PhiI3Y0') * v**2 / 6
    assert _same(_matrix(names, 'UMixCW1'), [[m2]])
    assert sp.simplify(_matrix(names, 'UMixRZ1Z2Z3')[2, 2] - m2) == 0


def _majorana_values(**changes):
    values = {'MChiI1Y0': 100.0, 'MChiI2Y3ChiI2Ym3': 200.0,
              'RyChiI1Y0ChiI2Y3': 0.5, 'IyChiI1Y0ChiI2Y3': 0.0,
              'RyChiI1Y0ChiI2Ym3': 0.3, 'IyChiI1Y0ChiI2Ym3': 0.0}
    return {**values, **changes}


def test_numeric_majorana_matrices():
    evaluate = numeric_mass_matrices(_sddm_model(MAJORANA))
    w = VEV
    matrices = evaluate(_majorana_values())
    assert np.allclose(matrices['UMixMZ1Z2Z3'], [[100, -w * 0.5 / np.sqrt(2), w * 0.3 / np.sqrt(2)],
                                                [-w * 0.5 / np.sqrt(2), 0, 200],
                                                [w * 0.3 / np.sqrt(2), 200, 0]])
    assert np.allclose(matrices['UMixDW1'], [[-200]])
    complex_y = evaluate(_majorana_values(IyChiI1Y0ChiI2Y3=0.1))['UMixMZ1Z2Z3']
    assert np.isclose(complex_y[0, 1], -w * (0.5 + 0.1j) / np.sqrt(2))


def test_numeric_real_block_is_real_for_complex_couplings():
    frm = _sddm_model({'PhiI1Y0', 'PhiI2Y3'})
    values = {p.name: 0.7 for p in input_parameters(frm)}
    values.update({'RmuPhiI1Y0PhiI2Y3': 30.0, 'ImuPhiI1Y0PhiI2Y3': 20.0,
                   'IlambdaPhiI2Y3PhiI2Y3': 0.4})
    real_block = numeric_mass_matrices(frm)(values)['UMixRZ1Z2Z3']
    assert np.allclose(real_block.imag, 0)
    assert np.isclose(real_block[0, 2].real, VEV * 20.0)              # S - phi_I = v Im(mu)


def test_numeric_matrices_from_the_input_file(tmp_path):
    path = tmp_path / 'input.txt'
    frm = _sddm_model(MAJORANA)
    write_input_template(frm, path)
    matrices = numeric_mass_matrices(frm)(read_input(path))
    assert np.allclose(matrices['UMixDW1'], [[-1]])


# top-down point MB0 = 250, M = 200, y(B0 Hd) = 1.5, y(B0 Hu) = 0.2: here ChiI2Ym3 = Hd,
# ChiI2Y3 = Hu, rows (B0, Hu0, Hd0)
SDDM_POINT = {'MChiI1Y0': 250.0, 'MChiI2Y3ChiI2Ym3': 200.0,
              'RyChiI1Y0ChiI2Ym3': 1.5, 'IyChiI1Y0ChiI2Ym3': 0.0,
              'RyChiI1Y0ChiI2Y3': 0.2, 'IyChiI1Y0ChiI2Y3': 0.0}


def test_spectrum_reproduces_the_top_down_sddm_point():
    masses, mixings = numeric_spectrum(_sddm_model(MAJORANA))(SDDM_POINT)
    assert masses.keys() == {1000022, 1000023, 1000024, 1000025}
    assert np.allclose([masses[p] for p in (1000022, 1000023, 1000025, 1000024)],
                       [109.28073066, 289.83963497, 430.55890431, 200.0])
    assert set(mixings) == {'UMixMZ1Z2Z3'}
    U = mixings['UMixMZ1Z2Z3']
    N, signs = U.conj().T, np.ones(3)              # SLHA1 form: real N, signed masses
    for k in range(3):
        if np.allclose(U[:, k].real, 0):
            N[k] *= -1j
            signs[k] = -1
    top_down = np.array([[-0.46777138, 0.35969618, 0.80734664],     # columns (B0, Hd0, Hu0)
                         [0.39122139, -0.73482453, 0.55405661],
                         [0.79255016, 0.5750231, 0.20300905]])
    assert np.allclose(N.real, top_down[:, [0, 2, 1]], atol=1e-7)
    assert np.allclose(N.imag, 0) and list(signs) == [1, -1, 1]


def test_spectrum_of_dirac_blocks_has_u_and_v():
    names = {'ChiI1Y0D1F1Lw', 'ChiI1Y0D1F2Lw', 'ChiI1Y0Dm1F1Lw', 'ChiI1Y0Dm1F2Lw',
             'ChiI2Y3D1Lw', 'ChiI2Ym3Dm1Lw'}
    frm = _sddm_model(names)
    values = {p.name: 0.0 for p in input_parameters(frm)}
    values.update({'MChiI1Y0D1F1ChiI1Y0Dm1F1': 100.0, 'MChiI1Y0D1F2ChiI1Y0Dm1F2': 300.0,
                   'MChiI2Y3D1ChiI2Ym3Dm1': 200.0, 'RyChiI1Y0D1F1ChiI2Ym3Dm1': 0.3,
                   'RyChiI2Y3D1ChiI1Y0Dm1F2': 0.5})
    masses, mixings = numeric_spectrum(frm)(values)
    assert set(mixings) == {'UMixDZ1Z2Z3', 'VMixDZ1Z2Z3'}
    M = numeric_mass_matrices(frm)(values)['UMixDZ1Z2Z3']
    U, V = mixings['UMixDZ1Z2Z3'], mixings['VMixDZ1Z2Z3']
    assert np.allclose(U.T @ M @ V, np.diag([masses[p] for p in (1000022, 1000023, 1000025)]))


def test_spectrum_of_scalar_blocks():
    frm = _sddm_model({'PhiI1Y0', 'PhiI2Y3'})
    values = {p.name: 0.0 for p in input_parameters(frm)}
    values.update({'MPhiI1Y0': 100.0, 'MPhiI2Y3': 300.0, 'RmuPhiI1Y0PhiI2Y3': 50.0,
                   'RlambdaPhiI2Y3PhiI2Y3': 0.5})
    masses, mixings = numeric_spectrum(frm)(values)
    assert set(mixings) == {'UMixRZ1Z2Z3'} and np.isclose(masses[1000024], 300.0)
    O = mixings['UMixRZ1Z2Z3']
    assert np.isrealobj(O) and np.allclose(O.T @ O, np.eye(3))
    values['RmuPhiI1Y0PhiI2Y3'] = 1000.0
    with pytest.raises(ValueError, match='UMixRZ1Z2Z3: negative physical mass'):
        numeric_spectrum(frm)(values)


def _template(frm, path, gf='1.166370e-05', skip=None):
    '''MG5-like param_card of the model: SM entries, every BSM entry at 1, a BR line.'''
    blocks = {}
    for p in input_parameters(frm):
        blocks.setdefault(p.block, []).append(f'    {p.order_block} 1.000000e+00 # {p.name}')
    matrices = {m.standard_symbol: m for m in frm.mixings.get_mixing_matrices()}
    labels = {q.pdg: frm.renames['particle'].get(q.label, q.label)
              for m in matrices.values() for q in m.particles}
    pdgs = sorted(labels)
    blocks['MASS'] = ['   23 9.118760e+01 # MZ'] + [f'  {pdg} 1.000000e+00 # M{labels[pdg]}'
                                                   for pdg in pdgs]
    for p in frm.parameters():
        if p.sector == 'Mixing' and p.kind == 'External':
            n = matrices[p.base].dim
            blocks[p.block] = [f'    {i}   {j}  0.0e+00 # {p.name}{i}x{j}'
                               for i in range(1, n + 1) for j in range(1, n + 1)]
    lines = ['Block SMINPUTS ', '    1 1.279000e+02 # aEWM1 ', f'    2 {gf} # Gf ']
    for name, entries in blocks.items():
        lines += ['', f'Block {name.lower()} '] + [e for e in entries if e.split('# ')[1] != skip]
    lines += ['DECAY  23 2.441404e+00 # WZ']
    lines += [f'DECAY  {pdg} 1.000000e+00 # W{pdg}' for pdg in pdgs]
    lines += ['   1.000000e+00   2   1000022   23 # BR']
    path.write_text('\n'.join(lines) + '\n')
    return path


def _read_card(path):
    '''{(BLOCK, *indices): value} and {('DECAY', pdg): width} of a card.'''
    data, block = {}, None
    for line in path.read_text().splitlines():
        fields = line.split('#')[0].split()
        if not fields:
            continue
        if fields[0].upper() in ('BLOCK', 'DECAY'):
            block = fields[1].upper() if fields[0].upper() == 'BLOCK' else None
            if fields[0].upper() == 'DECAY':
                data[('DECAY', int(fields[1]))] = fields[2]
            continue
        data[(block, *(int(f) for f in fields[:-1]))] = float(fields[-1])
    return data


def test_card_of_the_top_down_sddm_point(tmp_path):
    frm = _sddm_model(MAJORANA)
    template = _template(frm, tmp_path / 'template.dat')
    before = template.read_text()
    card = _read_card(write_card(frm, template, tmp_path / 'card.dat', SDDM_POINT))
    assert template.read_text() == before
    assert card[('BSMMASS', 1)] == 250 and card[('BSMMASS', 2)] == 200
    assert card[('BSMYUKAWA', 2)] == 1.5 and card[('IMBSMYUKAWA', 2)] == 0
    assert np.allclose([card[('MASS', p)] for p in (1000022, 1000023, 1000025, 1000024)],
                       [109.28073066, 289.83963497, 430.55890431, 200.0])
    U = numeric_spectrum(frm)(SDDM_POINT)[1]['UMixMZ1Z2Z3']
    for i, j in np.ndindex(3, 3):
        assert card[('UMIXMZ1Z2Z3', i + 1, j + 1)] == U[i, j].real
        assert card[('IMUMIXMZ1Z2Z3', i + 1, j + 1)] == U[i, j].imag
    assert np.allclose([card[('IMUMIXMZ1Z2Z3', i, 2)] for i in (1, 2, 3)],    # rows B0 Hu0 Hd0
                       [-0.39122139, -0.55405661, 0.73482453])
    assert card[('MASS', 23)] == 91.1876 and card[('DECAY', 23)] == '2.441404e+00'
    assert all(card[('DECAY', p)] == 'Auto' for p in (1000022, 1000023, 1000024, 1000025))
    text = (tmp_path / 'card.dat').read_text()
    assert ('Block umixmz1z2z3 \n# rows (gauge): ChiI1Y0q0 ChiI2Y3q0 ChiI2Ym3q0; '
            'columns (mass): Z1 Z2 Z3\n') in text
    assert '# BR' not in text and '# RUMixMZ1Z2Z32x2' in text


def test_card_of_a_real_scalar_block(tmp_path):
    frm = _sddm_model({'PhiI1Y0', 'PhiI2Y3'})
    values = {p.name: 0.0 for p in input_parameters(frm)}
    values.update({'MPhiI1Y0': 100.0, 'MPhiI2Y3': 300.0, 'RmuPhiI1Y0PhiI2Y3': 50.0,
                   'RlambdaPhiI2Y3PhiI2Y3': 0.5})
    path = write_card(frm, _template(frm, tmp_path / 'template.dat'), tmp_path / 'card.dat', values)
    O = numeric_spectrum(frm)(values)[1]['UMixRZ1Z2Z3']
    card = _read_card(path)
    assert all(card[('UMIXRZ1Z2Z3', i + 1, j + 1)] == O[i, j] for i, j in np.ndindex(3, 3))
    assert ('# rows (gauge): PhiI1Y0q0 PhiI2Y3q0_R PhiI2Y3q0_I; columns (mass): Z1 Z2 Z3'
            in path.read_text())


def test_card_rejects_a_template_of_another_model_or_vev(tmp_path):
    frm = _sddm_model(MAJORANA)
    template = _template(frm, tmp_path / 'template.dat', skip='IUMixMZ1Z2Z32x3')
    with pytest.raises(ValueError, match='entries missing: IMUMIXMZ1Z2Z3 2 3$'):
        write_card(frm, template, tmp_path / 'card.dat', SDDM_POINT)
    template = _template(frm, tmp_path / 'template.dat')
    text = template.read_text().replace('# RyChiI1Y0ChiI2Y3', '# RyB0Hd')
    template.write_text(text.replace('# MZ2', '# MZ3'))
    with pytest.raises(ValueError, match='BSMYUKAWA 1 is RyB0Hd, expected RyChiI1Y0ChiI2Y3\n'
                                         '  MASS 1000023 is MZ3, expected MZ2'):
        write_card(frm, template, tmp_path / 'card.dat', SDDM_POINT)
    template = _template(frm, tmp_path / 'template.dat', gf='1.166380e-05')
    with pytest.raises(ValueError, match='does not give the vev'):
        write_card(frm, template, tmp_path / 'card.dat', SDDM_POINT)
    assert not (tmp_path / 'card.dat').exists()


def test_card_filled_in_place_from_its_own_inputs(tmp_path):
    frm = _sddm_model(MAJORANA)
    reference = write_card(frm, _template(frm, tmp_path / 'template.dat'),
                           tmp_path / 'reference.dat', SDDM_POINT).read_text()
    card = _template(frm, tmp_path / 'card.dat')
    text = card.read_text()
    for name, value in SDDM_POINT.items():
        text = text.replace(f'1.000000e+00 # {name}\n', f'{value:.16e} # {name}\n')
        text = text.replace(f'0.000000e+00 # {name}\n', f'{value:.16e} # {name}\n')
    card.write_text(text)
    assert write_card(frm, card) == card
    assert card.read_text() == reference
    write_card(frm, card)                                   # idempotent
    assert card.read_text() == reference


def test_card_uses_the_output_names_of_renamed_particles(tmp_path):
    from tests.test_golden import STOP, _stop_model
    frm = _stop_model(STOP['stop_weylstop_scalargluino'])
    assert frm.renames['particle'] == {'B': 'B1'}               # B is reserved in SMMod.fr
    values, k = {}, 0
    for p in input_parameters(frm):
        values[p.name] = 300.0 + 150.0 * k if p.sector == 'Mass' else 0.1
        k += p.sector == 'Mass'
    path = write_card(frm, _template(frm, tmp_path / 'template.dat'), tmp_path / 'card.dat', values)
    (pdg,) = [q.pdg for m in frm.mixings.get_mixing_matrices() for q in m.particles if q.label == 'B']
    assert f'  {pdg} ' in path.read_text() and '# MB1' in path.read_text()
    assert _read_card(path)[('MASS', pdg)] == numeric_spectrum(frm)(values)[0][pdg]
