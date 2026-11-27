"""
.. module:: spectrum
   :synopsis: Diagonalizes the mass matrices of the mixing blocks.

.. moduleauthor:: Lucas Heck <lucasdealexandri.heck@gmail.com>
"""

import warnings

import numpy as np

TOLERANCE = 1e-8


def _largest(column) -> int:
    '''Index of the entry of largest modulus (first one among numerically equal).'''
    size = np.round(np.abs(column), 10)
    return int(np.argmax(size))


def _phases(U, sign_only=False) -> np.ndarray:
    '''Column factors that apply the phase convention to U: the entry of largest modulus of
    each column becomes real and positive; with sign_only (real scalars, Majorana
    fermions) only its sign is fixed: a positive real part, or a positive imaginary part
    if the real part vanishes.'''
    factors = np.ones(U.shape[1], dtype=complex)
    for j in range(U.shape[1]):
        z = U[_largest(U[:, j]), j]
        if sign_only:
            reference = z.real if abs(z.real) > TOLERANCE * abs(z) else z.imag
            factors[j] = -1 if reference < 0 else 1
        else:
            factors[j] = np.conj(z) / abs(z)
    return factors


def _check_masses(masses, block) -> None:
    '''Positive masses (relative tolerance TOLERANCE); degenerate ones are only warned
    about, since their mixing matrix is then not unique.'''
    scale = max(float(np.max(masses)), 1.0)
    if masses[0] <= TOLERANCE * scale:
        raise ValueError(f'{block}: zero mass ({masses[0]:.3e} GeV)')
    for low, high in zip(masses[:-1], masses[1:]):
        if high - low <= TOLERANCE * scale:
            warnings.warn(f'{block}: degenerate masses ({low:.10g}, {high:.10g} GeV): '
                          'the mixing matrix is not unique')


def _scalar_masses(squared, block) -> np.ndarray:
    '''Masses from the squared masses; ValueError for a negative physical mass squared
    (and the checks of _check_masses).'''
    scale = max(float(np.max(np.abs(squared))), 1.0)
    if squared[0] < -TOLERANCE * scale:
        raise ValueError(f'{block}: negative physical mass squared ({squared[0]:.6g} GeV^2)')
    masses = np.sqrt(np.clip(squared, 0, None))
    _check_masses(masses, block)
    return masses


def diagonalize_real_scalar(M2, block='real scalar block'):
    '''(m, O) with O^T M2 O = diag(m^2), m ascending, O real orthogonal. M2 may be
    given complex: its imaginary part must vanish.'''
    M2 = np.asarray(M2, dtype=complex)
    if np.max(np.abs(M2.imag)) > TOLERANCE * max(np.max(np.abs(M2)), 1.0):
        raise ValueError(f'{block}: real scalar mass matrix with an imaginary part')
    squared, O = np.linalg.eigh(M2.real)
    masses = _scalar_masses(squared, block)
    return masses, O * _phases(O, sign_only=True).real


def diagonalize_complex_scalar(M2, block='complex scalar block'):
    '''(m, U) with U^dag M2 U = diag(m^2), m ascending, U unitary.'''
    squared, U = np.linalg.eigh(np.asarray(M2, dtype=complex))
    masses = _scalar_masses(squared, block)
    return masses, U * _phases(U)


def diagonalize_dirac(M, block='Dirac block'):
    '''(m, U, V) with U^T M V = diag(m), m ascending, U and V unitary. From the SVD
    M = A diag(m) B^dag: U = A^*, V = B. The phase convention is applied to U, and V
    takes the opposite phases so that U^T M V is unchanged.'''
    A, singular, Bh = np.linalg.svd(np.asarray(M, dtype=complex))
    order = np.argsort(singular)
    masses, U, V = singular[order], A.conj()[:, order], Bh.conj().T[:, order]
    _check_masses(masses, block)
    factors = _phases(U)
    return masses, U * factors, V * factors.conj()


def _takagi_rotation(a11, a22, a12):
    '''2x2 Jacobi rotation of Hahn, physics/0607103, eqs. (18)-(31): returns (c, s, e) of
    R = [[c, s e], [-s e^*, c]] with R^* A R^dag diagonal, s = t c, e = e^{i phi}.'''
    numerator = np.conj(a11) * a12 + a22 * np.conj(a12)
    e = numerator / abs(numerator) if abs(numerator) > 0 else 1.0
    delta = (e * a11 - np.conj(e) * a22) / 2          # Delta-tilde, proportional to a12
    ratio = (delta / a12).real
    t = 1.0 / (ratio + np.copysign(np.hypot(ratio, 1.0), ratio))   # larger denominator
    c = 1.0 / np.sqrt(1.0 + t * t)
    return c, t * c, e


def takagi(M, block='Majorana block', max_sweeps=50):
    '''(m, U) with U^T M U = diag(m), m ascending and positive, U unitary, for a complex
    symmetric M. Cyclic Jacobi sweeps with the 2x2 rotations of Hahn (physics/0607103,
    sec. 4.2): R^* A R^dag zeroes A_pq; the product W of the rotations gives
    W^* M W^dag = diag(sigma); the phases of sigma are moved into W, and U = W^dag.'''
    A = np.array(M, dtype=complex)
    if np.max(np.abs(A - A.T)) > TOLERANCE * max(np.max(np.abs(A)), 1.0):
        raise ValueError(f'{block}: Majorana mass matrix is not symmetric')
    n = len(A)
    W = np.eye(n, dtype=complex)
    norm = np.linalg.norm(A)
    for _ in range(max_sweeps):
        off = np.linalg.norm(A - np.diag(np.diag(A)))
        if off <= 1e-15 * norm:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(A[p, q]) <= 1e-18 * norm:
                    continue
                c, s, e = _takagi_rotation(A[p, p], A[q, q], A[p, q])
                R = np.eye(n, dtype=complex)
                R[p, p] = R[q, q] = c
                R[p, q], R[q, p] = s * e, -s * np.conj(e)
                A = R.conj() @ A @ R.conj().T
                A = (A + A.T) / 2
                W = R @ W
    else:
        raise ValueError(f'{block}: Takagi factorization did not converge')
    sigma = np.diag(A)
    W = np.exp(1j * np.angle(sigma) / 2)[:, None] * W  # (P W)^* M (P W)^dag = diag(|sigma|)
    masses, U = np.abs(sigma), W.conj().T
    order = np.argsort(masses)
    masses, U = masses[order], U[:, order]
    _check_masses(masses, block)
    return masses, U * _phases(U, sign_only=True)
