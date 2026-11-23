"""Diagonalisation of the mass matrices (tools/spectrum.py)."""

import numpy as np
import pytest

from uvbuilder.tools.spectrum import (_largest, diagonalize_complex_scalar, diagonalize_dirac,
                            diagonalize_real_scalar, takagi)

RNG = np.random.default_rng(20261118)


def _complex(n):
    return RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))


def _largest_entries(U):
    return np.array([U[_largest(U[:, j]), j] for j in range(U.shape[1])])


def _unitary(U):
    return np.allclose(U.conj().T @ U, np.eye(len(U)))


@pytest.mark.parametrize('n', [1, 2, 3, 6])
def test_real_scalar(n):
    X = RNG.normal(size=(n, n))
    M2 = X @ X.T + n * np.eye(n)
    m, O = diagonalize_real_scalar(M2 + 0j)
    assert np.isrealobj(O) and np.allclose(O.T @ O, np.eye(n))
    assert np.allclose(O.T @ M2 @ O, np.diag(m**2)) and np.all(np.diff(m) > 0)
    assert np.all(_largest_entries(O) > 0)


@pytest.mark.parametrize('n', [1, 2, 3, 6])
def test_complex_scalar(n):
    X = _complex(n)
    M2 = X @ X.conj().T + n * np.eye(n)
    m, U = diagonalize_complex_scalar(M2)
    assert _unitary(U) and np.allclose(U.conj().T @ M2 @ U, np.diag(m**2))
    largest = _largest_entries(U)
    assert np.allclose(largest.imag, 0) and np.all(largest.real > 0)


@pytest.mark.parametrize('n', [1, 2, 3, 6])
def test_dirac(n):
    M = _complex(n)
    m, U, V = diagonalize_dirac(M)
    assert _unitary(U) and _unitary(V)
    assert np.allclose(U.T @ M @ V, np.diag(m)) and np.all(np.diff(m) > 0)
    largest = _largest_entries(U)
    assert np.allclose(largest.imag, 0) and np.all(largest.real > 0)


@pytest.mark.parametrize('n', [1, 2, 3, 4, 6, 8])
def test_takagi(n):
    X = _complex(n)
    M = X + X.T
    m, U = takagi(M)
    assert _unitary(U)
    assert np.allclose(U.T @ M @ U, np.diag(m), atol=1e-10)
    assert np.allclose(m, np.sort(np.linalg.svd(M, compute_uv=False)))
    largest = _largest_entries(U)
    assert np.all(np.where(np.abs(largest.real) > 1e-8, largest.real, largest.imag) > 0)


def test_takagi_negative_diagonal_gets_a_phase():
    m, U = takagi(np.diag([-3.0, 5.0]))
    assert np.allclose(m, [3, 5]) and np.allclose(U, np.diag([1j, 1]))


def test_takagi_near_degenerate_is_accepted():
    # higgsino-like pair split by a small mixing
    M = np.array([[100.0, 1.0, 1.0], [1.0, 0.0, 200.0], [1.0, 200.0, 0.0]])
    m, U = takagi(M)
    assert np.allclose(U.T @ M @ U, np.diag(m), atol=1e-9) and np.all(np.diff(m) > 1e-3)


@pytest.mark.parametrize('M, message', [
    ([[0, 1], [1, 0]], 'degenerate'),          # Haber, eq. (167)
    ([[1, 1j], [1j, -1]], 'zero mass'),        # Haber, eqs. (157)-(158)
    ([[1, 2], [3, 4]], 'not symmetric'),
])
def test_takagi_rejects(M, message):
    with pytest.raises(ValueError, match=message):
        takagi(np.array(M, dtype=complex))


def test_scalar_rejections():
    with pytest.raises(ValueError, match='negative physical mass'):
        diagonalize_real_scalar(np.diag([-4.0, 9.0]))
    with pytest.raises(ValueError, match='imaginary part'):
        diagonalize_real_scalar(np.array([[4, 1j], [-1j, 9]]))
    with pytest.raises(ValueError, match='degenerate'):
        diagonalize_complex_scalar(np.diag([4.0, 4.0]))
    with pytest.raises(ValueError, match='zero mass'):
        diagonalize_dirac(np.diag([0.0, 3.0]))
