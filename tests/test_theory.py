"""The signal model: closed forms against their series and Fourier definitions."""

import numpy as np
import pytest

from mess import complex_sum, magnitude_sum, theory

CASES = [(50, 10e-3, 0.83, 0.075, 0.18), (20, 9.83e-3, 4.16, 1.65, 0.059), (80, 5e-3, 1.56, 0.083, 0.32)]


@pytest.mark.parametrize("fa, tr, t1, t2, t2p", CASES)
def test_fourier_coefficients(fa, tr, t1, t2, t2p):
    """F_k^+ (Eq. 1) are the Fourier coefficients of M_xy^+(theta) (Eq. S1.1)."""
    n = 4096
    theta = 2 * np.pi * np.arange(n) / n
    coeffs = np.fft.fft(theory.bssfp_profile(theta, fa, tr, t1, t2)) / n
    k = np.arange(-6, 7)
    assert np.allclose(coeffs[k % n], theory.epg_coefficient(k, fa, tr, t1, t2), atol=1e-12)
    assert np.all(theory.epg_coefficient(k, fa, tr, t1, t2)[k >= 0] > 0)      # Eq. 2
    assert np.all(theory.epg_coefficient(k, fa, tr, t1, t2)[k < 0] < 0)


@pytest.mark.parametrize("fa, tr, t1, t2, t2p", CASES)
def test_closed_forms(fa, tr, t1, t2, t2p):
    """Series (Eqs. 5-7) and closed forms (Eqs. S4.2-S4.4) agree."""
    te, dte = 0.45 * tr, 0.075 * tr
    dw0 = 2 * np.pi * np.linspace(-150, 150, 7)
    dpsi, dphi = 0.7, 1.9
    b = theory.exponential_coefficients(te, 0.0, fa, tr, t1, t2, t2p)
    assert np.allclose(theory.m_bssfp(te, fa, tr, t1, t2, t2p, dw0, dpsi),
                       theory.m_infinity(dw0 * te, dw0 * tr + dpsi, b), rtol=1e-8)
    orders = np.arange(60, -61, -1)
    s = theory.mess_echoes(orders, te, dte, fa, tr, t1, t2, t2p, dw0, dpsi)
    c = theory.exponential_coefficients(te, dte, fa, tr, t1, t2, t2p)
    assert np.allclose(complex_sum(s, orders, dphi),
                       theory.m_infinity(dw0 * te, dw0 * (tr - dte) + dpsi + dphi, c), rtol=1e-6)
    assert np.allclose(magnitude_sum(s, orders, dphi), theory.m_infinity(0.0, dphi, c), rtol=1e-6)


def test_magnitude_sum_ignores_off_resonance():
    orders = np.arange(2, -4, -1)
    s = theory.mess_echoes(orders, 4.5e-3, 0.75e-3, 50, 9.83e-3, 1.56, 0.083, 0.32,
                           dw0=2 * np.pi * np.linspace(-300, 300, 11))
    m = magnitude_sum(s, orders, np.pi)
    assert np.allclose(m, m[0])
