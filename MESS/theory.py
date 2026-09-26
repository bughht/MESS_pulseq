"""Steady-state signal model of bSSFP and MESS, in the notation of the paper.

=================  ============================================================
``alpha``          flip angle (degrees in this module)
``TR``, ``TE_k``   repetition time, echo time of dephasing order ``k``
``E1``, ``E2``     ``exp(-TR/T1)``, ``exp(-TR/T2)``
``theta``          ``Δω_0 TR + Δψ``, phase accrued between two RF pulses
``M_xy^+(θ)``      complex bSSFP profile right after the RF pulse   (SI S1)
``F_k^+``          dephasing orders = Fourier coefficients of M_xy^+ (Eqs. 1, 2)
``S_k``            MESS echo of order ``k``                        (Eq. 3)
``M_bSSFP``        coherent sum of all orders at one TE            (Eq. 5)
``M_Σ``            complex sum of the acquired echoes              (Eq. 6)
``M_Σ|·|``         magnitude sum of the acquired echoes            (Eq. 7)
``λ±, μ±``         piecewise-exponential echoes                    (SI S3)
``M_∞``            shorthand for the closed forms of SI S4.2-S4.4
=================  ============================================================

``Δω_0`` is in rad/s (``2π`` times the off-resonance in Hz); the RF and virtual
phase-cycle increments ``Δψ`` and ``Δφ`` are in radians.  ``M_Σ`` and
``M_Σ|·|`` of model echoes are obtained with :func:`mess.recon.complex_sum` and
:func:`mess.recon.magnitude_sum`, exactly as for measured echoes.
"""

import numpy as np


def relaxation(tr, t1, t2):
    """``E1 = exp(-TR/T1)`` and ``E2 = exp(-TR/T2)``."""
    return np.exp(-tr / t1), np.exp(-tr / t2)


def bssfp_profile(theta, flip_angle, tr, t1, t2, m0=1.0):
    """Complex bSSFP profile ``M_xy^+(θ)`` immediately after the RF pulse (SI S1).

    ``M_xy^+(θ) = Q sin α (1 - E2 cos θ + i E2 sin θ)``.  Its magnitude is
    smallest at ``θ = 0`` (stop band) and largest at ``θ = π`` (pass band).
    """
    alpha = np.deg2rad(flip_angle)
    e1, e2 = relaxation(tr, t1, t2)
    cos_a, cos_t = np.cos(alpha), np.cos(theta)
    q = m0 * (1 - e1) / ((1 - e1 * cos_a) * (1 - e2 * cos_t) - e2 * (e1 - cos_a) * (e2 - cos_t))
    return q * np.sin(alpha) * (1 - e2 * cos_t + 1j * e2 * np.sin(theta))


def _b_and_c(flip_angle, tr, t1, t2, m0):
    """Auxiliary parameters ``b`` and ``c`` of Eq. S2.2."""
    alpha = np.deg2rad(flip_angle)
    e1, e2 = relaxation(tr, t1, t2)
    big_a = 1 - e1 * np.cos(alpha)
    big_b = e1 - np.cos(alpha)
    a = e2 * (big_b - big_a) / (big_a - big_b * e2 ** 2)
    root = np.sqrt(1 - a ** 2)
    b = -a / (1 + root)                        # = (sqrt(1 - a^2) - 1) / a, stable for a -> 0
    c = m0 * (1 - e1) * np.sin(alpha) / ((big_a - big_b * e2 ** 2) * root)
    return b, c, e2


def epg_coefficient(k, flip_angle, tr, t1, t2, m0=1.0):
    """Steady-state dephasing order ``F_k^+ = c (b^|k| - E2 b^|k+1|)`` (Eq. 1).

    ``F_k^+`` is real, positive for ``k >= 0`` and negative for ``k < 0``, i.e.
    ``F_k^+ = |F_k^+| exp(i [u(k) - 1] π)`` (Eq. 2).  The ``F_k^+`` are the
    Fourier coefficients of the bSSFP profile,
    ``M_xy^+(θ) = Σ_k F_k^+ exp(i k θ)`` (Eq. S1.1).
    """
    b, c, e2 = _b_and_c(flip_angle, tr, t1, t2, m0)
    k = np.asarray(k)
    return c * (b ** np.abs(k) - e2 * b ** np.abs(k + 1))


def polarity(k):
    """``exp(i [u(k) - 1] π)``: +1 for ``k >= 0`` and -1 for ``k < 0``."""
    return np.where(np.asarray(k) >= 0, 1.0, -1.0)


def echo_signal(k, te_k, flip_angle, tr, t1, t2, t2p=np.inf, dw0=0.0, dpsi=0.0, m0=1.0):
    """MESS echo ``S_k(TE_k, Δω_0, Δψ)`` (Eq. 3).

    ``S_k = F_k^+ exp(-TE_k/T2) exp(-|TE_k + k TR| / T2') exp(i Δω_0 (TE_k + k TR)) exp(i k Δψ)``

    ``dw0`` (rad/s) and ``dpsi`` (rad) broadcast against each other.
    """
    k = np.asarray(k)
    tau = te_k + k * tr                        # time since this pathway was refocused by the RF
    f = epg_coefficient(k, flip_angle, tr, t1, t2, m0)
    decay = np.exp(-te_k / t2) * np.exp(-np.abs(tau) / t2p)
    return f * decay * np.exp(1j * np.multiply(dw0, tau)) * np.exp(1j * k * np.asarray(dpsi))


def mess_echoes(orders, te0, delta_te, flip_angle, tr, t1, t2, t2p=np.inf,
                dw0=0.0, dpsi=0.0, m0=1.0):
    """``S_k`` of every acquired order, with ``TE_k = TE_0 - k ΔTE``.

    Returns an array of shape ``(len(orders),) + broadcast(dw0, dpsi).shape``.
    """
    dw0, dpsi = np.broadcast_arrays(np.asarray(dw0, float), np.asarray(dpsi, float))
    return np.stack([echo_signal(k, te0 - k * delta_te, flip_angle, tr, t1, t2, t2p,
                                 dw0, dpsi, m0) for k in orders])


def m_bssfp(te, flip_angle, tr, t1, t2, t2p=np.inf, dw0=0.0, dpsi=0.0, m0=1.0, tol=1e-12):
    """bSSFP signal ``M_bSSFP(Δω_0, Δψ) = Σ_k S_k`` with ``TE_k = TE`` (Eq. 5).

    The series is summed over all orders whose weight exceeds ``tol``; the
    closed form of Eq. S4.2 is :func:`m_infinity`.
    """
    b, _, _ = _b_and_c(flip_angle, tr, t1, t2, m0)
    k_max = int(np.ceil(np.log(tol) / np.log(max(abs(b), 1e-300)))) + 2
    k_max = min(max(k_max, 4), 5000)
    orders = np.arange(-k_max, k_max + 1)
    s = mess_echoes(orders, te, 0.0, flip_angle, tr, t1, t2, t2p, dw0, dpsi, m0)
    return s.sum(axis=0)


def exponential_coefficients(te0, delta_te, flip_angle, tr, t1, t2, t2p=np.inf, m0=1.0):
    """Coefficients of the piecewise-exponential echo model (Eqs. S3.2, S3.3).

    ``S_k = exp(λ_+ k + μ_+)`` for ``k >= 0`` and ``exp(λ_- k + μ_-)`` for
    ``k < 0`` (on resonance, without phase cycling).  ``μ_-`` is complex: it is
    the logarithm of the negative ``F_k^+`` of the negative orders.
    """
    b, c, e2 = _b_and_c(flip_angle, tr, t1, t2, m0)
    lam_p = np.log(b) + delta_te / t2 - (tr - delta_te) / t2p
    lam_m = -np.log(b) + delta_te / t2 + (tr - delta_te) / t2p
    mu_p = np.log(c * (1 - e2 * b)) - te0 / t2 - te0 / t2p
    mu_m = np.log(complex(c * (1 - e2 / b))) - te0 / t2 + te0 / t2p
    return lam_p, lam_m, mu_p, mu_m


def m_infinity(theta_te, theta_tr, coefficients):
    """Closed form of the infinite-order sums of SI Eqs. S4.2-S4.4.

    ``M_∞(θ_TE, θ_TR) = e^{μ_+ + iθ_TE} / (1 - e^{λ_+ + iθ_TR}) - e^{μ_- + iθ_TE} / (1 - e^{λ_- + iθ_TR})``
    is the geometric-series form the three closed forms share (the name is
    shorthand used in this package).  ``coefficients`` are
    ``(λ_+, λ_-, μ_+, μ_-)`` from :func:`exponential_coefficients`.  With them

    * ``M_bSSFP        = M_∞(Δω_0 TE,  Δω_0 TR + Δψ)``            (ΔTE = 0, TE_0 = TE)
    * ``M_Σ   (N -> ∞) = M_∞(Δω_0 TE_0, Δω_0 (TR - ΔTE) + Δψ + Δφ)``
    * ``M_Σ|·|(N -> ∞) = M_∞(0, Δφ)``
    """
    lam_p, lam_m, mu_p, mu_m = coefficients
    theta_te, theta_tr = np.asarray(theta_te), np.asarray(theta_tr)
    return (np.exp(mu_p + 1j * theta_te) / (1 - np.exp(lam_p + 1j * theta_tr))
            - np.exp(mu_m + 1j * theta_te) / (1 - np.exp(lam_m + 1j * theta_tr)))
