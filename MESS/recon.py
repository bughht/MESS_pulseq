"""Echo separation and the two MESS reconstructions.

A MESS readout samples the echoes ``S_k`` of the orders ``[k_1, ..., k_N]``
back to back inside one ADC.  After :func:`split_echoes` and a Fourier
transform per echo, the images are combined with a *virtual phase-cycle*
increment ``Δφ`` chosen at reconstruction time:

``complex_sum``    ``M_Σ(Δφ) = Σ_k S_k exp(i k Δφ)``                       (Eq. 6)
                   reproduces phase-cycled bSSFP, banding included.
``magnitude_sum``  ``M_Σ|·|(Δφ) = Σ_k |S_k| exp(i {k Δφ + [u(k) - 1] π})``  (Eq. 7)
                   drops the off-resonance phase of every echo, so the image
                   has bSSFP-like contrast and no banding.

``Δφ`` enters both sums exactly like the RF phase-cycle increment ``Δψ`` of
bSSFP.  All functions accept an array of ``Δφ`` values; the ``Δφ`` axes then
come first in the output.
"""

import numpy as np

from .theory import polarity


def split_echoes(adc, n_echoes, axis=-1):
    """Split concatenated MESS readouts into their echoes.

    Parameters
    ----------
    adc : array_like
        Raw data whose ``axis`` holds ``n_echoes * N_RO`` samples per readout.
    n_echoes : int
        Number of acquired orders ``N``.

    Returns
    -------
    ndarray
        Shape ``(n_echoes, ..., N_RO)``: entry ``j`` holds the ``j``-th echo in
        acquisition order, i.e. order ``k_{j+1}``.
    """
    data = np.moveaxis(np.asarray(adc), axis, -1)
    n_samples = data.shape[-1]
    if n_samples % n_echoes:
        raise ValueError(f"{n_samples} samples cannot be split into {n_echoes} echoes")
    data = data.reshape(*data.shape[:-1], n_echoes, n_samples // n_echoes)
    return np.moveaxis(data, -2, 0)


def ifft_centered(kspace, axes=(-2, -1)):
    """Centred inverse FFT, k = 0 at index ``N // 2`` of every transformed axis."""
    shifted = np.fft.ifftshift(kspace, axes=axes)
    return np.fft.fftshift(np.fft.ifftn(shifted, axes=axes), axes=axes)


def _weights(orders, dphi):
    k = np.asarray(orders)
    dphi = np.asarray(dphi, dtype=float)
    return np.exp(1j * np.multiply.outer(dphi, k))      # shape dphi.shape + (N,)


def complex_sum(echoes, orders, dphi=0.0):
    """Complex-sum reconstruction ``M_Σ(Δφ) = Σ_k S_k exp(i k Δφ)`` (Eq. 6).

    Parameters
    ----------
    echoes : array_like, shape (N, ...)
        Complex echo images ``S_k`` in the order given by ``orders``.
    orders : sequence of int
        Dephasing order of every echo.
    dphi : float or array_like
        Virtual phase-cycle increment(s) ``Δφ`` [rad].

    Returns
    -------
    ndarray, shape ``np.shape(dphi) + echoes.shape[1:]`` (complex)
    """
    echoes = np.asarray(echoes)
    return np.tensordot(_weights(orders, dphi), echoes, axes=(-1, 0))


def magnitude_sum(echoes, orders, dphi=np.pi):
    """Magnitude-sum reconstruction ``M_Σ|·|(Δφ)`` (Eq. 7).

    ``M_Σ|·|(Δφ) = Σ_k |S_k| exp(i {k Δφ + [u(k) - 1] π})``: the echo
    magnitudes carry no off-resonance phase, and the polarity term restores
    the sign of the negative orders (Eq. 2).  The default ``Δφ = π`` gives the
    bSSFP pass-band contrast.  Arguments as in :func:`complex_sum`.
    """
    echoes = np.abs(np.asarray(echoes))
    weights = _weights(orders, dphi) * polarity(orders)
    return np.tensordot(weights, echoes, axes=(-1, 0))
