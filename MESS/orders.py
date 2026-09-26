"""Dephasing orders and the MESS readout-gradient design (Algorithm 1).

Every SSFP coherence pathway is labelled by an integer *dephasing order* ``k``
(extended-phase-graph notation of Hong et al., MRM 2026).  ``F_k^+`` is the
transverse magnetisation of order ``k`` immediately after an RF pulse.

* In **bSSFP** the readout gradient is fully rewound (``G_1 + G_2 + G_3 = 0``),
  so every order refocuses at the same echo time and the image is the coherent
  sum of all of them.
* An **unbalanced** readout refocuses the orders one after the other, so one
  long ADC samples the echoes ``S_k`` of a contiguous list of orders
  ``[k_1, ..., k_N]``.  DESS (``[0, -1]``), TESS (``[1, 0, -1]``) and MESS
  (for example ``[2, 1, 0, -1, -2, -3]``) are all members of this family.

The readout axis of one TR carries three trapezoids, a prephaser, the readout
lobe and a rewinder, with zeroth moments ``G_1``, ``G_2`` and ``G_3``.  They are
given in units of

    u = N_RO / (2 FOV_RO),

half of the k-space extent that is traversed while one echo is sampled.  With
this unit one echo occupies ``2u`` of the readout lobe, and the net moment per
TR is ``G_1 + G_2 + G_3 = +-2``: one full k-space width, i.e. exactly one cycle
of dephasing across a voxel.  Order ``k`` refocuses where the moment
accumulated since the RF pulse equals ``-k (G_1 + G_2 + G_3)``.
"""

BSSFP_MOMENTS = (-1, 2, -1)
"""Readout moments ``(G_1, G_2, G_3)`` of the fully balanced bSSFP readout."""


def dephasing_orders(k1, kN):
    """Contiguous dephasing orders ``[k_1, ..., k_N]`` in acquisition order.

    >>> dephasing_orders(2, -3)
    [2, 1, 0, -1, -2, -3]
    >>> dephasing_orders(-1, 0)
    [-1, 0]
    """
    k1, kN = int(k1), int(kN)
    step = 1 if kN >= k1 else -1
    return list(range(k1, kN + step, step))


def readout_moments(k1, kN):
    """Readout moments that acquire the orders ``k_1 .. k_N`` (Algorithm 1).

    Parameters
    ----------
    k1, kN : int
        First and last acquired dephasing order.  The echoes are acquired in
        the order ``[k_1, ..., k_N]``; ``k_1 == k_N`` gives a single echo.

    Returns
    -------
    (G1, G2, G3) : tuple of int
        Zeroth moments of the prephaser, readout and rewinder lobes in units of
        ``u = N_RO / (2 FOV_RO)``.

    Examples
    --------
    >>> readout_moments(2, -3)      # six-echo MESS used in the paper
    (-5, 12, -5)
    >>> readout_moments(0, -1)      # DESS
    (-1, 4, -1)
    >>> readout_moments(-1, -1)     # PSIF / SSFP-echo
    (1, 2, -1)
    """
    orders = dephasing_orders(k1, kN)
    n = len(orders)
    if orders[0] < orders[-1]:
        i0, s = -orders[0], -2
    else:
        i0, s = orders[0], 2
    g1 = -2 * i0 - 1
    g2 = 2 * n
    g3 = s + 2 * (i0 - n) + 1
    return g1, g2, g3


def refocusing_moment(k, moments):
    """Readout moment (in units of ``u``) at which order ``k`` refocuses.

    The echo of order ``k`` forms when the moment accumulated since the RF
    pulse equals ``-k (G_1 + G_2 + G_3)``.

    >>> [refocusing_moment(k, (-5, 12, -5)) for k in dephasing_orders(2, -3)]
    [-4, -2, 0, 2, 4, 6]
    """
    return -k * sum(moments)
