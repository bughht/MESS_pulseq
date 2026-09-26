"""Algorithm 1 against every (G1, G2, G3) listed in the paper (Fig. 1B, Fig. S3)."""

import pytest

from mess import BSSFP_MOMENTS, dephasing_orders, readout_moments, refocusing_moment

PAPER = {
    (2, -3): (-5, 12, -5),    # six-echo MESS, Fig. 1B
    (0, -1): (-1, 4, -1),     # DESS, Fig. S3B
    (1, -1): (-3, 6, -1),     # TESS, Fig. S3C
    (1, -2): (-3, 8, -3),     # QESS, Fig. S3D
    (2, -2): (-5, 10, -3),    # five-echo MESS, Fig. S3E
    (-1, 0): (-3, 4, -3),     # DESS variants of the original Fig. 1
    (0, 1): (-1, 4, -5),
    (1, 0): (-3, 4, 1),
}


@pytest.mark.parametrize("orders, moments", PAPER.items())
def test_paper_examples(orders, moments):
    assert readout_moments(*orders) == moments


@pytest.mark.parametrize("k1, kN", [(k1, kN) for k1 in range(-4, 5) for kN in range(-4, 5)])
def test_echo_condition(k1, kN):
    """Every order refocuses in the middle of its own 2u-wide readout window."""
    g1, g2, g3 = readout_moments(k1, kN)
    orders = dephasing_orders(k1, kN)
    assert g2 == 2 * len(orders)
    assert abs(g1 + g2 + g3) == 2                    # one k-space width per TR
    for j, k in enumerate(orders):
        assert refocusing_moment(k, (g1, g2, g3)) == g1 + 2 * j + 1


def test_balanced():
    assert sum(BSSFP_MOMENTS) == 0
