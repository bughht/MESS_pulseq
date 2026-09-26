"""The readout diagrams draw every echo where Algorithm 1 refocuses it."""

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import mess  # noqa: E402
from mess import plotting  # noqa: E402

ORDERS = [(0, 0), (-1, -1), (0, -1), (-1, 0), (0, 1), (1, 0), (1, -1), (1, -2), (2, -3), (3, 0),
          (-3, 2)]


@pytest.fixture
def ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


@pytest.mark.parametrize("orders", ORDERS)
def test_echo_in_the_centre_of_its_window(ax, orders):
    """The j-th acquired order refocuses in the middle of the j-th 2u of the readout lobe."""
    acquired = mess.dephasing_orders(*orders)
    d = plotting.readout_diagram(ax, mess.readout_moments(*orders), acquired)
    for j, k in enumerate(acquired):
        assert d["echoes"][k] == pytest.approx(d["edges"][1] + 2 * j + 1)


@pytest.mark.parametrize("orders", ORDERS)
def test_pathway_crosses_zero_at_its_echo(ax, orders):
    """Every drawn pathway line passes through zero dephasing at its echo."""
    acquired = mess.dephasing_orders(*orders)
    d = plotting.readout_diagram(ax, mess.readout_moments(*orders), acquired)
    paths = [a for a in d["artists"] if isinstance(a, matplotlib.lines.Line2D)
             and len(a.get_xdata()) == 4]
    assert len(paths) == len(acquired)
    for k, line in zip(acquired, paths):
        assert np.interp(d["echoes"][k], line.get_xdata(), line.get_ydata()) == pytest.approx(d["zero"])


def test_bssfp_echo_in_the_middle_of_the_readout(ax):
    d = plotting.readout_diagram(ax, mess.BSSFP_MOMENTS)
    assert d["echoes"][0] == pytest.approx(0.5 * (d["edges"][1] + d["edges"][2]))


def test_figures_build():
    fig = plotting.readout_figure([[("bSSFP", None), ("DESS", (0, -1))], [("MESS", (2, -3))]])
    assert len(fig.axes[0].patches) == 3 * 3                 # three lobes per panel
    plt.close(fig)

    protocol = mess.MESS2D(orders=(1, -1), matrix=(16, 16))
    echoes = np.random.default_rng(0).normal(size=(3, 16, 16)) * (1 + 1j)
    fig = plotting.echo_figure(protocol, echoes, np.ones((16, 16), bool))
    assert len(fig.axes) == 1 + 2 * 3 + 1                    # diagram, |S_k| and arg S_k, colour bar
    plt.close(fig)
