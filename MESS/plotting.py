"""Figures and animations used by the notebooks (matplotlib only).

Tissue colours are the first three slots of a colour-blind-safe categorical
palette; curves always carry a legend or a direct label as well.  Magnitude
images use grey, phase images ``jet`` from -π to π, off-resonance maps the
diverging ``RdBu_r`` map centred on 0 Hz.  Readout diagrams follow Fig. 1 of
the paper: prephaser orange, readout green, rewinder blue.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

from .orders import BSSFP_MOMENTS, dephasing_orders, readout_moments, refocusing_moment

TISSUE_COLORS = {"WM": "#2a78d6", "GM": "#eb6834", "CSF": "#1baf7a"}
LOBE_COLORS = ("#f4b183", "#a9d18e", "#92cddc")      # G1, G2, G3 as in Fig. 1 of the paper
INK, INK_SECONDARY, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"


def use_style():
    """Quiet axes and a light surface for all notebook figures."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
        "axes.labelcolor": INK_SECONDARY, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.6,
        "font.size": 9, "axes.titlesize": 9.5, "legend.frameon": False,
        "lines.linewidth": 1.5, "image.interpolation": "nearest",
        "figure.dpi": 100,
    })


def show(ax, image, vmin=0.0, vmax=None, cmap="gray", title=None):
    """Display a map indexed ``[x, y]`` with anterior at the top."""
    handle = ax.imshow(np.asarray(image).T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_axis_off()
    if title:
        ax.set_title(title)
    return handle


def show_phase(ax, image, title=None):
    """Display the phase of a complex map, -π (blue) to π (red)."""
    return show(ax, np.angle(image), vmin=-np.pi, vmax=np.pi, cmap="jet", title=title)


def phase_colorbar(fig, handle, cax=None, ax=None, **kw):
    """Colour bar of a :func:`show_phase` map with ticks at -π, 0 and π."""
    bar = fig.colorbar(handle, cax=cax, ax=ax, ticks=[-np.pi, 0, np.pi], **kw)
    bar.ax.set_yticklabels(["\N{MINUS SIGN}π", "0", "π"])
    bar.outline.set_visible(False)
    return bar


# ---------------------------------------------------------------------- #
# sequence design
# ---------------------------------------------------------------------- #
def _signed(n, plus=True):
    """An integer with a typographic minus sign (and a plus sign if ``plus``)."""
    return (f"{n:+d}" if plus else f"{n:d}").replace("-", "\N{MINUS SIGN}")


def _lobes(moments, max_amplitude):
    """Lobe edges and amplitudes of one TR; the readout lobe has unit amplitude."""
    amplitude = np.array([np.sign(g) * (1.0 if i == 1 else min(abs(g), max_amplitude))
                          for i, g in enumerate(moments)])
    edges = np.concatenate([[0.0], np.cumsum(np.abs(moments) / np.abs(amplitude))])
    return edges, amplitude


def readout_diagram(ax, moments, orders=None, name=None, origin=(0.0, 0.0), max_amplitude=2.25,
                    height=1.0, epg_scale=0.3, gap=0.8, fontsize=9.0):
    """Readout gradient and extended phase graph of one TR, drawn like Fig. 1 of the paper.

    Top: prephaser, readout and rewinder, labelled with their moments
    ``G_1, G_2, G_3`` in units of ``u``.  The readout lobe has unit amplitude,
    so the moment grows by one ``u`` per unit of time; prephaser and rewinder
    are at most ``max_amplitude`` high and as long as their moments require.
    An amplitude of 1 is drawn ``height`` data units high.

    Bottom: pathway ``F_k`` enters the TR with the dephasing ``k (G_1+G_2+G_3)``
    (drawn ``epg_scale`` data units per ``u``), follows the accumulated moment
    and leaves the TR as ``F_{k+1}``.  It refocuses (●, "Echo k") where it
    crosses zero inside the readout lobe.  ``orders=None`` draws bSSFP, whose
    pathways all coincide.

    ``name`` (e.g. "DESS") is written to the left, above the list of acquired
    pathways.  ``origin`` is the start of the prephaser on the gradient
    baseline.  Returns a dict with the lobe ``edges``, the ``echoes``
    (order -> time), the height ``zero`` of the zero-dephasing line and the
    ``artists`` drawn.
    """
    x0, y0 = origin
    edges, amplitude = _lobes(moments, max_amplitude)
    x = x0 + edges
    s = sum(moments)
    ks = [0] if orders is None else list(orders)
    paths = np.array([np.cumsum([0, *moments]) + k * s for k in ks], float)
    artists = []

    def note(text, xy, offset, **kw):
        artists.append(ax.annotate(text, xy, xytext=offset, textcoords="offset points",
                                   annotation_clip=False, **kw))

    for i, (g, a, color) in enumerate(zip(moments, amplitude, LOBE_COLORS)):
        artists.append(ax.add_patch(Rectangle((x[i], y0), x[i + 1] - x[i], a * height,
                                              facecolor=color, edgecolor=INK, lw=0.9, zorder=2)))
        note(_signed(g), (0.5 * (x[i] + x[i + 1]), y0 + 0.5 * a * height), (0, 0),
             ha="center", va="center", fontsize=fontsize + 1, zorder=3)
    top = y0 + max(amplitude.max(), 0) * height
    zero = y0 + min(amplitude.min(), 0) * height - gap - paths.max() * epg_scale
    bottom = zero + paths.min() * epg_scale

    artists += ax.plot([x[0], x[-1]], [zero, zero], color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    for k, p in zip(ks, paths):
        y = zero + p * epg_scale
        artists += ax.plot(x, y, color=INK, lw=1.1, zorder=3)
        ends = (r"$\Sigma F_n$",) * 2 if orders is None else (f"$F_{{{k}}}$", f"$F_{{{k + 1}}}$")
        note(ends[0], (x[0], y[0]), (-4, 0), ha="right", va="center", fontsize=fontsize - 0.5)
        note(ends[1], (x[-1], y[-1]), (4, 0), ha="left", va="center", fontsize=fontsize - 0.5)

    echoes = {}
    for k in ks:
        t = x[1] + refocusing_moment(k, moments) - moments[0]
        echoes[k] = t
        artists += ax.plot([t, t], [bottom - 2 * epg_scale, top], color=INK, lw=0.7, ls=":",
                           zorder=1)
        artists += ax.plot([t], [zero], "o", ms=4, color=INK, zorder=4)
        note("Echo\nbSSFP" if orders is None else f"Echo {_signed(k, plus=False)}", (t, top),
             (0, 3), ha="center", va="bottom", fontsize=fontsize - 2, linespacing=1.0)

    if name:
        acquired = (r"$\left[\sum_{-\infty}^{+\infty} F_n\right]$" if orders is None else
                    "$[" + ",\\,".join(f"F_{{{k}}}" for k in ks) + "]$")
        note(name, (x[0], y0), (-30, 2), ha="right", va="bottom", fontsize=fontsize + 2.5)
        note(acquired, (x[0], y0), (-30, -2), ha="right", va="top", fontsize=fontsize + 2.5)
    return dict(edges=x, echoes=echoes, zero=zero, artists=artists)


def _diagram_extent(unit, draw):
    """Bounding box, in data units, of what ``draw(ax)`` puts on an axes of ``unit`` inches per unit."""
    scratch = Figure(figsize=(80 * unit, 80 * unit))
    renderer = FigureCanvasAgg(scratch).get_renderer()
    ax = scratch.add_axes([0, 0, 1, 1])
    ax.set_xlim(-40, 40)
    ax.set_ylim(-40, 40)
    artists = draw(ax)["artists"]
    box = Bbox.union([a.get_window_extent(renderer) for a in artists])
    return box.transformed(ax.transData.inverted())


def _diagram_of(name, orders, **kw):
    """``draw(ax, origin)`` for the readout diagram of ``orders = (k1, kN)``, or bSSFP if None."""
    if orders is None:
        return lambda ax, origin=(0.0, 0.0): readout_diagram(ax, BSSFP_MOMENTS, None, name, origin, **kw)
    moments, acquired = readout_moments(*orders), dephasing_orders(*orders)
    return lambda ax, origin=(0.0, 0.0): readout_diagram(ax, moments, acquired, name, origin, **kw)


def readout_figure(rows, unit=0.25, col_gap=1.5, row_gap=1.2, pad=0.3, **kw):
    """Readout diagrams of several sequences on one scale, arranged like Fig. 1 of the paper.

    ``rows`` is a list of rows, each a list of ``(name, orders)`` with
    ``orders = (k1, kN)``, or ``None`` for bSSFP.  ``unit`` is the length in
    inches of one unit of time (and of moment).  The panels of every row are
    spread over the full width; the rows with the most panels share columns,
    in which the prephasers start at the same time.  ``kw`` goes to
    :func:`readout_diagram`.
    """
    draws = [[_diagram_of(name, orders, **kw) for name, orders in row] for row in rows]
    boxes = [[_diagram_extent(unit, d) for d in row] for row in draws]
    n = max(len(row) for row in boxes)
    full = [row for row in boxes if len(row) == n]
    left = [max(-row[c].x0 for row in full) for c in range(n)]      # origin to column edges
    right = [max(row[c].x1 for row in full) for c in range(n)]
    used = [sum(left) + sum(right) if len(row) == n else sum(b.width for b in row) for row in boxes]
    width = max(u + col_gap * (len(row) - 1) for u, row in zip(used, boxes)) + 2 * pad

    placed, y = [], -pad
    for row_draws, row_boxes, row_used in zip(draws, boxes, used):
        baseline = y - max(b.y1 for b in row_boxes)
        m = len(row_boxes)
        step = (width - 2 * pad - row_used) / (m - 1) if m > 1 else 0.0
        x = pad if m > 1 else (width - row_used) / 2
        for c, (d, b) in enumerate(zip(row_draws, row_boxes)):
            if m == n:
                placed.append((d, (x + left[c], baseline)))
                x += left[c] + right[c] + step
            else:
                placed.append((d, (x - b.x0, baseline)))
                x += b.width + step
        y = baseline + min(b.y0 for b in row_boxes) - row_gap
    height = -(y + row_gap) + pad

    fig = plt.figure(figsize=(width * unit, height * unit))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(-height, 0)
    ax.set_axis_off()
    for d, origin in placed:
        d(ax, origin)
    return fig


def echo_figure(protocol, echoes, mask, image_size=1.35, spacing=0.1, fontsize=10.0, pad=0.12,
                **kw):
    """The echo images of one MESS acquisition under the readout diagram of its TR.

    Column ``j`` shows ``|S_k|`` (each on its own grey scale, maximum printed)
    and ``arg S_k`` of the ``j``-th acquired order ``k``.  Every echo occupies
    ``2 u`` of the readout lobe, which fixes the horizontal scale: each column
    sits under the echo that it was reconstructed from.

    ``echoes`` has the orders along the first axis (:func:`mess.simulation.reconstruct`);
    ``mask`` selects the object for the grey scales and the phase maps.  ``kw``
    goes to :func:`readout_diagram`.
    """
    unit = (image_size + spacing) / 2                           # inches per u
    kw = dict(dict(height=0.3, epg_scale=0.11, gap=0.35, fontsize=fontsize), **kw)
    draw = lambda ax, origin=(0.0, 0.0): readout_diagram(ax, protocol.moments, protocol.orders,
                                                         origin=origin, **kw)
    box = _diagram_extent(unit, draw)
    edges, _ = _lobes(protocol.moments, kw.get("max_amplitude", 2.25))
    title, bar = 0.32, 0.75                                     # inches
    width = (max(box.x1, edges[2] + bar / unit) - box.x0) * unit + 2 * pad
    strip = box.height * unit
    height = pad + strip + 2 * (title + image_size) + pad
    x_left = box.x0 - pad / unit                                # data x of the left figure edge

    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0, (height - pad - strip) / height, 1, strip / height])
    ax.set_xlim(x_left, x_left + width / unit)
    ax.set_ylim(box.y0, box.y1)
    ax.set_axis_off()
    draw(ax)

    def panel(j, row, text):
        left = (edges[1] + 2 * j - x_left) * unit + spacing / 2
        bottom = height - pad - strip - (row + 1) * (title + image_size)
        a = fig.add_axes([left / width, bottom / height, image_size / width, image_size / height])
        a.set_title(text, fontsize=fontsize - 1, pad=4)
        return a

    for j, k in enumerate(protocol.orders):
        top = np.percentile(np.abs(echoes[j])[mask], 99.5)
        ax_m = panel(j, 0, f"$|S_{{{k:+d}}}|$, TE = {protocol.echo_times[k] * 1e3:.2f} ms")
        show(ax_m, np.abs(echoes[j]), vmax=top)
        ax_m.text(0.03, 0.03, f"max {top:.3f}", transform=ax_m.transAxes, color="w",
                  fontsize=fontsize - 2)
        ax_p = panel(j, 1, f"$\\arg S_{{{k:+d}}}$")
        handle = show_phase(ax_p, np.where(mask, echoes[j], np.nan))
    pos = ax_p.get_position()
    cax = fig.add_axes([pos.x1 + 0.2 / width, pos.y0, 0.11 / width, pos.height])
    phase_colorbar(fig, handle, cax=cax).ax.tick_params(labelsize=fontsize - 1)
    return fig


def sequence_diagram(protocol, seq, i_tr=None, axes=None):
    """RF, gradients and ADC of one TR of ``seq``, with every ``TE_k`` marked.

    ``i_tr`` defaults to the first encoded TR, whose phase-encoding lobes are
    the largest.
    """
    if i_tr is None:
        i_tr = 0
    t0, t1 = i_tr * protocol.tr, (i_tr + 1) * protocol.tr
    waves, t_exc, _, t_adc, _, _ = seq.waveforms_and_times(append_RF=True, time_range=[t0, t1])
    if axes is None:
        _, axes = plt.subplots(4, 1, figsize=(8, 4.2), sharex=True,
                               gridspec_kw=dict(height_ratios=[1, 1, 1, 1]))
    names = ["RF", "$G_x$ (readout)", "$G_y$ (phase)", "$G_z$ (slice)"]
    channels = [3, 0, 1, 2]
    t_rf = t_exc[0][(t_exc[0] >= t0) & (t_exc[0] < t1)][0]
    for ax, name, ch in zip(axes, names, channels):
        w = waves[ch]
        if w.size:
            t = (np.real(w[0]) - t_rf) * 1e3
            y = np.abs(w[1]) if ch == 3 else w[1] / 1e3
            ax.plot(t, y, color=INK, lw=1.0)
            ax.fill_between(t, 0, y, color="#b7d3f6", lw=0)
        ax.axhline(0, color=AXIS, lw=0.6)
        ax.set_ylabel(name, rotation=0, ha="right", va="center")
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
    adc = t_adc[(t_adc >= t0) & (t_adc < t1)]
    axes[1].axvspan((adc[0] - t_rf) * 1e3, (adc[-1] - t_rf) * 1e3, color="#f0efec", lw=0, zorder=0)
    if protocol.echo_times:
        for k, te in protocol.echo_times.items():
            axes[1].axvline(te * 1e3, color=MUTED, lw=0.7, ls=":")
            axes[1].text(te * 1e3, 1.0, f"{k:+d}", transform=axes[1].get_xaxis_transform(),
                         ha="center", va="bottom", fontsize=7.5, color=INK_SECONDARY)
    axes[-1].set_xlabel("time after the RF centre [ms]")
    axes[-1].set_xlim((t0 - t_rf) * 1e3, (t1 - t_rf) * 1e3)
    return axes


# ---------------------------------------------------------------------- #
# phase-cycle animation
# ---------------------------------------------------------------------- #
def phase_cycle_animation(phases, left, right, left_title, right_title, curves=(),
                          rois=None, roi_values=None, legend=None, vmax=None,
                          ylabel="signal / $M_0$", xlabel="phase-cycle increment", path=None,
                          fps=5, dpi=90):
    """Two image series side by side under a signal profile with a moving cursor.

    Parameters
    ----------
    phases : (F,) array
        Phase-cycle increment of every frame [rad].
    left, right : (F, nx, ny) arrays
        Magnitude images shown in frame ``f``.
    left_title, right_title : callable
        ``f(phase_deg) -> str``, the title of each image panel.
    curves : sequence of dict
        Static profile curves: ``dict(x=rad, y=values, color=..., ls=..., label=...)``.
    rois : dict, optional
        ``name -> (ix, iy)`` voxel positions, marked on both images.
    roi_values : dict, optional
        ``name -> (left_values, right_values)``, each of shape (F,): the
        simulated ROI values drawn as markers that accumulate frame by frame.
    legend : (str, str, str, str), optional
        Labels of the solid curves, the dashed curves, the filled markers (left
        images) and the open markers (right images).  Tissues are identified
        by colour and by the curves' ``direct_label``.
    path : str, optional
        Write the animation to this GIF file.
    """
    use_style()
    phases = np.asarray(phases)
    deg = np.rad2deg(phases)
    vmax = np.percentile(np.asarray(left), 99.5) if vmax is None else vmax
    fig = plt.figure(figsize=(7.6, 6.9))
    grid = GridSpec(2, 2, height_ratios=[1.0, 1.55], hspace=0.28, wspace=0.04,
                    left=0.09, right=0.98, top=0.93, bottom=0.03)
    ax_p = fig.add_subplot(grid[0, :])
    ax_l, ax_r = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])

    for c in curves:
        ax_p.plot(np.rad2deg(c["x"]), c["y"], color=c.get("color", INK),
                  ls=c.get("ls", "-"), lw=c.get("lw", 1.5), label=c.get("label"))
        if c.get("direct_label"):
            i = int(np.argmax(c["y"]))
            ax_p.annotate(c["direct_label"], (np.rad2deg(c["x"][i]), c["y"][i]), xytext=(4, 4),
                          textcoords="offset points", color=INK_SECONDARY, fontsize=8)
    ax_p.set_xlim(0, 360)
    ax_p.set_xticks(np.arange(0, 361, 45))
    ax_p.set_xticklabels([f"{d}°" for d in range(0, 361, 45)])
    ax_p.set_xlabel(xlabel)
    ax_p.set_ylabel(ylabel)
    ax_p.set_ylim(bottom=0)
    ax_p.grid(True, axis="y")
    cursor = ax_p.axvline(deg[0], color=INK, lw=1.0)

    trails = {}
    if roi_values:
        for name, (vl, vr) in roi_values.items():
            color = TISSUE_COLORS.get(name, INK)
            trails[name] = (
                ax_p.plot([], [], "o", ms=5.5, color=color, mec=SURFACE, mew=0.8, zorder=5)[0],
                ax_p.plot([], [], "s", ms=5.5, mfc="none", mec=color, mew=1.3, zorder=6)[0])
    if legend:
        solid, dashed, filled, hollow = legend
        handles = [Line2D([], [], color=INK, lw=1.5, label=solid),
                   Line2D([], [], color=INK, lw=1.1, ls="--", label=dashed),
                   Line2D([], [], ls="none", marker="o", ms=5.5, color=INK, label=filled),
                   Line2D([], [], ls="none", marker="s", ms=5.5, mfc="none", mec=INK, mew=1.3,
                          label=hollow)]
        ax_p.legend(handles=handles, loc="center left", bbox_to_anchor=(0.38, 0.68), ncol=2,
                    fontsize=7.5, handlelength=2.2)

    im_l = show(ax_l, left[0], vmax=vmax)
    im_r = show(ax_r, right[0], vmax=vmax)
    t_l = ax_l.set_title(left_title(deg[0]))
    t_r = ax_r.set_title(right_title(deg[0]))
    if rois:
        for name, (ix, iy) in rois.items():
            for ax in (ax_l, ax_r):
                ax.plot([ix], [iy], "s", ms=7, mfc="none", mec=TISSUE_COLORS.get(name, "w"), mew=1.5)

    def update(f):
        im_l.set_data(np.asarray(left[f]).T)
        im_r.set_data(np.asarray(right[f]).T)
        t_l.set_text(left_title(deg[f]))
        t_r.set_text(right_title(deg[f]))
        cursor.set_xdata([deg[f], deg[f]])
        for name, (m_l, m_r) in trails.items():
            vl, vr = roi_values[name]
            m_l.set_data(deg[: f + 1], np.asarray(vl)[: f + 1])
            m_r.set_data(deg[: f + 1], np.asarray(vr)[: f + 1])
        return [im_l, im_r, t_l, t_r, cursor]

    anim = animation.FuncAnimation(fig, update, frames=len(phases), interval=1000 / fps, blit=False)
    if path is not None:
        anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    return fig, anim
