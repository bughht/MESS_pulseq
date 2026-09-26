"""Pulseq builders for MESS and for its matched bSSFP reference.

Every repetition (TR) built here has the same block structure::

    | RF + slice/slab | prephaser | TE fill | readout + ADC | rewinder | TR fill |
      gz                G1 (x)                G2 (x)          G3 (x)
                        +PE (y)                               -PE (y)
                        +rephase/PAR (z)                      -PAR/rephase (z)

* The readout (x) axis carries the moments ``G_1, G_2, G_3`` returned by
  :func:`mess.orders.readout_moments`; it is the only unbalanced axis.
* Phase encoding (y) and, in 3D, partition encoding (z) are rewound every TR.
* The slice/slab gradient is rewound exactly: the moment played after the RF
  centre is undone by the prephaser, the moment played before the RF centre is
  undone by the rewinder, so the z axis is balanced even for asymmetric pulses.
* The RF phase advances by the phase-cycle increment ``Δψ`` every TR and the
  ADC phase follows it.

A protocol object fixes all timing when it is created, so ``TR``, ``TE_0``,
``ΔTE`` and every ``TE_k`` can be inspected before the (possibly long) sequence
is built with :meth:`make_sequence`.  :meth:`bssfp` returns the balanced
reference that shares TR, TE, RF pulse and encoding with the MESS protocol.

Units: seconds, metres, flip angle in degrees, phases in radians.
"""

import copy
import warnings
from types import SimpleNamespace

import numpy as np
import pypulseq as pp

from .orders import BSSFP_MOMENTS, dephasing_orders, readout_moments


def default_system():
    """Hardware limits of the 5T scanner used for the experiments in the paper.

    Its ADC accepts any number of samples (``adc_samples_divisor=1``); PyPulseq's
    default divisor of 4 is a Siemens requirement and is kept for any
    ``pp.Opts`` you pass yourself.
    """
    return pp.Opts(max_grad=80, grad_unit="mT/m", max_slew=80, slew_unit="T/m/s",
                   rf_ringdown_time=10e-6, rf_dead_time=100e-6,
                   adc_dead_time=20e-6, grad_raster_time=10e-6, adc_samples_divisor=1)


def flip_angle_ramp(n):
    """Flip-angle scale factors of a linear ramp over ``n`` TRs.

    The ``n`` ramp pulses have ``alpha * (i + 1) / (n + 1)``, i = 0 .. n-1, and
    the pulse after the ramp has the full ``alpha``: a linear flip-angle series
    that damps the oscillating transient of SSFP (Deshpande et al., MRM
    2003;49:151-157).
    """
    return (np.arange(int(n)) + 1.0) / (int(n) + 1.0)


def adc_kspace(seq):
    """k-space position and time of every ADC sample, and the excitation times.

    Returns ``(k_adc, t_adc, t_excitation)`` with ``k_adc`` of shape (3, n).
    Uses PyPulseq's exact piecewise-polynomial integration where available.
    """
    if hasattr(seq, "calculate_kspacePP"):
        k_adc, t_adc, _, _, t_excitation, *_ = seq.calculate_kspacePP()
    else:
        k_adc, _, t_excitation, _, t_adc = seq.calculate_kspace()
    return np.asarray(k_adc), np.asarray(t_adc), np.asarray(t_excitation)


def _fresh(event):
    """A copy of a PyPulseq event that carries no registration of an earlier sequence.

    PyPulseq-Matlab-like caches the library id of an event on the event object,
    keyed by ``id(sequence)``.  Python reuses the id of a garbage-collected
    sequence, so an event object shared between sequences can carry stale ids
    into a new one and silently corrupt it.  Every event is therefore added to a
    sequence as a fresh copy.
    """
    event = copy.copy(event)
    event.__dict__.pop("_pypulseq_sequence_event_cache", None)
    return event


def _contiguous_orders(orders):
    """Validate ``orders`` (``(k_1, k_N)`` or the full contiguous list)."""
    orders = [int(k) for k in orders]
    if len(orders) == 0:
        raise ValueError("orders must not be empty")
    full = dephasing_orders(orders[0], orders[-1])
    if len(orders) > 2 and orders != full:
        raise ValueError(f"orders {orders} are not contiguous; expected {full}")
    return full


class _SSFPProtocol:
    """Timing and gradient design shared by :class:`MESS2D` and :class:`MESS3D`."""

    _name = "SSFP"

    def __init__(self, orders, flip_angle, fov, matrix, dwell, tr, te,
                 rf_duration, time_bw_product, rf_type, thickness, system):
        self.balanced = orders is None
        self.orders = None if self.balanced else _contiguous_orders(orders)
        self.flip_angle = float(flip_angle)
        self.fov = tuple(float(f) for f in fov)
        self.matrix = tuple(int(n) for n in matrix)
        self.dwell = float(dwell)
        self.rf_duration = float(rf_duration)
        self.time_bw_product = float(time_bw_product)
        self.rf_type = rf_type
        self.system = default_system() if system is None else system
        self._thickness = float(thickness)
        self._requested = dict(tr=tr, te=te)
        self._lobe_cache = {}
        self._design(tr, te)

    # ------------------------------------------------------------------ #
    # design: gradients, RF, block durations and echo times
    # ------------------------------------------------------------------ #
    def _design(self, tr, te):
        system = self.system
        raster = system.grad_raster_time
        n_ro, fov_ro = self.matrix[0], self.fov[0]

        # One echo is sampled in n_ro * dwell; it has to be a whole number of
        # gradient rasters so that the flat top stays on the raster for any N.
        self.echo_spacing = n_ro * self.dwell
        if abs(self.echo_spacing / raster - round(self.echo_spacing / raster)) > 1e-6:
            raise ValueError(
                f"matrix[0] * dwell = {self.echo_spacing * 1e6:.2f} us is not a multiple "
                f"of the gradient raster ({raster * 1e6:.2f} us)")

        # Readout axis: G_1, G_2, G_3 in units of u = N_RO / (2 FOV_RO).
        self.u = n_ro / (2 * fov_ro)
        if self.balanced:
            self.moments = BSSFP_MOMENTS
        else:
            self.moments = readout_moments(self.orders[0], self.orders[-1])
        self.n_echoes = 1 if self.balanced else len(self.orders)
        g1, g2, g3 = self.moments

        self.gx_readout = pp.make_trapezoid(
            "x", flat_area=g2 * self.u, flat_time=self.n_echoes * self.echo_spacing,
            system=system)
        rise_area = 0.5 * self.gx_readout.amplitude * self.gx_readout.rise_time
        fall_area = 0.5 * self.gx_readout.amplitude * self.gx_readout.fall_time
        # the ADC samples the flat top only, so the ramps count towards G_1 / G_3
        self.gx_prephaser = pp.make_trapezoid("x", area=g1 * self.u - rise_area, system=system)
        self.gx_rewinder = pp.make_trapezoid("x", area=g3 * self.u - fall_area, system=system)
        self.adc = pp.make_adc(num_samples=self.n_echoes * n_ro, dwell=self.dwell,
                               delay=self.gx_readout.rise_time, system=system)

        # Excitation and exact slice/slab rewinding (see module docstring).
        self.rf, self.gz, gz_rephaser = self._make_excitation()
        self._gz_after_rf = gz_rephaser.area
        self._gz_before_rf = -(self.gz.area + gz_rephaser.area)

        # Cartesian encoding tables, k = (i - N//2) / FOV.
        self.k_pe = (np.arange(self.matrix[1]) - self.matrix[1] // 2) / self.fov[1]
        self.k_par = self._partition_table()

        def lobe(channel, area):
            if abs(area) < 1e-12:
                return 0.0
            return pp.calc_duration(pp.make_trapezoid(channel, area=area, system=system))

        k_pe_max = np.max(np.abs(self.k_pe))
        t_pre = max(pp.calc_duration(self.gx_prephaser), lobe("y", k_pe_max),
                    *(lobe("z", self._gz_after_rf + kz) for kz in self.k_par[[0, -1]]))
        t_rew = max(pp.calc_duration(self.gx_rewinder), lobe("y", k_pe_max),
                    *(lobe("z", self._gz_before_rf - kz) for kz in self.k_par[[0, -1]]))
        t_exc = pp.calc_duration(self.rf, self.gz)
        t_rf_centre = self.rf.delay + pp.calc_rf_center(self.rf)[0]
        t_ro = pp.calc_duration(self.gx_readout, self.adc)

        # Echo times.  Echo j (acquisition order) is centred (j + 1/2) echo
        # spacings after the start of the flat top.  TE_0 is the echo time of
        # order 0 (extrapolated if order 0 is not acquired), so that
        # TE_k = TE_0 - k * delta_te for every acquired order.
        te_first = (t_exc - t_rf_centre) + t_pre + self.adc.delay + 0.5 * self.echo_spacing
        if self.balanced:
            j0, self.delta_te = 0, 0.0
        elif self.orders[0] >= self.orders[-1]:          # descending, as in the paper
            j0, self.delta_te = self.orders[0], self.echo_spacing
        else:                                            # ascending order list
            j0, self.delta_te = -self.orders[0], -self.echo_spacing
        te0_min = te_first + j0 * self.echo_spacing

        d_te = 0.0 if te is None else self._to_raster(te - te0_min)
        if d_te < -1e-9:
            raise ValueError(f"TE_0 = {te * 1e3:.3f} ms is below the minimum "
                             f"{te0_min * 1e3:.3f} ms of this protocol")
        d_te = max(d_te, 0.0)
        self.te0 = te0_min + d_te
        # shortest TR that fits this TE_0 (the TE fill included)
        self.min_tr = t_exc + t_pre + d_te + t_ro + t_rew
        d_tr = 0.0 if tr is None else self._to_raster(tr - self.min_tr)
        if d_tr < -1e-9:
            raise ValueError(f"TR = {tr * 1e3:.3f} ms is below the minimum "
                             f"{self.min_tr * 1e3:.3f} ms of this protocol")
        d_tr = max(d_tr, 0.0)
        self.tr = self.min_tr + d_tr

        if self.balanced:
            self.echo_times = None
        else:
            self.echo_times = {k: self.te0 - k * self.delta_te for k in self.orders}
        self._timing = SimpleNamespace(excitation=t_exc, rf_centre=t_rf_centre,
                                       prephaser=t_pre, te_fill=d_te, readout=t_ro,
                                       rewinder=t_rew, tr_fill=d_tr)

    def _to_raster(self, t):
        raster = self.system.grad_raster_time
        return float(np.round(t / raster) * raster)

    def _make_excitation(self):
        kwargs = dict(flip_angle=np.deg2rad(self.flip_angle), duration=self.rf_duration,
                      slice_thickness=self._thickness, time_bw_product=self.time_bw_product,
                      delay=self.system.rf_dead_time, system=self.system,
                      use="excitation", return_gz=True)
        rf_type = self.rf_type
        if rf_type == "slr" and not hasattr(pp, "make_slr_pulse"):
            # make_slr_pulse ships with the MATLAB-like PyPulseq used for the
            # experiments; stock PyPulseq falls back to an apodised sinc pulse.
            warnings.warn("this PyPulseq has no make_slr_pulse; using a sinc pulse",
                          stacklevel=3)
            rf_type = "sinc"
        if rf_type == "slr":
            return pp.make_slr_pulse(**kwargs, recenter_on_sample=True)
        if rf_type == "sinc":
            return pp.make_sinc_pulse(**kwargs, apodization=0.5)
        raise ValueError("rf_type must be 'slr' or 'sinc'")

    # ------------------------------------------------------------------ #
    # geometry-specific parts
    # ------------------------------------------------------------------ #
    def _partition_table(self):
        return np.zeros(1)

    def _encoding_steps(self):
        """Yield ``(labels, k_pe, k_par)`` for every phase-encoding step."""
        raise NotImplementedError

    def _geometry(self):
        raise NotImplementedError

    def _params(self):
        raise NotImplementedError

    # ------------------------------------------------------------------ #
    # sequence construction
    # ------------------------------------------------------------------ #
    def _lobe(self, channel, area, duration):
        """Trapezoid with the given area and duration (cached), or None for zero area."""
        key = (channel, round(area, 9), duration)
        if key not in self._lobe_cache:
            self._lobe_cache[key] = (None if abs(area) < 1e-12 else
                                     pp.make_trapezoid(channel, area=area, duration=duration,
                                                       system=self.system))
        return self._lobe_cache[key]

    def add_tr(self, seq, i_tr, phase_cycle=0.0, k_pe=0.0, k_par=0.0, adc=True, labels=None,
               flip_scale=1.0, trigger=False):
        """Append one repetition to ``seq``.

        :meth:`make_sequence` is a loop over this method; call it directly to
        build other schedules (for example the steady-state probe of notebook 2).

        Parameters
        ----------
        seq : pypulseq.Sequence
        i_tr : int
            Index of the TR in the sequence; RF and ADC phase are ``i_tr * phase_cycle``.
        phase_cycle : float
            RF phase-cycle increment ``Δψ`` [rad].
        k_pe, k_par : float
            Phase- and partition-encoding moment [1/m], rewound in the same TR.
        adc : bool
            Sample the echoes (False for a dummy TR).
        labels : dict, optional
            Pulseq labels of the ADC block, e.g. ``{"LIN": 3}``.
        flip_scale : float
            Scales the RF amplitude, i.e. the flip angle (flip-angle ramps).
        trigger : bool
            Add a 20 us digital output pulse on ``ext1`` to the ADC block.
        """
        t = self._timing
        # Rounded so that equal phases give identical RF events in the .seq file.
        phase = float(np.round(np.mod(i_tr * phase_cycle, 2 * np.pi), 9))

        rf = _fresh(self.rf)
        rf.phase_offset = phase
        if flip_scale != 1.0:
            rf.signal = self.rf.signal * flip_scale
        seq.add_block(rf, _fresh(self.gz))

        others = [_fresh(g) for g in (self._lobe("y", k_pe, t.prephaser),
                                      self._lobe("z", self._gz_after_rf + k_par, t.prephaser)) if g]
        seq.add_block(*pp.align(right=[_fresh(self.gx_prephaser)], left=others))
        if t.te_fill > 0:
            seq.add_block(pp.make_delay(t.te_fill))

        if adc:
            sampling = _fresh(self.adc)
            sampling.phase_offset = phase
            events = [_fresh(self.gx_readout), sampling]
            events += [pp.make_label(name, "SET", int(value))
                       for name, value in (labels or {}).items()]
            if trigger:
                events.append(pp.make_digital_output_pulse("ext1", duration=20e-6,
                                                           system=self.system))
            seq.add_block(*events)
        else:
            seq.add_block(_fresh(self.gx_readout), pp.make_delay(t.readout))

        others = [_fresh(g) for g in (self._lobe("y", -k_pe, t.rewinder),
                                      self._lobe("z", self._gz_before_rf - k_par, t.rewinder)) if g]
        seq.add_block(*pp.align(left=[_fresh(self.gx_rewinder)], right=others))
        if t.tr_fill > 0:
            seq.add_block(pp.make_delay(t.tr_fill))

    def make_sequence(self, phase_cycle=0.0, n_dummy=0, n_repeats=1, ramp=0, labels=True,
                      trigger=False):
        """Build the Pulseq sequence.

        Parameters
        ----------
        phase_cycle : float
            RF phase-cycle increment ``Δψ`` [rad] added every TR; the ADC phase
            follows the RF phase.
        n_dummy : int
            Dummy TRs played before the first encoded one (identical TRs with
            the ADC off and no phase encoding) to reach the steady state.
        n_repeats : int
            Number of passes through the phase-encoding table.
        ramp : int
            Length of a linear flip-angle ramp played on the first dummy TRs
            (:func:`flip_angle_ramp`); it damps the oscillating part of the
            transient.  Needs ``ramp <= n_dummy``.
        labels : bool
            Add ``LIN`` (and ``PAR`` / ``REP``) labels to the ADC blocks.
        trigger : bool
            Add a 20 us digital output pulse on ``ext1`` to every ADC block.

        Returns
        -------
        pypulseq.Sequence
        """
        if ramp > n_dummy:
            raise ValueError("the flip-angle ramp is played on dummy TRs: need ramp <= n_dummy")
        # Round away last-bit noise (np.deg2rad of an array and of a scalar can
        # differ), so that equal phase cycles always give identical .seq files.
        phase_cycle = float(np.round(phase_cycle, 12))
        seq = pp.Sequence(self.system)
        scales = flip_angle_ramp(ramp)
        i_tr = 0
        for n in range(int(n_dummy)):
            self.add_tr(seq, i_tr, phase_cycle, adc=False,
                        flip_scale=scales[n] if n < ramp else 1.0)
            i_tr += 1
        for rep in range(int(n_repeats)):
            for step_labels, k_pe, k_par in self._encoding_steps():
                if n_repeats > 1:
                    step_labels = dict(step_labels, REP=rep)
                self.add_tr(seq, i_tr, phase_cycle, k_pe, k_par, adc=True,
                            labels=step_labels if labels else None, trigger=trigger)
                i_tr += 1

        seq.set_definition("Name", "bssfp" if self.balanced else "mess")
        seq.set_definition("FOV", list(self._geometry()))
        seq.set_definition("DephasingOrders", "balanced" if self.balanced else self.orders)
        seq.set_definition("ReadoutMoments", list(self.moments))
        seq.set_definition("PhaseCycle", float(phase_cycle))
        seq.set_definition("DummyScans", int(n_dummy))
        if ramp:
            seq.set_definition("FlipAngleRamp", int(ramp))
        ok, errors = seq.check_timing()
        if not ok:
            raise RuntimeError("timing check failed:\n" + "\n".join(map(str, errors[:10])))
        return seq

    def bssfp(self):
        """Balanced SSFP reference with the same TR, TE (= TE_0), RF and encoding."""
        return type(self)(**dict(self._params(), orders=None, tr=self.tr, te=self.te0))

    # ------------------------------------------------------------------ #
    # reporting
    # ------------------------------------------------------------------ #
    @property
    def n_steps(self):
        """Number of encoded TRs in one pass through the encoding table."""
        return self.matrix[1] * len(self.k_par)

    def scan_time(self, n_dummy=0, n_repeats=1):
        """Duration of :meth:`make_sequence` with the same arguments [s]."""
        return (n_dummy + n_repeats * self.n_steps) * self.tr

    def summary(self):
        """Human-readable description of the protocol."""
        g1, g2, g3 = self.moments
        if self.balanced:
            head = f"{self._name} (bSSFP, balanced readout)"
        else:
            head = f"{self._name}  orders {self.orders}"
        lines = [
            head,
            f"  readout moments (G1, G2, G3) = ({g1:+d}, {g2:+d}, {g3:+d}) x u,"
            f"  u = {self.u:.1f} 1/m",
            f"  flip angle {self.flip_angle:g} deg   TR {self.tr * 1e3:.3f} ms"
            f"   TE_0 {self.te0 * 1e3:.3f} ms   delta TE {self.delta_te * 1e3:.3f} ms",
        ]
        if not self.balanced:
            lines.append("  TE_k [ms]: " + "  ".join(
                f"k={k:+d}: {te * 1e3:.3f}" for k, te in self.echo_times.items()))
        geom = " x ".join(f"{f * 1e3:.0f}" for f in self._geometry())
        lines.append(f"  matrix {' x '.join(map(str, self.matrix))}   FOV {geom} mm"
                     f"   scan time {self.scan_time():.1f} s (+ dummies)")
        return "\n".join(lines)

    def __repr__(self):
        return self.summary()


class MESS2D(_SSFPProtocol):
    """Slice-selective 2D MESS (or bSSFP) protocol.

    Parameters
    ----------
    orders : (k_1, k_N) or list of int or None
        Contiguous dephasing orders to acquire, e.g. ``(2, -3)`` for
        ``[2, 1, 0, -1, -2, -3]``.  ``None`` gives the balanced bSSFP readout.
    flip_angle : float
        Flip angle ``alpha`` [deg].
    fov : (float, float)
        Field of view (readout, phase encoding) [m].
    matrix : (int, int)
        Samples per echo ``N_RO`` and phase-encoding lines ``N_PE``.
    slice_thickness : float
        Excited slice thickness [m].
    dwell : float
        ADC dwell time [s]; ``matrix[0] * dwell`` is the echo spacing.
    tr : float or None
        Repetition time [s]; ``None`` uses the minimum.
    te : float or None
        Echo time ``TE_0`` of order 0 [s]; ``None`` uses the minimum.
    rf_duration, time_bw_product : float
        Duration [s] and time-bandwidth product of the excitation pulse.
    rf_type : {'slr', 'sinc'}
        Excitation pulse shape.
    system : pypulseq.Opts or None
        Hardware limits, :func:`default_system` if None.
    """

    _name = "MESS2D"

    def __init__(self, orders=(2, -3), flip_angle=50.0, fov=(200e-3, 200e-3),
                 matrix=(150, 150), slice_thickness=5e-3, dwell=5e-6, tr=None, te=None,
                 rf_duration=1.5e-3, time_bw_product=4.0, rf_type="slr", system=None):
        self.slice_thickness = float(slice_thickness)
        super().__init__(orders, flip_angle, fov, matrix, dwell, tr, te, rf_duration,
                         time_bw_product, rf_type, slice_thickness, system)

    def _encoding_steps(self):
        for i, k in enumerate(self.k_pe):
            yield {"LIN": i}, k, 0.0

    def _geometry(self):
        return self.fov[0], self.fov[1], self.slice_thickness

    def _params(self):
        return dict(orders=self.orders, flip_angle=self.flip_angle, fov=self.fov,
                    matrix=self.matrix, slice_thickness=self.slice_thickness,
                    dwell=self.dwell, tr=self._requested["tr"], te=self._requested["te"],
                    rf_duration=self.rf_duration, time_bw_product=self.time_bw_product,
                    rf_type=self.rf_type, system=self.system)


class MESS3D(_SSFPProtocol):
    """Slab-selective 3D MESS (or bSSFP) protocol, as used for the scans in the paper.

    Partitions are encoded along the slab axis (z) in the outer loop and phase
    encoding (y) in the inner loop; ADC blocks carry ``LIN`` and ``PAR`` labels.

    Parameters
    ----------
    fov : (float, float, float)
        Field of view (readout, phase, partition) [m].
    matrix : (int, int, int)
        ``N_RO`` samples per echo, ``N_PE`` lines and ``N_PAR`` partitions.
    slab_thickness : float or None
        Excited slab [m]; by default 80 % of the partition FOV, which keeps the
        slab-profile transition bands inside the encoded FOV.

    The remaining parameters are those of :class:`MESS2D`.
    """

    _name = "MESS3D"

    def __init__(self, orders=(2, -3), flip_angle=50.0, fov=(200e-3, 200e-3, 100e-3),
                 matrix=(150, 150, 75), slab_thickness=None, dwell=5e-6, tr=None, te=None,
                 rf_duration=1.5e-3, time_bw_product=4.0, rf_type="slr", system=None):
        self.slab_thickness = 0.8 * fov[2] if slab_thickness is None else float(slab_thickness)
        super().__init__(orders, flip_angle, fov, matrix, dwell, tr, te, rf_duration,
                         time_bw_product, rf_type, self.slab_thickness, system)

    def _partition_table(self):
        n_par = self.matrix[2]
        return (np.arange(n_par) - n_par // 2) / self.fov[2]

    def _encoding_steps(self):
        for p, k_par in enumerate(self.k_par):
            for i, k_pe in enumerate(self.k_pe):
                yield {"LIN": i, "PAR": p}, k_pe, k_par

    def _geometry(self):
        return self.fov

    def _params(self):
        return dict(orders=self.orders, flip_angle=self.flip_angle, fov=self.fov,
                    matrix=self.matrix, slab_thickness=self.slab_thickness,
                    dwell=self.dwell, tr=self._requested["tr"], te=self._requested["te"],
                    rf_duration=self.rf_duration, time_bw_product=self.time_bw_product,
                    rf_type=self.rf_type, system=self.system)
