"""The Pulseq output realises the designed echo positions, echo times and balancing."""

import numpy as np
import pytest

import mess

SMALL = dict(matrix=(32, 8), fov=(0.2, 0.2), dwell=10e-6, rf_type="sinc")


def kspace(seq):
    k_adc, t_adc, t_exc = mess.adc_kspace(seq)
    return k_adc, t_exc, t_adc


@pytest.mark.parametrize("orders", [(0, 0), (-1, -1), (0, -1), (-1, 0), (1, -1), (2, -3), (-3, 2)])
def test_echo_positions_and_times(orders):
    p = mess.MESS2D(orders=orders, **SMALL)
    k_adc, t_exc, t_adc = kspace(p.make_sequence(n_dummy=3))
    n = p.matrix[0]
    kx = k_adc[0].reshape(p.matrix[1], p.n_echoes, n)
    t = t_adc.reshape(p.matrix[1], p.n_echoes, n)
    centre_k = 0.5 * (kx[..., n // 2 - 1] + kx[..., n // 2])
    centre_t = 0.5 * (t[..., n // 2 - 1] + t[..., n // 2]) - t_exc[3:, None]
    for j, k in enumerate(p.orders):
        assert np.allclose(centre_k[:, j], mess.refocusing_moment(k, p.moments) * p.u, atol=1e-6)
        assert np.allclose(centre_t[:, j], p.echo_times[k], atol=1e-9)
    assert np.allclose(np.diff(t_exc), p.tr)


@pytest.mark.parametrize("cls, extra", [(mess.MESS2D, {}),
                                        (mess.MESS3D, dict(matrix=(32, 6, 4), fov=(0.2, 0.2, 0.1)))])
def test_net_moments(cls, extra):
    """Only the readout axis is unbalanced: (G1 + G2 + G3) u per TR on x, zero on y and z."""
    p = cls(orders=(2, -3), **{**SMALL, **extra})
    seq = p.make_sequence(n_dummy=2)
    waves = seq.waveforms_and_times()[0]
    t_exc = seq.waveforms_and_times()[1][0]
    for axis, expected in zip(range(3), (sum(p.moments) * p.u, 0.0, 0.0)):
        t, g = np.real(waves[axis][0]), np.real(waves[axis][1])
        k = np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(t))])
        k_at_rf = np.interp(t_exc, t, k)
        assert np.allclose(np.diff(k_at_rf), expected, atol=1e-6), axis


def test_bssfp_reference_matches_timing():
    p = mess.MESS2D(orders=(2, -3), **SMALL)
    b = p.bssfp()
    assert b.balanced and b.tr == pytest.approx(p.tr) and b.te0 == pytest.approx(p.te0)
    k_adc, t_exc, t_adc = kspace(b.make_sequence())
    n = b.matrix[0]
    t = t_adc.reshape(b.matrix[1], n)
    assert np.allclose(0.5 * (t[:, n // 2 - 1] + t[:, n // 2]) - t_exc, p.te0, atol=1e-9)


def test_te_and_tr_requests():
    p = mess.MESS2D(orders=(2, -3), **SMALL)
    q = mess.MESS2D(orders=(2, -3), tr=p.tr + 1e-3, te=p.te0 + 0.5e-3, **SMALL)
    assert q.tr == pytest.approx(p.tr + 1e-3) and q.te0 == pytest.approx(p.te0 + 0.5e-3)
    with pytest.raises(ValueError):
        mess.MESS2D(orders=(2, -3), tr=p.tr - 1e-3, **SMALL)
    with pytest.raises(ValueError):
        mess.MESS2D(orders=(2, -3), te=p.te0 - 1e-3, **SMALL)


def test_phase_cycling():
    p = mess.MESS2D(orders=(0, -1), **SMALL)
    seq = p.make_sequence(phase_cycle=np.pi / 2, n_dummy=1)
    phases = [seq.get_block(i).rf.phase_offset for i in range(1, len(seq.block_events) + 1)
              if seq.get_block(i).rf is not None]
    assert np.allclose(np.diff(np.unwrap(phases)), np.pi / 2)


def test_paper_protocol_timing():
    """The 3D protocol of the paper: TR 9.83 ms with this PyPulseq, TE_0 = 4.5 ms, 0.75 ms spacing."""
    p = mess.MESS3D(orders=(2, -3), flip_angle=50)
    assert p.moments == (-5, 12, -5)
    assert p.delta_te == pytest.approx(0.75e-3)
    assert p.te0 == pytest.approx(4.5e-3, abs=0.01e-3)
    assert p.scan_time() == pytest.approx(150 * 75 * p.tr)


def test_flip_angle_ramp():
    """A linear ramp over the first dummies; encoded TRs keep the full flip angle."""
    p = mess.MESS2D(orders=(2, -3), **SMALL)
    seq = p.make_sequence(n_dummy=6, ramp=4)
    peaks = [np.max(np.abs(seq.get_block(i).rf.signal)) for i in range(1, len(seq.block_events) + 1)
             if seq.get_block(i).rf is not None]
    assert np.allclose(np.array(peaks[:5]) / peaks[-1], [0.2, 0.4, 0.6, 0.8, 1.0])
    assert np.allclose(peaks[4:], peaks[-1])
    with pytest.raises(ValueError):
        p.make_sequence(n_dummy=2, ramp=3)


def test_equal_phase_cycles_give_identical_files(tmp_path):
    """Last-bit differences of the phase-cycle increment do not change the .seq file."""
    p = mess.MESS2D(orders=(0, -1), **SMALL)
    a, b = np.deg2rad(np.arange(0, 360, 15))[16], np.deg2rad(240)
    p.make_sequence(phase_cycle=a, n_dummy=200).write(str(tmp_path / "a.seq"))
    p.make_sequence(phase_cycle=b, n_dummy=200).write(str(tmp_path / "b.seq"))
    assert (tmp_path / "a.seq").read_text() == (tmp_path / "b.seq").read_text()


def test_sequences_built_in_a_loop_are_identical(tmp_path):
    """Sequences built one after another from one protocol match fresh-object builds.

    Guards against event objects carrying stale PyPulseq registrations from a
    garbage-collected sequence (whose id() Python reuses) into a new one.
    """
    import gc

    def text(seq):
        path = tmp_path / "s.seq"
        seq.write(str(path))
        return path.read_text()

    cycles = np.deg2rad(np.arange(0, 360, 30))
    shared = mess.MESS2D(orders=(2, -3), **SMALL).bssfp()
    for _ in range(2):
        for pc in cycles:
            fresh = mess.MESS2D(orders=(2, -3), **SMALL).bssfp()
            reference = text(fresh.make_sequence(phase_cycle=pc, n_dummy=10))
            assert text(shared.make_sequence(phase_cycle=pc, n_dummy=10)) == reference
            gc.collect()
