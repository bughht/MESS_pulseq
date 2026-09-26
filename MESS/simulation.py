"""MR-zero helpers used by the tutorial notebooks.

Requires ``MRzeroCore`` (and PyTorch): ``pip install MRzeroCore``.  MR-zero
simulates the very same ``.seq`` files that are played out on the scanner.

Conventions
-----------
* MR-zero writes free precession as ``exp(-i ω t)``, the paper (Eq. 3) as
  ``exp(+i ω t)``.  :func:`reconstruct` therefore returns the complex conjugate
  of the images, so that echo phases follow Eq. 4b exactly.  Magnitudes, and
  hence every banding pattern, are unaffected.
* Maps are indexed ``[x, y]`` (x: readout, y: phase encoding) with voxel ``i``
  at ``(i - N/2) * FOV/N``, which is also the pixel grid of :func:`reconstruct`.
* ``off_resonance`` maps are in Hz: ``Δω_0 = 2π * off_resonance``.
"""

import hashlib
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._lock import simulation_lock
from .recon import ifft_centered, split_echoes

TABLE_S1 = {
    # tissue: (PD [a.u.], T1 [s], T2 [s], T2' [s]) -- Table S1 of the paper
    "GM": (0.800, 1.560, 0.083, 0.320),
    "WM": (0.700, 0.830, 0.075, 0.180),
    "Fat": (1.000, 0.370, 0.125, 0.012),
    "CSF": (1.000, 4.160, 1.650, 0.059),
}

BRAIN_URL = ("https://github.com/MRsources/MRzero-Core/raw/main/documentation/"
             "playground_mr0/numerical_brain_cropped.mat")



@dataclass
class Phantom:
    """A 2D numerical phantom; every map is indexed ``[x, y]``."""

    pd: np.ndarray
    t1: np.ndarray
    t2: np.ndarray
    t2p: np.ndarray
    off_resonance: np.ndarray          # Hz
    tissue: np.ndarray                 # tissue name of every voxel ("" = background)
    fov: tuple = (0.2, 0.2)            # m
    thickness: float = 5e-3            # m
    diffusion: np.ndarray = None       # 1e-9 m^2/s (MR-zero units); None = no diffusion
    positions: np.ndarray = None       # x index of every probe voxel (probe_phantom)

    def mask(self, name):
        """Boolean mask of one tissue."""
        return self.tissue == name

    @property
    def shape(self):
        return self.pd.shape

    def to_mrzero(self):
        """The phantom as an ``MRzeroCore.VoxelGridPhantom``."""
        import MRzeroCore as mr0
        import torch

        nx, ny = self.shape

        def t(a):
            return torch.as_tensor(np.asarray(a, dtype=np.float32))[:, :, None]

        d = np.zeros(self.shape) if self.diffusion is None else self.diffusion
        dx, dy = self.fov[0] / nx, self.fov[1] / ny
        affine = torch.tensor([[dx * 1e3, 0, 0, -self.fov[0] / 2 * 1e3],
                               [0, dy * 1e3, 0, -self.fov[1] / 2 * 1e3],
                               [0, 0, self.thickness * 1e3, -self.thickness / 2 * 1e3]],
                              dtype=torch.float32)
        ones = torch.ones(1, nx, ny, 1)
        return mr0.VoxelGridPhantom(
            t(self.pd), t(self.t1), t(self.t2), t(self.t2p), t(d),
            t(self.off_resonance), ones.clone(), ones.clone(),
            torch.tensor([self.fov[0], self.fov[1], self.thickness]), affine)

    def fingerprint(self):
        h = hashlib.sha1()
        for a in (self.pd, self.t1, self.t2, self.t2p, self.off_resonance):
            h.update(np.ascontiguousarray(a, dtype=np.float32).tobytes())
        if self.diffusion is not None:
            h.update(np.ascontiguousarray(self.diffusion, dtype=np.float32).tobytes())
        h.update(repr((self.fov, self.thickness, self.shape)).encode())
        return h.hexdigest()


def _tissue_maps(tissue, tissues=TABLE_S1):
    t1, t2, t2p = (np.ones(tissue.shape) for _ in range(3))
    for name, (_, T1, T2, T2p) in tissues.items():
        m = tissue == name
        t1[m], t2[m], t2p[m] = T1, T2, T2p
    return t1, t2, t2p


def brain_phantom(n=150, off_resonance_max=250.0, cache_dir=".cache"):
    """MR-zero's numerical brain with the tissue parameters of Table S1.

    Voxels are labelled WM / GM / CSF from the phantom's own T1 map and receive
    T1, T2 and T2' of Table S1; the proton density map is kept for anatomical
    detail.  The phantom's measured B0 map (a frontal susceptibility hot spot)
    is scaled to a peak of ``off_resonance_max`` Hz to emulate high-field
    inhomogeneity.  Diffusion is off, as in the signal model of the paper.
    """
    import MRzeroCore as mr0

    path = Path(cache_dir) / "numerical_brain_cropped.mat"
    if not path.exists():
        from urllib.request import urlretrieve
        path.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(BRAIN_URL, path)
    ph = mr0.VoxelGridPhantom.load_mat(str(path)).interpolate(n, n, 1)
    pd = ph.PD[:, :, 0].numpy().astype(float)
    t1_mr0 = ph.T1[:, :, 0].numpy()
    b0 = ph.B0[:, :, 0].numpy().astype(float)

    brain = pd > 0.05
    tissue = np.full(pd.shape, "", dtype=object)
    tissue[brain & (t1_mr0 < 1.0)] = "WM"
    tissue[brain & (t1_mr0 >= 1.0) & (t1_mr0 < 2.5)] = "GM"
    tissue[brain & (t1_mr0 >= 2.5)] = "CSF"
    pd = np.where(brain, pd, 0.0)
    t1, t2, t2p = _tissue_maps(tissue)
    off_resonance = np.where(brain, b0 * off_resonance_max / np.abs(b0[brain]).max(), 0.0)
    return Phantom(pd, t1, t2, t2p, off_resonance, tissue, fov=(0.2, 0.2))


def probe_phantom(entries, n=150, fov=0.2):
    """A single row of isolated voxels, one per ``(tissue, off_resonance_Hz)`` entry.

    Acquired with a one-line protocol (``matrix=(n, 1)``), a 1D FFT of every
    readout returns the signal of each voxel separately, TR by TR, which is
    how the notebooks follow the approach to the steady state.
    """
    shape = (n, 1)
    tissue = np.full(shape, "", dtype=object)
    pd, df = np.zeros(shape), np.zeros(shape)
    positions = np.linspace(0, n, len(entries) + 2)[1:-1].round().astype(int)
    for x, (name, f) in zip(positions, entries):
        tissue[x, 0], pd[x, 0], df[x, 0] = name, TABLE_S1[name][0], f
    t1, t2, t2p = _tissue_maps(tissue)
    return Phantom(pd, t1, t2, t2p, df, tissue, fov=(fov, fov / n), positions=positions)


def simulate(seq, phantom, cache_dir=None, device=None, max_states=200, min_signal=1e-5,
             verbose=True):
    """Simulate a PyPulseq sequence on a :class:`Phantom` with MR-zero.

    Returns the complex ADC samples of the whole sequence, in acquisition
    order, as a 1D array.  With ``cache_dir`` the result is stored under a hash
    of the ``.seq`` file, the phantom and the simulation settings, and later
    calls with the same inputs load it instead of simulating again.

    ``max_states`` and ``min_signal`` control MR-zero's pruning of the phase
    distribution graph.  ``min_signal = 1e-5`` keeps the weak outer pathways of
    long-T2 tissue (CSF); with 1e-3 the orders ``|k| >= 2`` are lost entirely.

    Only one simulation runs at a time on this machine (:func:`simulation_lock`);
    cache hits do not wait for it.
    """
    import MRzeroCore as mr0
    import torch

    with tempfile.TemporaryDirectory() as tmp:
        seq_file = Path(tmp) / "sequence.seq"
        seq.write(str(seq_file))
        text = seq_file.read_text()
        text = "\n".join(line for line in text.splitlines() if not line.startswith("# Created"))
        key = hashlib.sha1((text + phantom.fingerprint()
                            + repr((max_states, min_signal))).encode()).hexdigest()[:16]
        cached = None if cache_dir is None else Path(cache_dir) / f"sim_{key}.npy"
        if cached is not None and cached.exists():
            return np.load(cached)
        s = mr0.Sequence.import_file(str(seq_file))

    with simulation_lock(verbose=verbose):
        if cached is not None and cached.exists():              # simulated while we waited
            return np.load(cached)
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        start = time.time()
        data = phantom.to_mrzero().build()
        graph = mr0.compute_graph(s, data, max_states, min_signal)
        if device == "cuda":
            s, data = s.cuda(), data.cuda()
        signal = mr0.execute_graph(graph, s, data, min_signal, min_signal, print_progress=False)
        signal = signal[:, 0].detach().cpu().numpy().astype(np.complex64)
        if not np.all(np.isfinite(signal)) or not np.any(signal):
            # e.g. a corrupted sequence or a failed simulation: never cache that
            raise RuntimeError("MR-zero returned an empty or non-finite signal")
        if verbose:
            print(f"simulated {len(s)} TRs, {signal.size} samples in {time.time() - start:.1f} s")
        if cached is not None:
            cached.parent.mkdir(parents=True, exist_ok=True)
            np.save(cached, signal)
    return signal


def reconstruct(signal, protocol, n_repeats=1):
    """Echo images from simulated MESS or bSSFP data.

    Returns an array of shape ``(N, nx, ny)`` (``N`` echoes in acquisition
    order, maps indexed ``[x, y]``) in the phase convention of the paper.
    For a one-line protocol (``matrix[1] == 1``) the result has shape
    ``(N, n_repeats, nx)``: one 1D profile per TR.
    """
    n_ro, n_pe = protocol.matrix[0], protocol.matrix[1]
    lines = np.asarray(signal).reshape(-1, protocol.n_echoes * n_ro)
    echoes = split_echoes(lines, protocol.n_echoes)          # (N, lines, n_ro)
    # The ADC samples sit half a k-space step off the echo centre, which the
    # centred FFT turns into a linear phase across x; remove it.
    ramp = np.exp(1j * np.pi * (np.arange(n_ro) - n_ro // 2) / n_ro)
    if n_pe == 1:
        profiles = ifft_centered(echoes, axes=(-1,)) * ramp
        return np.conj(profiles)
    kspace = echoes.reshape(protocol.n_echoes, n_repeats, n_pe, n_ro).mean(axis=1)
    images = ifft_centered(kspace, axes=(-2, -1)) * ramp       # (N, y, x)
    return np.conj(np.swapaxes(images, -1, -2))
