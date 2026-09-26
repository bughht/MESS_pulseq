"""MESS: multi-echo steady-state free precession with Pulseq.

Reference implementation of

    H. Hong, M. N. Hoff, W. Yang, Y. Dong, Z. Zhou, P. Hu.
    Multi-Echo SSFP: A Method for Synthetic bSSFP MRI Contrast Without Banding
    Artifacts at High Field. Magn Reson Med. 2026;96(4):1696-1708.
    https://doi.org/10.1002/mrm.70472

Modules
-------
``mess.orders``      Algorithm 1: readout moments for any contiguous set of dephasing orders
``mess.sequence``    Pulseq protocols :class:`MESS2D`, :class:`MESS3D` and their bSSFP reference
``mess.recon``       echo splitting, complex sum (Eq. 6) and magnitude sum (Eq. 7)
``mess.theory``      steady-state signal model (Eqs. 1-5, Supporting Information S1-S4)
``mess.simulation``  MR-zero helpers used by the notebooks (needs ``MRzeroCore``)
``mess.plotting``    figures and animations used by the notebooks (needs ``matplotlib``)
"""

from . import recon, theory
from .orders import BSSFP_MOMENTS, dephasing_orders, readout_moments, refocusing_moment
from .recon import complex_sum, magnitude_sum, split_echoes
from .sequence import MESS2D, MESS3D, adc_kspace, default_system, flip_angle_ramp

__version__ = "1.0.0"

__all__ = [
    "BSSFP_MOMENTS", "MESS2D", "MESS3D", "adc_kspace", "complex_sum", "default_system",
    "dephasing_orders", "flip_angle_ramp", "magnitude_sum", "readout_moments", "recon",
    "refocusing_moment", "split_echoes", "theory",
]
