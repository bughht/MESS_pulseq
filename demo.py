"""
Author: Haotian Hong
Date: 2025/4/22
Description: This script demonstrates the usage of the MESS_3D class to generate and save 3D MESS and bSSFP sequences with varying flip angles and phase cycling configurations.
"""

from MESS import MESS_3D
import pypulseq as pp
import os

save_path = "seq"
sys = pp.Opts(max_grad=80, grad_unit='mT/m', max_slew=80, slew_unit='T/m/s',
              rf_ringdown_time=10e-6, rf_dead_time=100e-6, adc_dead_time=20e-6, grad_raster_time=10e-6)

for FlipAngle in [20, 50, 80]:
    mess = MESS_3D(FA=FlipAngle, TR=10e-3, dwell=5e-6, rf_duration=1.5e-3, FLAG_min_TR=True,
                   num_PE=150, num_RO=150, num_SPE=75, fov_PE=200e-3, fov_RO=200e-3, fov_SPE=100e-3, system=sys)
    for PhaseCycle in [0, 180]:
        seq_bssfp = mess.make_sequence(
            0, 0, None, balance=True, delay_pre=2.26e-3, delay_post=3.01e-3, phase_cycle=PhaseCycle)
        seq_bssfp.write(os.path.join(save_path, "{}deg".format(
            FlipAngle), "seq_bssfp_{}pc".format(PhaseCycle)))

        seq_mess = mess.make_sequence(+2, -3, phase_cycle=PhaseCycle)
        seq_mess.write(os.path.join(save_path, "{}deg".format(
            FlipAngle), "seq_mess_+2_-3_{}pc".format(PhaseCycle)))
