"""
Author: Haotian Hong
Date: 2025/4/22
Description: This module implements the MESS_3D class for generating 3D MESS sequences. The class supports flexible configuration of imaging parameters and sequence generation.
"""

import numpy as np
import pypulseq as pp
from pypulseq.make_sigpy_pulse import sigpy_n_seq
from . import arbitarySSFP


class MESS_3D:
    """
    The MESS_3D class is designed to create 3D MRI sequences with customizable parameters such as flip angle, repetition time, field of view, and resolution.
    """

    def __init__(self,
                 FA=35,
                 TR=10e-3,
                 num_RO=120,
                 num_PE=120,
                 num_SPE=30,
                 fov_RO=200e-3,
                 fov_PE=200e-3,
                 fov_SPE=50e-3,
                 dwell=1e-05,
                 rf_duration=1e-3,
                 FLAG_min_TR=False,
                 system=None,
                 ):
        """
        Initialize the MESS_3D object with imaging parameters.

        Args:
            FA (float): Flip angle in degrees.
            TR (float): Repetition time in seconds.
            num_RO (int): Number of readout points.
            num_PE (int): Number of phase encoding steps.
            num_SPE (int): Number of slice phase encoding steps.
            fov_RO (float): Field of view in the readout direction (meters).
            fov_PE (float): Field of view in the phase encoding direction (meters).
            fov_SPE (float): Field of view in the slice phase encoding direction (meters).
            dwell (float): Dwell time for ADC sampling (seconds).
            rf_duration (float): Duration of the RF pulse (seconds).
            FLAG_min_TR (bool): Flag to automatically adjust TR to the minimum possible value.
            system (pp.Opts): System hardware constraints. Defaults to a predefined configuration.
        """
        self.FA = FA
        self.TR = TR
        self.num_RO = num_RO
        self.num_PE = num_PE
        self.num_SPE = num_SPE
        self.fov_RO = fov_RO
        self.fov_PE = fov_PE
        self.fov_SPE = fov_SPE
        self.dwell = dwell
        self.rf_duration = rf_duration
        self.FLAG_min_TR = FLAG_min_TR
        if system is None:
            self.system = pp.Opts(max_grad=40, grad_unit='mT/m', max_slew=150, slew_unit='T/m/s',
                                  rf_ringdown_time=10e-6, rf_dead_time=400e-6, adc_dead_time=70e-6, grad_raster_time=10e-6)
        else:
            self.system = system
        self.prepare_sequence()

    def prepare_sequence(self):
        """
        Prepare the basic sequence components, including RF pulses, gradients, and phase encoding steps.

        This method calculates the gradient areas, RF pulse configurations, and phase encoding tables 
        required for the 3D MRI sequence.
        """
        gx_single = pp.make_trapezoid(channel='x', flat_area=self.num_RO /
                                      self.fov_RO, flat_time=self.dwell*self.num_RO, system=self.system)
        self.ramp_area = (gx_single.area-gx_single.flat_area)/2
        self.half_grad_area = gx_single.flat_area/2

        self.rf_phase = 0

        pulse_cfg = pp.SigpyPulseOpts(
            pulse_type='slr',
            ptype='st',
            ftype='ls',
            d1=0.01,
            d2=0.01,
            cancel_alpha_phs=True,
            n_bands=4,
            band_sep=20,
            phs_0_pt='None',
        )

        self.rf, self.gz, self.gzr, _ = sigpy_n_seq(
            flip_angle=np.deg2rad(self.FA),
            duration=self.rf_duration,
            return_gz=True,
            slice_thickness=self.fov_SPE*0.5,
            system=self.system,
            # time_bw_product=4,
            time_bw_product=2,  # High Flip Angle
            pulse_cfg=pulse_cfg,
            plot=False,
        )

        # sinc_rf, self.gz, self.gzr = pp.make_sinc_pulse(
        #     flip_angle=self.FA * np.pi / 180, duration=self.rf_duration,
        #     # slice_thickness=self.fov_SPE*0.6, apodization=0.5, time_bw_product=4,
        #     slice_thickness=self.fov_SPE*0.25, apodization=0.2, time_bw_product=2.7,
        #     system=self.system, use="excitation", return_gz=True)
        # self.rf, self.gz, self.gzr = pp.make_gauss_pulse(
        #     flip_angle=self.FA * np.pi / 180, duration=self.rf_duration,
        #     # slice_thickness=self.fov_SPE*0.6, apodization=0.5, time_bw_product=4,
        #     slice_thickness=self.fov_SPE*0.25, apodization=0.2, time_bw_product=2.7,
        #     system=self.system, use="excitation", return_gz=True)

        self.rf.delay = self.gz.rise_time
        self.gz.delay = 0
        self.gzr_area = self.gzr.area
        # print("gzr_area:{}".format(self.gzr_area),
        #       "gz_area:{}".format(self.gz.area))

        self.SPE = np.linspace(-self.num_SPE / 2,
                               self.num_SPE / 2, self.num_SPE, endpoint=True) / self.fov_SPE
        self.PE = np.linspace(-self.num_PE / 2,
                              self.num_PE / 2, self.num_PE, endpoint=True) / self.fov_PE

    def make_sequence(self, start_echo=0, end_echo=0, ascend=None, balance=False, phase_cycle=180, delay_pre=None, delay_post=None):
        """
        Generate the 3D MRI sequence based on the specified parameters.

        Args:
            start_echo (int): Starting echo index.
            end_echo (int): Ending echo index.
            ascend (bool, optional): Direction of the echo sequence.
            balance (bool, optional): Whether to use a balanced SSFP sequence.
            phase_cycle (int): Phase cycling increment in degrees.
            delay_pre (float, optional): Pre-readout delay (seconds).
            delay_post (float, optional): Post-readout delay (seconds).

        Returns:
            pp.Sequence: The generated MRI sequence.
        """
        seq = pp.Sequence(self.system)
        a, b, c = arbitarySSFP(start_echo, end_echo, ascend, balance)
        # print("a:{} b:{} c:{}".format(a, b, c))
        num_echo = abs(start_echo - end_echo) + 1

        gx_pre = pp.make_trapezoid(
            channel='x', area=-self.ramp_area+self.half_grad_area*a, system=self.system)
        gx_rep = pp.make_trapezoid(
            channel='x', area=-self.ramp_area+self.half_grad_area*c, system=self.system)
        gx_ro = pp.make_trapezoid(channel='x', flat_area=self.half_grad_area*b,
                                  flat_time=self.dwell*self.num_RO*num_echo, system=self.system)
        adc = pp.make_adc(num_samples=self.num_RO*num_echo, duration=self.dwell *
                          self.num_RO*num_echo, system=self.system)
        adc.delay = gx_ro.rise_time
        # print("Delta TE: ", pp.calc_duration(adc))

        gpe_max = pp.make_trapezoid(channel='y', area=np.abs(
            self.PE).max(), system=self.system)
        gspe_max = pp.make_trapezoid(channel='z', area=np.abs(
            self.SPE+self.gzr_area).max(), system=self.system)
        gspe_rep_max = pp.make_trapezoid(channel='z', area=np.abs(-self.SPE+self.gzr_area).max(),
                                         system=self.system)
        min_pre_duration = max(pp.calc_duration(
            gpe_max), pp.calc_duration(gx_pre), pp.calc_duration(gspe_max))
        min_rep_duration = max(pp.calc_duration(
            gpe_max), pp.calc_duration(gx_rep), pp.calc_duration(gspe_rep_max))

        min_TR = min_pre_duration+min_rep_duration + \
            pp.calc_duration(gx_ro)+pp.calc_duration(self.gz)
        if not self.FLAG_min_TR:
            assert min_TR < self.TR, "TR {} is too short, minimum TR is {}".format(
                self.TR, min_TR)
        else:
            self.TR = min_TR
        # print("min_TR:{}".format(min_TR))
        delay_time = (self.TR-min_TR)
        for spe_idx, spe_area in enumerate(self.SPE):
            gspe_pre = pp.make_trapezoid(
                channel='z', area=spe_area+self.gzr_area, duration=min_pre_duration, system=self.system)
            gspe_rep = pp.make_trapezoid(
                channel='z', area=-spe_area+self.gzr_area, duration=min_rep_duration, system=self.system)
            for pe_idx, pe_area in enumerate(self.PE):
                # for pe_idx, pe_area in enumerate(self.PE[:5]):
                shot_idx = spe_idx*len(self.PE)+pe_idx
                self.rf.phase_offset = shot_idx*np.deg2rad(phase_cycle)
                adc.phase_offset = shot_idx*np.deg2rad(phase_cycle)
                gpe_pre = pp.make_trapezoid(
                    channel='y', area=pe_area, duration=min_pre_duration, system=self.system)
                gpe_rep = pp.make_trapezoid(
                    channel='y', area=-pe_area, duration=min_rep_duration, system=self.system)

                seq.add_block(self.rf, self.gz)
                seq.add_block(
                    *pp.align(right=gx_pre, left=[gpe_pre, gspe_pre]))
                if delay_pre is not None:
                    seq.add_block(pp.make_delay(delay_pre))
                seq.add_block(gx_ro, adc)
                if delay_post is not None:
                    seq.add_block(pp.make_delay(delay_post))
                seq.add_block(
                    *pp.align(left=[gx_rep, pp.make_delay(delay_time+min_rep_duration)], right=[gpe_rep, gspe_rep]))

        return seq
