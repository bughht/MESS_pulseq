"""Write the 2D sequences that the notebooks simulate, ready for the scanner.

    python examples/write_2d_sequences.py
    python examples/write_2d_sequences.py --flip-angle 35 --dummies 500

With the defaults the files are exactly the ones simulated in the notebooks::

    seq/2d/<FA>deg/seq_mess2d_+2_-3_0pc.seq     six-echo MESS                     (notebooks 3-5)
    seq/2d/<FA>deg/seq_bssfp2d_<pc>pc.seq       bSSFP, pc = 0, 15, ..., 345 deg   (notebooks 4-5)
    seq/2d/<FA>deg/seq_mess2d_<k1>_<kN>_0pc.seq FISP [0], PSIF [-1], DESS [0, -1] and
                                                TESS [1, 0, -1] with the TR and TE_0
                                                of the MESS protocol              (notebook 3)

A single slice has only 150 encoding steps, so the dummy TRs that bring the
magnetisation into the steady state dominate the scan time; notebook 2 shows
how many are needed for which tissue.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from a clone without installing
import mess  # noqa: E402

FAMILY = {"FISP": (0, 0), "PSIF": (-1, -1), "DESS": (0, -1), "TESS": (1, -1)}
BSSFP_PHASE_CYCLES = np.arange(0, 360, 15)            # degrees


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--flip-angle", type=float, default=50.0)
    parser.add_argument("--dummies", type=int, default=2000, help="dummy TRs before the first line")
    parser.add_argument("--ramp", type=int, default=0, help="linear flip-angle ramp on the first dummies")
    parser.add_argument("--out", default="seq/2d", help="output folder (default: seq/2d)")
    args = parser.parse_args()

    protocol = mess.MESS2D(orders=(2, -3), flip_angle=args.flip_angle)   # 150 x 150, 200 mm, 5 mm
    print(protocol)
    print(f"  with {args.dummies} dummies: {protocol.scan_time(args.dummies):.1f} s per image\n")
    folder = Path(args.out) / f"{args.flip_angle:g}deg"
    folder.mkdir(parents=True, exist_ok=True)

    def write(p, phase_cycle_deg, stem):
        seq = p.make_sequence(phase_cycle=np.deg2rad(phase_cycle_deg), n_dummy=args.dummies,
                              ramp=args.ramp)
        seq.write(str(folder / f"{stem}.seq"))
        print("wrote", folder / f"{stem}.seq")

    write(protocol, 0, "seq_mess2d_+2_-3_0pc")
    for pc in BSSFP_PHASE_CYCLES:
        write(protocol.bssfp(), pc, f"seq_bssfp2d_{pc}pc")
    for name, (k1, kN) in FAMILY.items():
        member = mess.MESS2D(orders=(k1, kN), flip_angle=args.flip_angle, tr=protocol.tr,
                             te=protocol.te0)
        write(member, 0, f"seq_mess2d_{k1:+d}_{kN:+d}_0pc")


if __name__ == "__main__":
    main()
