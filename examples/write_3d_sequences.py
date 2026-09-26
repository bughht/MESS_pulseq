"""Write the 3D MESS sequences and their matched bSSFP references for the scanner.

    python examples/write_3d_sequences.py              # every protocol below
    python examples/write_3d_sequences.py paper        # only the protocol of the paper
    python examples/write_3d_sequences.py --list

Each protocol is written for every flip angle and phase-cycle increment as::

    seq/3d/<protocol>/<FA>deg/seq_mess_<k1>_<kN>_<pc>pc.seq
    seq/3d/<protocol>/<FA>deg/seq_bssfp_<pc>pc.seq

The bSSFP reference shares TR, TE (= TE_0 of MESS), RF pulse and encoding with
the MESS sequence, so the two can be compared directly.  The gradient limits
are those of the 5T scanner of the paper (``mess.default_system()``); pass your
own ``pypulseq.Opts`` to the protocol classes for another scanner.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from a clone without installing
import mess  # noqa: E402

PROTOCOLS = {
    # Hong et al., MRM 2026: six orders, 1.3 mm isotropic over 200 x 200 x 100 mm (Methods).
    "paper": dict(orders=(2, -3), matrix=(150, 150, 75), fov=(200e-3, 200e-3, 100e-3),
                  flip_angles=(20, 50, 80), phase_cycles=(0, 180), trigger=False),
    # small field of view, with an external trigger on every readout (June 2026 scans).
    "small_fov": dict(orders=(2, -3), matrix=(80, 80, 24), fov=(100e-3, 100e-3, 30e-3),
                      flip_angles=(10, 30, 50), phase_cycles=(0,), trigger=True),
}


def write(name, out):
    spec = PROTOCOLS[name]
    k1, kN = spec["orders"]
    for flip_angle in spec["flip_angles"]:
        protocol = mess.MESS3D(orders=(k1, kN), flip_angle=flip_angle, fov=spec["fov"],
                               matrix=spec["matrix"], dwell=5e-6, rf_duration=1.5e-3)
        print(protocol, "\n")
        folder = Path(out) / name / f"{flip_angle}deg"
        folder.mkdir(parents=True, exist_ok=True)
        for pc in spec["phase_cycles"]:
            for p, stem in ((protocol, f"seq_mess_{k1:+d}_{kN:+d}_{pc}pc"),
                            (protocol.bssfp(), f"seq_bssfp_{pc}pc")):
                seq = p.make_sequence(phase_cycle=np.deg2rad(pc), trigger=spec["trigger"])
                seq.write(str(folder / f"{stem}.seq"))
                print("wrote", folder / f"{stem}.seq")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("protocols", nargs="*",
                        help=f"protocols to write, from {', '.join(PROTOCOLS)} (default: all)")
    parser.add_argument("--out", default="seq/3d", help="output folder (default: seq/3d)")
    parser.add_argument("--list", action="store_true", help="print the protocols and exit")
    args = parser.parse_args()
    if args.list:
        for name, spec in PROTOCOLS.items():
            print(f"{name:10s} {spec}")
        return
    unknown = sorted(set(args.protocols) - set(PROTOCOLS))
    if unknown:
        parser.error(f"unknown protocol(s): {', '.join(unknown)}")
    for name in args.protocols or PROTOCOLS:
        write(name, args.out)


if __name__ == "__main__":
    main()
