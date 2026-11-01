#!/usr/bin/env python3
"""Print SDSS spec FITS metadata (redshift, vdisp) for verification."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from astropy.io import fits


def inspect(path: Path) -> None:
    with fits.open(path) as hdul:
        print(f"\n=== {path.name} ===")
        hdul.info()
        pri = hdul[0].header
        for key in ("z", "Z", "REDSHIFT", "RA", "DEC", "VDISP"):
            if key in pri:
                print(f"  HDU0 {key} = {pri[key]}")

        if len(hdul) > 2 and hdul[2].data is not None:
            row = hdul[2].data
            names = row.names
            print(f"  HDU2 columns: {len(names)}")
            for col in ("Z", "Z_NOQSO", "Z_ERR", "ZWARNING", "VDISP", "VDISPN", "PLUG_RA", "PLUG_DEC"):
                if col in names:
                    print(f"  SPALL {col} = {row[col][0]}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fits", nargs="+", type=Path)
    args = p.parse_args()
    for path in args.fits:
        if not path.is_file():
            print(f"Missing: {path}", file=sys.stderr)
            return 2
        inspect(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
