#!/usr/bin/env python3
"""Print firefly spectral_resolution.dv for the MaStar model wavelength grid.

Usage (from firefly repo with STELLARPOPMODELS_DIR set):
  python ../FALCON/Scripts/Analysis/reference_dv.py

Or standalone with numpy only (matches firefly_instrument.spectrum_velocity_scale):
  python Scripts/Analysis/reference_dv.py --numpy-only --wave-file path/to/wavelength.npy
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_FALCON_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, ".."))


def dv_from_wave_numpy(wave: np.ndarray) -> float:
    """firefly: c_kms * mean(diff(ln(wave))) with natural log."""
    c_kms = 299792.458
    return float(c_kms * np.mean(np.diff(np.log(wave))))


def load_mastar_wavelength() -> np.ndarray:
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise SystemExit("astropy required: pip install astropy") from exc

    models_dir = os.environ.get(
        "FALCON_MODELS",
        os.path.join(_FALCON_ROOT, "Assets", "StellarPopulationModels"),
    )
    fits_path = os.path.join(models_dir, "MaStar_SSP_v0.2.fits.gz")
    if not os.path.isfile(fits_path):
        raise SystemExit(f"MaStar model not found: {fits_path}")

    with fits.open(fits_path, memmap=False) as hdul:
        return np.asarray(hdul[2].data[0], dtype=np.float64)


def main() -> int:
    p = argparse.ArgumentParser(description="Reference dv for MaStar wavelength grid")
    p.add_argument("--numpy-only", action="store_true", help="Skip firefly import")
    p.add_argument("--wave-file", type=str, default="", help="Optional .npy wavelength array")
    p.add_argument("--json", action="store_true", help="Emit JSON for golden tests")
    args = p.parse_args()

    if args.wave_file:
        wave = np.load(args.wave_file)
    else:
        wave = load_mastar_wavelength()

    dv_numpy = dv_from_wave_numpy(wave)
    result = {"dv_numpy": dv_numpy, "n_wave": int(wave.size)}

    if not args.numpy_only:
        ff_dir = os.environ.get("FF_DIR", os.path.join(_FALCON_ROOT, "..", "firefly"))
        if os.path.isdir(ff_dir) and ff_dir not in sys.path:
            sys.path.insert(0, ff_dir)
        try:
            from firefly_instrument import spectral_resolution

            res = spectral_resolution(wave.astype(np.float64), np.ones_like(wave), log10=True)
            result["dv_firefly"] = float(res.dv)
        except ImportError:
            result["dv_firefly"] = None

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"n_wave={result['n_wave']}")
        print(f"dv_numpy={dv_numpy:.17g}")
        if result.get("dv_firefly") is not None:
            print(f"dv_firefly={result['dv_firefly']:.17g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

