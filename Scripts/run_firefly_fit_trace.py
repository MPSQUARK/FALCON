#!/usr/bin/env python3
"""Run firefly on one spectrum with FIT_TRACE checkpoints (MaStar-E, Kroupa).

Usage (from firefly repo):
  set FIREFLY_FIT_TRACE=1
  set STELLARPOPMODELS_DIR=..\\FALCON\\assets\\stellar-population-models
  python ..\\FALCON\\Scripts\\run_firefly_fit_trace.py ..\\FALCON\\assets\\sample-data\\spec-0266-51630-0016.fits

Requires MaStar_SSP_v0.2.fits.gz under STELLARPOPMODELS_DIR (firefly native format).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

from fit_trace import emit, stats  # noqa: E402

os.environ.setdefault("FF_DIR", os.getcwd())
if "STELLARPOPMODELS_DIR" not in os.environ:
    falcon_assets = os.path.normpath(os.path.join(_script_dir, "..", "FALCON", "assets", "stellar-population-models"))
    if os.path.isdir(falcon_assets):
        os.environ["STELLARPOPMODELS_DIR"] = falcon_assets

sys.path.insert(0, os.environ["FF_DIR"])

import astropy.cosmology as co  # noqa: E402
import astropy.io.fits as pyfits  # noqa: E402

import firefly_models as fm  # noqa: E402
import firefly_setup as fs  # noqa: E402


def load_sdss_fits(path: str, redshift_override: float | None = None) -> fs.firefly_setup:
    with pyfits.open(path) as hdul:
        table = hdul[1].data
        wave = 10 ** np.asarray(table["loglam"], dtype=float)
        flux = np.asarray(table["flux"], dtype=float)
        ivar = np.asarray(table["ivar"], dtype=float)
        err = np.zeros_like(flux)
        good = ivar > 0
        err[good] = np.sqrt(1.0 / ivar[good])

        pri = hdul[0].header
        z = float(pri.get("z", pri.get("REDSHIFT", 0)))
        ra = float(pri.get("ra", 0))
        dec = float(pri.get("dec", 0))
        vdisp = float(pri.get("vdisp", 100))

        if len(hdul) > 2 and hdul[2].data is not None:
            row = hdul[2].data
            if "Z" in row.names and float(row["Z"][0]) > 1e-6:
                z = float(row["Z"][0])
            elif "Z_NOQSO" in row.names and float(row["Z_NOQSO"][0]) > 1e-6:
                z = float(row["Z_NOQSO"][0])
            if "PLUG_RA" in row.names:
                ra = float(row["PLUG_RA"][0])
            if "PLUG_DEC" in row.names:
                dec = float(row["PLUG_DEC"][0])
            if "VDISP" in row.names and float(row["VDISP"][0]) > 0:
                vdisp = float(row["VDISP"][0])

    if redshift_override is not None:
        z = float(redshift_override)

    r_inst = np.full_like(wave, 2000.0, dtype=float)
    emlines = [
        "H-α", "N-II", "H-β", "O-III", "H-γ", "H-δ", "H-ε",
    ]
    spec = fs.firefly_setup(path, milky_way_reddening=True, hpf_mode="on", N_angstrom_masked=20)
    spec.openSingleSpectrum(wave, flux, err, z, ra, dec, vdisp, emlines, r_inst)
    return spec


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: run_firefly_fit_trace.py <spectrum.fits> [redshift_override]", file=sys.stderr)
        return 2

    path = os.path.abspath(sys.argv[1])
    z_override = float(sys.argv[2]) if len(sys.argv) > 2 else None
    if not os.path.isfile(path):
        print(f"spectrum not found: {path}", file=sys.stderr)
        return 2

    mastar = os.path.join(os.environ["STELLARPOPMODELS_DIR"], "MaStar_SSP_v0.2.fits.gz")
    if not os.path.isfile(mastar):
        print(f"MaStar models not found: {mastar}", file=sys.stderr)
        print("Place MaStar_SSP_v0.2.fits.gz in STELLARPOPMODELS_DIR or set the env var.", file=sys.stderr)
        return 3

    t0 = time.time()
    cosmo = co.Planck15
    spec = load_sdss_fits(path, z_override)
    out_base = os.path.join(os.environ["FF_DIR"], "output", "trace-" + os.path.basename(path).replace(".fits", ""))

    emit("input", {
        "source": os.path.basename(path),
        "redshift": float(spec.redshift),
        "n_wave": int(len(spec.restframe_wavelength)),
        **stats("data", spec.flux),
        **stats("error", spec.error),
    })

    model = fm.StellarPopulationModel(
        spec, out_base, cosmo, t0,
        models="MaStar",
        model_libs=["E-MaStar"],
        imfs=["kr"],
        hpf_mode="on",
        age_limits=[0, 14.7],
        downgrade_models=True,
        dust_law="calzetti",
        max_ebv=1.0,
        num_dust_vals=200,
        dust_smoothing_length=100,
        max_iterations=10,
        fit_per_iteration_cap=1000,
        pdf_sampling=200,
        Z_limits=[-3, 3],
        wave_limits=[0, 99999990],
        write_results=False,
        flux_units=1e-17,
    )
    model.fit_models_to_data(t0)
    emit("done", {"elapsed_s": time.time() - t0})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
