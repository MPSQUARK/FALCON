#!/usr/bin/env python3
"""Convert FITS/HDF5 spectra to FALCON-readable ASCII triplets (wavelength, flux, error)."""
import argparse
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output", help="Output .dat path")
    args = p.parse_args()

    if args.input.endswith(".fits") or args.input.endswith(".fit"):
        from astropy.io import fits
        with fits.open(args.input) as hdul:
            table = hdul[1].data if len(hdul) > 1 else hdul[0].data
            flux = np.asarray(table["flux"] if "flux" in table.names else table[0], dtype=float)
            wave = np.asarray(table["loglam"] if "loglam" in table.names else table[1], dtype=float)
            if wave.mean() < 100:
                wave = 10 ** wave
            err = np.sqrt(np.maximum(flux, 1e-30))
    else:
        data = np.loadtxt(args.input)
        wave, flux, err = data[:, 0], data[:, 1], data[:, 2]

    out = np.column_stack([wave, flux, err])
    np.savetxt(args.output, out, fmt="%.8e")

if __name__ == "__main__":
    main()
