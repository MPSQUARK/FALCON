#!/usr/bin/env python3
"""Convert FALCON text output to firefly-compatible FITS."""
import argparse
import numpy as np
from astropy.io import fits

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_txt", help="FALCON output .txt from SpectralFitPipeline")
    p.add_argument("output_fits")
    args = p.parse_args()

    # Minimal placeholder FITS — extend with HIERARCH headers from parsed output
    wave = np.linspace(3500, 9000, 100)
    flux = np.ones_like(wave)
    err = np.ones_like(wave) * 0.01
    model = flux.copy()

    tbl = fits.BinTableHDU.from_columns([
        fits.Column(name="wavelength", format="D", array=wave),
        fits.Column(name="flux", format="D", array=flux),
        fits.Column(name="error", format="D", array=err),
        fits.Column(name="firefly_model", format="D", array=model),
    ])
    tbl.header["HIERARCH E(B-V)"] = (0.0, "dust attenuation")
    fits.HDUList([fits.PrimaryHDU(), tbl]).writeto(args.output_fits, overwrite=True)

if __name__ == "__main__":
    main()
