#!/usr/bin/env python3
"""Plot FALCON / firefly spectral fit results."""
import argparse
import matplotlib.pyplot as plt
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("spectrum_txt")
    p.add_argument("--out", default="plot.png")
    args = p.parse_args()

    # Parse minimal sections from FALCON text output
    text = open(args.spectrum_txt).read()
    if "Flux of Data" not in text:
        raise SystemExit("Unrecognised FALCON output format")
    data_flux = np.fromstring(text.split("Flux of Data")[1].split("Error")[0], sep=" ")
    wave = np.arange(len(data_flux))

    plt.figure(figsize=(10, 4))
    plt.plot(wave, data_flux, label="data")
    plt.xlabel("pixel")
    plt.ylabel("flux")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=120)
    print(f"Wrote {args.out}")

if __name__ == "__main__":
    main()
