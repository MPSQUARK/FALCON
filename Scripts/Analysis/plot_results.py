#!/usr/bin/env python3
"""Plot FALCON .results.json output: restframe data vs SSP model blend."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_COLOR = "#5B9BD5"
ERROR_FILL = "#B8D4EE"
MODEL_COLOR = "#E0634A"


def load_results(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def align_flux(wave: np.ndarray, model: np.ndarray) -> np.ndarray:
    n = min(len(wave), len(model))
    out = np.full(len(wave), np.nan, dtype=float)
    out[:n] = model[:n]
    return out


def finite_data_extent(data: np.ndarray) -> slice:
    finite = np.isfinite(data)
    if not np.any(finite):
        return slice(0, len(data))
    idx = np.where(finite)[0]
    return slice(int(idx[0]), int(idx[-1]) + 1)


def select_range(wave, data, model, doc, snip: bool):
    n = len(wave)
    data_sl = finite_data_extent(data)
    if not snip:
        return data_sl

    w_min = doc.get("model_wave_min") or 0
    w_max = doc.get("model_wave_max") or 0
    start = data_sl.start
    end = data_sl.stop
    if w_min > 0 and w_max > w_min:
        model_start = int(np.searchsorted(wave, w_min, side="left"))
        model_end = int(np.searchsorted(wave, w_max, side="right"))
        if model_end > model_start:
            start = max(start, model_start)
            end = min(end, model_end)

    if model.size == 0:
        return slice(start, end)

    max_abs = float(np.max(np.abs(model[:n]))) if n else 0.0
    eps = max_abs * 1e-12
    active = np.abs(model[start:end]) > eps
    if not np.any(active):
        return slice(start, end)
    rel = np.where(active)[0]
    return slice(start + int(rel[0]), start + int(rel[-1]) + 1)


def compute_display_limits(wave, data, model):
    x_min = float(wave[0])
    x_max = float(wave[-1])
    if len(wave) > 1:
        x_pad = (x_max - x_min) * 0.01
        x_min -= x_pad
        x_max += x_pad

    model_finite = model[np.isfinite(model)] if model.size else np.array([])
    if model_finite.size == 0:
        data_only = data[np.isfinite(data)]
        if data_only.size == 0:
            return x_min, x_max, 0.0, 1.0
        data_only = np.sort(data_only)
        lo = float(data_only[int(len(data_only) * 0.02)])
        hi = float(data_only[int(len(data_only) * 0.98)])
        pad = max((hi - lo) * 0.08, 0.1)
        return x_min, x_max, lo - pad, hi + pad

    model_vals = np.sort(model_finite)
    model_lo = float(model_vals[0])
    model_hi = float(model_vals[-1])
    model_med = float(model_vals[len(model_vals) // 2])
    model_span = max(model_hi - model_lo, max(model_med * 0.05, 0.5))
    band = model_span * 3 + max(model_med * 0.5, 1.0)

    near = []
    n = min(len(wave), len(data), len(model) if model.size else len(data))
    for i in range(n):
        if np.isfinite(data[i]):
            m = float(model[i]) if model.size > i and np.isfinite(model[i]) else model_med
            if abs(float(data[i]) - m) <= band:
                near.append(float(data[i]))
        if model.size > i and np.isfinite(model[i]):
            near.append(float(model[i]))

    if not near:
        near = model_vals.tolist()

    near = np.sort(np.asarray(near, dtype=float))
    count = len(near)
    p_lo = float(near[int(count * 0.02)])
    p_hi = float(near[int(count * 0.98)])
    p_hi = min(p_hi, model_hi + model_span * 0.35)
    p_lo = max(p_lo, model_lo - model_span * 0.15)
    if p_hi <= p_lo:
        p_lo, p_hi = model_lo, model_hi

    y_span = p_hi - p_lo
    y_pad = y_span * 0.06 if y_span > 0 else max(abs(p_hi), 1.0) * 0.06
    y_pad = max(y_pad, model_span * 0.08)
    return x_min, x_max, p_lo - y_pad, p_hi + y_pad


def plot_results(doc: dict, out_path: Path, snip: bool = True, show: bool = False) -> None:
    wave = np.asarray(doc["wavelength_angstrom"], dtype=float)
    data = np.asarray(doc["data_flux"], dtype=float)
    err = np.abs(np.asarray(doc["error_flux"], dtype=float))
    top = [s for s in (doc.get("top_solutions") or []) if len(s.get("model_flux") or []) > 0]
    best = top[0] if top else {}
    model = align_flux(wave, np.asarray(best.get("model_flux") or [], dtype=float))

    sl = select_range(wave, data, model, doc, snip)
    wave, data, err, model = wave[sl], data[sl], err[sl], model[sl]

    fig, ax = plt.subplots(figsize=(20, 10), facecolor="white")
    ax.set_facecolor("white")

    ax.fill_between(wave, data - err, data + err, color=ERROR_FILL, alpha=0.55, label="±1σ error", linewidth=0)
    ax.plot(wave, data, color=DATA_COLOR, lw=0.9, label="Data")

    if model.size > 0 and np.any(np.isfinite(model)):
        ax.plot(wave, model, color=MODEL_COLOR, lw=0.9, label=f"FALCON (#{best.get('rank', 1)})")

    x_min, x_max, y_min, y_max = compute_display_limits(wave, data, model)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    chi = best.get("chi_squared_per_dof") or 0
    if chi <= 0:
        dof = doc.get("degrees_of_freedom") or 1
        min_chi = doc.get("min_chi_squared") or 0
        chi = min_chi / dof if dof else 0

    src = Path(doc.get("source_file", "")).stem
    range_note = "model overlap" if snip else "full spectrum"
    ax.set_xlabel("Rest-frame wavelength (Å)")
    ax.set_ylabel("Flux (erg/s/Å/cm²)")
    ax.set_title(f"χ²/ν: {chi:.2f}  |  {src}  |  {range_note}")
    ax.grid(True, color="#E8EDF4", linewidth=0.8)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.92)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor="white")
    if show:
        plt.show()
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description="Plot FALCON results JSON")
    p.add_argument("results_json", type=Path, help="Path to *.results.json")
    p.add_argument("--out", type=Path, default=None, help="Output PNG (default: alongside JSON)")
    p.add_argument("--snip", action="store_true", help="Crop to model wavelength overlap (default)")
    p.add_argument("--full", action="store_true", help="Show full data wavelength range")
    p.add_argument("--show", action="store_true")
    args = p.parse_args()
    snip = not args.full
    doc = load_results(args.results_json)
    out = args.out or args.results_json.with_suffix(".png")
    plot_results(doc, out, snip=snip, show=args.show)


if __name__ == "__main__":
    main()
