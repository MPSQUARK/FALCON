#!/usr/bin/env python3
"""Compare FIT_TRACE checkpoint files from FALCON and firefly.

Usage:
  python Scripts/compare_fit_traces.py falcon.trace firefly.trace
  python Scripts/compare_fit_traces.py --rtol 0.05 falcon.trace firefly.trace

Parses lines: FIT_TRACE\\t{stage}\\t{json}
Reports the first stage where numeric fields diverge beyond tolerance.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

LINE_RE = re.compile(r"^FIT_TRACE\t([^\t]+)\t(.+)$")

# Keys that are counts / indices — compared exactly or with small abs tol
INTEGER_KEYS = {"n", "finite", "positive", "n_models", "n_solutions", "dof"}


def load_trace(path: Path) -> dict[str, dict]:
    stages: dict[str, dict] = {}
    for i, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        m = LINE_RE.match(raw.strip())
        if not m:
            continue
        stage, payload = m.group(1), m.group(2)
        try:
            stages[stage] = json.loads(payload)
        except json.JSONDecodeError as ex:
            print(f"WARN {path}:{i} bad JSON for stage {stage}: {ex}", file=sys.stderr)
    return stages


def numeric_diff(key: str, a, b, rtol: float, atol: float) -> str | None:
    if a is None and b is None:
        return None
    if a is None or b is None:
        return f"{key}: {a!r} vs {b!r}"
    if isinstance(a, bool) or isinstance(b, bool):
        return None if a == b else f"{key}: {a!r} vs {b!r}"
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        return None if a == b else f"{key}: {a!r} vs {b!r}"

    fa, fb = float(a), float(b)
    if not math.isfinite(fa) and not math.isfinite(fb):
        return None
    if not math.isfinite(fa) or not math.isfinite(fb):
        return f"{key}: {fa!r} vs {fb!r}"

    suffix = key.split("_")[-1]
    if suffix in INTEGER_KEYS or key.endswith("_n") or key.endswith("_finite"):
        if abs(fa - fb) <= max(atol, 1):
            return None
        return f"{key}: {fa} vs {fb} (Δ={fa-fb:g})"

    if abs(fa - fb) <= atol + rtol * max(abs(fa), abs(fb), 1e-30):
        return None
    rel = abs(fa - fb) / max(abs(fb), 1e-30)
    return f"{key}: {fa:.6g} vs {fb:.6g} (rel Δ={rel:.3g})"


def compare_stage(stage: str, left: dict, right: dict, rtol: float, atol: float) -> list[str]:
    diffs: list[str] = []
    keys = sorted(set(left) | set(right))
    for key in keys:
        if key not in left:
            diffs.append(f"missing in FALCON: {key}={right[key]!r}")
            continue
        if key not in right:
            diffs.append(f"missing in firefly: {key}={left[key]!r}")
            continue
        msg = numeric_diff(key, left[key], right[key], rtol, atol)
        if msg:
            diffs.append(msg)
    return diffs


def main() -> int:
    p = argparse.ArgumentParser(description="Compare FIT_TRACE files from FALCON and firefly")
    p.add_argument("falcon_trace", type=Path)
    p.add_argument("firefly_trace", type=Path)
    p.add_argument("--rtol", type=float, default=0.02, help="relative tolerance for floats")
    p.add_argument("--atol", type=float, default=1e-6, help="absolute tolerance")
    args = p.parse_args()

    falcon = load_trace(args.falcon_trace)
    firefly = load_trace(args.firefly_trace)

    if not falcon:
        print(f"No FIT_TRACE lines in {args.falcon_trace}", file=sys.stderr)
        return 2
    if not firefly:
        print(f"No FIT_TRACE lines in {args.firefly_trace}", file=sys.stderr)
        return 2

    all_stages = sorted(set(falcon) | set(firefly))
    first_divergence: str | None = None

    print(f"{'STAGE':<24} {'STATUS':<8} detail")
    print("-" * 72)
    for stage in all_stages:
        if stage not in falcon:
            print(f"{stage:<24} {'MISS-F':<8} only in firefly")
            if first_divergence is None:
                first_divergence = stage
            continue
        if stage not in firefly:
            print(f"{stage:<24} {'MISS-P':<8} only in FALCON")
            if first_divergence is None:
                first_divergence = stage
            continue

        diffs = compare_stage(stage, falcon[stage], firefly[stage], args.rtol, args.atol)
        if diffs:
            print(f"{stage:<24} {'DIFF':<8} {diffs[0]}")
            if len(diffs) > 1:
                for d in diffs[1:4]:
                    print(f"{'':24} {'':8} {d}")
                if len(diffs) > 4:
                    print(f"{'':24} {'':8} ... +{len(diffs)-4} more")
            if first_divergence is None:
                first_divergence = stage
        else:
            print(f"{stage:<24} {'OK':<8}")

    print()
    if first_divergence:
        print(f"First divergence at stage: {first_divergence}")
        return 1
    print("All shared stages match within tolerance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
