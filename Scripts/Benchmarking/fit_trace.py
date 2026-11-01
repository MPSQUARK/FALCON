#!/usr/bin/env python3
"""Shared FIT_TRACE emitter for firefly parity debugging.

Set FIREFLY_FIT_TRACE=1 (or FALCON_FIT_TRACE=1) to enable.
Output format: FIT_TRACE\\t{stage}\\t{json}

Used by Scripts/Benchmarking/compare_fit_traces.py to diff FALCON vs firefly checkpoints.
"""
from __future__ import annotations

import json
import os
import sys
from typing import Any, Mapping

import numpy as np


def enabled() -> bool:
    v = os.environ.get("FIREFLY_FIT_TRACE") or os.environ.get("FALCON_FIT_TRACE") or ""
    return v.strip().lower() in {"1", "true", "yes", "on"}


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        v = float(obj)
        if not np.isfinite(v):
            return None
        return v
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(type(obj))


def emit(stage: str, data: Mapping[str, Any]) -> None:
    if not enabled():
        return
    line = "FIT_TRACE\t" + stage + "\t" + json.dumps(dict(data), default=_json_default, separators=(",", ":"))
    print(line, flush=True)


def stats(label: str, values: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(values, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    out: dict[str, Any] = {
        f"{label}_n": int(arr.size),
        f"{label}_finite": int(finite.size),
        f"{label}_positive": int(np.sum(finite > 0)),
    }
    if finite.size:
        out[f"{label}_min"] = float(np.min(finite))
        out[f"{label}_max"] = float(np.max(finite))
        out[f"{label}_median"] = float(np.median(finite))
        out[f"{label}_mean"] = float(np.mean(finite))
    else:
        out[f"{label}_min"] = 0.0
        out[f"{label}_max"] = 0.0
        out[f"{label}_median"] = 0.0
        out[f"{label}_mean"] = 0.0
    return out


def model_grid_stats(models: np.ndarray, prefix: str = "model0") -> dict[str, Any]:
    models = np.asarray(models)
    if models.ndim == 1:
        models = models.reshape(1, -1)
    out = stats(prefix, models[0])
    out[f"{prefix}_n_models"] = int(models.shape[0])
    return out

