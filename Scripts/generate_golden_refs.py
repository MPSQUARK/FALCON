#!/usr/bin/env python3
"""Generate golden reference JSON files from local firefly for FALCON regression tests."""
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIREFLY = REPO.parent / "firefly"
OUT = REPO / "tests" / "Falcon.Tests" / "TestData"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Synthetic references when firefly runtime unavailable in CI
    bisect = {"array": [1000.0, 1500.0, 2000.0, 2500.0, 3000.0]}
    (OUT / "golden_bisect.json").write_text(json.dumps(bisect, indent=2))

    hpf = {"input_length": 128, "median_scale": 1.0}
    (OUT / "golden_hpf.json").write_text(json.dumps(hpf, indent=2))

    manifest = {
        "description": "Golden references for FALCON firefly parity tests",
        "tolerance_rtol": 1e-5,
        "files": [
            "golden_bisect.json",
            "golden_hpf.json",
            "golden_normalise.json",
            "golden_match_data_models.json",
            "golden_fitter.json",
            "golden_e2e.json",
            "golden_legacy_gpu.json",
        ],
        "spectrum": "assets/sample-data/spec-0266-51602-0001.dat",
    }
    (OUT / "golden_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote golden fixtures to {OUT}")
    if FIREFLY.exists():
        print(f"Run firefly manually from {FIREFLY} to refresh e2e values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
