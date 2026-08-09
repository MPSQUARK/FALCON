#!/usr/bin/env python3
"""Convert MaStar FITS / dust FITS to HDF5 cache for FALCON."""
import argparse
import os
import sys

# Delegate to existing dust converter when available
HERE = os.path.dirname(__file__)
LEGACY = os.path.join(HERE, "ConvertDustMapToHDF5.py")
if os.path.exists(LEGACY):
    import runpy
    sys.argv = ["convert_models.py", *sys.argv[1:]]
    runpy.run_path(LEGACY, run_name="__main__")
else:
    print("Place ConvertDustMapToHDF5.py in scripts/ or implement MaStar conversion here.")
