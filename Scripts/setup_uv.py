#!/usr/bin/env python3
"""Create or update the project-local uv environment."""
from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
subprocess.run([sys.executable, "-m", "uv", "sync"], cwd=SCRIPT_DIR, check=True)
print(f"FALCON Python environment ready at {SCRIPT_DIR / '.venv'}")
