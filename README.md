# FALCON

GPU-accelerated C# port of the firefly stellar population fitting algorithm.

## Quick start

```bash
# Build
dotnet build FALCON.sln

# Run CLI (requires MaStar FITS and SFD FITS assets)
dotnet run --project Falcon.Cli -- Configs/example.yaml

# Run tests
dotnet test Falcon.Tests

# Run Blazor UI
dotnet run --project Falcon.UI
```

## Layout

See [Documentation/DESIGN.md](Documentation/DESIGN.md) and [Documentation/ROADMAP.md](Documentation/ROADMAP.md).

## Python scripts

FALCON ships helper scripts under `Scripts/`. Use a project-local virtual environment:

```powershell
# Windows
python Scripts/setup_uv.py
.\Scripts\.venv\Scripts\Activate.ps1
python Scripts\plot_results.py output\spec-0266-51630-0016.results.json
```

```bash
# Linux / macOS
python Scripts/setup_uv.py
source Scripts/.venv/bin/activate
python Scripts/Analysis/plot_results.py output/spec-0266-51630-0016.results.json
```

Dependencies are listed in `Scripts/pyproject.toml`.

## Data assets

Place large files (gitignored):

- `Assets/StellarPopulationModels/MaStar_SSP_v0.2.fits.gz`
- `Assets/DustMaps/SFD_dust_4096_ngp.fits`
- `Assets/DustMaps/SFD_dust_4096_sgp.fits`

These files are the native FITS assets from Firefly. Falcon reads them directly with CFITS.

## Configuration

Copy and edit `Configs/example.yaml` — mirrors firefly.py inline options.

Logging is configured in `Falcon.Cli/appsettings.json`.

