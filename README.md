# FALCON

GPU-accelerated C# port of the firefly stellar population fitting algorithm.

## Quick start

```bash
# Build
dotnet build FALCON.sln

# Run CLI (requires MaStar FITS and SFD FITS assets)
dotnet run --project src/Falcon.Cli -- configs/example.yaml

# Run tests
dotnet test tests/Falcon.Tests

# Run Blazor UI
dotnet run --project src/Falcon.UI
```

## Layout

See [docs/DESIGN.md](docs/DESIGN.md) and [docs/ROADMAP.md](docs/ROADMAP.md).

## Python scripts

FALCON ships helper scripts under `Scripts/`. Use a project-local virtual environment:

```powershell
# Windows
.\Scripts\setup_venv.ps1
.\.venv\Scripts\Activate.ps1
python Scripts\plot_results.py output\spec-0266-51630-0016.results.json
```

```bash
# Linux / macOS
./Scripts/setup_venv.sh
source .venv/bin/activate
python Scripts/plot_results.py output/spec-0266-51630-0016.results.json
```

Dependencies are listed in `Scripts/requirements.txt`.

## Data assets

Place large files (gitignored):

- `assets/stellar-population-models/MaStar_SSP_v0.2.fits.gz`
- `assets/dust-maps/SFD_dust_4096_ngp.fits`
- `assets/dust-maps/SFD_dust_4096_sgp.fits`

These files are the native FITS assets from Firefly. Falcon reads them directly with CFITS.

## Configuration

Copy and edit `configs/example.yaml` — mirrors firefly.py inline options.

Logging is configured in `src/Falcon.Cli/appsettings.json`.
