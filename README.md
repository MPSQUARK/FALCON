# FALCON

GPU-accelerated C# port of the firefly stellar population fitting algorithm.

## Quick start

```bash
# Build
dotnet build FALCON.sln

# Run CLI (requires MaStar_SSP.h5 and dust.h5 in assets/)
dotnet run --project src/Falcon.Cli -- configs/example.yaml

# Run tests
dotnet test tests/Falcon.Tests

# Run Blazor UI
dotnet run --project src/Falcon.UI
```

## Layout

See [docs/DESIGN.md](docs/DESIGN.md) and [docs/ROADMAP.md](docs/ROADMAP.md).

## Data assets

Place large files (gitignored):

- `assets/stellar-population-models/MaStar_SSP.h5`
- `assets/dust-maps/dust.h5`

Convert dust maps: `python scripts/ConvertDustMapToHDF5.py`

## Configuration

Copy and edit `configs/example.yaml` — mirrors firefly.py inline options.

Logging is configured in `src/Falcon.Cli/appsettings.json`.
