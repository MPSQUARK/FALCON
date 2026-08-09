# FALCON Design

FALCON is a C# GPU-accelerated port of the [firefly](https://github.com/jjrbotteon/firefly) stellar population fitting algorithm, using [BAVCL](https://github.com/marcelbresink/BAVCL) for vectorised GPU numerics.

## Goals

- **Parity** with the local reference firefly copy at `C:\Users\marce\Repos\firefly` (MaStar models only)
- **Modern C#** — DI, interfaces, multi-project layout, YAML configuration
- **GPU performance** via BAVCL + ILGPU custom kernels
- **Testability** — component tests + bundled golden files for strict numeric regression

## Solution layout

```
src/Falcon.Core/     Domain logic, GPU kernels, firefly algorithm ports
src/Falcon.Cli/      Hosted console entry (`falcon` executable)
src/Falcon.UI/       Blazor Web App (config, run, plot viewer)
tests/Falcon.Tests/  xUnit component + golden tests
configs/             YAML firefly-equivalent configuration
assets/              Runtime data (models, dust maps, sample spectra)
scripts/             Python prep/post tools (manual invocation)
docs/                Design, roadmap, BAVCL feature requests
```

## Architecture

Hosts (`Falcon.Cli`, `Falcon.UI`) depend on `Falcon.Core` services registered via `ServiceCollectionExtensions.AddFalconCore()`.

Pipeline flow:

1. Load spectrum (FITS via `FitsSpectrumReader`)
2. Setup masks, MW dust (`DustService`)
3. Load/downgrade MaStar models (`SpectralFitPipeline` + `ResolutionService`)
4. Normalise (`StatisticsService`)
5. HPF + dust attenuation (`DustService` + `IterativeSpectralFitter`)
6. Full iterative fit (`IterativeSpectralFitter`)
7. PDF statistics + output (`StatisticsService`, FITS/text writers)

## Configuration

- **YAML** (`configs/*.yaml`) — all firefly inline options + data paths
- **appsettings.json** — logging only (Microsoft.Extensions.Logging)

## Python bridge

Scripts in `scripts/` are run manually:

- `convert_input.py` — FITS/HDF5 → FALCON-native formats
- `convert_models.py` — model asset conversion (evolves `ConvertDustMapToHDF5.py`)
- `convert_output.py` — results → firefly-compatible FITS
- `plot_results.py` — visualisation for UI

## BAVCL relationship

FALCON uses BAVCL as the GPU vector platform. Domain-specific spectral code (convolution, HPF, fitter) lives in FALCON until BAVCL adds APIs — see `docs/BAVCL_FEATURE_REQUESTS.md`.

## Legacy mapping

| Old root file | New location |
|---------------|--------------|
| `Spectral_Model.cs` | `Pipeline/SpectralFitPipeline.cs` |
| `Config.cs` | `Configuration/LegacyConfig.cs` + `FireflyOptions` |
| `UtilityMethods.cs` | `Services/IO/FitsSpectrumReader.cs` |
| `Spectral_resolution.cs` | `Services/Instrument/ResolutionService.cs` |

## Test strategy

- **Component tests** per service (fitter, dust, statistics)
- **Golden files** in `tests/Falcon.Tests/TestData/` — strict rtol ~1e-5 vs reference Python outputs
- GPU tests use BAVCL default backend; CPU fallback where GPU unavailable in CI
