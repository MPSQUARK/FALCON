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

- `convert_input.py` — optional FITS → FALCON-native formats
- `convert_models.py` — model asset inspection/conversion helpers
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

- **Unit tests** — one class per firefly function (`StatisticsMath`, fitter, dust, mass-loss ISD)
- **Golden JSON** in `tests/Falcon.Tests/TestData/` — strict rtol ~1e-5
- **Tiers:** `dotnet test` (CPU, always) vs `--filter Category=Integration` (GPU + MaStar FITS assets)
- **Target:** 100+ tests (currently 180+)
- Regenerate refs: `python Scripts/generate_golden_refs.py`

## Firefly parity (north star)

Canonical reference: `C:\Users\marce\Repos\firefly\firefly_fitter.py` and `firefly_models.py`.

| Behavior | Firefly | FALCON |
|----------|---------|--------|
| χ² per eval | Precomputed `chi_models` + `np.dot(weights, chi_models)` | `FloatChiSquaredEvaluator` (default) or batched `GpuChiSquaredEvaluator` |
| Mix | top-10 × full `final_fit_list` × 8 ratios | Same (`IterativeSpectralFitter.Mix`) |
| Iterate cap | `fit_cap` commented out (no cap) | `fit_per_iteration_cap` enforced (default 1000) |
| Dust | `determine_attenuation` → `fitter`, then main `fitter` | `DustService.DetermineAttenuation` + orchestrator fit |

Parity checks: set `FALCON_FIT_TRACE=1`, run FALCON and `Scripts/run_firefly_fit_trace.py`, then `Scripts/compare_fit_traces.py`. FIT_TRACE emits `fitter_summary` with `index_count`, `iterate_count`, `mix_eval_count`.

Config: `fitting.use_gpu_evaluator` (default `false`), `fitting.experimental_mix_search` (default `false`).

## Fitting engines

- `fitting.engine: firefly` (default) — `FireflyFitOrchestrator` full pipeline
- `fitting.engine: legacy_gpu` — original `GenerateInterpSSP` GPU chi² path
- `fitting.legacy_gpu_fallback: true` — fall back to GPU path if firefly path fails
