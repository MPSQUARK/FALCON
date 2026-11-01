# FALCON Roadmap

## M0 — Documentation & foundation (done)
Repo restructure, Falcon.Core/Cli/Tests/UI, BAVCL reference, DI + logging.

## M1 — BAVCL modernisation & resolution service (done)
Global BAVCL usings, `ResolutionService`, downgrade kernels.

## M2 — Spectrum setup & YAML config (done)
`FireflyOptions` with `engine`, `legacy_gpu_fallback`, `hpf_mode` (on/off/hpf_only), emission-line list.

## M3 — Iterative fitter (done)
`IterativeSpectralFitter` as primary fit path via `FireflyFitOrchestrator`.

## M4 — Dust, HPF, statistics, output (done)
- Full `StatisticsMath` / `StatisticsService` PDF pipeline
- `FireflyFitOrchestrator` wires match → normalise → dust → fit → stats
- `FireflyFitsWriter` JSON/HIERARCH sidecar output
- `GenerateInterpSSP` retained as `legacy_gpu` engine + fallback

## M5 — CLI polish & batch (done)
`Falcon.Cli` hosted entry, `SpectralFitRunner` with full DI.

## M6 — Falcon.UI (scaffolded)
Blazor Web App; config/run/plot pages stubbed.

## M7 — MaNGA + SFR (deferred)
Stubs remain until integration golden tests pass with the real MaStar FITS asset.

## M8 — Upstream firefly diff (deferred)

## Testing
- **177+ unit/golden tests** (`dotnet test tests/Falcon.Tests`)
- CI: all CPU + golden JSON tests
- Integration: `dotnet test --filter Category=Integration` (requires MaStar FITS + GPU)
