# FALCON Roadmap

## M0 — Documentation & foundation ✅
Repo restructure, Falcon.Core/Cli/Tests/UI, BAVCL reference fix, DI + logging scaffold.

## M1 — BAVCL modernisation & resolution service ✅
Global BAVCL usings, Statistics kernel workload, `ResolutionService` in `Services/Instrument/`.

## M2 — Spectrum setup & YAML config ✅
`FireflyOptions`, YAML loader, `LegacyConfig.ApplyFrom`, asset path updates.

## M3 — Iterative fitter ✅
`IterativeSpectralFitter` porting firefly tree expansion + BIC gate.

## M4 — Dust, HPF, statistics, output 🔄
`DustService`, `StatisticsService`, SFD dust wired; FITS output + full HPF integration in pipeline ongoing.

## M5 — CLI polish & batch ✅
`Falcon.Cli` hosted entry, parallel batch via `SpectralFitRunner`.

## M6 — Falcon.UI 🔄
Blazor Web App scaffolded; wire config upload, run trigger, plot viewer.

## M7 — MaNGA + SFR 🔄
`MangaBatchReader`, `StarFormationRateService` stubs after core golden tests.

## M8 — Upstream firefly diff (deferred)
Compare local reference vs current upstream; apply incremental changes.
