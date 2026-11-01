# Run the same spectrum fit with FIT_TRACE enabled in FALCON and firefly, then diff.
param(
    [string]$Spectrum = "assets/sample-data/spec-0266-51630-0016.fits",
    [string]$Config = "configs/example.yaml",
    [string]$FireflyDir = "",
    [string]$OutDir = "output/fit-trace"
)

$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot
if (-not $FireflyDir) {
    $FireflyDir = Join-Path (Split-Path -Parent $RepoRoot) "firefly"
}

New-Item -ItemType Directory -Force -Path (Join-Path $RepoRoot $OutDir) | Out-Null
$FalconTrace = Join-Path $RepoRoot "$OutDir/falcon.trace"
$FireflyTrace = Join-Path $RepoRoot "$OutDir/firefly.trace"

Write-Host "=== FALCON fit trace ===" -ForegroundColor Cyan
$env:FALCON_FIT_TRACE = "1"
Push-Location $RepoRoot
try {
    dotnet run --project src/Falcon.Cli/Falcon.Cli.csproj -c Debug -- $Config $Spectrum 2>&1 |
        Tee-Object -FilePath $FalconTrace
}
finally {
    Pop-Location
    Remove-Item Env:FALCON_FIT_TRACE -ErrorAction SilentlyContinue
}

if (-not (Test-Path $FireflyDir)) {
    Write-Warning "Firefly repo not found at $FireflyDir — skipping firefly trace."
    Write-Host "FALCON trace saved to $FalconTrace"
    exit 0
}

Write-Host "`n=== firefly fit trace ===" -ForegroundColor Cyan
$env:FIREFLY_FIT_TRACE = "1"
$env:FF_DIR = $FireflyDir
$env:STELLARPOPMODELS_DIR = Join-Path $RepoRoot "assets/stellar-population-models"
$SpecPath = if ([System.IO.Path]::IsPathRooted($Spectrum)) { $Spectrum } else { Join-Path $RepoRoot $Spectrum }

Push-Location $FireflyDir
try {
    python (Join-Path $RepoRoot "Scripts/run_firefly_fit_trace.py") $SpecPath 2>&1 |
        Tee-Object -FilePath $FireflyTrace
}
catch {
    Write-Warning "firefly trace failed: $_"
    Write-Host "FALCON trace saved to $FalconTrace"
    exit 0
}
finally {
    Pop-Location
    Remove-Item Env:FIREFLY_FIT_TRACE -ErrorAction SilentlyContinue
}

Write-Host "`n=== compare ===" -ForegroundColor Cyan
python (Join-Path $RepoRoot "Scripts/compare_fit_traces.py") $FalconTrace $FireflyTrace
