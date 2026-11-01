# Create and populate the FALCON Python virtual environment for Scripts/.
$ErrorActionPreference = "Stop"
$RepoRoot = Split-Path -Parent $PSScriptRoot
$VenvDir = Join-Path $RepoRoot ".venv"
$Python = Get-Command python -ErrorAction SilentlyContinue
if (-not $Python) {
    $Python = Get-Command py -ErrorAction SilentlyContinue
    if ($Python) {
        & py -3 -m venv $VenvDir
    } else {
        throw "Python 3 not found on PATH. Install Python 3 and retry."
    }
} else {
    & python -m venv $VenvDir
}

$Pip = Join-Path $VenvDir "Scripts\python.exe"
& $Pip -m pip install --upgrade pip
& $Pip -m pip install -r (Join-Path $RepoRoot "Scripts\requirements.txt")
Write-Host "Virtual environment ready at $VenvDir"
Write-Host "Activate with: .\.venv\Scripts\Activate.ps1"
