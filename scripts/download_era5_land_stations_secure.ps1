param(
    [switch]$DryRun,
    [int]$Limit = 0,
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$juliaScript = Join-Path $PSScriptRoot "download_era5_land_stations.jl"
$logDir = Join-Path $projectRoot "output\logs"
$logPath = Join-Path $logDir "era5_land_station_download.log"
New-Item -ItemType Directory -Path $logDir -Force | Out-Null
$arguments = @("--project=$projectRoot", $juliaScript, "--python", $Python)

if ($DryRun) {
    $arguments += "--dry-run"
}
if ($Limit -gt 0) {
    $arguments += @("--limit", $Limit.ToString())
}

if (-not $DryRun) {
    & $Python -c "import cdsapi"
    if ($LASTEXITCODE -ne 0) {
        throw 'Python package cdsapi is missing. Run: python -m pip install "cdsapi>=0.7.7"'
    }
}

$secureToken = $null
$tokenPointer = [IntPtr]::Zero
$transcriptStarted = $false
try {
    if (-not $DryRun) {
        $secureToken = Read-Host "Copernicus CDS Personal Access Token" -AsSecureString
        $tokenPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secureToken)
        $plainToken = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($tokenPointer)
        $env:CDSAPI_TOKEN = $plainToken
    }

    Start-Transcript -Path $logPath -Append | Out-Null
    $transcriptStarted = $true
    Write-Output "`n[$(Get-Date -Format s)] Starting ERA5-Land station download"

    $previousErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        & julia @arguments
        $juliaExitCode = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $previousErrorActionPreference
    }
    if ($juliaExitCode -ne 0) {
        throw "ERA5-Land download failed with Julia exit code $juliaExitCode"
    }
}
finally {
    if ($transcriptStarted) {
        Stop-Transcript | Out-Null
    }
    Remove-Item Env:CDSAPI_TOKEN -ErrorAction SilentlyContinue
    $plainToken = $null
    if ($tokenPointer -ne [IntPtr]::Zero) {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($tokenPointer)
    }
    if ($null -ne $secureToken) {
        $secureToken.Dispose()
    }
}
