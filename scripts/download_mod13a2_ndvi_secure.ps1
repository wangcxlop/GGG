param(
    [switch]$SubmitOnly,
    [string]$TaskId = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$juliaScript = Join-Path $PSScriptRoot "download_mod13a2_ndvi.jl"
$logDir = Join-Path $projectRoot "output\logs"
$logPath = Join-Path $logDir "mod13a2_ndvi_download.log"
New-Item -ItemType Directory -Path $logDir -Force | Out-Null

$username = Read-Host "NASA Earthdata username"
$securePassword = Read-Host "NASA Earthdata password" -AsSecureString
$passwordPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($securePassword)
$transcriptStarted = $false

try {
    $plainPassword = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($passwordPointer)
    $env:EARTHDATA_USERNAME = $username
    $env:EARTHDATA_PASSWORD = $plainPassword
    Start-Transcript -Path $logPath -Append | Out-Null
    $transcriptStarted = $true

    $arguments = @("--project=$projectRoot", $juliaScript)
    if ($SubmitOnly) {
        $arguments += "--submit-only"
    }
    if (-not [string]::IsNullOrWhiteSpace($TaskId)) {
        $arguments += @("--task-id", $TaskId)
    }

    Write-Output "`n[$(Get-Date -Format s)] Starting secure AppEEARS download"

    # Windows PowerShell converts native stderr into error records. Keep those
    # messages visible/logged without allowing retry notices to terminate Julia.
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
        throw "NDVI download failed with Julia exit code $juliaExitCode"
    }
}
finally {
    if ($transcriptStarted) {
        Stop-Transcript | Out-Null
    }
    Remove-Item Env:EARTHDATA_USERNAME -ErrorAction SilentlyContinue
    Remove-Item Env:EARTHDATA_PASSWORD -ErrorAction SilentlyContinue
    $plainPassword = $null
    if ($passwordPointer -ne [IntPtr]::Zero) {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($passwordPointer)
    }
    $securePassword.Dispose()
}
