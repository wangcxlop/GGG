param(
    [switch]$Overwrite,
    [string]$OutputDir = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$juliaScript = Join-Path $PSScriptRoot "download_copernicus_glo30.jl"
$logDir = Join-Path $projectRoot "output\logs"
$logPath = Join-Path $logDir "copernicus_glo30_download.log"
New-Item -ItemType Directory -Path $logDir -Force | Out-Null

$accessKey = Read-Host "CDSE S3 access key"
$secureSecretKey = Read-Host "CDSE S3 secret key" -AsSecureString
$secretKeyPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secureSecretKey)
$transcriptStarted = $false

try {
    $plainSecretKey = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($secretKeyPointer)
    $env:CDSE_S3_ACCESS_KEY = $accessKey
    $env:CDSE_S3_SECRET_KEY = $plainSecretKey
    Start-Transcript -Path $logPath -Append | Out-Null
    $transcriptStarted = $true

    $arguments = @("--project=$projectRoot", $juliaScript)
    if ($Overwrite) {
        $arguments += "--overwrite"
    }
    if (-not [string]::IsNullOrWhiteSpace($OutputDir)) {
        $arguments += @("--output", $OutputDir)
    }

    Write-Output "`n[$(Get-Date -Format s)] Starting secure Copernicus GLO-30 download"

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
        throw "Copernicus GLO-30 download failed with Julia exit code $juliaExitCode"
    }
}
finally {
    if ($transcriptStarted) {
        Stop-Transcript | Out-Null
    }
    Remove-Item Env:CDSE_S3_ACCESS_KEY -ErrorAction SilentlyContinue
    Remove-Item Env:CDSE_S3_SECRET_KEY -ErrorAction SilentlyContinue
    $plainSecretKey = $null
    if ($secretKeyPointer -ne [IntPtr]::Zero) {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($secretKeyPointer)
    }
    $secureSecretKey.Dispose()
}
