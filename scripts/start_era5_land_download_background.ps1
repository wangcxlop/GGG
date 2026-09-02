#Requires -Version 7

[CmdletBinding()]
param(
    [string]$Python = '',
    [switch]$LocalValidationOnly,
    [switch]$Worker,
    [string]$ElevationProbePath = ''
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$DownloadScript = Join-Path $PSScriptRoot 'download_era5_land_stations.jl'
$ProtectScript = Join-Path $PSScriptRoot 'protect_data.ps1'
$StationMeta = Join-Path $ProjectRoot 'data\processed\study_area\station_meta.csv'
$RawDir = Join-Path $ProjectRoot 'data\raw\era5_land\stations_2022_2024'
$AuditDir = Join-Path $ProjectRoot 'output\input_audit\era5_land_stations_2022_2024'
$LogDir = Join-Path $ProjectRoot 'output\logs'
$LogPath = Join-Path $LogDir 'era5_land_station_download.log'

function Test-IsElevated {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = [Security.Principal.WindowsPrincipal]::new($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Start-UnelevatedLauncher {
    $arguments = [Collections.Generic.List[string]]::new()
    $arguments.Add('pwsh.exe')
    $arguments.Add('-NoProfile')
    $arguments.Add('-NoExit')
    $arguments.Add('-File')
    $arguments.Add(('"{0}"' -f $PSCommandPath))
    if (-not [string]::IsNullOrWhiteSpace($Python)) {
        $arguments.Add('-Python')
        $arguments.Add(('"{0}"' -f $Python))
    }
    if ($LocalValidationOnly) {
        $arguments.Add('-LocalValidationOnly')
    }

    & runas.exe /trustlevel:0x20000 ($arguments -join ' ')
    if ($LASTEXITCODE -ne 0) {
        throw "Could not relaunch with a Basic User token (runas exit code $LASTEXITCODE)."
    }
    Write-Host 'Relaunched the ERA5 downloader with a Basic User token.'
}

function Resolve-PythonPath {
    param([string]$Requested)

    $candidates = [Collections.Generic.List[string]]::new()
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        $candidates.Add($Requested)
    }
    if (-not [string]::IsNullOrWhiteSpace($env:LOCALAPPDATA)) {
        $candidates.Add((Join-Path $env:LOCALAPPDATA 'Programs\Python\Python313\python.exe'))
        $candidates.Add((Join-Path $env:LOCALAPPDATA 'Programs\Python\Python312\python.exe'))
        $candidates.Add((Join-Path $env:LOCALAPPDATA 'Programs\Python\Python311\python.exe'))
    }
    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if ($null -ne $pythonCommand) {
        $candidates.Add($pythonCommand.Source)
    }

    foreach ($candidate in $candidates | Select-Object -Unique) {
        if (-not (Test-Path -LiteralPath $candidate -PathType Leaf)) {
            continue
        }
        & $candidate -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)'
        if ($LASTEXITCODE -ne 0) {
            continue
        }
        & $candidate -c 'import cdsapi; assert hasattr(cdsapi.Client, "status")'
        if ($LASTEXITCODE -eq 0) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    throw 'No Python >= 3.9 installation with cdsapi is available. Pass -Python with a suitable executable.'
}

function Assert-FreeSpace {
    $root = [IO.Path]::GetPathRoot($ProjectRoot)
    $drive = Get-PSDrive -Name $root.Substring(0, 1)
    $minimum = 4GB
    if ($drive.Free -lt $minimum) {
        throw ('At least 4 GiB free space is required; only {0:N2} GiB is available.' -f ($drive.Free / 1GB))
    }
    Write-Host ('Disk space: {0:N2} GiB free' -f ($drive.Free / 1GB))
}

function Invoke-LocalPreflight {
    param([string]$PythonPath)

    Test-IsElevated | Where-Object { $_ } | ForEach-Object {
        throw 'Run this launcher from a normal, non-administrator PowerShell window.'
    }
    foreach ($path in @($DownloadScript, $ProtectScript, $StationMeta)) {
        Test-Path -LiteralPath $path -PathType Leaf |
            Where-Object { -not $_ } |
            ForEach-Object { throw "Required file is missing: $path" }
    }
    $null = Get-Command julia -ErrorAction Stop
    Assert-FreeSpace

    & $PythonPath -c 'import sys, cdsapi; print(f"Python {sys.version.split()[0]}; cdsapi import OK")'
    if ($LASTEXITCODE -ne 0) {
        throw 'Python/cdsapi preflight failed.'
    }

    $tempRoot = Join-Path ([IO.Path]::GetTempPath()) ("era5-local-preflight-{0}-{1}" -f $PID, [Guid]::NewGuid().ToString('N'))
    $tempAudit = Join-Path $tempRoot 'audit'
    $tempRaw = Join-Path $tempRoot 'raw'
    try {
        New-Item -ItemType Directory -Path $tempRoot | Out-Null
        & julia "--project=$ProjectRoot" $DownloadScript --dry-run `
            --python $PythonPath --station-meta $StationMeta `
            --raw-dir $tempRaw --audit-dir $tempAudit
        if ($LASTEXITCODE -ne 0) {
            throw "ERA5 Julia dry-run failed with exit code $LASTEXITCODE."
        }
        $manifestPath = Join-Path $tempAudit 'era5_land_station_manifest.csv'
        $manifest = @(Import-Csv -LiteralPath $manifestPath)
        $requests = @(Get-ChildItem -LiteralPath (Join-Path $tempAudit 'requests') -File -Filter '*.json')
        if ($manifest.Count -ne 237 -or $requests.Count -ne 237) {
            throw "Expected 237 manifest rows and requests; found $($manifest.Count) and $($requests.Count)."
        }
        Write-Host 'Local dry-run: 237 stations and 237 request files validated'
    }
    finally {
        $resolvedTemp = [IO.Path]::GetFullPath($tempRoot)
        $resolvedSystemTemp = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
        if ($resolvedTemp.StartsWith($resolvedSystemTemp, [StringComparison]::OrdinalIgnoreCase) -and
            (Test-Path -LiteralPath $resolvedTemp)) {
            Remove-Item -LiteralPath $resolvedTemp -Recurse -Force
        }
    }

    $complete = 0
    if (Test-Path -LiteralPath $RawDir) {
        $complete = @(Get-ChildItem -LiteralPath $RawDir -Directory | Where-Object {
            (Test-Path -LiteralPath (Join-Path $_.FullName '.complete')) -and
            @(Get-ChildItem -LiteralPath $_.FullName -File -Filter '*.zip').Count -eq 1
        }).Count
    }
    Write-Host "Resume state: $complete/237 station downloads already complete"
    Write-Host 'Local preflight passed; no ERA5 data was downloaded or modified'
}

function Assert-DataProtectionAvailable {
    & pwsh -NoProfile -File $ProtectScript -Lock
    if ($LASTEXITCODE -ne 0) {
        throw 'The data delete-lock could not be applied. Formal download is blocked.'
    }
    Write-Host 'Data delete-lock: available and active'
}

function Test-CdsCredentialAndLicense {
    param([string]$PythonPath)

    $tempRoot = Join-Path ([IO.Path]::GetTempPath()) ("era5-cds-check-{0}-{1}" -f $PID, [Guid]::NewGuid().ToString('N'))
    try {
        New-Item -ItemType Directory -Path $tempRoot | Out-Null
        $code = @'
import os
import pathlib
import sys
import cdsapi

target = pathlib.Path(sys.argv[1]).resolve()
client = cdsapi.Client(
    url="https://cds.climate.copernicus.eu/api",
    key=os.environ["CDSAPI_TOKEN"],
    quiet=True,
)
client.status()
previous = pathlib.Path.cwd()
try:
    os.chdir(target)
    client.retrieve(
        "reanalysis-era5-land-timeseries",
        {
            "variable": ["2m_temperature"],
            "location": {"longitude": 110.312, "latitude": 31.25},
            "date": ["2022-01-01/2022-01-01"],
            "data_format": "csv",
        },
    ).download()
finally:
    os.chdir(previous)

assets = [path for path in target.rglob("*") if path.is_file() and path.stat().st_size > 0]
if not assets:
    raise RuntimeError("CDS validation returned no non-empty asset")
print("CDS credential, licence, API and sample retrieval validated")
'@
        & $PythonPath -c $code $tempRoot
        if ($LASTEXITCODE -ne 0) {
            throw "CDS online validation failed with exit code $LASTEXITCODE. Formal download was not started."
        }
    }
    finally {
        $resolvedTemp = [IO.Path]::GetFullPath($tempRoot)
        $resolvedSystemTemp = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
        if ($resolvedTemp.StartsWith($resolvedSystemTemp, [StringComparison]::OrdinalIgnoreCase) -and
            (Test-Path -LiteralPath $resolvedTemp)) {
            Remove-Item -LiteralPath $resolvedTemp -Recurse -Force
        }
    }
}

function Test-AllRawDownloads {
    param([string]$PythonPath)

    $code = @'
import csv
import pathlib
import sys
import zipfile

station_meta = pathlib.Path(sys.argv[1])
raw_root = pathlib.Path(sys.argv[2])
with station_meta.open("r", encoding="utf-8-sig", newline="") as stream:
    station_ids = [row["station_id"] for row in csv.DictReader(stream)]

expected_headers = {
    ("valid_time", "d2m", "t2m", "latitude", "longitude"),
    ("valid_time", "sp", "latitude", "longitude"),
    ("valid_time", "u10", "v10", "latitude", "longitude"),
}
errors = []
for station_id in station_ids:
    station_dir = raw_root / station_id
    archives = list(station_dir.glob("*.zip")) if station_dir.is_dir() else []
    if not (station_dir / ".complete").is_file() or len(archives) != 1:
        errors.append(f"{station_id}: missing marker or single ZIP")
        continue
    try:
        with zipfile.ZipFile(archives[0]) as archive:
            bad_member = archive.testzip()
            if bad_member is not None:
                errors.append(f"{station_id}: CRC failure in {bad_member}")
                continue
            csv_names = [name for name in archive.namelist() if name.lower().endswith(".csv")]
            if len(csv_names) != 3:
                errors.append(f"{station_id}: expected 3 CSV files, found {len(csv_names)}")
                continue
            headers = set()
            for name in csv_names:
                with archive.open(name) as member:
                    header = member.readline().decode("utf-8-sig").strip().split(",")
                    headers.add(tuple(header))
            if headers != expected_headers:
                errors.append(f"{station_id}: unexpected CSV headers")
    except Exception as exc:
        errors.append(f"{station_id}: {exc}")

if errors:
    print("\n".join(errors[:25]), file=sys.stderr)
    raise SystemExit(f"Raw ERA5 validation failed for {len(errors)} station(s)")
print(f"Raw ERA5 validation passed for {len(station_ids)} stations")
'@
    & $PythonPath -c $code $StationMeta $RawDir
    if ($LASTEXITCODE -ne 0) {
        throw "Final raw-download validation failed with exit code $LASTEXITCODE."
    }
}

function Invoke-Worker {
    param([string]$PythonPath)

    if (Test-IsElevated) {
        throw 'The ERA5 worker refuses to run elevated because the data delete-lock would be ineffective.'
    }
    if ([string]::IsNullOrWhiteSpace($env:CDSAPI_TOKEN)) {
        throw 'The background worker did not inherit CDSAPI_TOKEN.'
    }

    New-Item -ItemType Directory -Path $LogDir -Force | Out-Null
    $transcriptStarted = $false
    try {
        Start-Transcript -Path $LogPath -Append | Out-Null
        $transcriptStarted = $true
        Write-Output "`n[$(Get-Date -Format s)] ERA5 background worker starting"

        & pwsh -NoProfile -File $ProtectScript -Unlock -RawOnly
        if ($LASTEXITCODE -ne 0) {
            throw "Could not unlock protected data paths (exit code $LASTEXITCODE)."
        }

        & julia "--project=$ProjectRoot" $DownloadScript `
            --python $PythonPath --station-meta $StationMeta `
            --raw-dir $RawDir --audit-dir $AuditDir
        if ($LASTEXITCODE -ne 0) {
            throw "ERA5 download failed with Julia exit code $LASTEXITCODE."
        }

        Test-AllRawDownloads -PythonPath $PythonPath
        Write-Output "[$(Get-Date -Format s)] ERA5 background download and validation completed"
    }
    catch {
        Write-Error $_
        throw
    }
    finally {
        Remove-Item Env:CDSAPI_TOKEN -ErrorAction SilentlyContinue
        & pwsh -NoProfile -File $ProtectScript -Lock -RawOnly
        if ($transcriptStarted) {
            Stop-Transcript | Out-Null
        }
    }
}

if (-not [string]::IsNullOrWhiteSpace($ElevationProbePath)) {
    $probeParent = Split-Path -Parent $ElevationProbePath
    if (-not (Test-Path -LiteralPath $probeParent -PathType Container)) {
        throw "Elevation probe parent does not exist: $probeParent"
    }
    Set-Content -LiteralPath $ElevationProbePath -Value (Test-IsElevated) -Encoding ascii
    exit 0
}

if (-not $Worker -and (Test-IsElevated)) {
    Start-UnelevatedLauncher
    exit 0
}

$PythonPath = Resolve-PythonPath -Requested $Python
if ($Worker) {
    Invoke-Worker -PythonPath $PythonPath
    exit 0
}

Invoke-LocalPreflight -PythonPath $PythonPath
Assert-DataProtectionAvailable
if ($LocalValidationOnly) {
    exit 0
}

$secureToken = $null
$tokenPointer = [IntPtr]::Zero
$plainToken = $null
try {
    $secureToken = Read-Host 'Copernicus CDS Personal Access Token' -AsSecureString
    $tokenPointer = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($secureToken)
    $plainToken = [Runtime.InteropServices.Marshal]::PtrToStringBSTR($tokenPointer)
    if ([string]::IsNullOrWhiteSpace($plainToken)) {
        throw 'The CDS token is empty.'
    }
    $env:CDSAPI_TOKEN = $plainToken

    Test-CdsCredentialAndLicense -PythonPath $PythonPath

    New-Item -ItemType Directory -Path $LogDir -Force | Out-Null
    $argumentList = @(
        '-NoProfile',
        '-File', $PSCommandPath,
        '-Worker',
        '-Python', $PythonPath
    )
    $process = Start-Process -FilePath 'pwsh' -ArgumentList $argumentList `
        -WindowStyle Hidden -PassThru
    Write-Host "Formal download started in the background (PID $($process.Id))."
    Write-Host "Progress log: $LogPath"
}
finally {
    Remove-Item Env:CDSAPI_TOKEN -ErrorAction SilentlyContinue
    $plainToken = $null
    if ($tokenPointer -ne [IntPtr]::Zero) {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($tokenPointer)
    }
    if ($null -ne $secureToken) {
        $secureToken.Dispose()
    }
}
