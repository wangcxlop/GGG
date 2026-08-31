#Requires -Version 7
<#
.SYNOPSIS
    Lock / unlock the irreplaceable parts of data/ against deletion.

.DESCRIPTION
    Adds a deny-DELETE ACE for the current user to the data files that cannot be
    regenerated. Creating and overwriting files stays allowed; only deletion is
    blocked. See CLAUDE.md, "Data safety", for the incident this guards against.

    IMPORTANT -- this only works from a NON-ELEVATED session. A process running as
    Administrator holds SeBackupPrivilege/SeRestorePrivilege, which bypass the DACL
    entirely: measured on this machine, a delete succeeded against a deny ACE for
    Everyone AND Administrators. Run -Status to see whether the current session can
    be stopped by the lock at all.

    Protected: every entry directly under data/ except data/processed -- i.e.
    data/raw/, data/FY4B/ and the top-level source files -- plus a non-inherited
    lock on the data/ and data/processed/ directory objects themselves, so a
    recursive delete cannot remove the folders.

    NOT protected: the contents of data/processed/. Those files are regenerated
    from raw by the prepare_* pipelines, which rewrite them via
    write-temp -> mv(force) -> rm(temp) and therefore need DELETE rights.

.EXAMPLE
    pwsh -NoProfile -File scripts/protect_data.ps1 -Status
    pwsh -NoProfile -File scripts/protect_data.ps1 -Unlock   # before a re-download
    pwsh -NoProfile -File scripts/protect_data.ps1 -Lock     # after it finishes
#>
[CmdletBinding(DefaultParameterSetName = 'Status')]
param(
    [Parameter(ParameterSetName = 'Lock')]   [switch]$Lock,
    [Parameter(ParameterSetName = 'Unlock')] [switch]$Unlock,
    [Parameter(ParameterSetName = 'Status')] [switch]$Status
)

$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
$DataDir = Join-Path $RepoRoot 'data'
$ProcessedDir = Join-Path $DataDir 'processed'

if (-not (Test-Path -LiteralPath $DataDir)) {
    Write-Error "No data directory at $DataDir"
    exit 1
}

$WinIdentity = [System.Security.Principal.WindowsIdentity]::GetCurrent()
$Sid = $WinIdentity.User.Value
$Identity = "*$Sid"
$IsElevated = (New-Object System.Security.Principal.WindowsPrincipal($WinIdentity)).IsInRole(
    [System.Security.Principal.WindowsBuiltInRole]::Administrator)

function Write-ElevationWarning {
    if (-not $IsElevated) { return }
    Write-Host ''
    Write-Warning 'This session is ELEVATED (Administrator).'
    Write-Warning 'SeBackupPrivilege/SeRestorePrivilege bypass file ACLs, so the delete-lock'
    Write-Warning 'does NOT stop deletions made from this session -- verified by measurement,'
    Write-Warning 'not assumed. The lock only bites for a normal, non-elevated session.'
    Write-Warning 'Launch the terminal (and Claude Code) WITHOUT "Run as administrator" for'
    Write-Warning 'this layer to mean anything.'
    Write-Host ''
}

# Everything under data/ except processed/, which the pipelines must be able to rewrite.
function Get-ProtectedEntry {
    $entries = @(Get-ChildItem -LiteralPath $DataDir -Force |
        Where-Object { $_.Name -ne 'processed' })
    $entries | ForEach-Object {
        [pscustomobject]@{
            Path        = $_.FullName
            IsContainer = $_.PSIsContainer
            Inherited   = $true
        }
    }
}

# The two directory objects themselves: deny deleting the folder, without inheriting
# into the files the pipelines churn.
function Get-ShellEntry {
    $shells = @($DataDir)
    if (Test-Path -LiteralPath $ProcessedDir) { $shells += $ProcessedDir }
    $shells | ForEach-Object {
        [pscustomobject]@{ Path = $_; IsContainer = $true; Inherited = $false }
    }
}

function Invoke-Icacls {
    param([string[]]$IcaclsArgs)
    $output = & icacls @IcaclsArgs 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Warning "icacls failed for: $($IcaclsArgs -join ' ')"
        $output | ForEach-Object { Write-Warning "  $_" }
        return $false
    }
    return $true
}

function Set-DenyDelete {
    param([pscustomobject]$Entry)
    # DE = delete this object. DC = delete a child through the parent's right.
    # (OI)(CI) makes the ACE inherit to every descendant, which is what actually
    # saves the files inside a directory.
    $rights = if ($Entry.IsContainer -and $Entry.Inherited) { '(OI)(CI)(DE,DC)' } else { '(DE)' }
    Invoke-Icacls @($Entry.Path, '/deny', "${Identity}:$rights")
}

function Remove-DenyDelete {
    param([pscustomobject]$Entry)
    Invoke-Icacls @($Entry.Path, '/remove:d', $Identity)
}

function Test-Locked {
    param([string]$Path)
    $acl = & icacls $Path 2>&1
    return [bool]($acl | Select-String -SimpleMatch '(DENY)' -Quiet)
}

$all = @(Get-ProtectedEntry) + @(Get-ShellEntry)

switch ($PSCmdlet.ParameterSetName) {
    'Lock' {
        Write-Host "Locking $($all.Count) path(s) under $DataDir against deletion..."
        $failed = 0
        foreach ($e in $all) {
            if (-not (Set-DenyDelete $e)) { $failed++ }
        }
        Write-Host ''
        if ($failed -gt 0) {
            Write-Warning "$failed path(s) could not be locked -- see the warnings above."
        }
        Write-Host "Locked. Run with -Unlock before any data re-download: the download"
        Write-Host "scripts delete their own .partial temporaries under data/raw/."
        Write-ElevationWarning
    }
    'Unlock' {
        Write-Host "Unlocking $($all.Count) path(s) under $DataDir..."
        $failed = 0
        foreach ($e in $all) {
            if (-not (Remove-DenyDelete $e)) { $failed++ }
        }
        Write-Host ''
        if ($failed -gt 0) {
            Write-Warning "$failed path(s) could not be unlocked -- see the warnings above."
        }
        Write-Host "Unlocked. data/ is now deletable -- re-run with -Lock as soon as the"
        Write-Host "download or cleanup you needed this for is finished."
    }
    default {
        Write-Host "Delete-lock status for $DataDir"
        Write-Host "  identity : $Sid"
        Write-Host "  elevated : $IsElevated$(if ($IsElevated) { '   <-- lock is NOT enforced for this session' })"
        Write-Host ''
        $all | ForEach-Object {
            $state = if (Test-Locked $_.Path) { 'LOCKED  ' } else { 'unlocked' }
            $scope = if (-not $_.IsContainer) { 'file' }
                     elseif ($_.Inherited) { 'recursive' } else { 'folder only' }
            $rel = [IO.Path]::GetRelativePath($RepoRoot, $_.Path)
            Write-Host ("  {0}  {1,-11}  {2}" -f $state, $scope, $rel)
        }
        Write-Host ''
        Write-Host "data/processed/ contents are deliberately left deletable (regenerable"
        Write-Host "from raw, and the prepare_* pipelines rewrite them in place)."
        Write-ElevationWarning
    }
}
