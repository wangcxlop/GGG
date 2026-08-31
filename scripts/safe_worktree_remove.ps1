#Requires -Version 7
<#
.SYNOPSIS
    Remove a git worktree without following junctions out of it.

.DESCRIPTION
    'git worktree remove' recurses into NTFS junctions on Windows. A junction
    inside a worktree pointing at the main checkout's data/ is how the dataset was
    destroyed once before (see CLAUDE.md, "Data safety").

    This wrapper enumerates every reparse point inside the worktree WITHOUT
    following it, unlinks each one with 'cmd /c rmdir' (which removes the link and
    never touches the target), verifies none remain, and only then hands off to
    'git worktree remove'.

.EXAMPLE
    pwsh -NoProfile -File scripts/safe_worktree_remove.ps1 .claude/worktrees/my-branch
    pwsh -NoProfile -File scripts/safe_worktree_remove.ps1 .claude/worktrees/my-branch -Force
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory, Position = 0)]
    [string]$Path,

    # Passed through to 'git worktree remove' for a worktree with local changes.
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot

if (-not (Test-Path -LiteralPath $Path)) {
    Write-Error "No such path: $Path"
    exit 1
}
$Full = (Resolve-Path -LiteralPath $Path).Path

# Refuse anything git does not know as a worktree, and refuse the main checkout.
$registered = @(& git -C $RepoRoot worktree list --porcelain |
    Select-String -Pattern '^worktree (.+)$' |
    ForEach-Object { $_.Matches[0].Groups[1].Value -replace '/', '\' })

if ($registered.Count -eq 0) {
    Write-Error "git worktree list returned nothing -- is $RepoRoot a git repository?"
    exit 1
}
$mainWorktree = $registered[0]

if ($Full -ieq $mainWorktree) {
    Write-Error "Refusing: $Full is the MAIN checkout, not a worktree. This is where data/ lives."
    exit 1
}
if ($registered -notcontains $Full) {
    Write-Error "Refusing: $Full is not a registered worktree of $RepoRoot.`nRegistered: $($registered -join ', ')"
    exit 1
}

# -Recurse does NOT follow junctions/symlinks without -FollowSymlink, so this walk
# stays inside the worktree.
function Get-Link {
    @(Get-ChildItem -LiteralPath $Full -Recurse -Force -Attributes ReparsePoint -ErrorAction SilentlyContinue)
}

$links = Get-Link
if ($links.Count -eq 0) {
    Write-Host "No reparse points inside $Full."
}
else {
    Write-Host "Found $($links.Count) reparse point(s) inside the worktree. Unlinking (targets are NOT touched):"
    foreach ($link in $links) {
        $target = if ($link.LinkTarget) { $link.LinkTarget } else { '<unresolved>' }
        Write-Host "  $($link.FullName)  ->  $target"
        if ($link.PSIsContainer) {
            # rmdir on a junction/directory symlink removes the link only.
            & cmd /c rmdir "$($link.FullName)"
        }
        else {
            & cmd /c del /q "$($link.FullName)"
        }
        if (Test-Path -LiteralPath $link.FullName) {
            Write-Error "Failed to unlink $($link.FullName) -- aborting before git touches this tree."
            exit 1
        }
        if ($target -ne '<unresolved>' -and -not (Test-Path -LiteralPath $target)) {
            Write-Error "The link target $target disappeared while unlinking. STOP and check it before continuing."
            exit 1
        }
        Write-Host "    unlinked; target still present"
    }
}

$remaining = Get-Link
if ($remaining.Count -gt 0) {
    Write-Error "Reparse points still present after unlinking -- refusing to run git worktree remove."
    exit 1
}

Write-Host ''
Write-Host "Tree is link-free. Running git worktree remove..."
$gitArgs = @('-C', $RepoRoot, 'worktree', 'remove')
if ($Force) { $gitArgs += '--force' }
$gitArgs += $Full
& git @gitArgs
if ($LASTEXITCODE -ne 0) {
    Write-Error "git worktree remove failed with exit code $LASTEXITCODE."
    exit $LASTEXITCODE
}
Write-Host "Removed $Full."
