#Requires -Version 7
<#
.SYNOPSIS
    Blocks the command shapes that destroyed data/ once before.

.DESCRIPTION
    Wired up from .claude/settings.json as a PreToolUse hook (and a SessionStart
    tripwire). See the "Data safety" section of CLAUDE.md for the incident this
    guards against. Reads the hook payload as JSON on stdin; prints a deny
    decision on stdout and exits 0, or prints nothing and exits 0 to allow.
#>
[CmdletBinding()]
param(
    [ValidateSet('PreToolUse', 'Tripwire')]
    [string]$Mode = 'PreToolUse'
)

$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)

function Get-ReparsePointList {
    param([string]$Root)
    if ([string]::IsNullOrWhiteSpace($Root) -or -not (Test-Path -LiteralPath $Root)) { return @() }
    # -Recurse does NOT follow junctions/symlinks without -FollowSymlink, so this
    # enumeration is itself safe to run over a tree that contains one.
    @(Get-ChildItem -LiteralPath $Root -Recurse -Force -Attributes ReparsePoint -ErrorAction SilentlyContinue)
}

function Deny {
    param([string]$Reason)
    @{
        hookSpecificOutput = @{
            hookEventName            = 'PreToolUse'
            permissionDecision       = 'deny'
            permissionDecisionReason = $Reason
        }
    } | ConvertTo-Json -Depth 5 -Compress
    exit 0
}

if ($Mode -eq 'Tripwire') {
    $warnings = @()

    $links = Get-ReparsePointList $RepoRoot
    if ($links.Count -gt 0) {
        $list = ($links | ForEach-Object { $_.FullName }) -join '; '
        $warnings += "reparse point(s) present in this repo: $list -- a junction/symlink inside the repo is exactly what destroyed data/ before. Unlink it with 'cmd /c rmdir <link>' (never Remove-Item) before any worktree cleanup."
    }

    # An elevated session holds SeBackup/SeRestorePrivilege, which bypass file ACLs,
    # so scripts/protect_data.ps1's delete-lock cannot stop anything here.
    $winId = [System.Security.Principal.WindowsIdentity]::GetCurrent()
    $isElevated = (New-Object System.Security.Principal.WindowsPrincipal($winId)).IsInRole(
        [System.Security.Principal.WindowsBuiltInRole]::Administrator)
    if ($isElevated) {
        $warnings += "this session is ELEVATED (Administrator), so the data/ delete-lock is bypassed by SeRestorePrivilege and protects nothing. Restart Claude Code from a normal, non-elevated terminal before doing anything that deletes files."
    }

    if ($warnings.Count -gt 0) {
        @{ systemMessage = "DATA SAFETY: " + ($warnings -join ' || ') } | ConvertTo-Json -Compress
    }
    exit 0
}

$raw = [Console]::In.ReadToEnd()
if ([string]::IsNullOrWhiteSpace($raw)) { exit 0 }
try { $payload = $raw | ConvertFrom-Json } catch { exit 0 }

$tool = [string]$payload.tool_name

# The harness's own worktree teardown: refuse while a link is still inside the tree.
if ($tool -eq 'ExitWorktree') {
    $cwd = if ($payload.cwd) { [string]$payload.cwd } else { (Get-Location).Path }
    $links = Get-ReparsePointList $cwd
    if ($links.Count -gt 0) {
        $list = ($links | ForEach-Object { $_.FullName }) -join '; '
        Deny "Worktree teardown blocked: this worktree still contains reparse point(s): $list. Removing the worktree would follow the link and delete whatever it points at (this is how data/ was wiped before). Unlink each one with 'cmd /c rmdir <link>' first, or run scripts/safe_worktree_remove.ps1."
    }
    exit 0
}

$cmd = [string]$payload.tool_input.command
if ([string]::IsNullOrWhiteSpace($cmd)) { exit 0 }

# 1. Link creation -- the exact mechanism of the incident.
$linkPatterns = @(
    '\bmklink\b',
    '-ItemType\s+["'']?(Junction|SymbolicLink|HardLink)',
    '\bln\s+-s\b'
)
foreach ($p in $linkPatterns) {
    if ($cmd -imatch $p) {
        Deny "Blocked: this command creates a junction/symlink/hardlink. Links inside this repo are forbidden -- a junction pointing at data/ was followed by 'git worktree remove' and destroyed the dataset. If a worktree needs the data, do not link it: run the data-touching script in the main checkout at $RepoRoot instead. See CLAUDE.md, 'Data safety'."
    }
}

# 2. Worktree teardown must go through the de-linking wrapper.
if ($cmd -imatch '\bgit\s+worktree\s+(remove|prune)\b') {
    Deny "Blocked: 'git worktree remove/prune' recurses into junctions on Windows and once deleted the contents of data/ through one. Use 'pwsh -NoProfile -File scripts/safe_worktree_remove.ps1 <worktree-path>' instead -- it unlinks any reparse point first, then removes the worktree. See CLAUDE.md, 'Data safety'."
}

# 3. Deletes that reach gitignored files. data/ and output/ are gitignored, so
# 'git clean -x' and 'git stash --all' remove them from the working tree even
# though no path in the command names them.
if ($cmd -imatch '\bgit\s+clean\b([^;|&]*)') {
    $flags = $Matches[1]
    $dryRun = ($flags -cmatch '(^|\s)-[a-zA-Z]*n[a-zA-Z]*(\s|$)') -or ($flags -imatch '--dry-run\b')
    $ignoredToo = $flags -cmatch '(^|\s)-[a-zA-Z]*[xX][a-zA-Z]*(\s|$)'
    if ($ignoredToo -and -not $dryRun) {
        Deny "Blocked: 'git clean -x' removes IGNORED files, and data/ and output/ are gitignored -- this deletes the entire dataset, with no stash and no recycle bin. If you only meant to drop untracked build artefacts, use 'git clean -fd' (without -x). To see what it would touch, 'git clean -xdn' is a dry run and is allowed. See CLAUDE.md, 'Data safety'."
    }
}
if ($cmd -imatch '\bgit\s+stash\b([^;|&]*)') {
    $flags = $Matches[1]
    if (($flags -imatch '(^|\s)--all\b') -or ($flags -imatch '(^|\s)-[a-zA-Z]*a[a-zA-Z]*(\s|$)')) {
        Deny "Blocked: 'git stash --all' stashes IGNORED files too, which strips data/ and output/ off disk. Restoring them depends on the stash surviving intact -- do not risk the dataset on that. Use 'git stash -u' (untracked, not ignored) if you need a clean tree. See CLAUDE.md, 'Data safety'."
    }
}

# 4. Recursive deletes aimed at anything called data.
# 'data' as a path segment only: matches data, ./data, data/raw, "data", data\raw --
# but not guard-data.ps1, metadata, or data-vault.
if ($cmd -imatch '(^|[\s"''/\\;|&(=])data([\s"''/\\;|&)]|$)') {
    $recursiveDelete = @(
        '\brm\s+-\S*[rR]',
        'Remove-Item\b[^;|]*-Recurse',
        '\brmdir\s+/[sS]\b',
        '\brd\s+/[sS]\b',
        '\bdel\s+/[sS]\b'
    )
    foreach ($p in $recursiveDelete) {
        if ($cmd -imatch $p) {
            Deny "Blocked: recursive delete naming 'data'. data/ is gitignored, expensive to re-download, and was lost once already. If you genuinely need to remove something under data/, the user must do it by hand after running 'scripts/protect_data.ps1 -Unlock'. See CLAUDE.md, 'Data safety'."
        }
    }
}

exit 0
