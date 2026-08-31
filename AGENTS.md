# AGENTS.md

## Project Description

This is a satellite precipitation correction project written in Julia, mainly involving:

- GWR / MGWR
- Geographically weighted regression
- FY4B, GPM, GSMaP satellite precipitation
- Rain gauge observation data

## Data safety (read before touching worktrees)

`data/` was destroyed once. A git worktree under `.claude/worktrees/` had no `data/` of its own
(it is gitignored), so a **junction** was created inside the worktree pointing at the main
checkout's `data/`. When `git worktree remove` later cleaned that worktree up, its recursive
delete followed the junction and emptied the real `data/`. This was re-confirmed on
git 2.52.0.windows.1 in a controlled drill: `git worktree remove --force` deleted every file
inside the junction's target and left the empty directory behind.

The rules that follow are not superstition — that is the actual failure.

- **Never create a junction, symlink, or hardlink anywhere in this repo or its worktrees**, above
  all not one pointing into `data/`. `mklink`, `New-Item -ItemType Junction|SymbolicLink|HardLink`
  and `ln -s` are blocked by a PreToolUse hook (`.claude/hooks/guard-data.ps1`).
- **Worktrees are for code-only edits.** Anything that reads or writes `data/` runs in the main
  checkout at `C:\Users\qw123\GG\GeoWeightedRegression.jl`. That is what removes the reason the
  junction existed in the first place. If a worktree seems to need the data, that is the signal to
  move the work back to the main checkout, not to link anything.
- **Never run `git worktree remove` or `git worktree prune` directly.** Use
  `pwsh -NoProfile -File scripts/safe_worktree_remove.ps1 <worktree-path>`, which unlinks every
  reparse point in the tree (with `cmd /c rmdir`, which removes the link and never the target),
  verifies none remain, and only then calls git. The bare commands are blocked by both a deny rule
  and the hook.
- **Never add `data` to `worktree.symlinkDirectories`** in any settings file — that setting
  recreates the exact hazard, automatically, for every new worktree.
- If you ever do need to remove a link by hand, use `cmd /c rmdir "<link>"`. Never
  `Remove-Item -Recurse` on a junction.
- **Never run `git clean -x` or `git stash --all`.** Neither names `data/`, but `data/` and
  `output/` are gitignored and both commands reach ignored files — `git clean -xfd` deletes the
  whole dataset outright, and `git stash --all` strips it off disk and makes recovery depend on a
  stash surviving intact. Use `git clean -fd` (no `-x`) or `git stash -u` instead; `git clean -xdn`
  is a dry run and is allowed. Both are blocked by the hook.

### The delete-lock, and its one hard limit

`scripts/protect_data.ps1` puts a deny-DELETE ACL on the irreplaceable parts of `data/`
(`data/raw/`, `data/FY4B/`, the top-level source files). Writing and overwriting stay allowed; only
deletion is blocked. `data/processed/` is left deletable on purpose — it is regenerated from raw,
and the `prepare_*` pipelines rewrite it via write-temp → `mv(force)` → `rm(temp)`.

```sh
pwsh -NoProfile -File scripts/protect_data.ps1 -Status
pwsh -NoProfile -File scripts/protect_data.ps1 -Unlock   # before a data re-download
pwsh -NoProfile -File scripts/protect_data.ps1 -Lock     # immediately after it finishes
```

**The lock does nothing for an elevated session.** A process running as Administrator holds
`SeBackupPrivilege`/`SeRestorePrivilege`, which bypass the DACL: measured on this machine, a delete
succeeded against a deny ACE for Everyone *and* Administrators. So: **start Claude Code and the
terminal from a normal window, never "Run as administrator"** — that is what makes this layer real.
Every session start warns when the session is elevated (the `SessionStart` tripwire), and
`-Status` prints it too.

Residual risk, accepted deliberately: there is no second copy of `data/`. A recursive delete run
by an elevated process outside a Claude Code session is not defended against by anything here, and
re-downloading is the recovery path.

A "permission denied" on a delete under `data/` is the guard working. Do not route around it with
`takeown`, `icacls`, or `-Force`; run `-Unlock`, do the thing deliberately, then `-Lock` again.

### What each layer does and does not cover

| Layer | Covers | Does not cover |
|---|---|---|
| `.claude/hooks/guard-data.ps1` (PreToolUse) | link creation, `git worktree remove/prune`, `git clean -x`, `git stash --all`, recursive deletes naming `data`, and `ExitWorktree` while a link is present | commands issued outside this project's Claude Code sessions; a link created from *inside* a script, since the hook reads the command text |
| `.claude/settings.json` deny rules | the bare worktree commands, and `Edit`/`Write` into `data/` | anything phrased differently — the hook is the real check |
| `scripts/protect_data.ps1` | any deleter, including git itself | elevated sessions (see above) |
| Session-exit "keep or remove worktree?" prompt | — | not a tool call, so no hook sees it; answer **keep** if the worktree ever held a link |

## Coding Requirements

- Use Julia.
- Avoid irrelevant refactoring.
- Do not change existing function names or parameters unless necessary for the current task.

## Modification Principles

- Only modify code relevant to the current task.
- If you find other problems, point them out first, do not modify them directly.
- Do not delete existing functionality.
- Prioritize simple and easy-to-understand implementation methods.

## Code Style

- Strive for conciseness and clarity; refer to existing code.

## Project Directory and File Conventions

When creating, modifying, or moving files, the following directory structure must be followed. Do not create temporary scripts, data files, images, or output results arbitrarily in the project root directory.

### `assets/` Stores static resources for the project

- Images and illustrations
- Project diagrams

### `data/` Stores research data

data/
├── raw/
└── processed/

Where:

- `data/raw/`: Raw data, modification is generally prohibited
- `data/processed/`: Data that has been cleaned, transformed, matched, or preprocessed

### `output/` Stores all results generated by the program; writing the results back to `data` is prohibited

## `scripts/` Stores executable task scripts

- `scripts/` is responsible for "calling" the core code
- Core algorithms should not be written directly in the scripts

### `scr/` Stores the core Julia source code of the project

For example:

- Distance calculation
- Spatial weights
- Bandwidth selection
- Weighted regression solution
- GWR / MGWR
- Accuracy evaluation tools

Rules:

- Reusable core functions must be placed in `src` first
- Do not write one-off experimental workflows in src/.
- Avoid reading fixed local absolute paths in `src/`
- Core functions should receive data and configuration via parameters whenever possible

## New File Placement Rules

When creating a new file, first determine its purpose

- Program output → `output/`
- Core algorithm → `src/`
- Directly runnable experiments or processing procedures → `scripts/`

`scripts/` handles the process flow, `src/` handles reusable core algorithms