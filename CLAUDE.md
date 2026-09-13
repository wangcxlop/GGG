# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Description

This is a satellite precipitation correction project written in Julia, mainly involving:

- GWR / MGWR (Geographically Weighted Regression / Mixed GWR)
- FY4B, GPM, GSMaP satellite precipitation
- Rain gauge observation data

The Julia package itself is named `MixedGWR` (see `Project.toml`) — `using MixedGWR` is the correct import, not the repo directory name.

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

## Commands

```sh
# Install/resolve dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the full test suite (test/runtests.jl)
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file directly (most test files are self-contained:
# they `include` the src file they test and can run standalone)
julia --project=. test/test-study-area.jl
julia --project=. test/test-mger-five-kernels.jl
```

`test/test-speed.jl` is a standalone benchmark, not wired into `test/runtests.jl` — run it directly when checking performance.

The first `@testset "GWR"` block in `test/runtests.jl`, plus `test-solver.jl`, `test-ST_GWR.jl` and `test-GWR_mixed.jl`, validate results against the R `GWmodel` package via `RCall` using the `data/prcp_st174_shiyan.csv` fixture. `data/` is gitignored, so that fixture is absent from a fresh checkout and from CI; those four are skipped with a warning when it is missing (`HAS_SHIYAN_DATA` in `test/main_pkgs.jl`) and the other 14 test files still run. Restoring the file re-enables them. Running them needs a working R installation with `GWmodel` installed (`test = ["Test", "Distances", "RCall", "RTableTools"]` in `Project.toml`). Individual test files that don't touch `RCall`/`RTableTools` (e.g. `test-study-area.jl`, `test-era5-*.jl`) can be run without R.

Scripts under `scripts/run_*.jl` and `scripts/verify_*.jl` are the actual entry points for producing results (e.g. `scripts/run_mger_smoke_202206.jl`, `scripts/run_interpolation_benchmark.jl`). They set `LOAD_PATH` to `src/` and are run as plain Julia scripts:

```sh
julia --project=. scripts/run_mger_smoke_202206.jl
```

Run the interpolation benchmark multithreaded, but **not** with `-t auto`. Its dominant loop is
`dynamic_covariate_predict`, which fits one hour per task, and each of those tasks allocates
~1.6 MB — two n×n local-hat matrices among it. Past a handful of threads the large-object
allocation contends worse than the extra cores help, so the loop gets *slower*. Measured on this
24-core box, one full-mode call over 334 tuning hours:

| Julia threads | 4 | 8 | 24 |
|---|---|---|---|
| loop wall time | **167 ms** | 201 ms | 274 ms |

`-t 4` also beat `-t auto` end to end on the smoke benchmark (314 s vs 393 s), with byte-identical
output — the results do not depend on the thread count, only the wall clock does.

```sh
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates
```

Four is not a magic number: it is where this machine's allocation contention starts to bite.
Re-derive it on new hardware with `scripts/profile_hour_fit.jl`, which prints the per-hour cost
and the loop's speedup at the thread count it is given. Do not run single-threaded — that is
~3x slower than the best setting.

Do **not** set `BLAS.set_num_threads`. It looks like free speed — the GWR hot path is gemv and
p×p solves, where BLAS threading only contends with the Julia-level threading above — but `tps`
factorizes a dense (n+3)×(n+3) system that is well past OpenBLAS's threading threshold, and a
threaded LU accumulates in a different order. Pinning BLAS to one thread was measured to move
`tps`'s RMSE in the last three or four digits, which propagates into `metrics_pooled.csv`,
`paired_comparisons.csv` and `claim_assessment.csv`. Every other method was byte-identical. The
run's BLAS thread count is therefore part of what makes published numbers reproducible.

The same class of drift showed up again, independent of BLAS, when the full-mode invariance gate
(`scripts/verify_perf_invariance.jl`) was run against the hour-fit perf pass (commit `b822d78`):
`paired_comparisons.csv`'s `ci_high` column moved in its last one or two digits for most rows,
while `ci_low`, `delta_RMSE`, `relative_improvement`, and every other output file (105 of 106)
stayed byte-identical. `_daily_bootstrap_delta`/`paired_bootstrap_rows`
(`src/benchmark/InterpolationBenchmarkBootstrap.jl`) are themselves single-threaded and fully seeded, so the
difference traces to a sub-ULP perturbation in the `residual_gwr` prediction matrix upstream —
the perf pass's threaded prediction path reordering some sum. It is only visible in
`paired_comparisons.csv` because that file's RMSE is summed per-day (≤526 terms) before
bootstrapping, where a 1-ULP shift is a meaningfully large fraction of the sum; everywhere else
RMSE is pooled over the full ~2.6M-cell dataset, where the same shift is far below print
precision and rounds away. Treat a last-digit-only `ci_high` (never `ci_low`, never the
underlying RMSE/delta columns) as this same benign non-associativity, not a regression.

## Architecture

### The layout of `src/`

`src/` holds two files and seven directories. The directory a file is in says which of the three
loading disciplines applies to it — that used to be prose, and the filename could not carry it:
seven core-tier files are lowercase and four PascalCase, while the top-level fragments look
exactly like standalone modules.

```
src/
  MixedGWR.jl      package entry point (path fixed by Project.toml)
  load_modules.jl  the loader every script and test goes through
  core/            the MixedGWR module's own algorithm files
  sources/         one module per data source
  selection/       the variable-selection paths
  models/          the estimators
  benchmark/       the interpolation benchmark, and the diagnostics that read its output
  mger/            the MGER pipeline and its data prep
  util/            plumbing with no opinion about precipitation, shared across tiers
```

`util/` is deliberately narrow: it holds `TableIO` (the atomic CSV write and the case-insensitive
column lookup) because those callers span `sources/` and `mger/` and the helpers belong to neither.
It is not a drawer for anything that happens to be shared - domain vocabulary goes with its domain,
which is why `CovariateGroups` sits in `models/` rather than here.

Do not create `src/data/` or `src/docs/`: `.gitignore` carries bare `data/` and `docs` patterns,
which git matches at any depth, so files there would be silently untracked.

1. **`src/core/` — inside the `MixedGWR` module.** `fitted.jl`, `metrics.jl`, `kernel.jl`,
   `gw_weight.jl`, `PrecipitationCorrection.jl`, `solve_chol.jl`, `solve_reg.jl`, `GWR.jl`,
   `GWR_calib.jl`, `ST_GWR.jl`, `deprecated.jl`, each `include`d by `src/MixedGWR.jl` inside
   `module MixedGWR ... end`. Reached normally via `using MixedGWR`; they export the
   regression/kernel primitives (`GWR`, `ST_GWR`, `ST_GWR_fast`, kernel constants
   `GAUSSIAN`/`EXPONENTIAL`/`BISQUARE`/`TRICUBE`/`BOXCAR`, etc). These files define no module of
   their own, so they must never be `include`d directly by anything else.

2. **Standalone modules** — each defines its *own* `module X ... end` and is loaded through
   `src/load_modules.jl`, never through the `MixedGWR` module. Twenty of them:
   - `sources/`: `StudyArea.jl`, `FY4BPreprocessing.jl`, `ERA5LandStations.jl`,
     `ERA5LandProcessing.jl`, `ERA5LandCovariates.jl`, `MOD13A2NDVIProcessing.jl`,
     `AppEEARSNDVI.jl`, `TerrainFeatures.jl`, `LandformClassification.jl` (DEM relief regions for
     the gauges, window chosen by mean change-point) — one data source or ingest stage each.
   - `selection/`: `SelectionScaffolding.jl` (bookkeeping shared by the searches —
     `annotate_selection!`, `append_selection!`, `selection_schemes`), `ERA5VariableSelection.jl`,
     `NDVIVariableSelection.jl`, `JointVariableSelection.jl`.
   - `models/`: `TraditionalInterpolation.jl`, `DEMTerrainExperiment.jl`,
     `JointCovariateModels.jl`, `SatelliteFusion.jl` (least-squares or equal-weight merge of
     several products into one anchor) — everything that fits something — plus
     `CovariateGroups.jl`, the group/column/family vocabulary the fitting and the selection sides
     share.
   - `util/TableIO.jl`: `write_csv_atomic` and the case-insensitive `column` lookup, shared by
     seven modules across `sources/` and `mger/`.
   - `benchmark/BenchmarkDiagnostics.jl`, which reads a finished run's artefacts. It shares that
     directory with the fragments below but not their discipline: it is a real module and takes no
     part in their include chain.
   - `benchmark/SatelliteTemporalEvaluation.jl`: hourly, event, diurnal and regional-series scores
     of the satellite products against the gauges, stratified by intensity class, season and
     landform region. Same standalone discipline as `BenchmarkDiagnostics`.
   - `mger/MGERDataPrep.jl`, beside the pipeline it serves.

3. **`MGERPipeline.jl` and `InterpolationBenchmark.jl` are *not* modules** — they are top-level
   fragments (`using MixedGWR` + struct/function definitions) evaluated straight into `Main` by
   `load_pipeline`, after `using MixedGWR` is already active. They tie the core algorithms and the
   standalone modules together into full run/evaluate pipelines (e.g. `MGERConfig`,
   `run_multikernel_spatial_kfold_pipeline`).

`InterpolationBenchmark.jl` is a thin loader: it pulls in the modules the benchmark needs and then
includes ten concern-specific fragments, in this order — `Config`, `Folds`, `DEM`, `Joint`,
`Predictors`, `Hurdle`, `Tuning`, `Metrics`, `Bootstrap`, `Run`. They are plain top-level
fragments sharing one namespace, not modules, so a name defined in a later file may be called from
an earlier one; include order only has to put shared consts and structs first.
`InterpolationBenchmarkHurdle.jl` holds the deliberately-disabled `hurdle_gwr` model, which is
absent from `BENCHMARK_RUNS` and unreachable in a normal run.

Three opt-in flags widen what a run reports without moving any pre-existing number (the end-to-end
test asserts every incumbent metric row is identical with them on):

- `--satellite-wet-blend` adds `blend_residual_gwr`, `blend_mixed_gwr` and `blend_mgwr`: where the
  satellite reports rain, the anchored prediction is blended toward `adw`, with the weight chosen
  per fold on the inner selection split (`blend_selection.csv`).
- `--blend-axis agreement_envelope` adds `blend_agrenv_*` beside them, with one weight per band of
  (how many products call the cell wet) x (their maximum).
- `--fused-anchor` adds the derived products `MERGED_OLS` and `MERGED_MEAN`, whose anchor is a
  combination of FY4B, GPM and GSMaP. The OLS coefficients are fitted inside every fold on its
  training stations only (`fused_anchor_selection.csv`), so a held-out gauge never reaches its own
  anchor; it needs nested covariate selection. `run_claim_reassessment.jl` assesses every product
  the run scored, derived ones included.

Blended methods are scored on the shared mask but never join `MASK_METHODS`, so they cannot move
another method's denominator.

### Loading `src/` from a script or test

There is exactly one idiom. Never `include` a standalone module file directly: `include`ing a file
that defines `module X` a second time compiles a second, type-incompatible copy of it rather than
reusing the first.

```julia
include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")            # MGERPipeline / InterpolationBenchmark
load_standalone_modules("BenchmarkDiagnostics")    # any module; dependencies load first
using Main.BenchmarkDiagnostics
```

`load_standalone_modules` loads each module into `Main` at most once, pulling in that module's own
sibling dependencies first; `load_pipeline` does the same for the two top-level fragments, which
are detected by a sentinel struct rather than a module name. `load_modules.jl` itself is safe to
include more than once.

Names, not paths: `locate_src_file` walks `src/` for `<name>.jl`, so a caller never says which tier
a module lives in and moving a file between directories needs no edit anywhere. It refuses rather
than guesses in two cases - a name that matches nothing, and a name that matches twice. The second
is the one to know about: two files with the same basename in different directories (a leftover
copy from a half-finished move) would otherwise be resolved by walk order.

Do **not** add `pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))`. It makes `src/` an implicit
environment and causes `using MixedGWR` to load a second copy of the package, which is a hard
error on Julia 1.11.

Inside a standalone module, reach a sibling with `using Main.X` (add `using Main: X` as well if
the module name itself is used for qualified calls — `using Main.X: a, b` does not bind `X`).

When adding a new file to `src/`, the directory is the decision: a reusable regression/kernel
primitive goes in `core/` and needs an `include("core/...")` line in `src/MixedGWR.jl`; anything
else is its own standalone module under `sources/`, `selection/`, `models/` or `benchmark/`
following convention (2) above, and needs no registration at all — the loader finds it by name.

Two modules are themselves split, both by `include`ing fragments into the module rather than by
making sub-modules, so every name stays exactly where callers reach it today:
`benchmark/diagnostics/` holds one file per diagnostic, keeping the `D<n>` labels the sections
carried as banner comments, and `models/dem/` holds the terrain/GWR library with
`run_dem_experiment`, the study orchestrator that writes ~15 fixed-name CSVs, alone in
`experiment.jl`.

### Calling convention

`scripts/` calls into `src/` — algorithms and reusable logic live in `src/`, not in scripts. A script typically: sets `LOAD_PATH`, does `using MixedGWR`, `include`s any standalone modules/pipeline files it needs, then builds a config struct and calls a pipeline function.

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

```
data/
├── raw/
└── processed/
```

Where:

- `data/raw/`: Raw data, modification is generally prohibited
- `data/processed/`: Data that has been cleaned, transformed, matched, or preprocessed

### `output/` Stores all results generated by the program; writing the results back to `data` is prohibited

### `scripts/` Stores executable task scripts

- `scripts/` is responsible for "calling" the core code
- Core algorithms should not be written directly in the scripts

### `src/` Stores the core Julia source code of the project

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
