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

Every option that changes what is fitted also changes the output directory name, so a run can
never overwrite the baseline it is measured against. Three recent ones:
`--satellite-wet-blend` (`_satwetblend`) blends the anchored GWR-family predictions toward `adw`
where the satellite reports rain — the only variant that beats `adw`, see F4;
`--free-satellite-coefficient` (`_freesat`) fits the satellite's coefficient locally instead of
forcing it to 1, which F4 records as measured and negative; and `--equal-grids` (`_equalgrids`)
widens the GWR *and* the IDW/ADW/TPS search grids together, unlike `--local-grid`, which widened
only the GWR family and so left the baselines pinned against their own ceilings.

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
(`src/InterpolationBenchmarkBootstrap.jl`) are themselves single-threaded and fully seeded, so the
difference traces to a sub-ULP perturbation in the `residual_gwr` prediction matrix upstream —
the perf pass's threaded prediction path reordering some sum. It is only visible in
`paired_comparisons.csv` because that file's RMSE is summed per-day (≤526 terms) before
bootstrapping, where a 1-ULP shift is a meaningfully large fraction of the sum; everywhere else
RMSE is pooled over the full ~2.6M-cell dataset, where the same shift is far below print
precision and rounds away. Treat a last-digit-only `ci_high` (never `ci_low`, never the
underlying RMSE/delta columns) as this same benign non-associativity, not a regression.

## Architecture

### Two tiers of `src/`

1. **Core `MixedGWR` module** — algorithms included inside `module MixedGWR ... end` in `src/MixedGWR.jl` (the package entry point): `fitted.jl`, `metrics.jl`, `kernel.jl`, `gw_weight.jl`, `PrecipitationCorrection.jl`, `solve_chol.jl`, `solve_reg.jl`, `GWR.jl`, `GWR_calib.jl`, `ST_GWR.jl`, `deprecated.jl`. These are reached normally via `using MixedGWR` and export the regression/kernel primitives (`GWR`, `ST_GWR`, `ST_GWR_fast`, kernel constants `GAUSSIAN`/`EXPONENTIAL`/`BISQUARE`/`TRICUBE`/`BOXCAR`, etc).

2. **Standalone data-pipeline modules** — each of these files defines its *own* `module X ... end` and is loaded through `src/load_modules.jl`, not through the `MixedGWR` module: `StudyArea.jl`, `ERA5LandStations.jl`, `ERA5LandProcessing.jl`, `ERA5LandCovariates.jl`, `ERA5VariableSelection.jl`, `MOD13A2NDVIProcessing.jl`, `NDVIVariableSelection.jl`, `AppEEARSNDVI.jl`, `FY4BPreprocessing.jl`, `TerrainFeatures.jl`, `TraditionalInterpolation.jl`, `DEMTerrainExperiment.jl`, `JointCovariateModels.jl`, `JointVariableSelection.jl`, `MGERDataPrep.jl`, `BenchmarkDiagnostics.jl`. Each handles one data source or processing stage (ERA5-Land, MOD13A2 NDVI, FY4B, terrain/DEM, variable selection, benchmark diagnostics, etc).

3. **`SelectionScaffolding.jl`** — bookkeeping shared by the four variable-selection paths
   (`annotate_selection!`, `append_selection!`, `selection_schemes`). A standalone module like
   those in (2), loaded the same way.

4. `MGERPipeline.jl` and `InterpolationBenchmark.jl` are *not* modules — they are top-level scripts (`using MixedGWR` + struct/function definitions) meant to be `include`d directly by a script or test after `using MixedGWR` is already active. They tie the core GWR algorithms and the data-pipeline modules together into full run/evaluate pipelines (e.g. `MGERConfig`, `run_multikernel_spatial_kfold_pipeline`).

`InterpolationBenchmark.jl` is a thin loader: it pulls in the modules the benchmark needs and then
includes ten concern-specific fragments, in this order — `Config`, `Folds`, `DEM`, `Joint`,
`Predictors`, `Hurdle`, `Tuning`, `Metrics`, `Bootstrap`, `Run`. They are plain top-level
fragments sharing one namespace, not modules, so a name defined in a later file may be called from
an earlier one; include order only has to put shared consts and structs first.
`InterpolationBenchmarkHurdle.jl` holds the deliberately-disabled `hurdle_gwr` model, which is
absent from `BENCHMARK_RUNS` and unreachable in a normal run.

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

Do **not** add `pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))`. It makes `src/` an implicit
environment and causes `using MixedGWR` to load a second copy of the package, which is a hard
error on Julia 1.11.

Inside a standalone module, reach a sibling with `using Main.X` (add `using Main: X` as well if
the module name itself is used for qualified calls — `using Main.X: a, b` does not bind `X`).

When adding a new file to `src/`, follow the existing pattern: if it's a reusable regression/kernel primitive it belongs inside the `MixedGWR` module (add an `include(...)` line in `src/MixedGWR.jl`); if it's a data-source-specific processing step it should be its own standalone module following convention (2) above.

### Calling convention

`scripts/` calls into `src/` — algorithms and reusable logic live in `src/`, not in scripts. A script typically: sets `LOAD_PATH`, does `using MixedGWR`, `include`s any standalone modules/pipeline files it needs, then builds a config struct and calls a pipeline function.

## GWR-family design review — findings not yet acted on

A design review of the GWR-family methods (2026-08-31) produced nine findings. Five were fixed on
branch `fix/local-hat-unsupported-target-nan`, each commit carrying its own analysis: `2996b07`
F1 (unsupported local targets), `8ff27c6` F2 (duplicate methods), `f7967d8` E1+E2 (evaluation
mask), `81b8602` F3 (direct vs residual gwr), `6ee1497` F5 (mgwr nesting), `fb92dd8` F6 (role test
vs bandwidth search).

The four below were **not** acted on. They are recorded here rather than in an issue tracker
because each is a statement about how the benchmark is put together, and the evidence for it is
worth more than the one-line summary.

### F4 — residual mode is a property of the family, not a factor

`BENCHMARK_METHODS` (`src/InterpolationBenchmarkConfig.jl:1`) offers both `gwr` and
`residual_gwr`, but `TRADITIONAL_METHODS = ["idw", "adw", "tps"]` (line 6) has no residual
counterpart — there is no `residual_idw` or `residual_tps`. Every method that corrects satellite
residuals is therefore also a GWR, and every non-GWR method interpolates the gauge field directly.
"`residual_gwr` beats `tps`" consequently varies the framing and the estimator at the same time,
and the design cannot separate them. Making residual mode an orthogonal factor would need a
residual variant of at least one traditional method.

**Re-measured on the canonical baseline, 2026-09-02.** The framing is not just unseparated from
the estimator — it is where the whole loss comes from. `satellite_quadrant.csv` splits the MSE gap
against `adw` on the joint wet/dry state of gauge and satellite (balanced_spatial, `mgwr`; shares
are FY4B's, the other two differ by a point or two):

| quadrant | share | `mgwr` MSE FY4B / GPM / GSMaP | `adw` MSE | gap contribution |
|---|---|---|---|---|
| gauge dry, satellite dry | 87.4% | 0.047 / 0.015 / 0.039 | 0.049 / 0.039 / 0.040 | **−0.001 / −0.021 / −0.001** |
| gauge dry, satellite wet | 5.7% | 3.424 / 1.895 / 1.867 | 0.631 / 0.680 / 0.647 | **+0.158 / +0.079 / +0.081** |
| gauge wet, satellite dry | 4.6% | 5.567 / 5.670 / 4.333 | 5.340 / 5.162 / 4.191 | +0.010 / +0.012 / +0.003 |
| gauge wet, satellite wet | 2.4% | 18.220 / 12.468 / 14.376 | 18.296 / 12.187 / 12.665 | −0.002 / +0.013 / **+0.078** |

The core finding survives the rebuild and the IDW/ADW correction. The family still *beats* `adw`
in the dominant dry/dry quadrant, and the dry-gauge/wet-satellite quadrant is 5.7% of cells but
**95.7% of the gap** for FY4B and 94.7% for GPM. A false alarm is patchy at the satellite's own
error scale, so the training stations' residuals carry no information about it and a spatially
smooth correction cannot cancel a value the model was handed.

One thing did change: **GSMaP is no longer a one-quadrant story.** Its gap is now split roughly
evenly between dry/wet (50.2%) and wet/wet (48.5%), where it previously lost almost nothing on
heavy rain. Anything inferred from GSMaP alone needs re-checking against that.

Three things follow, the first two unchanged in direction and the third **reversed**:

- Gating the *correction* still cannot help — it only reaches the satellite-dry quadrants, where
  nothing is wrong. Measured on `mgwr`: the best dry-cell gate buys 0.1–1.0% RMSE and costs
  2.3–6.1% relative POD.
- The loss is reachable, but only by acting on satellite-wet cells. `satellite_wet_blend->adw` —
  blend toward `adw` where the satellite reports rain — is the **only** counterfactual in
  `anchor_discount_bounds.csv` that beats `adw` at all, and it lifts RMSE and POD together:
  `mgwr` goes **+0.27% / +2.98% / +2.21%** vs `adw` with POD rising 0.802→0.822, 0.787→0.811,
  0.792→0.846. It beats even the `oracle_dry_wet` variants, which are allowed to see the true
  gauge state — so the win is not about detecting false alarms, it is about not trusting the
  satellite's magnitude where it claims rain.
- **A plain anchor discount no longer works.** On the legacy run it lifted RMSE and POD together;
  on this baseline it loses to `adw` on all three products (−2.03% / −0.13% / −2.16%) *and* costs
  5–9 points of POD (0.802→0.759, 0.787→0.739, 0.792→0.703). Do not cite the old
  "0.910 against `adw`'s 0.944 with POD rising 0.813 → 0.831" result; it was measured against
  weaker `idw`/`adw` baselines, before those formed their weights per availability group.

**Both levers have now been run, and the blend is the one that works** (2026-09-03).

`--satellite-wet-blend` (suffix `_satwetblend`, commit `8f36c43`) turns the winning counterfactual
into a tuned method: it blends each anchored GWR-family prediction toward `adw` where the
satellite reports rain, with the blending weight chosen on the fold's inner selection split rather
than assumed. Scored on `balanced_spatial`, it is the first GWR-family method to beat `adw`:

| product | `adw` RMSE | `blend_mgwr` RMSE | vs `adw` | replay ceiling | tuning cost |
|---|---|---|---|---|---|
| FY4B | 0.874121 | 0.875248 | −0.13% | +0.27% | 0.40 pts |
| GPM | 0.870897 | 0.845295 | **+2.94%** | +2.98% | 0.04 pts |
| GSMaP | 0.870897 | 0.853139 | **+2.04%** | +2.21% | 0.17 pts |

GPM and GSMaP are significant under the paired daily bootstrap against **all three** traditional
baselines (2000 reps, 618 days, `ci_low > 0` and `pvalue_holm < 0.05`); FY4B is not (p = 0.752
raw against `adw`, 1.0 after Holm), and never could be — a +0.27% replay ceiling does not survive honest tuning. The
"tuning cost" column is the gap between the ceiling the replay promised and what selecting the
weight out of sample actually delivered: 0.04 and 0.17 points on GPM and GSMaP, 0.40 on FY4B.

**The pre-registered claim still fails, and on one gate only.** `product_supported` is `false` for
all three products because `assess_gwr_claim` (`src/InterpolationBenchmarkMetrics.jl:296`) demands
a **≥5%** RMSE improvement on the *heavy* stratum against every baseline. The blend delivers
+0.43% / +0.81% / +2.22% there — positive, and significant against `idw`/`adw` at p = 0.0, but an
order of magnitude short of that bar, so `heavy_win_count = 0` everywhere. Every other gate passes
on GPM and GSMaP: significance 3/3, moderate non-inferiority (it *gains* 2.4% / 1.7% rather than
degrading), 3 of 4 years, CSI/FAR not degraded, own coverage 0.9805. So "the blend beats every
traditional baseline on GPM and GSMaP" is defensible; "the GWR claim is supported" is not.

Reproduce with `scripts/run_claim_reassessment.jl <run dir>`, which is also the only place the
blended methods are assessed at all — `claim_assessment.csv` is written for `DEFAULT_CLAIM_METHOD`
alone. Heed its own warning: reporting the best of the eight assessed methods is a post-hoc
maximum over correlated tests, and that bites here. On GSMaP only `blend_mgwr` clears
significance; `blend_residual_gwr` (+0.08%) and `blend_mixed_gwr` (+0.10%) do not.

`--free-satellite-coefficient` acts on the anchor instead: it puts the satellite into the local
design as `JointCovariateModels.SATELLITE_GROUP`, so the effective coefficient becomes
`1 + b_sat(u)` rather than a forced 1. Off by default, suffix `_freesat`. It is **answered and
negative**, and it broke on the way. Where it completed (FY4B) it improved the family a great deal
and still lost: −3.25% against `adw`. On GPM it **crashed** — every kernel/bandwidth-family
combination failed to converge within `mgwr_max_tuning_iterations = 5`, `converged || continue`
dropped all of them, and the run died with `no common valid OOF samples across all methods` on 4
of 5 folds. That is F8 firing for real, on precisely the harder configuration that section names,
since the flag adds one covariate group. Its run directory and log therefore cover FY4B only, and
its 91 non-convergence warnings are a truncated count, not a census.

The blend avoids all of this because it never refits anything: it combines two already-fitted
predictions, adds no design column, and adds no back-fit failures beyond the pre-existing floor
recorded in F7 below. That robustness difference is independent of the RMSE result.

#### The two gates are opposed along the blending weight, and no banding fixes it

Measured 2026-09-04 with `scripts/verify_banded_blend_bounds.jl`, and this is the thing to know
before proposing another blend variant.

The family's error is regime-split. On `balanced_spatial`/FY4B the unblended `mgwr` is 73% worse
than `adw` on the 93% of cells where the gauge is dry, and **4.2% better than `tps`** on the 0.25%
where it is heavy — the satellite recovers part of the heavy-rain underestimation gauge-only
interpolation smooths away (heavy bias −10.63 against `adw`'s −11.53). Unblended `mgwr` on GSMaP
already clears the ≥5% heavy gate against all three baselines (5.78%) and fails only the overall
gate; the blend passes overall and fails heavy. It is tempting to read that as two configurations
that just need combining.

They cannot be combined. Pooled RMSE and the heavy stratum move in **opposite directions,
monotonically**, along the blending weight — `mgwr`, `balanced_spatial`, heavy quoted against
`tps`, which is the baseline that binds:

| λ | 0.0 | 0.3 | 0.5 | 0.7 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|
| FY4B pooled RMSE | 0.96391 | 0.90781 | 0.88421 | 0.87275 | **0.87176** | 0.87926 |
| FY4B heavy vs `tps` | **+4.23%** | +2.77% | +1.56% | +0.15% | −0.62% | −2.28% |
| GSMaP pooled RMSE | 0.95902 | 0.88824 | 0.86101 | **0.85165** | 0.85391 | 0.87211 |
| GSMaP heavy vs `tps` | **+5.78%** | +4.96% | +3.68% | +1.85% | +0.73% | −1.88% |

Pooled is minimised at λ ≈ 0.7–0.8; heavy is maximised at λ = 0 and falls from there. GSMaP needs
λ ≤ 0.3 to hold heavy ≥ 5%, and at λ = 0.3 pooled is 0.888 against `adw`'s 0.871. FY4B and GPM
never reach 5% against `tps` at **any** λ, λ = 0 included.

**Giving the weight a band per satellite-intensity class does not break the trade-off**, and the
reason generalises past this particular fix. Bands `[0.1, 2.5) / [2.5, 8) / [8, ∞)` on `y_sat`,
each band's λ minimised independently (exact — the bands partition the cells and blending never
changes which cells are finite, so pooled SSE is additive; the script audits that against a joint
enumeration rather than asserting it). The result is a ~0.1-point pooled gain and no heavy
recovery at all: FY4B `mgwr` +0.35% pooled with heavy −0.30%, GPM +3.10% / +0.67%, GSMaP +2.39% /
+1.52%. A second axis, `nearest_train_km`, does the same. Every band picks λ ≈ 0.6–0.9.

The mechanism is that **the satellite does not know the rain is heavy**, so conditioning on its
value cannot find the cells worth protecting:

| product | median `y_sat` on gauge-heavy cells | gauge-heavy cells the satellite calls dry | gauge-heavy share of the `y_sat ≥ 8` band |
|---|---|---|---|
| FY4B | 0.85 mm/h | 36.7% | 7.4% |
| GPM | 2.35 mm/h | 16.8% | 19.4% |
| GSMaP | 1.89 mm/h | 14.5% | 12.9% |

The median gauge value on those cells is 11.5 mm/h. Over a third of FY4B's gauge-heavy cells sit
below the wet threshold, where the blend never reaches them at all, and the satellite-heavy band
is 81–93% *not* gauge-heavy, so that band's own error is still minimised by blending hard. Any
rule keyed on `y_sat` inherits this; the useful conditioning variable would have to be one the
satellite's magnitude does not already fail at.

So the blend as shipped is at its ceiling, and "beats `adw` on GPM and GSMaP" and "supports the
pre-registered claim" are not two steps along one path. The heavy gate needs a method that
improves heavy rain *without* trading it against the dry cells — not a better weight.

### F7 — two back-fits, two stopping rules, one nominal tolerance

`DEMTerrainExperiment._backfit_components` (`src/DEMTerrainExperiment.jl:597` and `:708`) stops on
the relative change in **RSS**. `JointCovariateModels._multiscale_predict_damped`
(`src/JointCovariateModels.jl:607`) stops on the relative **L2 change in the fitted vector**. Both
default to `tolerance = 1e-5`, so the same configured number means two different things depending
on which path a model takes. They also disagree on failure: the joint one returns all-NaN with
`converged=false`, discarding the hour; the DEM one returns its last iterate.

**The joint one fires on the canonical baseline**, which was not previously written down here
(measured 2026-09-03). It discards ~228 of the 13471 hours (~1.7%) for FY4B, logging
`joint dynamic model did not converge for some hours` from
`src/InterpolationBenchmarkPredictors.jl:140` — 187 warnings in the baseline run, 168 of them at
exactly `failed_hours = 228`, spread over `residual_gwr`, `mixed_gwr` and `mgwr` alike. GPM (6)
and GSMaP (1) are almost untouched. Those hours are already excluded from the shared evaluation
mask, which is why FY4B's mask holds 3,020,184 cells against GPM/GSMaP's 3,072,559, so no
published number is wrong because of it — but a cross-*product* comparison is not on the same
cells, and the discard is silent apart from the warning.

Census it with `grep -c "did not converge"` on a run's stderr log plus
`grep -o "failed_hours = [0-9]*" | sort | uniq -c`. The completed `_satwetblend` run reproduces
the baseline census **exactly** — 187 warnings, 180 FY4B / 6 GPM / 1 GSMaP, identical triple by
triple on (product, method, failed_hours) — which is the sharpest available evidence that blending
adds no fitting failures, since it refits nothing. `_freesat`, by contrast, pushed the count to
236–419 on a third of its warnings. Adding the satellite to the local design does make the
back-fit harder; the 228-hour core is not its doing.

### F8 — non-convergence removes a candidate instead of penalising it

`converged || continue` at `src/InterpolationBenchmarkTuning.jl:199` drops an entire
kernel/bandwidth-family combination from the mgwr comparison when its coordinate descent has not
stabilised within `mgwr_max_tuning_iterations`; `src/InterpolationBenchmarkJoint.jl:236` does the
same one level down. A combination that converges slowly is not scored worse — it is not scored at
all, so convergence speed can decide which kernel wins.

Not currently firing, and re-checked on the canonical baseline (2026-09-02). Of the 81 selected
mgwr descent rows in `parameter_scan.csv`, 70 converge at iteration 2, 7 at 3 and 4 at 4; the
deepest is 4 and **none reaches the 5-iteration cap**, so `converged || continue` discards nothing
and no published number depends on it. The legacy run was 74 / 3 / 4 over the same 81 rows — the
descent got marginally slower with the larger dataset without approaching the cap.

The margin is thinner than "none at the cap" suggests, though. The cap is 5 and the deepest
observed descent is 4, so a single extra sweep separates the current numbers from the regime where
combinations start being dropped silently. The finding is about fragility under a harder
configuration — more covariate groups, a wider bandwidth grid, `--equal-grids`.

**And it has since fired** (2026-09-03). `--free-satellite-coefficient` adds one covariate group,
the first item on that list, and on GPM it exhausted the cap for *every* kernel/bandwidth-family
combination: the tuner dropped them all and the run died with `no common valid OOF samples across
all methods` on 4 of 5 folds (see F4). The fragility is demonstrated, not hypothetical, and one
extra covariate group was enough. Do not confuse it with F7's per-hour discard — a different
mechanism at a different level, which fires on the canonical baseline where this one does not.

Reproduce with: selected `method == "mgwr"` rows of `parameter_scan.csv`, counted by `iteration`.

### F9 — minor items

- `any_converged .|= group_converged` (`src/InterpolationBenchmarkJoint.jl:76`) should be `.&=` if
  the flag is meant to say "every group converged". Inert while the value is only consumed as
  "any", but the variable name already claims the stronger reading.
- `_predict_time` (`src/JointCovariateModels.jl:540`) is dead — zero callers — and carries a
  *third* back-fit with its own stopping rule, so F7's divergence is latent in triplicate.
  `_mixed_predict_complete` and `_multiscale_predict_complete` in `src/DEMTerrainExperiment.jl`
  are dead as well.
- The ridge is an absolute constant applied to designs whose columns differ in scale by orders of
  magnitude, so it does not regularise them comparably.
- GPM and GSMaP emit bit-identical baseline rows, which is worth confirming is intended.
- The bandwidth grid still saturates at its endpoints in some fold-cells, and the partial
  `--equal-grids` run below says which methods this actually binds: widening the search left
  `idw` and `tps` **byte-identical**, so their selections were interior, while `adw` moved.

### Outstanding verification: F1's benchmark-level gate

F1 changed unfittable local targets from a fabricated all-zero hat row to NaN. Unit tests cover
it; the benchmark-level gate has never run.

**The protocol this section used to describe is now unreachable, and must not be followed.** It
asked for a `--legacy-unsupported-zero` re-run to come out *byte-identical* to the baseline
directory on disk. Three separate things now make that impossible, and none of them is F1: the
data was rebuilt (11426 → 13471 full-mode hours, and 8067 → 8069 even in the Jun–Sep window), the
IDW/ADW kernel was rewritten to form its weights per availability group, and the old baseline is
archived under `output/_legacy_11426h_20260829/`. A byte comparison against it would fail for
reasons that have nothing to do with unsupported local targets.

The achievable replacement is a **paired run on one code and one dataset**, where the flag is the
only thing that differs (~10 h total):

```sh
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
      --legacy-unsupported-zero
julia --project=. scripts/compare_benchmark_runs.jl \
      output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only \
      output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_legacyzero
```

No snapshot step is needed any more: the two runs land in different directories by construction,
because the suffix is keyed on the legacy flag (see below). Watch the shared evaluation mask. F1
is *expected* to shrink it — targets that were silently counted as fitted now drop out honestly —
so a mask that does not move is the surprising outcome, not a passing one.

### Output directory naming: the clean name is the corrected default

Every option that changes what is fitted adds a suffix, so a run can never overwrite the baseline
it is measured against. The unsupported-target suffix was keyed the other way round until
2026-09-02: the *corrected* default carried `_nanunsupported` and the pre-fix path took the clean
name, so that a corrected run could not overwrite the pre-fix baseline it was being compared with.
That baseline is now archived and the comparison retired, so the key was inverted — the default
owns the clean name and `--legacy-unsupported-zero` produces `_legacyzero`. Left as it was, the
clean name would have sat empty and free for a future legacy run to claim, which is exactly the
collision the suffixes exist to prevent.

### The canonical baseline, and what is legacy

**Canonical:** `output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only`,
produced by `run_interpolation_benchmark.jl full --nested-covariates` on the dataset rebuilt
2026-09-02 — 13471 common hours spanning 2022-06-01T09 → 2025-01-01T08. Which commit produced it
is recorded *in the run*: `benchmark_scope.csv` now carries `git_commit`, `git_branch` and
`git_dirty`. Read those rather than inferring from a date, and treat `git_dirty = true` as saying
the run is not citable.

**Legacy:** every result computed on the 11426-hour grid ending 2024-10-01 is retired under
`output/_legacy_11426h_20260829/`, together with the matching `benchmark_diagnostics/`
subdirectories. See the `README.txt` there for what each one was.

Do **not** run `scripts/compare_benchmark_runs.jl` against anything in that archive. It rebuilds
the time grid from current data (13471 hours) and would intersect it with stored `oof_*.csv`
covering 11426, so it would print a delta between two different cell populations and nothing would
mark it as such.

### Partial: `--equal-grids` widened ADW's search and made ADW worse

Run `full --nested-covariates --equal-grids --satellite-wet-blend` on 2026-09-04, stopped
deliberately after `balanced_spatial/FY4B`. The directory
`output/interpolation_benchmark_full_joint_covariates_nested_equalgrids_mgwrintercept_only_satwetblend/`
holds that one cell and **no aggregates or `benchmark_scope.csv`**, so it is not citable; the
figures below come from its stored `oof_*.csv` scored over its own
`common_evaluation_mask.csv`. That reading path was validated first by reproducing the completed
`_satwetblend` run's `metrics_pooled.csv` exactly — all five methods, all five digits.

| method | normal grid | equal grids | change |
|---|---|---|---|
| `adw` | 0.87412 | 0.87766 | **+0.41% worse** |
| `idw` | 0.87688 | 0.87692 | unchanged (`oof_idw.csv` byte-identical) |
| `tps` | 0.93436 | 0.93439 | unchanged (`oof_tps.csv` byte-identical) |
| `mgwr` | 0.96391 | 0.97093 | +0.73% worse |
| `residual_gwr` | 1.00697 | 1.00166 | 0.53% better |
| `blend_mgwr` | 0.87525 | 0.87918 | +0.45% worse |
| `blend_residual_gwr` | 0.88913 | 0.88565 | 0.39% better |

The experiment was queued in case the GWR family's margin came from searching a wider grid than
the baselines were allowed. On this cell the opposite happened: **a wider search made the binding
baseline worse.** `adw` selected a parameter that scores better on the inner split and worse on
held-out cells — more search budget bought selection variance, not accuracy. `idw` and `tps` did
not move at all, their selections being interior to the old grid, which is the useful negative:
only `adw` was ever pinned against its endpoint, so "the baselines were handicapped by a narrow
grid" was true of exactly one of the three and cost it nothing.

The mask lost 237 cells of 3,020,184 (0.008%), which cannot account for the movement: dropping
cells can only lower pooled SSE, and the RMSE rose.

Two cautions on how far this reaches. It is **one cell, and the null product** — FY4B is where the
blend was flat either way, and the margin being tested (+2.94% GPM, +2.04% GSMaP) lives in cells
this run never reached. And the effect is not uniform across the family: the `mgwr` pair got worse
while the `residual_gwr` pair got better, so this is selection noise moving several ways at once,
not a single mechanism. Within the run the relative standing was unchanged — `blend_mgwr` sat
0.17% behind `adw`, against 0.13% on the normal grid.

### Queued experiments

Three remain. Each needs no new code and lands in its own suffixed directory, so none can
overwrite the baseline. Budget ~5 h each on the measured hardware, except `--equal-grids` at
~26 h. Run them only once the canonical baseline above stands.

- `--mgwr-grouping shared` settles two open questions at once: F5's decomposition (is
  `:intercept_only` mgwr a distinct model, or a nested extension of `mixed_gwr`?) and F6's
  corollary (how much of mgwr's advantage over `mixed_gwr` is its ability to demote an over-eager
  `local` role assignment to `bw = Inf`, rather than multiscale resolution?).
- ~~`--free-satellite-coefficient` (`_freesat`)~~ — **done, negative, abandoned** on 2026-09-03:
  −3.25% against `adw` where it completed, and a hard F8 crash on GPM. `--satellite-wet-blend`
  (`_satwetblend`) is the lever that worked, and is now a tuned method rather than an experiment.
  See F4 for both.
- `--equal-grids` (`_equalgrids`) — whether the margin survives giving IDW/ADW/TPS the same
  search budget. **Started 2026-09-04 and stopped after one of six cells**; the one cell is
  recorded just above and points the opposite way to the worry that motivated it. Needs re-running
  to GPM and GSMaP to be decisive, at ~26 h rather than ~5: the wider grid costs about 37% more per
  cell (FY4B 4 h 07 against 3 h 00).
- The paired F1 gate described under "Outstanding verification" above.

### Re-running the review's diagnostics

`scripts/verify_shared_mask_composition.jl` and `scripts/verify_role_vs_bandwidth.jl` reproduce
the measurements behind E1/E2 and F6 from a completed run directory.

`scripts/verify_dry_gate_replay.jl` and `scripts/verify_anchor_discount_bounds.jl` reproduce F4's
measurement. Both rescore a completed run's stored `oof_*.csv` under a counterfactual instead of
refitting, so they need only the run directory and `data/processed/study_area/` — they run while
the covariate inputs are missing. Each asserts an identity that must return the stored prediction
unchanged (`gate = 1`, `a = 1`) before reporting anything; treat a failure there as a broken
reading path, not a finding. Outputs land beside the other diagnostics as `dry_gate_replay.csv`,
`satellite_quadrant.csv` and `anchor_discount_bounds.csv`.

`scripts/verify_banded_blend_bounds.jl` is the fourth of these, and the only one that scores the
pooled gate and the heavy gate together — which is what showed they are opposed (see F4). It
sweeps a blending weight per band of `y_sat`, and per band of `nearest_train_km`, over the same
stored matrices, and writes `banded_blend_bounds.csv` and `banded_blend_sweep.csv`. It runs in
~4 minutes. Its checks are worth keeping if it is ever extended: a zero weight must reproduce
every stored prediction bit for bit; a constant weight must equal an independently written
`satellite_wet_blend`; the per-band minimisation must match a joint enumeration over a coarse
grid; and its constant-weight rows must agree with `anchor_discount_bounds.csv`, which on the
canonical baseline they do to 0.0 RMSE over 198 rows.

All four hold the run's selected bandwidth, kernel and shrinkage fixed, so none can see how the
tuner would re-select them, and all four pick their parameters on the held-out cells. They bound
what is worth building; they are not results. `_satwetblend` gave up 0.04–0.40 points between one
of these ceilings and honest per-fold tuning with a single parameter, so a bound that only just
clears a gate has not cleared it.

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
