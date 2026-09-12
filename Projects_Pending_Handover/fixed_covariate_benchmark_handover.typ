#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2.0cm),
  numbering: "1",
)

#set text(
  font: ("Segoe UI", "Arial", "Libertinus Serif"),
  size: 10pt,
  lang: "en",
)

#set par(justify: true, leading: 0.65em)

#show raw.where(block: true): it => block(
  fill: rgb("f4f4f4"),
  inset: 8pt,
  radius: 3pt,
  width: 100%,
  text(size: 8.5pt, font: ("Consolas", "DejaVu Sans Mono"), it),
)
#show raw.where(block: false): it => text(
  size: 9pt, font: ("Consolas", "DejaVu Sans Mono"), fill: rgb("1a4d80"), it,
)

#set heading(numbering: "1.1")
#show heading.where(level: 1): it => {
  set text(size: 13pt, weight: "bold")
  block(above: 1.2em, below: 0.5em)[#it]
}
#show heading.where(level: 2): it => {
  set text(size: 11pt, weight: "bold")
  block(above: 0.9em, below: 0.35em)[#it]
}

#let note = box(inset: (x: 4pt, y: 1pt), fill: rgb("fff3cd"), radius: 2pt)[
  #text(size: 8.5pt, weight: "bold")[NOTE]
]
#let bad(t) = text(fill: rgb("b31d28"), weight: "bold")[#t]

#align(center)[
  #text(size: 17pt, weight: "bold")[Fixed-Covariate Benchmark — Post-Run Verification Handover]

  #v(0.3em)
  #text(size: 9.5pt, fill: rgb("555555"))[
    Branch `diag/fixed-covariate-benchmark` · Commit `81f1b5a` · 2026-09-12
  ]
]

#v(0.6em)

#block(fill: rgb("eef4fb"), inset: 9pt, radius: 3pt, width: 100%)[
  *Purpose.* This document carries everything needed to verify and interpret the
  `run_interpolation_benchmark.jl full --no-nested-covariates` run once it completes, without
  reference to any prior conversation. Run section 3, verify with section 4, compare with
  section 5, read out with sections 6--7.
]

= Background and the question

In the current benchmark the joint covariates (DEM / ERA5 / NDVI) are #strong[re-selected inside
every outer training fold, separately for each satellite product] (`--nested-covariates`, the
default whenever `--legacy-dem` is not given).

The question to answer: #strong[if the covariates are fixed instead — one set per product — and
the benchmark is re-run, do the GWR-family products get better?]

#strong[No code changes are required.] `scripts/run_interpolation_benchmark.jl:272-278` already
parses a tri-state flag, and `--no-nested-covariates` reads a fixed full-data specification file.
The only difference between the two paths is where the `Dict{group => "local"|"global"}` reaching
`build_joint_fold_context` comes from (`src/benchmark/InterpolationBenchmarkRun.jl:514-530`).
Everything downstream — scaling, bandwidth selection, coefficients, predictions — is identical,
and still uses training stations only.

== What the fixed set actually differs from (verified)

The 7 rows with `final_included = true` in
`output/joint_variable_selection/full/joint_final_full_data_spec.csv`:

#table(
  columns: (auto, auto, auto),
  align: (left, left, left),
  [*Product*], [*Fixed covariates*], [*Roles*],
  [FY4B], [`d2m_c`, `u10`], [both local],
  [GPM], [`elevation`, `u10`, `sp_hpa`], [elevation global, rest local],
  [GSMaP], [`u10`, `sp_hpa`], [both local],
)

Against the per-fold selections under `balanced_spatial` (`joint_role_stability.csv`; the
majority set is defined as variables chosen in at least 3 of 5 folds):

#table(
  columns: (auto, auto, auto, auto),
  align: (left, left, left, center),
  [*Product*], [*Per-fold majority set (>=3/5)*], [*Fixed set*], [*Agree?*],
  [FY4B], [`d2m_c` (4/5), `u10` (4/5)], [`d2m_c`, `u10`], [yes],
  [GPM], [`elevation` (4/5), `u10` (3/5), `sp_hpa` (3/5)], [`elevation`, `u10`, `sp_hpa`], [yes],
  [GSMaP], [`d2m_c` (3/5)], [`u10`, `sp_hpa`], [#bad[no]],
)

#note GSMaP is the only substantive disagreement. `sp_hpa` is in the fixed set but was selected
in #strong[0 of 5 folds] under `balanced_spatial`, while `d2m_c` — the most frequently selected
variable there (3/5) — is absent from the fixed set. Expect the movement in this experiment to
concentrate on GSMaP.

The actual per-fold selections (`joint_fold_quality_control.csv`, `balanced_spatial`):

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (left, left, left, left, left, left),
  [*Product*], [*fold 1*], [*fold 2*], [*fold 3*], [*fold 4*], [*fold 5*],
  [FY4B], [d2m_c, u10], [t2m_c, d2m_c, u10], [d2m_c, u10], [elevation, ndvi], [d2m_c, u10],
  [GPM], [u10], [elevation, sp_hpa], [elevation, t2m_c, u10], [elevation, v10, sp_hpa, ndvi], [elevation, u10, v10, sp_hpa],
  [GSMaP], [d2m_c], [d2m_c, u10], [#bad[(empty)]], [elevation, slope, ndvi], [d2m_c, u10],
)

Selection is markedly unstable under spatial blocking: GSMaP fold 3 selects nothing at all, and
its fold 4 set (`elevation+slope+ndvi`) is disjoint from every other fold. No GPM variable is
selected in all five folds. Under the `:random` scheme the same selections are near-perfectly
stable, so the instability is a property of spatial blocking rather than of the search algorithm.

= Current state (already done — do not repeat)

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Item*], [*Status*],
  [Branch], [`diag/fixed-covariate-benchmark` at `81f1b5a`, the same commit as the nested baseline. Working tree clean apart from untracked plotting scripts.],
  [Specification], [Regenerated at HEAD. #strong[Byte-identical] to the 2026-09-02 version — all 15 files in the directory, `diff -rq` clean.],
  [Spec md5], [`d5e78608b15024a0780d0e8b11d394c9`],
  [Spec sha256], [`358959a159a9b6893159cada2b474a19d3454cfb51d14e74e19181e691941095`],
  [Backup of old spec], [`output/_joint_variable_selection_full_20260902/`],
  [Prerequisites], [DEM / ERA5 / NDVI `output/*_variable_selection/full/` all record `time_count=8069` with status `ok`; `prerequisite_audit()` passes.],
  [The benchmark itself], [#bad[NOT YET RUN.] Three launch attempts were killed by a session memory watchdog. All partial output directories were removed; no stray Julia processes remain.],
)

== What the specification regeneration established

The regenerated output is byte-identical to the 2026-09-02 version, which means:

+ `d63e485` (Jun--Sep grid re-pinned from 8067 to 8069) and `81f1b5a` (NDVI requires an
  `ndvi_land_qc` land pixel) #strong[did not change the selection outcome]. Both were
  re-derivations onto master of commits that already existed on
  `fix/local-hat-unsupported-target-nan`, and the 2026-09-02 selection run was made on that
  branch — so it already had both behaviours.
+ As a by-product, this is a reproducibility check: a 999-permutation, 8069-hour selection
  pipeline produced bit-identical output across two runs ten days apart.

The fixed-vs-nested comparison is therefore #strong[free of any code-version confound] and can
proceed.

== Why the three launch attempts failed (does not apply to a manual run)

Sampler measurements, 20 seconds apart — only two points were captured before the sampler itself
was killed:

#table(
  columns: (auto, auto, auto),
  align: (left, right, right),
  [*Time*], [*Julia RSS*], [*System free memory*],
  [03:39:52], [674 MB], [3883 MB],
  [03:40:13], [1115 MB], [3442 MB],
)

Julia was killed at 1.1 GB RSS with 3.4 GB still free — nowhere near any heap limit. The cause is
that the Claude Code session's low-memory watchdog triggers on #strong[system] free memory, and
this machine idles at roughly 4.3 GB free across 376 processes, so Julia taking even 1 GB crosses
the threshold. Adjusting `--heap-size-hint` or the thread count does not help, because the
problem is not Julia's. #strong[A manual run in an ordinary terminal is not subject to this
watchdog], and the machine has a 27 GB pagefile whose peak usage has only ever reached 1.89 GB.

= The command to run

Use an ordinary, #strong[non-elevated] terminal:

```powershell
cd C:\Users\qw123\GG\GeoWeightedRegression.jl
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --no-nested-covariates --with-random 2>&1 | Tee-Object output\log_benchmark_full_fixedcov_20260912.txt
```

Output lands in `output/interpolation_benchmark_full_joint_covariates_mgwrintercept_only/`, a
fresh directory. Neither the nested baseline nor the `_nocov` run is touched.

== Four constraints that must not be changed

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Constraint*], [*Reason*],
  [`-t 4`, never `-t auto`], [The dominant loop `dynamic_covariate_predict` allocates ~1.6 MB per task; past a handful of threads, large-object allocation contention outweighs the extra cores. Measured on this box: 4 threads 167 ms, 8 threads 201 ms, 24 threads 274 ms. On the smoke benchmark `-t 4` took 314 s against 393 s for `-t auto`, with byte-identical output.],
  [No `--heap-size-hint`], [That was a workaround for the session watchdog. It is unnecessary in a manual run; let Julia manage its own heap.],
  [Never set `BLAS.set_num_threads`], [`tps` factorizes a dense (n+3)x(n+3) system, and a threaded LU accumulates in a different order. Pinning BLAS to one thread was measured to move `tps`'s RMSE in the last three or four digits, propagating into `metrics_pooled.csv` and others. The run's BLAS thread count is part of what makes published numbers reproducible.],
  [Never "Run as administrator"], [An elevated process holds `SeBackupPrivilege`/`SeRestorePrivilege`, which bypass the deny-DELETE ACL on `data/` and silently remove that protection layer entirely.],
)

This run should be #strong[faster] than the nested baseline: it skips 999-permutation per-fold
selection across 30 fold-cells, and skips 2000 bootstrap replicates.

= Post-run verification checklist

Work through these in order. Resolve any failure before moving on to section 5.

== Check 1 — run integrity

```powershell
echo $LASTEXITCODE          # must be 0
Get-Content output\log_benchmark_full_fixedcov_20260912.txt -Tail 30
```

The log tail should show normal completion with no `ERROR`, `StackOverflow`, or `OutOfMemory`.

#strong[Freshness check (important).] The output directory's newest mtime must be later than
`src/`'s. This is exactly what `assert_run_is_fresh` in `scripts/verify_perf_invariance.jl:64`
exists for: if the pipeline errors partway, the previous run's output is still on disk and a
naive diff passes while telling you nothing.

```powershell
$run = Get-ChildItem output\interpolation_benchmark_full_joint_covariates_mgwrintercept_only -Recurse | Sort-Object LastWriteTime -Desc | Select-Object -First 1
$src = Get-ChildItem src -Recurse -File | Sort-Object LastWriteTime -Desc | Select-Object -First 1
"run=$($run.LastWriteTime)  src=$($src.LastWriteTime)  fresh=$($run.LastWriteTime -gt $src.LastWriteTime)"
```

Expect `fresh=True`. #note This command was trialled against the old nested baseline and returned
`fresh=False` (run at 09-11 01:49, `src/` last modified 09-11 03:49) — that is the check working,
not a broken command. A new run is necessarily later than `src/`, so it must report `True`. If
branches are switched in between, `git checkout` updates mtimes under `src/`; in that case re-run
the benchmark rather than waiving the check.

== Check 2 — provenance

```powershell
cd output\interpolation_benchmark_full_joint_covariates_mgwrintercept_only
Get-Content joint_spec_provenance.csv
Select-String -Path benchmark_scope.csv -Pattern "git_commit|git_branch|git_dirty|results_admissible"
```

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Field*], [*Expected*],
  [`joint_spec_provenance.csv` -> `selection_mode`], [`fixed_full_data`],
  [same -> `confirmatory`], [`false`],
  [same -> `sha256`], [`358959a159a9b6893159cada2b474a19d3454cfb51d14e74e19181e691941095`],
  [same -> `source_path`], [points at `output\joint_variable_selection\full\joint_final_full_data_spec.csv`],
  [`benchmark_scope.csv` -> `git_commit`], [`81f1b5a4c65e7b7dfc64f86ae9b6e6d9c0ba7d90`, matching the nested baseline],
  [same -> `git_dirty`], [`false`],
  [same -> `git_branch`], [`diag/fixed-covariate-benchmark`. This differs from the baseline's `fix/estimator-correctness`, which is expected and harmless — `git_commit` is the field that matters.],
)

Confirm the fixed specification was genuinely consumed:

```powershell
fc.exe joint_variable_spec_used.csv ..\joint_variable_selection\full\joint_final_full_data_spec.csv
```

== Check 3 — the run declares its own results inadmissible

```powershell
Select-String -Path benchmark_scope.csv -Pattern "results_admissible"
```

Expect a value beginning `NO - exploratory_only=true:`
(source: `InterpolationBenchmarkRun.jl:262-263`), against the nested baseline's
`yes - every selection step ran inside its training fold`.

#strong[If this reads "yes", the wrong mode was run] and the experiment is void.

== Check 4 — files that should be absent (absence is the check)

These will #strong[not] exist. That is by design, not failure, and their absence is itself a
verification:

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Should not exist*], [*Why*],
  [`paired_comparisons.csv`], [`bootstrap_reps=0` (`run_interpolation_benchmark.jl:199`); the write is guarded by `nrow(bootstrap) > 0` (`InterpolationBenchmarkRun.jl:196`).],
  [`claim_assessment.csv`, `claim_agreement.csv`], [`claim` requires `bootstrap_reps > 0` (`:171`); the write is guarded by `ncol(claim) > 0` (`:220-226`).],
  [`joint_fold_roles.csv`, `joint_fold_candidates.csv`, `joint_fold_vif.csv`, `joint_fold_spatial_variability.csv`, `joint_role_stability.csv`], [`_write_joint_selection_outputs` is called only when `nested_joint` is true (`:218`). In fixed mode there is no per-fold selection process to record.],
)

```powershell
foreach ($f in "paired_comparisons.csv","claim_assessment.csv","claim_agreement.csv","joint_fold_roles.csv","joint_role_stability.csv") {
  "{0,-38} {1}" -f $f, $(if (Test-Path $f) { "present -> UNEXPECTED" } else { "absent  -> correct" })
}
```

Files that #strong[should] be present: `metrics_pooled.csv`, `metrics_folds.csv`,
`metrics_stratified.csv`, `metrics_fold_summary.csv`, `method_rank_stability.csv`,
`parameter_scan.csv`, `run_status.csv`, `covariate_model_status.csv`, `benchmark_scope.csv`,
`auto_selection.csv`, `joint_variable_spec_used.csv`, `joint_spec_provenance.csv`,
`joint_bandwidths.csv`, `joint_fold_scaling.csv`, `joint_fold_quality_control.csv`,
`joint_era5_input_qc.csv`, `joint_ndvi_alignment_qc.csv`, plus the `balanced_spatial/` and
`random/` subdirectories.

== Check 5 — the fixed set actually took effect

Every residual-model row in `run_status.csv` / `covariate_model_status.csv` must carry
`covariate_selection_mode = fixed_full_data`, and `covariate_variables` must be
#strong[constant across all 5 folds] within each product:

```powershell
Import-Csv covariate_model_status.csv |
  Where-Object { $_.scheme -eq "balanced_spatial" } |
  Select-Object product, fold, method, covariate_selection_mode, covariate_variables |
  Sort-Object product, fold, method | Format-Table -AutoSize
```

#table(
  columns: (auto, auto),
  align: (left, left),
  [*Product*], [*Expected `covariate_variables` (identical in all 5 folds)*],
  [FY4B], [`d2m_c,u10`],
  [GPM], [`elevation,u10,sp_hpa`],
  [GSMaP], [`u10,sp_hpa`],
)

If that column varies between folds, the fixed path did not take and the experiment is void.

#note `covariate_effective_roles` and `covariate_selected_roles` may legitimately differ.
`residual_gwr` forces every selected group to `local` (`joint_effective_roles`,
`JointCovariateModels.jl:318`), while `mixed_gwr` and `mgwr` honour the specification's roles.
So GPM's `elevation=global` affects only the latter two. A consequence worth anticipating: when a
product's role map contains no `global` entry, `residual_gwr` and `mixed_gwr` are the same model
end to end, and the `duplicate_of` column in `run_status.csv` will say so. Both FY4B and GSMaP
fixed sets are all-local, so expect that pairing there; GPM has a global role and should not
collapse.

== Check 6 — sanity check (the most important one)

`idw`, `adw`, `tps` and direct `gwr` #strong[never touch the covariate path]. The gate at
`InterpolationBenchmarkRun.jl:12-13` requires `mode == "residual"` and a method in
`gwr`/`mixed_gwr`/`mgwr` for `uses_joint` to be true. On the intersected mask, therefore, these
four must come out #strong[unchanged] relative to the nested baseline.

In the section 5 comparison table their `delta_paired` should be zero or last-digit noise.
#strong[If they move materially, something other than the covariate set changed and the entire
comparison is invalid.]

= Comparison

```powershell
cd C:\Users\qw123\GG\GeoWeightedRegression.jl
julia -t 4 --project=. scripts/compare_benchmark_runs.jl output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only output/interpolation_benchmark_full_joint_covariates_mgwrintercept_only fixed_vs_nested_covariates
```

Why this script rather than subtracting two `metrics_pooled.csv` files: the common evaluation
mask keeps only cells where every method in `MASK_METHODS` is finite, so any change in one
method's coverage resizes the denominator for all of them — and the cells that move are
systematically the hard ones. This script intersects the two masks and compares pairwise.

Results are written to
`output/benchmark_diagnostics/fixed_vs_nested_covariates/run_comparison.csv`. Key columns
(defined in `src/benchmark/diagnostics/run_comparison.jl:17`):

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Column*], [*Meaning*],
  [`n_paired`], [cells finite in both runs],
  [`RMSE_paired_before`], [nested baseline on the intersection],
  [`RMSE_paired_after`], [fixed-spec run on the intersection],
  [`delta_paired`], [`after - before`; #strong[negative means the fixed set is better]],
  [`relative_paired`], [the same as a fraction],
  [`mask_cells_before` / `_after` / `_shared`], [mask sizes, for judging comparability],
)

== Do not compare against the `_nocov` run

#block(fill: rgb("fdecea"), inset: 9pt, radius: 3pt, width: 100%)[
  `output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_nocov/`
  records `git_commit = 8a83df3` on branch `feat/fused-anchor-agreement-blend` in its
  `benchmark_scope.csv`. That commit is #strong[not an ancestor of `81f1b5a`] and predates all
  five of the following:

  #table(
    columns: (auto, 1fr),
    align: (left, left),
    [*Commit*], [*Change*],
    [`35324db`], [Form IDW/ADW weights over the stations that actually report],
    [`ffd5c6d`], [Fix TPS's lambda scale once per call, not per missing-value group],
    [`8f83d13`], [Record an unfittable local target as missing, not as a zero correction — directly affects the joint/residual path],
    [`81f1b5a`], [Require a land pixel for the NDVI covariate],
    [`d63e485`], [Re-pin the Jun--Sep grid, gate auto on coverage, record the duplicate model],
  )

  So the difference in `idw`/`adw`/`tps` between those two runs comes from the estimator fixes
  themselves, not from a different number of evaluated cells; and the GWR-family difference
  cannot be separated from `8f83d13`. `compare_benchmark_runs.jl` aligns #strong[masks]; it
  cannot correct a #strong[code version] difference.

  The claim "removing covariates is better" therefore #strong[has no supporting evidence at
  present] and must not appear in any conclusion. Settling it requires a same-commit
  no-covariate run — see section 8.
]

= Baseline numbers

Nested baseline pooled RMSE (`balanced_spatial`, `group=overall`, `level=all`, empty `fold`).
Extraction command:

```powershell
$base = "output\interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only\metrics_pooled.csv"
Import-Csv $base | Where-Object { $_.scheme -eq "balanced_spatial" -and $_.group -eq "overall" -and $_.level -eq "all" -and $_.fold -eq "" } |
  Select-Object product, method, n, RMSE, MAE, Bias, r | Format-Table -AutoSize
```

#table(
  columns: (auto, auto, auto, auto),
  align: (left, right, right, right),
  [*Method*], [*FY4B*], [*GPM*], [*GSMaP*],
  [`raw` (bare satellite)], [1.3165], [0.9616], [1.2192],
  [`zero`], [1.0045], [1.0011], [1.0011],
  [`train_clim`], [0.9954], [0.9920], [0.9920],
  [`hour_field_mean`], [0.9427], [0.9395], [0.9395],
  [`idw`], [0.8769], [0.8737], [0.8737],
  [*`adw`*], [*0.8741*], [*0.8709*], [*0.8709*],
  [`tps`], [0.9344], [0.9309], [0.9309],
  [`gwr` (direct)], [0.9256], [0.9223], [0.9223],
  [`residual_gwr`], [1.0070], [0.9332], [1.0065],
  [`mixed_gwr`], [1.0061], [0.9324], [1.0060],
  [`mgwr`], [0.9639], [0.9173], [0.9590],
  [`auto`], [0.9639], [0.9173], [0.9621],
  [cells `n`], [3,020,184], [3,072,559], [3,072,559],
)

Current position: #strong[`adw` wins on all three products, and no GWR-family method beats `idw`
or `adw`.] On FY4B and GSMaP, `residual_gwr` and `mixed_gwr` are worse than the `zero` predictor.
The best GWR variant is always `mgwr`.

Fold spread for context (`metrics_fold_summary.csv`, `balanced_spatial` / `overall` / `all`,
FY4B): `adw` is 0.8742 +/- 0.0318 (range 0.8504--0.9260) and `mgwr` is 0.9634 +/- 0.0245
(0.9416--0.9915). The 0.089 gap is about 2.8 fold standard deviations, so it is not fold noise.

#note `RMSE_std` is a dispersion diagnostic, #strong[not a standard error] (see the docstring on
`summarize_fold_spread`). Do not use it for significance.

= Expected results and how to read them

== Expectations

+ #strong[FY4B and GPM should barely move]: the fixed set equals the per-fold majority set
  (section 1.1), so only individual folds change — FY4B fold 4, GPM folds 1 and 4.
+ #strong[Movement should concentrate on GSMaP]: its fixed set contains `sp_hpa` (selected 0/5
  per fold) and omits `d2m_c` (selected 3/5).
+ #strong[Method ranking is not expected to change]: `adw` should still lead, with the GWR family
  behind.

== The boundary on interpretation (must accompany any conclusion)

The fixed specification was selected using #strong[all 237 stations], including the ones each
fold holds out and then predicts. The source records this itself
(`InterpolationBenchmarkRun.jl:296-315`, written into `benchmark_scope.csv` in fixed mode):

#block(inset: (left: 12pt), stroke: (left: 2pt + rgb("cccccc")))[
  `covariate_selection_leakage`: non-nested variable selection used all 237 stations before
  spatial cross-validation

  `inference_policy`: exploratory performance metrics only; no paired significance or claim
  assessment
]

Therefore:

- Any RMSE improvement this run shows is #strong[partly selection leakage, not a genuine
  out-of-sample gain].
- This run #strong[cannot replace] the reported benchmark. The nested run remains the headline
  result; this one answers a diagnostic question about it.
- That is also why `bootstrap_reps` is zero and no significance test or claim assessment is
  produced.

== Read-out template

The final deliverable is a table of this shape, per product and per GWR-family method:

#table(
  columns: (auto, auto, auto, auto, auto),
  align: (left, left, right, right, right),
  [*Product*], [*Method*], [*nested*], [*fixed*], [*delta*],
  [FY4B], [`residual_gwr` / `mixed_gwr` / `mgwr`], [tbd], [tbd], [tbd],
  [GPM], [same], [tbd], [tbd], [tbd],
  [GSMaP], [same], [tbd], [tbd], [tbd],
)

Plus one sentence: whether fixing the covariates improved the GWR family, and whether the
movement was large enough to change any ranking.

= Follow-up work, by priority

== Priority 1 — a same-commit no-covariate run

This is the experiment that would actually answer whether the covariates help at all. The
`--no-covariates` flag was added in `8a83df3` (present on
`feat/fused-anchor-agreement-blend`, `cleanup/dead-code-and-exports` and `refactor/src-layout`)
and is #strong[not on master]. Cherry-pick `8a83df3` onto the current branch, re-run the
no-covariate benchmark on `81f1b5a` code, and compare against the nested baseline. Only that
difference is free of a code-version confound.

== Priority 2 — the heavy-rain lead

This is probably more promising than the covariate set. `paired_comparisons.csv` shows
`residual_gwr` significantly #strong[ahead] of all three baselines in the heavy stratum
(>= 8 mm/h):

#table(
  columns: (auto, auto, auto, auto),
  align: (left, right, right, right),
  [*Product*], [*vs idw*], [*vs adw*], [*vs tps*],
  [FY4B], [+4.67%], [+4.67%], [+3.01%],
  [GPM], [+4.70%], [+4.69%], [+2.99%],
  [GSMaP], [+6.71%], [+6.71%], [+5.05%],
)

All with `pvalue_holm <= 0.003`. It is significantly worse in the overall and moderate strata.
The reason is that pooled RMSE is dominated by the roughly 93% of cells with no rain — for FY4B
the `no_rain` stratum alone is 2,808,815 of 3,020,184 cells. In other words, #strong[the GWR
family already wins where it should, and loses on the average over dry cells.] Worth considering:
stratified reporting, or explicit treatment of dry cells (`hurdle_gwr` is implemented in
`InterpolationBenchmarkHurdle.jl` but deliberately excluded from `BENCHMARK_RUNS`).

== Priority 3 — two known stale items

+ `scripts/verify_interpolation_benchmark_nested_covariates.jl:12-17` carries a self-documented
  STALE note: it looks for a `..._nested` directory, but the runner now appends
  `_mgwrintercept_only`, so it reports every output as missing on a default run.
+ `output/ndvi_variable_selection/full/` dates from 2026-09-02, so `prerequisite_full_audit.csv`
  records stale provenance. This is #strong[harmless] — the joint selection run recomputes NDVI
  alignment at HEAD, and the audit gate only checks `time_count` and status, both of which pass.
  Refreshing it for the record costs about 45 minutes.

#v(1.5em)
#line(length: 100%, stroke: 0.5pt + rgb("cccccc"))
#align(right)[
  #text(size: 8.5pt, fill: rgb("777777"))[
    Every assertion here traces to a verified source line or an output file already on disk.\
    No numbers are predicted for the run that has not yet happened.
  ]
]
