# Handover — grid support diagnostics (point-vs-pixel cost of the benchmark's satellite inputs)

**Date:** 2026-09-19 · **State:** complete and verified, **uncommitted** · **Next decision:** open

---

## 1. Why this exists

An advisor objected that the satellite precipitation products are inherently **gridded**, and asked
whether the workflow extracts grid values at gauge points and then interpolates those point values
again — remedy being a downscaling step, or at minimum a resample of the coarse grid onto a finer
one, **scoped to the benchmark**.

**Half the premise is wrong and the correction matters.** The satellite is *not* re-interpolated
from station values. It enters as a point-extracted **anchor** at the target station:
`residual_*` methods fit `y_obs − y_sat` on training gauges and predict `y_sat_target + residual`
(`src/benchmark/InterpolationBenchmarkPredictors.jl:130,211`). What is interpolated is the **gauge
residual field**, not the satellite.

**The other half is right:** validation is only at gauge points, so the deliverable is a point
correction rather than a corrected field; and the anchor carries no sub-cell information, because
gauges sharing a cell are handed identical values.

**Decision already taken (by you, 2026-09-19):** diagnose before building. No downscaling
formulation chosen yet, FY4B only for anything grid-side, benchmark output plus scope rows for
reporting. This handover covers the diagnostic stage only.

Approved plan: `C:\Users\qw123\.claude\plans\my-advisor-told-me-resilient-kernighan.md`

---

## 2. What was delivered

| File | Status |
|---|---|
| `src/benchmark/GridSupportDiagnostics.jl` | new, untracked |
| `scripts/run_grid_support_diagnostics.jl` | new, untracked |
| `test/test-grid-support-diagnostics.jl` | new, untracked |
| `src/load_modules.jl` | +2 lines (dependency entry) |
| `test/runtests.jl` | +1 line (include) |
| `CLAUDE.md` | module list + count — see §7 |

A standalone module in the CLAUDE.md sense (2): own `module … end`, loaded by name, registered
only in `load_modules.jl`'s dependency list. It **measures inputs and changes no benchmark
output** — no config field, no CLI flag on the benchmark, no new method or product, `MASK_METHODS`
untouched, no published number moved.

### Verification state

- Full suite green: `julia --project=. -e 'using Pkg; Pkg.test()'` — 11 new testsets, 117
  assertions, plus every pre-existing test.
- `data/` untouched: a recursive `find data -newermt "-12 hours"` returned **nothing**.
- Only `output/grid_support_diagnostics/` was written under `output/`.
- D7 reproduces the shipped FY4B column **bit-for-bit** (`reference_mismatch = 0` over 2022-06),
  which is what makes the nearest-vs-bilinear comparison controlled.

### How to re-run

```sh
# main checkout only - it reads data/, so never in a worktree (CLAUDE.md data-safety rules)
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/run_grid_support_diagnostics.jl              # D7 over 2022-06, ~10 min
julia --project=. scripts/run_grid_support_diagnostics.jl --window none  # D1-D6 only, ~4 min
julia --project=. scripts/run_grid_support_diagnostics.jl --window full  # D7 over the record, hours
```

`--window` bounds **D7 only** (the one step that re-reads the NetCDF archive: ~2,900 files for one
month, ~160,000 for the record). `output/grid_support_diagnostics/run_scope.csv` records which
window produced the tables on disk.

---

## 3. Findings

### 3a. Two structural surprises — these change how the inputs must be described

**FY4B is not on one projection.** The archive holds **two sub-satellite longitudes**: 65,909 files
at 105.0°E and 59,334 at 133.0°E. `extract_precipitation` reads each file's own subpoint
(`dataset_satellite_lon`), so a gauge has *two* pixels, and two gauges receive the same value only
where they share a pixel in **both**. That intersection (`fy4b_pixel_effective`) explains the data
exactly — 12 cell groups, 12 series groups, 25 gauges either way. A single-projection view
overstates the merging (24 groups / 49 gauges at 133°E alone).

Consequence for resolution claims: the effective footprint **differs by view**. Over the study box
it is **21.7 km² from 105°E (1.35× the 4 km nadir spec)** and **24.6 km² from 133°E (1.54×)** —
never the 16 km² that "4 km" implies. `FY4BPreprocessing.jl:28`'s `RESOLUTION = 4000.0` is a nadir
figure and no report has ever stated what it becomes here.

**No lat/lon grid explains GPM or GSMaP.** The exports ask Earth Engine for `scale: 11132` m
(0.1° at the equator), but `reduceRegions` reprojects rather than reading native cells. Tested
against four candidate grids:

| Product | best candidate | cell groups explained |
|---|---|---|
| FY4B | `fy4b_pixel_effective` | **12 of 12** ✅ |
| GPM | `coarse_0p1deg` | 15 of 56 ❌ |
| GSMaP | `coarse_0p1deg` | 38 of 56 ❌ |

Gauges come out identical *across* assumed cells in every coarse case, so the grid is simply not
recoverable from this repository. **Use the empirical grouping** (`identical_series.csv`), which
assumes no geometry at all:

| Product | gauges it cannot tell apart | share of 237 | largest merged set |
|---|---|---|---|
| FY4B | 25 | 11% | 3 |
| GPM | 70 | 30% | 3 |
| GSMaP | **126** | **53%** | 5 |

### 3b. The numbers that decide the downscaling question

Gauge-vs-gauge disagreement — the floor on point-to-pixel error (a satellite cannot match both
gauges) and the ceiling on what any downscaling could recover:

| pairs | n | RMSE all hours | RMSE wet hours |
|---|---|---|---|
| share an FY4B pixel (median 1.9 km apart) | 14 | **0.520** | **1.831** |
| ≤5 km apart, different pixels | 56 | 0.649 | 2.255 |
| 5–10 km apart | 247 | 0.823 | 2.773 |
| 20–40 km apart | 3,149 | 1.128 | 3.626 |

Nearest-vs-bilinear re-extraction of the same files and hours:

| quantity | value |
|---|---|
| mean absolute difference | 0.013 mm/h |
| RMSE | 0.093 mm/h |
| **RMSE over wet hours** | **0.361 mm/h** |
| largest single-hour difference | 5.32 mm/h |
| wet hours flipping the rain call | 4.9% |
| cells failing to reproduce the shipped column | **0** |

**Read it this way:** resampling moves the anchor by roughly **a fifth** of the sub-pixel
disagreement it can never resolve (0.36 vs 1.83 mm/h on wet hours). Not nothing; not decisive.

**Important caveat before committing to anything.** The same-cell/different-cell contrast is nearly
flat once distance is held fixed — 0.520 vs 0.649 mm/h for pairs *both* under 5 km. **Proximity,
not cell-sharing, is what makes gauges agree.** The sub-pixel floor is essentially the ≤2.5 km
gauge-pair error, which means most of the point-vs-pixel gap is irreducible representativeness
error rather than anything the workflow can fix. A pure resample adds no information and this says
its ceiling is low.

### 3c. Gauge-network bias (relevant if terrain covariates are used to downscale)

The "fraction of the range no gauge occupies" figure (4.7%) badly understates it — the gauges sit
*inside* the terrain's range but low in it:

- gauge mean elevation **622 m** vs domain mean **850 m**
- gauge median **240 m below** the domain median (the gauge median sits near the domain's 25th pct)
- highest gauge 1,894 m vs highest terrain 2,984 m — **1,089 m of relief with no gauge above it**
- slope: gauge mean 21.2° vs domain 25.3°
- only **9–10%** of the FY4B pixels covering the study box contain a gauge at all

Domain distribution is the 264×264 (69,696-cell) 1/120° grid over `STUDY_BOUNDS`, aggregated from
the same Copernicus rasters the station covariates were sampled from, so it is one dataset at two
supports rather than two sources.

---

## 4. Output files

`output/grid_support_diagnostics/`

| File | Contents |
|---|---|
| `summary.csv` | **start here** — 49 headline numbers, one per row |
| `run_scope.csv` | what produced the directory (window, hours, stations, subpoints) |
| `pixel_assignment.csv` | per gauge per projection: pixel, fractional offset, footprint |
| `cell_membership.csv` / `cell_multiplicity.csv` | who shares a cell, per candidate grid |
| `identical_series.csv` | **the authoritative grouping** — no geometry assumed |
| `cell_grouping_agreement.csv` | which candidate grid the data endorses |
| `gauge_pair_disagreement.csv` | every gauge pair under 60 km, with share flags |
| `subcell_disagreement_summary.csv` | the table in §3b |
| `fy4b_footprint.csv` | footprint per projection, at gauges and box corners |
| `gauge_terrain_representativeness.csv` | deciles of gauge vs domain terrain |
| `extraction_sensitivity.csv` | nearest vs bilinear, overall and per gauge |
| `global_common_time_qc.csv` | written by the shared loader (not a diagnostic) |

---

## 5. The open decision, and the design already worked out for it

Nothing about the downscaling formulation has been decided. When you pick one, these were the three
arms on the table:

1. **Resample only** (bilinear/cubic coarse→fine). Satisfies the letter of the request. §3b says
   its ceiling is low.
2. **Resample + GWR weight field.** Fit the precipitation–terrain relation at coarse scale on
   monthly aggregates, apply at 1 km, residual-correct, normalise to a weight field with mean 1
   inside each coarse cell, then apply hour-by-hour. Mass-conserving, so it adds sub-cell structure
   without smuggling in bias correction. §3c says the terrain relation would be fitted on a biased
   gauge sample — worth stating.
3. **Per-hour GWR downscaling.** 13,471 hours × 69,696 cells; hourly precipitation has a weak,
   noisy elevation signal with ~93% dry cells, so most hours would be near-zero-skill fits
   dominated by the residual term.

### Assets already on disk

| Asset | Use |
|---|---|
| `data/FY4B/{2022,2023,2024}/**/*.NC` | native 4 km QPE — the only product with a local grid |
| `data/dem_ShiYan_1km.tif` | **264×264 at 1/120°, extent exactly `STUDY_BOUNDS`** — ready-made fine grid |
| `data/processed/dem/copernicus_glo30_utm49n_30m.tif` (+ slope, aspect) | gridded terrain predictors |
| `data/raw/gpm_imerg_v07_uncal_hubei_2022_2024/YYYY/DDD/*.nc4` | IMERG 0.1°, 30-min, **unread by any code** |
| `FY4BPreprocessing.jl:189` `latlon_to_xy` | lat/lon → FY4B pixel, so regridding needs no new geometry |
| `TraditionalInterpolation.jl:232,247,318` | `idw/adw/tps_predict` already accept arbitrary targets |
| `GridSupportDiagnostics.bilinear_neighbours` / `bilinear_combine` | QC-aware bilinear kernel, tested |

Gaps: **GSMaP has no gridded source on disk** (GEE point exports only). The local IMERG variable is
`precipitationUncal`, a *different product* from the benchmark's GEE-calibrated `GPM` column — a
field built from it is `GPM_UNCAL`, not a drop-in. Gridded NDVI and gridded ERA5-Land do not exist
locally, so a full *joint-covariate* gridded prediction is blocked; coordinate-only methods
(`idw`, `adw`, `tps`, residual GWR on `[1, lon, lat]`) are not.

### ⚠ Hard constraint for the next stage — read before adding any product

`load_global_common_product_data` (`src/mger/MGERPipeline.jl:145-191`) intersects **stations and
hours across all products**: `common_product_times` keeps only hours where *every* product is
available, asserted against `expected_common_time_count = 13471`. **Adding an entry to
`cfg.mger.sat_paths` would change the common grid and move every incumbent metric.**

So downscaled variants must enter as **extra products reindexed onto the already-fixed grid**
(NaN where absent), never as grid-defining products — the principle that lets `--fused-anchor` add
`MERGED_OLS`/`MERGED_MEAN` without disturbing anything.

The opt-in-flag idiom to copy: CLI parse (`run_interpolation_benchmark.jl:348-355`) → config field
(`InterpolationBenchmarkConfig.jl:151-164`) → **a function of the config, not a longer const**
(`benchmark_methods` / `benchmark_products`, `Config.jl:100-115`) → scope rows (`Run.jl:328-367`) →
a directory suffix (`run_interpolation_benchmark.jl:64-71`). The invariance guard to copy verbatim
is `test/test-interpolation-benchmark.jl:1687-1703`, which asserts every incumbent metric row is
byte-identical with the flag on.

Also note `scripts/verify_perf_invariance.jl` compares byte-for-byte and will report `EXTRA` /
`CHANGED` for any new method even when no incumbent number moved. It is transient by design
(its own docstring says to delete it once the perf pass lands) — re-baseline rather than fight it.

---

## 6. Problems found in existing code — reported, **not** modified

Per CLAUDE.md's "point them out first, do not modify them directly".

1. **`latlon_to_scan_angles`'s limb test never fires** (`src/sources/FY4BPreprocessing.jl:66-76`).
   It compares the angle *at the satellite*, which never exceeds ~8.7° for any point on Earth, so
   `cos_angle < 0.156` is unreachable. Measured: the antipode (−47°E, 0°N) maps to the nadir pixel,
   and (−40°E, 20°N) maps to a valid in-grid pixel. **Harmless for every published number** — all
   237 gauges are well inside the disk, and `extract_precipitation`'s bounds check is the real
   guard — but `latlon_to_xy` cannot be used to flag off-disk points. A test in
   `test-grid-support-diagnostics.jl` pins the current behaviour so a future fix announces itself
   rather than silently changing meaning.
2. **`CLAUDE.md:149`'s "~2.6M-cell dataset" is stale.** The current grid is 13,471 × 237 =
   3,192,627 cells, of which 3,020,184 are scored at 0.946 coverage (checked against the shipped
   `metrics_pooled.csv`). 2.6M belongs to the retired 11,426-hour grid, which
   `run_interpolation_benchmark.jl:144-145` already records as superseded.
3. **`scripts/verify_interpolation_benchmark_nested_covariates.jl` is already red**, unrelated to
   this work: it hard-codes eight method names and a directory name without the
   `_mgwrintercept_only` suffix (`:6-8`, `:69`). Its own header comment (`:12-17`) flags it stale.

---

## 7. Things to look at before committing

- **`CLAUDE.md` edit went slightly beyond the plan.** The file asserts a count of its
  standalone-module list, so adding `GridSupportDiagnostics` required the count to be right — and
  the list was already missing `benchmark/HeavyRainEvents.jl`. Both lines were added. Back out the
  `HeavyRainEvents` line if you'd rather keep that change separate; the count would then be wrong,
  which is why it went in.
  The count now reads **twenty-five and is correct** (verified: 25 files defining `module` under
  `src/`, excluding the `MixedGWR.jl` package entry).
- **`reference_code/`** (an uninitialized submodule — `git submodule status` shows `-`, no
  `.git/modules`) went missing from the working tree during the session. It held nothing, and the
  disappearance could not be attributed to any command run. The empty placeholder was recreated
  with `mkdir`, so `git status` matches the session-start state.
  `git submodule update --init reference_code` would populate it if ever wanted.
- **Parallel work already landed.** A concurrent session's no-rain evaluation was committed as
  `3a5bccf "Add the no-rain (gauge-dry) evaluation: module, drivers, plots, report draft"`
  (9 files: `src/benchmark/NoRainEvaluation.jl`, two `run_no_rain_*.jl`,
  `scripts/plot_no_rain_evaluation.py`, `test/test-no-rain-evaluation.jl`,
  `Interim_results/no_rain_evaluation.typ`, plus `CLAUDE.md` / `load_modules.jl` /
  `runtests.jl` lines). **Nothing of it is yours to commit** — the working tree is now clean of it.
- `scripts/__pycache__/` and `tmp/` were already untracked clutter before this work.
- **`final_plan.typ` shows as deleted** in `git status`. That predates this work — it was already
  staged-deleted at session start. Leave it alone unless you know why it went.
- **`.gitignore:3` is `*.pdf`**, so the existing `fixed_covariate_benchmark_handover.pdf` is
  *ignored* while this Markdown handover is not. Committing it is a choice: it makes the handover
  reviewable in history, but it is working notes, not a deliverable.

### Exactly what is yours to commit

```
new:  src/benchmark/GridSupportDiagnostics.jl
      scripts/run_grid_support_diagnostics.jl
      test/test-grid-support-diagnostics.jl
mod:  src/load_modules.jl   (+2: the dependency entry)
      test/runtests.jl      (+1: the include)
      CLAUDE.md             (count -> twenty-five, + HeavyRainEvents and
                             GridSupportDiagnostics bullets)
opt:  Projects_Pending_Handover/grid_support_diagnostics_handover.md  (this file)
```

Nothing else in `git status` belongs to this work.

---

## 8. What to say to the advisor

1. Correct the premise precisely: the satellite is a point-extracted **anchor**; what is
   interpolated is the gauge residual field. Do not concede a mechanism that is not there.
2. Concede the two real defects: validation is only at gauge points, and the anchor carries no
   sub-cell information — **53% of gauges cannot be told apart by GSMaP, 30% by GPM, 11% by FY4B**.
3. State the resolution facts no report has ever carried: FY4B's 4 km is a **nadir** figure that
   becomes **1.35–1.54× in area** over the study area, from **two different view geometries**; and
   the grid GPM/GSMaP were actually sampled on is not recoverable from the exports.
4. Put the ceiling on the table honestly: gauges inside one pixel disagree with each other by
   **1.83 mm/h on wet hours**, and that is essentially the ≤2.5 km gauge-pair error — so most of
   the point-vs-pixel gap is representativeness error, not retrieval error, and not something a
   resample fixes. Resampling moves the anchor by ~a fifth of that.
5. Then ask which they want: a fine-grid **anchor** (cheap, defuses the scale-mismatch objection),
   or a fine-grid **output** (the corrected field as the deliverable, which is the deeper ask and
   is blocked on gridded NDVI/ERA5 for the joint models but open for coordinate-only methods).
