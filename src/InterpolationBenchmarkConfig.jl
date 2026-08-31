const BENCHMARK_METHODS = [
    "raw", "zero", "train_clim", "hour_field_mean",
    "idw", "adw", "tps", "gwr",
    "residual_gwr", "mixed_gwr", "mgwr", "auto",
]
const TRADITIONAL_METHODS = ["idw", "adw", "tps"]
# Hyperparameter-free reference predictors, estimated from training stations only. They are
# reported so a method's RMSE can be read as skill rather than as a bare number: on hourly
# precipitation ~93% of station-hours are dry, so a large part of any RMSE is simply how hard
# a method shrinks toward zero. Deliberately not in TRADITIONAL_METHODS - they are context for
# the reader, not baselines the GWR claim is assessed against.
const NULL_METHODS = ["zero", "train_clim", "hour_field_mean"]
# The common evaluation mask stays pinned to the original comparison set. Every method must be
# finite for a cell to count, so letting a newly added method into the mask would silently
# re-score all the others and break comparability with earlier runs. This used to be a real
# confound - `mgwr`'s back-fit non-convergence dropped whole hours and dragged the shared mask
# down to ~0.76 - but that was a fit-level bug in mgwr's back-fit (fixed by raising
# max_iterations and defaulting mgwr_spatial_grouping to :intercept_only), not a `gwr` coverage
# problem: direct `gwr` has always had ~1.0 own coverage. Post-fix, the shared mask sits at
# ~0.96-0.98 and mgwr's remaining dropouts cost <1e-4 relative RMSE (see
# scripts/run_mgwr_diagnostics.jl). Diagnostic and reference methods are therefore scored *on*
# the mask without being allowed to *define* it.
const MASK_METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]
# `hurdle_gwr` is deliberately absent. The model and its benchmark adapters
# (`build_hurdle_context`, `_hurdle_predict`, and the branches in
# `select_interpolation_parameter!` / `predict_selected`) are kept so re-enabling it is a
# one-line change here, but it is not part of the reported comparison.
const BENCHMARK_RUNS = [
    ("direct", "idw"), ("direct", "adw"), ("direct", "tps"), ("direct", "gwr"),
    ("residual", "gwr"), ("residual", "mixed_gwr"), ("residual", "mgwr"),
]

"""Method name a `(mode, method)` run reports under."""
_output_method(mode::AbstractString, method::AbstractString) =
    method in ("mixed_gwr", "mgwr") ? String(method) :
        (mode == "direct" ? String(method) : "residual_$(method)")

const AUTO_METHOD = "auto"
# The runs `auto` may choose between: the GWR family only.
#
# `auto` exists because nothing else in the benchmark selects a *model*. Every method is fitted
# and reported, which is honest right up until a conclusion names a winner — at which point the
# winner was picked by reading held-out scores, the same selection-on-test that
# `run_multikernel_spatial_kfold_pipeline`'s own docstring identifies on the MGER side. `auto`
# makes that choice inside the training fold instead, so its held-out numbers estimate "use the
# best GWR-family model" rather than "use whichever one won once we had looked".
#
# The traditional baselines are deliberately excluded. They are what the GWR claim is assessed
# *against*, so an `auto` free to pick `adw` would leave `assess_gwr_claim`'s "beat every
# baseline" gate comparing `adw` with itself.
const AUTO_CANDIDATE_RUNS = [
    ("direct", "gwr"), ("residual", "gwr"), ("residual", "mixed_gwr"), ("residual", "mgwr"),
]

"""Convert the benchmark's `Int` kernel index into the weight-formula `Function` the
`mixed_gwr`/`mgwr` local-hat machinery (`DEMTerrainExperiment`/`JointCovariateModels`) expects.
`GWR_KERNELS` isn't exported, so it needs `MixedGWR.` qualification here."""
_kernel_function(kernel::Int) = MixedGWR.GWR_KERNELS[kernel + 1]

Base.@kwdef struct InterpolationBenchmarkConfig
    mger::MGERConfig
    terrain_path::Union{Nothing,String} = nothing
    dem::Union{Nothing,DEMExperimentConfig} = nothing
    joint_covariates::Union{Nothing,JointCovariateBenchmarkConfig} = nothing
    # Set to run joint-covariate variable selection fresh inside every training fold instead of
    # once on the full station set (`joint_covariates.spec_path`). Exactly one of the two must
    # be set when `joint_covariates` is enabled.
    joint_selection::Union{Nothing,JointSelectionConfig} = nothing
    # Acknowledges that this run's numbers are not admissible as a result.
    #
    # `joint_covariates.spec_path` names a variable set and local/global role map screened over
    # every station (`JointVariableSelection`'s `"full_data"` scheme), which is then reused inside
    # every training fold - so the outer held-out stations helped choose the covariates their own
    # predictions are built from. The run already labels the mode and zeroes `bootstrap_reps`, but
    # labelling is not a guard: `metrics_pooled.csv` still comes out looking like every other run's.
    #
    # `_validate_benchmark_config` therefore refuses `spec_path` unless this is set, so the leaky
    # path stays reachable for exploration and cannot be reached by accident.
    exploratory_only::Bool = false
    k::Int = 5
    seed::Int = 20260627
    # Repeated cross-validation: one independent fold partition per seed. Empty means a
    # single repeat using `seed`, which reproduces the historical single-split behaviour.
    seeds::Vector{Int} = Int[]
    # How `balanced_spatial` picks its initial fold centers. `:hilbert` is seed-free: the stations
    # are ordered along a Hilbert curve of the local-km frame, cut into `k` equal runs, and the run
    # centroids seed the capacity-constrained Lloyd loop. Repeat `i` uses rotation `i-1` of that
    # frame, so repeated cross-validation is indexed by a declared rotation rather than by an
    # arbitrary RNG draw, and rotation 0 is the canonical partition. `:kmeanspp` and `:farthest`
    # remain reachable and still read `seed`, so earlier runs reproduce exactly.
    fold_center_init::Symbol = :hilbert
    cv_schemes::Vector{Symbol} = [:balanced_spatial, :random]
    idw_powers::Vector{Float64} = [1.0, 1.5, 2.0, 2.5, 3.0]
    neighbor_candidates::Vector{Union{Nothing,Int}} = Union{Nothing,Int}[8, 16, 32, 64, nothing]
    tps_smooth_candidates::Vector{Float64} = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    # Gates two different quantities under one threshold. In `_select_candidate!` it is the
    # share of *inner-tuning* cells a candidate scored, and rejects the candidate; in
    # `_benchmark_status_row` the same number is compared against a fold's *outer prediction*
    # coverage to mark a run "partial". Those are different estimands over different samples and
    # there is no reason they should share a value - separating them is a config change with a
    # numbers impact, so it is recorded rather than done.
    min_tuning_coverage::Float64 = 0.95
    tuning_max_times::Int = 336
    # `:stratified` weights the tuning subsample so its RMSE estimates the pooled RMSE the
    # benchmark reports, instead of the wet-hour RMSE the unweighted subsample scores.
    #
    # It is NOT the default, despite being the more correct estimator, because measurement says
    # the level error it removes was very nearly a constant multiplier and therefore cancelled
    # in the ranking: over 30 product/fold/method cells it costs 0.586% mean regret against the
    # exact full-record curve where `:uniform` costs 0.096%, and on the joint-covariate path the
    # two are a wash (63.4% vs 63.7%). Switching the default would move published numbers for no
    # measured gain. See `scripts/verify_tuning_time_weighting.jl`, and the header of
    # `_tuning_time_sample` for what actually dominates selection error.
    tuning_time_weighting::Symbol = :uniform
    # How candidates are scored during selection. `:inner_spatial` predicts out-of-fold onto an
    # inner split of the training stations, built by the same splitter and scheme as the outer
    # partition, so the selection criterion is the same estimand the benchmark reports.
    # `:loocv` is the historical leave-one-out-at-training-stations criterion, which scores
    # interpolation skill (median 5.5 km to the nearest remaining gauge) while the results table
    # reports extrapolation skill (median 23.9 km under `balanced_spatial`). See
    # `selection_folds` for the measurements.
    tuning_geometry::Symbol = :inner_spatial
    # Groups in the inner selection split; 0 means "use `cfg.k`", which is what reproduces the
    # outer fold geometry. Ignored when `tuning_geometry === :loocv`.
    tuning_inner_k::Int = 0
    mgwr_max_tuning_iterations::Int = 5
    # Adds an explicit global (unweighted) candidate to the GWR family's bandwidth search, so
    # "this coefficient wants to be global" is a value the grid can express rather than
    # something that only happens by accident when an adaptive bandwidth overruns the training
    # set. See `_bandwidth_families`. `false` reproduces the historical purely-local search.
    bw_include_global::Bool = true
    # Scale applied to the joint dynamic models' residual correction before it is added back to
    # the satellite field. GWR is unbiased and never shrinks, so the correction it produces is
    # the right shape at the wrong magnitude: `satellite_offset.csv` measures the optimal rescale
    # at 0.84 / 0.47 / 0.53 for FY4B / GPM / GSMaP, i.e. the models over-correct by up to 2x.
    #
    # The scale is selected jointly with the bandwidth on the inner spatial split, and costs
    # nothing to add: `_joint_candidate_metrics` rescores the residual prediction it has already
    # computed, so the ladder needs no extra model fits. Scoring the clipped objective rather than
    # the closed-form `E[g*r]/E[g^2]` keeps selection consistent with the reported metric, which
    # the `max(., 0)` clip would otherwise break.
    #
    # The ladder starts at 0.1 rather than 0.3 because a smoke run with a 0.3 floor put 31 of 75
    # selected joint rows *on* that floor - the grid, not the data, would have been picking the
    # scale, which is the same failure the bandwidth grids hit before `--local-grid`.
    #
    # `[1.0]` disables shrinkage and reproduces the historical unshrunk behaviour exactly.
    residual_shrinkage_candidates::Vector{Float64} =
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    event_thresholds::Vector{Float64} = [0.1, 2.5, 8.0, 16.0]
    bootstrap_reps::Int = 2000
end

"""
Stamp identifying columns onto the front of a diagnostic table, in the order given.

The DEM and joint screening paths both write per-fold diagnostic tables that have to be traceable
back to the cell that produced them; they differ only in which columns that takes (`phase` on the
DEM side, `seed` on the joint side, and the two order `repeat`/`fold` oppositely). Column order is
the argument order, and it is part of the CSV contract, so it stays with the caller.

`keep_existing` names columns the caller is only supplying a default for: `product` is already a
column of some screening outputs, and there the table's own value wins.
"""
function _tag_selection_table!(
    table::DataFrame, columns::Pair{Symbol}...; keep_existing=(:product,),
)
    for (column, value) in columns
        column in keep_existing && column in propertynames(table) && continue
        table[!, column] = fill(value, nrow(table))
    end
    ordered = [column for (column, _) in columns]
    select!(table, ordered..., Not(ordered))
    return table
end

"""Fold seeds for the run: one repeat per seed, defaulting to a single repeat on `cfg.seed`."""
benchmark_seeds(cfg::InterpolationBenchmarkConfig) =
    isempty(cfg.seeds) ? [cfg.seed] : cfg.seeds

"""
Bandwidth families every GWR-family method searches: adaptive neighbour counts, fixed km, and -
when `cfg.bw_include_global` is set - an explicit global candidate. `adaptive_candidates` is
passed in because the joint path filters its own grid against the training size first.

Global is `(adaptive=false, bw=Inf)`. Every kernel in `GWR_KERNELS` returns exactly 1.0 at
`bw=Inf` (see `src/kernel.jl`), so the fit is unweighted OLS over the training stations.
Expressing it as a bandwidth rather than as its own code path means nothing downstream -
`predict_selected`, the joint back-fit, the `_scan_row` schema - needs to know about it, and
`gw_weight`'s `bw > n` branch is never reached during selection or at the final fit.

Routing every method through one helper is also what keeps the family comparable: direct `gwr`
used to search `cfg.mger.bw_adaptive` while the joint models searched a wider grid of their own,
so a direct-vs-residual comparison was not at equal search budget. The legacy `--legacy-dem`
path (`InterpolationBenchmarkDEM.jl`) is the one exception and keeps its own grid, because it is
mutually exclusive with the joint path and is not part of the reported comparison.
"""
_bandwidth_families(cfg::InterpolationBenchmarkConfig, adaptive_candidates) = (
    (true, collect(Float64, adaptive_candidates)),
    (false, cfg.bw_include_global ? vcat(cfg.mger.bw_fixed_km, Inf) : cfg.mger.bw_fixed_km),
)

"""
Smallest training-set size a candidate is actually fitted on during selection.

Candidates are scored on the inner selection split, which holds out a whole spatial group, so
the fit behind a scan row sees fewer stations than the outer training set the winner is later
refitted on. An adaptive bandwidth above this count falls through `gw_weight`'s `dn > 1` branch
to a near-global fit, which would make the scanned model and the refitted model different
things under the same `bw`. `:loocv` geometry holds out one station instead of a group.
"""
_selection_train_minimum(n_train::Int, selection_groups) =
    selection_groups === nothing ? n_train - 1 :
        n_train - maximum(length(group) for group in selection_groups)

"""Adaptive candidates that both the inner selection fit and the final fit can honour."""
_usable_adaptive_candidates(candidates, n_train::Int, selection_groups) =
    filter(<=(_selection_train_minimum(n_train, selection_groups)), candidates)

function _validate_benchmark_config(cfg::InterpolationBenchmarkConfig, n_station::Int)
    2 <= cfg.k <= n_station || throw(ArgumentError("k must be between 2 and station count"))
    length(unique(benchmark_seeds(cfg))) == length(benchmark_seeds(cfg)) ||
        throw(ArgumentError("seeds must be unique"))
    cfg.fold_center_init in (:hilbert, :kmeanspp, :farthest) ||
        throw(ArgumentError("fold_center_init must be :hilbert, :kmeanspp, or :farthest"))
    all(s -> s in (:balanced_spatial, :random, :strip), cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must contain :balanced_spatial, :random, or :strip"))
    length(unique(cfg.cv_schemes)) == length(cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must be unique"))
    all(>(0), cfg.idw_powers) || throw(ArgumentError("IDW/ADW powers must be positive"))
    all(x -> x === nothing || x > 0, cfg.neighbor_candidates) ||
        throw(ArgumentError("neighbor candidates must be positive or nothing"))
    all(>(0), cfg.tps_smooth_candidates) ||
        throw(ArgumentError("TPS smoothing candidates must be positive"))
    0 < cfg.min_tuning_coverage <= 1 ||
        throw(ArgumentError("min_tuning_coverage must be in (0, 1]"))
    cfg.bootstrap_reps >= 0 || throw(ArgumentError("bootstrap_reps must be non-negative"))
    cfg.tuning_max_times >= 0 || throw(ArgumentError("tuning_max_times must be non-negative"))
    cfg.tuning_time_weighting in (:stratified, :uniform) ||
        throw(ArgumentError("tuning_time_weighting must be :stratified or :uniform"))
    cfg.tuning_geometry in (:inner_spatial, :loocv) ||
        throw(ArgumentError("tuning_geometry must be :inner_spatial or :loocv"))
    cfg.tuning_inner_k == 0 || cfg.tuning_inner_k >= 2 ||
        throw(ArgumentError("tuning_inner_k must be 0 (use cfg.k) or at least 2"))
    cfg.mgwr_max_tuning_iterations > 0 ||
        throw(ArgumentError("mgwr_max_tuning_iterations must be positive"))
    isempty(cfg.residual_shrinkage_candidates) &&
        throw(ArgumentError("residual_shrinkage_candidates must not be empty"))
    all(s -> 0 < s <= 1, cfg.residual_shrinkage_candidates) ||
        throw(ArgumentError("residual shrinkage candidates must be in (0, 1]"))
    xor(cfg.terrain_path === nothing, cfg.dem === nothing) && throw(ArgumentError(
        "terrain_path and dem must either both be configured or both be omitted",
    ))
    cfg.joint_covariates !== nothing && cfg.terrain_path !== nothing && throw(ArgumentError(
        "joint covariates and legacy DEM screening cannot be enabled together",
    ))
    if cfg.dem !== nothing
        cfg.dem.min_wet_hours > 0 || throw(ArgumentError("DEM min_wet_hours must be positive"))
        cfg.dem.screen_permutations > 0 ||
            throw(ArgumentError("DEM screen_permutations must be positive"))
        cfg.dem.spatial_permutations > 0 ||
            throw(ArgumentError("DEM spatial_permutations must be positive"))
        all(>(1), cfg.dem.bandwidth_candidates) ||
            throw(ArgumentError("DEM bandwidth candidates must exceed one neighbor"))
    end
    if cfg.joint_covariates !== nothing
        joint = cfg.joint_covariates
        xor(joint.spec_path === nothing, cfg.joint_selection === nothing) || throw(ArgumentError(
            "joint covariates require exactly one of spec_path (fixed full-data selection) or " *
            "joint_selection (nested per-fold selection)",
        ))
        joint.spec_path === nothing || cfg.exploratory_only || throw(ArgumentError(
            "joint covariates selected once on the full station set (spec_path) let the outer " *
            "held-out fold influence which covariates and local/global roles its own predictions " *
            "use. Use joint_selection for nested per-fold selection, or set exploratory_only=true " *
            "to acknowledge that this run's metrics are not a reportable result",
        ))
        joint.spec_path === nothing || isfile(joint.spec_path) ||
            throw(ArgumentError("joint specification does not exist"))
        isfile(joint.terrain_path) || throw(ArgumentError("joint terrain table does not exist"))
        all(>(1), joint.bandwidth_candidates) ||
            throw(ArgumentError("joint bandwidth candidates must exceed one neighbor"))
        joint.max_iterations > 0 || throw(ArgumentError("joint max_iterations must be positive"))
        joint.tolerance > 0 || throw(ArgumentError("joint tolerance must be positive"))
        # Successive over-relaxation diverges outside (0, 2); see the field comment in
        # `JointCovariateBenchmarkConfig` for why the default sits above 1.
        0 < joint.relaxation < 2 ||
            throw(ArgumentError("joint relaxation must be in (0, 2)"))
        joint.mgwr_spatial_grouping in (:split, :shared, :intercept_only) || throw(ArgumentError(
            "mgwr_spatial_grouping must be :split, :shared, or :intercept_only",
        ))
        joint.unsupported_local_target in (:missing, :zero) || throw(ArgumentError(
            "unsupported_local_target must be :missing or :zero",
        ))
    end
    return nothing
end

"""
Per-repeat output root. A single-repeat run keeps the historical `outdir/<scheme>/...` layout so
existing verify scripts and downstream readers are unaffected; repeated runs nest under
`outdir/repeat_<i>/`.
"""
_repeat_dir(outdir::String, repeat_index::Int, n_repeats::Int) =
    n_repeats == 1 ? outdir : joinpath(outdir, "repeat_$(lpad(repeat_index, 2, '0'))")

_dem_enabled(cfg::InterpolationBenchmarkConfig) = cfg.terrain_path !== nothing
_joint_enabled(cfg::InterpolationBenchmarkConfig) = cfg.joint_covariates !== nothing
