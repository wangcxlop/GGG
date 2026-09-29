"""
Diagnostic fields describing whether a fold/method's prediction used a DEM or joint-covariate
model, and (when it used joint covariates) which variable roles were selected vs. actually used.
Factored out because the per-method fold loop in `run_interpolation_benchmark` needs the exact
same computation whether the method succeeded or raised — see the `try`/`catch` call sites below.
"""
function _run_status_role_fields(
    dem_context, joint_context, mode::String, method::String, output_method::String,
)
    uses_dem = dem_context !== nothing && mode == "residual" &&
        method in ("gwr", "mixed_gwr", "mgwr")
    uses_joint = joint_context !== nothing && mode == "residual" &&
        method in ("gwr", "mixed_gwr", "mgwr")
    selected_roles = uses_joint ? join([
        "$group=$(joint_context.roles[group])" for group in joint_context.variables
    ], ";") : ""
    effective_role_map = uses_joint ? joint_effective_roles(
        joint_context, output_method,
    ) : Dict{String,String}()
    effective_roles = uses_joint ? join([
        "$group=$(effective_role_map[group])" for group in joint_context.variables
    ], ";") : ""
    # `mixed_gwr` collapses onto `residual_gwr` whenever this fold's role map has no "global":
    # identical design, identical grid, identical numbers (see `joint_models_coincide`). Recorded
    # as a column so a reader of `run_status.csv` can see the two rows are one model, instead of
    # having to notice that their RMSEs match to the last digit.
    duplicate_of = uses_joint && output_method == "mixed_gwr" &&
        joint_models_coincide(joint_context, "mixed_gwr", "residual_gwr") ? "residual_gwr" : ""
    return (;
        uses_dem, uses_joint, selected_roles, effective_role_map, effective_roles, duplicate_of,
    )
end

"""Non-NaN share of the held-out cells a method could have predicted."""
_prediction_coverage(eligible, prediction) =
    count(eligible) == 0 ? 0.0 :
        count(eligible .& .!isnan.(prediction)) / count(eligible)

"""
`auto`'s `(status, error)` for one fold.

"skipped" rather than "failed" where `auto` was never applicable, so a genuine breakage stays
visible instead of being lost among expected non-runs.

A fold that did choose a contender is held to the same coverage gate as the seven fitted runs,
with the same reason string. Without that the `status` column meant two different things in one
table: `auto` read "success" at the very coverage that makes the model it selected read
"partial", so filtering `run_status.csv` on `status == "success"` kept the `auto` row and dropped
the joint-model rows describing the same fold at the same coverage. Measured on the 2026-09-02
full nested run, that was three rows — balanced_spatial/fold 4 for each product, at coverage
0.9028 against a 0.95 bar.
"""
function _auto_status(
    chose_contender::Bool, coverage::Float64, min_coverage::Float64,
    dem_enabled::Bool, has_selection_split::Bool, failures::Vector{String},
)
    chose_contender && return coverage >= min_coverage ? ("success", "") :
        ("partial", "prediction coverage below minimum")
    dem_enabled && return ("skipped", "auto is not run on the legacy DEM path")
    has_selection_split || return ("skipped", "fold has no inner selection split to choose on")
    return ("failed", "no auto contender scored: $(join(failures, " | "))")
end

"""
One row of `run_status.csv`.

The fold loop assembles this row in three places - a method that predicted, a method that raised,
and `auto` - which differ only in `status`, `error` and `prediction_coverage`. Passing `nothing`
for both contexts gives the all-"not_applicable" provenance `auto` reports: it fits no covariate
model of its own and delegates provenance to whichever method it chose, whose own row carries it.
Keeping that a `nothing` argument rather than a separate literal is what stops the three copies
drifting apart, which is also why `covariate_selection_mode` still reads "not_applicable" there
and `auto` stays out of `covariate_model_status.csv`.
"""
function _benchmark_status_row(
    dem_context, joint_context, joint_inputs, nested_joint::Bool;
    scheme, product, fold, repeat::Int, seed::Int, mode::String, method::String,
    output_method::String, status::String, error::String, prediction_coverage::Float64,
)
    (; uses_dem, uses_joint, selected_roles, effective_roles, duplicate_of) =
        _run_status_role_fields(dem_context, joint_context, mode, method, output_method)
    return (;
        scheme, product, fold, repeat, seed, method=output_method, status, error,
        prediction_coverage,
        dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
        dem_selection_status=uses_dem ? dem_context.selection_status : "not_applicable",
        dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
        dem_roles=uses_dem ? dem_context.dem_roles : "",
        covariate_selection_mode=uses_joint ?
            (nested_joint ? "nested_per_fold" : "fixed_full_data") : "not_applicable",
        covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
        covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
        covariate_selected_roles=selected_roles,
        covariate_effective_roles=effective_roles,
        covariate_spec_sha256=uses_joint ? something(joint_inputs.spec_sha256, "") : "",
        # Appended rather than grouped with the other status fields: existing readers of this
        # table index it by name, but the ad-hoc awk/pandas kind does not, and a new column at
        # the end cannot shift anything that already exists.
        duplicate_of,
    )
end

"""
Provenance and input-QC tables for a joint-covariate run.

`joint_spec_provenance.csv` is written on both paths so a reader can always tell which selection
mode produced the run: nested per-fold selection is confirmatory, a fixed full-data spec is not.
"""
function _write_joint_provenance(cfg::InterpolationBenchmarkConfig, joint_inputs, nested_joint::Bool)
    if nested_joint
        CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
            source_path="", sha256="", selection_mode="nested_per_fold", confirmatory=true,
        )]))
    else
        CSV.write(
            joinpath(cfg.mger.outdir, "joint_variable_spec_used.csv"),
            joint_inputs.specification.table,
        )
        CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
            source_path=abspath(something(cfg.joint_covariates).spec_path),
            sha256=joint_inputs.spec_sha256,
            selection_mode="fixed_full_data",
            confirmatory=false,
        )]))
    end
    CSV.write(joinpath(cfg.mger.outdir, "joint_era5_input_qc.csv"), joint_inputs.era5.qc)
    if joint_inputs.ndvi !== nothing
        CSV.write(joinpath(cfg.mger.outdir, "joint_ndvi_alignment_qc.csv"),
            joint_inputs.ndvi.aligned.qc)
    end
    return nothing
end

"""
Screen the DEM covariates once over every station, before any fold is drawn.

Reported alongside the per-fold screens so a reader can see which variables the full station set
supports; the `"full_data"` scheme label keeps it out of the cross-fold stability tables.
"""
function _screen_dem_full_data!(dem_store, cfg::InterpolationBenchmarkConfig, terrain, lonlat, products, product_data)
    dem = something(cfg.dem)
    for product in products
        data = product_data[product]
        selection = screen_dem_subset(
            terrain, lonlat, data.times, data.Y_obs, data.Y_sat, dem;
            scheme="full_data", product, fold=0, phase="full_data",
            seed=cfg.seed + 100_000 + Int(sum(codeunits(product))),
        )
        _store_dem_selection!(dem_store, selection)
    end
    return dem_store
end

"""
Assemble every result table from the accumulated rows and write the run's CSVs.

Split out of `run_interpolation_benchmark` so the orchestration reads as
setup -> repeat/scheme/product/fold loops -> outputs. Parameters are named for the accumulators
they receive so the body is unchanged from when it was inline.
"""
function _write_benchmark_outputs(
    cfg::InterpolationBenchmarkConfig, products, seeds, nested_joint::Bool, joint_inputs,
    dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
    all_metric_rows, all_scan_rows, all_bootstrap_rows, run_status_rows,
    auto_selection_rows, blend_selection_rows, fused_anchor_rows, hurdle_rows,
    selection_fallback_cells,
)
    metrics = DataFrame(all_metric_rows)
    scans = DataFrame(all_scan_rows)
    bootstrap = DataFrame(all_bootstrap_rows)
    status = DataFrame(run_status_rows)
    claim = if :balanced_spatial in cfg.cv_schemes && cfg.bootstrap_reps > 0
        assess_gwr_claim(metrics, bootstrap, products)
    else
        DataFrame()
    end
    repeat_summary = summarize_repeats(metrics)
    fold_summary = summarize_fold_spread(metrics)
    rank_stability = method_rank_stability(metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_stratified.csv"), metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_folds.csv"), filter(:fold => (x -> !ismissing(x)), metrics))
    CSV.write(joinpath(cfg.mger.outdir, "metrics_pooled.csv"), filter(:fold => ismissing, metrics))
    nrow(repeat_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_repeat_summary.csv"), repeat_summary)
    # Dispersion beside the pooled point estimate, so a headline gap can be read against how much
    # the methods move between folds. See `summarize_fold_spread` for why `_std` is not a standard
    # error.
    nrow(fold_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_fold_summary.csv"), fold_summary)
    nrow(rank_stability) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "method_rank_stability.csv"), rank_stability)
    CSV.write(joinpath(cfg.mger.outdir, "parameter_scan.csv"), scans)
    # How much of each hurdle prediction actually came from a local fit. Without this a
    # globally-degenerate model is indistinguishable from a local one in the metrics.
    isempty(hurdle_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "hurdle_diagnostics.csv"), DataFrame(hurdle_rows))
    nrow(bootstrap) > 0 && CSV.write(joinpath(cfg.mger.outdir, "paired_comparisons.csv"), bootstrap)
    CSV.write(joinpath(cfg.mger.outdir, "run_status.csv"), status)
    # Which method `auto` picked in each fold, and by how much. A `chosen` column that changes
    # from fold to fold is the honest reading of "no single GWR variant is best here", and the
    # gap to `runner_up_rmse` says whether the choice was decisive or a coin flip.
    isempty(auto_selection_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "auto_selection.csv"), DataFrame(auto_selection_rows))
    isempty(blend_selection_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "blend_selection.csv"), DataFrame(blend_selection_rows))
    # The fusion coefficients each derived product's anchor was built from, per fold. Worth
    # reporting for the same reason `auto_selection.csv` is: they are fitted quantities the
    # published numbers rest on, and their stability across folds is the evidence that the anchor
    # is a property of the products rather than of one partition. `fell_back` marks a fold whose
    # normal equations were singular and which therefore used equal weights instead.
    isempty(fused_anchor_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "fused_anchor_selection.csv"),
            DataFrame(fused_anchor_rows))
    _dem_enabled(cfg) && _write_dem_outputs(cfg.mger.outdir, dem_store, scans, status)
    if joint_inputs !== nothing
        scaling = isempty(joint_scaling_tables) ? DataFrame() :
            vcat(joint_scaling_tables...; cols=:union)
        quality = isempty(joint_qc_tables) ? DataFrame() :
            vcat(joint_qc_tables...; cols=:union)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_scaling.csv"), scaling)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_quality_control.csv"), quality)
        bandwidths = filter([:mode, :method, :selected] =>
            (mode, method, selected) -> mode == "residual" &&
                method in ("gwr", "mixed_gwr", "mgwr") && selected, scans)
        CSV.write(joinpath(cfg.mger.outdir, "joint_bandwidths.csv"), bandwidths)
        CSV.write(joinpath(cfg.mger.outdir, "covariate_model_status.csv"), filter(
            :covariate_selection_mode => in(("fixed_full_data", "nested_per_fold")), status,
        ))
        nested_joint && _write_joint_selection_outputs(cfg.mger.outdir, joint_store)
    end
    if ncol(claim) > 0
        CSV.write(joinpath(cfg.mger.outdir, "claim_assessment.csv"), claim)
        # One row per (repeat, product) above; this collapses them so "does the claim hold" can be
        # separated from "did it hold in one partition". Trivial for a single-partition run.
        agreement = claim_agreement(claim)
        nrow(agreement) > 0 &&
            CSV.write(joinpath(cfg.mger.outdir, "claim_agreement.csv"), agreement)
    end
    scope = DataFrame(
        key=[
            "repeated_cv_partitions", "repeated_cv_seeds", "fold_center_init",
            "fold_rotations",
            "validation_target", "training_signal", "temporal_holdout", "primary_cv",
            "secondary_cv", "common_evaluation_mask", "tuning_time_limit",
            "tuning_time_weighting", "tuning_geometry",
            "model_selection", "per_repeat_oof_tables", "results_admissible",
            "dem_variable_selection", "dem_role_assignment", "dem_leakage_control",
            "dem_empty_selection", "supported_claim", "unsupported_claim",
        ],
        value=[
            string(length(seeds)), join(seeds, ","), string(cfg.fold_center_init),
            cfg.fold_center_init === :hilbert ?
                "seed-free: partition i is rotation i-1 of the Hilbert frame (" *
                    join(0:(length(seeds) - 1), ",") * "); rotation 0 is the canonical one" :
                "not applicable: $(cfg.fold_center_init) draws its initial centers from the seed",
            "held-out stations at matched observation/satellite timestamps",
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus fold-selected DEM variables for residual_gwr, mixed_gwr, and mgwr",
            "false", "balanced two-dimensional spatial 5-fold", "random station 5-fold",
            "true across all eight methods", string(cfg.tuning_max_times),
            cfg.tuning_time_weighting === :stratified ?
                "stratified: wettest eighth taken with certainty, remaining hours systematically subsampled and inverse-probability weighted so the tuning RMSE estimates the reported pooled RMSE" :
                "uniform (default): unweighted RMSE over a wet-oversampled subsample, so the tuning RMSE runs ~3x the reported pooled RMSE it stands in for; the level error is close to a constant multiplier and cancels in the candidate ranking",
            cfg.tuning_geometry === :inner_spatial ?
                "inner_spatial (default): candidates scored by out-of-fold prediction onto an inner $(cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k)-group split of the training stations, built with the same splitter and scheme as the outer partition, so selection and reporting are the same estimand" *
                (isempty(selection_fallback_cells) ? "" :
                    "; FELL BACK to leave-one-out in $(length(selection_fallback_cells)) cell(s) too small to split: $(join(selection_fallback_cells, ", "))") :
                "loocv (legacy): candidates scored by leave-one-out at training stations, which measures interpolation next to a retained gauge while the reported metric measures extrapolation to a station 20+ km from any gauge",
            "auto: one GWR-family method chosen per fold on the inner selection split ($(join([_output_method(mode, method) for (mode, method) in AUTO_CANDIDATE_RUNS], ", "))), contenders re-predicted and compared on the intersection of their masks; every method is also reported on its own, but naming a per-method winner from those columns is selection on the held-out fold, which is what auto exists to avoid",
            # `oof_*.csv` and `common_evaluation_mask.csv` are large, so only the first partition
            # writes them. Stated here because the claim path no longer assumes repeat 1 exists.
            length(seeds) == 1 ? "written for the single partition" :
                "written for repeat_01 only; later partitions report metrics but not per-station OOF tables",
            cfg.exploratory_only ?
                "NO - exploratory_only=true: joint covariates and roles were selected once over every station, so the outer held-out fold helped choose the covariates its own predictions use" :
                "yes - every selection step ran inside its training fold",
            _dem_enabled(cfg) ? "training-fold Pearson/Spearman direction check, joint aspect F test, BH q<0.05, and grouped VIF<5" : "disabled",
            _dem_enabled(cfg) ? "training-fold GWR Monte Carlo spatial nonstationarity test with BH q<0.05" : "disabled",
            _dem_enabled(cfg) ? "all DEM screening, scaling, role tests, and bandwidth selection use training stations only" : "not applicable",
            _dem_enabled(cfg) ? "fall back to the original spatial-only residual design and record dem_variable_count=0" : "not applicable",
            "gauge-network-assisted interpolation/correction at stations not used for fitting",
            "temporal forecast or correction without concurrent gauge observations",
        ],
    )
    if joint_inputs !== nothing
        joint = something(cfg.joint_covariates)
        append!(scope, DataFrame(
            key=["mgwr_spatial_grouping", "backfit_relaxation", "backfit_max_iterations",
                "residual_shrinkage"],
            value=[
                joint.mgwr_spatial_grouping === :split ?
                    "split (default): the intercept, longitude and latitude each form their own single-column back-fitting group" :
                    joint.mgwr_spatial_grouping === :shared ?
                        "shared: intercept, longitude and latitude solved as one weighted least squares, so only the covariates carry separate bandwidths" :
                        "intercept_only: the coordinate columns are dropped; a locally varying intercept plus one group per local covariate",
                "$(joint.relaxation) (successive over-relaxation on the back-fitting sweeps; converges to the same fixed point for any admissible value, so this trades rate only — plain Gauss-Seidel at 1.0 fails to converge on most hours)",
                string(joint.max_iterations),
                length(cfg.residual_shrinkage_candidates) == 1 &&
                    only(cfg.residual_shrinkage_candidates) == 1.0 ?
                    "disabled: the residual correction is added back unshrunk, as GWR produces it" :
                    "enabled: a scale in $(minimum(cfg.residual_shrinkage_candidates))-$(maximum(cfg.residual_shrinkage_candidates)) is selected with the bandwidth on the inner spatial split and applied to the joint dynamic models' residual correction",
            ],
        ))
        training_row = findfirst(==("training_signal"), scope.key)
        scope.value[training_row] = nested_joint ?
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a per-fold, training-station-only product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr" :
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a fixed full-data product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr"
        append!(scope, DataFrame(
            key=[
                "covariate_selection_mode", "covariate_selection_leakage",
                "covariate_parameter_leakage", "covariate_temporal_fit",
                "inference_policy",
            ],
            value=nested_joint ? [
                "nested_per_fold",
                "variable and role selection re-run inside every training fold; validation stations never seen by selection",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "eligible for the same paired significance test and claim assessment as the legacy DEM path",
            ] : [
                "fixed_full_data",
                "non-nested variable selection used all 237 stations before spatial cross-validation",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "exploratory performance metrics only; no paired significance or claim assessment",
            ],
        ))
    end
    if cfg.satellite_wet_blend
        append!(scope, DataFrame(
            key=["satellite_wet_blend", "satellite_wet_blend_selection", "blend_axes"],
            value=[
                "blended counterparts of $(join(BLEND_SOURCE_METHODS, ", ")) reported as " *
                "$(join([_blend_method(m, a) for a in cfg.blend_axes
                         for m in BLEND_SOURCE_METHODS], ", ")): where the " *
                "satellite reports at least $(cfg.mger.rain_threshold) mm, the prediction is " *
                "blended toward $(BLEND_FALLBACK_METHOD)",
                "blending weight chosen per fold on the inner selection split over " *
                "$(BLEND_LAMBDAS), never on held-out cells; see blend_selection.csv. The " *
                "blended methods are scored on the shared mask but do not define it, so adding " *
                "them leaves every other method's denominator unchanged",
                "$(join(cfg.blend_axes, ", ")). `constant` is one weight for every " *
                "satellite-wet cell. `agreement_envelope` gives one weight per band of (how " *
                "many source products call the cell wet) x (their maximum, on " *
                "$(BLEND_ENVELOPE_EDGES)), which widens the intervention set: a cell this " *
                "product calls dry but another calls wet is blended, where the constant axis " *
                "would not have touched it",
            ],
        ))
    end
    if cfg.stacked_blend
        append!(scope, DataFrame(
            key=["stacked_blend", "stacked_blend_selection"],
            value=[
                "$(STACK_METHOD): non-negative combination of $(join(STACK_FEATURES, ", ")), " *
                "one weight vector per band of (agreement-envelope band including no-product-wet) " *
                "x (nearest training gauge at $(STACK_DISTANCE_EDGES) km)",
                "weights fitted per fold on the inner selection split over the full training " *
                "record, never on held-out cells; see stack_selection.csv. Scored on the shared " *
                "mask without defining it",
            ],
        ))
    end
    if !isempty(cfg.fused_anchor_variants)
        append!(scope, DataFrame(
            key=["fused_anchor", "fused_anchor_leakage"],
            value=[
                "derived products $(join([_fused_product_name(v) for v in
                    sort(cfg.fused_anchor_variants)], ", ")): the anchor is a combination of the " *
                "products that have files, not a satellite field as shipped. Their `raw` row is " *
                "therefore a gauge-trained predictor that uses only satellites at prediction " *
                "time, and is reported so the fusion's contribution can be read separately from " *
                "the correction built on top of it",
                "fusion coefficients fitted inside every fold on that fold's training stations " *
                "alone and applied to every station, so a held-out gauge never reaches the " *
                "anchor its own prediction is built on; see fused_anchor_selection.csv. A " *
                "derived product's evaluation mask needs every source product finite, so it is " *
                "a different cell population from any of them",
            ],
        ))
    end
    grouped = [v for v in sort(cfg.fused_anchor_variants) if v in GROUPED_FUSION_VARIANTS]
    if !isempty(grouped)
        append!(scope, DataFrame(
            key=["fused_anchor_grouped"],
            value=[
                "$(join([_fused_product_name(v) for v in grouped], ", ")): regressed on every " *
                "source product at t, t-1 and t+1 and on each product's mean over the 8 nearest " *
                "other stations, one coefficient set per agreement-envelope band; features read " *
                "satellite values only. Coefficients per fold and band in " *
                "fused_anchor_grouped_coefficients.csv, not fused_anchor_selection.csv",
            ],
        ))
    end
    append!(scope, _git_provenance())
    CSV.write(joinpath(cfg.mger.outdir, "benchmark_scope.csv"), scope)
    return (; metrics, scans, bootstrap, status, claim, repeat_summary, fold_summary,
        rank_stability, auto_selection=DataFrame(auto_selection_rows))
end

"""
Which code produced this run, as `scope` rows.

`benchmark_scope.csv` recorded the CV scope and nothing about provenance, so a published run
could be dated but not attributed. `git_dirty` is the load-bearing one: a clean tree means
`git_commit` fully determines the code, and a dirty tree means it does not, which is the
difference between a citable baseline and a plausible one.

Every call is wrapped. A missing git, a checkout that is not a repository, or an ownership check
that refuses records `"unavailable"` rather than aborting a run that takes hours.
"""
function _git_provenance()
    # Two levels up: this file lives in `src/benchmark/`. Derived rather than configured —
    # `src/` must not read fixed absolute paths, and the process working directory is whatever
    # the calling script was launched from.
    root = normpath(joinpath(@__DIR__, "..", ".."))
    run_git(args::Vector{String}) = try
        String(strip(read(pipeline(`git -C $root $args`; stderr=devnull), String)))
    catch
        nothing
    end
    commit = run_git(["rev-parse", "HEAD"])
    branch = run_git(["rev-parse", "--abbrev-ref", "HEAD"])
    porcelain = run_git(["status", "--porcelain"])
    return DataFrame(
        key=["git_commit", "git_branch", "git_dirty"],
        value=[
            something(commit, "unavailable"),
            something(branch, "unavailable"),
            porcelain === nothing ? "unavailable" : string(!isempty(porcelain)),
        ],
    )
end
