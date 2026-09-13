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

"""
`auto` for one fold: choose a single GWR-family method using only the inner selection split, and
record the choice.

Split out of `_run_benchmark_fold!`, which was 262 lines. Mutates `fold_predictions` (writing the
`auto` entry, NaN when no contender could be scored), `auto_selection_rows` and `run_status_rows`.
Parameters keep the names the enclosing locals had, so the body is unchanged from when it was
inline.

Returns the inner-split contenders and the tuning hours they were scored over, which
`_run_fold_blend!` reuses.
"""
function _run_fold_auto!(
    cfg::InterpolationBenchmarkConfig, fold::Int, scheme, product, repeat_index::Int,
    repeat_seed::Int, val_idx, y_obs, y_sat_val, y_obs_train, y_sat_train, train_lonlat,
    selection_groups, joint_selection_contexts, joint_inputs, nested_joint::Bool,
    fold_selected, fold_predictions, predictions, auto_selection_rows, run_status_rows,
)

    # `auto`: choose one GWR-family method for this fold on the inner split alone.
    #
    # Contenders are re-predicted rather than compared on their scan rows, because the
    # scan scores each method over its own non-NaN cells and the winner would otherwise
    # partly be whoever failed on the hardest ones. Same tuning hours as the scan, so
    # the choice is made against the criterion the candidates were tuned on.
    auto_time_indices, auto_time_weights =
        cfg.tuning_max_times > 0 && size(y_obs_train, 2) > cfg.tuning_max_times ?
            _tuning_time_sample(
                y_obs_train, cfg.tuning_max_times, cfg.tuning_time_weighting,
            ) : (collect(axes(y_obs_train, 2)), nothing)
    auto_contenders = NamedTuple[]
    auto_failures = String[]
    # Not run on the legacy DEM path: `inner_selection_prediction` deliberately carries
    # no `dem_context`, and that path is mutually exclusive with the joint one and is
    # not part of the reported comparison.
    auto_applicable = !_dem_enabled(cfg) && selection_groups !== nothing
    if auto_applicable
        for (mode, method) in AUTO_CANDIDATE_RUNS
            output_method = _output_method(mode, method)
            haskey(fold_selected, output_method) || continue
            try
                push!(auto_contenders, (; method=output_method,
                    prediction=inner_selection_prediction(
                        fold_selected[output_method], method, mode, train_lonlat,
                        y_obs_train, y_sat_train, selection_groups;
                        joint_contexts=joint_selection_contexts,
                        time_indices=auto_time_indices,
                    )))
            catch e
                push!(auto_failures, "$output_method: $(sprint(showerror, e))")
            end
        end
    end
    auto_choice = select_auto_method(
        auto_contenders,
        Matrix{Float64}(y_obs_train[:, auto_time_indices]),
        Matrix{Float64}(y_sat_train[:, auto_time_indices]);
        time_weights=auto_time_weights,
    )
    if auto_choice !== nothing && haskey(fold_predictions, auto_choice.chosen)
        fold_predictions[AUTO_METHOD] = fold_predictions[auto_choice.chosen]
        predictions[AUTO_METHOD][val_idx, :] = fold_predictions[AUTO_METHOD]
        push!(auto_selection_rows, merge(
            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed),
            auto_choice,
            (; skipped=join(auto_failures, " | ")),
        ))
    else
        # No inner split (leave-one-out geometry, or a fold too small to split), or no
        # contender survived. Left unpredicted rather than defaulted to a method,
        # which would make `auto` mean something different in different folds.
        fold_predictions[AUTO_METHOD] = fill(NaN, length(val_idx), size(y_obs, 2))
    end
    # Measured the same way as every other method's row - non-NaN share of the
    # held-out cells it could have predicted - so the column means one thing across the
    # table. `auto_selection.csv` carries the inner-split mask coverage separately.
    auto_eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
    auto_coverage = _prediction_coverage(auto_eligible, fold_predictions[AUTO_METHOD])
    auto_status, auto_error = _auto_status(
        auto_choice !== nothing, auto_coverage, cfg.min_tuning_coverage,
        _dem_enabled(cfg), selection_groups !== nothing, auto_failures,
    )
    push!(run_status_rows, _benchmark_status_row(
        nothing, nothing, joint_inputs, nested_joint;
        scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
        mode="", method=AUTO_METHOD, output_method=AUTO_METHOD,
        status=auto_status, error=auto_error, prediction_coverage=auto_coverage,
    ))
    # Handed back rather than discarded: `_run_fold_blend!` needs the same inner-split predictions
    # for the same methods over the same tuning hours, and re-deriving them would double the most
    # expensive thing this function does.
    return (; contenders=auto_contenders, time_indices=auto_time_indices,
        time_weights=auto_time_weights, applicable=auto_applicable)
end

"""
One `blend_selection.csv` row shape for both weight rules.

The constant axis chooses a scalar and the banded one a vector, and a table whose column count
depends on the axis would be a nuisance to read and impossible to concatenate across runs. The
weights therefore go into one `lambdas` string with `n_bands` beside it, and `lambda` is the single
weight for the constant axis while carrying the first band's weight for the others.
"""
function _blend_selection_weights(choice, axis::Symbol)
    weights = hasproperty(choice, :lambdas) ? choice.lambdas : [choice.lambda]
    return (;
        lambda=choice.lambda, lambdas=join(weights, "|"), n_bands=length(weights),
        inner_RMSE=choice.inner_RMSE, inner_MAE=choice.inner_MAE,
        unblended_RMSE=choice.unblended_RMSE, n=choice.n, coverage=choice.coverage,
        wet_cells=choice.wet_cells,
    )
end

"""
Blended counterparts of the anchored methods: on satellite-wet cells, blend toward `adw`.

Runs after `_run_fold_auto!` and reuses its inner-split contenders, so the only fit this adds is
the fallback's - `adw` is deliberately not an `auto` candidate, so its inner prediction is the one
thing not already on hand. Everything else here is arithmetic over matrices that exist.

The weight is chosen on the inner selection split, never on the held-out cells. A replay that
sweeps the weight over the test set bounds what is reachable; a method has to pick one without
looking and will therefore score worse than that bound.

A fold with no inner split, or with no fallback prediction, leaves the blended methods NaN rather
than falling back to the unblended source. Silently reporting the source under a blended name
would make the method mean different things in different folds - the same reason `auto` leaves
itself unpredicted rather than defaulting.

Mutates `predictions`, `fold_predictions`, `blend_selection_rows` and `run_status_rows`.
"""
function _run_fold_blend!(
    cfg::InterpolationBenchmarkConfig, fold::Int, scheme, product, repeat_index::Int,
    repeat_seed::Int, val_idx, y_obs, y_sat_val, y_obs_train, y_sat_train, train_lonlat,
    selection_groups, auto_inner, fold_selected, fold_predictions, predictions,
    blend_band_train, blend_band_val, blend_selection_rows, run_status_rows,
)
    threshold = cfg.mger.rain_threshold
    eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)

    # One `run_status.csv` row for a blended method, which carries no covariate model.
    status_row(output_method, status, message, coverage) = _benchmark_status_row(
        nothing, nothing, nothing, false;
        scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
        mode="", method=output_method, output_method,
        status, error=message, prediction_coverage=coverage,
    )

    # Leave every blended method unpredicted, with one status row each saying why.
    function skip_all(status::String, message::String)
        for axis in cfg.blend_axes, source in BLEND_SOURCE_METHODS
            output_method = _blend_method(source, axis)
            fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
            push!(run_status_rows, status_row(output_method, status, message, 0.0))
        end
        return nothing
    end

    auto_inner.applicable || return skip_all(
        "skipped", "fold has no inner selection split to choose a blending weight on",
    )
    haskey(fold_selected, BLEND_FALLBACK_METHOD) || return skip_all(
        "skipped", "$BLEND_FALLBACK_METHOD did not fit on this fold, so there is nothing to " *
        "blend toward",
    )

    # The one fit this adds. `adw` is not an `auto` candidate - the traditional baselines are what
    # the GWR claim is assessed against - so its inner-split prediction is not already computed.
    fallback_inner = try
        inner_selection_prediction(
            fold_selected[BLEND_FALLBACK_METHOD], BLEND_FALLBACK_METHOD, "direct",
            train_lonlat, y_obs_train, y_sat_train, selection_groups;
            time_indices=auto_inner.time_indices,
        )
    catch e
        return skip_all("failed",
            "$BLEND_FALLBACK_METHOD inner prediction failed: $(sprint(showerror, e))")
    end

    inner_obs = Matrix{Float64}(y_obs_train[:, auto_inner.time_indices])
    inner_sat = Matrix{Float64}(y_sat_train[:, auto_inner.time_indices])

    for axis in cfg.blend_axes, source in BLEND_SOURCE_METHODS
        output_method = _blend_method(source, axis)
        contender = findfirst(entry -> entry.method == source, auto_inner.contenders)
        if contender === nothing || !haskey(fold_predictions, source) ||
                !haskey(fold_predictions, BLEND_FALLBACK_METHOD)
            fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
            push!(run_status_rows, status_row(output_method, "skipped",
                "$source has no inner-split prediction to choose a weight on", 0.0))
            continue
        end
        # The constant axis keeps its scalar code path rather than being expressed as a one-band
        # case of the banded one. Both give the same answer, and only one of them cannot regress
        # the constant-weight method by accident.
        inner_prediction = auto_inner.contenders[contender].prediction
        choice = if axis === :constant
            select_blend_lambda(
                inner_obs, inner_sat, inner_prediction, fallback_inner,
                BLEND_LAMBDAS, threshold; time_weights=auto_inner.time_weights,
            )
        else
            select_blend_lambdas(
                inner_obs, inner_sat, blend_band_train[axis][:, auto_inner.time_indices],
                inner_prediction, fallback_inner, BLEND_LAMBDAS, blend_band_count(axis);
                time_weights=auto_inner.time_weights,
            )
        end
        if choice === nothing
            fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
            push!(run_status_rows, status_row(output_method, "failed",
                "no inner-split cell was scorable for $source and $BLEND_FALLBACK_METHOD", 0.0))
            continue
        end
        blended = if axis === :constant
            satellite_wet_blend_prediction(
                y_sat_val, fold_predictions[source], fold_predictions[BLEND_FALLBACK_METHOD],
                choice.lambda, threshold,
            )
        else
            banded_blend_prediction(
                blend_band_val[axis], fold_predictions[source],
                fold_predictions[BLEND_FALLBACK_METHOD], choice.lambdas,
            )
        end
        fold_predictions[output_method] = blended
        predictions[output_method][val_idx, :] = blended
        push!(blend_selection_rows, merge(
            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                method=output_method, source, fallback=BLEND_FALLBACK_METHOD, threshold,
                axis=String(axis)),
            _blend_selection_weights(choice, axis),
        ))
        coverage = _prediction_coverage(eligible, blended)
        push!(run_status_rows, status_row(output_method,
            coverage >= cfg.min_tuning_coverage ? "success" : "partial",
            coverage >= cfg.min_tuning_coverage ? "" : "prediction coverage below minimum",
            coverage))
    end
    return nothing
end

"""
Everything one cross-validation fold does: split the stations, build the fold's DEM/joint/hurdle
contexts, tune and predict each `BENCHMARK_RUNS` method, run `auto` and the blends, and append the
fold's metric, scan and status rows.

Extracted verbatim from `run_interpolation_benchmark`, whose body was a single 580-line function
with the fold loop nested four deep. Parameters are named for the values they receive so the body
is unchanged from when it was inline; the length of this argument list is the honest measure of
how much per-run state a fold touches, and is a fair target for a later pass.

Mutates `predictions`, `nearest_train_distance`, the DEM/joint stores and the row accumulators.
"""
function _run_benchmark_fold!(
    cfg::InterpolationBenchmarkConfig, fold::Int, folds, id_map, ids, products, product,
    data, lonlat, y_obs, y_sat, terrain, joint_inputs, nested_joint::Bool,
    scheme, scheme_symbol, repeat_index::Int, repeat_seed::Int,
    predictions, nearest_train_distance, blend_bands,
    dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
    all_metric_rows, all_scan_rows, run_status_rows, auto_selection_rows,
    blend_selection_rows, hurdle_rows, selection_fallback_cells,
)

    train_ids, val_ids, train_idx, val_idx = _fold_station_indices(folds, id_map, cfg.k, fold)
    train_lonlat = Matrix{Float64}(lonlat[train_idx, :])
    val_lonlat = Matrix{Float64}(lonlat[val_idx, :])
    y_obs_train = Matrix{Float64}(y_obs[train_idx, :])
    y_sat_train = Matrix{Float64}(y_sat[train_idx, :])
    y_sat_val = Matrix{Float64}(y_sat[val_idx, :])
    distance_train_val = haversine_distance_matrix(train_lonlat, val_lonlat)
    nearest_train_distance[val_idx] = vec(minimum(distance_train_val, dims=1))
    null_predictions = _null_fold_predictions(y_obs_train, length(val_idx))
    for (null_method, null_prediction) in null_predictions
        predictions[null_method][val_idx, :] = null_prediction
    end
    # One inner split for the whole fold, so every method's candidates are scored
    # against the same held-out stations, and so the joint contexts below agree with
    # what the spatial-only methods use.
    selection_groups = cfg.tuning_geometry === :inner_spatial ?
        selection_folds(cfg, scheme_symbol, train_ids, train_lonlat, fold,
            repeat_seed; repeat_index) :
        nothing
    if cfg.tuning_geometry === :inner_spatial && selection_groups === nothing
        push!(selection_fallback_cells, "$scheme/$product/fold$fold")
    end

    dem_context = nothing
    if _dem_enabled(cfg)
        dem = something(cfg.dem)
        selection = screen_dem_subset(
            terrain[train_idx, :], train_lonlat, data.times,
            y_obs_train, y_sat_train, dem; scheme, product, fold, phase="cv",
            seed=repeat_seed + 10_000 * fold + Int(sum(codeunits(scheme * product))),
            repeat=repeat_index,
        )
        _store_dem_selection!(dem_store, selection)
        dem_context = build_dem_fold_context(
            selection, terrain[train_idx, :], terrain[val_idx, :],
            train_lonlat, val_lonlat, dem,
        )
    end
    joint_context = nothing
    joint_selection_contexts = nothing
    if joint_inputs !== nothing
        joint = something(cfg.joint_covariates)
        role_map = if nested_joint
            product_index = findfirst(==(product), products)
            run_seed = repeat_seed + 1000 * product_index + 10 * fold
            dem_seed = repeat_seed + 20260815 + Int(sum(codeunits(product))) +
                10_000 + 10 * fold
            joint_selection = _screen_joint_subset(
                # `select_joint_covariates` expects population-wide matrices and
                # slices them via `train_idx` itself (unlike the DEM path's
                # pre-sliced `terrain[train_idx,:]` convention above).
                product, y_obs, y_sat, joint_inputs.terrain,
                joint_inputs.era5.values,
                joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
                train_idx, ids, data.times, lonlat, something(cfg.joint_selection);
                scheme, fold, repeat=repeat_index, seed=repeat_seed, run_seed, dem_seed,
            )
            _store_joint_selection!(joint_store, joint_selection)
            joint_selection.role_map
        else
            joint_inputs.specification.role_maps[product]
        end
        joint_context = build_joint_fold_context(
            product, role_map,
            train_idx, val_idx, lonlat, y_obs, y_sat,
            joint_inputs.terrain, joint_inputs.era5.values,
            joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
            joint,
        )
        # `repeat`/`seed` alongside scheme/fold: without them a repeated run's rows are
        # indistinguishable on disk, since (scheme, fold, product, variable_group)
        # repeats once per partition.
        scaling = copy(joint_context.scaling)
        insertcols!(scaling, 1,
            :scheme => fill(scheme, nrow(scaling)),
            :repeat => fill(repeat_index, nrow(scaling)),
            :seed => fill(repeat_seed, nrow(scaling)),
            :fold => fill(fold, nrow(scaling)))
        push!(joint_scaling_tables, scaling)
        quality = copy(joint_context.quality_control)
        insertcols!(quality, 1,
            :scheme => fill(scheme, nrow(quality)),
            :repeat => fill(repeat_index, nrow(quality)),
            :seed => fill(repeat_seed, nrow(quality)),
            :fold => fill(fold, nrow(quality)))
        push!(joint_qc_tables, quality)
        # One context per inner selection group, so joint candidates can be scored by
        # predicting onto held-out stations instead of leave-one-out. Same builder,
        # same role map, inner indices — all scaling refits on the inner training set.
        if selection_groups !== nothing
            joint_selection_contexts = [(
                target_positions=group,
                context=build_joint_fold_context(
                    product, role_map,
                    train_idx[setdiff(1:length(train_idx), group)],
                    train_idx[group], lonlat, y_obs, y_sat,
                    joint_inputs.terrain, joint_inputs.era5.values,
                    joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
                    joint,
                ),
            ) for group in selection_groups]
        end
    end

    hurdle_context = build_hurdle_context(
        data.times, cfg, hurdle_rows; scheme, product, fold, repeat=repeat_index,
    )

    fold_predictions = merge(
        Dict{String,Matrix{Float64}}("raw" => y_sat_val), null_predictions,
    )
    scan_start = length(all_scan_rows) + 1
    # Winning hyperparameters per method, kept so `auto` can re-predict them across the
    # inner split and choose between the methods without seeing a held-out station.
    fold_selected = Dict{String,Any}()
    for (mode, method) in BENCHMARK_RUNS
        output_method = _output_method(mode, method)
        try
            selected = select_interpolation_parameter!(
                all_scan_rows, cfg, method, mode, scheme_symbol, product, fold,
                train_lonlat, y_obs_train, y_sat_train;
                dem_context, joint_context, hurdle_context,
                selection_groups, joint_selection_contexts, repeat_seed,
                repeat_index,
            )
            fold_predictions[output_method] = predict_selected(
                selected, method, mode, train_lonlat, val_lonlat,
                y_obs_train, y_sat_train, y_sat_val;
                dem_context, joint_context, hurdle_context,
            )
            predictions[output_method][val_idx, :] = fold_predictions[output_method]
            fold_selected[output_method] = selected
            eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
            coverage = _prediction_coverage(eligible, fold_predictions[output_method])
            push!(run_status_rows, _benchmark_status_row(
                dem_context, joint_context, joint_inputs, nested_joint;
                scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                mode, method, output_method,
                status=coverage >= cfg.min_tuning_coverage ? "success" : "partial",
                error=coverage >= cfg.min_tuning_coverage ? "" :
                    "prediction coverage below minimum",
                prediction_coverage=coverage,
            ))
        catch e
            fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
            push!(run_status_rows, _benchmark_status_row(
                dem_context, joint_context, joint_inputs, nested_joint;
                scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                mode, method, output_method, status="failed",
                error=sprint(showerror, e), prediction_coverage=0.0,
            ))
        end
    end

    # `mixed_gwr` and `residual_gwr` are the same model whenever this fold's role map has no
    # "global" role, so the second of the two scans is redundant work. It is kept rather than
    # skipped, and turned into a check: the two search the same grid through the same scorer over
    # the same designs, so if they ever disagree here, the designs or the grids have drifted
    # apart and the `duplicate_of` column in `run_status.csv` is lying.
    if joint_context !== nothing &&
            haskey(fold_predictions, "mixed_gwr") && haskey(fold_predictions, "residual_gwr") &&
            joint_models_coincide(joint_context, "mixed_gwr", "residual_gwr") &&
            !isequal(fold_predictions["mixed_gwr"], fold_predictions["residual_gwr"])
        @warn(
            "mixed_gwr and residual_gwr have identical designs on this fold but produced " *
            "different predictions",
            scheme, product, fold, repeat=repeat_index,
        )
    end

    for index in scan_start:length(all_scan_rows)
        all_scan_rows[index] = merge(
            all_scan_rows[index], (; repeat=repeat_index, seed=repeat_seed),
        )
    end

    auto_inner = _run_fold_auto!(
        cfg, fold, scheme, product, repeat_index, repeat_seed, val_idx, y_obs, y_sat_val,
        y_obs_train, y_sat_train, train_lonlat, selection_groups, joint_selection_contexts,
        joint_inputs, nested_joint, fold_selected, fold_predictions, predictions,
        auto_selection_rows, run_status_rows,
    )
    if cfg.satellite_wet_blend
        blend_band_train = Dict(axis => band[train_idx, :] for (axis, band) in blend_bands)
        blend_band_val = Dict(axis => band[val_idx, :] for (axis, band) in blend_bands)
        _run_fold_blend!(
            cfg, fold, scheme, product, repeat_index, repeat_seed, val_idx, y_obs, y_sat_val,
            y_obs_train, y_sat_train, train_lonlat, selection_groups, auto_inner,
            fold_selected, fold_predictions, predictions, blend_band_train, blend_band_val,
            blend_selection_rows, run_status_rows,
        )
    end
    fold_mask = _common_method_mask(Matrix{Float64}(y_obs[val_idx, :]), fold_predictions)
    if any(fold_mask)
        for method in benchmark_methods(cfg)
            append_stratified_metrics!(
                all_metric_rows, scheme, product, method, data.times,
                Matrix{Float64}(y_obs[val_idx, :]), fold_predictions[method], fold_mask,
                nearest_train_distance[val_idx], cfg.event_thresholds; fold,
                repeat=repeat_index, seed=repeat_seed,
            )
        end
    end
    return nothing
end

"""
The station ids and row indices a fold trains on and validates on.

Lifted out of `_run_benchmark_fold!` so the product loop can build a fold's fused anchor against
exactly the training stations the fold will then use, rather than deriving the split twice.
"""
function _fold_station_indices(folds, id_map, k::Int, fold::Int)
    val_ids = folds[fold]
    train_ids = reduce(vcat, (folds[index] for index in 1:k if index != fold))
    return train_ids, val_ids, [id_map[id] for id in train_ids], [id_map[id] for id in val_ids]
end

"""
Products derived by fusing the ones that have files, keyed by product name.

Each entry carries the variant and the source matrices, which are the real products' own `Y_sat`
in `products` order - so the fusion's design is the same one everywhere, and
`fused_anchor_selection.csv` can name its coefficients.

The derived product shares the source products' `times`, `ids` and `Y_obs`:
`load_global_common_product_data` has already put every product on one common station set and one
common time grid, so there is nothing to align. Its `Y_sat` starts as NaN and is filled per fold,
because the anchor depends on which stations were held out - which is exactly what keeps a held-out
gauge out of the anchor its own prediction is built on.
"""
function _build_fused_products!(
    cfg::InterpolationBenchmarkConfig, products::Vector{String}, product_data::Dict{String,Any},
)
    isempty(cfg.fused_anchor_variants) && return String[], Dict{String,Any}()
    source_names = copy(products)
    sources = [Matrix{Float64}(product_data[name].Y_sat) for name in source_names]
    template = product_data[first(source_names)]
    derived = String[]
    specifications = Dict{String,Any}()
    for variant in sort(cfg.fused_anchor_variants)
        name = _fused_product_name(variant)
        haskey(product_data, name) &&
            throw(ArgumentError("derived product $name collides with a configured product"))
        product_data[name] = (;
            template.times, template.ids, template.Y_obs,
            Y_sat=fill(NaN, size(template.Y_obs)),
        )
        specifications[name] = (; variant, source_names, sources)
        push!(derived, name)
    end
    return derived, specifications
end

function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    source_products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))
    # Derived products are appended after validation, and the joint inputs below are loaded for the
    # source products only: `load_joint_benchmark_inputs` reads a fixed specification keyed by
    # product, which a product with no file cannot appear in, and which the config validation has
    # already ruled out alongside a fused anchor.
    _, fused_specifications = _build_fused_products!(cfg, source_products, product_data)
    products = benchmark_products(cfg, source_products)
    terrain = _dem_enabled(cfg) ? load_aligned_terrain(something(cfg.terrain_path), ids) : nothing
    dem_store = _empty_dem_store()
    common_times = product_data[first(products)].times
    joint_inputs = _joint_enabled(cfg) ? load_joint_benchmark_inputs(
        something(cfg.joint_covariates), source_products, ids, common_times,
    ) : nothing
    nested_joint = _joint_enabled(cfg) && cfg.joint_selection !== nothing
    # Band matrices for the non-constant blend axes. They are read off the source products alone,
    # so they depend on neither the fold nor the scheme nor which product is being corrected, and
    # are built once here rather than per fold. The constant axis keeps its own code path in
    # `_run_fold_blend!` so that turning another axis on cannot perturb what it reports.
    blend_bands = Dict{Symbol,Matrix{Int}}()
    if cfg.satellite_wet_blend
        band_sources = [Matrix{Float64}(product_data[name].Y_sat) for name in source_products]
        for axis in cfg.blend_axes
            axis === :constant && continue
            blend_bands[axis] = blend_band_matrix(
                axis, first(band_sources), band_sources, cfg.mger.rain_threshold,
            )
        end
    end
    joint_store = _empty_joint_store()
    joint_scaling_tables = DataFrame[]
    joint_qc_tables = DataFrame[]
    joint_inputs === nothing || _write_joint_provenance(cfg, joint_inputs, nested_joint)

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]
    auto_selection_rows = NamedTuple[]
    blend_selection_rows = NamedTuple[]
    fused_anchor_rows = NamedTuple[]
    hurdle_rows = NamedTuple[]
    # Cells where the fold was too small for an inner selection split and fell back to
    # leave-one-out. Reported in `benchmark_scope.csv` so the fallback is never silent.
    selection_fallback_cells = String[]

    _dem_enabled(cfg) && _screen_dem_full_data!(dem_store, cfg, terrain, lonlat, products, product_data)

    seeds = benchmark_seeds(cfg)
    # Repeated cross-validation: one independent fold partition per repeat.
    #
    # `repeat_seed` and the partition index are deliberately separate. Under the default
    # `:hilbert` initialisation the partition comes from `rotation = repeat_index - 1` and no RNG
    # is involved at all; `repeat_seed` still seeds the genuinely stochastic sub-processes below
    # (paired bootstrap, DEM permutation tests, joint variable selection) and the `:random`
    # scheme, which is random by definition.
    for (repeat_index, repeat_seed) in enumerate(seeds)
        repeat_root = _repeat_dir(cfg.mger.outdir, repeat_index, length(seeds))
        for scheme_symbol in cfg.cv_schemes
            scheme = string(scheme_symbol)
            scheme_dir = joinpath(repeat_root, scheme)
            mkpath(scheme_dir)
            folds = benchmark_folds(
                scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
                center_init=cfg.fold_center_init, rotation=repeat_index - 1,
            )
            _write_split(joinpath(scheme_dir, "split_common.csv"), ids, folds, scheme_symbol)
            id_map = Dict(id => index for (index, id) in enumerate(ids))

            for product in products
                data = product_data[product]
                y_obs = data.Y_obs
                y_sat = data.Y_sat
                fusion = get(fused_specifications, product, nothing)
                predictions = Dict(
                    method => fill(NaN, size(y_obs)) for method in benchmark_methods(cfg)
                )
                # A derived product has no anchor until a fold says which stations may be used to
                # fit one, so `raw` is filled per fold below instead of in one assignment here.
                fusion === nothing && (predictions["raw"] .= y_sat)
                nearest_train_distance = fill(NaN, length(ids))

                for fold in 1:cfg.k
                    fold_sat = y_sat
                    if fusion !== nothing
                        _, _, fusion_train_idx, fusion_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        coefficients, used, fell_back = fusion_coefficients(
                            fusion.variant, y_obs, fusion.sources, fusion_train_idx,
                        )
                        # Fitted on the training stations, applied everywhere: the fold's own
                        # training rows need the same anchor its held-out rows get, or the
                        # residual target would mean two different things inside one fit.
                        fold_sat = apply_satellite_fusion(fusion.sources, coefficients)
                        predictions["raw"][fusion_val_idx, :] = fold_sat[fusion_val_idx, :]
                        push!(fused_anchor_rows, merge(
                            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                                variant=String(fusion.variant),
                                n_train_station=length(fusion_train_idx), n_train_cell=used,
                                fell_back, intercept=coefficients[1]),
                            NamedTuple{Tuple(Symbol("beta_", lowercase(name))
                                             for name in fusion.source_names)}(
                                Tuple(coefficients[2:end])),
                        ))
                    end
                    _run_benchmark_fold!(
                        cfg, fold, folds, id_map, ids, products, product,
                        data, lonlat, y_obs, fold_sat, terrain, joint_inputs, nested_joint,
                        scheme, scheme_symbol, repeat_index, repeat_seed,
                        predictions, nearest_train_distance, blend_bands,
                        dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
                        all_metric_rows, all_scan_rows, run_status_rows, auto_selection_rows,
                        blend_selection_rows, hurdle_rows, selection_fallback_cells,
                    )
                end

                common_mask = _common_method_mask(y_obs, predictions)
                if !any(common_mask)
                    failures = ["fold=$(row.fold) method=$(row.method): $(row.error)" for
                        row in run_status_rows if row.scheme == scheme && row.product == product &&
                        row.status != "success"]
                    error("[$scheme/$product] no common valid OOF samples across all methods; " *
                        join(failures, " | "))
                end
                product_dir = joinpath(scheme_dir, lowercase(product))
                mkpath(product_dir)
                for method in benchmark_methods(cfg)
                    # The per-station OOF tables are large; only the first repeat writes them.
                    repeat_index == 1 && write_wide(
                        joinpath(product_dir, "oof_$(method).csv"), data.times, ids, predictions[method],
                    )
                    append_stratified_metrics!(
                        all_metric_rows, scheme, product, method, data.times, y_obs,
                        predictions[method], common_mask, nearest_train_distance, cfg.event_thresholds;
                        repeat=repeat_index, seed=repeat_seed,
                    )
                end
                if repeat_index == 1
                    mask_df = DataFrame(time=Dates.format.(data.times, dateformat"yyyy-mm-ddTHH:MM:SS"))
                    for (station_index, station_id) in enumerate(ids)
                        mask_df[!, Symbol(station_id)] = common_mask[station_index, :]
                    end
                    CSV.write(joinpath(product_dir, "common_evaluation_mask.csv"), mask_df)
                end
                if scheme_symbol == :balanced_spatial && cfg.bootstrap_reps > 0
                    append!(all_bootstrap_rows, paired_bootstrap_rows(
                        cfg, scheme, product, data.times, y_obs, predictions, common_mask;
                        repeat=repeat_index, seed=repeat_seed,
                    ))
                end
            end
        end
    end # repeat loop

    return _write_benchmark_outputs(
        cfg, products, seeds, nested_joint, joint_inputs,
        dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
        all_metric_rows, all_scan_rows, all_bootstrap_rows, run_status_rows,
        auto_selection_rows, blend_selection_rows, fused_anchor_rows, hurdle_rows,
        selection_fallback_cells,
    )
end
