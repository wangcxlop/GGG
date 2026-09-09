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
    return (; uses_dem, uses_joint, selected_roles, effective_role_map, effective_roles)
end

"""Non-NaN share of the held-out cells a method could have predicted."""
_prediction_coverage(eligible, prediction) =
    count(eligible) == 0 ? 0.0 :
        count(eligible .& .!isnan.(prediction)) / count(eligible)

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
    (; uses_dem, uses_joint, selected_roles, effective_roles) =
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
    auto_selection_rows, hurdle_rows, selection_fallback_cells,
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
    CSV.write(joinpath(cfg.mger.outdir, "benchmark_scope.csv"), scope)
    return (; metrics, scans, bootstrap, status, claim, repeat_summary, fold_summary,
        rank_stability, auto_selection=DataFrame(auto_selection_rows))
end

"""
`auto` for one fold: choose a single GWR-family method using only the inner selection split, and
record the choice.

Split out of `_run_benchmark_fold!`, which was 262 lines. Mutates `fold_predictions` (writing the
`auto` entry, NaN when no contender could be scored), `auto_selection_rows` and `run_status_rows`.
Parameters keep the names the enclosing locals had, so the body is unchanged from when it was
inline.
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
    # "skipped" rather than "failed" where `auto` was never applicable, so a genuine
    # breakage stays visible instead of being lost among expected non-runs.
    auto_status, auto_error = if auto_choice !== nothing
        ("success", "")
    elseif _dem_enabled(cfg)
        ("skipped", "auto is not run on the legacy DEM path")
    elseif selection_groups === nothing
        ("skipped", "fold has no inner selection split to choose on")
    else
        ("failed", "no auto contender scored: $(join(auto_failures, " | "))")
    end
    # Measured the same way as every other method's row - non-NaN share of the
    # held-out cells it could have predicted - so the column means one thing across the
    # table. `auto_selection.csv` carries the inner-split mask coverage separately.
    auto_eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
    auto_coverage = _prediction_coverage(auto_eligible, fold_predictions[AUTO_METHOD])
    push!(run_status_rows, _benchmark_status_row(
        nothing, nothing, joint_inputs, nested_joint;
        scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
        mode="", method=AUTO_METHOD, output_method=AUTO_METHOD,
        status=auto_status, error=auto_error, prediction_coverage=auto_coverage,
    ))
    return nothing
end

"""
Everything one cross-validation fold does: split the stations, build the fold's DEM/joint/hurdle
contexts, tune and predict each `BENCHMARK_RUNS` method, run `auto`, and append the fold's metric,
scan and status rows.

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
    predictions, nearest_train_distance,
    dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
    all_metric_rows, all_scan_rows, run_status_rows, auto_selection_rows, hurdle_rows,
    selection_fallback_cells,
)

    val_ids = folds[fold]
    train_ids = reduce(vcat, (folds[index] for index in 1:cfg.k if index != fold))
    train_idx = [id_map[id] for id in train_ids]
    val_idx = [id_map[id] for id in val_ids]
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

    for index in scan_start:length(all_scan_rows)
        all_scan_rows[index] = merge(
            all_scan_rows[index], (; repeat=repeat_index, seed=repeat_seed),
        )
    end

    _run_fold_auto!(
        cfg, fold, scheme, product, repeat_index, repeat_seed, val_idx, y_obs, y_sat_val,
        y_obs_train, y_sat_train, train_lonlat, selection_groups, joint_selection_contexts,
        joint_inputs, nested_joint, fold_selected, fold_predictions, predictions,
        auto_selection_rows, run_status_rows,
    )
    fold_mask = _common_method_mask(Matrix{Float64}(y_obs[val_idx, :]), fold_predictions)
    if any(fold_mask)
        for method in BENCHMARK_METHODS
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

function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))
    terrain = _dem_enabled(cfg) ? load_aligned_terrain(something(cfg.terrain_path), ids) : nothing
    dem_store = _empty_dem_store()
    common_times = product_data[first(products)].times
    joint_inputs = _joint_enabled(cfg) ? load_joint_benchmark_inputs(
        something(cfg.joint_covariates), products, ids, common_times,
    ) : nothing
    nested_joint = _joint_enabled(cfg) && cfg.joint_selection !== nothing
    joint_store = _empty_joint_store()
    joint_scaling_tables = DataFrame[]
    joint_qc_tables = DataFrame[]
    joint_inputs === nothing || _write_joint_provenance(cfg, joint_inputs, nested_joint)

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]
    auto_selection_rows = NamedTuple[]
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
                predictions = Dict(method => fill(NaN, size(y_obs)) for method in BENCHMARK_METHODS)
                predictions["raw"] .= y_sat
                nearest_train_distance = fill(NaN, length(ids))

                for fold in 1:cfg.k
                    _run_benchmark_fold!(
                        cfg, fold, folds, id_map, ids, products, product,
                        data, lonlat, y_obs, y_sat, terrain, joint_inputs, nested_joint,
                        scheme, scheme_symbol, repeat_index, repeat_seed,
                        predictions, nearest_train_distance,
                        dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
                        all_metric_rows, all_scan_rows, run_status_rows, auto_selection_rows,
                        hurdle_rows, selection_fallback_cells,
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
                for method in BENCHMARK_METHODS
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
        auto_selection_rows, hurdle_rows, selection_fallback_cells,
    )
end
