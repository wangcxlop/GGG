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
    if joint_inputs !== nothing
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
    end

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]
    hurdle_rows = NamedTuple[]
    # Cells where the fold was too small for an inner selection split and fell back to
    # leave-one-out. Reported in `benchmark_scope.csv` so the fallback is never silent.
    selection_fallback_cells = String[]

    if _dem_enabled(cfg)
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
    end

    seeds = benchmark_seeds(cfg)
    # Repeated cross-validation: one independent fold partition per seed. The body is left at its
    # original indentation to keep the diff reviewable.
    for (repeat_index, repeat_seed) in enumerate(seeds)
    repeat_root = _repeat_dir(cfg.mger.outdir, repeat_index, length(seeds))
    for scheme_symbol in cfg.cv_schemes
        scheme = string(scheme_symbol)
        scheme_dir = joinpath(repeat_root, scheme)
        mkpath(scheme_dir)
        folds = benchmark_folds(
            scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
            center_init=cfg.fold_center_init,
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
                    selection_folds(cfg, scheme_symbol, train_ids, train_lonlat, fold, repeat_seed) :
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
                    scaling = copy(joint_context.scaling)
                    insertcols!(scaling, 1,
                        :scheme => fill(scheme, nrow(scaling)),
                        :fold => fill(fold, nrow(scaling)))
                    push!(joint_scaling_tables, scaling)
                    quality = copy(joint_context.quality_control)
                    insertcols!(quality, 1,
                        :scheme => fill(scheme, nrow(quality)),
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
                    data.times, cfg, hurdle_rows; scheme, product, fold,
                )

                fold_predictions = merge(
                    Dict{String,Matrix{Float64}}("raw" => y_sat_val), null_predictions,
                )
                scan_start = length(all_scan_rows) + 1
                for (mode, method) in BENCHMARK_RUNS
                    output_method = method in ("mixed_gwr", "mgwr") ? method :
                        (mode == "direct" ? method : "residual_$(method)")
                    try
                        selected = select_interpolation_parameter!(
                            all_scan_rows, cfg, method, mode, scheme_symbol, product, fold,
                            train_lonlat, y_obs_train, y_sat_train;
                            dem_context, joint_context, hurdle_context,
                            selection_groups, joint_selection_contexts, repeat_seed,
                        )
                        fold_predictions[output_method] = predict_selected(
                            selected, method, mode, train_lonlat, val_lonlat,
                            y_obs_train, y_sat_train, y_sat_val;
                            dem_context, joint_context, hurdle_context,
                        )
                        predictions[output_method][val_idx, :] = fold_predictions[output_method]
                        (; uses_dem, uses_joint, selected_roles, effective_role_map, effective_roles) =
                            _run_status_role_fields(dem_context, joint_context, mode, method, output_method)
                        eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
                        coverage = count(eligible) == 0 ? 0.0 :
                            count(eligible .& .!isnan.(fold_predictions[output_method])) / count(eligible)
                        push!(run_status_rows, (;
                            scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            method=output_method,
                            status=coverage >= cfg.min_tuning_coverage ? "success" : "partial",
                            error=coverage >= cfg.min_tuning_coverage ? "" :
                                "prediction coverage below minimum",
                            prediction_coverage=coverage,
                            dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
                            dem_selection_status=uses_dem ? dem_context.selection_status :
                                "not_applicable",
                            dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
                            dem_roles=uses_dem ? dem_context.dem_roles : "",
                            covariate_selection_mode=uses_joint ?
                                (nested_joint ? "nested_per_fold" : "fixed_full_data") :
                                "not_applicable",
                            covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
                            covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
                            covariate_selected_roles=selected_roles,
                            covariate_effective_roles=effective_roles,
                            covariate_spec_sha256=uses_joint ?
                                something(joint_inputs.spec_sha256, "") : "",
                        ))
                    catch e
                        fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
                        (; uses_dem, uses_joint, selected_roles, effective_role_map, effective_roles) =
                            _run_status_role_fields(dem_context, joint_context, mode, method, output_method)
                        push!(run_status_rows, (;
                            scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            method=output_method, status="failed",
                            error=sprint(showerror, e), prediction_coverage=0.0,
                            dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
                            dem_selection_status=uses_dem ? dem_context.selection_status :
                                "not_applicable",
                            dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
                            dem_roles=uses_dem ? dem_context.dem_roles : "",
                            covariate_selection_mode=uses_joint ?
                                (nested_joint ? "nested_per_fold" : "fixed_full_data") :
                                "not_applicable",
                            covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
                            covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
                            covariate_selected_roles=selected_roles,
                            covariate_effective_roles=effective_roles,
                            covariate_spec_sha256=uses_joint ?
                                something(joint_inputs.spec_sha256, "") : "",
                        ))
                    end
                end

                for index in scan_start:length(all_scan_rows)
                    all_scan_rows[index] = merge(
                        all_scan_rows[index], (; repeat=repeat_index, seed=repeat_seed),
                    )
                end

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
    rank_stability = method_rank_stability(metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_stratified.csv"), metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_folds.csv"), filter(:fold => (x -> !ismissing(x)), metrics))
    CSV.write(joinpath(cfg.mger.outdir, "metrics_pooled.csv"), filter(:fold => ismissing, metrics))
    nrow(repeat_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_repeat_summary.csv"), repeat_summary)
    nrow(rank_stability) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "method_rank_stability.csv"), rank_stability)
    CSV.write(joinpath(cfg.mger.outdir, "parameter_scan.csv"), scans)
    # How much of each hurdle prediction actually came from a local fit. Without this a
    # globally-degenerate model is indistinguishable from a local one in the metrics.
    isempty(hurdle_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "hurdle_diagnostics.csv"), DataFrame(hurdle_rows))
    nrow(bootstrap) > 0 && CSV.write(joinpath(cfg.mger.outdir, "paired_comparisons.csv"), bootstrap)
    CSV.write(joinpath(cfg.mger.outdir, "run_status.csv"), status)
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
    end
    scope = DataFrame(
        key=[
            "repeated_cv_partitions", "repeated_cv_seeds", "fold_center_init",
            "validation_target", "training_signal", "temporal_holdout", "primary_cv",
            "secondary_cv", "common_evaluation_mask", "tuning_time_limit",
            "tuning_time_weighting", "tuning_geometry",
            "dem_variable_selection", "dem_role_assignment", "dem_leakage_control",
            "dem_empty_selection", "supported_claim", "unsupported_claim",
        ],
        value=[
            string(length(seeds)), join(seeds, ","), string(cfg.fold_center_init),
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
    return (; metrics, scans, bootstrap, status, claim, repeat_summary, rank_stability)
end
