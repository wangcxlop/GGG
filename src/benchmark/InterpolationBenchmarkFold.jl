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
`STACK_METHOD` for one fold: fit the banded non-negative stack on the inner selection split, then
apply it to the fold's held-out stations.

The inner design is each member re-predicted across the inner split over the full training record
(`time_indices=nothing`), not `auto`'s tuning hours - the lagged features need every hour's
neighbours, and the tuning sample's wet oversampling would skew the dry band's weights. `raw` is the
anchor itself. Distances on the inner split are to the nearest *inner*-training station, so the
distance bands mean the same thing there as on the fold.

Skips, leaving the method NaN with a status row saying why, when there is no inner split or a
member did not fit - the same reason `_run_fold_blend!` gives for not defaulting to the source.

Mutates `fold_predictions`, `predictions`, `stack_selection_rows` and `run_status_rows`.
"""
function _run_fold_stack!(
    cfg::InterpolationBenchmarkConfig, fold::Int, scheme, product, repeat_index::Int,
    repeat_seed::Int, times, val_idx, y_obs, y_sat_val, y_obs_train, y_sat_train, train_lonlat,
    selection_groups, joint_selection_contexts, fold_selected, fold_predictions, predictions,
    agreement_train, agreement_val, val_distance, stack_selection_rows, run_status_rows,
)
    eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
    status_row(status, message, coverage) = _benchmark_status_row(
        nothing, nothing, nothing, false;
        scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
        mode="", method=STACK_METHOD, output_method=STACK_METHOD,
        status, error=message, prediction_coverage=coverage,
    )
    function skip(status::String, message::String)
        fold_predictions[STACK_METHOD] = fill(NaN, length(val_idx), size(y_obs, 2))
        push!(run_status_rows, status_row(status, message, 0.0))
        return nothing
    end

    (_dem_enabled(cfg) || selection_groups === nothing) &&
        return skip("skipped", "fold has no inner selection split to fit the stack on")
    for (mode, method) in STACK_MEMBER_RUNS
        output_method = _output_method(mode, method)
        haskey(fold_selected, output_method) && haskey(fold_predictions, output_method) ||
            return skip("skipped", "$output_method did not fit on this fold")
    end

    inner = Dict{String,Matrix{Float64}}()
    for (mode, method) in STACK_MEMBER_RUNS
        output_method = _output_method(mode, method)
        inner[output_method] = try
            inner_selection_prediction(
                fold_selected[output_method], method, mode, train_lonlat,
                y_obs_train, y_sat_train, selection_groups;
                joint_contexts=joint_selection_contexts,
            )
        catch e
            return skip("failed",
                "$output_method inner prediction failed: $(sprint(showerror, e))")
        end
    end

    distance_train = haversine_distance_matrix(train_lonlat, train_lonlat)
    inner_distance = fill(NaN, size(train_lonlat, 1))
    for group in selection_groups
        inner_train = setdiff(axes(train_lonlat, 1), group)
        isempty(inner_train) && continue
        inner_distance[group] = vec(minimum(distance_train[inner_train, group]; dims=1))
    end

    n_bands = stack_band_count()
    inner_features = stack_feature_matrices(
        inner["mgwr"], inner["adw"], y_sat_train, inner["tps"], times,
    )
    choice = select_stack_weights(
        y_obs_train, inner_features, stack_band_matrix(agreement_train, inner_distance), n_bands,
    )
    choice === nothing &&
        return skip("failed", "no inner-split cell was scorable for every stack member")

    stacked = stacked_prediction(
        stack_feature_matrices(
            fold_predictions["mgwr"], fold_predictions["adw"], y_sat_val,
            fold_predictions["tps"], times,
        ),
        stack_band_matrix(agreement_val, val_distance), choice.weights,
    )
    fold_predictions[STACK_METHOD] = stacked
    predictions[STACK_METHOD][val_idx, :] = stacked
    for band in 1:n_bands
        push!(stack_selection_rows, merge(
            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                method=STACK_METHOD, band, n_cell=choice.used[band],
                fell_back=choice.fell_back[band]),
            NamedTuple{Tuple(Symbol.("w_" .* STACK_FEATURES))}(Tuple(choice.weights[:, band])),
            (; inner_RMSE=choice.inner_RMSE, source_RMSE=choice.source_RMSE),
        ))
    end
    coverage = _prediction_coverage(eligible, stacked)
    push!(run_status_rows, status_row(
        coverage >= cfg.min_tuning_coverage ? "success" : "partial",
        coverage >= cfg.min_tuning_coverage ? "" : "prediction coverage below minimum",
        coverage))
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
    stack_agreement=nothing, stack_selection_rows=NamedTuple[],
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
    if cfg.stacked_blend
        _run_fold_stack!(
            cfg, fold, scheme, product, repeat_index, repeat_seed, data.times, val_idx, y_obs,
            y_sat_val, y_obs_train, y_sat_train, train_lonlat, selection_groups,
            joint_selection_contexts, fold_selected, fold_predictions, predictions,
            stack_agreement[train_idx, :], stack_agreement[val_idx, :],
            nearest_train_distance[val_idx], stack_selection_rows, run_status_rows,
        )
    end
    fold_mask =_common_method_mask(Matrix{Float64}(y_obs[val_idx, :]), fold_predictions)
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
