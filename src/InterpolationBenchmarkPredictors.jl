function _gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64}, target_lonlat::Matrix{Float64};
    kernel::Int, adaptive::Bool, bw::Float64, exclude_self::Bool=false,
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center=center)
    X_target = build_X_intercept_centered(target_lonlat; center=center)
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("GWR exclude_self requires matching rows"))
        for i in 1:size(distances, 1)
            distances[i, i] = Inf
        end
    end
    weights = gw_weight(distances, bw; kernel=kernel, adaptive=adaptive)
    if exclude_self
        # The `Inf` distance above is what keeps self out of the adaptive neighbour ranking, but
        # it does not by itself produce a zero weight: at the global candidate (`bw = Inf`)
        # `kernel(Inf, Inf)` is NaN for four of the five kernels and 1.0 for boxcar, which would
        # leak the held-out station straight back into its own fit. Clear it explicitly.
        for i in axes(weights, 1)
            weights[i, i] = 0.0
        end
    end
    return st_gwr_predict_nanaware(X_train, values, weights; Xpred=X_target, min_obs=3)
end

"""
Per-fold context for `hurdle_gwr`.

`fit_hurdle_gwr`/`predict_hurdle_gwr` need the fold's `times` vector, which the rest of the
benchmark's method dispatch never has to carry; this mirrors the existing `dem_context` /
`joint_context` keyword pattern rather than widening every signature. `diagnostics` is a shared
accumulator (same pattern as `scan_rows`) for the degeneracy report described at
`_hurdle_predict`.
"""
function build_hurdle_context(
    times::Vector{DateTime}, cfg::InterpolationBenchmarkConfig,
    diagnostics::Vector{NamedTuple}; scheme::String, product::String, fold::Int,
    repeat::Int=1,
)
    wet_threshold = cfg.mger.rain_threshold
    # `HurdleGWRConfig` requires `first(intensity_breaks) == wet_threshold`, sorted and unique.
    # Deriving the breaks from the benchmark's own event thresholds keeps the occurrence table
    # binned the same way the categorical scores are.
    breaks = sort(unique(vcat(wet_threshold, filter(>(wet_threshold), cfg.event_thresholds))))
    return (; times, wet_threshold, intensity_breaks=breaks, diagnostics,
        scheme, product, fold, repeat)
end

"""Adapt the benchmark's `(adaptive, bw)` convention to HurdleGWR's `SpatialBandwidth`."""
_hurdle_bandwidth(adaptive::Bool, bw::Float64) =
    adaptive ? AdaptiveBandwidth(Int(round(bw))) : FixedBandwidth(bw)

_hurdle_config(context, kernel::Int, adaptive::Bool, bw::Float64) = HurdleGWRConfig(
    wet_threshold=context.wet_threshold,
    intensity_breaks=context.intensity_breaks,
    kernel=kernel,
    bandwidth=_hurdle_bandwidth(adaptive, bw),
)

"""
Fit the hurdle model on the training fold and predict `P(wet) * E[amount | wet]`.

Unlike every other method here this is a time-pooled calibration, not a per-hour spatial fit:
one local regression per target station, with time entering through annual/diurnal harmonics and
a calendar-month occurrence bin. The satellite enters the design as a *fitted* `log1p_satellite`
coefficient rather than as an offset with its coefficient forced to 1.

`predict_hurdle_gwr` takes no gauge observations at the target at all, so held-out leakage is
impossible by construction. `exclude_self` drives its `exclude_train_indices`, which zeroes a
station's own spatial weight - the same leave-one-out convention `_gwr_predict` uses.

Returns the full prediction object: the caller needs `used_global_amount` and
`occurrence_fallback` as well as `corrected`, because when a local fit misses
`min_amount_samples` / `min_effective_stations` the model *silently* falls back to global
coefficients, and an undiagnosed fallback would report a global model as a local one.
"""
function _hurdle_predict(
    train_lonlat::Matrix{Float64}, y_obs_train::Matrix{Float64},
    y_sat_train::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_sat_target::Matrix{Float64}, times::Vector{DateTime};
    config::HurdleGWRConfig, exclude_self::Bool=false,
)
    model = fit_hurdle_gwr(y_obs_train, y_sat_train, train_lonlat, times; config)
    exclusions = exclude_self ? collect(1:size(train_lonlat, 1)) : nothing
    return predict_hurdle_gwr(
        model, y_sat_target, target_lonlat, times; exclude_train_indices=exclusions,
    )
end

"""Summarise how much of a hurdle prediction actually came from a *local* fit."""
function _hurdle_diagnostic_row(prediction, context; kernel::Int, adaptive::Bool, bw::Float64)
    fallback = prediction.occurrence_fallback
    total = length(fallback)
    share(code) = total == 0 ? NaN : count(==(code), fallback) / total
    return (;
        scheme=context.scheme, product=context.product, fold=context.fold,
        repeat=context.repeat, kernel, adaptive, bw,
        n_target=length(prediction.used_global_amount),
        used_global_amount_fraction=mean(prediction.used_global_amount),
        mean_effective_stations=mean(prediction.effective_stations),
        occurrence_monthly=share(0x00), occurrence_all_month=share(0x01),
        occurrence_global=share(0x02), occurrence_skipped=share(0x03),
    )
end

"""
Hyperparameter-free reference predictors for one fold, estimated from training stations only.

- `zero`: constant 0. Hourly precipitation is ~93% dry, so this is a surprisingly strong
  RMSE competitor and the floor any method must clear.
- `train_clim`: the training stations' overall mean. A per-station climatology is not
  estimable for a held-out station, so the pooled mean is the station-free analogue.
- `hour_field_mean`: each hour's spatial mean over training stations - "is it raining
  anywhere in the domain right now", carrying no spatial structure at all. An interpolator
  that does not beat this is not doing spatial work.

Reported alongside the real methods so `metrics_*.csv` can be read as skill rather than as a
bare RMSE. See `NULL_METHODS`.
"""
function _null_fold_predictions(y_obs_train::Matrix{Float64}, n_val::Int)
    n_time = size(y_obs_train, 2)
    finite = filter(isfinite, vec(y_obs_train))
    climatology = isempty(finite) ? NaN : mean(finite)
    field_mean = Matrix{Float64}(undef, n_val, n_time)
    @inbounds for time in 1:n_time
        total = 0.0
        count = 0
        for station in axes(y_obs_train, 1)
            value = y_obs_train[station, time]
            if !isnan(value)
                total += value
                count += 1
            end
        end
        field_mean[:, time] .= count > 0 ? total / count : NaN
    end
    return Dict{String,Matrix{Float64}}(
        "zero" => zeros(Float64, n_val, n_time),
        "train_clim" => fill(climatology, n_val, n_time),
        "hour_field_mean" => field_mean,
    )
end

"""Build the provisional mixed-GWR design; replace this when variable roles are finalized."""
function build_mixed_gwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    local_train = build_X_intercept_centered(train_lonlat; center)
    local_target = build_X_intercept_centered(target_lonlat; center)
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    return (; local_train, local_target, global_train, global_target)
end

function _mixed_gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bw::Float64, kernel::Int, adaptive::Bool=true,
    exclude_self::Bool=false,
)
    designs = build_mixed_gwr_designs(train_lonlat, target_lonlat)
    prediction, _ = mixed_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, bw,
        _kernel_function(kernel); adaptive, exclude_self,
    )
    return prediction
end

"""Build the provisional MGWR groups; replace these when covariates are finalized."""
function build_mgwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center)
    X_target = build_X_intercept_centered(target_lonlat; center)
    local_train = [X_train[:, index:index] for index in axes(X_train, 2)]
    local_target = [X_target[:, index:index] for index in axes(X_target, 2)]
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    group_names = ["intercept", "longitude", "latitude"]
    return (; local_train, local_target, global_train, global_target, group_names)
end

function _mgwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bandwidths::Vector{Float64}, kernel::Int,
    adaptive::Bool=true, exclude_self::Bool=false,
)
    designs = build_mgwr_designs(train_lonlat, target_lonlat)
    prediction, _ = multiscale_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, bandwidths,
        _kernel_function(kernel); adaptive, exclude_self,
    )
    return prediction
end

function predict_selected(
    selected, method::String, mode::String,
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64}, y_sat_target::Matrix{Float64};
    dem_context=nothing, joint_context=nothing, hurdle_context=nothing,
)
    (mode, method) in BENCHMARK_RUNS ||
        throw(ArgumentError("unsupported benchmark method/mode pair: $method/$mode"))
    values = mode == "direct" ? y_obs_train : y_obs_train .- y_sat_train
    is_dem_model = dem_context !== nothing && hasproperty(selected, :dem_model) && selected.dem_model
    is_joint_model = joint_context !== nothing &&
        hasproperty(selected, :joint_model) && selected.joint_model
    interpolated = if is_joint_model && mode == "residual"
        prediction, converged = dynamic_covariate_predict(
            joint_context, values, selected.joint_method,
            Float64.(selected.bandwidths), _kernel_function(selected.kernel);
            adaptive=selected.adaptive,
        )
        any(.!converged) && @warn(
            "joint dynamic model did not converge for some hours",
            product=joint_context.product, method=selected.joint_method,
            failed_hours=count(.!converged),
        )
        # Folded in here rather than at the shared `y_sat_target .+ interpolated` tail below,
        # which serves every other method and must keep its historical behaviour.
        shrink = hasproperty(selected, :shrink) ? selected.shrink : 1.0
        isfinite(shrink) ? prediction .* shrink : prediction
    elseif is_dem_model && method == "gwr" && mode == "residual"
        designs = dem_context.all_local
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, zeros(Float64, size(train_lonlat, 1), 0), values,
            train_lonlat, designs.mixed_local_target,
            zeros(Float64, size(target_lonlat, 1), 0), target_lonlat,
            selected.bw, _kernel_function(selected.kernel); adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mixed_gwr"
        designs = dem_context.mixed
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, designs.global_train, values, train_lonlat,
            designs.mixed_local_target, designs.global_target, target_lonlat,
            selected.bw, _kernel_function(selected.kernel); adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mgwr"
        designs = dem_context.mixed
        prediction, _ = multiscale_gwr_predict(
            designs.multiscale_train, designs.global_train, values, train_lonlat,
            designs.multiscale_target, designs.global_target, target_lonlat,
            selected.bandwidths, _kernel_function(selected.kernel);
            adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif method in ("idw", "adw")
        selected_neighbors = ismissing(selected.neighbors) ? nothing :
            (selected.neighbors == 0 ? nothing : Int(selected.neighbors))
        predictor = method == "idw" ? idw_predict : adw_predict
        predictor(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "tps"
        tps_predict(train_lonlat, values, target_lonlat; smooth=selected.smooth)
    elseif method == "gwr"
        _gwr_predict(train_lonlat, values, target_lonlat;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw)
    elseif method == "hurdle_gwr"
        hurdle_context === nothing &&
            throw(ArgumentError("hurdle_gwr requires a hurdle_context carrying the fold times"))
        prediction = _hurdle_predict(
            train_lonlat, y_obs_train, y_sat_train, target_lonlat, y_sat_target,
            hurdle_context.times;
            config=_hurdle_config(
                hurdle_context, selected.kernel, selected.adaptive, selected.bw,
            ),
        )
        push!(hurdle_context.diagnostics, _hurdle_diagnostic_row(
            prediction, hurdle_context;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw,
        ))
        prediction.corrected
    elseif method == "mixed_gwr"
        _mixed_gwr_predict(
            train_lonlat, values, target_lonlat;
            bw=selected.bw, kernel=selected.kernel, adaptive=selected.adaptive,
        )
    elseif method == "mgwr"
        _mgwr_predict(
            train_lonlat, values, target_lonlat; bandwidths=Float64.(selected.bandwidths),
            kernel=selected.kernel, adaptive=selected.adaptive,
        )
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat_target .+ interpolated, 0.0)
end

"""
Refit a method's already-selected hyperparameters across the inner selection split and stitch the
held-out inner predictions into one `n_train × n_time` matrix.

This is what makes `auto` a fair comparison rather than a scan-row lookup. Each method's
`selected` scan row already carries an inner-split RMSE, but those numbers are scored on slightly
different cell sets: a method is only required to reach `min_tuning_coverage`, and in the full run
the selected candidates' coverage ranges 0.98-1.00, with `mgwr` — the method most likely to win —
sitting lowest. Comparing RMSEs computed over different denominators would quietly reward whichever
method dropped the hardest cells, on margins of about a percent. Re-predicting here lets the caller
score every contender on the intersection of their masks.

Deliberately no `dem_context`: the legacy DEM path is mutually exclusive with the joint path and is
not part of the reported comparison, so `auto` does not run there. `joint_contexts` is
`joint_selection_contexts` — one `JointFoldContext` per inner group, already built for the scan, and
matched here by `target_positions` rather than by position so a reordering cannot silently misalign
a context with the wrong stations.

`time_indices` restricts the result to the tuning hours, so `auto` compares methods over exactly
the hours the scan scored them on. The two model families need it applied at different points: a
`JointFoldContext` indexes hours internally and `predict_selected` hands it the full residual
matrix, so the joint methods predict every hour and are subset afterwards, while everything else
can be given pre-sliced inputs and never touches the other hours at all. Predicting the full
record for the joint methods and slicing costs roughly 5-10% of their scan, which is why this
does not thread a `time_indices` keyword through `predict_selected` and risk the reported
prediction path for it.
"""
function inner_selection_prediction(
    selected, method::String, mode::String, train_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64},
    selection_groups::Vector{Vector{Int}}; joint_contexts=nothing,
    time_indices::Union{Nothing,Vector{Int}}=nothing,
)
    is_joint = joint_contexts !== nothing &&
        hasproperty(selected, :joint_model) && selected.joint_model
    columns = time_indices === nothing ? collect(axes(y_obs_train, 2)) : time_indices
    obs = is_joint ? y_obs_train : y_obs_train[:, columns]
    sat = is_joint ? y_sat_train : y_sat_train[:, columns]
    n_train = size(obs, 1)
    out_of_fold = _selection_oof(selection_groups, n_train, size(obs, 2), function (inner_train, group)
        joint_context = if joint_contexts === nothing
            nothing
        else
            entry = findfirst(candidate -> candidate.target_positions == group, joint_contexts)
            entry === nothing && throw(ArgumentError(
                "no joint selection context matches inner group $(group)",
            ))
            joint_contexts[entry].context
        end
        return predict_selected(
            selected, method, mode,
            train_lonlat[inner_train, :], train_lonlat[group, :],
            obs[inner_train, :], sat[inner_train, :], sat[group, :];
            joint_context,
        )
    end)
    return is_joint ? out_of_fold[:, columns] : out_of_fold
end
