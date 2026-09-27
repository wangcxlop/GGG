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
        _predict_hurdle_selected(selected, train_lonlat, target_lonlat,
            y_obs_train, y_sat_train, y_sat_target, hurdle_context)
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

"""
Blend a residual-family prediction toward a gauge-only one on the cells the satellite calls wet.

`lambda = 0` returns `prediction` untouched, `lambda = 1` replaces it with `fallback` on those
cells and leaves the satellite-dry cells alone either way. The test is on the satellite, which is
knowable at prediction time - this is a model, not an oracle.

A NaN in `prediction` or the satellite propagates. A NaN in `fallback` does not: the cell simply
keeps its unblended prediction, so a fallback that failed on some cells cannot shrink the blended
method's coverage below its source's. That matters because coverage feeds the run's status rows.
"""
function satellite_wet_blend_prediction(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64}, fallback::Matrix{Float64},
    lambda::Float64, threshold::Float64,
)
    size(prediction) == size(y_sat) == size(fallback) || throw(DimensionMismatch(
        "satellite, prediction and fallback must have the same shape",
    ))
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        satellite = y_sat[index]
        predicted = prediction[index]
        other = fallback[index]
        if isnan(satellite) || isnan(predicted)
            out[index] = NaN
        elseif satellite >= threshold && !isnan(other)
            out[index] = (1 - lambda) * predicted + lambda * other
        else
            out[index] = predicted
        end
    end
    return out
end

"""
Blend toward the fallback with one weight per band, the band read off a precomputed matrix.

The vector-weight twin of `satellite_wet_blend_prediction`, and it keeps that function's rules
exactly: band 0 is never touched, a NaN prediction propagates, and a NaN fallback leaves the cell
at its unblended value rather than poisoning it. That last one is what makes the per-band
selection in `select_blend_lambdas` exact - the set of finite cells does not move as the weights
move, so pooled SSE is additive over the bands and each band minimises independently.

`band` is `0` for a cell no weight reaches, and otherwise an index into `lambdas`.
"""
function banded_blend_prediction(
    band::Matrix{Int}, prediction::Matrix{Float64}, fallback::Matrix{Float64},
    lambdas::Vector{Float64},
)
    size(prediction) == size(band) == size(fallback) || throw(DimensionMismatch(
        "band, prediction and fallback must have the same shape",
    ))
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        predicted = prediction[index]
        other = fallback[index]
        position = band[index]
        if isnan(predicted)
            out[index] = NaN
        elseif position == 0 || isnan(other)
            out[index] = predicted
        else
            1 <= position <= length(lambdas) || throw(BoundsError(lambdas, position))
            lambda = lambdas[position]
            out[index] = (1 - lambda) * predicted + lambda * other
        end
    end
    return out
end

"""
Band index per cell for `axis`, from the source products the fusion also draws on.

`:constant` reproduces the scalar intervention set exactly - band 1 wherever the *target's* own
anchor calls the cell wet - so a one-band run of the banded path is the scalar path.

`:agreement_envelope` bands on the source products instead: `(agreement - 1) * 3 + envelope`,
where `agreement` counts how many products call the cell wet and `envelope` places `max` over them
in `BLEND_ENVELOPE_EDGES`. A NaN source counts as not wet rather than poisoning the cell, so a
cell can still be reached when one product is missing and the count is then a lower bound - which
biases the weight toward the unblended model rather than toward blending on absent evidence. Any
product being wet puts the envelope at or above the first edge, so a nonzero agreement always
pairs with a nonzero envelope band.
"""
function blend_band_matrix(
    axis::Symbol, y_sat::Matrix{Float64}, sources::Vector{Matrix{Float64}}, threshold::Float64,
)
    axis in BLEND_AXES || throw(ArgumentError("unknown blend axis: $axis"))
    out = zeros(Int, size(y_sat))
    if axis === :constant
        @inbounds for index in eachindex(out)
            value = y_sat[index]
            out[index] = (isnan(value) || value < threshold) ? 0 : 1
        end
        return out
    end
    isempty(sources) && throw(ArgumentError("$axis needs the source products to band on"))
    for source in sources
        size(source) == size(y_sat) ||
            throw(DimensionMismatch("every source product must have the anchor's shape"))
    end
    n_envelope = length(BLEND_ENVELOPE_EDGES)
    @inbounds for index in eachindex(out)
        agreement = 0
        envelope = -Inf
        for source in sources
            value = source[index]
            isnan(value) && continue
            envelope = max(envelope, value)
            value >= threshold && (agreement += 1)
        end
        agreement == 0 && continue
        agreement = min(agreement, 3)
        band = 1
        for (position, edge) in enumerate(BLEND_ENVELOPE_EDGES)
            envelope >= edge && (band = position)
        end
        out[index] = (agreement - 1) * n_envelope + band
    end
    return out
end

"""
Choose one blending weight per band on the inner selection split.

Scores the cells where the source and the fallback are both finite, as the scalar version does, so
the choice is not partly a question of which one dropped the harder cells. Each band's weight is
then minimised over `lambdas` independently, which is exact rather than a shortcut: see
`banded_blend_prediction` for why the bands separate.

Ties go to the smaller weight, matching `select_blend_lambda`: at equal error the answer that
departs less from the fitted model is the one to prefer, and it keeps the choice off the order of
the grid. A band with no scorable cell keeps weight 0, i.e. is left unblended - the same answer a
fold that cannot be helped gets.

Returns `nothing` when no cell is scorable at all, which the caller reports as a skipped fold.
"""
function select_blend_lambdas(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, band::Matrix{Int},
    prediction::Matrix{Float64}, fallback::Matrix{Float64}, lambdas::Vector{Float64},
    n_bands::Int; time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    shared = .!isnan.(y_obs) .& .!isnan.(y_sat) .&
        .!isnan.(prediction) .& .!isnan.(fallback)
    any(shared) || return nothing
    scored_band = ifelse.(shared, band, 0)
    chosen = zeros(Float64, n_bands)
    for target in 1:n_bands
        best_sse = Inf
        best_lambda = 0.0
        for lambda in lambdas
            sse = 0.0
            count = 0
            @inbounds for index in eachindex(scored_band)
                scored_band[index] == target || continue
                difference = (1 - lambda) * prediction[index] + lambda * fallback[index] -
                    y_obs[index]
                sse += difference * difference
                count += 1
            end
            count == 0 && continue
            if sse < best_sse - 1e-12
                best_sse = sse
                best_lambda = lambda
            end
        end
        chosen[target] = best_lambda
    end
    restricted_prediction = ifelse.(shared, prediction, NaN)
    restricted_fallback = ifelse.(shared, fallback, NaN)
    blended = banded_blend_prediction(
        scored_band, restricted_prediction, restricted_fallback, chosen,
    )
    best = _candidate_metrics(y_obs, y_sat, blended; time_weights)
    unblended = _candidate_metrics(y_obs, y_sat, restricted_prediction; time_weights)
    return (;
        lambda=chosen[1], lambdas=chosen, inner_RMSE=best.RMSE, inner_MAE=best.MAE,
        unblended_RMSE=unblended.RMSE, n=best.n, coverage=best.coverage,
        wet_cells=count(index -> scored_band[index] > 0, eachindex(scored_band)),
    )
end

"""
Choose the blending weight on the inner selection split.

Scores every candidate in `lambdas` on the cells where the source and the fallback are both finite,
so the choice is not partly a question of which one dropped the harder cells - the same discipline
`select_auto_method` applies to its contenders.

Ties go to the smaller `lambda`: at equal inner RMSE the answer that departs less from the fitted
model is the one to prefer, and it keeps the choice from depending on the order of the grid.

Returns `nothing` when no cell is scorable, which the caller reports as a skipped fold rather than
defaulting to a weight.
"""
function select_blend_lambda(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64},
    fallback::Matrix{Float64}, lambdas::Vector{Float64}, threshold::Float64;
    time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    shared = .!isnan.(y_obs) .& .!isnan.(y_sat) .&
        .!isnan.(prediction) .& .!isnan.(fallback)
    any(shared) || return nothing
    restricted_prediction = ifelse.(shared, prediction, NaN)
    restricted_fallback = ifelse.(shared, fallback, NaN)
    scored = NamedTuple[]
    for lambda in lambdas
        blended = satellite_wet_blend_prediction(
            y_sat, restricted_prediction, restricted_fallback, lambda, threshold,
        )
        push!(scored, merge((; lambda),
            _candidate_metrics(y_obs, y_sat, blended; time_weights)))
    end
    order = sortperm(scored; by=row -> (row.RMSE, row.MAE, row.lambda))
    best = scored[order[1]]
    unblended = scored[findfirst(row -> row.lambda == 0.0, scored)]
    return (;
        lambda=best.lambda, inner_RMSE=best.RMSE, inner_MAE=best.MAE,
        unblended_RMSE=unblended.RMSE, n=best.n, coverage=best.coverage,
        wet_cells=count(index -> shared[index] && y_sat[index] >= threshold,
            eachindex(shared)),
    )
end

"""
`matrix` moved one hour along `times`: `shift = -1` gives each cell its value at t-1, `+1` at t+1.

A cell whose neighbouring column is not exactly `shift` hours away (a gap in the common grid), or is
NaN there, keeps its own value, so a lagged stack feature degrades to the member it lags.
"""
function stack_hour_shift(matrix::AbstractMatrix{Float64}, times::AbstractVector, shift::Int)
    length(times) == size(matrix, 2) ||
        throw(DimensionMismatch("times must match the matrix's columns"))
    out = Matrix{Float64}(matrix)
    @inbounds for time in axes(matrix, 2)
        other = time + shift
        (1 <= other <= size(matrix, 2) && times[other] - times[time] == Hour(shift)) || continue
        for station in axes(matrix, 1)
            value = matrix[station, other]
            isnan(value) || (out[station, time] = value)
        end
    end
    return out
end

"""The stack's design, one matrix per entry of `STACK_FEATURES` and in that order."""
function stack_feature_matrices(
    mgwr::Matrix{Float64}, adw::Matrix{Float64}, raw::Matrix{Float64}, tps::Matrix{Float64},
    times::AbstractVector,
)
    size(mgwr) == size(adw) == size(raw) == size(tps) ||
        throw(DimensionMismatch("every stack member must have the same shape"))
    return [mgwr, adw, raw, tps,
        stack_hour_shift(mgwr, times, -1), stack_hour_shift(mgwr, times, 1),
        stack_hour_shift(adw, times, -1), stack_hour_shift(adw, times, 1)]
end

"""
Stack band per cell: `agreement + 1 + n_agreement * (distance_band - 1)`, in `1:stack_band_count()`.

`agreement` is `blend_band_matrix(:agreement_envelope, ...)` for these stations, whose band 0 (no
product wet) the stack fits like any other. `distance_km` is each station's distance to the nearest
station the predicting model was trained on, which is knowable at prediction time.
"""
function stack_band_matrix(agreement::AbstractMatrix{Int}, distance_km::AbstractVector{<:Real})
    length(distance_km) == size(agreement, 1) ||
        throw(DimensionMismatch("need one distance per station row"))
    n_agreement = blend_band_count(:agreement_envelope) + 1
    distance_band = [1 + count(edge -> distance >= edge, STACK_DISTANCE_EDGES)
                     for distance in distance_km]
    return agreement .+ 1 .+ n_agreement .* (distance_band .- 1)
end

"""
Non-negative least squares from normal equations, by enumerating the support.

Exact for the eight or so features the stack has (255 subsets of 8×8 solves). A singular subset is
skipped; if every subset is singular the zero vector comes back, which the caller never uses as
the fitted answer because the pooled system over every band is not singular in practice.
"""
function _nnls_normal(gram::Matrix{Float64}, rhs::Vector{Float64}, yy::Float64)
    p = length(rhs)
    best_sse, best = Inf, zeros(Float64, p)
    for support in 1:(2^p - 1)
        active = [j for j in 1:p if (support >> (j - 1)) & 1 == 1]
        solution = try
            cholesky(Symmetric(gram[active, active])) \ rhs[active]
        catch
            continue
        end
        all(>=(0), solution) || continue
        weights = zeros(Float64, p)
        weights[active] = solution
        sse = yy - 2 * dot(weights, rhs) + dot(weights, gram * weights)
        sse < best_sse && ((best_sse, best) = (sse, weights))
    end
    return best
end

"""
Fit one non-negative weight vector per stack band on the inner selection split.

Scores the cells where the gauge and every feature are finite, so no member's gaps decide which
cells a weight is fitted on. A band with fewer than `min_cells` such cells takes the pooled fit
over every band and is flagged in `fell_back`.

Returns `nothing` when no cell is usable, which the caller reports as a skipped fold. Otherwise
`weights` is `length(features) x n_bands`, beside the inner RMSE of the stack and of `features[1]`
alone on the same cells.
"""
function select_stack_weights(
    y_obs::Matrix{Float64}, features::Vector{Matrix{Float64}}, band::Matrix{Int}, n_bands::Int;
    min_cells::Int=STACK_MIN_BAND_CELLS,
)
    for feature in features
        size(feature) == size(y_obs) == size(band) || throw(DimensionMismatch(
            "every feature and the band matrix must have the observation's shape",
        ))
    end
    p = length(features)
    gram = zeros(Float64, p, p, n_bands)
    rhs = zeros(Float64, p, n_bands)
    yy = zeros(Float64, n_bands)
    used = zeros(Int, n_bands)
    row = zeros(Float64, p)
    @inbounds for index in eachindex(y_obs)
        observed = y_obs[index]
        isnan(observed) && continue
        usable = true
        for (j, feature) in enumerate(features)
            value = feature[index]
            isnan(value) && (usable = false; break)
            row[j] = value
        end
        usable || continue
        target = band[index]
        1 <= target <= n_bands || throw(ArgumentError("band $target is outside 1:$n_bands"))
        used[target] += 1
        yy[target] += observed * observed
        for i in 1:p
            rhs[i, target] += row[i] * observed
            for j in i:p
                gram[i, j, target] += row[i] * row[j]
            end
        end
    end
    sum(used) == 0 && return nothing
    for target in 1:n_bands, i in 1:p, j in 1:(i - 1)
        gram[i, j, target] = gram[j, i, target]
    end
    pooled = _nnls_normal(dropdims(sum(gram; dims=3); dims=3), vec(sum(rhs; dims=2)), sum(yy))
    weights = zeros(Float64, p, n_bands)
    fell_back = falses(n_bands)
    sse_stack, sse_source = 0.0, 0.0
    for target in 1:n_bands
        fell_back[target] = used[target] < min_cells
        weights[:, target] = fell_back[target] ? pooled :
            _nnls_normal(gram[:, :, target], rhs[:, target], yy[target])
        used[target] == 0 && continue
        w = weights[:, target]
        sse_stack += yy[target] - 2 * dot(w, rhs[:, target]) + dot(w, gram[:, :, target] * w)
        sse_source += yy[target] - 2 * rhs[1, target] + gram[1, 1, target]
    end
    n = sum(used)
    return (; weights, used, fell_back, n,
        inner_RMSE=sqrt(max(sse_stack, 0.0) / n), source_RMSE=sqrt(max(sse_source, 0.0) / n))
end

"""
Apply per-band stack weights. A cell with any feature missing keeps `features[1]`'s value - the
source prediction, NaN if the source itself is - so a gap in a secondary member can neither shrink
the stack's coverage below its source's nor poison a cell, the rule the blends follow too.
"""
function stacked_prediction(
    features::Vector{Matrix{Float64}}, band::Matrix{Int}, weights::Matrix{Float64},
)
    size(weights, 1) == length(features) ||
        throw(DimensionMismatch("need one weight per feature"))
    source = first(features)
    out = similar(source)
    @inbounds for index in eachindex(source)
        total = 0.0
        complete = true
        for (j, feature) in enumerate(features)
            value = feature[index]
            isnan(value) && (complete = false; break)
            total += weights[j, band[index]] * value
        end
        out[index] = complete ? total : source[index]
    end
    return out
end
