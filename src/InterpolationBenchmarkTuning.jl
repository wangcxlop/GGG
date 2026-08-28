"""
Score one tuning candidate.

`time_weights`, when given, holds one weight per column of `y_obs` and makes RMSE/MAE weighted
means, so a stratified tuning subsample (see [`_tuning_time_sample`](@ref)) estimates the metric
over the full hour set rather than over the subsample's own wet-heavy mixture. `n` and
`coverage` stay raw cell counts either way: they only gate whether the candidate produced
predictions at all (`min_tuning_coverage`), which reweighting would silently redefine.

`require_satellite` is now `true` for direct candidates too, where it used to be
`mode == "residual"`. The reported metric is computed on `_common_method_mask`, which includes
`raw` (the satellite field) among `MASK_METHODS` — so a cell with a missing satellite value is
never scored, whatever method produced it. Tuning a direct method over the wider obs-only mask
therefore selected against a population the results table does not contain, the same
selection/reporting estimand mismatch that `tuning_geometry` fixed in the spatial dimension.
Holding it fixed across modes is also what lets `auto` compare a direct candidate against a
residual one at all (see `AUTO_CANDIDATE_RUNS`).
"""
function _candidate_metrics(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64};
    require_satellite::Bool=true, time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    base_mask = require_satellite ? (.!isnan.(y_obs) .& .!isnan.(y_sat)) : .!isnan.(y_obs)
    mask = base_mask .& .!isnan.(prediction)
    base_n = count(base_mask)
    n = count(mask)
    n == 0 && return (; n=0, coverage=0.0, RMSE=Inf, MAE=Inf)
    if time_weights === nothing
        metric = metric_continuous(y_obs, prediction; mask=mask)
        return (; n, coverage=n / base_n, RMSE=metric.RMSE, MAE=metric.MAE)
    end
    length(time_weights) == size(y_obs, 2) || throw(ArgumentError(
        "time_weights must have one entry per tuning hour",
    ))
    weighted_square = 0.0
    weighted_absolute = 0.0
    weight_total = 0.0
    for time in axes(y_obs, 2)
        weight = time_weights[time]
        for station in axes(y_obs, 1)
            mask[station, time] || continue
            error = prediction[station, time] - y_obs[station, time]
            weighted_square += weight * error^2
            weighted_absolute += weight * abs(error)
            weight_total += weight
        end
    end
    weight_total > 0 || return (; n, coverage=n / base_n, RMSE=Inf, MAE=Inf)
    return (; n, coverage=n / base_n,
        RMSE=sqrt(weighted_square / weight_total), MAE=weighted_absolute / weight_total)
end

function _scan_row(; scheme, product, fold, mode, method, group="all", iteration=0,
    repeat=1, seed=0,
    power=NaN, neighbors=missing,
    smooth=NaN, kernel=missing, adaptive=missing, bw=NaN, shrink=NaN, n=0, coverage=0.0,
    RMSE=Inf, MAE=Inf, status="success", error="", selected=false)
    stored_neighbors = neighbors === nothing ? 0 : neighbors
    return (;
        scheme=string(scheme), product=String(product), fold=Int(fold), mode=String(mode),
        method=String(method), group=String(group), iteration=Int(iteration),
        repeat=Int(repeat), seed=Int(seed),
        power=Float64(power), neighbors=stored_neighbors, smooth=Float64(smooth),
        kernel, adaptive, bw=Float64(bw), shrink=Float64(shrink),
        n=Int(n), coverage=Float64(coverage),
        RMSE=Float64(RMSE), MAE=Float64(MAE), status=String(status), error=String(error),
        selected=Bool(selected),
    )
end

function _select_candidate!(rows::Vector{NamedTuple}, first_row::Int, min_coverage::Float64)
    valid = [i for i in first_row:length(rows) if
        rows[i].status == "success" && rows[i].coverage >= min_coverage && isfinite(rows[i].RMSE)]
    isempty(valid) && throw(ArgumentError("all parameter candidates failed or had insufficient coverage"))
    best_index = sort(valid; by=i -> (rows[i].RMSE, rows[i].MAE, -rows[i].coverage))[1]
    rows[best_index] = merge(rows[best_index], (; selected=true))
    return rows[best_index]
end

function select_mgwr_bandwidths!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    scheme::Symbol, product::String, fold::Int, train_lonlat::Matrix{Float64},
    residuals::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    time_weights::Union{Nothing,Vector{Float64}}=nothing, selection_groups=nothing,
)
    n_train, n_time = size(residuals)
    designs = build_mgwr_designs(train_lonlat, train_lonlat)
    first_row = length(scan_rows) + 1
    # Each (kernel, adaptive/fixed bandwidth-family) combination runs its own independent
    # per-group coordinate descent to convergence (identical logic to the historical
    # bisquare-only, adaptive-only search, just parameterised); combinations are then compared
    # by the RMSE of the last group scored in their converged sweep.
    #
    # That is a whole-configuration score, not a per-group one, so the comparison is like for
    # like. The descent only stops once a full sweep leaves `bandwidths` unchanged, which means
    # every group in that sweep was scored with all the others already sitting at their converged
    # values — so all of a combination's winning rows in the final iteration carry the identical
    # RMSE, and taking the last is the same as taking any. (An earlier version of this comment
    # claimed the descent "never computes one joint RMSE" and that this was merely the best number
    # available; that was wrong, and `select_joint_parameter!` states the correct reasoning.)
    #
    # A combination that has no candidates, or fails to converge, is skipped rather than aborting
    # the search — and a non-converged sweep is exactly the case where the last-group row would
    # *not* be a joint score, which is why `converged || continue` guards the comparison.
    best_result = nothing
    bandwidth_families = _bandwidth_families(
        cfg, _usable_adaptive_candidates(cfg.mger.bw_adaptive, n_train, selection_groups),
    )
    for kernel in cfg.mger.kernels
        for (adaptive, raw_candidates) in bandwidth_families
            candidates = adaptive ? sort(unique(round.(raw_candidates))) : sort(unique(raw_candidates))
            isempty(candidates) && continue
            # Start every group at the widest candidate - the global one, when it is in the
            # family. A back-fit that starts wide and tightens is the conventional direction,
            # and the descent runs to convergence either way.
            bandwidths = fill(last(candidates), length(designs.group_names))
            final_iteration = 0
            last_index = 0
            converged = false
            try
                for iteration in 1:cfg.mgwr_max_tuning_iterations
                    previous = copy(bandwidths)
                    for (group_index, group) in enumerate(designs.group_names)
                        candidate_rows = Int[]
                        for bw in candidates
                            trial = copy(bandwidths)
                            trial[group_index] = bw
                            try
                                interpolated = selection_groups === nothing ?
                                    _mgwr_predict(
                                        train_lonlat, residuals, train_lonlat;
                                        bandwidths=trial, kernel, adaptive, exclude_self=true,
                                    ) :
                                    _selection_oof(selection_groups, n_train, n_time, (tr, va) ->
                                        _mgwr_predict(
                                            train_lonlat[tr, :], residuals[tr, :],
                                            train_lonlat[va, :]; bandwidths=trial, kernel, adaptive,
                                        ))
                                prediction = max.(y_sat .+ interpolated, 0.0)
                                metrics = _candidate_metrics(y_obs, y_sat, prediction; time_weights)
                                push!(scan_rows, _scan_row(;
                                    scheme, product, fold, mode="residual", method="mgwr",
                                    group, iteration, kernel, adaptive, bw, metrics...,
                                ))
                            catch e
                                push!(scan_rows, _scan_row(;
                                    scheme, product, fold, mode="residual", method="mgwr",
                                    group, iteration, kernel, adaptive, bw,
                                    status="failed", error=sprint(showerror, e),
                                ))
                            end
                            push!(candidate_rows, length(scan_rows))
                        end
                        valid = filter(candidate_rows) do index
                            row = scan_rows[index]
                            row.status == "success" && row.coverage >= cfg.min_tuning_coverage &&
                                isfinite(row.RMSE)
                        end
                        isempty(valid) && throw(ArgumentError(
                            "all MGWR candidates failed for group $group " *
                            "(kernel=$kernel, adaptive=$adaptive)",
                        ))
                        best_index = sort(valid; by=index -> (
                            scan_rows[index].RMSE, scan_rows[index].MAE, -scan_rows[index].coverage,
                        ))[1]
                        bandwidths[group_index] = scan_rows[best_index].bw
                        last_index = best_index
                    end
                    final_iteration = iteration
                    if bandwidths == previous
                        converged = true
                        break
                    end
                end
            catch error
                error isa ArgumentError || rethrow()
                continue
            end
            converged || continue
            rmse = scan_rows[last_index].RMSE
            if best_result === nothing || rmse < best_result.rmse
                best_result = (; kernel, adaptive, bandwidths=copy(bandwidths), final_iteration, rmse)
            end
        end
    end
    best_result === nothing && throw(ArgumentError(
        "MGWR bandwidth search did not converge for any kernel/bandwidth-family combination " *
        "in $(cfg.mgwr_max_tuning_iterations) iterations",
    ))
    kernel = best_result.kernel
    adaptive = best_result.adaptive
    bandwidths = best_result.bandwidths
    final_iteration = best_result.final_iteration

    for (group_index, group) in enumerate(designs.group_names)
        matching_rows = filter(first_row:length(scan_rows)) do index
            row = scan_rows[index]
            row.method == "mgwr" && row.group == group && row.kernel == kernel &&
                row.adaptive == adaptive && row.iteration == final_iteration &&
                row.status == "success" && row.bw == bandwidths[group_index]
        end
        isempty(matching_rows) && error("MGWR selected bandwidth row is missing for $group")
        selected_index = last(matching_rows)
        scan_rows[selected_index] = merge(scan_rows[selected_index], (; selected=true))
    end
    return (; bandwidths, kernel, adaptive)
end

"""
Two-stratum tuning-hour sample and its Horvitz-Thompson weights.

Hyperparameters are chosen on this subsample but the benchmark reports pooled RMSE over every
hour, and the two populations are nothing alike: on the 11426 common hours the wettest half of
the subsample drives the tuning set to 40% wet cells against 7% in the reported population, and
13x its mean squared observation. An unweighted RMSE over the subsample therefore scores
wet-hour skill while the results table scores mostly-dry-hour skill.

So the wettest `maximum_times ÷ 8` hours are taken with certainty and the rest of the hours are
systematically subsampled and up-weighted by the inverse of their inclusion probability. The
weighted MSE is then an unbiased estimator of the pooled MSE over all `size(y_obs, 2)` hours.

The certainty stratum is an eighth of the budget rather than the historical half because under
weighting it only carries its true population share of the weight (~1.5%); its job is to cover
the extreme-error tail, and every hour beyond that is budget taken away from the weighted
remainder, which is where almost all of the weight — and therefore the estimator's variance —
actually sits. Measured over 12 product/fold/method cells against the exact full-record curve,
an eighth reproduces the legacy objective's candidate choice in 12/12 cells at 0.0% regret while
cutting the estimator's bias from +197% to -2.5%; a half costs 0.7-0.9% mean regret, and
dropping the certainty stratum entirely costs 0.4%. See
`scripts/verify_tuning_time_weighting.jl`.

`weighting=:uniform` returns the historical wet-heavy index set with no weights, and is the
default. What this sampler corrects is real but second-order:

!!! note "The temporal mismatch is not what dominates selection error"
    Scored against the exact full-record curve, the unweighted subsample's RMSE is +201% too
    high — yet that error is nearly a constant multiplier, so it cancels in the *ranking* and
    costs only ~0.1% regret. Correcting the level trades that harmless bias for real variance
    at a fixed hour budget and costs ~0.6%.

    What dominates is a different mismatch entirely, and it is spatial. Candidates are scored by
    leave-one-out at *training* stations, where the held-out station sits inside its own
    neighbourhood, but the benchmark reports RMSE at *fold* stations 20+ km from any training
    gauge. On the joint-covariate path the held-out RMSE curve falls monotonically from bw=8 to
    bw=160 while both tuning curves are U-shaped with a minimum at 12-30: over most of the grid
    the two criteria are anti-correlated, and the tuner's pick costs +19% to +101%. This is why
    widening the grid (`--local-grid`, floor 30 -> 8) made `residual_gwr`/`mixed_gwr` worse — the
    old floor was accidentally shielding the tuner from its own preference. Fixing that needs an
    inner *spatial* split of the training fold, not a reweighted hour sample.
"""
function _tuning_time_sample(
    y_obs::Matrix{Float64}, maximum_times::Int, weighting::Symbol=:stratified,
)
    total = size(y_obs, 2)
    (maximum_times <= 0 || total <= maximum_times) && return collect(1:total), nothing
    wetness = [let values = filter(isfinite, @view(y_obs[:, time]))
        isempty(values) ? -Inf : mean(values)
    end for time in 1:total]
    if weighting === :uniform
        wet_count = div(maximum_times, 2)
        wet = partialsortperm(wetness, 1:wet_count; rev=true)
        spaced = unique(round.(Int, range(1, total; length=maximum_times - wet_count)))
        return sort(unique(vcat(wet, spaced))), nothing
    end
    wet_count = max(1, div(maximum_times, 8))
    wet = partialsortperm(wetness, 1:wet_count; rev=true)
    # `rest` excludes the certainty stratum, so the two strata are disjoint and the sample size
    # is exactly `maximum_times` (the historical union could collapse to fewer hours).
    rest = setdiff(1:total, wet)
    dry_count = min(maximum_times - wet_count, length(rest))
    # Sample the remainder systematically along the wetness ordering rather than the calendar
    # ordering. Inclusion probability is uniform either way, so the weights stay valid, but
    # spanning the wetness distribution stops the estimate from swinging on how many wet hours
    # that fell outside the certainty stratum happen to be picked up at ~67x weight. Measured
    # over 40 synthetic replicates this cuts the mean error against the full-sample RMSE from
    # 6.3% to 2.2%, and the worst case from 19.8% to 2.2%.
    frame = rest[sortperm(wetness[rest])]
    dry = frame[unique(round.(Int, range(1, length(frame); length=dry_count)))]
    indices = vcat(wet, dry)
    weights = vcat(ones(Float64, length(wet)), fill(length(rest) / length(dry), length(dry)))
    order = sortperm(indices)
    return indices[order], weights[order]
end

_tuning_time_indices(y_obs::Matrix{Float64}, maximum_times::Int) =
    first(_tuning_time_sample(y_obs, maximum_times, :uniform))

function select_interpolation_parameter!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    method::String, mode::String, scheme::Symbol, product::String, fold::Int,
    train_lonlat::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    dem_context=nothing, joint_context=nothing, hurdle_context=nothing,
    selection_groups=nothing, joint_selection_contexts=nothing, repeat_seed::Int=cfg.seed,
    repeat_index::Int=1,
)
    (mode, method) in BENCHMARK_RUNS ||
        throw(ArgumentError("unsupported benchmark method/mode pair: $method/$mode"))
    if dem_context !== nothing && mode == "residual" && method in ("gwr", "mixed_gwr", "mgwr")
        return select_dem_parameter!(
            scan_rows, cfg, method, mode, scheme, product, fold, dem_context,
        )
    end
    # One inner split per fold, shared by every method, so candidates from different methods are
    # scored against the same held-out stations. The caller normally supplies it (it also needs
    # it to build the joint inner contexts); computing it here keeps direct callers working.
    if selection_groups === nothing && cfg.tuning_geometry === :inner_spatial
        selection_groups = selection_folds(
            cfg, scheme, string.(1:size(train_lonlat, 1)), train_lonlat, fold, repeat_seed;
            repeat_index,
        )
    end
    if joint_context !== nothing && mode == "residual" && method in ("gwr", "mixed_gwr", "mgwr")
        return select_joint_parameter!(
            scan_rows, cfg, method, mode, scheme, product, fold,
            joint_context, y_obs, y_sat; joint_selection_contexts,
        )
    end
    # `hurdle_gwr` is the one method that needs the timestamps themselves, so the tuning subset
    # has to be applied to them as well as to the data matrices.
    tuning_times = hurdle_context === nothing ? nothing : hurdle_context.times
    time_weights = nothing
    if cfg.tuning_max_times > 0 && size(y_obs, 2) > cfg.tuning_max_times
        time_idx, time_weights = _tuning_time_sample(
            y_obs, cfg.tuning_max_times, cfg.tuning_time_weighting,
        )
        y_obs = y_obs[:, time_idx]
        y_sat = y_sat[:, time_idx]
        tuning_times = tuning_times === nothing ? nothing : tuning_times[time_idx]
    end
    target_values = mode == "direct" ? y_obs : y_obs .- y_sat
    n_train, n_time = size(target_values)
    # `:loocv` keeps the historical criterion — predict at every training station with only that
    # station excluded. `:inner_spatial` predicts out-of-fold onto whole spatial groups instead,
    # which is the geometry the benchmark reports. Both produce an `n_train × n_time` matrix, so
    # everything downstream is identical.
    scored(predict_pair, predict_loocv) = selection_groups === nothing ?
        predict_loocv() : _selection_oof(selection_groups, n_train, n_time, predict_pair)
    # Shared by every GWR-family branch below, so `gwr`, `hurdle_gwr` and `mixed_gwr` search the
    # same candidates as each other and as the joint models in `select_joint_parameter!`.
    bandwidth_families = _bandwidth_families(
        cfg, _usable_adaptive_candidates(cfg.mger.bw_adaptive, n_train, selection_groups),
    )
    first_row = length(scan_rows) + 1
    if method in ("idw", "adw")
        predictor = method == "idw" ? idw_predict : adw_predict
        for power in cfg.idw_powers, neighbors in cfg.neighbor_candidates
            try
                interpolated = scored(
                    (tr, va) -> predictor(
                        train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                        power=power, neighbors=neighbors,
                    ),
                    () -> predictor(
                        train_lonlat, target_values, train_lonlat;
                        power=power, neighbors=neighbors, exclude_self=true,
                    ),
                )
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(
                    y_obs, y_sat, prediction;
                    require_satellite=true, time_weights,
                )
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, power, neighbors, metrics...,
                ))
            catch e
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, power, neighbors,
                    status="failed", error=sprint(showerror, e),
                ))
            end
        end
    elseif method == "tps"
        for smooth in cfg.tps_smooth_candidates
            try
                interpolated = scored(
                    (tr, va) -> tps_predict(
                        train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                        smooth=smooth,
                    ),
                    () -> tps_loo_predict(train_lonlat, target_values; smooth=smooth),
                )
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(
                    y_obs, y_sat, prediction;
                    require_satellite=true, time_weights,
                )
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, smooth, metrics...,
                ))
            catch e
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, smooth,
                    status="failed", error=sprint(showerror, e),
                ))
            end
        end
    elseif method == "gwr"
        for kernel in cfg.mger.kernels
            for (adaptive, candidates) in bandwidth_families
                for bw in candidates
                    try
                        interpolated = scored(
                            (tr, va) -> _gwr_predict(
                                train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                                kernel, adaptive, bw,
                            ),
                            () -> _gwr_predict(
                                train_lonlat, target_values, train_lonlat;
                                kernel, adaptive, bw, exclude_self=true,
                            ),
                        )
                        prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                        metrics = _candidate_metrics(
                            y_obs, y_sat, prediction;
                            require_satellite=true, time_weights,
                        )
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw, metrics...,
                        ))
                    catch e
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw,
                            status="failed", error=sprint(showerror, e),
                        ))
                    end
                end
            end
        end
    elseif method == "hurdle_gwr"
        mode == "direct" || throw(ArgumentError("hurdle_gwr only supports direct mode"))
        hurdle_context === nothing &&
            throw(ArgumentError("hurdle_gwr requires a hurdle_context carrying the fold times"))
        for kernel in cfg.mger.kernels
            for (adaptive, candidates) in bandwidth_families
                for bw in candidates
                    # `FixedBandwidth` rejects a non-finite bandwidth by construction, so the
                    # family's global candidate cannot be expressed as a `HurdleGWRConfig`.
                    # Skip it rather than logging a guaranteed-failure row per kernel; giving
                    # `hurdle_gwr` a global option means teaching `SpatialBandwidth` about one.
                    isfinite(bw) || continue
                    try
                        hurdle_config = _hurdle_config(hurdle_context, kernel, adaptive, bw)
                        corrected = scored(
                            (tr, va) -> _hurdle_predict(
                                train_lonlat[tr, :], y_obs[tr, :], y_sat[tr, :],
                                train_lonlat[va, :], y_sat[va, :], something(tuning_times);
                                config=hurdle_config,
                            ).corrected,
                            () -> _hurdle_predict(
                                train_lonlat, y_obs, y_sat, train_lonlat, y_sat,
                                something(tuning_times);
                                config=hurdle_config, exclude_self=true,
                            ).corrected,
                        )
                        metrics = _candidate_metrics(y_obs, y_sat, corrected; time_weights)
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw, metrics...,
                        ))
                    catch e
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw,
                            status="failed", error=sprint(showerror, e),
                        ))
                    end
                end
            end
        end
    elseif method == "mixed_gwr"
        mode == "residual" || throw(ArgumentError("mixed_gwr only supports residual mode"))
        for kernel in cfg.mger.kernels
            for (adaptive, candidates) in bandwidth_families
                for bw in candidates
                    try
                        interpolated = scored(
                            (tr, va) -> _mixed_gwr_predict(
                                train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                                bw, kernel, adaptive,
                            ),
                            () -> _mixed_gwr_predict(
                                train_lonlat, target_values, train_lonlat;
                                bw, kernel, adaptive, exclude_self=true,
                            ),
                        )
                        prediction = max.(y_sat .+ interpolated, 0.0)
                        metrics = _candidate_metrics(y_obs, y_sat, prediction; time_weights)
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw, metrics...,
                        ))
                    catch e
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw,
                            status="failed", error=sprint(showerror, e),
                        ))
                    end
                end
            end
        end
    elseif method == "mgwr"
        mode == "residual" || throw(ArgumentError("mgwr only supports residual mode"))
        return select_mgwr_bandwidths!(
            scan_rows, cfg, scheme, product, fold, train_lonlat,
            target_values, y_obs, y_sat; time_weights, selection_groups,
        )
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return _select_candidate!(scan_rows, first_row, cfg.min_tuning_coverage)
end

"""
Pick one GWR-family method for a fold, using only the inner selection split.

`contenders` is one `(; method, prediction)` per candidate, where `prediction` is that method's
already-selected hyperparameters re-predicted across the inner split by
[`inner_selection_prediction`](@ref) — so every entry is an out-of-fold prediction over the same
training stations and the same tuning hours.

The comparison is made on the intersection of the contenders' masks. Each is scored by the same
`_candidate_metrics` the scan uses, after blanking every cell outside the shared mask, so a method
cannot come out ahead by having failed on the cells it found hardest. `shared_mask_coverage` on
the returned row is the share of scoreable cells that survived that intersection; a low value
means some contender was patchy and the comparison rests on less than the full record. It is
not the same quantity as `run_status.csv`'s `prediction_coverage`, which is measured on the
outer held-out stations.

Returns `nothing` when there is no shared mask to score on — the caller then leaves `auto`
unpredicted for the fold rather than guessing, and `run_status.csv` records it.
"""
function select_auto_method(
    contenders::Vector, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    isempty(contenders) && return nothing
    shared = .!isnan.(y_obs) .& .!isnan.(y_sat)
    for contender in contenders
        shared = shared .& .!isnan.(contender.prediction)
    end
    any(shared) || return nothing
    scored = [merge((; contender.method), _candidate_metrics(
        y_obs, y_sat, ifelse.(shared, contender.prediction, NaN); time_weights,
    )) for contender in contenders]
    # Method name breaks ties so a fold's choice does not depend on `BENCHMARK_RUNS` ordering.
    order = sortperm(scored; by=row -> (row.RMSE, row.MAE, row.method))
    best = scored[order[1]]
    runner_up = length(order) > 1 ? scored[order[2]] : nothing
    return (;
        chosen=best.method, chosen_rmse=best.RMSE, chosen_mae=best.MAE,
        runner_up=runner_up === nothing ? "" : runner_up.method,
        runner_up_rmse=runner_up === nothing ? NaN : runner_up.RMSE,
        n=best.n, shared_mask_coverage=best.coverage, n_contenders=length(contenders),
    )
end
