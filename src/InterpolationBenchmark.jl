using CSV, DataFrames, Dates, LinearAlgebra, Random, Statistics
using MixedGWR

if !isdefined(@__MODULE__, :MGERConfig)
    include(joinpath(@__DIR__, "MGERPipeline.jl"))
end
if !isdefined(@__MODULE__, :TraditionalInterpolation)
    include(joinpath(@__DIR__, "TraditionalInterpolation.jl"))
end
using .TraditionalInterpolation

const BENCHMARK_METHODS = [
    "raw", "idw", "adw", "tps", "gwr",
    "residual_idw", "residual_adw", "residual_tps", "residual_gwr",
]
const TRADITIONAL_RESIDUAL_METHODS = ["residual_idw", "residual_adw", "residual_tps"]

Base.@kwdef struct InterpolationBenchmarkConfig
    mger::MGERConfig
    k::Int = 5
    seed::Int = 20260627
    cv_schemes::Vector{Symbol} = [:balanced_spatial, :random]
    idw_powers::Vector{Float64} = [1.0, 1.5, 2.0, 2.5, 3.0]
    neighbor_candidates::Vector{Union{Nothing,Int}} = Union{Nothing,Int}[8, 16, 32, 64, nothing]
    tps_smooth_candidates::Vector{Float64} = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    min_tuning_coverage::Float64 = 0.95
    tuning_max_times::Int = 336
    event_thresholds::Vector{Float64} = [0.1, 2.5, 8.0, 16.0]
    bootstrap_reps::Int = 2000
end

function _validate_benchmark_config(cfg::InterpolationBenchmarkConfig, n_station::Int)
    2 <= cfg.k <= n_station || throw(ArgumentError("k must be between 2 and station count"))
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
    return nothing
end

function _initial_spatial_centers(xy::Matrix{Float64}, k::Int, rng::AbstractRNG)
    n = size(xy, 1)
    centers = Matrix{Float64}(undef, k, 2)
    first_index = rand(rng, 1:n)
    centers[1, :] = xy[first_index, :]
    nearest2 = fill(Inf, n)
    for center_index in 2:k
        previous = @view centers[center_index - 1, :]
        for station in 1:n
            d2 = sum(abs2, @view(xy[station, :]) .- previous)
            nearest2[station] = min(nearest2[station], d2)
        end
        next_index = argmax(nearest2)
        centers[center_index, :] = xy[next_index, :]
    end
    return centers
end

function _capacity_assignment(xy::Matrix{Float64}, centers::Matrix{Float64}, capacities::Vector{Int})
    n, k = size(xy, 1), size(centers, 1)
    distances = Matrix{Float64}(undef, n, k)
    for station in 1:n, cluster in 1:k
        distances[station, cluster] = sum(abs2, @view(xy[station, :]) .- @view(centers[cluster, :]))
    end
    certainty = Vector{Float64}(undef, n)
    for station in 1:n
        ordered = partialsort(@view(distances[station, :]), 1:min(2, k))
        certainty[station] = length(ordered) == 1 ? Inf : ordered[2] - ordered[1]
    end
    order = sortperm(1:n; by=station -> (-certainty[station], station))
    remaining = copy(capacities)
    assignments = zeros(Int, n)
    for station in order
        cluster_order = sortperm(1:k; by=cluster -> (distances[station, cluster], cluster))
        chosen = findfirst(cluster -> remaining[cluster] > 0, cluster_order)
        chosen === nothing && error("capacity assignment failed")
        cluster = cluster_order[chosen]
        assignments[station] = cluster
        remaining[cluster] -= 1
    end
    return assignments
end

"""Deterministic, compact, capacity-balanced two-dimensional spatial folds."""
function split_stations_balanced_spatial_kfold(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, seed::Int=20260627,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat, 1) == n || throw(DimensionMismatch("lonlat rows must match station_ids"))
    xy = local_km_coordinates(lonlat)
    capacities = fill(div(n, k), k)
    capacities[1:rem(n, k)] .+= 1
    centers = _initial_spatial_centers(xy, k, MersenneTwister(seed))
    assignments = zeros(Int, n)
    for _ in 1:100
        updated = _capacity_assignment(xy, centers, capacities)
        updated == assignments && break
        assignments = updated
        for cluster in 1:k
            idx = findall(==(cluster), assignments)
            centers[cluster, :] = vec(mean(xy[idx, :], dims=1))
        end
    end
    folds = [String[] for _ in 1:k]
    for station in 1:n
        push!(folds[assignments[station]], station_ids[station])
    end
    return folds
end

function benchmark_folds(
    scheme::Symbol, ids::Vector{String}, lonlat::Matrix{Float64}; k::Int, seed::Int,
)
    if scheme == :balanced_spatial
        return split_stations_balanced_spatial_kfold(ids, lonlat; k=k, seed=seed)
    elseif scheme == :random
        return split_stations_kfold(ids; k=k, rng=MersenneTwister(seed))
    elseif scheme == :strip
        return split_stations_spatial_block_kfold(ids, lonlat; k=k)
    end
    throw(ArgumentError("unsupported CV scheme: $scheme"))
end

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
    return st_gwr_predict_nanaware(X_train, values, weights; Xpred=X_target, min_obs=3)
end

function _candidate_metrics(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64};
    require_satellite::Bool=true,
)
    base_mask = require_satellite ? (.!isnan.(y_obs) .& .!isnan.(y_sat)) : .!isnan.(y_obs)
    mask = base_mask .& .!isnan.(prediction)
    base_n = count(base_mask)
    n = count(mask)
    n == 0 && return (; n=0, coverage=0.0, RMSE=Inf, MAE=Inf)
    metric = metric_continuous(y_obs, prediction; mask=mask)
    return (; n, coverage=n / base_n, RMSE=metric.RMSE, MAE=metric.MAE)
end

function _scan_row(; scheme, product, fold, mode, method, power=NaN, neighbors=missing,
    smooth=NaN, kernel=missing, adaptive=missing, bw=NaN, n=0, coverage=0.0,
    RMSE=Inf, MAE=Inf, status="success", error="")
    stored_neighbors = neighbors === nothing ? 0 : neighbors
    return (;
        scheme=string(scheme), product=String(product), fold=Int(fold), mode=String(mode),
        method=String(method), power=Float64(power), neighbors=stored_neighbors, smooth=Float64(smooth),
        kernel, adaptive, bw=Float64(bw), n=Int(n), coverage=Float64(coverage),
        RMSE=Float64(RMSE), MAE=Float64(MAE), status=String(status), error=String(error),
        selected=false,
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

function select_interpolation_parameter!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    method::String, mode::String, scheme::Symbol, product::String, fold::Int,
    train_lonlat::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64},
)
    if cfg.tuning_max_times > 0 && size(y_obs, 2) > cfg.tuning_max_times
        wetness = [let values = filter(isfinite, @view(y_obs[:, time]))
            isempty(values) ? -Inf : mean(values)
        end for time in axes(y_obs, 2)]
        n_wet = div(cfg.tuning_max_times, 2)
        wet_idx = partialsortperm(wetness, 1:n_wet; rev=true)
        evenly_spaced = unique(round.(Int, range(1, size(y_obs, 2); length=cfg.tuning_max_times - n_wet)))
        time_idx = sort(unique(vcat(wet_idx, evenly_spaced)))
        y_obs = y_obs[:, time_idx]
        y_sat = y_sat[:, time_idx]
    end
    target_values = mode == "direct" ? y_obs : y_obs .- y_sat
    first_row = length(scan_rows) + 1
    if method in ("idw", "adw")
        predictor = method == "idw" ? idw_predict : adw_predict
        for power in cfg.idw_powers, neighbors in cfg.neighbor_candidates
            try
                interpolated = predictor(
                    train_lonlat, target_values, train_lonlat;
                    power=power, neighbors=neighbors, exclude_self=true,
                )
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(y_obs, y_sat, prediction; require_satellite=mode == "residual")
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
                interpolated = tps_loo_predict(train_lonlat, target_values; smooth=smooth)
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(y_obs, y_sat, prediction; require_satellite=mode == "residual")
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
            for (adaptive, candidates) in ((true, cfg.mger.bw_adaptive), (false, cfg.mger.bw_fixed_km))
                for bw in candidates
                    try
                        interpolated = _gwr_predict(
                            train_lonlat, target_values, train_lonlat;
                            kernel, adaptive, bw, exclude_self=true,
                        )
                        prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                        metrics = _candidate_metrics(y_obs, y_sat, prediction; require_satellite=mode == "residual")
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
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return _select_candidate!(scan_rows, first_row, cfg.min_tuning_coverage)
end

function predict_selected(
    selected, method::String, mode::String,
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64}, y_sat_target::Matrix{Float64},
)
    values = mode == "direct" ? y_obs_train : y_obs_train .- y_sat_train
    selected_neighbors = ismissing(selected.neighbors) ? nothing :
        (selected.neighbors == 0 ? nothing : Int(selected.neighbors))
    interpolated = if method == "idw"
        idw_predict(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "adw"
        adw_predict(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "tps"
        tps_predict(train_lonlat, values, target_lonlat; smooth=selected.smooth)
    elseif method == "gwr"
        _gwr_predict(train_lonlat, values, target_lonlat;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw)
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat_target .+ interpolated, 0.0)
end

function _empty_metric_row(; scheme, product, method, fold=missing, group="overall", level="all",
    threshold=NaN, n=0, coverage=0.0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN,
    POD=NaN, FAR=NaN, CSI=NaN)
    return (;
        scheme=String(scheme), product=String(product), method=String(method), fold,
        group=String(group), level=String(level), threshold=Float64(threshold), n=Int(n),
        coverage=Float64(coverage), RMSE=Float64(RMSE), MAE=Float64(MAE),
        Bias=Float64(Bias), r=Float64(r), POD=Float64(POD), FAR=Float64(FAR), CSI=Float64(CSI),
    )
end

function _continuous_row(y_obs, prediction, mask; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs...)
    metric = metric_continuous(y_obs, prediction; mask=mask)
    return _empty_metric_row(;
        kwargs..., n=metric.n, coverage=metric.n / length(mask), RMSE=metric.RMSE,
        MAE=metric.MAE, Bias=metric.Bias, r=metric.r,
    )
end

function _event_row(y_obs, prediction, mask, threshold; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs..., threshold)
    metric = metric_event(y_obs, prediction; mask=mask, thr=threshold)
    return _empty_metric_row(;
        kwargs..., threshold, n, coverage=n / length(mask),
        POD=metric.POD, FAR=metric.FAR, CSI=metric.CSI,
    )
end

function append_stratified_metrics!(
    rows::Vector{NamedTuple}, scheme::String, product::String, method::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, prediction::Matrix{Float64},
    common_mask::BitMatrix, nearest_distance::Vector{Float64}, thresholds::Vector{Float64};
    fold=missing,
)
    base = (; scheme, product, method, fold)
    push!(rows, _continuous_row(y_obs, prediction, common_mask; base..., group="overall", level="all"))

    rain_groups = (
        ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
        ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
    )
    for (name, lower, upper) in rain_groups
        stratum = common_mask .& (y_obs .>= lower) .& (y_obs .< upper)
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="rain_intensity", level=name))
    end

    for year_value in sort(unique(year.(times)))
        time_mask = reshape(year.(times) .== year_value, 1, :)
        stratum = common_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="year", level=string(year_value)))
    end
    for month_value in sort(unique(month.(times)))
        time_mask = reshape(month.(times) .== month_value, 1, :)
        stratum = common_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="month", level=lpad(month_value, 2, '0')))
    end

    distance_groups = (
        ("0_20", -Inf, 20.0), ("20_50", 20.0, 50.0),
        ("50_100", 50.0, 100.0), ("100_plus", 100.0, Inf),
    )
    for (name, lower, upper) in distance_groups
        station_mask = reshape((nearest_distance .>= lower) .& (nearest_distance .< upper), :, 1)
        stratum = common_mask .& station_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="nearest_train_km", level=name))
    end

    for threshold in thresholds
        push!(rows, _event_row(y_obs, prediction, common_mask, threshold;
            base..., group="event_threshold", level=string(threshold)))
    end
    return rows
end

function _common_method_mask(y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}})
    mask = .!isnan.(y_obs)
    for method in BENCHMARK_METHODS
        mask .&= .!isnan.(predictions[method])
    end
    return BitMatrix(mask)
end

function _daily_bootstrap_delta(
    rng::AbstractRNG, times::Vector{DateTime}, y_obs::Matrix{Float64},
    baseline::Matrix{Float64}, gwr::Matrix{Float64}, mask::BitMatrix, reps::Int,
)
    days = Date.(times)
    unique_days = sort(unique(days))
    daily_n = zeros(Int, length(unique_days))
    daily_sse_baseline = zeros(Float64, length(unique_days))
    daily_sse_gwr = zeros(Float64, length(unique_days))
    for (day_index, day) in enumerate(unique_days)
        time_idx = findall(==(day), days)
        day_mask = mask[:, time_idx]
        daily_n[day_index] = count(day_mask)
        daily_sse_baseline[day_index] = sum(abs2, (baseline[:, time_idx] .- y_obs[:, time_idx])[day_mask])
        daily_sse_gwr[day_index] = sum(abs2, (gwr[:, time_idx] .- y_obs[:, time_idx])[day_mask])
    end
    keep = findall(>(0), daily_n)
    isempty(keep) && return Float64[]
    deltas = Vector{Float64}(undef, reps)
    for rep in 1:reps
        sampled = rand(rng, keep, length(keep))
        n = sum(daily_n[sampled])
        rmse_baseline = sqrt(sum(daily_sse_baseline[sampled]) / n)
        rmse_gwr = sqrt(sum(daily_sse_gwr[sampled]) / n)
        deltas[rep] = rmse_baseline - rmse_gwr
    end
    return deltas
end

function _holm_adjust(pvalues::AbstractVector)
    pvalues = Float64.(pvalues)
    m = length(pvalues)
    order = sortperm(pvalues)
    adjusted = fill(NaN, m)
    running = 0.0
    for (rank_index, original_index) in enumerate(order)
        value = min(1.0, (m - rank_index + 1) * pvalues[original_index])
        running = max(running, value)
        adjusted[original_index] = running
    end
    return adjusted
end

function paired_bootstrap_rows(
    cfg::InterpolationBenchmarkConfig, scheme::String, product::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    common_mask::BitMatrix,
)
    rows = NamedTuple[]
    cfg.bootstrap_reps == 0 && return rows
    strata = (
        ("overall", common_mask),
        ("moderate", common_mask .& (y_obs .>= 2.5) .& (y_obs .< 8.0)),
        ("heavy", common_mask .& (y_obs .>= 8.0)),
    )
    for (stratum, mask) in strata
        local_rows = NamedTuple[]
        for (baseline_index, baseline_method) in enumerate(TRADITIONAL_RESIDUAL_METHODS)
            deltas = _daily_bootstrap_delta(
                MersenneTwister(cfg.seed + 1000 * baseline_index + sum(codeunits(product))),
                times, y_obs, predictions[baseline_method], predictions["residual_gwr"],
                BitMatrix(mask), cfg.bootstrap_reps,
            )
            isempty(deltas) && continue
            observed_baseline = metric_continuous(y_obs, predictions[baseline_method]; mask=mask).RMSE
            observed_gwr = metric_continuous(y_obs, predictions["residual_gwr"]; mask=mask).RMSE
            delta = observed_baseline - observed_gwr
            pvalue = min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0)))
            push!(local_rows, (;
                scheme, product, stratum, baseline=baseline_method,
                n=count(mask), n_day=length(unique(Date.(times))), reps=cfg.bootstrap_reps,
                RMSE_baseline=observed_baseline, RMSE_gwr=observed_gwr,
                delta_RMSE=delta, relative_improvement=delta / observed_baseline,
                ci_low=quantile(deltas, 0.025), ci_high=quantile(deltas, 0.975),
                pvalue, pvalue_holm=NaN,
            ))
        end
        adjusted = _holm_adjust([row.pvalue for row in local_rows])
        for (index, row) in enumerate(local_rows)
            push!(rows, merge(row, (; pvalue_holm=adjusted[index])))
        end
    end
    return rows
end

function _write_split(path::String, ids::Vector{String}, folds::Vector{Vector{String}}, scheme::Symbol)
    fold_map = Dict(id => fold for (fold, fold_ids) in enumerate(folds) for id in fold_ids)
    CSV.write(path, DataFrame(
        station_id=ids, fold=[fold_map[id] for id in ids], scheme=fill(string(scheme), length(ids)),
    ))
end

function _one_metric(metrics::DataFrame, product::String, method::String, group::String, level::String)
    rows = filter(row ->
        row.scheme == "balanced_spatial" && row.product == product && row.method == method &&
        ismissing(row.fold) && row.group == group && row.level == level,
        metrics,
    )
    nrow(rows) == 1 || throw(ArgumentError(
        "expected one pooled metric for $product/$method/$group/$level, got $(nrow(rows))",
    ))
    return rows[1, :]
end

function assess_gwr_claim(metrics::DataFrame, bootstrap::DataFrame, products::Vector{String})
    rows = NamedTuple[]
    for product in products
        overall_traditional = [_one_metric(metrics, product, method, "overall", "all")
            for method in TRADITIONAL_RESIDUAL_METHODS]
        best_index = argmin([row.RMSE for row in overall_traditional])
        best_method = TRADITIONAL_RESIDUAL_METHODS[best_index]
        best_overall = overall_traditional[best_index]
        gwr_overall = _one_metric(metrics, product, "residual_gwr", "overall", "all")

        bootstrap_match = filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.stratum == "overall" && row.baseline == best_method,
            bootstrap,
        )
        significant = nrow(bootstrap_match) == 1 &&
            bootstrap_match.ci_low[1] > 0 && bootstrap_match.pvalue_holm[1] < 0.05

        best_heavy = minimum(
            _one_metric(metrics, product, method, "rain_intensity", "heavy").RMSE
            for method in TRADITIONAL_RESIDUAL_METHODS
        )
        gwr_heavy = _one_metric(metrics, product, "residual_gwr", "rain_intensity", "heavy").RMSE
        heavy_improvement = isfinite(best_heavy) && best_heavy > 0 ?
            (best_heavy - gwr_heavy) / best_heavy : NaN

        best_moderate = minimum(
            _one_metric(metrics, product, method, "rain_intensity", "moderate").RMSE
            for method in TRADITIONAL_RESIDUAL_METHODS
        )
        gwr_moderate = _one_metric(metrics, product, "residual_gwr", "rain_intensity", "moderate").RMSE
        moderate_degradation = isfinite(best_moderate) && best_moderate > 0 ?
            (gwr_moderate - best_moderate) / best_moderate : NaN

        year_levels = unique(filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.method == "residual_gwr" && ismissing(row.fold) && row.group == "year",
            metrics,
        ).level)
        year_wins = 0
        for year_level in year_levels
            gwr_year = _one_metric(metrics, product, "residual_gwr", "year", year_level).RMSE
            best_year = minimum(
                _one_metric(metrics, product, method, "year", year_level).RMSE
                for method in TRADITIONAL_RESIDUAL_METHODS
            )
            year_wins += gwr_year < best_year
        end
        majority_years = isempty(year_levels) ? false : year_wins > length(year_levels) / 2

        best_event_rows = [_one_metric(metrics, product, method, "event_threshold", "0.1")
            for method in TRADITIONAL_RESIDUAL_METHODS]
        event_best_index = argmax([row.CSI for row in best_event_rows])
        best_event = best_event_rows[event_best_index]
        gwr_event = _one_metric(metrics, product, "residual_gwr", "event_threshold", "0.1")
        event_not_degraded = gwr_event.CSI >= best_event.CSI - 0.02 &&
            gwr_event.FAR <= best_event.FAR + 0.02
        coverage_acceptable = gwr_overall.coverage >= 0.95

        direction_improved = gwr_overall.RMSE < best_overall.RMSE
        product_supported = significant && heavy_improvement >= 0.05 &&
            moderate_degradation <= 0.02 && majority_years && event_not_degraded &&
            coverage_acceptable
        push!(rows, (;
            product, best_traditional=best_method,
            RMSE_best_traditional=best_overall.RMSE, RMSE_residual_gwr=gwr_overall.RMSE,
            overall_relative_improvement=(best_overall.RMSE - gwr_overall.RMSE) / best_overall.RMSE,
            paired_significant=significant, heavy_relative_improvement=heavy_improvement,
            moderate_relative_degradation=moderate_degradation,
            year_win_count=year_wins, year_count=length(year_levels), majority_years,
            event_not_degraded, common_coverage=gwr_overall.coverage, coverage_acceptable,
            direction_improved, product_supported,
        ))
    end
    supported_products = count(row -> row.product_supported, rows)
    improved_products = count(row -> row.direction_improved, rows)
    overall_supported = supported_products >= 2 && improved_products >= 2
    return DataFrame([merge(row, (;
        supported_product_count=supported_products,
        improved_product_count=improved_products,
        overall_claim_supported=overall_supported,
    )) for row in rows])
end

function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]

    for scheme_symbol in cfg.cv_schemes
        scheme = string(scheme_symbol)
        scheme_dir = joinpath(cfg.mger.outdir, scheme)
        mkpath(scheme_dir)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=cfg.k, seed=cfg.seed)
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

                fold_predictions = Dict{String,Matrix{Float64}}("raw" => y_sat_val)
                for mode in ("direct", "residual"), method in ("idw", "adw", "tps", "gwr")
                    output_method = mode == "direct" ? method : "residual_$(method)"
                    try
                        selected = select_interpolation_parameter!(
                            all_scan_rows, cfg, method, mode, scheme_symbol, product, fold,
                            train_lonlat, y_obs_train, y_sat_train,
                        )
                        fold_predictions[output_method] = predict_selected(
                            selected, method, mode, train_lonlat, val_lonlat,
                            y_obs_train, y_sat_train, y_sat_val,
                        )
                        predictions[output_method][val_idx, :] = fold_predictions[output_method]
                        push!(run_status_rows, (;
                            scheme, product, fold, method=output_method, status="success", error="",
                        ))
                    catch e
                        fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
                        push!(run_status_rows, (;
                            scheme, product, fold, method=output_method, status="failed",
                            error=sprint(showerror, e),
                        ))
                    end
                end

                fold_mask = _common_method_mask(Matrix{Float64}(y_obs[val_idx, :]), fold_predictions)
                if any(fold_mask)
                    for method in BENCHMARK_METHODS
                        append_stratified_metrics!(
                            all_metric_rows, scheme, product, method, data.times,
                            Matrix{Float64}(y_obs[val_idx, :]), fold_predictions[method], fold_mask,
                            nearest_train_distance[val_idx], cfg.event_thresholds; fold,
                        )
                    end
                end
            end

            common_mask = _common_method_mask(y_obs, predictions)
            any(common_mask) || error("[$scheme/$product] no common valid OOF samples across all methods")
            product_dir = joinpath(scheme_dir, lowercase(product))
            mkpath(product_dir)
            for method in BENCHMARK_METHODS
                write_wide(
                    joinpath(product_dir, "oof_$(method).csv"), data.times, ids, predictions[method],
                )
                append_stratified_metrics!(
                    all_metric_rows, scheme, product, method, data.times, y_obs,
                    predictions[method], common_mask, nearest_train_distance, cfg.event_thresholds,
                )
            end
            mask_df = DataFrame(time=Dates.format.(data.times, dateformat"yyyy-mm-ddTHH:MM:SS"))
            for (station_index, station_id) in enumerate(ids)
                mask_df[!, Symbol(station_id)] = common_mask[station_index, :]
            end
            CSV.write(joinpath(product_dir, "common_evaluation_mask.csv"), mask_df)
            if scheme_symbol == :balanced_spatial
                append!(all_bootstrap_rows, paired_bootstrap_rows(
                    cfg, scheme, product, data.times, y_obs, predictions, common_mask,
                ))
            end
        end
    end

    metrics = DataFrame(all_metric_rows)
    scans = DataFrame(all_scan_rows)
    bootstrap = DataFrame(all_bootstrap_rows)
    status = DataFrame(run_status_rows)
    claim = if :balanced_spatial in cfg.cv_schemes && cfg.bootstrap_reps > 0
        assess_gwr_claim(metrics, bootstrap, products)
    else
        DataFrame()
    end
    CSV.write(joinpath(cfg.mger.outdir, "metrics_stratified.csv"), metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_folds.csv"), filter(:fold => (x -> !ismissing(x)), metrics))
    CSV.write(joinpath(cfg.mger.outdir, "metrics_pooled.csv"), filter(:fold => ismissing, metrics))
    CSV.write(joinpath(cfg.mger.outdir, "parameter_scan.csv"), scans)
    CSV.write(joinpath(cfg.mger.outdir, "paired_comparisons.csv"), bootstrap)
    CSV.write(joinpath(cfg.mger.outdir, "run_status.csv"), status)
    if ncol(claim) > 0
        CSV.write(joinpath(cfg.mger.outdir, "claim_assessment.csv"), claim)
    end
    scope = DataFrame(
        key=[
            "validation_target", "training_signal", "temporal_holdout", "primary_cv",
            "secondary_cv", "common_evaluation_mask", "tuning_time_limit",
            "supported_claim", "unsupported_claim",
        ],
        value=[
            "held-out stations at matched observation/satellite timestamps",
            "concurrent training-station observations or observation-minus-satellite residuals",
            "false", "balanced two-dimensional spatial 5-fold", "random station 5-fold",
            "true across all nine methods", string(cfg.tuning_max_times),
            "gauge-network-assisted interpolation/correction at stations not used for fitting",
            "temporal forecast or correction without concurrent gauge observations",
        ],
    )
    CSV.write(joinpath(cfg.mger.outdir, "benchmark_scope.csv"), scope)
    return (; metrics, scans, bootstrap, status, claim)
end
