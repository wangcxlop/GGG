using Dates

abstract type SpatialBandwidth end

"""A fixed spatial bandwidth expressed in kilometres."""
struct FixedBandwidth <: SpatialBandwidth
    km::Float64

    function FixedBandwidth(km::Real)
        isfinite(km) && km > 0 || throw(ArgumentError("fixed bandwidth must be positive and finite"))
        new(Float64(km))
    end
end

"""An adaptive spatial bandwidth defined by the number of nearest stations."""
struct AdaptiveBandwidth <: SpatialBandwidth
    neighbors::Int

    function AdaptiveBandwidth(neighbors::Integer)
        neighbors > 0 || throw(ArgumentError("adaptive bandwidth must use at least one neighbour"))
        new(Int(neighbors))
    end
end

"""
Configuration for a historical, two-part satellite-precipitation calibration model.

The occurrence component estimates a spatially weighted wet probability for each
satellite-intensity category. The positive-amount component fits a local ridge
regression to `log1p(observed)` using `log1p(satellite)`, periodic time terms, and
optional station covariates. Coordinates define spatial weights; they are not used
as substitutes for satellite predictors.
"""
Base.@kwdef struct HurdleGWRConfig
    wet_threshold::Float64 = 0.1
    intensity_breaks::Vector{Float64} = [0.1, 2.5, 8.0, 16.0]
    kernel::Int = BISQUARE
    bandwidth::SpatialBandwidth = AdaptiveBandwidth(50)

    annual_cycle::Bool = true
    diurnal_cycle::Bool = true
    occurrence_by_month::Bool = true

    amount_ridge::Float64 = 1e-4
    min_effective_stations::Float64 = 8.0
    min_amount_samples::Float64 = 100.0
    min_occurrence_samples::Float64 = 100.0
    occurrence_prior_strength::Float64 = 20.0
    condition_limit::Float64 = 1e10
    use_lognormal_smearing::Bool = true
    max_positive_prediction::Float64 = Inf
end

"""Fitted historical hurdle-GWR satellite calibration model."""
struct HurdleGWRCalibrator
    config::HurdleGWRConfig
    train_lonlat::Matrix{Float64}
    feature_names::Vector{Symbol}
    feature_mean::Vector{Float64}
    feature_scale::Vector{Float64}
    station_xtx::Array{Float64,3}
    station_xty::Matrix{Float64}
    station_yty::Vector{Float64}
    station_amount_count::Vector{Int}
    occurrence_total::Array{Int,3}
    occurrence_wet::Array{Int,3}
    global_occurrence_total::Vector{Int}
    global_occurrence_wet::Vector{Int}
    global_amount_beta::Vector{Float64}
    global_amount_sigma2::Float64
end

"""Predictions and diagnostics returned by [`predict_hurdle_gwr`](@ref)."""
struct HurdleGWRPrediction
    corrected::Matrix{Float64}
    wet_probability::Matrix{Float64}
    positive_mean::Matrix{Float64}
    amount_coefficients::Matrix{Float64}
    amount_sigma2::Vector{Float64}
    effective_stations::Vector{Float64}
    used_global_amount::BitVector
    occurrence_fallback::Matrix{UInt8}
end

function _validate_hurdle_config(config::HurdleGWRConfig)
    isfinite(config.wet_threshold) && config.wet_threshold > 0 ||
        throw(ArgumentError("wet_threshold must be positive and finite"))
    isempty(config.intensity_breaks) && throw(ArgumentError("intensity_breaks must not be empty"))
    all(isfinite, config.intensity_breaks) ||
        throw(ArgumentError("intensity_breaks must be finite"))
    issorted(config.intensity_breaks) ||
        throw(ArgumentError("intensity_breaks must be sorted"))
    length(unique(config.intensity_breaks)) == length(config.intensity_breaks) ||
        throw(ArgumentError("intensity_breaks must be unique"))
    isapprox(first(config.intensity_breaks), config.wet_threshold; atol=0.0, rtol=1e-12) ||
        throw(ArgumentError("the first intensity break must equal wet_threshold"))
    0 <= config.kernel < length(GWR_KERNELS) ||
        throw(ArgumentError("kernel must be one of 0:$(length(GWR_KERNELS) - 1)"))
    config.amount_ridge >= 0 && isfinite(config.amount_ridge) ||
        throw(ArgumentError("amount_ridge must be finite and non-negative"))
    config.min_effective_stations > 0 && isfinite(config.min_effective_stations) ||
        throw(ArgumentError("min_effective_stations must be positive and finite"))
    config.min_amount_samples > 0 && isfinite(config.min_amount_samples) ||
        throw(ArgumentError("min_amount_samples must be positive and finite"))
    config.min_occurrence_samples > 0 && isfinite(config.min_occurrence_samples) ||
        throw(ArgumentError("min_occurrence_samples must be positive and finite"))
    config.occurrence_prior_strength > 0 && isfinite(config.occurrence_prior_strength) ||
        throw(ArgumentError("occurrence_prior_strength must be positive and finite"))
    config.condition_limit > 1 && isfinite(config.condition_limit) ||
        throw(ArgumentError("condition_limit must be finite and greater than one"))
    (isinf(config.max_positive_prediction) ||
        (isfinite(config.max_positive_prediction) &&
         config.max_positive_prediction >= config.wet_threshold)) ||
        throw(ArgumentError("max_positive_prediction must be Inf or at least wet_threshold"))
    return nothing
end

function _validate_lonlat(lonlat::AbstractMatrix{<:Real}, name::AbstractString)
    size(lonlat, 2) == 2 || throw(DimensionMismatch("$name must have two columns: lon and lat"))
    all(isfinite, lonlat) || throw(ArgumentError("$name contains non-finite coordinates"))
    all(x -> -180 <= x <= 180, @view(lonlat[:, 1])) ||
        throw(ArgumentError("$name contains longitude outside [-180, 180]"))
    all(x -> -90 <= x <= 90, @view(lonlat[:, 2])) ||
        throw(ArgumentError("$name contains latitude outside [-90, 90]"))
    return Matrix{Float64}(lonlat)
end

function _materialize_covariates(covariates, nrow_expected::Int, name::AbstractString)
    if covariates === nothing
        return zeros(Float64, nrow_expected, 0)
    end
    size(covariates, 1) == nrow_expected ||
        throw(DimensionMismatch("$name rows must match station rows"))
    values = Matrix{Float64}(covariates)
    all(isfinite, values) || throw(ArgumentError("$name contains non-finite values"))
    return values
end

function _feature_names(config::HurdleGWRConfig, n_covariate::Int)
    names = Symbol[:intercept, :log1p_satellite]
    if config.annual_cycle
        append!(names, [:annual_sin, :annual_cos])
    end
    if config.diurnal_cycle
        append!(names, [:diurnal_sin, :diurnal_cos])
    end
    append!(names, [Symbol("covariate_", i) for i in 1:n_covariate])
    return names
end

function _raw_amount_features!(
    x::Vector{Float64}, satellite::Float64, time::TimeType,
    covariates::AbstractVector{<:Real}, config::HurdleGWRConfig,
)
    index = 1
    x[index] = 1.0
    index += 1
    x[index] = log1p(satellite)
    index += 1

    if config.annual_cycle
        days_in_year = isleapyear(year(time)) ? 366.0 : 365.0
        day_fraction = time isa Date ? 0.0 :
            (hour(time) + minute(time) / 60 + second(time) / 3600) / 24
        phase = 2pi * ((dayofyear(time) - 1 + day_fraction) / days_in_year)
        x[index] = sin(phase)
        x[index + 1] = cos(phase)
        index += 2
    end

    if config.diurnal_cycle
        decimal_hour = time isa Date ? 0.0 :
            hour(time) + minute(time) / 60 + second(time) / 3600
        phase = 2pi * decimal_hour / 24
        x[index] = sin(phase)
        x[index + 1] = cos(phase)
        index += 2
    end

    @inbounds for value in covariates
        x[index] = Float64(value)
        index += 1
    end
    return x
end

function _standardize_features!(
    x::Vector{Float64}, feature_mean::Vector{Float64}, feature_scale::Vector{Float64},
)
    @inbounds for k in 2:length(x)
        x[k] = (x[k] - feature_mean[k]) / feature_scale[k]
    end
    return x
end

@inline function _valid_precipitation(value::Real)
    return isfinite(value) && value >= 0
end

@inline function _satellite_category(value::Float64, breaks::Vector{Float64})
    return searchsortedlast(breaks, value) + 1
end

function _haversine_km(lon1::Float64, lat1::Float64, lon2::Float64, lat2::Float64)
    dlon = deg2rad(lon2 - lon1)
    dlat = deg2rad(lat2 - lat1)
    a = sin(dlat / 2)^2 + cosd(lat1) * cosd(lat2) * sin(dlon / 2)^2
    return 2 * 6378.388 * asin(sqrt(clamp(a, 0.0, 1.0)))
end

function _spatial_weights(
    train_lonlat::Matrix{Float64}, target_lon::Float64, target_lat::Float64,
    config::HurdleGWRConfig; exclude_train_index::Int=0,
)
    n_station = size(train_lonlat, 1)
    0 <= exclude_train_index <= n_station ||
        throw(ArgumentError("exclude_train_index is outside the training-station range"))

    distances = Vector{Float64}(undef, n_station)
    @inbounds for station in 1:n_station
        distances[station] = station == exclude_train_index ? Inf :
            _haversine_km(
                train_lonlat[station, 1], train_lonlat[station, 2],
                target_lon, target_lat,
            )
    end

    bandwidth = if config.bandwidth isa FixedBandwidth
        config.bandwidth.km
    else
        k = config.bandwidth.neighbors
        finite_distances = filter(isfinite, distances)
        length(finite_distances) >= k ||
            throw(ArgumentError("adaptive bandwidth requests $k neighbours but only $(length(finite_distances)) are available"))
        partialsort!(finite_distances, k)
        finite_distances[k]
    end
    bandwidth = max(bandwidth, sqrt(eps(Float64)))

    kernel_function = GWR_KERNELS[config.kernel + 1]
    weights = zeros(Float64, n_station)
    @inbounds for station in 1:n_station
        distance = distances[station]
        if isfinite(distance)
            weight = kernel_function(distance, bandwidth)
            weights[station] = isfinite(weight) && weight > 0 ? weight : 0.0
        end
    end
    return weights
end

function _aggregate_amount_statistics(
    station_xtx::Array{Float64,3}, station_xty::Matrix{Float64},
    station_yty::Vector{Float64}, station_count::Vector{Int},
    weights::Vector{Float64},
)
    p = size(station_xty, 1)
    A = zeros(Float64, p, p)
    b = zeros(Float64, p)
    yty = 0.0
    weighted_samples = 0.0
    sum_station_weight = 0.0
    sum_station_weight2 = 0.0

    @inbounds for station in eachindex(weights)
        weight = weights[station]
        count = station_count[station]
        if weight > 0 && count > 0
            for col in 1:p, row in 1:p
                A[row, col] += weight * station_xtx[row, col, station]
            end
            for row in 1:p
                b[row] += weight * station_xty[row, station]
            end
            yty += weight * station_yty[station]
            weighted_samples += weight * count
            sum_station_weight += weight
            sum_station_weight2 += weight^2
        end
    end
    effective_stations = sum_station_weight2 > 0 ?
        sum_station_weight^2 / sum_station_weight2 : 0.0
    return A, b, yty, weighted_samples, effective_stations
end

function _solve_amount_model(
    A::Matrix{Float64}, b::Vector{Float64}, yty::Float64,
    weighted_samples::Float64, config::HurdleGWRConfig,
)
    p = length(b)
    A_regularized = copy(A)
    penalty = config.amount_ridge * max(weighted_samples, 1.0)
    @inbounds for k in 2:p
        A_regularized[k, k] += penalty
    end

    condition_number = try
        cond(Symmetric(A_regularized))
    catch
        Inf
    end
    isfinite(condition_number) && condition_number <= config.condition_limit ||
        return nothing

    beta = try
        cholesky(Symmetric(A_regularized)) \ b
    catch
        return nothing
    end
    all(isfinite, beta) || return nothing

    residual_ss = max(yty - 2dot(beta, b) + dot(beta, A * beta), 0.0)
    sigma2 = residual_ss / max(weighted_samples - p, 1.0)
    isfinite(sigma2) || return nothing
    return beta, sigma2, condition_number
end

function _occurrence_counts(
    occurrence_total::Array{Int,3}, occurrence_wet::Array{Int,3},
    weights::Vector{Float64}, months, category::Int,
)
    weighted_total = 0.0
    weighted_wet = 0.0
    sum_station_weight = 0.0
    sum_station_weight2 = 0.0

    @inbounds for station in eachindex(weights)
        station_total = 0
        station_wet = 0
        for month_index in months
            station_total += occurrence_total[station, month_index, category]
            station_wet += occurrence_wet[station, month_index, category]
        end
        weight = weights[station]
        if weight > 0 && station_total > 0
            weighted_total += weight * station_total
            weighted_wet += weight * station_wet
            sum_station_weight += weight
            sum_station_weight2 += weight^2
        end
    end
    effective_stations = sum_station_weight2 > 0 ?
        sum_station_weight^2 / sum_station_weight2 : 0.0
    return weighted_wet, weighted_total, effective_stations
end

function _global_occurrence_rate(model::HurdleGWRCalibrator, category::Int)
    total = model.global_occurrence_total[category]
    wet = model.global_occurrence_wet[category]
    if total == 0
        total = sum(model.global_occurrence_total)
        wet = sum(model.global_occurrence_wet)
    end
    return total > 0 ? clamp(wet / total, 1e-6, 1 - 1e-6) : 0.5
end

function _wet_probability(
    model::HurdleGWRCalibrator, weights::Vector{Float64},
    month_index::Int, category::Int,
)
    config = model.config
    prior_rate = _global_occurrence_rate(model, category)
    prior_wet = config.occurrence_prior_strength * prior_rate
    prior_dry = config.occurrence_prior_strength * (1 - prior_rate)

    if config.occurrence_by_month
        wet, total, effective = _occurrence_counts(
            model.occurrence_total, model.occurrence_wet,
            weights, month_index:month_index, category,
        )
        if total >= config.min_occurrence_samples &&
           effective >= config.min_effective_stations
            probability = (wet + prior_wet) / (total + prior_wet + prior_dry)
            return clamp(probability, 0.0, 1.0), UInt8(0)
        end
    end

    wet, total, effective = _occurrence_counts(
        model.occurrence_total, model.occurrence_wet,
        weights, 1:12, category,
    )
    if total >= config.min_occurrence_samples &&
       effective >= config.min_effective_stations
        probability = (wet + prior_wet) / (total + prior_wet + prior_dry)
        return clamp(probability, 0.0, 1.0), UInt8(1)
    end

    total = model.global_occurrence_total[category]
    wet = model.global_occurrence_wet[category]
    probability = (wet + prior_wet) / (total + prior_wet + prior_dry)
    return clamp(probability, 0.0, 1.0), UInt8(2)
end

"""
    fit_hurdle_gwr(Y_obs, Y_sat, train_lonlat, times; station_covariates=nothing, config=HurdleGWRConfig())

Fit a historical two-part spatially varying satellite calibration model.

Rows of `Y_obs` and `Y_sat` are stations and columns are times. Only historical
training observations are accepted here. Prediction has a separate API that does
not accept target observations, preventing accidental use of concurrent validation
gauges in the standalone-calibration path.
"""
function fit_hurdle_gwr(
    Y_obs::AbstractMatrix{<:Real}, Y_sat::AbstractMatrix{<:Real},
    train_lonlat::AbstractMatrix{<:Real}, times::AbstractVector{<:TimeType};
    station_covariates=nothing, config::HurdleGWRConfig=HurdleGWRConfig(),
)
    _validate_hurdle_config(config)
    size(Y_obs) == size(Y_sat) ||
        throw(DimensionMismatch("Y_obs and Y_sat must have the same size"))
    n_station, n_time = size(Y_obs)
    n_station > 0 && n_time > 0 || throw(ArgumentError("training matrices must not be empty"))
    length(times) == n_time ||
        throw(DimensionMismatch("times length must match the number of matrix columns"))

    lonlat = _validate_lonlat(train_lonlat, "train_lonlat")
    size(lonlat, 1) == n_station ||
        throw(DimensionMismatch("train_lonlat rows must match training stations"))
    covariates = _materialize_covariates(
        station_covariates, n_station, "station_covariates",
    )
    if config.bandwidth isa AdaptiveBandwidth
        config.bandwidth.neighbors <= n_station ||
            throw(ArgumentError("adaptive bandwidth exceeds the number of training stations"))
    end

    feature_names = _feature_names(config, size(covariates, 2))
    p = length(feature_names)
    x = zeros(Float64, p)
    feature_sum = zeros(Float64, p)
    feature_sum2 = zeros(Float64, p)
    amount_count = 0

    @inbounds for station in 1:n_station, time_index in 1:n_time
        observed = Float64(Y_obs[station, time_index])
        satellite = Float64(Y_sat[station, time_index])
        if _valid_precipitation(observed) && observed >= config.wet_threshold &&
           _valid_precipitation(satellite)
            _raw_amount_features!(
                x, satellite, times[time_index], @view(covariates[station, :]), config,
            )
            for k in 2:p
                feature_sum[k] += x[k]
                feature_sum2[k] += x[k]^2
            end
            amount_count += 1
        end
    end
    amount_count > p ||
        throw(ArgumentError("insufficient positive-precipitation samples to fit the amount model"))

    feature_mean = zeros(Float64, p)
    feature_scale = ones(Float64, p)
    @inbounds for k in 2:p
        feature_mean[k] = feature_sum[k] / amount_count
        variance = max(feature_sum2[k] / amount_count - feature_mean[k]^2, 0.0)
        scale = sqrt(variance)
        feature_scale[k] = scale > sqrt(eps(Float64)) ? scale : 1.0
    end

    station_xtx = zeros(Float64, p, p, n_station)
    station_xty = zeros(Float64, p, n_station)
    station_yty = zeros(Float64, n_station)
    station_amount_count = zeros(Int, n_station)

    n_category = length(config.intensity_breaks) + 1
    occurrence_total = zeros(Int, n_station, 12, n_category)
    occurrence_wet = zeros(Int, n_station, 12, n_category)

    @inbounds for station in 1:n_station, time_index in 1:n_time
        observed = Float64(Y_obs[station, time_index])
        satellite = Float64(Y_sat[station, time_index])
        if !_valid_precipitation(observed) || !_valid_precipitation(satellite)
            continue
        end

        month_index = config.occurrence_by_month ? month(times[time_index]) : 1
        category = _satellite_category(satellite, config.intensity_breaks)
        occurrence_total[station, month_index, category] += 1
        occurrence_wet[station, month_index, category] += observed >= config.wet_threshold

        if observed >= config.wet_threshold
            _raw_amount_features!(
                x, satellite, times[time_index], @view(covariates[station, :]), config,
            )
            _standardize_features!(x, feature_mean, feature_scale)
            response = log1p(observed)
            for col in 1:p, row in 1:p
                station_xtx[row, col, station] += x[row] * x[col]
            end
            for row in 1:p
                station_xty[row, station] += x[row] * response
            end
            station_yty[station] += response^2
            station_amount_count[station] += 1
        end
    end

    global_occurrence_total = vec(sum(occurrence_total, dims=(1, 2)))
    global_occurrence_wet = vec(sum(occurrence_wet, dims=(1, 2)))
    sum(global_occurrence_total) > 0 ||
        throw(ArgumentError("no valid observation-satellite pairs are available"))

    global_weights = ones(Float64, n_station)
    A, b, yty, weighted_samples, _ = _aggregate_amount_statistics(
        station_xtx, station_xty, station_yty, station_amount_count, global_weights,
    )
    global_fit = _solve_amount_model(A, b, yty, weighted_samples, config)
    global_fit === nothing &&
        throw(ArgumentError("global positive-amount calibration is numerically singular"))
    global_beta, global_sigma2, _ = global_fit

    return HurdleGWRCalibrator(
        config, lonlat, feature_names, feature_mean, feature_scale,
        station_xtx, station_xty, station_yty, station_amount_count,
        occurrence_total, occurrence_wet,
        global_occurrence_total, global_occurrence_wet,
        global_beta, global_sigma2,
    )
end

function _local_amount_parameters(
    model::HurdleGWRCalibrator, weights::Vector{Float64},
)
    A, b, yty, weighted_samples, effective_stations = _aggregate_amount_statistics(
        model.station_xtx, model.station_xty, model.station_yty,
        model.station_amount_count, weights,
    )
    local_fit = if weighted_samples >= model.config.min_amount_samples &&
                   effective_stations >= model.config.min_effective_stations
        _solve_amount_model(A, b, yty, weighted_samples, model.config)
    else
        nothing
    end
    if local_fit === nothing
        return model.global_amount_beta, model.global_amount_sigma2,
            effective_stations, true
    end
    beta, sigma2, _ = local_fit
    return beta, sigma2, effective_stations, false
end

"""
    predict_hurdle_gwr(model, Y_sat_target, target_lonlat, times; target_covariates=nothing)

Predict standalone satellite calibration at target locations. This function never
accepts target rain-gauge observations. `exclude_train_indices[j]` may optionally
identify a collocated training station to exclude when producing leave-one-station-
out predictions; use zero for no exclusion.
"""
function predict_hurdle_gwr(
    model::HurdleGWRCalibrator,
    Y_sat_target::AbstractMatrix{<:Real}, target_lonlat::AbstractMatrix{<:Real},
    times::AbstractVector{<:TimeType}; target_covariates=nothing,
    exclude_train_indices=nothing,
)
    n_target, n_time = size(Y_sat_target)
    length(times) == n_time ||
        throw(DimensionMismatch("times length must match target matrix columns"))
    lonlat = _validate_lonlat(target_lonlat, "target_lonlat")
    size(lonlat, 1) == n_target ||
        throw(DimensionMismatch("target_lonlat rows must match target stations"))
    n_covariate = length(model.feature_names) - 2 -
        2Int(model.config.annual_cycle) - 2Int(model.config.diurnal_cycle)
    covariates = _materialize_covariates(
        target_covariates, n_target, "target_covariates",
    )
    size(covariates, 2) == n_covariate ||
        throw(DimensionMismatch("target covariate columns must match fitted covariates"))

    exclusions = if exclude_train_indices === nothing
        zeros(Int, n_target)
    else
        length(exclude_train_indices) == n_target ||
            throw(DimensionMismatch("exclude_train_indices length must match target stations"))
        Int.(exclude_train_indices)
    end

    corrected = fill(NaN, n_target, n_time)
    wet_probability = fill(NaN, n_target, n_time)
    positive_mean = fill(NaN, n_target, n_time)
    amount_coefficients = fill(NaN, n_target, length(model.feature_names))
    amount_sigma2 = fill(NaN, n_target)
    effective_stations = zeros(Float64, n_target)
    used_global_amount = falses(n_target)
    occurrence_fallback = fill(UInt8(3), n_target, n_time)

    x = zeros(Float64, length(model.feature_names))
    @inbounds for target in 1:n_target
        weights = _spatial_weights(
            model.train_lonlat, lonlat[target, 1], lonlat[target, 2], model.config;
            exclude_train_index=exclusions[target],
        )
        beta, sigma2, effective, used_global = _local_amount_parameters(model, weights)
        amount_coefficients[target, :] .= beta
        amount_sigma2[target] = sigma2
        effective_stations[target] = effective
        used_global_amount[target] = used_global

        for time_index in 1:n_time
            satellite = Float64(Y_sat_target[target, time_index])
            _valid_precipitation(satellite) || continue

            category = _satellite_category(satellite, model.config.intensity_breaks)
            month_index = model.config.occurrence_by_month ? month(times[time_index]) : 1
            probability, fallback = _wet_probability(
                model, weights, month_index, category,
            )

            _raw_amount_features!(
                x, satellite, times[time_index], @view(covariates[target, :]), model.config,
            )
            _standardize_features!(x, model.feature_mean, model.feature_scale)
            mean_log = dot(x, beta)
            if model.config.use_lognormal_smearing
                mean_log += 0.5 * sigma2
            end
            finite_log_limit = log(floatmax(Float64)) - 2
            conditional_mean = exp(clamp(mean_log, -finite_log_limit, finite_log_limit)) - 1
            conditional_mean = max(conditional_mean, model.config.wet_threshold)
            if isfinite(model.config.max_positive_prediction)
                conditional_mean = min(conditional_mean, model.config.max_positive_prediction)
            end

            wet_probability[target, time_index] = probability
            positive_mean[target, time_index] = conditional_mean
            corrected[target, time_index] = probability * conditional_mean
            occurrence_fallback[target, time_index] = fallback
        end
    end

    return HurdleGWRPrediction(
        corrected, wet_probability, positive_mean, amount_coefficients,
        amount_sigma2, effective_stations, used_global_amount, occurrence_fallback,
    )
end

export SpatialBandwidth, FixedBandwidth, AdaptiveBandwidth
export HurdleGWRConfig, HurdleGWRCalibrator, HurdleGWRPrediction
export fit_hurdle_gwr, predict_hurdle_gwr
