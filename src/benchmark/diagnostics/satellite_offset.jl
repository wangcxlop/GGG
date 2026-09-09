# --------------------------------------------------------------- D1: satellite offset

"""
Out-of-fold `a + b·satellite`, with `a, b` fit by OLS on training stations only.

`residual_gwr` forces `b = 1`; this measures what the same satellite field is worth when
that single coefficient is allowed to be estimated instead.
"""
function satellite_rescale_matrix(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, ids::Vector{String},
    fold_map::Dict{String,Int}; clip::Bool=true,
)
    n_station, n_time = size(y_obs)
    prediction = fill(NaN, n_station, n_time)
    fold_of = [fold_map[id] for id in ids]
    coefficients = NamedTuple[]
    for fold in sort(unique(fold_of))
        val = findall(==(fold), fold_of)
        train = findall(!=(fold), fold_of)
        (isempty(val) || isempty(train)) && continue

        x = Float64[]
        y = Float64[]
        @inbounds for station in train, time in 1:n_time
            observed = y_obs[station, time]
            satellite = y_sat[station, time]
            (isnan(observed) || isnan(satellite)) && continue
            push!(x, satellite)
            push!(y, observed)
        end
        mean_x = mean(x)
        mean_y = mean(y)
        variance_x = mean(x .^ 2) - mean_x^2
        slope = variance_x > 0 ? (mean(x .* y) - mean_x * mean_y) / variance_x : 0.0
        intercept = mean_y - slope * mean_x
        push!(coefficients, (; fold, intercept, slope, n_train=length(x)))

        @inbounds for station in val, time in 1:n_time
            satellite = y_sat[station, time]
            isnan(satellite) && continue
            value = intercept + slope * satellite
            prediction[station, time] = clip ? max(value, 0.0) : value
        end
    end
    return prediction, DataFrame(coefficients)
end

"""
Exact decomposition of what the residual GWR does to the satellite field.

With `ĝ = prediction − satellite` on cells the `max(·, 0)` clip did not touch,

    MSE(satellite + ĝ) − MSE(satellite) = E[ĝ²] − 2·E[ĝ·(obs − sat)]

The first term is the variance the correction injects; the second is the error it actually
removes. If `E[ĝ·(obs − sat)] ≈ 0` while `E[ĝ²] > 0`, the correction is provably adding
noise, and no bandwidth or covariate change can rescue the formulation.
"""
function satellite_offset_table(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64},
    predictions::Dict{String,Matrix{Float64}}, mask::AbstractMatrix;
    scheme::String, product::String, methods::Vector{String},
)
    rows = NamedTuple[]
    for method in methods
        haskey(predictions, method) || continue
        prediction = predictions[method]
        correction = Float64[]
        truth = Float64[]
        clipped = 0
        total = 0
        @inbounds for index in eachindex(mask)
            mask[index] || continue
            observed = y_obs[index]
            satellite = y_sat[index]
            predicted = prediction[index]
            (isnan(observed) || isnan(satellite) || isnan(predicted)) && continue
            total += 1
            # The clip is only active where it drove the prediction to exactly zero.
            if predicted <= 0 && satellite > 0
                clipped += 1
                continue
            end
            push!(correction, predicted - satellite)
            push!(truth, observed - satellite)
        end
        n = length(correction)
        n == 0 && continue
        injected = mean(correction .^ 2)
        removed = 2 * mean(correction .* truth)
        # Rescaling the whole correction by `s` gives MSE(s) = s^2 E[g^2] - 2s E[g*r],
        # minimised at s = E[g*r] / E[g^2]. s well below 1 means the correction has the
        # right idea but the wrong magnitude - GWR is unbiased and never shrinks.
        scale = injected > 0 ? (removed / 2) / injected : NaN
        push!(rows, (;
            scheme, product, method, n, clipped_fraction=total > 0 ? clipped / total : NaN,
            correction_sd=sqrt(max(injected - mean(correction)^2, 0.0)),
            residual_sd=sqrt(max(mean(truth .^ 2) - mean(truth)^2, 0.0)),
            corr_correction_residual=n > 1 ? cor(correction, truth) : NaN,
            variance_injected=injected, error_removed=removed,
            delta_mse_vs_raw=injected - removed,
            optimal_correction_scale=scale,
            delta_mse_if_shrunk=isnan(scale) ? NaN : -scale^2 * injected,
        ))
    end
    return DataFrame(rows)
end

