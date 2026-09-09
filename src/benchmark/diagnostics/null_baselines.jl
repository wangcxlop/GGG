# ------------------------------------------------------------------- D0: null baselines

"""
Out-of-fold null predictors, all estimated from training stations only.

- `zero`: constant 0 — the "it is almost always dry" predictor.
- `train_clim`: the training stations' overall mean, one scalar per fold. A *per-station*
  climatology is not estimable for a held-out station, so the global mean is the honest
  station-free analogue.
- `hour_field_mean`: each hour's spatial mean over training stations, i.e. "how much is it
  raining somewhere in Hubei right now" with no spatial structure at all. This is the
  baseline any spatial interpolator has to beat to be doing spatial work.
"""
function null_baseline_matrices(
    y_obs::Matrix{Float64}, ids::Vector{String}, fold_map::Dict{String,Int},
)
    n_station, n_time = size(y_obs)
    fold_of = [fold_map[id] for id in ids]
    predictions = Dict(
        "zero" => zeros(Float64, n_station, n_time),
        "train_clim" => fill(NaN, n_station, n_time),
        "hour_field_mean" => fill(NaN, n_station, n_time),
    )
    for fold in sort(unique(fold_of))
        val = findall(==(fold), fold_of)
        train = findall(!=(fold), fold_of)
        (isempty(val) || isempty(train)) && continue

        finite = filter(isfinite, vec(view(y_obs, train, :)))
        predictions["train_clim"][val, :] .= isempty(finite) ? NaN : mean(finite)

        field = predictions["hour_field_mean"]
        @inbounds for time in 1:n_time
            total = 0.0
            count = 0
            for station in train
                value = y_obs[station, time]
                if !isnan(value)
                    total += value
                    count += 1
                end
            end
            hour_mean = count > 0 ? total / count : NaN
            for station in val
                field[station, time] = hour_mean
            end
        end
    end
    return predictions
end

"""
Score every method and every null on the benchmark's own common mask, overall and by rain
class, and report each method's skill against each null.
"""
function null_baseline_table(
    y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    mask::AbstractMatrix; scheme::String, product::String, nulls::Vector{String},
)
    rows = NamedTuple[]
    strata = vcat([("overall", "all", mask)], [
        (("rain_intensity", name, _rain_mask(y_obs, mask, low, high)))
        for (name, low, high) in RAIN_CLASSES
    ])
    for (group, level, stratum_mask) in strata
        scores = Dict(
            method => _metrics(y_obs, prediction, stratum_mask)
            for (method, prediction) in predictions
        )
        total = count(stratum_mask)
        for (method, score) in sort(collect(scores); by=first)
            skill = Dict(
                "skill_vs_$(null)" => (isnan(score.RMSE) || isnan(scores[null].RMSE)) ? NaN :
                    1 - score.RMSE / scores[null].RMSE
                for null in nulls if haskey(scores, null)
            )
            push!(rows, merge((;
                scheme, product, group, level, method, n=score.n,
                sample_share=total > 0 ? score.n / count(mask) : NaN,
                score.RMSE, score.MAE, score.Bias, score.r,
            ), NamedTuple(Symbol(k) => v for (k, v) in sort(collect(skill); by=first))))
        end
    end
    return DataFrame(rows)
end

