"""
Post-hoc diagnostics for `run_interpolation_benchmark` results.

Everything here works off the artefacts an existing benchmark run already wrote
(`oof_*.csv`, `common_evaluation_mask.csv`, `split_common.csv`, `parameter_scan.csv`, ...),
so a diagnosis costs no benchmark re-run.

The scoring helper `_metrics` is a deliberate independent reimplementation of
`MGERPipeline.metric_continuous` rather than a call into it: reproducing `metrics_pooled.csv`
from these matrices is the correctness check for the whole reading path, and sharing the
scoring code would make that check vacuous. (`MGERPipeline.jl` is also a top-level script,
not a module, so a standalone module cannot import from it anyway.)
"""
module BenchmarkDiagnostics

using CSV, DataFrames, Dates, Statistics

export read_wide_matrix, read_mask_matrix, read_fold_map
export null_baseline_matrices, satellite_rescale_matrix
export null_baseline_table, satellite_offset_table, mse_decomposition_table
export bandwidth_saturation_table, covariate_contribution_table
export RAIN_CLASSES

"""Rain-intensity strata, mirroring `InterpolationBenchmark.append_stratified_metrics!`."""
const RAIN_CLASSES = [
    ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
    ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
]

# ---------------------------------------------------------------------------- reading

"""Read a wide `time × station` benchmark CSV into a `station × time` matrix ordered by `ids`."""
function read_wide_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = Matrix{Float64}(undef, length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        column = df[!, Symbol(id)]
        out[index, :] = [value === missing ? NaN : Float64(value) for value in column]
    end
    return out
end

"""Read `common_evaluation_mask.csv` into a `station × time` `BitMatrix` ordered by `ids`."""
function read_mask_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = falses(length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        out[index, :] = Bool.(df[!, Symbol(id)])
    end
    return out
end

"""Read `split_common.csv` into `station_id => fold`."""
function read_fold_map(path::AbstractString)
    df = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    return Dict(String(row.station_id) => Int(row.fold) for row in eachrow(df))
end

# ---------------------------------------------------------------------------- scoring

function _metrics(y_true::AbstractMatrix, y_pred::AbstractMatrix, mask::AbstractMatrix)
    a = Float64[]
    b = Float64[]
    @inbounds for index in eachindex(mask)
        mask[index] || continue
        observed = y_true[index]
        predicted = y_pred[index]
        (isnan(observed) || isnan(predicted)) && continue
        push!(a, observed)
        push!(b, predicted)
    end
    n = length(a)
    n == 0 && return (; n=0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN, MSE=NaN, variance=NaN)
    e = b .- a
    bias = mean(e)
    mse = mean(e .^ 2)
    return (;
        n, RMSE=sqrt(mse), MAE=mean(abs.(e)), Bias=bias,
        r=n > 1 ? cor(a, b) : NaN, MSE=mse, variance=mse - bias^2,
    )
end

"""`station × time` mask selecting the cells of `mask` whose gauge value falls in `[low, high)`."""
function _rain_mask(y_obs::AbstractMatrix, mask::AbstractMatrix, low::Float64, high::Float64)
    out = falses(size(mask))
    @inbounds for index in eachindex(mask)
        mask[index] || continue
        observed = y_obs[index]
        isnan(observed) && continue
        out[index] = low <= observed < high
    end
    return out
end

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

# ------------------------------------------------------- D3: where the MSE gap lives

"""
Per-method MSE split into bias² + variance within each rain class, plus each class's share
of the total MSE gap against `reference` (the best traditional interpolator).
"""
function mse_decomposition_table(
    y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    mask::AbstractMatrix; scheme::String, product::String, reference::String="adw",
)
    overall_n = count(mask)
    reference_overall = _metrics(y_obs, predictions[reference], mask)
    rows = NamedTuple[]
    for method in sort(collect(keys(predictions)))
        method == reference && continue
        overall = _metrics(y_obs, predictions[method], mask)
        # Total MSE gap, weighted the way the pooled RMSE weights it.
        total_gap = overall.MSE - reference_overall.MSE
        for (name, low, high) in RAIN_CLASSES
            stratum = _rain_mask(y_obs, mask, low, high)
            score = _metrics(y_obs, predictions[method], stratum)
            reference_score = _metrics(y_obs, predictions[reference], stratum)
            score.n == 0 && continue
            weight = score.n / overall_n
            contribution = weight * (score.MSE - reference_score.MSE)
            push!(rows, (;
                scheme, product, method, reference, rain_class=name, n=score.n,
                sample_share=weight, RMSE=score.RMSE, MSE=score.MSE,
                bias_squared=score.Bias^2, variance=score.variance,
                reference_MSE=reference_score.MSE,
                mse_gap_contribution=contribution,
                gap_share=total_gap != 0 ? contribution / total_gap : NaN,
            ))
        end
    end
    return DataFrame(rows)
end

# --------------------------------------------------- D2a: is the bandwidth grid binding?

"""
Summarise a `parameter_scan.csv` / `joint_bandwidths.csv` scan: for every selected
candidate, whether it sits on an endpoint of the grid it was chosen from, and whether the
CV curve was monotone up to that endpoint (which is what says the grid, not the data,
picked the value).
"""
function bandwidth_saturation_table(scan::DataFrame)
    rows = NamedTuple[]
    keys_of_interest = [:scheme, :product, :fold, :mode, :method, :group, :kernel, :adaptive]
    for subgroup in groupby(unique(scan), keys_of_interest)
        candidates = filter(row -> row.status == "success" && isfinite(row.bw), subgroup)
        nrow(candidates) < 2 && continue
        selected = filter(row -> row.selected === true, candidates)
        nrow(selected) == 1 || continue
        # MGWR tunes each group by coordinate descent, so a group spans several sweeps over the
        # same grid. Only the sweep that produced the winner is a comparable candidate set.
        candidates = filter(row -> row.iteration == selected.iteration[1], candidates)
        nrow(candidates) < 2 && continue
        order = sortperm(candidates.bw)
        widths = candidates.bw[order]
        errors = candidates.RMSE[order]
        chosen = selected.bw[1]
        at_min = chosen == first(widths)
        at_max = chosen == last(widths)
        # Monotone toward the chosen endpoint means the search was clipped, not resolved.
        increasing = all(diff(errors) .> 0)
        decreasing = all(diff(errors) .< 0)
        push!(rows, (;
            scheme=subgroup[1, :scheme], product=subgroup[1, :product],
            fold=subgroup[1, :fold], mode=subgroup[1, :mode],
            method=subgroup[1, :method], group=subgroup[1, :group],
            kernel=subgroup[1, :kernel], adaptive=subgroup[1, :adaptive],
            n_candidates=nrow(candidates), bw_min=first(widths), bw_max=last(widths),
            bw_selected=chosen, at_grid_min=at_min, at_grid_max=at_max,
            monotone_increasing=increasing, monotone_decreasing=decreasing,
            clipped_by_grid=(at_min && increasing) || (at_max && decreasing),
            RMSE_selected=selected.RMSE[1],
            RMSE_at_min=first(errors), RMSE_at_max=last(errors),
        ))
    end
    return DataFrame(rows)
end

# ------------------------------------------------- D4: do the covariates buy anything?

"""
Join each fold's selected-covariate count onto that fold's RMSE, so "more covariates" can
be checked against "better fold". Folds that selected nothing run intercept-only and are the
natural control.
"""
function covariate_contribution_table(status::DataFrame, folds::DataFrame; level::String="all")
    metrics = filter(row -> row.group == "overall" && row.level == level, folds)
    joined = innerjoin(
        select(status, [:scheme, :product, :fold, :method, :covariate_variable_count,
            :covariate_variables, :covariate_effective_roles, :prediction_coverage]),
        select(metrics, [:scheme, :product, :fold, :method, :RMSE, :MAE, :Bias, :r, :n]),
        on=[:scheme, :product, :fold, :method],
    )
    baseline = filter(row -> row.method == "adw", metrics)
    rename!(baseline, :RMSE => :RMSE_adw)
    joined = leftjoin(
        joined, select(baseline, [:scheme, :product, :fold, :RMSE_adw]),
        on=[:scheme, :product, :fold],
    )
    joined.rmse_gap_vs_adw = joined.RMSE .- joined.RMSE_adw
    return sort!(joined, [:scheme, :product, :method, :fold])
end

end # module
