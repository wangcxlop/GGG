"""
Post-hoc diagnostics for `run_interpolation_benchmark` results.

Everything here works off the artefacts an existing benchmark run already wrote
(`oof_*.csv`, `common_evaluation_mask.csv`, `split_common.csv`, `parameter_scan.csv`, ...),
so a diagnosis costs no benchmark re-run.

The scoring helper `_metrics` is a deliberate independent reimplementation of
`MixedGWR.metric_continuous` rather than a call into it: reproducing `metrics_pooled.csv`
from these matrices is the correctness check for the whole reading path, and sharing the
scoring code would make that check vacuous.

That is now the only reason. `metric_continuous` used to live in `MGERPipeline.jl`, a top-level
script a module could not import from at all; it has since moved into `MixedGWR`, so the
duplication here is a choice rather than a constraint. The two also differ at the edges -
`_metrics` returns a NaN row for an empty sample where `metric_continuous` asserts, and adds
`MSE`/`variance` - so they are not interchangeable as written.
"""
module BenchmarkDiagnostics

using CSV, DataFrames, Dates, Statistics
using MixedGWR
using Main: TraditionalInterpolation

export read_wide_matrix, read_mask_matrix, read_fold_map, load_prediction_matrices
export read_run_grid, load_gauge_matrix
export null_baseline_matrices, satellite_rescale_matrix
export null_baseline_table, satellite_offset_table, mse_decomposition_table
export satellite_quadrant_table
export station_weight_matrix, local_satellite_slope, local_anchor_prediction
export false_alarm_coherence_table
export bandwidth_saturation_table, covariate_contribution_table
export rebuild_common_mask, dropout_table, mask_cost_table, run_comparison_table
export JOINT_MASK_METHODS
export benchmark_diagnostics_outdir
export RAIN_CLASSES, MASK_METHODS

"""Rain-intensity strata, mirroring `InterpolationBenchmark.append_stratified_metrics!`."""
const RAIN_CLASSES = [
    ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
    ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
]

"""
Methods whose finiteness defines the shared evaluation mask, mirroring
`InterpolationBenchmark.MASK_METHODS`. Duplicated rather than imported because
`InterpolationBenchmark.jl` is a top-level script, not a module; `rebuild_common_mask` is
checked against the run's own `common_evaluation_mask.csv` so the copy cannot drift silently.
"""
const MASK_METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]

# ---------------------------------------------------------------------------- reading

"""
Where a diagnosis of `run_dir` is written: `output/benchmark_diagnostics/<run name>`.

Each diagnostic run gets its own subdirectory named for the benchmark run it read, so diagnosing a
second run never overwrites the first one's numbers. Repeated verbatim in
`run_benchmark_diagnostics.jl`, `run_mgwr_diagnostics.jl` and `run_claim_reassessment.jl`, which
also each checked `run_dir` exists first - done here so the error is worded once.
"""
function benchmark_diagnostics_outdir(root::AbstractString, run_dir::AbstractString)
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(root, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)
    return outdir
end


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

"""Read every `oof_<method>.csv` in one `(scheme, product)` directory into `method => matrix`."""
function load_prediction_matrices(product_dir::AbstractString, ids::Vector{String})
    predictions = Dict{String,Matrix{Float64}}()
    for path in readdir(product_dir; join=true)
        name = basename(path)
        startswith(name, "oof_") && endswith(name, ".csv") || continue
        predictions[name[5:end-4]] = read_wide_matrix(path, ids)
    end
    return predictions
end

"""Read `split_common.csv` into `station_id => fold`."""
function read_fold_map(path::AbstractString)
    df = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    return Dict(String(row.station_id) => Int(row.fold) for row in eachrow(df))
end

_as_datetime(value::DateTime) = value
# Gauge files stamp a trailing `Z` where the benchmark's own artefacts do not; stripped the same
# way `MGERDataPrep.jl:9` and `MGERPipeline.jl:62` already do.
_as_datetime(value::AbstractString) = DateTime(replace(strip(String(value)), "Z" => ""))

"""
Station ids and timestamps of a wide `time × station` benchmark CSV, in file order.

The run's own artefacts define the grid a diagnosis has to score on, so it is read off one of
them rather than rebuilt from a config - which also means a diagnosis needs none of the inputs
the run itself consumed.
"""
function read_run_grid(path::AbstractString)
    df = CSV.read(path, DataFrame)
    ids = String.(names(df)[2:end])
    times = [_as_datetime(value) for value in df[!, 1]]
    return ids, times
end

"""
Gauge observations as a `station × time` matrix on a run's own grid, plus the number of the
run's hours the file does not cover.

Callers are expected to treat a non-zero shortfall as an error: silently scoring against a
partly-absent truth is worse than not scoring at all.
"""
function load_gauge_matrix(path::AbstractString, ids::Vector{String}, times::Vector{DateTime})
    df = CSV.read(path, DataFrame)
    row_of = Dict(_as_datetime(value) => row for (row, value) in enumerate(df[!, 1]))
    columns = [df[!, Symbol(id)] for id in ids]
    out = fill(NaN, length(ids), length(times))
    unmatched = 0
    for (t, time) in enumerate(times)
        row = get(row_of, time, 0)
        if row == 0
            unmatched += 1
            continue
        end
        for s in eachindex(columns)
            value = columns[s][row]
            out[s, t] = value === missing ? NaN : Float64(value)
        end
    end
    return out, unmatched
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

"""
The same gap split, but by the *joint* wet/dry state of gauge and satellite rather than by
rain class alone.

`mse_decomposition_table` says the GWR family's loss lives in the dry gauge cells; it cannot say
*which* dry cells, and the two halves behave in opposite directions. The residual framing
predicts `y_sat + correction`, so a satellite false alarm (`obs` dry, `sat` wet) is an error the
correction has to cancel out of a value it was handed, while a cell both agree is dry costs the
correction nothing. Splitting on both is what separates "the estimator is worse" from "the
anchor is wrong", and those have different fixes.

Quadrants are named `<obs>_<sat>` over `{dry, wet}` at `threshold`, matching the benchmark's
`rain_threshold`. `mean_sat` is reported because the size of the anchor is the size of the
problem in the false-alarm quadrant.
"""
function satellite_quadrant_table(
    y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    mask::AbstractMatrix; scheme::String, product::String, reference::String="adw",
    satellite::String="raw", threshold::Float64=0.1,
)
    overall_n = count(mask)
    y_sat = predictions[satellite]
    reference_overall = _metrics(y_obs, predictions[reference], mask)
    quadrants = [
        ("dry_dry", false, false), ("dry_wet", false, true),
        ("wet_dry", true, false), ("wet_wet", true, true),
    ]
    rows = NamedTuple[]
    for method in sort(collect(keys(predictions)))
        method == reference && continue
        total_gap = _metrics(y_obs, predictions[method], mask).MSE - reference_overall.MSE
        for (name, obs_wet, sat_wet) in quadrants
            cell = falses(size(mask))
            satellite_sum = 0.0
            @inbounds for index in eachindex(mask)
                mask[index] || continue
                observed = y_obs[index]
                anchor = y_sat[index]
                (isnan(observed) || isnan(anchor)) && continue
                (observed >= threshold) == obs_wet && (anchor >= threshold) == sat_wet || continue
                cell[index] = true
                satellite_sum += anchor
            end
            score = _metrics(y_obs, predictions[method], cell)
            score.n == 0 && continue
            reference_score = _metrics(y_obs, predictions[reference], cell)
            weight = count(cell) / overall_n
            contribution = weight * (score.MSE - reference_score.MSE)
            push!(rows, (;
                scheme, product, method, reference, quadrant=name, n=score.n,
                sample_share=weight, mean_satellite=satellite_sum / count(cell),
                RMSE=score.RMSE, MSE=score.MSE,
                bias_squared=score.Bias^2, variance=score.variance,
                reference_MSE=reference_score.MSE,
                satellite_MSE=_metrics(y_obs, y_sat, cell).MSE,
                mse_gap_contribution=contribution,
                gap_share=total_gap != 0 ? contribution / total_gap : NaN,
            ))
        end
    end
    return DataFrame(rows)
end

# ------------------------------------------- D5: which cells does a method fail to predict?

"""
Rebuild `common_evaluation_mask.csv` from the out-of-fold matrices: a cell counts only when
the gauge value and every mask-defining method's prediction are all finite.

Mirrors `InterpolationBenchmark._common_method_mask`. Pass a shorter `methods` list to ask
what the mask would have been had a given method not been allowed to define it.
"""
function rebuild_common_mask(
    y_obs::Matrix{Float64}, predictions::AbstractDict{String,Matrix{Float64}};
    methods::Vector{String}=MASK_METHODS,
)
    mask = .!isnan.(y_obs)
    for method in methods
        haskey(predictions, method) ||
            throw(ArgumentError("mask-defining method is absent from the run: $method"))
        mask .&= .!isnan.(predictions[method])
    end
    return BitMatrix(mask)
end

"""
Per fold, how many cells a method failed to predict, and whether the failures come in whole
hours or scattered station-by-station.

The denominator is the *evaluable* cell count — gauge and satellite both finite — not the
shared mask, because the shared mask has already had these very failures removed from it.
A whole-hour dropout (every evaluable validation station NaN in the same hour) is the
signature of a fit-level failure such as back-fitting non-convergence, which discards the
hour for all targets at once; scattered dropouts instead mean individual local fits were
underdetermined. The split is what distinguishes the two, and they need different fixes.
"""
function dropout_table(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64},
    predictions::AbstractDict{String,Matrix{Float64}}, ids::Vector{String},
    fold_map::Dict{String,Int}; scheme::String, product::String,
    methods::Vector{String}=["mgwr", "residual_gwr", "mixed_gwr"],
)
    fold_of = [fold_map[id] for id in ids]
    n_time = size(y_obs, 2)
    rows = NamedTuple[]
    for method in methods
        haskey(predictions, method) || continue
        prediction = predictions[method]
        for fold in sort(unique(fold_of))
            val = findall(==(fold), fold_of)
            isempty(val) && continue
            evaluable = 0
            dropped = 0
            hours_evaluable = 0
            hours_all_dropped = 0
            cells_in_dropped_hours = 0
            @inbounds for time in 1:n_time
                hour_evaluable = 0
                hour_dropped = 0
                for station in val
                    (isnan(y_obs[station, time]) || isnan(y_sat[station, time])) && continue
                    hour_evaluable += 1
                    isnan(prediction[station, time]) && (hour_dropped += 1)
                end
                hour_evaluable == 0 && continue
                hours_evaluable += 1
                evaluable += hour_evaluable
                dropped += hour_dropped
                if hour_dropped == hour_evaluable
                    hours_all_dropped += 1
                    cells_in_dropped_hours += hour_dropped
                end
            end
            push!(rows, (;
                scheme, product, fold, method,
                n_evaluable=evaluable, n_dropped=dropped,
                own_coverage=evaluable > 0 ? 1 - dropped / evaluable : NaN,
                hours_evaluable, hours_all_dropped,
                hour_dropout_fraction=hours_evaluable > 0 ?
                    hours_all_dropped / hours_evaluable : NaN,
                cells_in_all_dropped_hours=cells_in_dropped_hours,
                cells_scattered=dropped - cells_in_dropped_hours,
            ))
        end
    end
    return DataFrame(rows)
end

"""
Methods whose mask failures are correlated, and must therefore be excluded together to be seen.

The joint-covariate models share `build_joint_fold_context`'s predictor matrices and
`dynamic_covariate_predict`'s `valid` guard, so a station-hour with a missing covariate is
dropped by all three at once. See [`mask_cost_table`](@ref) for why that defeats a one-at-a-time
exclusion.
"""
const JOINT_MASK_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]

"""
What the mask-defining methods in `excluded` cost every other method.

The shared mask keeps only cells where *every* method in `MASK_METHODS` is finite, so an
unstable method silently shrinks the denominator the whole field is scored on. This rebuilds
the mask with `excluded` removed from the defining set and rescores everyone on both, so the
cost is a number rather than an argument. Cells recovered are by construction the harder ones
(the excluded methods could not fit them), so a method's RMSE is expected to rise on the wider
mask — the comparison of interest is between methods, not against zero.

`excluded` takes a **set**, not just one name, and that is the point. Excluding one method at a
time is blind to correlated failure, which on this data is the only failure there is: measured on
the full nested run (`balanced_spatial`), every one of the 3,963 / 3,986 / 3,986 cells the mask
drops for FY4B / GPM / GSMaP is dropped by `residual_gwr`, `mixed_gwr` **and** `mgwr` together —
they share the same predictor matrices and the same `valid` guard, so a station-hour with a
missing covariate fails all three. `raw`, `idw`, `adw`, `tps` and direct `gwr` cost nothing at all.
Excluding `mgwr` alone therefore recovers 48 cells on FY4B and *zero* on GPM and GSMaP, because
the other two still mask them, and that zero reads as "the shared mask is harmless" when it only
means "one at a time cannot see this". Pass `JOINT_MASK_METHODS` to get the honest number.
"""
function mask_cost_table(
    y_obs::Matrix{Float64}, predictions::AbstractDict{String,Matrix{Float64}};
    scheme::String, product::String,
    excluded::Union{AbstractString,AbstractVector{<:AbstractString}}="mgwr",
    methods::Vector{String}=MASK_METHODS,
)
    dropped = excluded isa AbstractString ? [String(excluded)] : String.(excluded)
    isempty(dropped) &&
        throw(ArgumentError("excluded must name at least one mask-defining method"))
    allunique(dropped) || throw(ArgumentError("excluded names a method twice: $(dropped)"))
    for name in dropped
        name in methods || throw(ArgumentError(
            "$name does not define the mask, so removing it changes nothing",
        ))
    end
    length(dropped) < length(methods) || throw(ArgumentError(
        "excluding every mask-defining method leaves no mask to compare against",
    ))
    full = rebuild_common_mask(y_obs, predictions; methods)
    reduced = rebuild_common_mask(y_obs, predictions; methods=filter(!in(dropped), methods))
    # One row per method carries the same `excluded_method`; joined rather than one column per
    # name so the CSV keeps its shape, and a single-name call still writes exactly that name.
    excluded_label = join(dropped, ",")
    rows = NamedTuple[]
    for method in sort(collect(keys(predictions)))
        on_full = _metrics(y_obs, predictions[method], full)
        on_reduced = _metrics(y_obs, predictions[method], reduced)
        push!(rows, (;
            scheme, product, excluded_method=excluded_label, method,
            mask_cells_full=count(full), mask_cells_reduced=count(reduced),
            cells_recovered=count(reduced) - count(full),
            n_full=on_full.n, RMSE_full=on_full.RMSE,
            n_reduced=on_reduced.n, RMSE_reduced=on_reduced.RMSE,
            RMSE_delta=on_reduced.RMSE - on_full.RMSE,
            relative_change=(isnan(on_full.RMSE) || on_full.RMSE == 0) ? NaN :
                (on_reduced.RMSE - on_full.RMSE) / on_full.RMSE,
        ))
    end
    return DataFrame(rows)
end

# ------------------------------------------------------- D6: comparing two benchmark runs

"""
Score two runs of the benchmark against each other on cells they both evaluated.

A raw RMSE delta between two run directories is meaningless here. The shared evaluation mask
keeps only cells where every mask-defining method is finite, so any change that alters one
method's coverage resizes the denominator for all of them — and the cells that appear or vanish
are systematically the hard ones. Two masks are therefore intersected before anything is scored,
and each method is additionally restricted to cells where *both* runs produced a finite value,
which makes `delta_paired` a genuine paired difference rather than two numbers from two samples.

`RMSE_before` / `RMSE_after` are each run's own published-style number (its own mask, its own
NaNs dropped) and are reported alongside so the size of the mask effect is visible rather than
hidden. Trust `delta_paired`; read the own-mask columns only to see how much the mask moved.
"""
function run_comparison_table(
    y_obs::Matrix{Float64},
    before::AbstractDict{String,Matrix{Float64}}, before_mask::AbstractMatrix,
    after::AbstractDict{String,Matrix{Float64}}, after_mask::AbstractMatrix;
    scheme::String, product::String,
)
    shared = BitMatrix(before_mask .& after_mask)
    rows = NamedTuple[]
    for method in sort(collect(union(keys(before), keys(after))))
        in_before = haskey(before, method)
        in_after = haskey(after, method)
        own_before = in_before ? _metrics(y_obs, before[method], before_mask) : nothing
        own_after = in_after ? _metrics(y_obs, after[method], after_mask) : nothing
        paired_before = (; n=0, RMSE=NaN)
        paired_after = (; n=0, RMSE=NaN)
        if in_before && in_after
            pair = BitMatrix(shared .& .!isnan.(before[method]) .& .!isnan.(after[method]))
            paired_before = _metrics(y_obs, before[method], pair)
            paired_after = _metrics(y_obs, after[method], pair)
        end
        delta = paired_after.RMSE - paired_before.RMSE
        push!(rows, (;
            scheme, product, method,
            present_in=in_before && in_after ? "both" : (in_before ? "before" : "after"),
            mask_cells_before=count(before_mask), mask_cells_after=count(after_mask),
            mask_cells_shared=count(shared),
            n_before=own_before === nothing ? 0 : own_before.n,
            n_after=own_after === nothing ? 0 : own_after.n,
            RMSE_before=own_before === nothing ? NaN : own_before.RMSE,
            RMSE_after=own_after === nothing ? NaN : own_after.RMSE,
            n_paired=paired_before.n,
            RMSE_paired_before=paired_before.RMSE,
            RMSE_paired_after=paired_after.RMSE,
            delta_paired=delta,
            relative_paired=isnan(delta) || paired_before.RMSE == 0 ? NaN :
                delta / paired_before.RMSE,
        ))
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
        # `!isnan` rather than `isfinite`: the filter's job is excluding the `NaN` bw of the
        # idw/adw/tps rows, and the GWR family's explicit global candidate is `bw = Inf`, which
        # sorts last and so still reads as "chose the widest candidate offered".
        candidates = filter(row -> row.status == "success" && !isnan(row.bw), subgroup)
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
        # Landing on the widest candidate only counts as clipping if a wider one could exist.
        # The GWR family's top candidate is the explicit global fit (`bw = Inf`), which is the
        # end of the bandwidth continuum, not the end of an arbitrary grid - selecting it is a
        # resolved answer ("this coefficient wants to be global"), so it must not be reported
        # as a grid that needs extending.
        extensible_max = isfinite(last(widths))
        push!(rows, (;
            scheme=subgroup[1, :scheme], product=subgroup[1, :product],
            fold=subgroup[1, :fold], mode=subgroup[1, :mode],
            method=subgroup[1, :method], group=subgroup[1, :group],
            kernel=subgroup[1, :kernel], adaptive=subgroup[1, :adaptive],
            n_candidates=nrow(candidates), bw_min=first(widths), bw_max=last(widths),
            bw_selected=chosen, at_grid_min=at_min, at_grid_max=at_max,
            monotone_increasing=increasing, monotone_decreasing=decreasing,
            clipped_by_grid=(at_min && increasing) || (at_max && decreasing && extensible_max),
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

# ------------------------- D8: can a *locally* varying satellite coefficient close the gap?

"""
Mean-Earth radius, the metric `JointCovariateModels._design_at` measures its bandwidths on.

`TraditionalInterpolation.haversine_distance_matrix` uses the WGS72 equatorial radius instead (see
the comment at the top of that file: the two families are deliberately not on the same metric). A
bandwidth read out of `joint_bandwidths.csv` belongs to the joint path, so distances are rescaled
onto its radius rather than the trigonometry being written a third time - the haversine formula is
linear in the radius, so the rescaling is exact.
"""
const JOINT_EARTH_RADIUS_KM = 6371.0088

"""
Geographic weights between a run's stations, on the joint-covariate path's distance metric.

The diagonal is set to `Inf` before weighting rather than zeroed afterwards. Every kernel in
`MixedGWR.GWR_KERNELS` returns 0 at an infinite distance, and under `adaptive` the self-distance
would otherwise occupy one of the `bw` nearest slots and shift the bandwidth inward. Each row
therefore describes a leave-one-out fit, which is what the benchmark does: a held-out station is
never in its own training set.
"""
function station_weight_matrix(
    lonlat::Matrix{Float64}; kernel::Int, bw::Float64, adaptive::Bool,
)
    n = size(lonlat, 1)
    size(lonlat, 2) == 2 || throw(DimensionMismatch("lonlat must have two columns"))
    adaptive && bw >= n && throw(ArgumentError(
        "adaptive bandwidth $bw needs more than $n stations once the target is excluded",
    ))
    distances = TraditionalInterpolation.haversine_distance_matrix(lonlat, lonlat)
    distances .*= JOINT_EARTH_RADIUS_KM / TraditionalInterpolation.EARTH_RADIUS_KM
    @inbounds for index in 1:n
        distances[index, index] = Inf
    end
    return MixedGWR.gw_weight(distances, bw; kernel, adaptive)
end

"""
Geographically weighted slope of the satellite residual on the satellite itself, per cell.

This is the one-dimensional analogue of the coefficient `--free-satellite-coefficient` fits. With
that flag on, `JointCovariateModels.build_joint_fold_context` puts the satellite into the local
design, so the residual model becomes `r = b0(u) + b_sat(u)*y_sat + ...` and the effective anchor
is `1 + b_sat(u)` instead of a forced 1. Fitting the same weighted regression here, from an
existing run's stored matrices, bounds how far that coefficient could move the predictions without
paying for a re-run.

Returns `slope` (the weighted least-squares slope of `r = y_obs - y_sat` on `y_sat` over the
*other* stations), `satellite_mean` (the neighbourhood's weighted mean satellite value, which the
increment is measured against) and `degenerate_cells`.

`slope` is 0 wherever the local design is rank-deficient. That is overwhelmingly the dry/dry case:
every neighbour reports `y_sat = 0`, so there is no variation to regress on. Zero is exactly the
forced anchor, so a degenerate cell falls back to the model the run already fitted rather than to
an arbitrary number, and the count comes back so its share is visible rather than assumed small.

The per-hour weighted normal equations are accumulated as five matrix products rather than a solve
per cell: at 237 stations x 13471 hours the loop would be 3.2M two-by-two systems.
"""
function local_satellite_slope(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, weights::Matrix{Float64};
    variance_floor::Float64=1e-8, block::Int=512,
)
    size(y_obs) == size(y_sat) ||
        throw(DimensionMismatch("y_obs and y_sat must have the same shape"))
    n, nt = size(y_obs)
    size(weights) == (n, n) ||
        throw(DimensionMismatch("weights must be station x station"))
    slope = zeros(Float64, n, nt)
    satellite_mean = fill(NaN, n, nt)
    degenerate = 0
    transposed = permutedims(weights)
    V = Matrix{Float64}(undef, n, block)
    VX = similar(V); VX2 = similar(V); VR = similar(V); VXR = similar(V)
    for start in 1:block:nt
        stop = min(start + block - 1, nt)
        width = stop - start + 1
        v = view(V, :, 1:width); vx = view(VX, :, 1:width); vx2 = view(VX2, :, 1:width)
        vr = view(VR, :, 1:width); vxr = view(VXR, :, 1:width)
        fill!(v, 0.0); fill!(vx, 0.0); fill!(vx2, 0.0); fill!(vr, 0.0); fill!(vxr, 0.0)
        @inbounds for column in 1:width
            time = start + column - 1
            for i in 1:n
                observed = y_obs[i, time]
                x = y_sat[i, time]
                (isnan(observed) || isnan(x)) && continue
                r = observed - x
                v[i, column] = 1.0
                vx[i, column] = x
                vx2[i, column] = x * x
                vr[i, column] = r
                vxr[i, column] = x * r
            end
        end
        sum_w = transposed * v
        sum_wx = transposed * vx
        sum_wx2 = transposed * vx2
        sum_wr = transposed * vr
        sum_wxr = transposed * vxr
        @inbounds for column in 1:width
            time = start + column - 1
            for j in 1:n
                total = sum_w[j, column]
                if !(total > 0)
                    degenerate += 1
                    continue
                end
                mean_x = sum_wx[j, column] / total
                satellite_mean[j, time] = mean_x
                variance_x = sum_wx2[j, column] / total - mean_x^2
                if !(variance_x > variance_floor)
                    degenerate += 1
                    continue
                end
                mean_r = sum_wr[j, column] / total
                covariance = sum_wxr[j, column] / total - mean_x * mean_r
                slope[j, time] = covariance / variance_x
            end
        end
    end
    return (; slope, satellite_mean, degenerate_cells=degenerate)
end

"""
Replay a run's stored predictions with the satellite anchor freed locally.

The stored prediction is `p = max(y_sat + s*correction, 0)` with the anchor forced to 1. Freeing
the coefficient adds a satellite column to the local design, and for a weighted least squares
*that already carries an intercept* the fitted value moves by exactly
`b_sat(u) * (y_sat(u) - xbar(u))`, where `xbar` is the neighbourhood's weighted mean satellite
value - the intercept absorbs everything else. So the increment applied here is
`shrink * slope * (y_sat - satellite_mean)`.

That form, rather than the `a*y_sat + correction` of `verify_anchor_discount_bounds.jl`, is the
whole point: a per-cell `a` would double-count the neighbourhood mean, which the run's stored
correction already carries. It also makes the limit legible - the increment is 0 wherever the
target's satellite value is what its neighbourhood would have predicted, which is the honest
statement of what a coefficient varying in *location* can and cannot do about an error that varies
with the satellite's *state*.

`floor_at_zero` clamps the increment so the effective anchor never falls below 0: the satellite may
be removed entirely, never inverted.

`clipped_cells` counts where the replay is not exact. Where `p > 0` the applied correction is
recoverable as `p - y_sat`, so the replay is exact. Where `p == 0` the true correction is only
known to satisfy `y_sat + s*correction <= 0`; a negative increment still gives 0 and stays exact,
so only a *positive* increment on a clipped cell is uncertain, and that is what is counted. On such
a cell the freed model would return `max(y_sat + s*correction + increment, 0)`, which is at most
the `increment` this returns - so the default over-predicts there.

`suppress_clipped` takes the other end: it drops the increment on exactly those cells, returning
the stored 0. The truth lies cellwise between the two, so running both brackets what the stored
matrices can determine rather than reporting one side of it as though it were the answer.
"""
function local_anchor_prediction(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64},
    slope::Matrix{Float64}, satellite_mean::Matrix{Float64};
    shrink::Float64=1.0, floor_at_zero::Bool=false, suppress_clipped::Bool=false,
)
    out = similar(prediction)
    clipped = 0
    @inbounds for index in eachindex(prediction)
        predicted = prediction[index]
        satellite = y_sat[index]
        if isnan(predicted) || isnan(satellite)
            out[index] = NaN
            continue
        end
        mean_x = satellite_mean[index]
        increment = isnan(mean_x) ? 0.0 : shrink * slope[index] * (satellite - mean_x)
        floor_at_zero && (increment = max(increment, -satellite))
        if predicted == 0.0 && increment > 0
            clipped += 1
            suppress_clipped && (increment = 0.0)
        end
        out[index] = max(predicted + increment, 0.0)
    end
    return (; prediction=out, clipped_cells=clipped)
end

"""
Is the satellite's false alarm spatially coherent enough for a local coefficient to see it?

`satellite_quadrant_table` puts nearly the whole GWR-family gap in the gauge-dry/satellite-wet
quadrant. A coefficient that varies with *location* can only correct that quadrant where the
neighbouring stations are in it too, so this measures exactly that, per quadrant:

- `neighbour_share` - the geographically weighted share of a cell's neighbours in the same
  quadrant, averaged over the quadrant's own cells.
- `hour_share` - the share of the *other* valid stations in that quadrant in the same hour,
  averaged the same way. This is the null that matters. A false alarm is trivially more likely in
  an hour that has many of them, and a per-hour effect is one the model's local intercept already
  absorbs; only an excess over this hour share is spatial information a local coefficient could
  use.
- `lift` - `neighbour_share / hour_share`. At 1 there is no spatial structure beyond the hour.

`mean_slope`, `share_slope_negative` and `mean_increment` summarise what `local_satellite_slope`
actually produced on those cells, so the coherence measure and the replay can be read against each
other.
"""
function false_alarm_coherence_table(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, mask::AbstractMatrix,
    weights::Matrix{Float64}, slope::Matrix{Float64}, satellite_mean::Matrix{Float64};
    scheme::String, product::String, method::String, kernel::Int, bw::Float64,
    adaptive::Bool, threshold::Float64=0.1, block::Int=512,
)
    n, nt = size(y_obs)
    valid = falses(n, nt)
    @inbounds for index in eachindex(valid)
        valid[index] = !isnan(y_obs[index]) && !isnan(y_sat[index])
    end
    transposed = permutedims(weights)
    overall_n = count(mask)
    rows = NamedTuple[]
    for (name, obs_wet, sat_wet) in (
        ("dry_dry", false, false), ("dry_wet", false, true),
        ("wet_dry", true, false), ("wet_wet", true, true),
    )
        indicator = falses(n, nt)
        @inbounds for index in eachindex(valid)
            valid[index] || continue
            indicator[index] = (y_obs[index] >= threshold) == obs_wet &&
                (y_sat[index] >= threshold) == sat_wet
        end
        cell = indicator .& mask
        cell_count = count(cell)
        cell_count == 0 && continue
        neighbour_total = 0.0; hour_total = 0.0; scored = 0
        slope_total = 0.0; negative = 0; increment_total = 0.0
        indicator_block = Matrix{Float64}(undef, n, block)
        valid_block = Matrix{Float64}(undef, n, block)
        for start in 1:block:nt
            stop = min(start + block - 1, nt)
            width = stop - start + 1
            ib = view(indicator_block, :, 1:width); vb = view(valid_block, :, 1:width)
            @inbounds for column in 1:width, i in 1:n
                time = start + column - 1
                ib[i, column] = indicator[i, time] ? 1.0 : 0.0
                vb[i, column] = valid[i, time] ? 1.0 : 0.0
            end
            weighted_hits = transposed * ib
            weighted_valid = transposed * vb
            @inbounds for column in 1:width
                time = start + column - 1
                hour_valid = 0; hour_hits = 0
                for i in 1:n
                    valid[i, time] || continue
                    hour_valid += 1
                    indicator[i, time] && (hour_hits += 1)
                end
                hour_valid > 1 || continue
                for j in 1:n
                    cell[j, time] || continue
                    total = weighted_valid[j, column]
                    total > 0 || continue
                    scored += 1
                    neighbour_total += weighted_hits[j, column] / total
                    # The target is excluded from its own neighbourhood (the weight matrix has an
                    # infinite diagonal), so the hour null excludes it too - otherwise the cell
                    # would be compared against a rate it is itself counted in.
                    hour_total += (hour_hits - 1) / (hour_valid - 1)
                    slope_total += slope[j, time]
                    slope[j, time] < 0 && (negative += 1)
                    mean_x = satellite_mean[j, time]
                    isnan(mean_x) ||
                        (increment_total += slope[j, time] * (y_sat[j, time] - mean_x))
                end
            end
        end
        neighbour_share = scored > 0 ? neighbour_total / scored : NaN
        hour_share = scored > 0 ? hour_total / scored : NaN
        push!(rows, (;
            scheme, product, method, kernel, bw, adaptive, quadrant=name,
            n=cell_count, sample_share=cell_count / overall_n, scored,
            neighbour_share, hour_share,
            lift=hour_share > 0 ? neighbour_share / hour_share : NaN,
            mean_slope=scored > 0 ? slope_total / scored : NaN,
            share_slope_negative=scored > 0 ? negative / scored : NaN,
            mean_increment=scored > 0 ? increment_total / scored : NaN,
        ))
    end
    return DataFrame(rows)
end

end # module
