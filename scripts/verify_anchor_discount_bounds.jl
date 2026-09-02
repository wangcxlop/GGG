#!/usr/bin/env julia

# How much of the GWR family's loss is reachable at all by discounting the satellite anchor,
# bounded from an existing benchmark run without refitting anything.
#
#   julia --project=. scripts/verify_anchor_discount_bounds.jl
#   julia --project=. scripts/verify_anchor_discount_bounds.jl <benchmark_output_dir>
#
# `verify_dry_gate_replay.jl` established that the family's entire pooled-RMSE gap against the
# best traditional interpolator sits in one quadrant - the gauge is dry and the satellite falsely
# reports rain - and that no gate on the *correction* can reach it, because the correction is not
# what is wrong there. The residual framing predicts `y_sat + correction`; a false alarm is an
# error in the anchor the model was handed.
#
# So the question this script answers is: if the model could discount that anchor, how much would
# it get back? Three counterfactuals, in increasing order of how much they are allowed to know,
# so that a real method can be placed between them:
#
#  1. `anchor_discount` - FEASIBLE and EXACT. Predict `max(a*y_sat + correction, 0)` for a global
#     constant `a`. This is the simplest possible version of "let the satellite coefficient be
#     fitted rather than forced to 1", and a locally varying coefficient can only do better, so
#     it is a *lower* bound on that fix.
#  2. `satellite_wet_blend` - FEASIBLE. On the cells the satellite calls wet (which is knowable
#     at prediction time), blend the residual prediction toward a method that ignores the
#     satellite entirely. Bounds the more aggressive "fall back where the satellite is not
#     trusted" family of fixes.
#  3. `oracle_dry_wet` - NOT FEASIBLE, reported as a ceiling only. Replace the prediction with
#     the fallback in the false-alarm quadrant *as identified by the gauge*. No method can know
#     this, so it is the most any anchor-discounting fix could possibly achieve.
#
# Exactness of (1): the stored prediction is `p = max(y_sat + s*r, 0)`. Where `p > 0` the applied
# correction is `s*r = p - y_sat` exactly, so `max(a*y_sat + s*r, 0) = max(p - (1-a)*y_sat, 0)`.
# Where the clip bound (`p == 0`) we know `y_sat + s*r <= 0`, hence `a*y_sat + s*r <= -(1-a)*y_sat
# <= 0` for every `a <= 1` and `y_sat >= 0` - so the rescaled value is 0 there too, which is what
# the same formula returns. The identity is therefore exact over the whole grid for `a <= 1`, and
# the `a = 1` row must reproduce the stored metrics, which this script asserts.
#
# The standing caveat from the replay applies unchanged: bandwidth, kernel and the run's own
# shrinkage stay at their already-selected values, so none of this can see how the tuner would
# re-select them. These are bounds on what is worth building, not results.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("BenchmarkDiagnostics")
using Main.BenchmarkDiagnostics
using Main: BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const OBS_PATH = joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only",
)

"""Residual-mode methods: the only ones carrying a satellite anchor to discount."""
const ANCHORED_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]

"""
Methods that ignore the satellite, used as the thing to fall back toward.

`adw` is the best traditional interpolator and so the target to beat; `gwr` is the within-family
analogue, i.e. what the same estimator does when it is not handed an anchor at all.
"""
const FALLBACKS = ["adw", "gwr"]

const SWEEP = collect(0.0:0.1:1.0)

"""Matches `rain_threshold` / `wet_threshold` in the benchmark config."""
const WET_THRESHOLD = 0.1

"""Discount the satellite anchor by a global constant, keeping the applied correction."""
function anchor_discount(y_sat::Matrix{Float64}, prediction::Matrix{Float64}, a::Float64)
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        satellite = y_sat[index]
        predicted = prediction[index]
        out[index] = (isnan(satellite) || isnan(predicted)) ? NaN :
            max(predicted - (1 - a) * satellite, 0.0)
    end
    return out
end

"""Blend toward `fallback` on the cells the satellite calls wet. `lambda = 0` changes nothing."""
function satellite_wet_blend(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64}, fallback::Matrix{Float64},
    lambda::Float64, threshold::Float64,
)
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
Replace the prediction with `fallback` exactly where the gauge is dry and the satellite is wet.

Uses the gauge to choose, so it is unattainable by construction; reported as the ceiling any
anchor-discounting fix is working toward.
"""
function oracle_dry_wet(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64},
    fallback::Matrix{Float64}, threshold::Float64,
)
    out = copy(prediction)
    @inbounds for index in eachindex(prediction)
        observed = y_obs[index]
        satellite = y_sat[index]
        other = fallback[index]
        (isnan(observed) || isnan(satellite) || isnan(other)) && continue
        observed < threshold <= satellite && (out[index] = other)
    end
    return out
end

"""Pooled score plus detection at the wet threshold, as one row."""
function score_row(
    y_obs::Matrix{Float64}, prediction::Matrix{Float64}, mask::BitMatrix; base...,
)
    metrics = BenchmarkDiagnostics._metrics(y_obs, prediction, mask)
    scored = mask .& .!isnan.(prediction) .& .!isnan.(y_obs)
    event = metric_event(y_obs, prediction; mask=scored, thr=WET_THRESHOLD)
    return (; base..., metrics..., event...)
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    outdir = benchmark_diagnostics_outdir(ROOT, run_dir)
    isfile(OBS_PATH) || error("gauge observations not found: $OBS_PATH")

    println("Reading benchmark run: $run_dir")
    rows = NamedTuple[]

    for scheme in ("balanced_spatial", "random")
        scheme_dir = joinpath(run_dir, scheme)
        isdir(scheme_dir) || continue
        for product in sort(filter(name -> isdir(joinpath(scheme_dir, name)), readdir(scheme_dir)))
            product_dir = joinpath(scheme_dir, product)
            raw_path = joinpath(product_dir, "oof_raw.csv")
            isfile(raw_path) || continue

            ids, times = read_run_grid(raw_path)
            y_obs, unmatched = load_gauge_matrix(OBS_PATH, ids, times)
            unmatched == 0 || error(
                "$scheme/$product: $unmatched of $(length(times)) hours absent from $OBS_PATH",
            )
            predictions = load_prediction_matrices(product_dir, ids)
            mask = read_mask_matrix(joinpath(product_dir, "common_evaluation_mask.csv"), ids)
            y_sat = predictions["raw"]
            println("  $scheme/$product: $(count(mask)) evaluated cells")

            for method in ("raw", "idw", "adw", "tps", "gwr")
                haskey(predictions, method) || continue
                push!(rows, score_row(y_obs, predictions[method], mask;
                    scheme, product, method, counterfactual="none", fallback="",
                    parameter=NaN))
            end

            for method in ANCHORED_METHODS
                haskey(predictions, method) || continue
                prediction = predictions[method]

                # `a = 1` is the identity; check it before trusting anything the sweep says.
                identity = anchor_discount(y_sat, prediction, 1.0)
                for index in eachindex(prediction)
                    isnan(prediction[index]) && continue
                    identity[index] == prediction[index] || error(
                        "$scheme/$product/$method: a=1 did not reproduce the stored prediction " *
                        "at linear index $index ($(identity[index]) vs $(prediction[index]))",
                    )
                end

                for a in SWEEP
                    push!(rows, score_row(y_obs, anchor_discount(y_sat, prediction, a), mask;
                        scheme, product, method, counterfactual="anchor_discount",
                        fallback="", parameter=a))
                end

                for name in FALLBACKS
                    haskey(predictions, name) || continue
                    fallback = predictions[name]
                    for lambda in SWEEP
                        push!(rows, score_row(y_obs,
                            satellite_wet_blend(y_sat, prediction, fallback, lambda, WET_THRESHOLD),
                            mask; scheme, product, method,
                            counterfactual="satellite_wet_blend", fallback=name, parameter=lambda))
                    end
                    push!(rows, score_row(y_obs,
                        oracle_dry_wet(y_obs, y_sat, prediction, fallback, WET_THRESHOLD),
                        mask; scheme, product, method,
                        counterfactual="oracle_dry_wet", fallback=name, parameter=NaN))
                end
            end
        end
    end

    isempty(rows) && error("no scheme/product directories with stored predictions under $run_dir")
    table = DataFrame(rows)
    path = joinpath(outdir, "anchor_discount_bounds.csv")
    CSV.write(path, table)
    println("\nWrote $path ($(nrow(table)) rows)")

    report(table)
    return table
end

"""Print each counterfactual's best against the traditional method it has to beat."""
function report(table::DataFrame)
    for scheme in unique(table.scheme), product in unique(table.product)
        cell = filter(row -> row.scheme == scheme && row.product == product, table)
        isempty(cell) && continue
        traditional = filter(row -> row.method in ("idw", "adw", "tps"), cell)
        isempty(traditional) && continue
        target = traditional[argmin(traditional.RMSE), :]
        println("\n$scheme / $product - target: $(target.method) RMSE " *
            "$(round(target.RMSE; digits=5))")
        for method in ANCHORED_METHODS
            baseline = filter(row -> row.method == method &&
                row.counterfactual == "anchor_discount" && row.parameter == 1.0, cell)
            isempty(baseline) && continue
            println("  $method (as run): RMSE $(round(only(baseline).RMSE; digits=5)) " *
                "POD $(round(only(baseline).POD; digits=4))")
            for counterfactual in ("anchor_discount", "satellite_wet_blend", "oracle_dry_wet")
                for fallback in unique(filter(row -> row.method == method &&
                        row.counterfactual == counterfactual, cell).fallback)
                    variant = filter(row -> row.method == method &&
                        row.counterfactual == counterfactual && row.fallback == fallback, cell)
                    isempty(variant) && continue
                    best = variant[argmin(variant.RMSE), :]
                    gain = 100 * (target.RMSE - best.RMSE) / target.RMSE
                    label = isempty(fallback) ? counterfactual : "$counterfactual->$fallback"
                    parameter = isnan(best.parameter) ? "-" : string(best.parameter)
                    println("    $(rpad(label, 32)) best=$(rpad(parameter, 5)) " *
                        "RMSE $(round(best.RMSE; digits=5)) POD $(round(best.POD; digits=4)) " *
                        "[$(round(gain; digits=2))% vs $(target.method)]")
                end
            end
        end
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
