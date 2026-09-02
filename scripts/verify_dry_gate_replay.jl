#!/usr/bin/env julia

# What a dry-cell gate on the residual correction would be worth, measured from an existing
# benchmark run without refitting anything.
#
#   julia --project=. scripts/verify_dry_gate_replay.jl
#   julia --project=. scripts/verify_dry_gate_replay.jl <benchmark_output_dir>
#
# Why this is possible at all: a residual-mode prediction is `p = max(y_sat + s*r, 0)`, and the
# run stores both `p` (`oof_<method>.csv`) and `y_sat` (`oof_raw.csv`) per cell. So the applied
# correction `c = p - y_sat` is recoverable, and a gate that rescales it by `g` on the cells the
# satellite calls dry is just `max(y_sat + g*c, 0)`. At `g = 1` that identity must reproduce the
# stored prediction exactly, which is the self-check this script asserts before reporting
# anything.
#
# Two limits, both deliberate and neither repairable from stored artefacts:
#
#  1. This is a feasibility estimate, not a result. Bandwidth, kernel and the run's own
#     shrinkage stay at the values that were already selected, so the replay cannot see how the
#     tuner would re-select them once a gate exists. Only a real run can.
#  2. `g` is swept against out-of-fold predictions. Reading a winner off this sweep and
#     reporting it would be selection on the test set. The gate has to be selected on the inner
#     split inside the benchmark; this sweep only says whether it is worth building.
#
# One small inexactness is reported rather than hidden: where the `max(., 0)` clip bound
# (`p == 0`) the true correction is only known to be at most `-y_sat`, so `c` understates its
# magnitude. That cannot matter when `y_sat == 0` (every gate leaves 0 at 0), so it is confined
# to dry cells with `0 < y_sat < threshold` that were clipped. The `ambiguous_cells` column
# counts them.

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

"""Residual-mode methods: the only ones that carry a correction a gate could rescale."""
const GATED_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]

"""Scored unchanged, so the gated rows can be read against the target they have to beat."""
const REFERENCE_METHODS = ["raw", "idw", "adw", "tps", "gwr"]

const GATE_FACTORS = collect(0.0:0.1:1.0)

"""Matches `rain_threshold` / `wet_threshold` in the benchmark config."""
const DRY_THRESHOLD = 0.1

"""The run's own event thresholds, from `metrics_stratified.csv`."""
const EVENT_THRESHOLDS = [0.1, 2.5, 8.0, 16.0]

"""
`max(y_sat + g*c, 0)` with the recovered correction `c = prediction - y_sat`, and `g` applied
only where the satellite reports below `threshold`. `gate == 1` returns `prediction` unchanged.
"""
function gated_prediction(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64}, gate::Float64, threshold::Float64,
)
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        satellite = y_sat[index]
        predicted = prediction[index]
        if isnan(satellite) || isnan(predicted)
            out[index] = NaN
        else
            scale = satellite < threshold ? gate : 1.0
            out[index] = max(satellite + scale * (predicted - satellite), 0.0)
        end
    end
    return out
end

"""Dry cells inside `mask` whose correction the clip left only partially known."""
function ambiguous_cells(
    y_sat::Matrix{Float64}, prediction::Matrix{Float64}, mask::BitMatrix, threshold::Float64,
)
    total = 0
    @inbounds for index in eachindex(mask)
        mask[index] || continue
        satellite = y_sat[index]
        (isnan(satellite) || satellite <= 0 || satellite >= threshold) && continue
        prediction[index] == 0.0 && (total += 1)
    end
    return total
end

"""One method at one gate, scored overall, per rain class, and per event threshold."""
function score_rows!(
    rows::Vector{NamedTuple}, y_obs::Matrix{Float64}, prediction::Matrix{Float64},
    mask::BitMatrix; scheme, product, method, gate, ambiguous,
)
    base = (; scheme, product, method, gate, ambiguous_cells=ambiguous)
    overall = BenchmarkDiagnostics._metrics(y_obs, prediction, mask)
    push!(rows, (; base..., group="overall", level="all", threshold=NaN,
        overall..., POD=NaN, FAR=NaN, CSI=NaN))

    for (name, low, high) in RAIN_CLASSES
        stratum = BenchmarkDiagnostics._rain_mask(y_obs, mask, low, high)
        metrics = BenchmarkDiagnostics._metrics(y_obs, prediction, stratum)
        push!(rows, (; base..., group="rain_intensity", level=name, threshold=NaN,
            metrics..., POD=NaN, FAR=NaN, CSI=NaN))
    end

    scored = mask .& .!isnan.(prediction) .& .!isnan.(y_obs)
    for threshold in EVENT_THRESHOLDS
        event = metric_event(y_obs, prediction; mask=scored, thr=threshold)
        push!(rows, (; base..., group="event_threshold", level=string(threshold),
            threshold, n=count(scored), RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN,
            MSE=NaN, variance=NaN, event...))
    end
    return rows
end

"""Print the sweep as the decision it is meant to support: RMSE gain against detection loss."""
function report(table::DataFrame)
    overall = filter(row -> row.group == "overall", table)
    detection = filter(row -> row.group == "event_threshold" && row.level == "0.1", table)
    for scheme in unique(overall.scheme), product in unique(overall.product)
        cell = filter(row -> row.scheme == scheme && row.product == product, overall)
        isempty(cell) && continue
        reference = filter(row -> row.method in ("idw", "adw", "tps"), cell)
        isempty(reference) && continue
        best_reference = reference[argmin(reference.RMSE), :]
        println("\n$scheme / $product - best traditional: " *
            "$(best_reference.method) RMSE $(round(best_reference.RMSE; digits=5))")
        for method in GATED_METHODS
            gated = filter(row -> row.method == method && isfinite(row.gate), cell)
            isempty(gated) && continue
            ungated = only(filter(row -> row.gate == 1.0, gated))
            best = gated[argmin(gated.RMSE), :]
            pod_of(gate) = begin
                hit = filter(row -> row.scheme == scheme && row.product == product &&
                    row.method == method && row.gate == gate, detection)
                isempty(hit) ? NaN : only(hit).POD
            end
            improvement = 100 * (best_reference.RMSE - best.RMSE) / best_reference.RMSE
            println("  $(rpad(method, 13)) gate=1.0 RMSE $(round(ungated.RMSE; digits=5)) " *
                "POD $(round(pod_of(1.0); digits=4))  ->  " *
                "best gate=$(best.gate) RMSE $(round(best.RMSE; digits=5)) " *
                "POD $(round(pod_of(best.gate); digits=4))  " *
                "[$(round(improvement; digits=2))% vs $(best_reference.method), " *
                "$(best.ambiguous_cells) ambiguous cells]")
        end
    end
    return nothing
end

"""
Print where the gap against the best traditional interpolator actually lives.

The gate above can only touch `dry_dry` and `wet_dry` - the cells the satellite calls dry. If
the gap sits in `dry_wet` instead, no setting of the gate can reach it, and the finding is about
the residual framing rather than about the correction's magnitude.
"""
function report_quadrants(table::DataFrame)
    for scheme in unique(table.scheme), product in unique(table.product)
        cell = filter(row -> row.scheme == scheme && row.product == product &&
            row.method == "mgwr", table)
        isempty(cell) && continue
        println("\n$scheme / $product - mgwr vs $(first(cell.reference)), by (gauge, satellite):")
        for row in eachrow(cell)
            println("  $(rpad(row.quadrant, 8)) share=$(rpad(round(row.sample_share; digits=4), 7)) " *
                "mean_sat=$(rpad(round(row.mean_satellite; digits=3), 6)) " *
                "MSE $(rpad(round(row.MSE; digits=4), 8)) vs " *
                "$(rpad(round(row.reference_MSE; digits=4), 8)) " *
                "(satellite $(rpad(round(row.satellite_MSE; digits=4), 8)))  " *
                "gap contrib $(round(row.mse_gap_contribution; digits=5))")
        end
    end
    return nothing
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    outdir = benchmark_diagnostics_outdir(ROOT, run_dir)
    isfile(OBS_PATH) || error("gauge observations not found: $OBS_PATH")

    println("Reading benchmark run: $run_dir")
    rows = NamedTuple[]
    quadrants = DataFrame[]

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
            println("  $scheme/$product: $(length(ids)) stations x $(length(times)) hours, " *
                "$(count(mask)) evaluated cells")

            for method in REFERENCE_METHODS
                haskey(predictions, method) || continue
                score_rows!(rows, y_obs, predictions[method], mask;
                    scheme, product, method, gate=NaN, ambiguous=0)
            end

            # Written alongside the sweep because it is what explains the sweep's result: the
            # gate can only reach the `dry_dry` quadrant, and that is not where the loss is.
            push!(quadrants, satellite_quadrant_table(
                y_obs, predictions, mask; scheme, product, threshold=DRY_THRESHOLD,
            ))

            for method in GATED_METHODS
                haskey(predictions, method) || continue
                prediction = predictions[method]
                # The identity the whole replay rests on, checked before it is used.
                replayed = gated_prediction(y_sat, prediction, 1.0, DRY_THRESHOLD)
                for index in eachindex(prediction)
                    isnan(prediction[index]) && continue
                    replayed[index] == prediction[index] || error(
                        "$scheme/$product/$method: gate=1 did not reproduce the stored " *
                        "prediction at linear index $index " *
                        "($(replayed[index]) vs $(prediction[index]))",
                    )
                end
                ambiguous = ambiguous_cells(y_sat, prediction, mask, DRY_THRESHOLD)
                for gate in GATE_FACTORS
                    score_rows!(rows, y_obs,
                        gated_prediction(y_sat, prediction, gate, DRY_THRESHOLD), mask;
                        scheme, product, method, gate, ambiguous)
                end
            end
        end
    end

    isempty(rows) && error("no scheme/product directories with stored predictions under $run_dir")
    table = DataFrame(rows)
    path = joinpath(outdir, "dry_gate_replay.csv")
    CSV.write(path, table)
    println("\nWrote $path ($(nrow(table)) rows)")

    quadrant = vcat(quadrants...)
    quadrant_path = joinpath(outdir, "satellite_quadrant.csv")
    CSV.write(quadrant_path, quadrant)
    println("Wrote $quadrant_path ($(nrow(quadrant)) rows)")

    report(table)
    report_quadrants(quadrant)
    return (; sweep=table, quadrant)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
