#!/usr/bin/env julia

# What a *locally* varying satellite coefficient could buy, bounded from an existing benchmark run
# without refitting anything.
#
#   julia --project=. scripts/verify_local_anchor_bound.jl
#   julia --project=. scripts/verify_local_anchor_bound.jl <benchmark_output_dir>
#
# `verify_anchor_discount_bounds.jl` swept a *global* satellite discount and found it now loses to
# `adw` on all three products while costing 5-9 points of POD, and that the one counterfactual which
# does win - `satellite_wet_blend->adw` - conditions on the satellite's *state at the cell*. The
# `--free-satellite-coefficient` flag varies the coefficient by *location* instead. This script
# measures how much of the gap that third option can reach, before a 5-hour re-run is spent on it.
#
# The replay fits the same weighted regression the flag fits, in one dimension: at every cell, the
# geographically weighted slope of `r = y_obs - y_sat` on `y_sat` over the *other* stations. For a
# weighted least squares that already carries a local intercept - which the residual model does -
# freeing the satellite coefficient moves the fitted value by exactly `b_sat * (y_sat - xbar)`,
# where `xbar` is the neighbourhood's weighted mean satellite value. That increment, times the run's
# own shrinkage, is what is added to the stored prediction.
#
# Note the shape of that increment, because it is the finding as much as the numbers are: it
# vanishes wherever the target's satellite value is what its neighbourhood would have predicted. A
# location-varying coefficient can only correct a false alarm that its *neighbours* also see.
# `false_alarm_coherence.csv` measures whether they do, against the right null - the share of the
# same hour's other stations in the same quadrant, since a purely per-hour effect is one the local
# intercept already absorbs.
#
# Exactness, and why there are three variants rather than one number. Where the stored `p > 0` the
# applied correction is `p - y_sat` exactly and the replay is exact. Where `p == 0` the correction
# is only known to satisfy `y_sat + s*correction <= 0`. A negative increment still gives 0 there, so
# only a *positive* increment on a clipped cell is undetermined - and on those cells the freed model
# returns at most the increment. `local_freed` applies it (over-predicting), `local_freed_exact`
# drops it (under-predicting), and the freed model lies cellwise between them. `clipped_cells` says
# how many cells that gap spans; it is 10-25% of the grid here, which is too many to round away.
# Forcing the slope to 0 must reproduce the stored prediction bit for bit, which this script asserts
# before reporting anything.
#
# The remaining biases all run one way, toward flattering the flag:
#
#  1. Neighbours are drawn from all of the run's stations. A real 5-fold fit sees about 80% of them,
#     so the neighbourhood here is denser than the one the model actually gets.
#  2. Kernel, bandwidth and shrinkage are frozen at the values the run already selected, so nothing
#     here can see how the tuner would re-select them - the standing caveat from the other two
#     replays - and `report` prints the best RMSE over those geometries, which is a choice made on
#     the held-out cells.
#  3. The regression is on the satellite alone. The real fit shares the local design with the
#     screened covariates, which would explain away part of what is attributed to the satellite.
#
# So this bounds `--free-satellite-coefficient` from above on every axis except the clipped cells,
# where the two variants bracket it. It is a bound on what is worth building, not a result.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("BenchmarkDiagnostics")
using Main.BenchmarkDiagnostics
using Main: BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const OBS_PATH = joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv")
const META_PATH = joinpath(STUDY_DATA, "station_meta.csv")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only",
)

"""Matches `rain_threshold` / `wet_threshold` in the benchmark config."""
const WET_THRESHOLD = 0.1

"""
Residual-mode methods, and the `joint_bandwidths.csv` group whose bandwidth stands for the model's
spatial scale.

`residual_gwr` and `mixed_gwr` carry one shared bandwidth, so the choice is forced. Joint `mgwr`
fits one per group and has no satellite group in this run - the flag is what would add it - so the
`intercept` bandwidth is used: it is the spatial scale the model resolves the residual field at,
and it is the only group present in every fold.
"""
const ANCHORED_METHODS = [
    ("residual_gwr", "residual", "gwr", "shared_all_local"),
    ("mixed_gwr", "residual", "mixed_gwr", "shared_local"),
    ("mgwr", "residual", "mgwr", "intercept"),
]

"""Traditional baselines the family has to beat, scored unchanged for reference."""
const REFERENCE_METHODS = ["raw", "idw", "adw", "tps", "gwr"]

"""Station coordinates in the run's own column order."""
function load_lonlat(path::AbstractString, ids::Vector{String})
    meta = CSV.read(path, DataFrame)
    lookup = Dict(string(row.station_id) => (Float64(row.lon), Float64(row.lat))
                  for row in eachrow(meta))
    lonlat = Matrix{Float64}(undef, length(ids), 2)
    for (index, id) in enumerate(ids)
        haskey(lookup, id) || error("station $id is in the run but not in $path")
        lonlat[index, 1], lonlat[index, 2] = lookup[id]
    end
    return lonlat
end

"""
The distinct `(kernel, bw, adaptive)` combinations the run selected for one method, with the
shrinkage that went with them.

Read out of the run rather than chosen here: the replay has no business inventing a bandwidth. A
method whose folds disagree contributes several rows, so the bound is "the best of the geometries
the tuner actually picked" rather than the best of an arbitrary sweep.
"""
function selected_geometries(bandwidths::DataFrame, scheme, product, mode, method, group)
    rows = filter(row -> row.scheme == scheme &&
        lowercase(String(row.product)) == lowercase(product) &&
        row.mode == mode && row.method == method && row.group == group, bandwidths)
    isempty(rows) && return NamedTuple[]
    out = NamedTuple[]
    for key in unique(zip(rows.kernel, rows.bw, rows.adaptive))
        kernel, bw, adaptive = key
        matching = filter(row -> row.kernel == kernel && row.bw == bw &&
            row.adaptive == adaptive, rows)
        push!(out, (; kernel=Int(kernel), bw=Float64(bw), adaptive=Bool(adaptive),
            shrink=median(Float64.(matching.shrink)), folds=nrow(matching)))
    end
    return sort(out, by=entry -> (entry.kernel, entry.bw, entry.adaptive))
end

"""Pooled score plus detection at the wet threshold, as one row. Mirrors the other two replays."""
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
    isfile(META_PATH) || error("station metadata not found: $META_PATH")
    bandwidth_path = joinpath(run_dir, "joint_bandwidths.csv")
    isfile(bandwidth_path) || error("no joint_bandwidths.csv under $run_dir")
    bandwidths = CSV.read(bandwidth_path, DataFrame)

    println("Reading benchmark run: $run_dir")
    rows = NamedTuple[]
    coherence = DataFrame[]

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
            lonlat = load_lonlat(META_PATH, ids)
            println("  $scheme/$product: $(count(mask)) evaluated cells, $(length(ids)) stations")

            for method in REFERENCE_METHODS
                haskey(predictions, method) || continue
                push!(rows, score_row(y_obs, predictions[method], mask;
                    scheme, product, method, variant="none", kernel=-1, bw=NaN,
                    adaptive=false, shrink=NaN, folds=0,
                    clipped_cells=0, degenerate_cells=0))
            end

            # Keyed on the geometry alone: the slope field depends on the weights and the data,
            # not on which method selected them, and two methods often select the same geometry.
            slope_cache = Dict{Tuple{Int,Float64,Bool},Any}()

            for (output_method, mode, method, group) in ANCHORED_METHODS
                haskey(predictions, output_method) || continue
                prediction = predictions[output_method]
                geometries = selected_geometries(
                    bandwidths, scheme, product, mode, method, group,
                )
                isempty(geometries) && (@warn(
                    "no selected bandwidth rows", scheme, product, method, group,
                ); continue)

                for geometry in geometries
                    key = (geometry.kernel, geometry.bw, geometry.adaptive)
                    fitted = get!(slope_cache, key) do
                        weights = station_weight_matrix(
                            lonlat; kernel=geometry.kernel, bw=geometry.bw,
                            adaptive=geometry.adaptive,
                        )
                        merge(local_satellite_slope(y_obs, y_sat, weights), (; weights))
                    end

                    # The forced anchor is `slope = 0`; it must return the run unchanged.
                    zeroed = local_anchor_prediction(
                        y_sat, prediction, zeros(size(prediction)), fitted.satellite_mean;
                        shrink=geometry.shrink,
                    )
                    for index in eachindex(prediction)
                        isnan(prediction[index]) && continue
                        zeroed.prediction[index] == prediction[index] || error(
                            "$scheme/$product/$output_method: a zero slope did not reproduce " *
                            "the stored prediction at linear index $index " *
                            "($(zeroed.prediction[index]) vs $(prediction[index]))",
                        )
                    end

                    base = (; scheme, product, method=output_method,
                        kernel=geometry.kernel, bw=geometry.bw, adaptive=geometry.adaptive,
                        shrink=geometry.shrink, folds=geometry.folds,
                        degenerate_cells=fitted.degenerate_cells)
                    push!(rows, score_row(y_obs, zeroed.prediction, mask;
                        base..., variant="identity", clipped_cells=zeroed.clipped_cells))
                    # Three variants, because the stored matrices do not determine one number.
                    # `local_freed` applies the whole increment and so over-predicts on the
                    # clipped cells; `local_freed_exact` drops it there and so under-predicts. The
                    # freed model lies cellwise between them. `local_freed_floor` is the same as
                    # `local_freed` with the effective anchor kept at or above 0.
                    for (name, floor_at_zero, suppress_clipped) in (
                        ("local_freed", false, false),
                        ("local_freed_exact", false, true),
                        ("local_freed_floor", true, false),
                    )
                        replayed = local_anchor_prediction(
                            y_sat, prediction, fitted.slope, fitted.satellite_mean;
                            shrink=geometry.shrink, floor_at_zero, suppress_clipped,
                        )
                        push!(rows, score_row(y_obs, replayed.prediction, mask;
                            base..., variant=name, clipped_cells=replayed.clipped_cells))
                    end

                    push!(coherence, false_alarm_coherence_table(
                        y_obs, y_sat, mask, fitted.weights, fitted.slope, fitted.satellite_mean;
                        scheme, product, method=output_method, kernel=geometry.kernel,
                        bw=geometry.bw, adaptive=geometry.adaptive, threshold=WET_THRESHOLD,
                    ))
                end
            end
        end
    end

    isempty(rows) && error("no scheme/product directories with stored predictions under $run_dir")
    table = DataFrame(rows)
    path = joinpath(outdir, "local_anchor_bound.csv")
    CSV.write(path, table)
    println("\nWrote $path ($(nrow(table)) rows)")

    coherence_table = vcat(coherence...)
    coherence_path = joinpath(outdir, "false_alarm_coherence.csv")
    CSV.write(coherence_path, coherence_table)
    println("Wrote $coherence_path ($(nrow(coherence_table)) rows)")

    report(table, coherence_table)
    return (; table, coherence=coherence_table)
end

"""Print the freed-coefficient bound against the traditional method it has to beat."""
function report(table::DataFrame, coherence::DataFrame)
    for scheme in unique(table.scheme), product in unique(table.product)
        cell = filter(row -> row.scheme == scheme && row.product == product, table)
        isempty(cell) && continue
        traditional = filter(row -> row.method in ("idw", "adw", "tps"), cell)
        isempty(traditional) && continue
        target = traditional[argmin(traditional.RMSE), :]
        println("\n$scheme / $product - target: $(target.method) RMSE " *
            "$(round(target.RMSE; digits=5)) POD $(round(target.POD; digits=4))")
        for method in unique(row.method for row in eachrow(cell) if row.variant != "none")
            asrun = filter(row -> row.method == method && row.variant == "identity", cell)
            isempty(asrun) && continue
            reference = asrun[1, :]
            println("  $method (as run): RMSE $(round(reference.RMSE; digits=5)) " *
                "POD $(round(reference.POD; digits=4))")
            for variant in ("local_freed", "local_freed_exact", "local_freed_floor")
                rows = filter(row -> row.method == method && row.variant == variant, cell)
                isempty(rows) && continue
                best = rows[argmin(rows.RMSE), :]
                gain = 100 * (target.RMSE - best.RMSE) / target.RMSE
                geometry = "k$(best.kernel)/bw$(best.bw)/$(best.adaptive ? "adap" : "fixed")"
                println("    $(rpad(variant, 18)) $(rpad(geometry, 22)) " *
                    "RMSE $(round(best.RMSE; digits=5)) POD $(round(best.POD; digits=4)) " *
                    "[$(round(gain; digits=2))% vs $(target.method)] " *
                    "clipped=$(best.clipped_cells)")
            end
        end
        wet = filter(row -> row.scheme == scheme && row.product == product &&
            row.quadrant == "dry_wet", coherence)
        isempty(wet) && continue
        entry = wet[1, :]
        println("  dry_wet coherence: neighbour $(round(entry.neighbour_share; digits=4)) " *
            "vs same-hour $(round(entry.hour_share; digits=4)) " *
            "(lift $(round(entry.lift; digits=3))), " *
            "mean slope $(round(entry.mean_slope; digits=4)), " *
            "negative $(round(100 * entry.share_slope_negative; digits=1))%")
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
