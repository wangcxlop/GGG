#!/usr/bin/env julia

# What a blending weight that varies *by band* could buy, bounded from an existing benchmark run
# without refitting anything.
#
#   julia --project=. scripts/verify_banded_blend_bounds.jl
#   julia --project=. scripts/verify_banded_blend_bounds.jl <benchmark_output_dir>
#
# `verify_anchor_discount_bounds.jl` found `satellite_wet_blend->adw` to be the only counterfactual
# that beats `adw` at all, and `--satellite-wet-blend` then shipped it as a tuned method: one
# scalar `lambda` per fold, applied uniformly to every cell the satellite calls wet. That fixed the
# pooled gate and broke the other one. On the canonical baseline the unblended family beats the
# best traditional interpolator on the heavy stratum by 4.2% / 4.3% / 5.8% (FY4B / GPM / GSMaP),
# and blending at the selected lambda ~ 0.7-0.8 drops that to 0.4% / 0.8% / 2.2% - because the
# same weight that suppresses a false alarm also mixes 75% of `adw` into the heavy cells, which is
# exactly where the family was winning. `assess_gwr_claim` wants >= 5% on heavy against *every*
# baseline, so no run so far passes both gates at once.
#
# This script asks whether one weight per band closes that. Two conditioning axes, swept
# separately rather than crossed:
#
#  1. `satellite_intensity` - the band is read off `y_sat` at the cell. Knowable at prediction
#     time, so this stays a model and not an oracle, exactly as the wet/dry test already is.
#  2. `nearest_train_km` - the band is read off the target station's distance to the nearest
#     training station in its own fold. `metrics_pooled.csv` says the family only overtakes `adw`
#     in the sparse 50-100 km band while losing 9% in the dense 0-20 km band, so this is the other
#     axis with visible signal. Reproduced here the way the run computes it
#     (`InterpolationBenchmarkRun.jl:646`), from `split_common.csv` and `station_meta.csv`.
#
# Why the sweep is small, and why that is exact rather than a shortcut. The bands partition the
# cells, and blending never changes which cells are finite - a NaN fallback leaves the cell at its
# unblended value, and a NaN prediction or satellite propagates regardless of lambda - so the
# scored cell population does not move as lambda moves. Pooled SSE is therefore additive over the
# bands and each band's lambda minimises independently: 3 x 11 evaluations rather than 11^3, and
# that separable optimum *is* the joint optimum. `audit_separability` does not take this on trust -
# it enumerates a coarse grid over all bands at once, scores each combination through the same
# `_metrics` the results use, and checks the joint argmin against the separable pick.
#
# The standing caveat from the other three replays applies unchanged, and is sharper here.
# Bandwidth, kernel and shrinkage stay at the values the run already selected, so nothing here can
# see how the tuner would re-select them; and the lambdas are picked on the held-out cells, so
# every number below is a *ceiling*. `_satwetblend` gave up 0.04 (GPM) / 0.17 (GSMaP) / 0.40
# (FY4B) points between this kind of ceiling and honest out-of-sample tuning with ONE parameter.
# Three or four parameters will give up more. A bound that only just clears 5% has not cleared it.

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

"""Residual-mode methods: the only ones carrying a satellite anchor to blend away."""
const ANCHORED_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]

"""
The blend target.

`gwr` was the other fallback `anchor_discount_bounds.csv` swept and it lost at every lambda on
every cell, so there is no reason to re-sweep it here.
"""
const FALLBACK_METHOD = "adw"

"""Scored unchanged for reference. `assess_gwr_claim` grades against all three of the first ones."""
const TRADITIONAL_METHODS = ["idw", "adw", "tps"]
const REFERENCE_METHODS = ["raw", "idw", "adw", "tps", "gwr"]

"""Matches `rain_threshold` / `wet_threshold` in the benchmark config."""
const WET_THRESHOLD = 0.1

"""The grid `BLEND_LAMBDAS` offers, so the bound and the shipped method are comparable."""
const SWEEP = collect(0.0:0.1:1.0)

"""Coarse grid for the joint-vs-separable audit. 3^bands combinations, each scored for real."""
const COARSE_SWEEP = [0.0, 0.5, 1.0]

"""
Bands of the satellite value at the cell, above the wet threshold.

Edges are the ones `append_stratified_metrics!` (`InterpolationBenchmarkMetrics.jl:48`) already
stratifies the *observations* on, read off `y_sat` instead so the band is knowable at prediction
time. Cells the satellite calls dry are in no band and are never touched, as now.
"""
const INTENSITY_BANDS = [("light", 0.1, 2.5), ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf)]

"""
Bands of the target station's distance to its fold's nearest training station.

Same edges as the `nearest_train_km` strata in `metrics_pooled.csv`, so a band here lines up with
a row there.
"""
const DISTANCE_BANDS = [
    ("0_20", -Inf, 20.0), ("20_50", 20.0, 50.0),
    ("50_100", 50.0, 100.0), ("100_plus", 100.0, Inf),
]

const MAX_BANDS = 4

# ------------------------------------------------------------------------------ run inputs

"""Station coordinates in the run's own column order. Mirrors `verify_local_anchor_bound.jl`."""
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
Each station's distance to the nearest station outside its own fold.

This is what `_run_fold!` records as `nearest_train_distance` (`InterpolationBenchmarkRun.jl:646`):
a station is held out in exactly one fold, and its distance is measured to that fold's training
set. Computed on `TraditionalInterpolation`'s own radius, which is the one the run used.
"""
function nearest_train_km(lonlat::Matrix{Float64}, fold_of::Vector{Int})
    distances = Main.TraditionalInterpolation.haversine_distance_matrix(lonlat, lonlat)
    out = fill(NaN, length(fold_of))
    for target in eachindex(fold_of)
        best = Inf
        for other in eachindex(fold_of)
            fold_of[other] == fold_of[target] && continue
            best = min(best, distances[other, target])
        end
        out[target] = best
    end
    return out
end

# ------------------------------------------------------------------------------- blending

"""
Blend toward the fallback with one weight per band. A `band` of 0 means the cell is never touched.

The twin of `satellite_wet_blend_prediction` (`src/InterpolationBenchmarkPredictors.jl:291`) with
the scalar generalised to a vector, and it keeps that function's NaN rule: a NaN fallback leaves
the cell at its unblended value rather than poisoning it, so the scored population cannot move
with the weights. That invariance is what makes the per-band minimisation exact.
"""
function banded_blend(
    prediction::Vector{Float64}, fallback::Vector{Float64},
    band::Vector{Int}, lambdas::Vector{Float64},
)
    out = similar(prediction)
    @inbounds for index in eachindex(prediction)
        b = band[index]
        other = fallback[index]
        if b == 0 || isnan(other)
            out[index] = prediction[index]
        else
            lambda = lambdas[b]
            out[index] = (1 - lambda) * prediction[index] + lambda * other
        end
    end
    return out
end

"""
Matrix-form `satellite_wet_blend`, copied from `verify_anchor_discount_bounds.jl:88`.

Kept here as an independent implementation on purpose: the constant-weight identity check compares
the banded path against it, and a check against a shared helper would be no check at all.
"""
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

"""Band index per scored cell from the satellite value; 0 below the wet threshold."""
function intensity_bands(satellite::Vector{Float64})
    out = zeros(Int, length(satellite))
    @inbounds for index in eachindex(satellite)
        value = satellite[index]
        (isnan(value) || value < WET_THRESHOLD) && continue
        for (position, (_, low, high)) in enumerate(INTENSITY_BANDS)
            if low <= value < high
                out[index] = position
                break
            end
        end
    end
    return out
end

"""
Band index per scored cell from the station's nearest-training-station distance.

Still gated on the satellite calling the cell wet, so the two axes differ only in what the weight
is conditioned on and not in which cells it reaches.
"""
function distance_bands(satellite::Vector{Float64}, distance::Vector{Float64})
    out = zeros(Int, length(satellite))
    @inbounds for index in eachindex(satellite)
        value = satellite[index]
        (isnan(value) || value < WET_THRESHOLD) && continue
        km = distance[index]
        isnan(km) && continue
        for (position, (_, low, high)) in enumerate(DISTANCE_BANDS)
            if low <= km < high
                out[index] = position
                break
            end
        end
    end
    return out
end

# -------------------------------------------------------------------------------- scoring

"""
Pooled score, the four rain strata, and detection at the wet threshold, for one flat prediction.

Reuses `BenchmarkDiagnostics._metrics` rather than restating the arithmetic; the vectors are
reshaped to one-column matrices because that is the signature it advertises. `rain_masks` is
passed in because it depends only on the observations, which do not change across a sweep.
"""
function score_flat(
    observed::Vector{Float64}, prediction::Vector{Float64},
    rain_masks::Vector{<:AbstractMatrix},
)
    scored = reshape(trues(length(observed)), :, 1)
    obs = reshape(observed, :, 1)
    pred = reshape(prediction, :, 1)
    overall = BenchmarkDiagnostics._metrics(obs, pred, scored)
    event = metric_event(obs, pred; mask=scored, thr=WET_THRESHOLD)
    strata = NamedTuple()
    for (position, (name, _, _)) in enumerate(RAIN_CLASSES)
        stratum = BenchmarkDiagnostics._metrics(obs, pred, rain_masks[position])
        strata = merge(strata, (;
            Symbol("RMSE_", name) => stratum.RMSE, Symbol("n_", name) => stratum.n,
        ))
    end
    return merge((;
        n=overall.n, RMSE=overall.RMSE, MAE=overall.MAE, Bias=overall.Bias,
        r=overall.r, MSE=overall.MSE,
    ), strata, event)
end

"""Sum of squared error inside one band, which is what the per-band minimisation compares."""
function band_sse(
    observed::Vector{Float64}, prediction::Vector{Float64}, fallback::Vector{Float64},
    band::Vector{Int}, target::Int, lambda::Float64,
)
    total = 0.0
    count = 0
    @inbounds for index in eachindex(band)
        band[index] == target || continue
        other = fallback[index]
        value = isnan(other) ? prediction[index] :
            (1 - lambda) * prediction[index] + lambda * other
        difference = value - observed[index]
        total += difference * difference
        count += 1
    end
    return total, count
end

"""
The per-band minimising weights, and the sweep they were chosen from.

Ties go to the smaller weight, matching `select_blend_lambda`
(`src/InterpolationBenchmarkPredictors.jl:345`): at equal error the answer that departs less from
the fitted model is the one to prefer, and it keeps the choice off the order of the grid.
"""
function select_banded_lambdas(
    observed::Vector{Float64}, prediction::Vector{Float64}, fallback::Vector{Float64},
    band::Vector{Int}, n_bands::Int; sweep::Vector{Float64}=SWEEP,
)
    chosen = zeros(Float64, n_bands)
    rows = NamedTuple[]
    for target in 1:n_bands
        best_sse = Inf
        best_lambda = 0.0
        for lambda in sweep
            sse, count = band_sse(observed, prediction, fallback, band, target, lambda)
            push!(rows, (; band_index=target, lambda, band_sse=sse, band_n=count,
                band_RMSE=count > 0 ? sqrt(sse / count) : NaN))
            if count > 0 && sse < best_sse - 1e-12
                best_sse = sse
                best_lambda = lambda
            end
        end
        chosen[target] = best_lambda
    end
    return chosen, rows
end

# ------------------------------------------------------------------------------- auditing

"""
Check the separable pick against a joint enumeration, through the same scorer the results use.

The additivity argument is only as good as the invariance it rests on, so it is tested rather than
asserted: every combination of `COARSE_SWEEP` over the bands is materialised and scored, and the
joint minimum must be the combination the per-band minimisation picks out of the same coarse grid.
"""
function audit_separability(
    observed::Vector{Float64}, prediction::Vector{Float64}, fallback::Vector{Float64},
    band::Vector{Int}, n_bands::Int, label::AbstractString,
)
    scored = reshape(trues(length(observed)), :, 1)
    obs = reshape(observed, :, 1)
    separable, _ = select_banded_lambdas(
        observed, prediction, fallback, band, n_bands; sweep=COARSE_SWEEP,
    )
    best_rmse = Inf
    best_combo = Float64[]
    for combo in Iterators.product(ntuple(_ -> COARSE_SWEEP, n_bands)...)
        lambdas = collect(Float64, combo)
        blended = banded_blend(prediction, fallback, band, lambdas)
        rmse = BenchmarkDiagnostics._metrics(obs, reshape(blended, :, 1), scored).RMSE
        if rmse < best_rmse - 1e-12
            best_rmse = rmse
            best_combo = lambdas
        end
    end
    separable == best_combo || error(
        "separability audit failed for $label: per-band pick $separable but the joint minimum " *
        "over the same coarse grid is $best_combo (RMSE $best_rmse)",
    )
    println("  separability audit ($label): joint minimum over " *
            "$(length(COARSE_SWEEP)^n_bands) combinations matches the per-band pick $best_combo")
    return best_combo
end

"""
Compare the constant-weight rows against `anchor_discount_bounds.csv`, when that file is present.

Same harness, same stored matrices, same weights: the two must agree. A disagreement means one of
the two reading paths has drifted, which is worth knowing before any of this is believed.
"""
function cross_check_constant(stored::Union{Nothing,DataFrame}, table::DataFrame)
    stored === nothing && return nothing
    worst = 0.0
    compared = 0
    for row in eachrow(table)
        row.variant == "constant" || continue
        match = filter(other ->
            other.scheme == row.scheme && other.product == row.product &&
            other.method == row.method && other.counterfactual == "satellite_wet_blend" &&
            other.fallback == FALLBACK_METHOD && other.parameter == row.lambda_1, stored)
        nrow(match) == 1 || continue
        worst = max(worst, abs(only(match.RMSE) - row.RMSE))
        compared += 1
    end
    compared == 0 && return nothing
    if worst > 1e-9
        @warn "constant-weight rows disagree with anchor_discount_bounds.csv" compared worst
    else
        println("\nCross-check: $compared constant-weight rows agree with " *
                "anchor_discount_bounds.csv to $(round(worst; sigdigits=3)) RMSE")
    end
    return worst
end

# ------------------------------------------------------------------------------------ run

"""Pad a weight vector out to the widest axis, so every row carries the same columns."""
function lambda_columns(lambdas::Vector{Float64})
    padded = fill(NaN, MAX_BANDS)
    padded[1:length(lambdas)] = lambdas
    return (; lambda_1=padded[1], lambda_2=padded[2], lambda_3=padded[3], lambda_4=padded[4])
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    outdir = benchmark_diagnostics_outdir(ROOT, run_dir)
    isfile(OBS_PATH) || error("gauge observations not found: $OBS_PATH")
    isfile(META_PATH) || error("station metadata not found: $META_PATH")

    stored_path = joinpath(outdir, "anchor_discount_bounds.csv")
    stored = isfile(stored_path) ? CSV.read(stored_path, DataFrame) : nothing
    stored === nothing &&
        @warn "no anchor_discount_bounds.csv to cross-check against" stored_path

    println("Reading benchmark run: $run_dir")
    rows = NamedTuple[]
    sweep_rows = NamedTuple[]
    audited = false

    for scheme in ("balanced_spatial", "random")
        scheme_dir = joinpath(run_dir, scheme)
        isdir(scheme_dir) || continue
        split_path = joinpath(scheme_dir, "split_common.csv")
        isfile(split_path) || (@warn("no split_common.csv", scheme); continue)
        fold_map = read_fold_map(split_path)

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

            # One scored population for every method, so the heavy comparison the claim gate makes
            # is between numbers computed on identical cells. The shared mask already guarantees
            # each of MASK_METHODS is finite here; the assertion is what makes that a fact rather
            # than an assumption inherited from another file.
            scored_mask = mask .& .!isnan.(y_obs)
            indices = findall(scored_mask)
            for method in vcat(REFERENCE_METHODS, ANCHORED_METHODS)
                haskey(predictions, method) || continue
                bad = count(index -> isnan(predictions[method][index]), indices)
                bad == 0 || error(
                    "$scheme/$product: $method is NaN on $bad of $(length(indices)) cells inside " *
                    "the shared evaluation mask, which the mask is supposed to preclude",
                )
            end

            lonlat = load_lonlat(META_PATH, ids)
            fold_of = [get(fold_map, id, 0) for id in ids]
            any(==(0), fold_of) && error("$scheme: stations missing from $split_path")
            distance_of_station = nearest_train_km(lonlat, fold_of)

            observed = [y_obs[index] for index in indices]
            satellite = [y_sat[index] for index in indices]
            distance = [distance_of_station[index[1]] for index in indices]
            fallback = [predictions[FALLBACK_METHOD][index] for index in indices]

            all_scored = reshape(trues(length(observed)), :, 1)
            rain_masks = [BenchmarkDiagnostics._rain_mask(
                reshape(observed, :, 1), all_scored, low, high,
            ) for (_, low, high) in RAIN_CLASSES]

            axes_of_run = [
                ("satellite_intensity", intensity_bands(satellite),
                 [name for (name, _, _) in INTENSITY_BANDS]),
                ("nearest_train_km", distance_bands(satellite, distance),
                 [name for (name, _, _) in DISTANCE_BANDS]),
            ]
            println("  $scheme/$product: $(length(indices)) scored cells, " *
                    "$(count(>(0), axes_of_run[1][2])) satellite-wet")

            base = (; scheme, product)
            for method in REFERENCE_METHODS
                haskey(predictions, method) || continue
                flat = [predictions[method][index] for index in indices]
                push!(rows, merge((; base..., method, variant="reference", axis="",
                    band_labels="", lambda_columns(Float64[])...),
                    score_flat(observed, flat, rain_masks)))
            end

            for method in ANCHORED_METHODS
                haskey(predictions, method) || continue
                prediction = [predictions[method][index] for index in indices]
                intensity_band = axes_of_run[1][2]

                # A zero weight vector is the identity; check it before trusting the sweep.
                for (axis, band, labels) in axes_of_run
                    identity = banded_blend(prediction, fallback, band, zeros(length(labels)))
                    for position in eachindex(prediction)
                        identity[position] == prediction[position] || error(
                            "$scheme/$product/$method/$axis: a zero weight did not reproduce the " *
                            "stored prediction at scored index $position " *
                            "($(identity[position]) vs $(prediction[position]))",
                        )
                    end
                end

                # A constant weight vector must reproduce the shipped blend's own operation, which
                # is written independently above rather than shared with the banded path.
                for lambda in (0.3, 0.7)
                    banded = banded_blend(prediction, fallback, intensity_band,
                                          fill(lambda, length(INTENSITY_BANDS)))
                    reference = satellite_wet_blend(
                        y_sat, predictions[method], predictions[FALLBACK_METHOD],
                        lambda, WET_THRESHOLD,
                    )
                    for (position, index) in enumerate(indices)
                        banded[position] == reference[index] || error(
                            "$scheme/$product/$method: the banded path and satellite_wet_blend " *
                            "disagree at lambda $lambda, scored index $position " *
                            "($(banded[position]) vs $(reference[index]))",
                        )
                    end
                end

                push!(rows, merge((; base..., method, variant="unblended", axis="",
                    band_labels="", lambda_columns(Float64[])...),
                    score_flat(observed, prediction, rain_masks)))

                for lambda in SWEEP
                    blended = banded_blend(prediction, fallback, intensity_band,
                                           fill(lambda, length(INTENSITY_BANDS)))
                    push!(rows, merge((; base..., method, variant="constant",
                        axis="satellite_intensity", band_labels="all",
                        lambda_columns(fill(lambda, MAX_BANDS))...),
                        score_flat(observed, blended, rain_masks)))
                end

                for (axis, band, labels) in axes_of_run
                    n_bands = length(labels)
                    if !audited
                        audit_separability(observed, prediction, fallback, band, n_bands,
                                           "$scheme/$product/$method/$axis")
                        audited = true
                    end
                    chosen, band_rows = select_banded_lambdas(
                        observed, prediction, fallback, band, n_bands,
                    )
                    for row in band_rows
                        push!(sweep_rows, merge((; base..., method, axis,
                            band_label=labels[row.band_index]), row))
                    end
                    blended = banded_blend(prediction, fallback, band, chosen)
                    push!(rows, merge((; base..., method, variant="banded", axis,
                        band_labels=join(labels, "|"), lambda_columns(chosen)...),
                        score_flat(observed, blended, rain_masks)))
                end
            end
        end
    end

    isempty(rows) && error("no scheme/product directories with stored predictions under $run_dir")
    table = DataFrame(rows)
    path = joinpath(outdir, "banded_blend_bounds.csv")
    CSV.write(path, table)
    println("\nWrote $path ($(nrow(table)) rows)")

    sweep_table = DataFrame(sweep_rows)
    sweep_path = joinpath(outdir, "banded_blend_sweep.csv")
    CSV.write(sweep_path, sweep_table)
    println("Wrote $sweep_path ($(nrow(sweep_table)) rows)")

    cross_check_constant(stored, table)
    report(table)
    return (; table, sweep=sweep_table)
end

"""
Print each variant against both gates at once, which is the thing no existing artefact does.

`heavy` is reported as the *minimum* relative improvement over the three traditional baselines,
because that is what `assess_gwr_claim` gates on (`InterpolationBenchmarkMetrics.jl:297`). `tps` is
usually the binding one, and quoting the `adw` figure alone would flatter every row here.
"""
function report(table::DataFrame)
    for scheme in unique(table.scheme), product in unique(table.product)
        cell = filter(row -> row.scheme == scheme && row.product == product, table)
        isempty(cell) && continue
        traditional = filter(row -> row.method in TRADITIONAL_METHODS &&
            row.variant == "reference", cell)
        nrow(traditional) == length(TRADITIONAL_METHODS) || continue
        target = traditional[argmin(traditional.RMSE), :]
        heavy_baselines = traditional.RMSE_heavy

        println("\n$scheme / $product - pooled target: $(target.method) " *
                "RMSE $(round(target.RMSE; digits=5)); heavy baselines " *
                join([string(row.method, " ", round(row.RMSE_heavy; digits=3))
                      for row in eachrow(traditional)], ", "))

        for method in ANCHORED_METHODS
            rows = filter(row -> row.method == method, cell)
            isempty(rows) && continue
            for variant in ("unblended", "constant", "banded")
                subset = filter(row -> row.variant == variant, rows)
                isempty(subset) && continue
                axes_present = variant == "banded" ? unique(subset.axis) : [""]
                for axis in axes_present
                    candidates = variant == "banded" ?
                        filter(row -> row.axis == axis, subset) : subset
                    best = candidates[argmin(candidates.RMSE), :]
                    pooled = 100 * (target.RMSE - best.RMSE) / target.RMSE
                    heavy = 100 * minimum((heavy_baselines .- best.RMSE_heavy) ./ heavy_baselines)
                    weights = join([isnan(value) ? "-" : string(round(value; digits=1))
                                    for value in (best.lambda_1, best.lambda_2,
                                                  best.lambda_3, best.lambda_4)], "/")
                    label = variant == "banded" ? "banded:$axis" : variant
                    gate = heavy >= 5.0 ? "HEAVY OK" : "heavy short"
                    println("  $(rpad(method, 13)) $(rpad(label, 28)) " *
                            "lambda $(rpad(weights, 20)) " *
                            "RMSE $(round(best.RMSE; digits=5)) [$(round(pooled; digits=2))% " *
                            "vs $(target.method)]  heavy $(round(heavy; digits=2))% " *
                            "($gate)  POD $(round(best.POD; digits=4))")
                end
            end
        end
    end
    println("\nBoth gates must hold at once: pooled >= 0% against the best traditional method " *
            "and heavy >= 5% against the worst-case one. These are ceilings chosen on held-out " *
            "cells - leave margin for what honest per-fold tuning gives up.")
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
