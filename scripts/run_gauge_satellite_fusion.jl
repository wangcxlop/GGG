#!/usr/bin/env julia

# What does combining gauges with satellite precipitation gain, which interpolator and which fusion
# method do best, and which product - or merge of products - makes the best anchor?
#
#   julia --project=. scripts/run_gauge_satellite_fusion.jl
#   julia --project=. scripts/run_gauge_satellite_fusion.jl <benchmark_output_dir>
#
# Reads a finished `run_interpolation_benchmark.jl full --nested-covariates --satellite-wet-blend
# --blend-axis agreement_envelope --fused-anchor` run and the `run_claim_reassessment.jl` tables
# written for it; fits nothing. Writes CSVs to output/gauge_satellite_fusion/, which
# `py -3.13 scripts/plot_gauge_satellite_fusion.py` draws.
#
# Every "best method" read off these tables is a choice among ~50 configurations made on held-out
# scores. `auto` (a GWR-family method chosen inside each training fold) is reported beside them as
# the selection-free answer, and the Holm columns correct across the whole family of comparisons.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Random, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("BenchmarkDiagnostics", "TableIO")
using Main.BenchmarkDiagnostics
using Main.TableIO: write_csv_atomic

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_" *
        "blendagrenv_fusedanchor",
)
const OUTDIR = joinpath(ROOT, "output", "gauge_satellite_fusion")
const SCHEME = "balanced_spatial"
const SEED = 20260913
const BOOTSTRAP_REPS = 2000

const SINGLE_PRODUCTS = ["FY4B", "GPM", "GSMaP"]
const MERGED_PRODUCTS = ["MERGED_MEAN", "MERGED_OLS"]
const ANCHORS = vcat(SINGLE_PRODUCTS, MERGED_PRODUCTS)
const INTERPOLATORS = ["idw", "adw", "tps", "gwr"]
const FUSION_METHODS = [
    "residual_gwr", "mixed_gwr", "mgwr",
    "blend_residual_gwr", "blend_mixed_gwr", "blend_mgwr",
    "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr",
]
# Methods compared across anchors on one shared mask, and scored for detection skill. Fixed here
# rather than picked per anchor from its results, so the anchor comparison is not a comparison of
# each anchor's luckiest method.
const PANEL_METHODS = ["raw", "mgwr", "blend_mgwr", "blend_agrenv_mgwr", "auto"]
const EVENT_THRESHOLDS = [0.1, 2.5, 8.0, 16.0]
# The gauge-only reference every fusion method is measured against: the best interpolator on
# pooled RMSE in every earlier run, and the blends' fallback.
const REFERENCE = "adw"

"""The `MGERConfig` the full benchmark run used, so the common station/time grid matches."""
function study_config(outdir::AbstractString)
    return MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA, "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2024_full_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2024_full_aligned.csv"),
        ),
        outdir=outdir,
        rain_threshold=0.1,
        analysis_start=DateTime(2022, 1, 1, 9),
        analysis_end=DateTime(2025, 1, 1, 8),
        expected_common_time_count=13471,
    )
end

_pooled(metrics) = filter(row -> ismissing(row.fold), metrics)

"""One value from the pooled metrics table, NaN when the row is absent."""
function _value(index, scheme, product, method, group, level, column)
    row = get(index, (scheme, product, method, group, level), nothing)
    return row === nothing ? NaN : Float64(row[column])
end

"""
Overall accuracy and detection skill per (scheme, product, method), with RMSE also expressed as a
relative improvement over the product's own satellite field (`raw`) and over gauge-only `adw`.
"""
function method_summary(pooled::DataFrame)
    index = Dict((row.scheme, row.product, row.method, row.group, row.level) => row
                 for row in eachrow(pooled))
    rows = NamedTuple[]
    for key in unique(select(pooled, :scheme, :product, :method)) |> eachrow
        scheme, product, method = String(key.scheme), String(key.product), String(key.method)
        value(group, level, column) = _value(index, scheme, product, method, group, level, column)
        rmse = value("overall", "all", :RMSE)
        raw = _value(index, scheme, product, "raw", "overall", "all", :RMSE)
        reference = _value(index, scheme, product, REFERENCE, "overall", "all", :RMSE)
        events = (; (Symbol("$(metric)_$(threshold)") =>
            value("event_threshold", string(threshold), metric)
            for threshold in EVENT_THRESHOLDS for metric in (:POD, :FAR, :CSI))...)
        push!(rows, merge((;
            scheme, product, method,
            n=Int(value("overall", "all", :n)), coverage=value("overall", "all", :coverage),
            RMSE=rmse, MAE=value("overall", "all", :MAE), Bias=value("overall", "all", :Bias),
            r=value("overall", "all", :r),
            RMSE_improvement_vs_raw=(raw - rmse) / raw,
            RMSE_improvement_vs_adw=(reference - rmse) / reference,
        ), events))
    end
    return sort!(DataFrame(rows), [:scheme, :product, :method])
end

"""Accuracy by rain intensity, distance to the nearest training gauge and year, beside `adw`."""
function stratified_summary(pooled::DataFrame)
    groups = ("rain_intensity", "nearest_train_km", "year")
    table = filter(row -> row.group in groups && row.scheme == SCHEME && row.n > 0, pooled)
    reference = Dict((row.product, row.group, row.level) => row.RMSE
                     for row in eachrow(table) if row.method == REFERENCE)
    return sort!(DataFrame([(;
        scheme=String(row.scheme), product=String(row.product), method=String(row.method),
        group=String(row.group), level=String(row.level), n=Int(row.n),
        RMSE=row.RMSE, MAE=row.MAE, Bias=row.Bias, r=row.r,
        RMSE_improvement_vs_adw=let base = get(reference, (row.product, row.group, row.level), NaN)
            (base - row.RMSE) / base
        end,
    ) for row in eachrow(table)]), [:product, :method, :group, :level])
end

"""
Every fusion method on every anchor against `adw`, from the claim re-assessment's paired bootstrap.

`pvalue_holm_family` corrects across all (product, method) pairs within a stratum, which is the
family a "best fusion method" is actually chosen from; the re-assessment's own `pvalue_holm` only
corrects across the three baselines.
"""
function fusion_vs_gauge_only(comparisons::DataFrame)
    table = filter(row -> row.baseline == REFERENCE && row.scheme == SCHEME, comparisons)
    table = select(table, :product, :method, :stratum, :n, :RMSE_baseline, :RMSE_gwr,
        :delta_RMSE, :relative_improvement, :ci_low, :ci_high, :pvalue)
    rename!(table, :RMSE_gwr => :RMSE_method, :RMSE_baseline => :RMSE_adw)
    table.relative_ci_low = table.ci_low ./ table.RMSE_adw
    table.relative_ci_high = table.ci_high ./ table.RMSE_adw
    table.pvalue_holm_family = fill(NaN, nrow(table))
    for group in groupby(table, :stratum)
        group.pvalue_holm_family .= _holm_adjust(group.pvalue)
    end
    return sort!(table, [:stratum, :product, :method])
end

"""A `station x time` prediction matrix per requested method, for one product directory."""
function load_methods(product_dir::AbstractString, ids::Vector{String}, methods)
    return Dict(method => read_wide_matrix(joinpath(product_dir, "oof_$(method).csv"), ids)
                for method in methods if isfile(joinpath(product_dir, "oof_$(method).csv")))
end

"""
Paired day-block bootstrap of RMSE for the interpolators against each other, on one product's mask.

Interpolators never read the satellite, so their predictions are the same whichever product
directory they are read from; only the evaluation mask differs. GPM's is used because GPM and
GSMaP share it and it is the larger of the two masks.
"""
function interpolation_ranking(cfg, grid, predictions, mask)
    rows = NamedTuple[]
    for method in ("adw", "gwr")
        for row in paired_bootstrap_rows(cfg, SCHEME, "GPM", grid.times, grid.Y_obs,
                predictions, mask; method, pairwise_mask=true)
            row.baseline == row.method && continue
            push!(rows, row)
        end
    end
    return DataFrame(rows)
end

"""
The anchors compared with each other on the cells every one of them was scored on.

Each derived product's mask needs all three source products finite, and each single product's mask
needs its own, so the five `common_evaluation_mask.csv` files describe five different populations.
Their intersection is the one population on which "which anchor is better" has an answer. Within it,
every anchor is bootstrapped against GPM - the best single product - with the same method on both
sides.
"""
function anchor_comparison(cfg, grid, anchor_predictions, masks, nearest_distance)
    shared = reduce(.&, (masks[product] for product in ANCHORS))
    metric_rows = NamedTuple[]
    bootstrap_rows = NamedTuple[]
    reference = anchor_predictions["GPM"][REFERENCE]
    for method in vcat(PANEL_METHODS, [REFERENCE])
        for (product_index, product) in enumerate(ANCHORS)
            haskey(anchor_predictions[product], method) || continue
            method == REFERENCE && product != "GPM" && continue
            prediction = anchor_predictions[product][method]
            append_stratified_metrics!(metric_rows, SCHEME,
                method == REFERENCE ? "gauge_only" : product, method, grid.times, grid.Y_obs,
                prediction, BitMatrix(shared), nearest_distance, EVENT_THRESHOLDS)
            (method == REFERENCE || product == "GPM") && continue
            baseline = anchor_predictions["GPM"][method]
            test_mask = BitMatrix(shared .& .!isnan.(baseline) .& .!isnan.(prediction))
            deltas = _daily_bootstrap_delta(MersenneTwister(SEED + 100 * product_index),
                grid.times, grid.Y_obs, baseline, prediction, test_mask, BOOTSTRAP_REPS)
            base_rmse = metric_continuous(grid.Y_obs, baseline; mask=test_mask).RMSE
            method_rmse = metric_continuous(grid.Y_obs, prediction; mask=test_mask).RMSE
            push!(bootstrap_rows, (;
                method, baseline_product="GPM", product, n=count(test_mask),
                RMSE_baseline=base_rmse, RMSE_product=method_rmse,
                relative_improvement=(base_rmse - method_rmse) / base_rmse,
                relative_ci_low=quantile(deltas, 0.025) / base_rmse,
                relative_ci_high=quantile(deltas, 0.975) / base_rmse,
                pvalue=min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0))),
            ))
        end
    end
    comparisons = DataFrame(bootstrap_rows)
    comparisons.pvalue_holm = _holm_adjust(comparisons.pvalue)
    metrics = filter(row -> row.n > 0 && row.group in
        ("overall", "rain_intensity", "event_threshold"), DataFrame(metric_rows))
    select!(metrics, :product, :method, :group, :level, :n, :RMSE, :MAE, :Bias, :r,
        :POD, :FAR, :CSI)
    return metrics, comparisons, count(shared)
end

"""
POD, FAR and CSI of each panel method against `adw`, with day-block bootstrap intervals.

Detection skill is where the satellite is expected to pay: a gauge-only interpolator smooths a
storm's peak away, and `metric_event` counts that as a miss. The intervals come from resampling
whole days' contingency counts, so they carry the same weather-system dependence the RMSE
intervals do.
"""
function detection_skill(grid, product_predictions, masks)
    rows = NamedTuple[]
    for (product_index, product) in enumerate(ANCHORS)
        predictions = product_predictions[product]
        baseline = predictions[REFERENCE]
        methods = vcat(PANEL_METHODS, ["idw", "tps", "gwr"])
        for (method_index, method) in enumerate(methods)
            haskey(predictions, method) || continue
            treatment = predictions[method]
            test_mask = BitMatrix(masks[product] .& .!isnan.(baseline) .& .!isnan.(treatment))
            for (threshold_index, threshold) in enumerate(EVENT_THRESHOLDS)
                result = daily_bootstrap_event_delta(
                    MersenneTwister(SEED + 10_000 * product_index + 100 * method_index +
                                    threshold_index),
                    grid.times, grid.Y_obs, baseline, treatment, test_mask, threshold,
                    BOOTSTRAP_REPS,
                )
                for (metric_index, metric) in enumerate(("POD", "FAR", "CSI"))
                    deltas = filter(isfinite, result.deltas[:, metric_index])
                    push!(rows, (;
                        product, method, threshold, metric,
                        value_adw=result.baseline[metric_index],
                        value_method=result.treatment[metric_index],
                        delta=result.treatment[metric_index] - result.baseline[metric_index],
                        ci_low=isempty(deltas) ? NaN : quantile(deltas, 0.025),
                        ci_high=isempty(deltas) ? NaN : quantile(deltas, 0.975),
                        pvalue=isempty(deltas) ? NaN :
                            min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0))),
                    ))
                end
            end
        end
    end
    return DataFrame(rows)
end

"""Fold-to-fold range of the merged anchors' fusion weights."""
function merged_anchor_coefficients(selection::DataFrame)
    columns = [:intercept, :beta_fy4b, :beta_gpm, :beta_gsmap]
    return combine(groupby(selection, [:scheme, :product, :variant]),
        :fell_back => sum => :folds_fell_back,
        [column => mean => Symbol(column, "_mean") for column in columns]...,
        [column => minimum => Symbol(column, "_min") for column in columns]...,
        [column => maximum => Symbol(column, "_max") for column in columns]...)
end

"""The blending weights each fold chose, one row per (product, method, axis)."""
blend_weights(selection::DataFrame) = combine(
    groupby(filter(:scheme => ==(SCHEME), selection), [:product, :method, :axis]),
    :lambda => mean => :lambda_mean, :lambda => minimum => :lambda_min,
    :lambda => maximum => :lambda_max, :lambdas => (x -> join(x, " / ")) => :fold_lambdas,
    [:inner_RMSE, :unblended_RMSE] =>
        ((a, b) -> mean((b .- a) ./ b)) => :mean_inner_improvement,
)

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    diagnostics_dir = joinpath(ROOT, "output", "benchmark_diagnostics", basename(run_dir))
    comparisons_path = joinpath(diagnostics_dir, "paired_comparisons_by_method.csv")
    isfile(comparisons_path) || error(
        "no $comparisons_path; run scripts/run_claim_reassessment.jl $run_dir first")
    mkpath(OUTDIR)
    write(joinpath(OUTDIR, "source_run.txt"), run_dir * "\n")

    pooled = _pooled(CSV.read(joinpath(run_dir, "metrics_pooled.csv"), DataFrame;
        types=Dict(:level => String)))
    write_csv_atomic(joinpath(OUTDIR, "method_summary.csv"), method_summary(pooled))
    write_csv_atomic(joinpath(OUTDIR, "stratified_summary.csv"), stratified_summary(pooled))
    write_csv_atomic(joinpath(OUTDIR, "fusion_vs_gauge_only.csv"),
        fusion_vs_gauge_only(CSV.read(comparisons_path, DataFrame)))
    write_csv_atomic(joinpath(OUTDIR, "merged_anchor_coefficients.csv"), merged_anchor_coefficients(
        CSV.read(joinpath(run_dir, "fused_anchor_selection.csv"), DataFrame)))
    write_csv_atomic(joinpath(OUTDIR, "blend_weights.csv"), blend_weights(
        CSV.read(joinpath(run_dir, "blend_selection.csv"), DataFrame)))
    scope = CSV.read(joinpath(run_dir, "benchmark_scope.csv"), DataFrame)
    write_csv_atomic(joinpath(OUTDIR, "source_provenance.csv"),
        filter(:key => in(("git_commit", "git_branch", "git_dirty")), scope))
    println("Wrote the tables read off the run's own outputs")

    # The rest re-scores stored predictions, which needs the observations on the run's grid.
    # `load_global_common_product_data` writes its QC table into `outdir`, so it points here.
    mger = study_config(OUTDIR)
    _, ids, product_data = load_global_common_product_data(mger)
    grid = product_data["GPM"]
    cfg = InterpolationBenchmarkConfig(mger=mger, bootstrap_reps=BOOTSTRAP_REPS, seed=20260627)
    scheme_dir = joinpath(run_dir, SCHEME)
    station_meta = load_station_meta(mger.station_meta_path;
        station_id_col=mger.station_id_col, lon_col=mger.lon_col, lat_col=mger.lat_col)
    fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))
    nearest_distance = nearest_train_km([fold_map[id] for id in ids],
        build_X_lonlat(station_meta, ids))

    masks = Dict(product => read_mask_matrix(
        joinpath(scheme_dir, lowercase(product), "common_evaluation_mask.csv"), ids)
        for product in ANCHORS)
    methods = unique(vcat(PANEL_METHODS, INTERPOLATORS))
    predictions = Dict(product => load_methods(joinpath(scheme_dir, lowercase(product)), ids, methods)
        for product in ANCHORS)
    println("Loaded predictions for $(join(ANCHORS, ", "))")

    ranking = interpolation_ranking(cfg, grid, predictions["GPM"], masks["GPM"])
    write_csv_atomic(joinpath(OUTDIR, "interpolation_ranking.csv"), ranking)
    println("Wrote interpolation_ranking.csv")

    anchor_metrics, anchor_tests, shared_cells = anchor_comparison(
        cfg, grid, predictions, masks, nearest_distance)
    write_csv_atomic(joinpath(OUTDIR, "anchor_metrics.csv"), anchor_metrics)
    write_csv_atomic(joinpath(OUTDIR, "anchor_comparison.csv"), anchor_tests)
    println("Wrote the anchor comparison on $shared_cells shared cells")

    write_csv_atomic(joinpath(OUTDIR, "detection_skill_bootstrap.csv"),
        detection_skill(grid, predictions, masks))
    println("Wrote detection_skill_bootstrap.csv to $OUTDIR")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
