#!/usr/bin/env julia

# How much does MERGED_OLS_LAGNBR - the OLS fusion widened with each product's t-1/t+1 values and its
# nearest-station mean, fitted per agreement-envelope band - add over MERGED_OLS, and where does it sit
# among the other anchors and correction methods?
#
#   julia -t 4 --project=. scripts/run_fused_anchor_lagnbr_report.jl
#   julia -t 4 --project=. scripts/run_fused_anchor_lagnbr_report.jl <benchmark_output_dir>
#
# Reads a finished `run_interpolation_benchmark.jl full --nested-covariates --satellite-wet-blend
# --blend-axis agreement_envelope --fused-anchor --fused-anchor-lagnbr` run; fits nothing. Writes CSVs
# to output/fused_anchor_lagnbr/, which `py -3.13 scripts/plot_fused_anchor_lagnbr.py` draws.
#
# The run scores each product on its own common mask, and MERGED_OLS_LAGNBR's is 60k cells smaller
# than MERGED_OLS's (fold 4's covariate model left four stations unpredicted), so every comparison
# between products here is re-scored from the stored predictions on the cells both - or all - share.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Random, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("BenchmarkDiagnostics", "TableIO")
using Main.BenchmarkDiagnostics
using Main.TableIO: write_csv_atomic

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const RUN_PREFIX = "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_" *
    "satwetblend_blendagrenv_fusedanchor"
const DEFAULT_RUN = joinpath(ROOT, "output", RUN_PREFIX * "_fusedlagnbr")
# The same configuration without the new product: every row it wrote must be unchanged.
const BASELINE_RUN = joinpath(ROOT, "output", RUN_PREFIX)
const OUTDIR = joinpath(ROOT, "output", "fused_anchor_lagnbr")
const SCHEMES = ["balanced_spatial", "random"]
const SEED = 20260922
const BOOTSTRAP_REPS = 2000

const NEW = "MERGED_OLS_LAGNBR"
const OLD = "MERGED_OLS"
const ANCHORS = ["FY4B", "GPM", "GSMaP", "MERGED_MEAN", OLD, NEW]
const REFERENCE = "adw"
const CORRECTION_METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr", "auto",
    "blend_residual_gwr", "blend_mixed_gwr", "blend_mgwr",
    "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr",
]
# Fixed rather than picked from the results, so the anchor comparison is not each anchor's luckiest.
const PANEL_METHODS = ["raw", "residual_gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
const VS_ADW_METHODS = ["residual_gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
const STRATA = [("all", -Inf, Inf), ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
    ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf)]
const EVENT_THRESHOLDS = [0.1, 2.5, 8.0, 16.0]

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

"""
Pooled overall scores for every (scheme, product, method), as the run reported them.

Each product is on its own mask here, so these rank methods *within* a product; the cross-product
comparisons below re-score on shared cells instead.
"""
function method_summary(pooled::DataFrame)
    overall = filter(row -> row.group == "overall", pooled)
    reference(scheme, product, method) = only(filter(row -> row.scheme == scheme &&
        row.product == product && row.method == method, overall).RMSE)
    return sort!(DataFrame([(;
        scheme=String(row.scheme), product=String(row.product), method=String(row.method),
        n=Int(row.n), coverage=row.coverage, RMSE=row.RMSE, MAE=row.MAE, Bias=row.Bias, r=row.r,
        RMSE_improvement_vs_raw=1 - row.RMSE / reference(row.scheme, row.product, "raw"),
        RMSE_improvement_vs_adw=1 - row.RMSE / reference(row.scheme, row.product, REFERENCE),
    ) for row in eachrow(overall)]), [:scheme, :product, :method])
end

"""RMSE and bias of `prediction` on `mask`."""
function scores(y_obs, prediction, mask)
    errors = prediction[mask] .- y_obs[mask]
    return (; n=count(mask), RMSE=sqrt(mean(abs2, errors)), Bias=mean(errors))
end

"""
`treatment` against `baseline` on `mask`: both RMSEs and the day-block bootstrap interval of the
relative improvement. Positive means `treatment` has the lower RMSE.
"""
function paired_row(rng, times, y_obs, baseline, treatment, mask)
    base = scores(y_obs, baseline, mask)
    treat = scores(y_obs, treatment, mask)
    deltas = _daily_bootstrap_delta(rng, times, y_obs, baseline, treatment, mask, BOOTSTRAP_REPS)
    return (; n=base.n, RMSE_baseline=base.RMSE, RMSE_treatment=treat.RMSE,
        Bias_baseline=base.Bias, Bias_treatment=treat.Bias,
        relative_improvement=1 - treat.RMSE / base.RMSE,
        ci_low=quantile(deltas, 0.025) / base.RMSE, ci_high=quantile(deltas, 0.975) / base.RMSE,
        pvalue=min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0))))
end

stratum_mask(y_obs, mask, low, high) = BitMatrix(mask .& (y_obs .>= low) .& (y_obs .< high))

"""A `station x time` prediction matrix per method that the product directory holds."""
load_methods(product_dir, ids, methods) = Dict(
    method => read_wide_matrix(joinpath(product_dir, "oof_$(method).csv"), ids)
    for method in methods if isfile(joinpath(product_dir, "oof_$(method).csv")))

"""
MERGED_OLS_LAGNBR against MERGED_OLS for every correction method and stratum, and both against
`adw`, on the cells both products' common masks share and both methods predicted.
"""
function head_to_head(scheme, grid, scheme_dir, ids)
    masks = Dict(product => read_mask_matrix(
        joinpath(scheme_dir, lowercase(product), "common_evaluation_mask.csv"), ids)
        for product in (OLD, NEW))
    predictions = Dict(product => load_methods(joinpath(scheme_dir, lowercase(product)), ids,
        CORRECTION_METHODS) for product in (OLD, NEW))
    shared = BitMatrix(masks[OLD] .& masks[NEW] .& .!isnan.(grid.Y_obs))
    rows = NamedTuple[]
    adw_rows = NamedTuple[]
    for (method_index, method) in enumerate(CORRECTION_METHODS)
        haskey(predictions[OLD], method) && haskey(predictions[NEW], method) || continue
        old, new = predictions[OLD][method], predictions[NEW][method]
        cells = BitMatrix(shared .& .!isnan.(old) .& .!isnan.(new))
        for (stratum_index, (stratum, low, high)) in enumerate(STRATA)
            mask = stratum_mask(grid.Y_obs, cells, low, high)
            rng = MersenneTwister(SEED + 1000 * method_index + stratum_index)
            push!(rows, merge((; scheme, method, stratum),
                paired_row(rng, grid.times, grid.Y_obs, old, new, mask)))
        end
        method in VS_ADW_METHODS || continue
        for product in (OLD, NEW), (stratum_index, (stratum, low, high)) in enumerate(STRATA)
            reference = predictions[product][REFERENCE]
            mask = stratum_mask(grid.Y_obs,
                BitMatrix(cells .& .!isnan.(reference)), low, high)
            rng = MersenneTwister(SEED + 50_000 + 1000 * method_index + stratum_index)
            push!(adw_rows, merge((; scheme, product, method, stratum),
                paired_row(rng, grid.times, grid.Y_obs, reference, predictions[product][method],
                    mask)))
        end
    end
    comparison = DataFrame(rows)
    comparison.pvalue_holm = _holm_adjust(comparison.pvalue)
    vs_adw = DataFrame(adw_rows)
    vs_adw.pvalue_holm = _holm_adjust(vs_adw.pvalue)
    detection = detection_skill(grid, predictions, shared)
    detection.scheme .= scheme
    return comparison, vs_adw, detection, count(shared)
end

"""POD/FAR/CSI of the panel methods on both merged anchors against `adw`, with day-block CIs."""
function detection_skill(grid, predictions, shared)
    rows = NamedTuple[]
    for (product_index, product) in enumerate((OLD, NEW))
        reference = predictions[product][REFERENCE]
        for (method_index, method) in enumerate(PANEL_METHODS)
            treatment = predictions[product][method]
            mask = BitMatrix(shared .& .!isnan.(reference) .& .!isnan.(treatment))
            for (threshold_index, threshold) in enumerate(EVENT_THRESHOLDS)
                result = daily_bootstrap_event_delta(
                    MersenneTwister(SEED + 10_000 * product_index + 100 * method_index +
                                    threshold_index),
                    grid.times, grid.Y_obs, reference, treatment, mask, threshold, BOOTSTRAP_REPS)
                for (metric_index, metric) in enumerate(("POD", "FAR", "CSI"))
                    deltas = filter(isfinite, result.deltas[:, metric_index])
                    push!(rows, (;
                        product, method, threshold, metric,
                        value_adw=result.baseline[metric_index],
                        value_method=result.treatment[metric_index],
                        delta=result.treatment[metric_index] - result.baseline[metric_index],
                        ci_low=isempty(deltas) ? NaN : quantile(deltas, 0.025),
                        ci_high=isempty(deltas) ? NaN : quantile(deltas, 0.975),
                    ))
                end
            end
        end
    end
    return DataFrame(rows)
end

"""
All six anchors on the cells every one of them was scored on, with the panel methods, each
bootstrapped against MERGED_OLS with the same method on both sides.
"""
function anchor_comparison(grid, scheme, scheme_dir, ids)
    masks = [read_mask_matrix(
        joinpath(scheme_dir, lowercase(product), "common_evaluation_mask.csv"), ids)
        for product in ANCHORS]
    shared = BitMatrix(reduce(.&, masks) .& .!isnan.(grid.Y_obs))
    metric_rows = NamedTuple[]
    test_rows = NamedTuple[]
    for (method_index, method) in enumerate(vcat(PANEL_METHODS, [REFERENCE]))
        predictions = Dict(product => read_wide_matrix(
            joinpath(scheme_dir, lowercase(product), "oof_$(method).csv"), ids)
            for product in ANCHORS)
        cells = BitMatrix(reduce(.&, (.!isnan.(p) for p in values(predictions))) .& shared)
        for (product_index, product) in enumerate(ANCHORS)
            # `adw` never reads the satellite: one row, not six copies of it.
            method == REFERENCE && product != OLD && continue
            for (stratum, low, high) in STRATA
                mask = stratum_mask(grid.Y_obs, cells, low, high)
                push!(metric_rows, merge(
                    (; scheme, product=method == REFERENCE ? "gauge_only" : product, method, stratum),
                    scores(grid.Y_obs, predictions[product], mask)))
            end
            (method == REFERENCE || product == OLD) && continue
            rng = MersenneTwister(SEED + 100 * method_index + product_index +
                                  (scheme == "random" ? 100_000 : 0))
            push!(test_rows, merge((; scheme, product, method, baseline_product=OLD),
                paired_row(rng, grid.times, grid.Y_obs, predictions[OLD], predictions[product],
                    cells)))
        end
    end
    tests = DataFrame(test_rows)
    tests.pvalue_holm = _holm_adjust(tests.pvalue)
    return DataFrame(metric_rows), tests, count(shared)
end

"""Fold mean and range of every band x term coefficient, with the band's cells and fallbacks."""
function coefficient_summary(coefficients::DataFrame)
    table = filter(row -> row.product == NEW, coefficients)
    return combine(groupby(table, [:scheme, :group, :term]),
        :coefficient => mean => :coefficient_mean, :coefficient => minimum => :coefficient_min,
        :coefficient => maximum => :coefficient_max, :n_cell => mean => :n_cell_mean,
        :fell_back => sum => :folds_fell_back, nrow => :folds)
end

"""
How many of the baseline run's rows the new run reproduces byte for byte.

The new product is only appended, so every row of the run without it must be there unchanged;
`improved_product_count` in `claim_assessment.csv` is the one expected exception, because it counts
across products.
"""
function invariance(run_dir)
    rows = NamedTuple[]
    for file in ("metrics_pooled.csv", "paired_comparisons.csv", "claim_assessment.csv")
        old_lines = readlines(joinpath(BASELINE_RUN, file))[2:end]
        new_lines = Set(readlines(joinpath(run_dir, file)))
        push!(rows, (; file, baseline_rows=length(old_lines),
            identical=count(in(new_lines), old_lines),
            changed=count(!in(new_lines), old_lines)))
    end
    return DataFrame(rows)
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    mkpath(OUTDIR)
    write(joinpath(OUTDIR, "source_run.txt"), run_dir * "\n")

    pooled = filter(row -> ismissing(row.fold), CSV.read(
        joinpath(run_dir, "metrics_pooled.csv"), DataFrame; types=Dict(:level => String)))
    write_csv_atomic(joinpath(OUTDIR, "method_summary.csv"), method_summary(pooled))
    write_csv_atomic(joinpath(OUTDIR, "claim_assessment.csv"),
        CSV.read(joinpath(run_dir, "claim_assessment.csv"), DataFrame))
    write_csv_atomic(joinpath(OUTDIR, "grouped_coefficients_summary.csv"), coefficient_summary(
        CSV.read(joinpath(run_dir, "fused_anchor_grouped_coefficients.csv"), DataFrame)))
    status = CSV.read(joinpath(run_dir, "covariate_model_status.csv"), DataFrame)
    write_csv_atomic(joinpath(OUTDIR, "coverage_note.csv"), select(
        filter(row -> row.product in (OLD, NEW, "GPM"), status),
        :scheme, :product, :fold, :method, :status, :prediction_coverage, :covariate_variables))
    isdir(BASELINE_RUN) &&
        write_csv_atomic(joinpath(OUTDIR, "invariance.csv"), invariance(run_dir))
    println("Wrote the tables read off the run's own outputs")

    # The rest re-scores stored predictions, which needs the observations on the run's grid.
    # `load_global_common_product_data` writes its QC table into `outdir`, so it points here.
    _, ids, product_data = load_global_common_product_data(study_config(OUTDIR))
    ids = String.(ids)
    grid = (; times=product_data["GPM"].times, Y_obs=Matrix{Float64}(product_data["GPM"].Y_obs))

    comparisons, vs_adw, detections, cells = DataFrame(), DataFrame(), DataFrame(), NamedTuple[]
    metrics, tests = DataFrame(), DataFrame()
    for scheme in SCHEMES
        scheme_dir = joinpath(run_dir, scheme)
        comparison, adw, detection, shared = head_to_head(scheme, grid, scheme_dir, ids)
        append!(comparisons, comparison)
        append!(vs_adw, adw)
        append!(detections, detection)
        push!(cells, (; scheme, comparison="lagnbr_vs_ols", shared_cells=shared))
        println("$scheme: head to head on $shared shared cells")
        GC.gc()
        scheme_metrics, scheme_tests, all_shared = anchor_comparison(grid, scheme, scheme_dir, ids)
        append!(metrics, scheme_metrics)
        append!(tests, scheme_tests)
        push!(cells, (; scheme, comparison="all_anchors", shared_cells=all_shared))
        println("$scheme: six anchors on $all_shared shared cells")
        GC.gc()
    end
    write_csv_atomic(joinpath(OUTDIR, "lagnbr_vs_ols.csv"), comparisons)
    write_csv_atomic(joinpath(OUTDIR, "vs_adw.csv"), vs_adw)
    write_csv_atomic(joinpath(OUTDIR, "detection_skill.csv"), detections)
    write_csv_atomic(joinpath(OUTDIR, "anchor_metrics.csv"), metrics)
    write_csv_atomic(joinpath(OUTDIR, "anchor_comparison.csv"), tests)
    write_csv_atomic(joinpath(OUTDIR, "shared_cells.csv"), DataFrame(cells))
    println("Wrote the comparisons to $OUTDIR")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
