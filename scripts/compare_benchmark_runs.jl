#!/usr/bin/env julia

# Compare two interpolation-benchmark runs on the cells they both evaluated.
#
# A raw RMSE delta between two run directories is unreadable: the shared evaluation mask keeps
# only cells where every mask-defining method is finite, so any change to one method's coverage
# resizes the denominator for all fourteen, and the cells that move are systematically the hard
# ones. This intersects the two masks and reports a paired per-method difference.
#
#   julia --project=. scripts/compare_benchmark_runs.jl <before_dir> <after_dir> [out_name]
#   julia --project=. scripts/compare_benchmark_runs.jl --smoke <before_dir> <after_dir> [out_name]
#
# `--smoke` rebuilds the June-2022 grid that `run_interpolation_benchmark.jl smoke` scores on,
# rather than the full 2022-2024 one; the two runs must have been produced in the same mode.
#
# Writes `output/benchmark_diagnostics/<out_name>/run_comparison.csv`, defaulting `out_name` to
# `compare_<basename(before)>_vs_<basename(after)>`.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Printf

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("BenchmarkDiagnostics")
using Main.BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

"""
Rebuild the common station/time grid the full benchmark runs used.

Same study-area wiring as `run_benchmark_diagnostics.jl` and `run_mgwr_diagnostics.jl`; it stays
in the script because it hard-codes file locations, which do not belong in `src/`.
"""
function load_common_data(outdir::AbstractString; smoke::Bool=false)
    cfg = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
                    "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv" :
                    "hubei_gpm_hourly_2022_2024_full_aligned.csv",
            ),
            "GSMaP" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv" :
                    "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
            ),
        ),
        outdir=outdir,
        rain_threshold=0.1,
        analysis_start=smoke ? DateTime(2022, 6, 1, 9) : DateTime(2022, 1, 1, 9),
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2025, 1, 1, 8),
        expected_common_time_count=smoke ? nothing : 11426,
    )
    return load_global_common_product_data(cfg)
end

function main(args=ARGS)
    smoke = "--smoke" in args
    positional = filter(!startswith("--"), args)
    length(positional) >= 2 || error(
        "usage: compare_benchmark_runs.jl [--smoke] <before_dir> <after_dir> [out_name]",
    )
    before_dir = abspath(positional[1])
    after_dir = abspath(positional[2])
    for directory in (before_dir, after_dir)
        isdir(directory) || error("benchmark output directory not found: $directory")
    end
    out_name = length(positional) >= 3 ? positional[3] :
        "compare_$(basename(before_dir))_vs_$(basename(after_dir))"
    outdir = joinpath(ROOT, "output", "benchmark_diagnostics", out_name)
    mkpath(outdir)

    mask_moved = String[]
    println("before: $before_dir")
    println("after:  $after_dir")
    products, ids, product_data = load_common_data(outdir; smoke)

    comparison_rows = DataFrame[]
    for scheme in ["balanced_spatial", "random"]
        for product in products
            before_product = joinpath(before_dir, scheme, lowercase(product))
            after_product = joinpath(after_dir, scheme, lowercase(product))
            (isdir(before_product) && isdir(after_product)) || continue
            data = product_data[product]
            before = load_prediction_matrices(before_product, ids)
            after = load_prediction_matrices(after_product, ids)
            before_mask = read_mask_matrix(
                joinpath(before_product, "common_evaluation_mask.csv"), ids,
            )
            after_mask = read_mask_matrix(
                joinpath(after_product, "common_evaluation_mask.csv"), ids,
            )
            table = run_comparison_table(
                data.Y_obs, before, before_mask, after, after_mask; scheme, product,
            )
            push!(comparison_rows, table)
            shared = count(before_mask .& after_mask)
            @printf("  %-17s %-6s mask %d -> %d, shared %d (%.1f%% of the wider)\n",
                scheme, product, count(before_mask), count(after_mask), shared,
                100 * shared / max(count(before_mask), count(after_mask)))
            count(before_mask) == count(after_mask) || push!(
                mask_moved,
                "$scheme/$product $(count(before_mask))->$(count(after_mask))",
            )
        end
    end
    isempty(comparison_rows) && error("the two runs share no scheme/product directory")
    # The mask is an intersection over `MASK_METHODS`, so a change to any one of those
    # methods moves the denominator *every* method is scored on, the traditional
    # baselines included. When that happens the unpaired columns of
    # `run_comparison.csv` are not comparable between the two runs - only the
    # `*_paired` ones, which are restricted to cells both runs evaluated. Raised as a
    # warning rather than left as two numbers on a line, because reading a mask change
    # as a method result is the mistake this table exists to prevent.
    isempty(mask_moved) || @warn(
        "the evaluation mask changed between these runs; compare the *_paired columns " *
        "only, since RMSE_before and RMSE_after are computed over different cell sets",
        cells=join(mask_moved, ", "),
    )

    comparison = vcat(comparison_rows...)
    CSV.write(joinpath(outdir, "run_comparison.csv"), comparison)
    println("Wrote $(joinpath(outdir, "run_comparison.csv"))")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
