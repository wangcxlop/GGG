#!/usr/bin/env julia

# MGWR-specific post-hoc diagnostics for an existing interpolation-benchmark run.
#
# Answers two questions that the general diagnostics in `run_benchmark_diagnostics.jl` do not:
#   1. Which cells does MGWR fail to predict, and do the failures arrive in whole hours
#      (a fit-level failure such as back-fitting non-convergence) or scattered station-by-station?
#   2. MGWR is one of the methods whose finiteness defines the shared evaluation mask, so its
#      failures shrink the denominator every other method is scored on. What does that cost?
#
# Reads only artefacts the benchmark already wrote, so this needs no benchmark re-run.
#
#   julia --project=. scripts/run_mgwr_diagnostics.jl
#   julia --project=. scripts/run_mgwr_diagnostics.jl <benchmark_output_dir>

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates

include(joinpath(ROOT, "src", "MGERPipeline.jl"))
include(joinpath(ROOT, "src", "BenchmarkDiagnostics.jl"))
using .BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output", "interpolation_benchmark_full_joint_covariates_nested_localgrid",
)
const DROPOUT_METHODS = ["mgwr", "residual_gwr", "mixed_gwr", "gwr"]

"""
Rebuild the exact common station/time grid the full benchmark run used.

Deliberately a copy of the same helper in `run_benchmark_diagnostics.jl` rather than a shared
function: it hard-codes the study-area file wiring, which belongs to a script rather than to
`src/`. The generic parts (matrix/mask/fold readers, `load_prediction_matrices`) are shared.
"""
function load_common_data(outdir::AbstractString)
    cfg = MGERConfig(
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
        expected_common_time_count=11426,
    )
    return load_global_common_product_data(cfg)
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(ROOT, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)

    println("Reading benchmark run: $run_dir")
    products, ids, product_data = load_common_data(outdir)

    dropout_rows = DataFrame[]
    mask_cost_rows = DataFrame[]
    mask_checks = NamedTuple[]

    for scheme in ["balanced_spatial", "random"]
        scheme_dir = joinpath(run_dir, scheme)
        isdir(scheme_dir) || continue
        fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))

        for product in products
            product_dir = joinpath(scheme_dir, lowercase(product))
            isdir(product_dir) || continue
            data = product_data[product]
            predictions = load_prediction_matrices(product_dir, ids)

            # The mask rule is duplicated from `InterpolationBenchmark._common_method_mask`;
            # checking the rebuild against the run's own mask is what stops the copy drifting.
            stored = read_mask_matrix(
                joinpath(product_dir, "common_evaluation_mask.csv"), ids,
            )
            rebuilt = rebuild_common_mask(data.Y_obs, predictions)
            disagreements = count(stored .!= rebuilt)
            push!(mask_checks, (;
                scheme, product, stored_cells=count(stored),
                rebuilt_cells=count(rebuilt), disagreements,
            ))
            disagreements == 0 || @warn(
                "rebuilt common mask differs from the run's own mask — MASK_METHODS may have " *
                "changed since this run was produced; mask_cost numbers are not trustworthy",
                scheme, product, disagreements,
            )

            push!(dropout_rows, dropout_table(
                data.Y_obs, data.Y_sat, predictions, ids, fold_map;
                scheme, product, methods=DROPOUT_METHODS,
            ))
            push!(mask_cost_rows, mask_cost_table(
                data.Y_obs, predictions; scheme, product, excluded="mgwr",
            ))
            println("  $scheme / $product: mask $(count(stored)) cells, " *
                "$(disagreements) rebuild disagreements")
        end
    end

    CSV.write(joinpath(outdir, "mgwr_dropout.csv"), vcat(dropout_rows...))
    CSV.write(joinpath(outdir, "mgwr_mask_cost.csv"), vcat(mask_cost_rows...))
    CSV.write(joinpath(outdir, "mgwr_mask_rebuild_check.csv"), DataFrame(mask_checks))

    println("Wrote MGWR diagnostics to $outdir")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
