#!/usr/bin/env julia

# Post-hoc diagnosis of an existing interpolation-benchmark run: why do the traditional
# interpolators beat every GWR variant? Reads the artefacts the benchmark already wrote,
# so this needs no (expensive) benchmark re-run.
#
#   julia --project=. scripts/run_benchmark_diagnostics.jl
#   julia --project=. scripts/run_benchmark_diagnostics.jl <benchmark_output_dir>

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "MGERPipeline.jl"))
include(joinpath(ROOT, "src", "BenchmarkDiagnostics.jl"))
using .BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_localgrid_mgwrintercept_only",
)
const NULLS = ["zero", "train_clim", "hour_field_mean"]

"""Rebuild the exact common station/time grid the full benchmark run used."""
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
    # Diagnostics always go to their own subdirectory, so diagnosing a second
    # run never overwrites the first one's numbers.
    outdir = joinpath(ROOT, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)

    println("Reading benchmark run: $run_dir")
    products, ids, product_data = load_common_data(outdir)
    println("Common grid: $(length(ids)) stations x " *
        "$(length(product_data[first(products)].times)) hours")

    null_rows = DataFrame[]
    offset_rows = DataFrame[]
    coefficient_rows = DataFrame[]
    decomposition_rows = DataFrame[]

    for scheme in ["balanced_spatial", "random"]
        scheme_dir = joinpath(run_dir, scheme)
        isdir(scheme_dir) || continue
        fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))
        nulls = null_baseline_matrices(product_data[first(products)].Y_obs, ids, fold_map)

        for product in products
            product_dir = joinpath(scheme_dir, lowercase(product))
            isdir(product_dir) || continue
            data = product_data[product]
            mask = read_mask_matrix(
                joinpath(product_dir, "common_evaluation_mask.csv"), ids,
            )
            predictions = load_prediction_matrices(product_dir, ids)
            merge!(predictions, nulls)

            rescaled, coefficients = satellite_rescale_matrix(
                data.Y_obs, data.Y_sat, ids, fold_map,
            )
            predictions["sat_rescaled"] = rescaled
            coefficients.scheme .= scheme
            coefficients.product .= product
            push!(coefficient_rows, coefficients)

            println("  $scheme / $product: $(length(predictions)) prediction fields, " *
                "mask coverage $(round(count(mask) / length(mask); digits=4))")

            push!(null_rows, null_baseline_table(
                data.Y_obs, predictions, mask; scheme, product, nulls=NULLS,
            ))
            push!(offset_rows, satellite_offset_table(
                data.Y_obs, data.Y_sat, predictions, mask; scheme, product,
                methods=["residual_gwr", "mixed_gwr", "mgwr", "sat_rescaled"],
            ))
            push!(decomposition_rows, mse_decomposition_table(
                data.Y_obs, predictions, mask; scheme, product,
            ))
        end
    end

    CSV.write(joinpath(outdir, "null_baselines.csv"), vcat(null_rows...))
    CSV.write(joinpath(outdir, "satellite_offset.csv"), vcat(offset_rows...))
    CSV.write(joinpath(outdir, "satellite_rescale_coefficients.csv"), vcat(coefficient_rows...))
    CSV.write(joinpath(outdir, "mse_decomposition.csv"), vcat(decomposition_rows...))

    # `joint_bandwidths.csv` is a subset of `parameter_scan.csv` with the same schema, so the
    # scan file alone already covers the joint models - concatenating them would double every row.
    scan = CSV.read(joinpath(run_dir, "parameter_scan.csv"), DataFrame)
    CSV.write(joinpath(outdir, "bandwidth_saturation.csv"), bandwidth_saturation_table(scan))

    CSV.write(joinpath(outdir, "covariate_contribution.csv"), covariate_contribution_table(
        CSV.read(joinpath(run_dir, "covariate_model_status.csv"), DataFrame),
        CSV.read(joinpath(run_dir, "metrics_folds.csv"), DataFrame),
    ))

    println("Wrote diagnostics to $outdir")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
