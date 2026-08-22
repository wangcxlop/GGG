#!/usr/bin/env julia

# Re-run the formal GWR claim assessment for every GWR-family method, not only `residual_gwr`.
#
# The benchmark's own `claim_assessment.csv` is produced with the defaults of
# `assess_gwr_claim`, which assess `residual_gwr` and gate coverage on the *shared* evaluation
# mask. That mask is pinned to `MASK_METHODS` and sits at ~0.80 because direct `gwr` fails ~20%
# of cells, so the coverage gate can never pass for anyone, and the methods that actually win
# this run (`hurdle_gwr`, `residual_gwr_const`) are never assessed at all.
#
# This reads the artefacts the run already wrote, so it costs no benchmark re-run.
#
#   julia --project=. scripts/run_claim_reassessment.jl
#   julia --project=. scripts/run_claim_reassessment.jl <benchmark_output_dir>

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using CSV, DataFrames, Dates, Random, Statistics

include(joinpath(ROOT, "src", "InterpolationBenchmark.jl")) # also brings in MGERPipeline.jl
include(joinpath(ROOT, "src", "BenchmarkDiagnostics.jl"))
using .BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output", "interpolation_benchmark_full_joint_covariates_nested_localgrid",
)
# The claim is defined on the balanced-spatial partition only (see `_one_metric`).
const SCHEME = "balanced_spatial"
const ASSESSED_METHODS = [
    "gwr", "gwr_const", "residual_gwr", "residual_gwr_const",
    "mixed_gwr", "mgwr", "hurdle_gwr",
]

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
        expected_common_time_count=11426,
    )
end

"""Load every method's out-of-fold prediction matrix for one (scheme, product) cell."""
function load_predictions(product_dir::AbstractString, ids::Vector{String})
    predictions = Dict{String,Matrix{Float64}}()
    for path in readdir(product_dir; join=true)
        name = basename(path)
        startswith(name, "oof_") && endswith(name, ".csv") || continue
        predictions[name[5:end-4]] = read_wide_matrix(path, ids)
    end
    return predictions
end

"""
Cells the benchmark could score at all: a finite observation *and* a finite satellite value.

This is also what the shared evaluation mask requires, since `raw` is one of the `MASK_METHODS`,
so it is the common denominator for "did this method predict where it could?" — unlike the
`coverage` column in `metrics_pooled.csv`, which divides by the whole station x time grid and
which every method inherits from direct `gwr`'s failures.
"""
evaluable_cells(y_obs::Matrix{Float64}, y_sat::Matrix{Float64}) =
    .!isnan.(y_obs) .& .!isnan.(y_sat)

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(ROOT, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)

    metrics = CSV.read(joinpath(run_dir, "metrics_pooled.csv"), DataFrame)
    mger = study_config(outdir)
    products, ids, product_data = load_global_common_product_data(mger)
    # Only `bootstrap_reps` and `seed` are read out of the config here; both must match the run
    # that wrote the artefacts or the bootstrap draws differ.
    cfg = InterpolationBenchmarkConfig(mger=mger, bootstrap_reps=2000, seed=20260627)

    # Read each product's OOF matrices once. They are ~2.7M values apiece and every assessed
    # method needs the same three traditional baselines, so re-reading per method would dominate
    # the runtime.
    # Product strings stay UPPERCASE: the bootstrap RNG seed hashes the product name
    # (`seed + 1000 * baseline_index + sum(codeunits(product))`), so feeding it the lowercase
    # directory name would silently change every draw.
    predictions = Dict(product => load_predictions(
        joinpath(run_dir, SCHEME, lowercase(product)), ids) for product in products)
    masks = Dict(product => read_mask_matrix(
        joinpath(run_dir, SCHEME, lowercase(product), "common_evaluation_mask.csv"), ids)
        for product in products)
    evaluable = Dict(product => evaluable_cells(
        product_data[product].Y_obs, product_data[product].Y_sat) for product in products)

    bootstrap_rows = NamedTuple[]
    claim_frames = DataFrame[]
    coverage_rows = NamedTuple[]

    for method in ASSESSED_METHODS
        method_products = String[]
        own = Dict{String,Float64}()
        for product in products
            haskey(predictions[product], method) || continue
            data = product_data[product]
            mask = masks[product]
            n_evaluable = count(evaluable[product])
            own[product] = count(
                evaluable[product] .& .!isnan.(predictions[product][method])) / n_evaluable
            push!(method_products, product)
            push!(coverage_rows, (;
                scheme=SCHEME, product, method, n_evaluable, own_coverage=own[product],
                # Same denominator as `own_coverage`, so the two are directly comparable; the
                # gap between them is what the other MASK_METHODS cost this method.
                shared_mask_of_evaluable=count(mask) / n_evaluable,
                shared_mask_of_grid=count(mask) / length(mask),
            ))
            append!(bootstrap_rows, paired_bootstrap_rows(
                cfg, SCHEME, product, data.times, data.Y_obs, predictions[product], mask;
                method, pairwise_mask=true,
            ))
            println("  $method / $product: own coverage " *
                "$(round(own[product]; digits=4)), shared mask " *
                "$(round(count(mask) / n_evaluable; digits=4)) of evaluable")
        end
        isempty(method_products) && continue
        push!(claim_frames, assess_gwr_claim(
            metrics, DataFrame(bootstrap_rows), method_products;
            method, own_coverage=own,
        ))
    end

    # `assess_gwr_claim` reports `supported_product_count` / `overall_claim_supported` per call,
    # so each frame's summary columns already refer to that method's own three products.
    claims = vcat(claim_frames...)
    comparisons = DataFrame(bootstrap_rows)
    CSV.write(joinpath(outdir, "claim_assessment_by_method.csv"), claims)
    CSV.write(joinpath(outdir, "paired_comparisons_by_method.csv"), comparisons)
    CSV.write(joinpath(outdir, "method_coverage.csv"), DataFrame(coverage_rows))
    println("Wrote claim re-assessment to $outdir")
    return (; claims, comparisons)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
