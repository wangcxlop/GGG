#!/usr/bin/env julia

# Re-run the formal GWR claim assessment for every GWR-family method, not only `residual_gwr`.
#
# The benchmark's own `claim_assessment.csv` is produced with the defaults of
# `assess_gwr_claim`, which assess `residual_gwr` and gate coverage on the *shared* evaluation
# mask. That mask is pinned to `MASK_METHODS` and sits at ~0.80 because direct `gwr` fails ~20%
# of cells, so the coverage gate can never pass for anyone, and every GWR-family method other
# than `residual_gwr` is never assessed at all.
#
# This reads the artefacts the run already wrote, so it costs no benchmark re-run.
#
#   julia --project=. scripts/run_claim_reassessment.jl
#   julia --project=. scripts/run_claim_reassessment.jl <benchmark_output_dir>

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Random, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("BenchmarkDiagnostics")
using Main.BenchmarkDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_localgrid_mgwrintercept_only",
)
# The claim is defined on the balanced-spatial partition only (see `_one_metric`).
const SCHEME = "balanced_spatial"

# `hurdle_gwr` is no longer produced by the benchmark, but it is kept here so older run
# directories that still contain it can be re-assessed; methods absent from disk are skipped.
const ASSESSED_METHODS = [
    "gwr", "residual_gwr", "mixed_gwr", "mgwr", "hurdle_gwr",
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

"""
Locate a run's scheme directory under either output layout.

A single-partition run keeps `<run>/<scheme>/`; a repeated run nests under
`<run>/repeat_NN/<scheme>/` (`_repeat_dir`). Only the first partition writes `oof_*.csv`
(`run_interpolation_benchmark` guards them on `repeat_index == 1` because they are large), so a
re-assessment is inherently a first-partition exercise and says so rather than searching further.
"""
function scheme_directory(run_dir::AbstractString)
    flat = joinpath(run_dir, SCHEME)
    isdir(flat) && return flat
    nested = joinpath(run_dir, "repeat_01", SCHEME)
    isdir(nested) && return nested
    error("no $SCHEME directory in $run_dir (looked for $flat and $nested)")
end

"""
Per-station distance to the nearest station outside its own fold.

`append_stratified_metrics!` needs this for the `nearest_train_km` stratum. The benchmark computes
it in the fold loop and does not write it out, so it is rebuilt here from `split_common.csv`
rather than passing `NaN` and leaving a whole stratum empty in the rebuilt table.
"""
function nearest_train_km(fold_of::Vector{Int}, lonlat::Matrix{Float64})
    distance = fill(NaN, length(fold_of))
    for fold in sort(unique(fold_of))
        validation = findall(==(fold), fold_of)
        training = findall(!=(fold), fold_of)
        (isempty(validation) || isempty(training)) && continue
        distance[validation] = vec(minimum(
            haversine_distance_matrix(lonlat[training, :], lonlat[validation, :]), dims=1,
        ))
    end
    return distance
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(ROOT, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)

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
    scheme_dir = scheme_directory(run_dir)
    predictions = Dict(product => load_predictions(
        joinpath(scheme_dir, lowercase(product)), ids) for product in products)
    masks = Dict(product => read_mask_matrix(
        joinpath(scheme_dir, lowercase(product), "common_evaluation_mask.csv"), ids)
        for product in products)
    # Rebuilt for the `nearest_train_km` stratum of the metrics table assembled below.
    station_meta = load_station_meta(mger.station_meta_path;
        station_id_col=mger.station_id_col, lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))
    fold_of = [fold_map[id] for id in ids]
    nearest_distance = nearest_train_km(fold_of, lonlat)
    evaluable = Dict(product => evaluable_cells(
        product_data[product].Y_obs, product_data[product].Y_sat) for product in products)

    bootstrap_rows = NamedTuple[]
    claim_frames = DataFrame[]
    coverage_rows = NamedTuple[]

    for method in ASSESSED_METHODS
        method_products = String[]
        own = Dict{String,Float64}()
        # Rebuilt per method rather than read from `metrics_pooled.csv`, and kept separate from
        # every other method's rows. Both matter:
        #
        # 1. Estimand. The bootstrap gate is computed on a pairwise mask while the pooled table was
        #    scored on the shared `MASK_METHODS` mask, so a claim row used to mix two different
        #    sample masks - the significance test answering one question and the heavy/moderate/
        #    year/event gates beside it answering another. Everything now runs on `claim_mask`.
        # 2. Isolation. The accumulated `bootstrap_rows` used to be handed to every call, and only
        #    `assess_gwr_claim`'s internal `row.method == method` filter kept the earlier methods'
        #    rows out of the completeness check. Passing this method's rows alone removes the
        #    dependence on that filter holding.
        method_metric_rows = NamedTuple[]
        method_bootstrap_rows = NamedTuple[]
        for product in products
            haskey(predictions[product], method) || continue
            data = product_data[product]
            mask = masks[product]
            n_evaluable = count(evaluable[product])
            own[product] = count(
                evaluable[product] .& .!isnan.(predictions[product][method])) / n_evaluable
            push!(method_products, product)
            # One mask for the whole claim row: evaluable cells this method and all three
            # baselines predicted. Every gate in the row is then scored on the same cells, and
            # `pairwise_mask=true` below becomes a no-op restriction rather than a second mask.
            claim_mask = evaluable[product] .& .!isnan.(predictions[product][method])
            for baseline in TRADITIONAL_METHODS
                claim_mask = claim_mask .& .!isnan.(predictions[product][baseline])
            end
            claim_mask = BitMatrix(claim_mask)
            push!(coverage_rows, (;
                scheme=SCHEME, product, method, n_evaluable, own_coverage=own[product],
                # Same denominator as `own_coverage`, so the two are directly comparable; the
                # gap between them is what the other MASK_METHODS cost this method.
                shared_mask_of_evaluable=count(mask) / n_evaluable,
                shared_mask_of_grid=count(mask) / length(mask),
                claim_mask_of_evaluable=count(claim_mask) / n_evaluable,
            ))
            for scored_method in vcat(TRADITIONAL_METHODS, [method])
                append_stratified_metrics!(
                    method_metric_rows, SCHEME, product, scored_method, data.times, data.Y_obs,
                    predictions[product][scored_method], claim_mask, nearest_distance,
                    cfg.event_thresholds; repeat=1, seed=cfg.seed,
                )
            end
            append!(method_bootstrap_rows, paired_bootstrap_rows(
                cfg, SCHEME, product, data.times, data.Y_obs, predictions[product], claim_mask;
                method, pairwise_mask=true,
            ))
            println("  $method / $product: own coverage " *
                "$(round(own[product]; digits=4)), shared mask " *
                "$(round(count(mask) / n_evaluable; digits=4)) of evaluable, claim mask " *
                "$(round(count(claim_mask) / n_evaluable; digits=4))")
        end
        isempty(method_products) && continue
        push!(claim_frames, assess_gwr_claim(
            DataFrame(method_metric_rows), DataFrame(method_bootstrap_rows), method_products;
            method, own_coverage=own,
        ))
        append!(bootstrap_rows, method_bootstrap_rows)
    end

    # `assess_gwr_claim` reports `supported_product_count` / `overall_claim_supported` per call,
    # so each frame's summary columns already refer to that method's own three products.
    #
    # What it does NOT do is correct for testing five methods. Holm runs within a stratum across
    # the three baselines only (`paired_bootstrap_rows`), so reading down this table and reporting
    # whichever GWR variant came out supported is a post-hoc maximum over five correlated tests on
    # the held-out fold. The benchmark's `auto` method exists to answer "which model should be
    # used" without that; this table answers "how did each one do", which is a different question.
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
