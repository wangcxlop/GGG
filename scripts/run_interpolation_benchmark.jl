#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using Dates

include(joinpath(ROOT, "src", "InterpolationBenchmark.jl"))

function benchmark_config(mode::Symbol; with_random::Bool=false)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be :smoke or :full"))
    smoke = mode == :smoke
    outdir = joinpath(
        ROOT, "output", smoke ? "interpolation_benchmark_smoke_202206" : "interpolation_benchmark_full",
    )
    mger = MGERConfig(
        station_meta_path=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        obs_hourly_wide_path=joinpath(ROOT, "data", "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                ROOT, "data", "processed",
                smoke ? "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
                    "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                ROOT, "data", "processed", "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv",
            ),
            "GSMaP" => joinpath(
                ROOT, "data", "processed", "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv",
            ),
        ),
        outdir=outdir,
        kernels=smoke ? [GAUSSIAN, BISQUARE] :
            [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR],
        bw_adaptive=smoke ? [30.0, 80.0] : [30.0, 50.0, 80.0, 120.0],
        bw_fixed_km=smoke ? [30.0, 50.0] : [10.0, 20.0, 30.0, 50.0],
        rain_threshold=0.1,
        use_loocv_eval=true,
        analysis_start=DateTime(2022, 6, 1, 9),
        # The observation file has no June-September 2025 truth; full validation ends in 2024.
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2024, 10, 1, 8),
        expected_common_time_count=smoke ? nothing : 8067,
    )
    schemes = smoke && !with_random ? [:balanced_spatial] : [:balanced_spatial, :random]
    return InterpolationBenchmarkConfig(
        mger=mger,
        k=5,
        seed=20260627,
        cv_schemes=schemes,
        idw_powers=smoke ? [1.5, 2.0, 2.5] : [1.0, 1.5, 2.0, 2.5, 3.0],
        neighbor_candidates=smoke ? Union{Nothing,Int}[16, 32] :
            Union{Nothing,Int}[8, 16, 32, 64, nothing],
        tps_smooth_candidates=smoke ? [1e-3, 1e-2, 1e-1] :
            [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
        min_tuning_coverage=0.95,
        tuning_max_times=smoke ? 72 : 336,
        event_thresholds=[0.1, 2.5, 8.0, 16.0],
        bootstrap_reps=smoke ? 200 : 2000,
    )
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    with_random = "--with-random" in args
    cfg = benchmark_config(mode; with_random)
    println("Running interpolation benchmark: mode=$mode, output=$(cfg.mger.outdir)")
    result = run_interpolation_benchmark(cfg)
    println("Finished: $(nrow(result.metrics)) metric rows, $(nrow(result.scans)) scan rows")
    if nrow(result.claim) > 0
        println(result.claim)
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

