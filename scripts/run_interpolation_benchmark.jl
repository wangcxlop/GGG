#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using DataFrames, Dates

include(joinpath(ROOT, "src", "InterpolationBenchmark.jl"))

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

"""Deterministic seed list for repeated cross-validation, so the sweep itself is reproducible."""
repeat_seeds(repeats::Int) = repeats <= 1 ? Int[] : [20260627 + 1000 * i for i in 0:(repeats - 1)]

function benchmark_config(
    mode::Symbol; with_random::Bool=false, legacy_dem::Bool=false, repeats::Int=1,
    nested_covariates::Bool=false,
)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be :smoke or :full"))
    nested_covariates && legacy_dem && throw(ArgumentError(
        "--nested-covariates only applies to the joint-covariate path, not --legacy-dem",
    ))
    smoke = mode == :smoke
    # Repeated runs, and nested-selection runs, get their own directory so they never overwrite
    # the canonical single-partition/fixed-selection outputs that the verify scripts target.
    outdir = joinpath(
        ROOT, "output", (legacy_dem ?
            (smoke ? "interpolation_benchmark_smoke_dem" : "interpolation_benchmark_full_dem") :
            (smoke ? "interpolation_benchmark_smoke_joint_covariates" :
                "interpolation_benchmark_full_joint_covariates")) *
            (nested_covariates ? "_nested" : "") *
            (repeats > 1 ? "_repeats$(repeats)" : ""),
    )
    mkpath(outdir)
    mger = MGERConfig(
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
        kernels=smoke ? [GAUSSIAN, BISQUARE] :
            [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR],
        bw_adaptive=smoke ? [30.0, 80.0] : [30.0, 50.0, 80.0, 120.0],
        bw_fixed_km=smoke ? [30.0, 50.0] : [10.0, 20.0, 30.0, 50.0],
        rain_threshold=0.1,
        use_loocv_eval=true,
        # Full mode spans Jan 2022 - Dec 2024. FY4B's strict-completeness QC yields zero
        # usable hours in every November/December on record (2022-2024), so this window
        # is a no-op vs. a narrower Nov-ending window - stated explicitly rather than
        # relying on that gap.
        analysis_start=smoke ? DateTime(2022, 6, 1, 9) : DateTime(2022, 1, 1, 9),
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2025, 1, 1, 8),
        # FY4B's off-season NC files rarely form a strict-complete hour (near-zero in
        # Nov/Dec and Jan-Mar), so the common-time count grows from 8067 (Jun-Sep only)
        # to 11426, not to a full 3-year hourly count.
        expected_common_time_count=smoke ? nothing : 11426,
    )
    schemes = smoke && !with_random ? [:balanced_spatial] : [:balanced_spatial, :random]
    dem = legacy_dem ? DEMTerrainExperiment.DEMExperimentConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours=smoke ? 20 : 100,
        k=5,
        seed=20260815,
        bandwidth_candidates=[30, 50, 80, 120, 160],
        screen_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
        vif_threshold=5.0,
        tolerance=1e-5,
        max_iterations=200,
    ) : nothing
    joint = legacy_dem ? nothing : JointCovariateBenchmarkConfig(
        spec_path=nested_covariates ? nothing : joinpath(
            ROOT, "output", "joint_variable_selection", "full",
            "joint_final_full_data_spec.csv",
        ),
        terrain_path=joinpath(
            ROOT, "data", "processed", "covariates", "station_terrain.csv",
        ),
        era5_paths=Dict(year_value => joinpath(
            ROOT, "data", "processed", "covariates", "era5_land",
            "era5_land_station_hourly_utc_$(year_value).csv",
        ) for year_value in 2022:2024),
        ndvi_path=joinpath(
            ROOT, "data", "processed", "covariates",
            "station_ndvi_16day_2022_2024.csv",
        ),
        bandwidth_candidates=[30, 50, 80, 120, 160],
        feature_time_offset_hours=9,
        max_ndvi_age_days=32,
        wet_threshold=0.1,
        tolerance=1e-5,
        max_iterations=200,
    )
    # Same permutation rigor as the per-fold nested selection is asked to match the legacy DEM
    # path's precedent (999 full / 99 smoke), rather than reducing it for speed.
    joint_selection = nested_covariates ? JointSelectionConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours=smoke ? 20 : 100,
        k=5,
        bandwidth_candidates=[30, 50, 80, 120, 160],
        independent_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
        vif_threshold=5.0,
    ) : nothing
    # The joint path's significance test was disabled only because fixed full-data selection
    # made it unreliable; nested selection removes that leak, so it gets the same bootstrap
    # budget as the legacy DEM path.
    joint_bootstrap_reps = nested_covariates ? (smoke ? 200 : 2000) : 0
    return InterpolationBenchmarkConfig(
        mger=mger,
        terrain_path=legacy_dem ? joinpath(
            ROOT, "data", "processed", "covariates", "station_terrain.csv",
        ) : nothing,
        dem=dem,
        joint_covariates=joint,
        joint_selection=joint_selection,
        k=5,
        seed=20260627,
        seeds=repeat_seeds(repeats),
        cv_schemes=schemes,
        idw_powers=smoke ? [1.5, 2.0, 2.5] : [1.0, 1.5, 2.0, 2.5, 3.0],
        neighbor_candidates=smoke ? Union{Nothing,Int}[16, 32] :
            Union{Nothing,Int}[8, 16, 32, 64, nothing],
        tps_smooth_candidates=smoke ? [1e-3, 1e-2, 1e-1] :
            [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
        min_tuning_coverage=0.95,
        tuning_max_times=smoke ? 72 : 336,
        event_thresholds=[0.1, 2.5, 8.0, 16.0],
        bootstrap_reps=legacy_dem ? (smoke ? 200 : 2000) : joint_bootstrap_reps,
    )
end

"""Parse `--repeats N`; absent means a single partition (the historical behaviour)."""
function parse_repeats(args)
    index = findfirst(startswith("--repeats"), args)
    index === nothing && return 1
    text = args[index] == "--repeats" ? get(args, index + 1, "1") :
        split(args[index], "=", limit=2)[2]
    repeats = tryparse(Int, text)
    (repeats === nothing || repeats < 1) &&
        throw(ArgumentError("--repeats expects a positive integer, got $text"))
    return repeats
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    with_random = "--with-random" in args
    legacy_dem = "--legacy-dem" in args
    # :full defaults to the leak-free nested per-fold covariate selection (closes the
    # full-data-reused-across-every-fold leak); --no-nested-covariates reproduces the old
    # fixed full-data comparison path, and --legacy-dem never defaults into nested selection.
    nested_covariates = if "--nested-covariates" in args
        true
    elseif "--no-nested-covariates" in args
        false
    else
        mode == :full && !legacy_dem
    end
    repeats = parse_repeats(args)
    cfg = benchmark_config(mode; with_random, legacy_dem, repeats, nested_covariates)
    !legacy_dem && Threads.nthreads() == 1 && @warn(
        "Joint dynamic models are compute intensive; use julia -t auto for parallel hourly fits",
    )
    nested_covariates && @warn(
        "Nested per-fold covariate selection reruns variable selection for every " *
        "(scheme, product, fold, repeat) cell — expect substantially longer runtime than " *
        "the fixed full-data selection path.",
    )
    println("Running interpolation benchmark: mode=$mode, repeats=$repeats, " *
        "nested_covariates=$nested_covariates, output=$(cfg.mger.outdir)")
    result = run_interpolation_benchmark(cfg)
    println("Finished: $(nrow(result.metrics)) metric rows, $(nrow(result.scans)) scan rows")
    if nrow(result.claim) > 0
        println(result.claim)
    end
    if repeats > 1
        println("Across-partition spread: metrics_repeat_summary.csv, method_rank_stability.csv")
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
