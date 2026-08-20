#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using Dates

include(joinpath(ROOT, "src", "MGERPipeline.jl"))

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

function main(args=ARGS)
    nested_kernel = "--nested-kernel" in args
    # Nested selection gets its own directory so it never overwrites the canonical paired
    # per-kernel comparison output.
    config = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA,
                "hubei_fy4b_hourly_202206_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                STUDY_DATA,
                "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv",
            ),
            "GSMaP" => joinpath(
                STUDY_DATA,
                "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv",
            ),
        ),
        outdir=joinpath(
            ROOT, "output",
            nested_kernel ? "mger_nested_kernel_202206" : "mger_smoke_202206_5kernels_5fold",
        ),
        kernels=[
            MixedGWR.GAUSSIAN,
            MixedGWR.EXPONENTIAL,
            MixedGWR.BISQUARE,
            MixedGWR.TRICUBE,
            MixedGWR.BOXCAR,
        ],
        bw_adaptive=[30.0, 50.0, 80.0, 120.0],
        bw_fixed_km=[10.0, 20.0, 30.0, 50.0],
        rain_threshold=0.1,
        use_loocv_eval=true,
        analysis_start=DateTime(2022, 6, 1, 9),
        analysis_end=DateTime(2022, 7, 1, 8),
        expected_common_time_count=nothing,
    )
    pipeline = nested_kernel ? run_nested_kernel_spatial_kfold_pipeline :
        run_multikernel_spatial_kfold_pipeline
    return pipeline(
        config;
        k=5,
        seed=20260627,
        fold_scheme=:random,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
