#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using Dates

include(joinpath(ROOT, "src", "MGERPipeline.jl"))

function main()
    config = MGERConfig(
        station_meta_path=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        obs_hourly_wide_path=joinpath(ROOT, "data", "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                ROOT, "data", "processed",
                "hubei_fy4b_hourly_202206_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                ROOT, "data", "processed",
                "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv",
            ),
            "GSMaP" => joinpath(
                ROOT, "data", "processed",
                "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv",
            ),
        ),
        outdir=joinpath(ROOT, "output", "mger_smoke_202206_5kernels_5fold"),
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
    return run_multikernel_spatial_kfold_pipeline(
        config;
        k=5,
        seed=20260627,
        fold_scheme=:random,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
