#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using Dates

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

function main(args=ARGS)
    nested_kernel = "--nested-kernel" in args
    # Two full windows exist in this project and they are not interchangeable. The default is the
    # Jun-Sep window (8067 common hours) that `scripts/audit_mger_inputs.jl` and every
    # variable-selection script already treat as "full". `--full-year` switches to the Jan
    # 2022 - Dec 2024 window (11426 common hours) that `run_interpolation_benchmark.jl full`
    # uses, so the two pipelines can be read on the same hours when that comparison is wanted.
    full_year = "--full-year" in args
    # Nested selection gets its own directory so it never overwrites the canonical paired
    # per-kernel comparison output.
    config = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA,
                full_year ? "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv" :
                    "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                STUDY_DATA,
                full_year ? "hubei_gpm_hourly_2022_2024_full_aligned.csv" :
                    "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv",
            ),
            "GSMaP" => joinpath(
                STUDY_DATA,
                full_year ? "hubei_gsmap_hourly_2022_2024_full_aligned.csv" :
                    "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv",
            ),
        ),
        outdir=joinpath(
            ROOT, "output",
            (nested_kernel ? "mger_nested_kernel_full" : "mger_full_5kernels_5fold") *
                (full_year ? "_fullyear" : ""),
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
        analysis_start=full_year ? DateTime(2022, 1, 1, 9) : DateTime(2022, 6, 1, 9),
        analysis_end=full_year ? DateTime(2025, 1, 1, 8) : DateTime(2024, 10, 1, 8),
        # The smoke leaves this `nothing`; at full scale a silent drift in the common-time
        # intersection would quietly change every reported metric, so it is pinned. Both counts
        # were measured with `audit_mger_inputs`. A mismatch aborts during data load and points
        # at `global_common_time_qc.csv`.
        expected_common_time_count=full_year ? 11426 : 8067,
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
