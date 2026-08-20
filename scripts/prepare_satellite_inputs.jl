#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MGERDataPrep.jl"))

using .MGERDataPrep

const PROCESSED = joinpath(ROOT, "data", "processed")
const AUDIT = joinpath(ROOT, "output", "input_audit")

function main()
    mkpath(PROCESSED)
    mkpath(AUDIT)
    results = [
        prepare_satellite_wide(
            "gpm",
            joinpath(ROOT, "data", "GEE_GPM_HUBEI"),
            joinpath(PROCESSED, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv");
            value_column=:gpm_mm_h,
            duplicate_qc_file=joinpath(AUDIT, "gpm_duplicate_station_times.csv"),
        ),
        prepare_satellite_wide(
            "gsmap",
            joinpath(ROOT, "data", "GEE_GSMAP_HUBEI"),
            joinpath(PROCESSED, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv");
            value_column=:gsmap_mm_h,
            duplicate_qc_file=joinpath(AUDIT, "gsmap_duplicate_station_times.csv"),
        ),
        prepare_satellite_wide(
            "gpm",
            joinpath(ROOT, "data", "GEE_GPM_HUBEI"),
            joinpath(PROCESSED, "hubei_gpm_hourly_2022_2024_full_aligned.csv");
            value_column=:gpm_mm_h,
            years=2022:2024,
            months=1:12,
            duplicate_qc_file=joinpath(AUDIT, "gpm_full_duplicate_station_times.csv"),
        ),
        prepare_satellite_wide(
            "gsmap",
            joinpath(ROOT, "data", "GEE_GSMAP_HUBEI"),
            joinpath(PROCESSED, "hubei_gsmap_hourly_2022_2024_full_aligned.csv");
            value_column=:gsmap_mm_h,
            years=2022:2024,
            months=1:12,
            duplicate_qc_file=joinpath(AUDIT, "gsmap_full_duplicate_station_times.csv"),
        ),
    ]
    foreach(println, results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
