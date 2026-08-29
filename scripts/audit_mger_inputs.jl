#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("MGERDataPrep")

using Main.MGERDataPrep
using Dates

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

function main()
    smoke = "--smoke-202206" in ARGS
    fy4b_name = smoke ?
        "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
        "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"
    result = audit_mger_inputs(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        observation_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        satellite_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA,
                fy4b_name,
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
        output_dir=joinpath(ROOT, "output", smoke ? "input_audit_smoke_202206" : "input_audit"),
        analysis_start=DateTime(2022, 6, 1, 9),
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2024, 10, 1, 8),
        raw_station_meta_path=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        study_bounds=(109.4, 111.6, 31.2, 33.4),
    )
    println(result.summary)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
