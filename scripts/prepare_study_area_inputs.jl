#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "StudyArea.jl"))

using .StudyArea

const PROCESSED = joinpath(ROOT, "data", "processed")
const STUDY_DATA = joinpath(PROCESSED, "study_area")

function main()
    result = prepare_study_area_inputs(
        station_meta_path=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        wide_sources=Dict(
            "observation" => joinpath(ROOT, "data", "hubei_obs_hourly_2022_2025_JunSep.csv"),
            "FY4B_smoke" => joinpath(PROCESSED, "hubei_fy4b_hourly_202206_strict_navcorrected.csv"),
            "FY4B_full" => joinpath(PROCESSED, "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"),
            "GPM" => joinpath(PROCESSED, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
            "GSMaP" => joinpath(PROCESSED, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
            "FY4B_full_year" => joinpath(PROCESSED, "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv"),
            "GPM_full_year" => joinpath(PROCESSED, "hubei_gpm_hourly_2022_2024_full_aligned.csv"),
            "GSMaP_full_year" => joinpath(PROCESSED, "hubei_gsmap_hourly_2022_2024_full_aligned.csv"),
        ),
        output_dir=STUDY_DATA,
        differences_path=joinpath(
            ROOT, "output", "input_audit", "study_area_source_station_differences.csv",
        ),
    )
    println(
        "Study area $(STUDY_BOUNDS.west)-$(STUDY_BOUNDS.east) E, " *
        "$(STUDY_BOUNDS.south)-$(STUDY_BOUNDS.north) N: " *
        "$(result.station_count)/$(result.raw_station_count) stations",
    )
    println(
        "Actual station bounds: lon $(result.actual_lon_min)-$(result.actual_lon_max), " *
        "lat $(result.actual_lat_min)-$(result.actual_lat_max)",
    )
    foreach(println, values(result.results))
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
