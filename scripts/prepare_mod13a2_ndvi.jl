#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "MOD13A2NDVIProcessing.jl"))

using .MOD13A2NDVIProcessing

function main()
    result = prepare_mod13a2_ndvi(
        joinpath(
            ROOT, "data", "raw", "ndvi", "mod13a2_2022_2024",
            "MOD13A2-NDVI-Hubei-2022-2024-MOD13A2-061-results.csv",
        ),
        joinpath(
            ROOT, "data", "processed", "covariates",
            "station_ndvi_16day_2022_2024.csv",
        ),
        joinpath(
            ROOT, "output", "validation", "mod13a2_ndvi_station_qc.csv",
        ),
    )
    table = result.table
    quality_valid = count(!ismissing, table.ndvi_qc)
    land_valid = count(!ismissing, table.ndvi_land_qc)
    review_stations = count(result.audit.review_required)
    println("Prepared MOD13A2 NDVI: $(size(table, 1)) rows, " *
            "$(length(unique(table.station_id))) stations, " *
            "$(length(result.period_dates)) composite periods")
    println("Quality-valid NDVI: $quality_valid ($(round(100quality_valid / size(table, 1); digits=2))%)")
    println("Land-and-quality-valid NDVI: $land_valid ($(round(100land_valid / size(table, 1); digits=2))%)")
    println("Stations requiring review: $review_stations")
    println("Processed table: $(result.output_path)")
    println("Station audit: $(result.audit_path)")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
