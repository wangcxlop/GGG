#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "TerrainFeatures.jl"))

using .TerrainFeatures

function main()
    result = extract_station_terrain(
        joinpath(ROOT, "data", "processed", "study_area", "station_meta.csv"),
        joinpath(ROOT, "data", "raw", "copernicus_glo30"),
        joinpath(ROOT, "data", "processed", "dem"),
        joinpath(ROOT, "data", "processed", "covariates", "station_terrain.csv"),
    )
    table = result.table
    println("Extracted DEM features for $(size(table, 1)) stations")
    println("Elevation range: $(extrema(table.elevation_m)) m")
    println("Slope range: $(extrema(table.slope_deg)) degrees")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
