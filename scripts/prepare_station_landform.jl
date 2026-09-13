#!/usr/bin/env julia

# Landform regions for the gauge network, from Copernicus GLO-30 local relief.
#
#   julia --project=. scripts/prepare_station_landform.jl
#
# Picks the relief window by the mean change-point method, classifies each gauge, and writes
# data/processed/covariates/station_landform.csv plus the diagnostics that justify the window to
# output/satellite_temporal_evaluation/landform/. Needs GDAL and GMT on PATH, as
# prepare_dem_station_features.jl does.

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("LandformClassification")

using Main.LandformClassification

function main()
    result = prepare_station_landform(;
        dem_path=joinpath(ROOT, "data", "processed", "dem", "copernicus_glo30_utm49n_30m.tif"),
        slope_path=joinpath(ROOT, "data", "processed", "dem", "copernicus_glo30_slope_deg.tif"),
        terrain_path=joinpath(ROOT, "data", "processed", "covariates", "station_terrain.csv"),
        output_path=joinpath(ROOT, "data", "processed", "covariates", "station_landform.csv"),
        diagnostics_dir=joinpath(ROOT, "output", "satellite_temporal_evaluation", "landform"),
    )
    println("\nRelief window (mean change point): $(result.window_side_km) km")
    show(result.counts; allrows=true, allcols=true)
    stations = result.stations
    println("\n\nGauges whose region changes at half or double the window: $(count(.!stations.class_stable))")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
