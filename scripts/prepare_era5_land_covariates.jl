#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "ERA5LandCovariates.jl"))

using .ERA5LandCovariates

function main()
    result = prepare_era5_annual_covariates(
        joinpath(ROOT, "data", "processed", "era5_land", "stations_2022_2024"),
        joinpath(ROOT, "data", "processed", "covariates", "era5_land"),
        joinpath(ROOT, "data", "processed", "study_area", "station_meta.csv"),
        joinpath(ROOT, "output", "validation", "era5_land_annual_covariates.csv"),
    )
    println(result.audit)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
