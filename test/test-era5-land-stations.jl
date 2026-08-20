using CSV, DataFrames, Dates, JSON, Test

if !isdefined(@__MODULE__, :ERA5LandStations)
    include(joinpath(@__DIR__, "..", "src", "ERA5LandStations.jl"))
end
using .ERA5LandStations

@testset "ERA5-Land station requests" begin
    request = build_era5_land_request(110.312, 31.25)
    @test request["variable"] == ERA5_LAND_VARIABLES
    @test request["location"] == Dict(
        "longitude" => 110.312, "latitude" => 31.25,
    )
    @test request["date"] == ["2022-01-01/2024-12-31"]
    @test request["data_format"] == "csv"
    @test_throws ArgumentError build_era5_land_request(
        110.312, 31.25; start_date=Date(2024), end_date=Date(2022),
    )

    mktempdir() do temp_dir
        station_path = joinpath(temp_dir, "station_meta.csv")
        stations = DataFrame(
            station_id=[lpad(string(index), 8, '0') for index in 1:237],
            lon=fill(110.312, 237),
            lat=fill(31.25, 237),
        )
        CSV.write(station_path, stations)
        validated = read_era5_stations(station_path)
        @test nrow(validated) == 237

        audit_dir = joinpath(temp_dir, "audit")
        raw_dir = joinpath(temp_dir, "raw")
        manifest, manifest_path = write_station_requests(
            validated, audit_dir, raw_dir,
        )
        @test isfile(manifest_path)
        @test nrow(manifest) == 237
        request_path = joinpath(audit_dir, "requests", "00000001.json")
        @test isfile(request_path)
        saved = JSON.parsefile(request_path)
        @test saved["location"]["longitude"] == 110.312

        CSV.write(station_path, first(stations, 236))
        @test_throws ArgumentError read_era5_stations(station_path)
    end
end
