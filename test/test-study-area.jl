using Test
using CSV, DataFrames

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("StudyArea")
using Main.StudyArea

function write_study_wide(path, station_ids)
    table = DataFrame(time=["2022-06-01T09:00:00Z", "2022-06-01T10:00:00Z"])
    for (index, station_id) in enumerate(station_ids)
        table[!, Symbol(station_id)] = [index, index + 1]
    end
    CSV.write(path, table)
end

@testset "Study area input filtering" begin
    bounds = StudyBounds(109.4, 111.6, 31.2, 33.4)
    metadata = DataFrame(
        station_id=["west", "east", "south", "north", "inside", "outside"],
        lon=[109.4, 111.6, 110.0, 110.0, 110.5, 111.7],
        lat=[32.0, 32.0, 31.2, 33.4, 32.5, 32.0],
    )
    selected = filter_stations(metadata; bounds)
    @test selected.station_id == ["east", "inside", "north", "south", "west"]
    @test nrow(selected) == 5

    missing_lon = copy(metadata)
    allowmissing!(missing_lon, :lon)
    missing_lon.lon[1] = missing
    @test_throws ArgumentError filter_stations(missing_lon; bounds)

    duplicate = vcat(metadata, metadata[1:1, :])
    @test_throws ArgumentError filter_stations(duplicate; bounds)

    mktempdir() do temp_dir
        metadata_path = joinpath(temp_dir, "stations.csv")
        CSV.write(metadata_path, metadata)
        complete_path = joinpath(temp_dir, "complete.csv")
        write_study_wide(complete_path, ["outside", "west", "north", "inside", "east", "south"])
        output_dir = joinpath(temp_dir, "processed")
        differences_path = joinpath(temp_dir, "differences.csv")

        result = prepare_study_area_inputs(
            station_meta_path=metadata_path,
            wide_sources=Dict("observation" => complete_path),
            output_dir=output_dir,
            differences_path=differences_path,
            bounds=bounds,
        )
        @test result.station_count == 5
        output = CSV.read(joinpath(output_dir, "complete.csv"), DataFrame)
        @test names(output) == ["time", "east", "inside", "north", "south", "west"]
        @test nrow(output) == 2

        incomplete_path = joinpath(temp_dir, "incomplete.csv")
        write_study_wide(incomplete_path, ["west", "north", "inside", "east", "outside"])
        @test_throws ArgumentError prepare_study_area_inputs(
            station_meta_path=metadata_path,
            wide_sources=Dict("observation" => incomplete_path),
            output_dir=joinpath(temp_dir, "failed"),
            differences_path=joinpath(temp_dir, "failed_differences.csv"),
            bounds=bounds,
        )
        differences = CSV.read(joinpath(temp_dir, "failed_differences.csv"), DataFrame)
        @test any(
            (differences.station_id .== "south") .&
            (differences.issue .== "missing_required_station"),
        )
    end
end
