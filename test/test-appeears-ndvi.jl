using CSV, DataFrames, Downloads, Test

if !isdefined(@__MODULE__, :AppEEARSNDVI)
    include(joinpath(@__DIR__, "..", "src", "AppEEARSNDVI.jl"))
end

@testset "Downloads response compatibility" begin
    ok = Downloads.Response(nothing, "https://example.test", 200, "OK", Pair{String,String}[])
    failed_response = Downloads.Response(
        nothing, "https://example.test", 503, "Service Unavailable", Pair{String,String}[],
    )
    failed = Downloads.RequestError(
        "https://example.test", 22, "HTTP request returned an error", failed_response,
    )

    @test AppEEARSNDVI._response_status(ok) == 200
    @test AppEEARSNDVI._response_status(failed) == 503
    @test AppEEARSNDVI._retryable_status(0)
    @test AppEEARSNDVI._retryable_status(503)
    @test !AppEEARSNDVI._retryable_status(403)
end
using .AppEEARSNDVI

@testset "AppEEARS MOD13A2 request" begin
    mktempdir() do temp_dir
        station_path = joinpath(temp_dir, "stations.csv")
        CSV.write(station_path, DataFrame(
            station_id=["001", "002"],
            lon=[112.1, 113.2],
            lat=[30.5, 31.6],
        ))

        coordinates = read_station_coordinates(station_path)
        @test length(coordinates) == 2
        @test coordinates[1]["id"] == "001"
        @test coordinates[2]["longitude"] == 113.2
        @test validate_coordinate_bounds(
            coordinates; west=109.4, east=113.2, south=30.5, north=33.4,
        )
        @test_throws ArgumentError validate_coordinate_bounds(
            coordinates; west=109.4, east=111.6, south=31.2, north=33.4,
        )

        task = build_mod13a2_point_task(coordinates)
        @test task["task_type"] == "point"
        @test task["params"]["dates"][1]["startDate"] == "01-01-2022"
        @test task["params"]["dates"][1]["endDate"] == "12-31-2024"
        @test [layer["layer"] for layer in task["params"]["layers"]] == MOD13A2_LAYERS

        request_path = joinpath(temp_dir, "request.json")
        write_request_json(request_path, task)
        @test isfile(request_path)
        @test occursin("MOD13A2.061", read(request_path, String))

        duplicate_path = joinpath(temp_dir, "duplicate.csv")
        CSV.write(duplicate_path, DataFrame(
            station_id=["001", "001"], lon=[112.1, 112.2], lat=[30.5, 30.6],
        ))
        @test_throws ArgumentError read_station_coordinates(duplicate_path)
    end
end
