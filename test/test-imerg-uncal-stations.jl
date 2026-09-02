using CSV, DataFrames, Dates, JSON, NCDatasets, Test

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("IMERGUncalStations")
using Main.IMERGUncalStations

function write_imerg_fixture(path, spec, start_time)
    NCDataset(path, "c") do dataset
        defDim(dataset, "lat", spec.lat_last - spec.lat_first + 1)
        defDim(dataset, "lon", spec.lon_last - spec.lon_first + 1)
        defDim(dataset, "time", 1)
        latitude = defVar(dataset, "lat", Float32, ("lat",))
        longitude = defVar(dataset, "lon", Float32, ("lon",))
        time = defVar(dataset, "time", Int32, ("time",))
        precipitation = defVar(
            dataset,
            "precipitationUncal",
            Float32,
            ("lat", "lon", "time");
            attrib=Dict("units" => "mm/hr"),
            fillvalue=-9999.9f0,
        )
        latitude[:] = Float32.(-89.95 .+ 0.1 .* collect(spec.lat_first:spec.lat_last))
        longitude[:] = Float32.(-179.95 .+ 0.1 .* collect(spec.lon_first:spec.lon_last))
        time[:] = Int32[0]
        precipitation[:, :, 1] = reshape(
            Float32.(1:length(latitude) * length(longitude)),
            length(latitude), length(longitude),
        )
        dataset.attrib["FileHeader"] = "StartGranuleDateTime=" *
            Dates.format(start_time, dateformat"yyyy-mm-ddTHH:MM:SS") * ".000Z;"
    end
    return path
end

@testset "IMERG precipitationUncal station extraction" begin
    mktempdir() do temp_dir
        station_path = joinpath(temp_dir, "stations.csv")
        stations_source = DataFrame(
            station_id=[lpad(string(index), 8, '0') for index in 1:3],
            lon=[108.672, 110.312, 116.029],
            lat=[29.232, 31.25, 33.2],
        )
        CSV.write(station_path, stations_source)
        stations = read_imerg_stations(station_path; expected_count=3)
        @test stations.station_id == ["00000001", "00000002", "00000003"]
        @test_throws ArgumentError read_imerg_stations(station_path)

        spec = build_subset_spec(stations)
        @test spec.lon_first == 2886
        @test spec.lon_last == 2960
        @test spec.lat_first == 1192
        @test spec.lat_last == 1232
        constraint = build_subset_constraint(spec)
        @test occursin("/Grid/Intermediate/precipitationUncal", constraint)
        @test occursin("[2886:1:2960][1192:1:1232]", constraint)

        start_time = DateTime(2022, 1, 1)
        fixture = write_imerg_fixture(joinpath(temp_dir, "subset.nc4"), spec, start_time)
        @test validate_subset_file(fixture, spec; expected_time=start_time)
        values = read_station_values(fixture, stations, spec; expected_time=start_time)
        @test length(values) == 3
        @test all(!ismissing, values)
        @test Float64.(values) == [1.0, 718.0, 3075.0]
        hourly = aggregate_hour_values(values, values)
        @test hourly == values
        @test_throws ErrorException aggregate_hour_values([1.0], [missing])
        @test isequal(
            aggregate_hour_values([1.0], [missing]; allow_missing=true),
            [missing],
        )

        second_fixture = write_imerg_fixture(
            joinpath(temp_dir, "subset_0030.nc4"), spec, start_time + Minute(30),
        )
        granules = [
            GranuleRef(start_time, "first.HDF5", "https://example.test/first"),
            GranuleRef(start_time + Minute(30), "second.HDF5", "https://example.test/second"),
        ]
        monthly = IMERGUncalStations._month_table(
            granules, [fixture, second_fixture], stations, spec,
        )
        @test nrow(monthly) == 3
        @test unique(monthly.time) == ["2022-01-01T00:00:00Z"]
        month_path = joinpath(temp_dir, "gpm_hubei_hourly_long_202201.csv")
        CSV.write(month_path, monthly)
        summary = validate_month_csv(
            month_path, stations, start_time, start_time + Hour(1),
        )
        @test summary.rows == 3
        @test summary.hours == 1
        @test summary.stations == 3
    end
end

@testset "IMERG CMR catalog checks" begin
    ids = [
        "3B-HHR.MS.MRG.3IMERG.20220101-S000000-E002959.0000.V07B.HDF5",
        "3B-HHR.MS.MRG.3IMERG.20220101-S003000-E005959.0030.V07B.HDF5",
    ]
    entries = [
        Dict(
            "producer_granule_id" => id,
            "links" => [Dict(
                "href" => "https://opendap.earthdata.nasa.gov/collections/test/granules/" * id,
            )],
        ) for id in ids
    ]
    payload = Dict("feed" => Dict("entry" => entries))
    start_time = DateTime(2022, 1, 1)
    stop_time = start_time + Hour(1)
    catalog = parse_cmr_catalog(payload, start_time, stop_time)
    @test length(catalog) == 2
    @test parse_granule_start(ids[2]) == start_time + Minute(30)
    @test validate_catalog(catalog, start_time, stop_time)
    @test_throws ErrorException validate_catalog(first(catalog, 1), start_time, stop_time)

    spec = SubsetSpec(2886, 2960, 1192, 1232)
    url = build_subset_url(first(catalog), spec)
    @test occursin(".nc4?%2FGrid%2FIntermediate%2FprecipitationUncal", url)
    @test occursin("%5B2886%3A1%3A2960%5D", url)
    @test expected_half_hours(DateTime(2022), DateTime(2025)) == 52_608
    @test expected_hours(DateTime(2022), DateTime(2025)) == 26_304

    csv_body = "Producer Granule ID\n\"$(ids[1])\"\n\"$(ids[2])\"\n"
    compact_catalog = parse_cmr_csv(csv_body, start_time, stop_time)
    @test [item.producer_id for item in compact_catalog] == ids
    @test occursin("GPM_3IMERGHH.07%3A", first(compact_catalog).opendap_url)
    @test IMERGUncalStations._retryable_status(200)
    @test !IMERGUncalStations._retryable_status(403)
end
