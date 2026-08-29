#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("ERA5LandStations", "ERA5LandProcessing")

using Main.ERA5LandStations
using Main.ERA5LandProcessing
using CSV, DataFrames

function parse_args(args)
    options = Dict{Symbol,Any}(
        :limit => nothing,
        :overwrite => false,
        :station_meta => joinpath(
            ROOT, "data", "processed", "study_area", "station_meta.csv",
        ),
        :raw_dir => joinpath(
            ROOT, "data", "raw", "era5_land", "stations_2022_2024",
        ),
        :output_dir => joinpath(
            ROOT, "data", "processed", "era5_land", "stations_2022_2024",
        ),
        :audit_path => joinpath(
            ROOT, "output", "validation", "era5_land_station_processing.csv",
        ),
    )
    value_options = Dict(
        "--limit" => :limit,
        "--station-meta" => :station_meta,
        "--raw-dir" => :raw_dir,
        "--output-dir" => :output_dir,
        "--audit-path" => :audit_path,
    )
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--overwrite"
            options[:overwrite] = true
        elseif haskey(value_options, argument)
            index == length(args) && throw(ArgumentError("missing value after $argument"))
            index += 1
            key = value_options[argument]
            options[key] = key == :limit ? parse(Int, args[index]) : args[index]
        else
            throw(ArgumentError("unknown argument: $argument"))
        end
        index += 1
    end
    options[:limit] === nothing || options[:limit] > 0 ||
        throw(ArgumentError("--limit must be positive"))
    return options
end

function main(args=ARGS)
    options = parse_args(args)
    stations = read_era5_stations(options[:station_meta])
    selected = options[:limit] === nothing ? stations :
        first(stations, min(options[:limit], nrow(stations)))
    audit = DataFrame(
        station_id=selected.station_id,
        status=fill("pending", nrow(selected)),
        output_path=[joinpath(options[:output_dir], string(id, ".csv"))
            for id in selected.station_id],
    )

    for (index, row) in enumerate(eachrow(selected))
        output_path = audit.output_path[index]
        if isfile(output_path) && !options[:overwrite]
            audit.status[index] = "complete"
            println("[$index/$(nrow(selected))] Already processed: $(row.station_id)")
            continue
        end
        station_dir = joinpath(options[:raw_dir], row.station_id)
        if !isfile(joinpath(station_dir, ".complete"))
            audit.status[index] = "raw_pending"
            continue
        end
        println("[$index/$(nrow(selected))] Processing station $(row.station_id)")
        process_era5_station_directory(
            station_dir, output_path;
            station_id=row.station_id,
            station_lon=row.lon,
            station_lat=row.lat,
        )
        audit.status[index] = "complete"
    end

    mkpath(dirname(options[:audit_path]))
    CSV.write(options[:audit_path], audit)
    println("Processed $(count(==("complete"), audit.status)) stations; " *
        "$(count(==("raw_pending"), audit.status)) raw downloads pending.")
    return audit
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error_value
        showerror(stderr, error_value)
        println(stderr)
        exit(1)
    end
end
