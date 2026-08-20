#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "ERA5LandStations.jl"))

using .ERA5LandStations
using DataFrames

function usage()
    println("""
Download hourly ERA5-Land time series for the 237 study-area stations.

Usage:
  julia --project=. scripts/download_era5_land_stations.jl [options]

Options:
  --dry-run              Validate stations and write request JSON files only
  --limit N              Download only the first N stations (for testing)
  --python PATH          Python executable (default: python)
  --station-meta PATH    Station metadata CSV
  --raw-dir PATH         Raw station download directory
  --audit-dir PATH       Request and status audit directory
  -h, --help             Show this help

Authentication:
  Set CDSAPI_TOKEN in the environment, preferably through the secure wrapper.
""")
end

function parse_args(args)
    options = Dict{Symbol,Any}(
        :dry_run => false,
        :limit => nothing,
        :python => "python",
        :station_meta => joinpath(
            ROOT, "data", "processed", "study_area", "station_meta.csv",
        ),
        :raw_dir => joinpath(
            ROOT, "data", "raw", "era5_land", "stations_2022_2024",
        ),
        :audit_dir => joinpath(
            ROOT, "output", "input_audit", "era5_land_stations_2022_2024",
        ),
        :help => false,
    )
    value_options = Dict(
        "--limit" => :limit,
        "--python" => :python,
        "--station-meta" => :station_meta,
        "--raw-dir" => :raw_dir,
        "--audit-dir" => :audit_dir,
    )

    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("-h", "--help")
            options[:help] = true
        elseif arg == "--dry-run"
            options[:dry_run] = true
        elseif haskey(value_options, arg)
            index == length(args) && throw(ArgumentError("missing value after $arg"))
            index += 1
            key = value_options[arg]
            options[key] = key == :limit ? parse(Int, args[index]) : args[index]
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        index += 1
    end
    options[:limit] === nothing || options[:limit] > 0 ||
        throw(ArgumentError("--limit must be positive"))
    return options
end

function main(args=ARGS)
    options = parse_args(args)
    if options[:help]
        usage()
        return nothing
    end

    stations = read_era5_stations(options[:station_meta])
    println("Validated $(nrow(stations)) stations in 109.4-111.6 E, 31.2-33.4 N")
    println("Date range: 2022-01-01 through 2024-12-31 (hourly UTC)")
    println("Variables: $(join(ERA5_LAND_VARIABLES, ", "))")

    selected = options[:limit] === nothing ? stations :
        first(stations, min(options[:limit], nrow(stations)))
    if options[:dry_run]
        _, manifest_path = write_station_requests(
            selected, options[:audit_dir], options[:raw_dir],
        )
        println("Dry run complete: $(nrow(selected)) requests")
        println("Manifest: $manifest_path")
        return manifest_path
    end

    token = strip(get(ENV, "CDSAPI_TOKEN", ""))
    isempty(token) && error("missing CDSAPI_TOKEN; use the secure PowerShell wrapper")
    success(`$(options[:python]) -c "import cdsapi"`) || error(
        "Python package cdsapi is unavailable. Run: " *
        "python -m pip install \"cdsapi>=0.7.7\"",
    )
    manifest = download_era5_land_stations(
        stations;
        raw_dir=options[:raw_dir],
        audit_dir=options[:audit_dir],
        token,
        python=options[:python],
        limit=options[:limit],
    )
    println("Completed $(count(==("complete"), manifest.status)) station downloads.")
    return manifest
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
