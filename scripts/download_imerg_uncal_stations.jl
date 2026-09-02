#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("IMERGUncalStations")

using Main.IMERGUncalStations
using Dates
using DataFrames

function usage()
    println("""
Download NASA GPM IMERG V07 precipitationUncal for all 318 Hubei stations.

The script requests only the 0.1-degree Hubei grid rectangle from OPeNDAP,
retains every validated half-hour subset, and creates monthly hourly long CSVs.

Usage:
  julia --project=. scripts/download_imerg_uncal_stations.jl [options]

Modes:
  --dry-run                    Validate local inputs and print the complete plan
  --catalog-check              Verify all monthly CMR granules, then exit
  --validate-hour TIMESTAMP    Download and verify one UTC hour, then exit
  (no mode flag)               Download complete calendar months

Options:
  --start YYYY-MM-DD           First month, inclusive (default: 2022-01-01)
  --end YYYY-MM-DD             Month boundary, exclusive (default: 2025-01-01)
  --station-meta PATH          Station CSV (default: data/hubei_station_meta.csv)
  --raw-dir PATH               Retained NASA subset files
  --output-dir PATH            Monthly long CSV directory
  --audit-dir PATH             CMR catalogs, QC, and validation output
  --staging-dir PATH           Temporary download directory outside data/raw
  --concurrency N              Concurrent NASA requests (default: 4)
  --max-attempts N             Attempts per subset (default: 4)
  --overwrite-output           Replace an existing invalid/monthly output CSV
  -h, --help                   Show this help

Authentication:
  Set EARTHDATA_TOKEN in the environment. The token is never written to disk.

Examples:
  julia --project=. scripts/download_imerg_uncal_stations.jl --dry-run
  julia --project=. scripts/download_imerg_uncal_stations.jl --validate-hour 2022-01-01T00:00:00
  julia --project=. scripts/download_imerg_uncal_stations.jl
""")
end

function _date(value)
    try
        return Date(value, dateformat"yyyy-mm-dd")
    catch
        throw(ArgumentError("invalid date '$value'; expected YYYY-MM-DD"))
    end
end

function _datetime(value)
    text = replace(strip(value), "Z" => "")
    for format in (dateformat"yyyy-mm-ddTHH:MM:SS", dateformat"yyyy-mm-ddTHH:MM")
        try
            return DateTime(text, format)
        catch
        end
    end
    throw(ArgumentError(
        "invalid timestamp '$value'; expected YYYY-MM-DDTHH:MM[:SS]",
    ))
end

function parse_args(args)
    audit_default = joinpath(
        ROOT, "output", "input_audit", "gpm_imerg_uncal_stations_2022_2024",
    )
    options = Dict{Symbol,Any}(
        :dry_run => false,
        :catalog_check => false,
        :validate_hour => nothing,
        :start => Date(2022, 1, 1),
        :stop => Date(2025, 1, 1),
        :station_meta => joinpath(ROOT, "data", "hubei_station_meta.csv"),
        :raw_dir => joinpath(
            ROOT, "data", "raw", "gpm_imerg_v07_uncal_hubei_2022_2024",
        ),
        :output_dir => joinpath(
            ROOT, "output", "gpm_imerg_uncal_stations_2022_2024",
        ),
        :audit_dir => audit_default,
        :staging_dir => joinpath(audit_default, "staging"),
        :concurrency => 4,
        :max_attempts => 4,
        :overwrite_output => false,
        :help => false,
    )
    value_options = Dict(
        "--validate-hour" => :validate_hour,
        "--start" => :start,
        "--end" => :stop,
        "--station-meta" => :station_meta,
        "--raw-dir" => :raw_dir,
        "--output-dir" => :output_dir,
        "--audit-dir" => :audit_dir,
        "--staging-dir" => :staging_dir,
        "--concurrency" => :concurrency,
        "--max-attempts" => :max_attempts,
    )

    index = 1
    while index <= length(args)
        argument = args[index]
        if argument in ("-h", "--help")
            options[:help] = true
        elseif argument == "--dry-run"
            options[:dry_run] = true
        elseif argument == "--catalog-check"
            options[:catalog_check] = true
        elseif argument == "--overwrite-output"
            options[:overwrite_output] = true
        elseif haskey(value_options, argument)
            index == length(args) && throw(ArgumentError("missing value after $argument"))
            index += 1
            key = value_options[argument]
            value = args[index]
            options[key] = if key == :validate_hour
                _datetime(value)
            elseif key in (:start, :stop)
                _date(value)
            elseif key in (:concurrency, :max_attempts)
                parse(Int, value)
            else
                value
            end
        elseif startswith(argument, "--") && occursin('=', argument)
            name, value = split(argument, '='; limit=2)
            haskey(value_options, name) || throw(ArgumentError("unknown option: $name"))
            key = value_options[name]
            options[key] = if key == :validate_hour
                _datetime(value)
            elseif key in (:start, :stop)
                _date(value)
            elseif key in (:concurrency, :max_attempts)
                parse(Int, value)
            else
                value
            end
        else
            throw(ArgumentError("unknown argument: $argument"))
        end
        index += 1
    end

    options[:concurrency] > 0 || throw(ArgumentError("--concurrency must be positive"))
    options[:max_attempts] > 0 || throw(ArgumentError("--max-attempts must be positive"))
    options[:start] < options[:stop] || throw(ArgumentError("--start must precede --end"))
    for key in (:start, :stop)
        date = options[key]
        day(date) == 1 || throw(ArgumentError("--$(key) must be the first day of a month"))
    end
    options[:dry_run] && options[:validate_hour] !== nothing &&
        throw(ArgumentError("--dry-run and --validate-hour are separate modes"))
    options[:dry_run] && options[:catalog_check] &&
        throw(ArgumentError("--dry-run and --catalog-check are separate modes"))
    options[:catalog_check] && options[:validate_hour] !== nothing &&
        throw(ArgumentError("--catalog-check and --validate-hour are separate modes"))
    return options
end

function _month_starts(start_date::Date, stop_date::Date)
    starts = DateTime[]
    current = DateTime(start_date)
    boundary = DateTime(stop_date)
    while current < boundary
        push!(starts, current)
        current += Month(1)
    end
    current == boundary || error("date range does not end on a calendar-month boundary")
    return starts
end

function _preflight_catalogs(months)
    catalogs = Vector{Vector{GranuleRef}}(undef, length(months))
    total = 0
    for (index, month_start) in enumerate(months)
        month_stop = month_start + Month(1)
        catalog = fetch_cmr_catalog(month_start, month_stop)
        catalogs[index] = catalog
        total += length(catalog)
        println("  catalog $(Dates.format(month_start, dateformat"yyyy-mm")): " *
                "$(length(catalog)) granules OK")
    end
    expected = expected_half_hours(first(months), last(months) + Month(1))
    total == expected || error("catalog total $total does not equal expected $expected")
    println("Catalog preflight passed: $(length(months)) months, $total granules.")
    return catalogs
end

function _token()
    token = strip(get(ENV, "EARTHDATA_TOKEN", ""))
    isempty(token) && error(
        "EARTHDATA_TOKEN is missing. Configure it in the environment and restart the terminal/app.",
    )
    return token
end

function _print_plan(options, stations)
    start_time = DateTime(options[:start])
    stop_time = DateTime(options[:stop])
    spec = build_subset_spec(stations)
    half_hours = expected_half_hours(start_time, stop_time)
    hours = expected_hours(start_time, stop_time)
    rows = hours * nrow(stations)
    println("Validated stations: $(nrow(stations)) unique IDs")
    println("Station bounds: $(minimum(stations.lon))-$(maximum(stations.lon)) E, " *
            "$(minimum(stations.lat))-$(maximum(stations.lat)) N")
    println("IMERG subset: lon indices $(spec.lon_first):$(spec.lon_last), " *
            "lat indices $(spec.lat_first):$(spec.lat_last)")
    println("Date range: $start_time UTC (inclusive) to $stop_time UTC (exclusive)")
    println("Expected: $(length(_month_starts(options[:start], options[:stop]))) months, " *
            "$half_hours half-hour granules, $hours hours, $rows station rows")
    println("Variable: /Grid/Intermediate/precipitationUncal " *
            "($(IMERGUncalStations.IMERG_UNITS))")
    println("Hourly formula: 0.5 * minute-00 rate + 0.5 * minute-30 rate")
    println("Raw subsets: $(options[:raw_dir])")
    println("Monthly CSVs: $(options[:output_dir])")
    println("Audit/QC: $(options[:audit_dir])")
end

function main(args=ARGS)
    options = parse_args(args)
    if options[:help]
        usage()
        return nothing
    end
    stations = read_imerg_stations(options[:station_meta])
    _print_plan(options, stations)
    if options[:dry_run]
        println("Dry run complete; no network request or file download was made.")
        return nothing
    end

    months = _month_starts(options[:start], options[:stop])
    if options[:catalog_check]
        _preflight_catalogs(months)
        return nothing
    end

    token = _token()
    if options[:validate_hour] !== nothing
        output_path, audit_path, audit = validate_remote_hour(
            options[:validate_hour],
            stations;
            token,
            audit_dir=options[:audit_dir],
            concurrency=min(options[:concurrency], 2),
            max_attempts=options[:max_attempts],
        )
        println("Remote validation passed: $(audit["granules"]) granules, " *
                "$(audit["stations"]) stations, $(audit["missing_values"]) missing values")
        println("Validation CSV: $output_path")
        println("Validation audit: $audit_path")
        return output_path
    end

    println("Preflighting every monthly NASA CMR catalog before downloading...")
    catalogs = _preflight_catalogs(months)
    outputs = String[]
    for (index, month_start) in enumerate(months)
        println("[$index/$(length(months))] Processing $(Dates.format(month_start, dateformat"yyyy-mm"))")
        output_path, summary = process_imerg_month(
            month_start,
            stations;
            raw_dir=options[:raw_dir],
            output_dir=options[:output_dir],
            audit_dir=options[:audit_dir],
            staging_dir=options[:staging_dir],
            token,
            concurrency=options[:concurrency],
            max_attempts=options[:max_attempts],
            overwrite_output=options[:overwrite_output],
            catalog=catalogs[index],
        )
        push!(outputs, output_path)
        println("  valid: $(summary.hours) hours × $(summary.stations) stations")
    end
    println("Completed and validated $(length(outputs)) monthly CSV files.")
    return outputs
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error_value
        showerror(stderr, error_value, catch_backtrace())
        println(stderr)
        exit(1)
    end
end
