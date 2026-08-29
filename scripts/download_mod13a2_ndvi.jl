#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("AppEEARSNDVI")

using Main.AppEEARSNDVI

function usage()
    println("""
Download full-year 2022--2024 MOD13A2.061 NDVI at project station locations.

Usage:
  julia --project=. scripts/download_mod13a2_ndvi.jl [options]

Options:
  --dry-run                 Validate stations and write request JSON only
  --submit-only             Submit the request and print/save task id, then exit
  --task-id ID              Resume polling/downloading an existing AppEEARS task
  --station-meta PATH       Station CSV (default: processed study-area stations)
  --output-dir PATH         Raw bundle directory
  --request-json PATH       Local request/audit JSON path
  --task-id-file PATH       File used to save the submitted task id
  --poll-seconds N          Poll interval in seconds (default: 30)
  --overwrite               Replace existing bundle files
  -h, --help                Show this help

Authentication (environment variables; never pass passwords as arguments):
  APPEEARS_TOKEN            Existing AppEEARS bearer token, or
  EARTHDATA_USERNAME        NASA Earthdata Login username
  EARTHDATA_PASSWORD        NASA Earthdata Login password
""")
end

function parse_args(args)
    options = Dict{Symbol,Any}(
        :dry_run => false,
        :submit_only => false,
        :task_id => nothing,
        :station_meta => joinpath(
            ROOT, "data", "processed", "study_area", "station_meta.csv",
        ),
        :output_dir => joinpath(ROOT, "data", "raw", "ndvi", "mod13a2_2022_2024"),
        :request_json => joinpath(
            ROOT, "output", "input_audit", "mod13a2_ndvi_2022_2024_request.json",
        ),
        :task_id_file => joinpath(
            ROOT, "output", "input_audit", "mod13a2_ndvi_2022_2024_task_id.txt",
        ),
        :poll_seconds => 30.0,
        :overwrite => false,
        :help => false,
    )

    value_options = Dict(
        "--task-id" => :task_id,
        "--station-meta" => :station_meta,
        "--output-dir" => :output_dir,
        "--request-json" => :request_json,
        "--task-id-file" => :task_id_file,
        "--poll-seconds" => :poll_seconds,
    )

    index = 1
    while index <= length(args)
        arg = args[index]
        if arg in ("-h", "--help")
            options[:help] = true
        elseif arg == "--dry-run"
            options[:dry_run] = true
        elseif arg == "--submit-only"
            options[:submit_only] = true
        elseif arg == "--overwrite"
            options[:overwrite] = true
        elseif haskey(value_options, arg)
            index == length(args) && throw(ArgumentError("missing value after $arg"))
            index += 1
            key = value_options[arg]
            options[key] = key == :poll_seconds ? parse(Float64, args[index]) : args[index]
        elseif startswith(arg, "--") && occursin('=', arg)
            name, value = split(arg, '='; limit=2)
            haskey(value_options, name) || throw(ArgumentError("unknown option: $name"))
            key = value_options[name]
            options[key] = key == :poll_seconds ? parse(Float64, value) : value
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        index += 1
    end
    options[:poll_seconds] > 0 || throw(ArgumentError("--poll-seconds must be positive"))
    options[:dry_run] && options[:task_id] !== nothing &&
        throw(ArgumentError("--dry-run cannot be combined with --task-id"))
    return options
end

function write_task_id(path::AbstractString, task_id::AbstractString)
    mkpath(dirname(path))
    temp_path = string(path, ".tmp-", getpid())
    open(temp_path, "w") do io
        println(io, task_id)
    end
    mv(temp_path, path; force=true)
    return path
end

function main(args=ARGS)
    options = parse_args(args)
    if options[:help]
        usage()
        return nothing
    end

    task_id = options[:task_id]
    if task_id === nothing
        coordinates = read_station_coordinates(options[:station_meta])
        validate_coordinate_bounds(
            coordinates; west=109.4, east=111.6, south=31.2, north=33.4,
        )
        task = build_mod13a2_point_task(coordinates)
        write_request_json(options[:request_json], task)
        println("Validated $(length(coordinates)) stations in 109.4-111.6 E, 31.2-33.4 N")
        println("Request JSON: $(options[:request_json])")
        println("Date range: 2022-01-01 through 2024-12-31")
        println("Layers: $(join(MOD13A2_LAYERS, ", "))")

        if options[:dry_run]
            println("Dry run complete; no request was submitted.")
            return options[:request_json]
        end

        println("Checking current AppEEARS product layers...")
        validate_product_layers()
        token = appeears_token()
        task_id = submit_task(token, task)
        write_task_id(options[:task_id_file], task_id)
        println("Submitted AppEEARS task: $task_id")
        println("Task id saved to: $(options[:task_id_file])")
        if options[:submit_only]
            return task_id
        end
    else
        token = appeears_token()
        println("Resuming AppEEARS task: $task_id")
    end

    wait_for_task(token, task_id; poll_seconds=options[:poll_seconds])
    files = download_task_bundle(
        token, task_id, options[:output_dir]; overwrite=options[:overwrite],
    )
    println("Downloaded $(length(files)) bundle files to: $(options[:output_dir])")
    return files
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
