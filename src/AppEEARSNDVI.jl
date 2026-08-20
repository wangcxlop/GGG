module AppEEARSNDVI

using Base64, CSV, DataFrames, Downloads, JSON, SHA

export MOD13A2_PRODUCT, MOD13A2_LAYERS,
    read_station_coordinates, validate_coordinate_bounds, build_mod13a2_point_task,
    write_request_json, appeears_token, validate_product_layers,
    submit_task, wait_for_task, download_task_bundle

const API_BASE = "https://appeears.earthdatacloud.nasa.gov/api"
const MOD13A2_PRODUCT = "MOD13A2.061"
const MOD13A2_LAYERS = [
    "_1_km_16_days_NDVI",
    "_1_km_16_days_VI_Quality",
    "_1_km_16_days_pixel_reliability",
    "_1_km_16_days_composite_day_of_the_year",
]

function _column(df::DataFrame, requested::Symbol)
    columns = Dict(Symbol(lowercase(String(name))) => name for name in names(df))
    column = get(columns, Symbol(lowercase(String(requested))), nothing)
    column === nothing && throw(ArgumentError("missing column: $requested"))
    return column
end

"""Require all point coordinates to lie inside an inclusive lon/lat rectangle."""
function validate_coordinate_bounds(
    coordinates;
    west::Real,
    east::Real,
    south::Real,
    north::Real,
)
    west < east || throw(ArgumentError("west must be less than east"))
    south < north || throw(ArgumentError("south must be less than north"))
    outside = [
        String(coordinate["id"]) for coordinate in coordinates
        if !(west <= coordinate["longitude"] <= east &&
             south <= coordinate["latitude"] <= north)
    ]
    isempty(outside) || throw(ArgumentError(
        "$(length(outside)) stations fall outside [$west, $east] x [$south, $north]: " *
        join(first(outside, min(10, length(outside))), ", "),
    ))
    return true
end

"""Read and validate AppEEARS point coordinates from station metadata."""
function read_station_coordinates(
    path::AbstractString;
    station_id_col::Symbol=:station_id,
    lon_col::Symbol=:lon,
    lat_col::Symbol=:lat,
)
    isfile(path) || throw(ArgumentError("station metadata not found: $path"))
    stations = CSV.read(path, DataFrame; types=Dict(station_id_col => String))
    id_col = _column(stations, station_id_col)
    x_col = _column(stations, lon_col)
    y_col = _column(stations, lat_col)

    ids = strip.(string.(stations[!, id_col]))
    any(isempty, ids) && throw(ArgumentError("station metadata contains an empty station id"))
    allunique(ids) || throw(ArgumentError("station metadata contains duplicate station ids"))

    coordinates = Vector{Dict{String,Any}}(undef, nrow(stations))
    for row in 1:nrow(stations)
        lon = Float64(stations[row, x_col])
        lat = Float64(stations[row, y_col])
        isfinite(lon) && -180.0 <= lon <= 180.0 ||
            throw(ArgumentError("invalid longitude for station $(ids[row]): $lon"))
        isfinite(lat) && -90.0 <= lat <= 90.0 ||
            throw(ArgumentError("invalid latitude for station $(ids[row]): $lat"))
        coordinates[row] = Dict(
            "id" => ids[row],
            "category" => "hubei_rainfall_station",
            "longitude" => lon,
            "latitude" => lat,
        )
    end
    return coordinates
end

"""Build the full-year 2022--2024 MOD13A2 station point request."""
function build_mod13a2_point_task(
    coordinates;
    task_name::AbstractString="MOD13A2_NDVI_Hubei_2022_2024",
    start_date::AbstractString="01-01-2022",
    end_date::AbstractString="12-31-2024",
    layers::AbstractVector{<:AbstractString}=MOD13A2_LAYERS,
)
    isempty(coordinates) && throw(ArgumentError("at least one station coordinate is required"))
    isempty(layers) && throw(ArgumentError("at least one MOD13A2 layer is required"))
    return Dict(
        "task_type" => "point",
        "task_name" => String(task_name),
        "params" => Dict(
            "dates" => [Dict(
                "startDate" => String(start_date),
                "endDate" => String(end_date),
                "recurring" => false,
            )],
            "layers" => [
                Dict("product" => MOD13A2_PRODUCT, "layer" => String(layer))
                for layer in layers
            ],
            "coordinates" => coordinates,
        ),
    )
end

function write_request_json(path::AbstractString, task)
    mkpath(dirname(path))
    temp_path = string(path, ".tmp-", getpid())
    open(temp_path, "w") do io
        JSON.print(io, task, 2)
        write(io, '\n')
    end
    mv(temp_path, path; force=true)
    return path
end

_response_status(response::Downloads.Response) = response.status
_response_status(response::Downloads.RequestError) = response.response.status
_retryable_status(status::Integer) = status == 0 || status in (408, 425, 429, 500, 502, 503, 504)

function _request_json(
    endpoint::AbstractString;
    method::AbstractString="GET",
    token::Union{Nothing,AbstractString}=nothing,
    basic_auth::Union{Nothing,AbstractString}=nothing,
    body=nothing,
    max_attempts::Integer=6,
    retry_delay::Real=10,
)
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive"))
    retry_delay >= 0 || throw(ArgumentError("retry_delay must be nonnegative"))
    headers = Pair{String,String}["Accept" => "application/json"]
    token === nothing || push!(headers, "Authorization" => "Bearer $token")
    basic_auth === nothing || push!(headers, "Authorization" => "Basic $basic_auth")
    body_text = nothing
    if body !== nothing
        push!(headers, "Content-Type" => "application/json")
        body_text = JSON.json(body)
    end

    for attempt in 1:max_attempts
        input = body_text === nothing ? nothing : IOBuffer(body_text)
        output = IOBuffer()
        response = Downloads.request(
            string(API_BASE, endpoint);
            method,
            headers,
            input,
            output,
            timeout=60,
            throw=false,
        )
        response_body = String(take!(output))
        status = _response_status(response)
        retryable = _retryable_status(status)

        if response isa Downloads.RequestError || !(200 <= status < 300)
            detail = isempty(strip(response_body)) ?
                (response isa Downloads.RequestError ? response.message : "no response body") :
                strip(response_body)
            if retryable && attempt < max_attempts
                println(stderr,
                    "Transient AppEEARS request failure (attempt $attempt/$max_attempts, " *
                    "HTTP $status); retrying in $(retry_delay) seconds.",
                )
                sleep(Float64(retry_delay))
                continue
            end
            error("AppEEARS $method $endpoint failed with HTTP $status: $detail")
        end
        return isempty(strip(response_body)) ? Dict{String,Any}() : JSON.parse(response_body)
    end
    error("unreachable AppEEARS request state")
end

"""
Return an AppEEARS bearer token. `APPEEARS_TOKEN` takes precedence; otherwise
`EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD` are used without printing them.
"""
function appeears_token()
    existing = strip(get(ENV, "APPEEARS_TOKEN", ""))
    !isempty(existing) && return existing

    username = strip(get(ENV, "EARTHDATA_USERNAME", ""))
    password = get(ENV, "EARTHDATA_PASSWORD", "")
    isempty(username) && error("missing EARTHDATA_USERNAME (or APPEEARS_TOKEN)")
    isempty(password) && error("missing EARTHDATA_PASSWORD (or APPEEARS_TOKEN)")

    basic = base64encode(string(username, ":", password))
    response = _request_json("/login"; method="POST", basic_auth=basic)
    haskey(response, "token") || error("AppEEARS login response did not contain a token")
    return String(response["token"])
end

function validate_product_layers(layers=MOD13A2_LAYERS)
    product = _request_json("/product/$MOD13A2_PRODUCT")
    missing_layers = [String(layer) for layer in layers if !haskey(product, String(layer))]
    isempty(missing_layers) || error(
        "AppEEARS no longer exposes requested MOD13A2 layers: $(join(missing_layers, ", "))",
    )
    unavailable = [
        String(layer) for layer in layers
        if !Bool(get(product[String(layer)], "Available", false))
    ]
    isempty(unavailable) || error(
        "AppEEARS reports unavailable MOD13A2 layers: $(join(unavailable, ", "))",
    )
    return true
end

function submit_task(token::AbstractString, task)
    response = _request_json("/task"; method="POST", token, body=task)
    haskey(response, "task_id") || error("AppEEARS submit response did not contain task_id")
    return String(response["task_id"])
end

function wait_for_task(
    token::AbstractString,
    task_id::AbstractString;
    poll_seconds::Real=30,
)
    poll_seconds > 0 || throw(ArgumentError("poll_seconds must be positive"))
    last_status = ""
    while true
        task = _request_json("/task/$task_id"; token)
        status = String(get(task, "status", "unknown"))
        if status != last_status
            println("AppEEARS task $task_id: $status")
            last_status = status
        end
        status == "done" && return task
        if status == "error"
            detail = get(task, "error", "unknown AppEEARS processing error")
            error("AppEEARS task $task_id failed: $detail")
        end
        status in ("pending", "processing", "queued") ||
            error("AppEEARS task $task_id returned unexpected status: $status")
        sleep(Float64(poll_seconds))
    end
end

function _sha256(path::AbstractString)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function _download_bundle_file(
    token::AbstractString,
    task_id::AbstractString,
    file::AbstractDict,
    outdir::AbstractString;
    overwrite::Bool=false,
)
    remote_name = String(file["file_name"])
    file_name = basename(remote_name)
    file_name == remote_name || error("unsafe AppEEARS bundle filename: $remote_name")
    destination = joinpath(outdir, file_name)
    expected_sha = lowercase(String(get(file, "sha256", "")))

    if isfile(destination) && !overwrite
        if isempty(expected_sha) || _sha256(destination) == expected_sha
            println("Already downloaded: $destination")
            return destination
        end
        error("existing file failed checksum; rerun with --overwrite: $destination")
    end

    temp_path = string(destination, ".part-", getpid())
    headers = ["Authorization" => "Bearer $token"]
    for attempt in 1:6
        response = open(temp_path, "w") do io
            Downloads.request(
                "$API_BASE/bundle/$task_id/$(file["file_id"])";
                headers,
                output=io,
                timeout=600,
                throw=false,
            )
        end
        status = _response_status(response)
        if response isa Downloads.RequestError || !(200 <= status < 300)
            rm(temp_path; force=true)
            if _retryable_status(status) && attempt < 6
                println("Transient download failure for $remote_name " *
                    "(attempt $attempt/6, HTTP $status); retrying in 10 seconds.")
                sleep(10)
                continue
            end
            detail = response isa Downloads.RequestError ? response.message : "request failed"
            error("failed to download $remote_name: HTTP $status: $detail")
        end
        if !isempty(expected_sha) && _sha256(temp_path) != expected_sha
            rm(temp_path; force=true)
            error("SHA-256 verification failed for $remote_name")
        end
        mv(temp_path, destination; force=true)
        println("Downloaded: $destination")
        return destination
    end
    error("unreachable AppEEARS download state")
end

function download_task_bundle(
    token::AbstractString,
    task_id::AbstractString,
    outdir::AbstractString;
    overwrite::Bool=false,
)
    mkpath(outdir)
    bundle = _request_json("/bundle/$task_id"; token)
    files = get(bundle, "files", Any[])
    isempty(files) && error("AppEEARS task $task_id has no downloadable files")
    return [
        _download_bundle_file(token, task_id, file, outdir; overwrite)
        for file in files
    ]
end

end
