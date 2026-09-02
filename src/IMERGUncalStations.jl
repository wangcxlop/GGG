module IMERGUncalStations

using CSV
using DataFrames
using Dates
using Downloads
using JSON
using NCDatasets
using SHA

export GranuleRef, SubsetSpec,
    read_imerg_stations, build_subset_spec, build_subset_constraint,
    build_subset_url, parse_granule_start, parse_cmr_catalog, parse_cmr_csv,
    fetch_cmr_catalog, validate_catalog, validate_subset_file,
    read_station_values, aggregate_hour_values, process_imerg_month,
    validate_remote_hour, validate_month_csv, expected_half_hours,
    expected_hours

const IMERG_COLLECTION_ID = "C2723754847-GES_DISC"
const IMERG_SHORT_NAME = "GPM_3IMERGHH"
const IMERG_VERSION = "07"
const IMERG_VARIABLE = "precipitationUncal"
const IMERG_UNITS = "mm/hr"
const HALF_HOUR = Minute(30)
const GRID_STEP = 0.1
const GRID_LON_ORIGIN = -179.95
const GRID_LAT_ORIGIN = -89.95
const GRID_NLON = 3600
const GRID_NLAT = 1800
const DEFAULT_EXPECTED_STATIONS = 318

struct GranuleRef
    start_time::DateTime
    producer_id::String
    opendap_url::String
end

struct SubsetSpec
    lon_first::Int
    lon_last::Int
    lat_first::Int
    lat_last::Int
end

expected_half_hours(start_time::DateTime, stop_time::DateTime) =
    _expected_steps(start_time, stop_time, HALF_HOUR)

expected_hours(start_time::DateTime, stop_time::DateTime) =
    _expected_steps(start_time, stop_time, Hour(1))

function _expected_steps(start_time::DateTime, stop_time::DateTime, step::Period)
    start_time < stop_time || throw(ArgumentError("start_time must precede stop_time"))
    milliseconds = Dates.value(stop_time - start_time)
    step_milliseconds = Dates.value(Millisecond(step))
    milliseconds % step_milliseconds == 0 || throw(ArgumentError(
        "time range must be an exact multiple of $step",
    ))
    return milliseconds ÷ step_milliseconds
end

function _column(df::DataFrame, requested::Symbol)
    mapping = Dict(Symbol(lowercase(String(name))) => name for name in names(df))
    column = get(mapping, Symbol(lowercase(String(requested))), nothing)
    column === nothing && throw(ArgumentError("missing station column: $requested"))
    return column
end

function _station_id(value)
    id = strip(string(value))
    isempty(id) && throw(ArgumentError("station_id cannot be empty"))
    occursin(r"^[A-Za-z0-9_-]+$", id) ||
        throw(ArgumentError("unsafe station_id: $id"))
    return id
end

"""Read and validate the full 318-station Hubei metadata table."""
function read_imerg_stations(
    path::AbstractString;
    expected_count::Integer=DEFAULT_EXPECTED_STATIONS,
    west::Real=108.672,
    east::Real=116.029,
    south::Real=29.232,
    north::Real=33.200,
)
    isfile(path) || throw(ArgumentError("station metadata not found: $path"))
    source = CSV.read(
        path, DataFrame;
        stringtype=String,
        types=Dict(:station_id => String),
    )
    id_column = _column(source, :station_id)
    lon_column = _column(source, :lon)
    lat_column = _column(source, :lat)
    stations = DataFrame(
        station_id=[_station_id(value) for value in source[!, id_column]],
        lon=Float64.(source[!, lon_column]),
        lat=Float64.(source[!, lat_column]),
    )

    nrow(stations) == expected_count || throw(ArgumentError(
        "expected $expected_count stations, found $(nrow(stations))",
    ))
    allunique(stations.station_id) ||
        throw(ArgumentError("duplicate station_id values found"))
    all(isfinite, stations.lon) ||
        throw(ArgumentError("station longitude contains non-finite values"))
    all(isfinite, stations.lat) ||
        throw(ArgumentError("station latitude contains non-finite values"))
    outside = .!((west .<= stations.lon .<= east) .&
                 (south .<= stations.lat .<= north))
    any(outside) && throw(ArgumentError(
        "$(count(outside)) stations fall outside $west-$east E, $south-$north N",
    ))
    return stations
end

function _nearest_grid_index(value::Real, origin::Real, count::Integer)
    isfinite(value) || throw(ArgumentError("grid coordinate must be finite"))
    index = floor(Int, (Float64(value) - origin) / GRID_STEP + 0.5)
    0 <= index < count || throw(ArgumentError("coordinate $value is outside the IMERG grid"))
    return index
end

"""Return the smallest IMERG grid rectangle containing every nearest station cell."""
function build_subset_spec(stations::DataFrame)
    nrow(stations) > 0 || throw(ArgumentError("station table is empty"))
    lon_indices = [_nearest_grid_index(value, GRID_LON_ORIGIN, GRID_NLON)
                   for value in stations.lon]
    lat_indices = [_nearest_grid_index(value, GRID_LAT_ORIGIN, GRID_NLAT)
                   for value in stations.lat]
    return SubsetSpec(
        minimum(lon_indices), maximum(lon_indices),
        minimum(lat_indices), maximum(lat_indices),
    )
end

function build_subset_constraint(spec::SubsetSpec)
    return "/Grid/Intermediate/$IMERG_VARIABLE" *
        "[0:1:0][$(spec.lon_first):1:$(spec.lon_last)]" *
        "[$(spec.lat_first):1:$(spec.lat_last)]"
end

function _percent_encode(value::AbstractString)
    io = IOBuffer()
    for byte in codeunits(value)
        character = Char(byte)
        unreserved = (UInt8('A') <= byte <= UInt8('Z')) ||
            (UInt8('a') <= byte <= UInt8('z')) ||
            (UInt8('0') <= byte <= UInt8('9')) ||
            character in ('-', '_', '.', '~')
        if unreserved
            write(io, byte)
        else
            print(io, '%', uppercase(string(byte; base=16, pad=2)))
        end
    end
    return String(take!(io))
end

build_subset_url(granule::GranuleRef, spec::SubsetSpec) =
    string(granule.opendap_url, ".nc4?", _percent_encode(build_subset_constraint(spec)))

function parse_granule_start(producer_id::AbstractString)
    matched = match(r"3IMERG\.(\d{8})-S(\d{6})-E\d{6}", producer_id)
    matched === nothing && throw(ArgumentError(
        "cannot parse IMERG start time from producer id: $producer_id",
    ))
    return DateTime(string(matched.captures[1], matched.captures[2]), dateformat"yyyymmddHHMMSS")
end

function _opendap_link(entry::AbstractDict)
    links = get(entry, "links", Any[])
    for link in links
        href = string(get(link, "href", ""))
        startswith(href, "https://opendap.earthdata.nasa.gov/") && return href
    end
    throw(ArgumentError(
        "CMR entry has no NASA OPeNDAP link: $(get(entry, "producer_granule_id", "unknown"))",
    ))
end

"""Parse legacy CMR JSON and keep granules whose start lies in `[start, stop)`."""
function parse_cmr_catalog(payload, start_time::DateTime, stop_time::DateTime)
    entries = get(get(payload, "feed", Dict{String,Any}()), "entry", Any[])
    granules = GranuleRef[]
    for entry in entries
        producer_id = string(get(entry, "producer_granule_id", ""))
        isempty(producer_id) && continue
        granule_start = parse_granule_start(producer_id)
        start_time <= granule_start < stop_time || continue
        push!(granules, GranuleRef(granule_start, producer_id, _opendap_link(entry)))
    end
    sort!(granules; by=granule -> granule.start_time)
    return granules
end

function _opendap_url(producer_id::AbstractString)
    granule_ur = string(IMERG_SHORT_NAME, '.', IMERG_VERSION, ':', producer_id)
    return "https://opendap.earthdata.nasa.gov/collections/$IMERG_COLLECTION_ID/granules/" *
        _percent_encode(granule_ur)
end

"""Parse the compact CMR CSV response used by production catalog checks."""
function parse_cmr_csv(body::AbstractString, start_time::DateTime, stop_time::DateTime)
    table = CSV.read(
        IOBuffer(body), DataFrame;
        select=["Producer Granule ID"],
        types=String,
        stringtype=String,
    )
    producer_column = Symbol("Producer Granule ID")
    granules = GranuleRef[]
    for producer_id in table[!, producer_column]
        granule_start = parse_granule_start(producer_id)
        start_time <= granule_start < stop_time || continue
        push!(granules, GranuleRef(
            granule_start,
            producer_id,
            _opendap_url(producer_id),
        ))
    end
    sort!(granules; by=granule -> granule.start_time)
    return granules
end

function _cmr_url(start_time::DateTime, stop_time::DateTime)
    stop_inclusive = stop_time - Millisecond(1)
    timestamp(value) = Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
    parameters = [
        "short_name" => IMERG_SHORT_NAME,
        "version" => IMERG_VERSION,
        "temporal" => string(timestamp(start_time), ',', timestamp(stop_inclusive)),
        "page_size" => "2000",
        "sort_key" => "start_date",
    ]
    query = join((string(_percent_encode(key), '=', _percent_encode(value))
                  for (key, value) in parameters), '&')
    return "https://cmr.earthdata.nasa.gov/search/granules.csv?$query"
end

function _request_cmr_csv(
    url::AbstractString;
    timeout::Real=120,
    max_attempts::Integer=8,
)
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive"))
    for attempt in 1:max_attempts
        try
            output = IOBuffer()
            response = Downloads.request(
                url;
                output,
                timeout,
                headers=["Accept" => "text/csv"],
            )
            response.status == 200 || error("CMR request failed with HTTP $(response.status)")
            return String(take!(output))
        catch error_value
            status = _response_status(error_value)
            if attempt == max_attempts || !_retryable_status(status)
                rethrow(error_value)
            end
            delay = min(30, 5 * 2^(attempt - 1))
            println(stderr,
                "CMR catalog request failed (HTTP $status), attempt $attempt/$max_attempts; " *
                "retrying in $delay seconds",
            )
            sleep(delay)
        end
    end
    error("unreachable CMR request state")
end

"""Fetch and strictly validate one CMR catalog range (maximum 2000 half-hours)."""
function fetch_cmr_catalog(start_time::DateTime, stop_time::DateTime)
    expected_half_hours(start_time, stop_time) <= 2000 || throw(ArgumentError(
        "one CMR query is limited to 2000 entries; request at most one calendar month",
    ))
    body = _request_cmr_csv(_cmr_url(start_time, stop_time))
    granules = parse_cmr_csv(body, start_time, stop_time)
    validate_catalog(granules, start_time, stop_time)
    return granules
end

function validate_catalog(
    granules::AbstractVector{GranuleRef},
    start_time::DateTime,
    stop_time::DateTime,
)
    expected = expected_half_hours(start_time, stop_time)
    length(granules) == expected || error(
        "CMR returned $(length(granules)) granules; expected $expected for " *
        "$start_time to $stop_time",
    )
    starts = [granule.start_time for granule in granules]
    expected_starts = collect(start_time:HALF_HOUR:(stop_time - HALF_HOUR))
    starts == expected_starts || error(
        "CMR granules are duplicated, missing, or not on the 30-minute UTC grid",
    )
    length(unique(granule.producer_id for granule in granules)) == expected ||
        error("CMR catalog contains duplicate producer ids")
    return true
end

function _hdf5_signature_ok(path::AbstractString)
    isfile(path) || return false
    filesize(path) >= 8 || return false
    open(path, "r") do io
        return read(io, 8) == UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]
    end
end

_grid_coordinate(origin::Real, index::Integer) = origin + GRID_STEP * index

function _validate_coordinates(values, expected_first, expected_last, name)
    isempty(values) && error("$name coordinate is empty")
    issorted(values) || error("$name coordinate is not ascending")
    isapprox(first(values), expected_first; atol=2e-4, rtol=0) ||
        error("unexpected first $name coordinate: $(first(values))")
    isapprox(last(values), expected_last; atol=2e-4, rtol=0) ||
        error("unexpected last $name coordinate: $(last(values))")
    length(values) == 1 || all(isapprox.(diff(Float64.(values)), GRID_STEP; atol=2e-4, rtol=0)) ||
        error("$name coordinate spacing is not 0.1 degree")
    return true
end

function _expected_header_timestamp(time::DateTime)
    return "StartGranuleDateTime=" *
        Dates.format(time, dateformat"yyyy-mm-ddTHH:MM:SS") * ".000Z;"
end

"""Validate signature, variable, units, dimensions, coordinates, and granule header."""
function validate_subset_file(
    path::AbstractString,
    spec::SubsetSpec;
    expected_time::Union{Nothing,DateTime}=nothing,
)
    _hdf5_signature_ok(path) || error("not a valid HDF5/NetCDF4 subset: $path")
    NCDataset(path, "r") do dataset
        for variable_name in ("lat", "lon", "time", IMERG_VARIABLE)
            haskey(dataset, variable_name) ||
                error("subset is missing variable $variable_name: $path")
        end
        latitude = dataset["lat"][:]
        longitude = dataset["lon"][:]
        precipitation = dataset[IMERG_VARIABLE]
        NCDatasets.dimnames(precipitation) == ("lat", "lon", "time") || error(
            "unexpected $IMERG_VARIABLE dimension order: " *
            string(NCDatasets.dimnames(precipitation)),
        )
        expected_size = (
            spec.lat_last - spec.lat_first + 1,
            spec.lon_last - spec.lon_first + 1,
            1,
        )
        size(precipitation) == expected_size || error(
            "unexpected $IMERG_VARIABLE size $(size(precipitation)); expected $expected_size",
        )
        units = string(get(precipitation.attrib, "units", ""))
        units == IMERG_UNITS || error(
            "unexpected $IMERG_VARIABLE units '$units'; expected '$IMERG_UNITS'",
        )
        _validate_coordinates(
            longitude,
            _grid_coordinate(GRID_LON_ORIGIN, spec.lon_first),
            _grid_coordinate(GRID_LON_ORIGIN, spec.lon_last),
            "longitude",
        )
        _validate_coordinates(
            latitude,
            _grid_coordinate(GRID_LAT_ORIGIN, spec.lat_first),
            _grid_coordinate(GRID_LAT_ORIGIN, spec.lat_last),
            "latitude",
        )
        if expected_time !== nothing
            header = string(get(dataset.attrib, "FileHeader", ""))
            occursin(_expected_header_timestamp(expected_time), header) || error(
                "subset FileHeader does not match expected granule start $expected_time",
            )
        end
    end
    return true
end

function _local_grid_index(
    values,
    target::Real,
    origin::Real,
    global_count::Integer,
    subset_first::Integer,
)
    global_index = _nearest_grid_index(target, origin, global_count)
    local_index = global_index - subset_first + 1
    1 <= local_index <= length(values) ||
        error("station coordinate $target falls outside the requested subset")
    expected_coordinate = _grid_coordinate(origin, global_index)
    isapprox(Float64(values[local_index]), expected_coordinate; atol=2e-4, rtol=0) ||
        error("subset coordinate does not match the expected nearest station cell")
    return local_index
end

"""Read nearest-cell `precipitationUncal` for every station in input order."""
function read_station_values(
    path::AbstractString,
    stations::DataFrame,
    spec::SubsetSpec;
    expected_time::Union{Nothing,DateTime}=nothing,
)
    validate_subset_file(path, spec; expected_time)
    return NCDataset(path, "r") do dataset
        latitude = dataset["lat"][:]
        longitude = dataset["lon"][:]
        precipitation = dataset[IMERG_VARIABLE][:, :, :]
        values = Vector{Union{Missing,Float64}}(undef, nrow(stations))
        for (index, station) in enumerate(eachrow(stations))
            lon_index = _local_grid_index(
                longitude, station.lon, GRID_LON_ORIGIN, GRID_NLON, spec.lon_first,
            )
            lat_index = _local_grid_index(
                latitude, station.lat, GRID_LAT_ORIGIN, GRID_NLAT, spec.lat_first,
            )
            value = precipitation[lat_index, lon_index, 1]
            if ismissing(value)
                values[index] = missing
            else
                numeric = Float64(value)
                isfinite(numeric) || error("non-finite IMERG value at station $(station.station_id)")
                numeric >= 0 || error("negative IMERG value at station $(station.station_id)")
                values[index] = numeric
            end
        end
        values
    end
end

"""Convert two 30-minute rates (mm/hr) to one hourly accumulation (mm)."""
function aggregate_hour_values(first_half, second_half; allow_missing::Bool=false)
    length(first_half) == length(second_half) ||
        throw(ArgumentError("half-hour station vectors have different lengths"))
    result = Vector{Union{Missing,Float64}}(undef, length(first_half))
    for index in eachindex(first_half, second_half)
        first_value = first_half[index]
        second_value = second_half[index]
        if ismissing(first_value) || ismissing(second_value)
            allow_missing || error(
                "missing station value in one or both half-hours at vector index $index",
            )
            result[index] = missing
        else
            result[index] = 0.5 * Float64(first_value) + 0.5 * Float64(second_value)
        end
    end
    return result
end

function _response_status(error_value)
    error_value isa Downloads.RequestError || return 0
    response = error_value.response
    response === nothing && return 0
    return response.status
end

# Downloads.RequestError can carry HTTP 200 when libcurl aborts an incomplete/too-slow
# transfer after the response headers arrived. That is still a transport failure and must
# be retried; a genuinely complete HTTP 200 never reaches this catch path.
_retryable_status(status::Integer) =
    status == 0 || status == 200 || status == 408 || status == 429 || status >= 500

function _raw_subset_path(raw_dir::AbstractString, granule::GranuleRef)
    year = Dates.format(granule.start_time, dateformat"yyyy")
    day_of_year = lpad(string(Dates.dayofyear(Date(granule.start_time))), 3, '0')
    filename = replace(granule.producer_id, r"\.HDF5$" => ".subset.nc4")
    return joinpath(raw_dir, year, day_of_year, filename)
end

function _download_subset(
    granule::GranuleRef,
    spec::SubsetSpec;
    raw_dir::AbstractString,
    staging_dir::AbstractString,
    token::AbstractString,
    max_attempts::Integer=4,
)
    destination = _raw_subset_path(raw_dir, granule)
    if isfile(destination)
        validate_subset_file(destination, spec; expected_time=granule.start_time)
        return destination, "cached"
    end

    mkpath(dirname(destination))
    mkpath(staging_dir)
    staging_path = joinpath(staging_dir, string(basename(destination), ".part"))
    url = build_subset_url(granule, spec)
    for attempt in 1:max_attempts
        isfile(staging_path) && rm(staging_path; force=true)
        try
            Downloads.download(
                url,
                staging_path;
                headers=["Authorization" => "Bearer $token"],
                timeout=300,
            )
            validate_subset_file(staging_path, spec; expected_time=granule.start_time)
            isfile(destination) && error(
                "raw destination appeared during download and will not be overwritten: $destination",
            )
            mv(staging_path, destination)
            return destination, "downloaded"
        catch error_value
            isfile(staging_path) && rm(staging_path; force=true)
            status = _response_status(error_value)
            if attempt == max_attempts || !_retryable_status(status)
                rethrow(error_value)
            end
            delay = min(30, 5 * 2^(attempt - 1))
            println(stderr,
                "IMERG request failed (HTTP $status), attempt $attempt/$max_attempts; " *
                "retrying in $delay seconds: $(granule.producer_id)",
            )
            sleep(delay)
        end
    end
    error("unreachable IMERG download state")
end

function _write_atomic(writer::Function, path::AbstractString; overwrite::Bool=false)
    mkpath(dirname(path))
    temporary = string(path, ".tmp")
    isfile(temporary) && rm(temporary; force=true)
    try
        writer(temporary)
        if isfile(path) && !overwrite
            error("refusing to overwrite existing file: $path")
        end
        mv(temporary, path; force=overwrite)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

function _write_catalog(path::AbstractString, granules, raw_dir::AbstractString)
    table = DataFrame(
        time=[Dates.format(item.start_time, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"
              for item in granules],
        producer_id=[item.producer_id for item in granules],
        opendap_url=[item.opendap_url for item in granules],
        raw_subset=[_raw_subset_path(raw_dir, item) for item in granules],
    )
    return _write_atomic(path; overwrite=true) do temporary
        CSV.write(temporary, table)
    end
end

function _parse_utc_hour(value)
    text = replace(strip(string(value)), "Z" => "")
    return DateTime(text, dateformat"yyyy-mm-ddTHH:MM:SS")
end

function validate_month_csv(
    path::AbstractString,
    stations::DataFrame,
    month_start::DateTime,
    month_stop::DateTime,
)
    isfile(path) || error("monthly CSV does not exist: $path")
    table = CSV.read(
        path, DataFrame;
        select=["time", "station_id", "gpm_mm_h"],
        types=Dict(:station_id => String),
    )
    hour_count = expected_hours(month_start, month_stop)
    expected_rows = hour_count * nrow(stations)
    nrow(table) == expected_rows || error(
        "$(basename(path)) has $(nrow(table)) rows; expected $expected_rows",
    )
    any(ismissing, table.gpm_mm_h) &&
        error("$(basename(path)) contains missing gpm_mm_h values")
    values = Float64.(table.gpm_mm_h)
    all(isfinite, values) || error("$(basename(path)) contains non-finite values")
    all(>=(0), values) || error("$(basename(path)) contains negative precipitation")

    station_ids = stations.station_id
    timestamps = _parse_utc_hour.(table.time)
    for hour_index in 0:(hour_count - 1)
        rows = (hour_index * nrow(stations) + 1):((hour_index + 1) * nrow(stations))
        expected_time = month_start + Hour(hour_index)
        all(==(expected_time), @view timestamps[rows]) ||
            error("$(basename(path)) has a missing, duplicate, or out-of-order hour")
        table.station_id[rows] == station_ids ||
            error("$(basename(path)) station order or membership is invalid")
    end
    return (
        rows=nrow(table),
        hours=hour_count,
        stations=nrow(stations),
        missing_values=0,
        minimum=minimum(values),
        maximum=maximum(values),
        sha256=bytes2hex(open(SHA.sha256, path)),
    )
end

function _month_table(granules, paths, stations, spec; allow_missing::Bool=false)
    length(granules) == length(paths) || error("granule/path count mismatch")
    iseven(length(granules)) || error("monthly granule count is not even")
    station_count = nrow(stations)
    hour_count = length(granules) ÷ 2
    table = DataFrame(
        time=Vector{String}(undef, hour_count * station_count),
        station_id=Vector{String}(undef, hour_count * station_count),
        gpm_mm_h=Vector{Union{Missing,Float64}}(undef, hour_count * station_count),
    )
    for hour_index in 1:hour_count
        first_index = 2 * hour_index - 1
        second_index = first_index + 1
        first_granule = granules[first_index]
        second_granule = granules[second_index]
        minute(first_granule.start_time) == 0 ||
            error("first granule of hour is not at minute 00")
        second_granule.start_time == first_granule.start_time + HALF_HOUR ||
            error("hour does not contain exactly the 00 and 30 minute granules")
        first_values = read_station_values(
            paths[first_index], stations, spec; expected_time=first_granule.start_time,
        )
        second_values = read_station_values(
            paths[second_index], stations, spec; expected_time=second_granule.start_time,
        )
        hourly = aggregate_hour_values(first_values, second_values; allow_missing)
        rows = ((hour_index - 1) * station_count + 1):(hour_index * station_count)
        table.time[rows] .= Dates.format(
            first_granule.start_time, dateformat"yyyy-mm-ddTHH:MM:SS",
        ) * "Z"
        table.station_id[rows] .= stations.station_id
        table.gpm_mm_h[rows] .= hourly
    end
    return table
end

function _download_catalog(
    granules,
    spec;
    raw_dir,
    staging_dir,
    token,
    concurrency,
    max_attempts,
)
    concurrency > 0 || throw(ArgumentError("concurrency must be positive"))
    completed = Ref(0)
    statuses = asyncmap(granules; ntasks=concurrency) do granule
        result = _download_subset(
            granule,
            spec;
            raw_dir,
            staging_dir,
            token,
            max_attempts,
        )
        completed[] += 1
        if completed[] == 1 || completed[] % 50 == 0 || completed[] == length(granules)
            println("  subsets: $(completed[])/$(length(granules))")
        end
        result
    end
    return first.(statuses), last.(statuses)
end

function _write_audit_json(path, audit)
    return _write_atomic(path; overwrite=true) do temporary
        open(temporary, "w") do io
            JSON.print(io, audit, 2)
            println(io)
        end
    end
end

"""Download, validate, aggregate, and audit one complete UTC calendar month."""
function process_imerg_month(
    month_start::DateTime,
    stations::DataFrame;
    raw_dir::AbstractString,
    output_dir::AbstractString,
    audit_dir::AbstractString,
    staging_dir::AbstractString,
    token::AbstractString,
    concurrency::Integer=4,
    max_attempts::Integer=4,
    overwrite_output::Bool=false,
    catalog::Union{Nothing,AbstractVector{GranuleRef}}=nothing,
)
    day(month_start) == 1 && hour(month_start) == 0 && minute(month_start) == 0 ||
        throw(ArgumentError("month_start must be the first day at 00:00 UTC"))
    month_stop = month_start + Month(1)
    suffix = Dates.format(month_start, dateformat"yyyymm")
    output_path = joinpath(output_dir, "gpm_hubei_hourly_long_$suffix.csv")
    if isfile(output_path) && !overwrite_output
        summary = validate_month_csv(output_path, stations, month_start, month_stop)
        println("Already complete and valid: $output_path")
        return output_path, summary
    end

    granules = catalog === nothing ? fetch_cmr_catalog(month_start, month_stop) : catalog
    validate_catalog(granules, month_start, month_stop)
    spec = build_subset_spec(stations)
    _write_catalog(
        joinpath(audit_dir, "catalog_$suffix.csv"), granules, raw_dir,
    )
    paths, statuses = _download_catalog(
        granules,
        spec;
        raw_dir,
        staging_dir,
        token,
        concurrency,
        max_attempts,
    )
    table = _month_table(granules, paths, stations, spec)
    _write_atomic(output_path; overwrite=overwrite_output) do temporary
        CSV.write(temporary, table)
    end
    summary = validate_month_csv(output_path, stations, month_start, month_stop)
    audit = Dict(
        "month" => suffix,
        "source_product" => "$IMERG_SHORT_NAME.$IMERG_VERSION",
        "source_variable" => "/Grid/Intermediate/$IMERG_VARIABLE",
        "source_units" => IMERG_UNITS,
        "aggregation" => "0.5 * rate_at_minute_00 + 0.5 * rate_at_minute_30",
        "time_semantics" => "UTC hour start",
        "granules" => length(granules),
        "expected_granules" => expected_half_hours(month_start, month_stop),
        "downloaded_subsets" => count(==("downloaded"), statuses),
        "cached_subsets" => count(==("cached"), statuses),
        "rows" => summary.rows,
        "hours" => summary.hours,
        "stations" => summary.stations,
        "missing_values" => summary.missing_values,
        "minimum_gpm_mm_h" => summary.minimum,
        "maximum_gpm_mm_h" => summary.maximum,
        "csv_sha256" => summary.sha256,
    )
    _write_audit_json(joinpath(audit_dir, "qc_$suffix.json"), audit)
    return output_path, summary
end

"""Download two temporary subsets and write one validated 318-row hourly test CSV."""
function validate_remote_hour(
    hour_start::DateTime,
    stations::DataFrame;
    token::AbstractString,
    audit_dir::AbstractString,
    concurrency::Integer=2,
    max_attempts::Integer=4,
)
    minute(hour_start) == 0 && second(hour_start) == 0 ||
        throw(ArgumentError("validation hour must start at HH:00:00"))
    granules = fetch_cmr_catalog(hour_start, hour_start + Hour(1))
    spec = build_subset_spec(stations)
    mkpath(audit_dir)
    return mktempdir() do temporary_raw
        paths, statuses = _download_catalog(
            granules,
            spec;
            raw_dir=temporary_raw,
            staging_dir=joinpath(temporary_raw, "staging"),
            token,
            concurrency,
            max_attempts,
        )
        table = _month_table(granules, paths, stations, spec)
        suffix = Dates.format(hour_start, dateformat"yyyymmdd_HHMM")
        output_path = joinpath(audit_dir, "validation_hour_$suffix.csv")
        _write_atomic(output_path; overwrite=true) do temporary
            CSV.write(temporary, table)
        end
        values = Float64.(table.gpm_mm_h)
        audit = Dict(
            "validated" => true,
            "hour_start_utc" => Dates.format(
                hour_start, dateformat"yyyy-mm-ddTHH:MM:SS",
            ) * "Z",
            "granules" => length(granules),
            "downloaded_subsets" => count(==("downloaded"), statuses),
            "stations" => nrow(stations),
            "rows" => nrow(table),
            "missing_values" => count(ismissing, table.gpm_mm_h),
            "minimum_gpm_mm_h" => minimum(values),
            "maximum_gpm_mm_h" => maximum(values),
            "subset_lon_indices" => [spec.lon_first, spec.lon_last],
            "subset_lat_indices" => [spec.lat_first, spec.lat_last],
            "variable" => "/Grid/Intermediate/$IMERG_VARIABLE",
            "source_units" => IMERG_UNITS,
            "output_units" => "mm per hourly interval",
            "csv_sha256" => bytes2hex(open(SHA.sha256, output_path)),
        )
        audit_path = joinpath(audit_dir, "validation_hour_$suffix.json")
        _write_audit_json(audit_path, audit)
        return output_path, audit_path, audit
    end
end

end
