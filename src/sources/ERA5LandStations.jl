module ERA5LandStations

using CSV, DataFrames, Dates, JSON

export ERA5_LAND_DATASET, ERA5_LAND_VARIABLES,
    read_era5_stations, build_era5_land_request,
    write_station_requests, download_era5_land_stations

const ERA5_LAND_DATASET = "reanalysis-era5-land-timeseries"
const ERA5_LAND_VARIABLES = [
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_pressure",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
]
const START_DATE = Date(2022, 1, 1)
const END_DATE = Date(2024, 12, 31)

const PYTHON_DRIVER = raw"""
import json
import os
import pathlib
import sys

import cdsapi

request_path = pathlib.Path(sys.argv[1]).resolve()
output_dir = pathlib.Path(sys.argv[2]).resolve()
output_dir.mkdir(parents=True, exist_ok=True)

with request_path.open("r", encoding="utf-8") as stream:
    request = json.load(stream)

token = os.environ.get("CDSAPI_TOKEN", "").strip()
if not token:
    raise RuntimeError("CDSAPI_TOKEN is empty")

client = cdsapi.Client(
    url="https://cds.climate.copernicus.eu/api",
    key=token,
    quiet=False,
)

previous_dir = pathlib.Path.cwd()
try:
    os.chdir(output_dir)
    result = client.retrieve("reanalysis-era5-land-timeseries", request).download()
    print(f"CDS download result: {result}")
finally:
    os.chdir(previous_dir)
"""

function _column(df::DataFrame, requested::Symbol)
    mapping = Dict(Symbol(lowercase(String(name))) => name for name in names(df))
    column = get(mapping, Symbol(lowercase(String(requested))), nothing)
    column === nothing && throw(ArgumentError("missing column: $requested"))
    return column
end

function _safe_station_id(value)
    station_id = strip(string(value))
    isempty(station_id) && throw(ArgumentError("station id cannot be empty"))
    occursin(r"^[A-Za-z0-9_-]+$", station_id) ||
        throw(ArgumentError("unsafe station id: $station_id"))
    return station_id
end

"""Read and validate the 237 study-area stations used by this project."""
function read_era5_stations(
    path::AbstractString;
    expected_count::Integer=237,
    west::Real=109.4,
    east::Real=111.6,
    south::Real=31.2,
    north::Real=33.4,
)
    isfile(path) || throw(ArgumentError("station metadata not found: $path"))
    df = CSV.read(
        path, DataFrame;
        stringtype=String,
        types=Dict(:station_id => String),
    )
    id_column = _column(df, :station_id)
    lon_column = _column(df, :lon)
    lat_column = _column(df, :lat)

    stations = DataFrame(
        station_id=[_safe_station_id(value) for value in df[!, id_column]],
        lon=Float64.(df[!, lon_column]),
        lat=Float64.(df[!, lat_column]),
    )
    nrow(stations) == expected_count || throw(ArgumentError(
        "expected $expected_count stations, found $(nrow(stations))",
    ))
    allunique(stations.station_id) || throw(ArgumentError("duplicate station ids found"))

    outside = filter(row -> !(
        west <= row.lon <= east && south <= row.lat <= north
    ), stations)
    isempty(outside) || throw(ArgumentError(
        "$(nrow(outside)) stations fall outside $west-$east E, $south-$north N",
    ))
    return stations
end

function build_era5_land_request(
    longitude::Real,
    latitude::Real;
    start_date::Date=START_DATE,
    end_date::Date=END_DATE,
)
    start_date <= end_date || throw(ArgumentError("start_date must not exceed end_date"))
    return Dict(
        "variable" => copy(ERA5_LAND_VARIABLES),
        "location" => Dict(
            "longitude" => Float64(longitude),
            "latitude" => Float64(latitude),
        ),
        "date" => [
            string(Dates.format(start_date, dateformat"yyyy-mm-dd"), "/",
                Dates.format(end_date, dateformat"yyyy-mm-dd")),
        ],
        "data_format" => "csv",
    )
end

function _write_json(path::AbstractString, value)
    mkpath(dirname(path))
    temp_path = string(path, ".tmp-", getpid())
    open(temp_path, "w") do io
        JSON.print(io, value, 2)
        println(io)
    end
    mv(temp_path, path; force=true)
    return path
end

function _manifest(stations::DataFrame, raw_dir::AbstractString)
    return DataFrame(
        station_id=stations.station_id,
        lon=stations.lon,
        lat=stations.lat,
        status=fill("pending", nrow(stations)),
        raw_directory=[joinpath(raw_dir, id) for id in stations.station_id],
    )
end

"""Write one auditable CDS request JSON per station and return the manifest."""
function write_station_requests(
    stations::DataFrame,
    audit_dir::AbstractString,
    raw_dir::AbstractString,
)
    request_dir = joinpath(audit_dir, "requests")
    mkpath(request_dir)
    for row in eachrow(stations)
        request = build_era5_land_request(row.lon, row.lat)
        _write_json(joinpath(request_dir, string(row.station_id, ".json")), request)
    end
    manifest = _manifest(stations, raw_dir)
    manifest_path = joinpath(audit_dir, "era5_land_station_manifest.csv")
    mkpath(dirname(manifest_path))
    CSV.write(manifest_path, manifest)
    return manifest, manifest_path
end

function _downloaded_assets(station_dir::AbstractString)
    isdir(station_dir) || return String[]
    assets = String[]
    for (directory, _, files) in walkdir(station_dir)
        for file in files
            file == ".complete" && continue
            path = joinpath(directory, file)
            filesize(path) > 0 && push!(assets, path)
        end
    end
    return assets
end

function _write_complete_marker(path::AbstractString)
    open(path, "w") do io
        println(io, Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    end
end

function _run_station_download(
    python::AbstractString,
    driver_path::AbstractString,
    request_path::AbstractString,
    station_dir::AbstractString,
    token::AbstractString;
    max_attempts::Integer=3,
)
    for attempt in 1:max_attempts
        try
            command = addenv(
                `$python $driver_path $request_path $station_dir`,
                "CDSAPI_TOKEN" => token,
            )
            run(command)
            isempty(_downloaded_assets(station_dir)) &&
                error("CDS returned no non-empty files")
            _write_complete_marker(joinpath(station_dir, ".complete"))
            return true
        catch error_value
            attempt == max_attempts && rethrow(error_value)
            println(stderr,
                "Station download attempt $attempt/$max_attempts failed; " *
                "retrying in 30 seconds.",
            )
            sleep(30)
        end
    end
    return false
end

"""Download 2022--2024 hourly ERA5-Land data for each station, with resume markers."""
function download_era5_land_stations(
    stations::DataFrame;
    raw_dir::AbstractString,
    audit_dir::AbstractString,
    token::AbstractString,
    python::AbstractString="python",
    limit::Union{Nothing,Integer}=nothing,
)
    isempty(strip(token)) && error("CDS API token is empty")
    selected = limit === nothing ? stations : first(stations, min(limit, nrow(stations)))
    manifest, manifest_path = write_station_requests(selected, audit_dir, raw_dir)

    mktemp() do driver_path, driver_io
        write(driver_io, PYTHON_DRIVER)
        close(driver_io)

        for (index, row) in enumerate(eachrow(selected))
            station_dir = joinpath(raw_dir, row.station_id)
            marker_path = joinpath(station_dir, ".complete")
            if isfile(marker_path) && !isempty(_downloaded_assets(station_dir))
                manifest.status[index] = "complete"
                println("[$index/$(nrow(selected))] Already complete: $(row.station_id)")
                continue
            end

            mkpath(station_dir)
            request_path = joinpath(audit_dir, "requests", string(row.station_id, ".json"))
            println("[$index/$(nrow(selected))] Downloading station $(row.station_id)")
            _run_station_download(
                python, driver_path, request_path, station_dir, token,
            )
            manifest.status[index] = "complete"
            CSV.write(manifest_path, manifest)
        end
    end
    CSV.write(manifest_path, manifest)
    return manifest
end

end
