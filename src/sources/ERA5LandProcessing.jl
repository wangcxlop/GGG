module ERA5LandProcessing

using CSV, DataFrames, Dates, Statistics

export relative_humidity, wind_speed, wind_direction,
    process_era5_station_zip, process_era5_station_directory

const REQUIRED_VARIABLES = [:d2m, :t2m, :sp, :u10, :v10]
const EXPECTED_HOURLY_ROWS = 26_304

"""Relative humidity (%) from 2 m temperature and dewpoint in Kelvin."""
function relative_humidity(t2m_k::Real, d2m_k::Real)
    temperature_c = Float64(t2m_k) - 273.15
    dewpoint_c = Float64(d2m_k) - 273.15
    ratio = exp(
        17.625 * dewpoint_c / (243.04 + dewpoint_c) -
        17.625 * temperature_c / (243.04 + temperature_c),
    )
    return clamp(100 * ratio, 0.0, 100.0)
end

wind_speed(u10::Real, v10::Real) = hypot(Float64(u10), Float64(v10))

"""Meteorological wind direction in degrees clockwise from north."""
function wind_direction(u10::Real, v10::Real)
    return mod(rad2deg(atan(-Float64(u10), -Float64(v10))), 360.0)
end

function _extract_zip(zip_path::AbstractString, destination::AbstractString)
    run(`tar -xf $zip_path -C $destination`)
    return destination
end

function _parse_time(value)
    value isa DateTime && return value
    return DateTime(String(value), dateformat"yyyy-mm-dd HH:MM:SS")
end

function _load_zip_tables(zip_path::AbstractString)
    mktempdir() do temp_dir
        _extract_zip(zip_path, temp_dir)
        csv_paths = filter(
            path -> endswith(lowercase(path), ".csv"),
            readdir(temp_dir; join=true),
        )
        length(csv_paths) == 3 || error(
            "expected three ERA5-Land CSV groups, found $(length(csv_paths))",
        )

        merged = nothing
        grid_coordinates = Tuple{Float64,Float64}[]
        for path in csv_paths
            table = CSV.read(path, DataFrame; stringtype=String)
            :valid_time in propertynames(table) || error("missing valid_time in $path")
            variables = intersect(REQUIRED_VARIABLES, propertynames(table))
            isempty(variables) && error("no required ERA5-Land variables in $path")
            all(name -> name in propertynames(table), [:latitude, :longitude]) ||
                error("missing ERA5 grid coordinates in $path")
            push!(grid_coordinates, (
                Float64(first(table.latitude)), Float64(first(table.longitude)),
            ))
            part = select(table, :valid_time, variables...)
            merged = merged === nothing ? part :
                outerjoin(merged, part; on=:valid_time, validate=(true, true))
        end

        all(coordinate -> isapprox(coordinate[1], grid_coordinates[1][1]) &&
            isapprox(coordinate[2], grid_coordinates[1][2]), grid_coordinates) ||
            error("ERA5 variable groups use inconsistent grid coordinates")
        return merged, grid_coordinates[1]
    end
end

function _pressure_anomaly!(table::DataFrame)
    table.month_bjt = Dates.format.(table.time_bjt, dateformat"yyyy-mm")
    transform!(
        groupby(table, :month_bjt),
        :sp_hpa => (values -> values .- mean(values)) => :sp_anomaly_hpa,
    )
    return table
end

"""Convert one raw station ZIP into a model-ready hourly CSV."""
function process_era5_station_zip(
    zip_path::AbstractString,
    output_path::AbstractString;
    station_id::AbstractString,
    station_lon::Real,
    station_lat::Real,
)
    isfile(zip_path) || throw(ArgumentError("ERA5 station ZIP not found: $zip_path"))
    table, (grid_lat, grid_lon) = _load_zip_tables(zip_path)
    all(variable -> variable in propertynames(table), REQUIRED_VARIABLES) ||
        error("ERA5 ZIP does not contain all five required variables")

    table.time_utc = _parse_time.(table.valid_time)
    sort!(table, :time_utc)
    nrow(table) == EXPECTED_HOURLY_ROWS || error(
        "expected $EXPECTED_HOURLY_ROWS hourly rows, found $(nrow(table))",
    )
    allunique(table.time_utc) || error("duplicate ERA5 hourly timestamps found")
    first(table.time_utc) == DateTime(2022, 1, 1) ||
        error("unexpected first ERA5 timestamp: $(first(table.time_utc))")
    last(table.time_utc) == DateTime(2024, 12, 31, 23) ||
        error("unexpected last ERA5 timestamp: $(last(table.time_utc))")

    output = DataFrame(
        station_id=fill(String(station_id), nrow(table)),
        station_lon=fill(Float64(station_lon), nrow(table)),
        station_lat=fill(Float64(station_lat), nrow(table)),
        era5_lon=fill(grid_lon, nrow(table)),
        era5_lat=fill(grid_lat, nrow(table)),
        time_utc=table.time_utc,
        time_bjt=table.time_utc .+ Hour(8),
        t2m_c=Float64.(table.t2m) .- 273.15,
        d2m_c=Float64.(table.d2m) .- 273.15,
        relative_humidity=relative_humidity.(table.t2m, table.d2m),
        u10=Float64.(table.u10),
        v10=Float64.(table.v10),
        wind_speed=wind_speed.(table.u10, table.v10),
        wind_direction=wind_direction.(table.u10, table.v10),
        sp_hpa=Float64.(table.sp) ./ 100,
    )
    _pressure_anomaly!(output)

    mkpath(dirname(output_path))
    temp_path = string(output_path, ".tmp-", getpid())
    CSV.write(temp_path, output)
    mv(temp_path, output_path; force=true)
    return output_path
end

function process_era5_station_directory(
    station_dir::AbstractString,
    output_path::AbstractString;
    station_id::AbstractString,
    station_lon::Real,
    station_lat::Real,
)
    isfile(joinpath(station_dir, ".complete")) ||
        error("station download is not marked complete: $station_dir")
    zip_paths = filter(
        path -> endswith(lowercase(path), ".zip"),
        readdir(station_dir; join=true),
    )
    length(zip_paths) == 1 || error(
        "expected one station ZIP in $station_dir, found $(length(zip_paths))",
    )
    return process_era5_station_zip(
        only(zip_paths), output_path;
        station_id, station_lon, station_lat,
    )
end

end
