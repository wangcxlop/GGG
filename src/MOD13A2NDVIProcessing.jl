module MOD13A2NDVIProcessing

using CSV
using DataFrames
using Dates

export parse_mod13a2_date, composite_observation_date, quality_class,
    land_water_class, clean_ndvi, process_mod13a2_table, prepare_mod13a2_ndvi

const NDVI_COLUMN = Symbol("MOD13A2_061__1_km_16_days_NDVI")
const COMPOSITE_DOY_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_composite_day_of_the_year")
const RELIABILITY_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_pixel_reliability")
const LAND_WATER_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_VI_Quality_Land/Water_Mask_Description")

const REQUIRED_COLUMNS = [
    :ID, :Latitude, :Longitude, :Date,
    NDVI_COLUMN, COMPOSITE_DOY_COLUMN, RELIABILITY_COLUMN, LAND_WATER_COLUMN,
]

parse_mod13a2_date(value::Date) = value
parse_mod13a2_date(value::DateTime) = Date(value)
parse_mod13a2_date(value) = Date(strip(string(value)))

function _integer_code(value)
    ismissing(value) && return nothing
    number = tryparse(Float64, strip(string(value)))
    number === nothing && return nothing
    isfinite(number) && isinteger(number) || return nothing
    return Int(number)
end

function composite_observation_date(composite_start, day_of_year)
    start_date = parse_mod13a2_date(composite_start)
    doy = _integer_code(day_of_year)
    doy === nothing && return missing
    max_doy = isleapyear(year(start_date)) ? 366 : 365
    1 <= doy <= max_doy || return missing
    return Date(year(start_date), 1, 1) + Day(doy - 1)
end

function quality_class(pixel_reliability)
    code = _integer_code(pixel_reliability)
    code == -1 && return "fill"
    code == 0 && return "good"
    code == 1 && return "marginal"
    code == 2 && return "snow_ice"
    code == 3 && return "cloudy"
    return "unknown"
end

function land_water_class(description)
    ismissing(description) && return "unknown"
    text = lowercase(strip(string(description)))
    occursin("nothing else but land", text) && return "land"
    occursin("shoreline", text) && return "shoreline"
    occursin("shallow inland water", text) && return "shallow_inland_water"
    occursin("water", text) && return "water"
    return "other"
end

function clean_ndvi(
    raw_value,
    pixel_reliability;
    land_class::Union{Nothing,AbstractString}=nothing,
    require_land::Bool=false,
)
    ismissing(raw_value) && return missing
    value = tryparse(Float64, strip(string(raw_value)))
    value === nothing && return missing
    isfinite(value) && -0.2 <= value <= 1.0 || return missing
    quality_class(pixel_reliability) in ("good", "marginal") || return missing
    require_land && land_class != "land" && return missing
    return value
end

function _audit_stations(table::DataFrame; low_coverage_threshold::Float64=0.8)
    rows = NamedTuple[]
    for station in groupby(table, :station_id)
        periods = nrow(station)
        quality_valid = count(!ismissing, station.ndvi_qc)
        land_valid = count(!ismissing, station.ndvi_land_qc)
        quality_values = collect(skipmissing(station.ndvi_qc))
        classes = station.quality_class
        surfaces = station.land_water_class
        quality_fraction = quality_valid / periods
        land_fraction = land_valid / periods
        has_nonland = any(!=("land"), surfaces)
        low_quality = quality_fraction < low_coverage_threshold
        low_land = land_fraction < low_coverage_threshold
        push!(rows, (
            station_id=first(station.station_id),
            lon=first(station.lon),
            lat=first(station.lat),
            periods,
            quality_valid_count=quality_valid,
            quality_valid_fraction=quality_fraction,
            land_valid_count=land_valid,
            land_valid_fraction=land_fraction,
            good_count=count(==("good"), classes),
            marginal_count=count(==("marginal"), classes),
            snow_ice_count=count(==("snow_ice"), classes),
            cloudy_count=count(==("cloudy"), classes),
            fill_count=count(==("fill"), classes),
            unknown_quality_count=count(==("unknown"), classes),
            land_count=count(==("land"), surfaces),
            shoreline_count=count(==("shoreline"), surfaces),
            shallow_inland_water_count=count(==("shallow_inland_water"), surfaces),
            other_surface_count=count(x -> x in ("water", "other", "unknown"), surfaces),
            ndvi_qc_min=isempty(quality_values) ? missing : minimum(quality_values),
            ndvi_qc_max=isempty(quality_values) ? missing : maximum(quality_values),
            has_nonland_pixel=has_nonland,
            low_quality_coverage=low_quality,
            low_land_coverage=low_land,
            review_required=has_nonland || low_quality || low_land,
        ))
    end
    return sort!(DataFrame(rows), :station_id)
end

function process_mod13a2_table(
    raw::DataFrame;
    start_date::Date=Date(2022, 1, 1),
    end_date::Date=Date(2024, 12, 31),
    expected_station_count::Int=237,
    expected_period_count::Int=69,
    expected_first_date::Date=Date(2022, 1, 1),
    expected_last_date::Date=Date(2024, 12, 18),
)
    missing_columns = setdiff(REQUIRED_COLUMNS, propertynames(raw))
    isempty(missing_columns) || error(
        "MOD13A2 input missing columns: $(join(string.(missing_columns), ", "))",
    )

    dates = parse_mod13a2_date.(raw.Date)
    selected = start_date .<= dates .<= end_date
    source = raw[selected, :]
    composite_start = dates[selected]

    result = DataFrame(
        station_id=strip.(string.(source.ID)),
        lon=Float64.(source.Longitude),
        lat=Float64.(source.Latitude),
        composite_start=composite_start,
        observation_date=[
            composite_observation_date(date, doy)
            for (date, doy) in zip(composite_start, source[!, COMPOSITE_DOY_COLUMN])
        ],
        ndvi_raw=Float64.(source[!, NDVI_COLUMN]),
        pixel_reliability=[
            something(_integer_code(value), typemin(Int))
            for value in source[!, RELIABILITY_COLUMN]
        ],
    )
    result.quality_class = quality_class.(result.pixel_reliability)
    result.land_water_class = land_water_class.(source[!, LAND_WATER_COLUMN])
    result.ndvi_qc = [
        clean_ndvi(value, reliability)
        for (value, reliability) in zip(result.ndvi_raw, result.pixel_reliability)
    ]
    result.ndvi_land_qc = [
        clean_ndvi(value, reliability; land_class=surface, require_land=true)
        for (value, reliability, surface) in
            zip(result.ndvi_raw, result.pixel_reliability, result.land_water_class)
    ]
    select!(result, [
        :station_id, :lon, :lat, :composite_start, :observation_date,
        :ndvi_raw, :ndvi_qc, :ndvi_land_qc, :pixel_reliability,
        :quality_class, :land_water_class,
    ])
    sort!(result, [:station_id, :composite_start])

    station_ids = unique(result.station_id)
    length(station_ids) == expected_station_count || error(
        "Expected $expected_station_count stations, found $(length(station_ids))",
    )
    expected_rows = expected_station_count * expected_period_count
    nrow(result) == expected_rows || error(
        "Expected $expected_rows station-period rows, found $(nrow(result))",
    )
    allunique(result[:, [:station_id, :composite_start]]) ||
        error("Duplicate station/composite-date rows")

    period_dates = sort(unique(result.composite_start))
    length(period_dates) == expected_period_count || error(
        "Expected $expected_period_count composite dates, found $(length(period_dates))",
    )
    extrema(period_dates) == (expected_first_date, expected_last_date) || error(
        "Unexpected composite date range: $(extrema(period_dates))",
    )
    for station in groupby(result, :station_id)
        station.composite_start == period_dates ||
            error("Incomplete composite-date sequence for $(first(station.station_id))")
        all(==(first(station.lon)), station.lon) ||
            error("Longitude changes for $(first(station.station_id))")
        all(==(first(station.lat)), station.lat) ||
            error("Latitude changes for $(first(station.station_id))")
    end
    for column in (:ndvi_qc, :ndvi_land_qc)
        all(value -> ismissing(value) || -0.2 <= value <= 1.0, result[!, column]) ||
            error("$column contains a value outside [-0.2, 1.0]")
    end

    audit = _audit_stations(result)
    return (; table=result, audit, period_dates)
end

function _write_csv_atomic(path::AbstractString, table)
    mkpath(dirname(path))
    temporary = string(path, ".tmp-", getpid())
    try
        CSV.write(temporary, table)
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

function prepare_mod13a2_ndvi(
    input_path::AbstractString,
    output_path::AbstractString,
    audit_path::AbstractString;
    kwargs...,
)
    isfile(input_path) || error("MOD13A2 input not found: $input_path")
    raw = CSV.read(input_path, DataFrame; types=Dict(:ID => String))
    result = process_mod13a2_table(raw; kwargs...)
    _write_csv_atomic(output_path, result.table)
    _write_csv_atomic(audit_path, result.audit)
    return merge(result, (; output_path=String(output_path), audit_path=String(audit_path)))
end

end
