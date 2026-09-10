module ERA5LandCovariates

using CSV
using DataFrames
using Dates

export parse_era5_datetime, prepare_era5_annual_covariates

const OUTPUT_COLUMNS = [
    :station_id, :time_utc, :time_bjt,
    :t2m_c, :d2m_c, :relative_humidity,
    :u10, :v10, :wind_speed, :wind_direction,
    :wind_direction_sin, :wind_direction_cos,
    :sp_hpa, :sp_anomaly_hpa,
]

function parse_era5_datetime(value)
    value isa DateTime && return value
    text = replace(strip(string(value)), r"\.\d+$" => "")
    return DateTime(text, dateformat"yyyy-mm-ddTHH:MM:SS")
end

function validate_station_table(table::DataFrame, station_id::String)
    required = [
        :station_id, :time_utc, :time_bjt, :t2m_c, :d2m_c,
        :relative_humidity, :u10, :v10, :wind_speed, :wind_direction,
        :sp_hpa, :sp_anomaly_hpa,
    ]
    missing_columns = setdiff(required, propertynames(table))
    isempty(missing_columns) || error("$station_id missing columns: $(join(missing_columns, ", "))")
    nrow(table) == 26_304 || error("$station_id has $(nrow(table)) rows; expected 26304")
    all(string.(table.station_id) .== station_id) || error("Station ID mismatch in $station_id")

    utc = parse_era5_datetime.(table.time_utc)
    bjt = parse_era5_datetime.(table.time_bjt)
    all(bjt .== utc .+ Hour(8)) || error("UTC/BJT mismatch in $station_id")
    allunique(utc) || error("Duplicate UTC timestamps in $station_id")
    issorted(utc) || error("Unsorted UTC timestamps in $station_id")
    return utc, bjt
end

function prepare_partition(table::DataFrame, utc, bjt, year_value::Int)
    selected = year.(utc) .== year_value
    part = table[selected, [
        :station_id, :t2m_c, :d2m_c, :relative_humidity,
        :u10, :v10, :wind_speed, :wind_direction, :sp_hpa, :sp_anomaly_hpa,
    ]]
    directions = Float64.(part.wind_direction)
    insertcols!(part, 2,
        :time_utc => Dates.format.(utc[selected], dateformat"yyyy-mm-ddTHH:MM:SS"),
        :time_bjt => Dates.format.(bjt[selected], dateformat"yyyy-mm-ddTHH:MM:SS"),
    )
    insertcols!(part, 11,
        :wind_direction_sin => sind.(directions),
        :wind_direction_cos => cosd.(directions),
    )
    select!(part, OUTPUT_COLUMNS)
    return part
end

function prepare_era5_annual_covariates(
    input_dir::AbstractString,
    output_dir::AbstractString,
    station_meta_path::AbstractString,
    audit_path::AbstractString;
    years=2022:2024,
)
    stations = CSV.read(station_meta_path, DataFrame; types=Dict(:station_id => String))
    station_ids = sort(String.(stations.station_id))
    length(station_ids) == 237 || error("Expected 237 stations, found $(length(station_ids))")
    allunique(station_ids) || error("Duplicate station IDs in metadata")

    mkpath(output_dir)
    outputs = Dict(year_value => joinpath(output_dir, "era5_land_station_hourly_utc_$year_value.csv") for year_value in years)
    temporary = Dict(year_value => string(outputs[year_value], ".tmp-", getpid()) for year_value in years)
    foreach(path -> isfile(path) && rm(path; force=true), values(temporary))
    written = Dict(year_value => false for year_value in years)

    try
        for (index, station_id) in enumerate(station_ids)
            path = joinpath(input_dir, "$station_id.csv")
            isfile(path) || error("Missing ERA5 station file: $path")
            table = CSV.read(path, DataFrame; types=Dict(:station_id => String))
            utc, bjt = validate_station_table(table, station_id)
            for year_value in years
                part = prepare_partition(table, utc, bjt, year_value)
                expected = isleapyear(year_value) ? 8_784 : 8_760
                nrow(part) == expected || error("$station_id UTC $year_value has $(nrow(part)) rows")
                CSV.write(temporary[year_value], part; append=written[year_value], writeheader=!written[year_value])
                written[year_value] = true
            end
            index % 25 == 0 && println("Prepared ERA5 covariates for $index/$(length(station_ids)) stations")
        end
        for year_value in years
            mv(temporary[year_value], outputs[year_value]; force=true)
        end
    finally
        foreach(path -> isfile(path) && rm(path; force=true), values(temporary))
    end

    audit = DataFrame(
        utc_year=Int[], stations=Int[], hours_per_station=Int[], rows=Int[],
        first_utc=String[], last_utc=String[], first_bjt=String[], last_bjt=String[],
        output_path=String[],
    )
    for year_value in years
        hours = isleapyear(year_value) ? 8_784 : 8_760
        push!(audit, (
            year_value, length(station_ids), hours, length(station_ids) * hours,
            "$year_value-01-01T00:00:00", "$year_value-12-31T23:00:00",
            "$year_value-01-01T08:00:00", "$(year_value + 1)-01-01T07:00:00",
            outputs[year_value],
        ))
    end
    mkpath(dirname(audit_path))
    CSV.write(audit_path, audit)
    return (; outputs, audit)
end

end
