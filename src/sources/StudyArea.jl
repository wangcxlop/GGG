module StudyArea

using CSV, DataFrames

export StudyBounds, STUDY_BOUNDS, filter_stations, prepare_study_area_inputs

struct StudyBounds
    west::Float64
    east::Float64
    south::Float64
    north::Float64

    function StudyBounds(west::Real, east::Real, south::Real, north::Real)
        values = Float64.((west, east, south, north))
        all(isfinite, values) || throw(ArgumentError("study-area bounds must be finite"))
        values[1] < values[2] || throw(ArgumentError("west must be smaller than east"))
        values[3] < values[4] || throw(ArgumentError("south must be smaller than north"))
        return new(values...)
    end
end

const STUDY_BOUNDS = StudyBounds(109.4, 111.6, 31.2, 33.4)

function _write_csv_atomic(path::AbstractString, table)
    mkpath(dirname(path))
    temporary_path = string(path, ".tmp-", getpid())
    try
        CSV.write(temporary_path, table)
        mv(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
    return path
end

function filter_stations(metadata::DataFrame; bounds::StudyBounds=STUDY_BOUNDS)
    required = [:station_id, :lon, :lat]
    missing_columns = setdiff(required, propertynames(metadata))
    isempty(missing_columns) ||
        throw(ArgumentError("station metadata is missing columns: $(join(missing_columns, ", "))"))

    any(ismissing, metadata.station_id) && throw(ArgumentError("station_id contains missing values"))
    any(ismissing, metadata.lon) && throw(ArgumentError("lon contains missing values"))
    any(ismissing, metadata.lat) && throw(ArgumentError("lat contains missing values"))

    station_ids = string.(metadata.station_id)
    length(unique(station_ids)) == length(station_ids) ||
        throw(ArgumentError("station metadata contains duplicate station_id values"))

    lon = try
        Float64.(metadata.lon)
    catch
        throw(ArgumentError("lon must contain numeric values"))
    end
    lat = try
        Float64.(metadata.lat)
    catch
        throw(ArgumentError("lat must contain numeric values"))
    end
    all(isfinite, lon) || throw(ArgumentError("lon contains non-finite values"))
    all(isfinite, lat) || throw(ArgumentError("lat contains non-finite values"))

    inside = (bounds.west .<= lon .<= bounds.east) .&
        (bounds.south .<= lat .<= bounds.north)
    selected = copy(metadata[inside, :])
    isempty(selected) && throw(ArgumentError("no stations fall inside the study area"))
    selected.station_id = station_ids[inside]
    selected.lon = lon[inside]
    selected.lat = lat[inside]
    sort!(selected, :station_id)
    return selected
end

function _subset_wide_file(
    source_path::AbstractString,
    output_path::AbstractString,
    station_ids::Vector{String},
)
    table = CSV.read(source_path, DataFrame)
    isempty(names(table)) && throw(ArgumentError("wide table has no columns: $source_path"))
    time_column = first(names(table))
    available = Set(string.(names(table)[2:end]))
    missing_ids = [station_id for station_id in station_ids if !(station_id in available)]
    isempty(missing_ids) || throw(ArgumentError(
        "$(basename(source_path)) is missing $(length(missing_ids)) study-area stations: " *
        join(missing_ids, ", "),
    ))

    selected = select(table, time_column, station_ids...)
    _write_csv_atomic(output_path, selected)
    return (path=String(output_path), rows=nrow(selected), stations=length(station_ids))
end

function prepare_study_area_inputs(;
    station_meta_path::AbstractString,
    wide_sources::Dict{String,String},
    output_dir::AbstractString,
    differences_path::AbstractString,
    bounds::StudyBounds=STUDY_BOUNDS,
)
    metadata = CSV.read(station_meta_path, DataFrame)
    stations = filter_stations(metadata; bounds)
    station_ids = Vector{String}(stations.station_id)
    selected_set = Set(station_ids)

    differences = DataFrame(source=String[], station_id=String[], issue=String[])
    all_metadata_ids = Set(string.(metadata.station_id))
    for station_id in sort!(collect(setdiff(all_metadata_ids, selected_set)))
        push!(differences, ("station_meta", station_id, "excluded_outside_study_area"))
    end

    missing_by_source = Dict{String,Vector{String}}()
    for (source, path) in sort!(collect(wide_sources); by=first)
        header = CSV.File(path; limit=1)
        available = Set(string.(propertynames(header))[2:end])
        missing_ids = [station_id for station_id in station_ids if !(station_id in available)]
        missing_by_source[source] = missing_ids
        for station_id in missing_ids
            push!(differences, (source, station_id, "missing_required_station"))
        end
        for station_id in sort!(collect(setdiff(available, selected_set)))
            issue = station_id in all_metadata_ids ?
                "excluded_outside_study_area" : "not_in_station_metadata"
            push!(differences, (source, station_id, issue))
        end
    end
    _write_csv_atomic(differences_path, differences)

    failing_sources = [source for (source, ids) in missing_by_source if !isempty(ids)]
    isempty(failing_sources) || throw(ArgumentError(
        "study-area stations are missing from: $(join(sort(failing_sources), ", ")); " *
        "see $differences_path",
    ))

    mkpath(output_dir)
    station_output = joinpath(output_dir, "station_meta.csv")
    _write_csv_atomic(station_output, stations)
    results = Dict{String,Any}()
    for (source, path) in sort!(collect(wide_sources); by=first)
        output_path = joinpath(output_dir, basename(path))
        results[source] = _subset_wide_file(path, output_path, station_ids)
    end

    return (;
        station_output,
        station_count=length(station_ids),
        raw_station_count=nrow(metadata),
        actual_lon_min=minimum(stations.lon),
        actual_lon_max=maximum(stations.lon),
        actual_lat_min=minimum(stations.lat),
        actual_lat_max=maximum(stations.lat),
        results,
    )
end

end # module
