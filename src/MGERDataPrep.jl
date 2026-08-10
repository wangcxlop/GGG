module MGERDataPrep

using CSV, DataFrames, Dates, Statistics

export audit_mger_inputs, parse_time_utcish, prepare_satellite_wide

function parse_time_utcish(value)
    text = replace(strip(String(value)), "Z" => "")
    for format in (
        dateformat"yyyy-mm-ddTHH:MM:SS",
        dateformat"yyyy-mm-dd HH:MM:SS",
        dateformat"yyyy/mm/dd HH:MM:SS",
        dateformat"yyyy-mm-ddTHH:MM",
        dateformat"yyyy-mm-dd HH:MM",
    )
        try
            return DateTime(text, format)
        catch
        end
    end
    error("Cannot parse timestamp: $value")
end

function write_csv_atomic(path::AbstractString, table)
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

function monthly_files(input_dir::AbstractString, product::AbstractString; years=2022:2025, months=6:9)
    prefix = lowercase(product)
    pattern = Regex("^$(prefix)_hubei_hourly_long_(\\d{4})(\\d{2})\\.csv\$")
    selected = String[]
    for file in readdir(input_dir; join=true)
        match_value = match(pattern, lowercase(basename(file)))
        match_value === nothing && continue
        year_value = parse(Int, match_value.captures[1])
        month_value = parse(Int, match_value.captures[2])
        year_value in years && month_value in months && push!(selected, file)
    end
    sort!(selected)
    isempty(selected) && error("No monthly $(product) CSV files found in $input_dir")
    return selected
end

function prepare_satellite_wide(
    product::AbstractString,
    input_dir::AbstractString,
    output_file::AbstractString;
    value_column::Symbol,
    shift_hours::Int=9,
    years=2022:2025,
    months=6:9,
    duplicate_qc_file::Union{Nothing, AbstractString}=nothing,
)
    files = monthly_files(input_dir, product; years, months)
    tables = DataFrame[]
    for file in files
        source = CSV.read(file, DataFrame; select=["time", "station_id", String(value_column)])
        rename!(source, value_column => :value)
        source.time = parse_time_utcish.(source.time) .+ Hour(shift_hours)
        source.station_id = string.(source.station_id)
        source.value = Float64.(source.value)
        push!(tables, source)
    end

    long = vcat(tables...)
    grouped = combine(
        groupby(long, [:time, :station_id]),
        :value => mean => :value,
        nrow => :source_rows,
        :value => (values -> length(unique(values))) => :unique_value_count,
    )
    duplicate_qc = grouped[grouped.source_rows .> 1, :]
    if nrow(duplicate_qc) > 0
        all(duplicate_qc.unique_value_count .<= 1) ||
            error("$product contains duplicate station-time rows with conflicting values")
    end
    duplicate_qc_file === nothing || write_csv_atomic(duplicate_qc_file, duplicate_qc)

    wide = unstack(select(grouped, :time, :station_id, :value), :time, :station_id, :value)
    sort!(wide, :time)
    station_columns = sort(names(wide, Not(:time)); by=string)
    select!(wide, :time, station_columns...)
    wide.time = Dates.format.(wide.time, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z"
    write_csv_atomic(output_file, wide)

    return (
        product=String(product),
        output_file=String(output_file),
        source_file_count=length(files),
        source_row_count=nrow(long),
        output_time_count=nrow(wide),
        station_count=length(station_columns),
        duplicate_pair_count=nrow(duplicate_qc),
    )
end

function read_wide_signature(path::AbstractString)
    table = CSV.read(path, DataFrame)
    time_column = Symbol(first(names(table)))
    times = parse_time_utcish.(table[!, time_column])
    station_ids = string.(names(table, Not(time_column)))
    return (; table, times, station_ids)
end

function audit_mger_inputs(;
    station_meta_path::AbstractString,
    observation_path::AbstractString,
    satellite_paths::Dict{String, String},
    output_dir::AbstractString,
    analysis_start::DateTime=DateTime(2022, 6, 1, 9),
    analysis_end::DateTime=DateTime(2024, 10, 1, 8),
)
    mkpath(output_dir)
    metadata = CSV.read(station_meta_path, DataFrame)
    required_columns = [:station_id, :lon, :lat]
    all(column -> column in propertynames(metadata), required_columns) ||
        error("Station metadata must contain station_id, lon, and lat")
    metadata.station_id = string.(metadata.station_id)
    nrow(unique(metadata, :station_id)) == nrow(metadata) ||
        error("Station metadata contains duplicate station_id values")

    observation = read_wide_signature(observation_path)
    metadata_ids = Set(metadata.station_id)
    observation_ids = Set(observation.station_ids)
    common_ids = intersect(metadata_ids, observation_ids)
    common_times = Set(filter(time -> analysis_start <= time <= analysis_end, observation.times))

    source_rows = NamedTuple[]
    push!(source_rows, (
        source="station_meta", rows=nrow(metadata), station_count=length(metadata_ids),
        time_count=missing, start_time=missing, end_time=missing,
    ))
    push!(source_rows, (
        source="observation", rows=nrow(observation.table), station_count=length(observation_ids),
        time_count=length(observation.times), start_time=string(first(observation.times)),
        end_time=string(last(observation.times)),
    ))

    station_rows = NamedTuple[]
    all_sources = Dict("station_meta" => metadata_ids, "observation" => observation_ids)
    for product in sort(collect(keys(satellite_paths)))
        signature = read_wide_signature(satellite_paths[product])
        ids = Set(signature.station_ids)
        times = Set(filter(time -> analysis_start <= time <= analysis_end, signature.times))
        all_sources[product] = ids
        common_ids = intersect(common_ids, ids)
        common_times = intersect(common_times, times)
        push!(source_rows, (
            source=product, rows=nrow(signature.table), station_count=length(ids),
            time_count=length(signature.times), start_time=string(first(signature.times)),
            end_time=string(last(signature.times)),
        ))
    end

    union_ids = union(values(all_sources)...)
    for station_id in sort(collect(union_ids))
        for source in sort(collect(keys(all_sources)))
            station_id in all_sources[source] && continue
            push!(station_rows, (station_id=station_id, missing_from=source))
        end
    end

    summary = DataFrame(source_rows)
    summary[!, :analysis_start] = fill(string(analysis_start), nrow(summary))
    summary[!, :analysis_end] = fill(string(analysis_end), nrow(summary))
    summary[!, :global_common_station_count] = fill(length(common_ids), nrow(summary))
    summary[!, :global_common_timestamp_count] = fill(length(common_times), nrow(summary))
    write_csv_atomic(joinpath(output_dir, "input_audit.csv"), summary)
    write_csv_atomic(joinpath(output_dir, "station_id_differences.csv"), DataFrame(station_rows))
    write_csv_atomic(
        joinpath(output_dir, "global_common_stations.csv"),
        DataFrame(station_id=sort(collect(common_ids))),
    )
    return (; summary, common_ids, common_times)
end

end # module
