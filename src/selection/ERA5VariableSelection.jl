module ERA5VariableSelection

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

# Shared bookkeeping for every variable-selection path; see src/selection/SelectionScaffolding.jl.
using Main.SelectionScaffolding

export ERA5SelectionConfig, ERA5_VARIABLES, align_feature_time, load_era5_panel
export balanced_spatial_folds, prepare_dynamic_panel, station_block_permutation
export dynamic_panel_screen, panel_spatial_variability_test, run_era5_variable_selection

const ERA5_VARIABLES = (:t2m_c, :d2m_c, :u10, :v10, :sp_hpa)

Base.@kwdef struct ERA5SelectionConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    min_stations_per_time::Int = 12
    k::Int = 5
    seed::Int = 20260816
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    association_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
    time_offset_hours::Int = 9
end

align_feature_time(time_utc::DateTime; offset_hours::Int=9) = time_utc + Hour(offset_hours)

function _parse_datetime(value)
    value isa DateTime && return value
    text = String(value)
    try
        return DateTime(text)
    catch
        for format in (dateformat"yyyy-mm-dd HH:MM:SS", dateformat"yyyy-mm-ddTHH:MM:SS")
            try
                return DateTime(text, format)
            catch
            end
        end
    end
    throw(ArgumentError("invalid ERA5 timestamp: $text"))
end

function _finite_float(value)
    ismissing(value) && return NaN
    value isa Real && return Float64(value)
    parsed = tryparse(Float64, String(value))
    return parsed === nothing ? NaN : parsed
end

function _csv_field(text::AbstractString)
    value = strip(text)
    if ncodeunits(value) >= 2 && startswith(value, '"') && endswith(value, '"')
        return replace(value[2:end-1], "\"\"" => "\"")
    end
    return value
end

function _bh_adjust(pvalues::AbstractVector{<:Real})
    p = Float64.(pvalues)
    m = length(p)
    m == 0 && return Float64[]
    order = sortperm(p)
    adjusted = fill(NaN, m)
    running = 1.0
    for rank in m:-1:1
        index = order[rank]
        running = min(running, p[index] * m / rank)
        adjusted[index] = min(running, 1.0)
    end
    return adjusted
end

"""Load requested station-hour ERA5 values using `time_utc + offset`, never `time_bjt`.

The function streams annual files and returns station × feature-time matrices plus a
quality-control table. Missing, duplicate, non-finite, or incorrectly aligned records
are reported and cause an error by default; no temporal filling is performed.
"""
function load_era5_panel(
    annual_paths::AbstractDict{<:Integer,<:AbstractString}, station_ids::Vector{String},
    feature_times::Vector{DateTime}; variables=collect(ERA5_VARIABLES),
    offset_hours::Int=9, strict::Bool=true, validate_annual_complete::Bool=true,
)
    length(unique(station_ids)) == length(station_ids) ||
        throw(ArgumentError("station IDs must be unique"))
    length(unique(feature_times)) == length(feature_times) ||
        throw(ArgumentError("feature times must be unique"))
    vars = Symbol.(variables)
    all(v -> v in ERA5_VARIABLES, vars) || throw(ArgumentError("unsupported ERA5 variable"))
    id_index = Dict(id => i for (i, id) in enumerate(station_ids))
    time_index = Dict(time => j for (j, time) in enumerate(feature_times))
    matrices = Dict(v => fill(NaN, length(station_ids), length(feature_times)) for v in vars)
    seen = falses(length(station_ids), length(feature_times))
    duplicates = 0
    offset_mismatch = 0
    matched = 0
    annual_expected = 0
    annual_matched = 0
    annual_duplicates = 0
    annual_nonfinite = 0
    annual_invalid_times = 0
    annual_unknown_stations = 0
    needed = vcat([:station_id, :time_utc, :time_bjt], vars)
    target_years = Set(year(t - Hour(offset_hours)) for t in feature_times)
    missing_years = sort(collect(setdiff(target_years, Set(Int.(keys(annual_paths))))))
    isempty(missing_years) || throw(ArgumentError("missing ERA5 annual paths: $(join(missing_years, ", "))"))

    for yr in sort(collect(target_years))
        path = String(annual_paths[yr])
        isfile(path) || throw(ArgumentError("ERA5 annual file does not exist: $path"))
        year_start = DateTime(yr, 1, 1)
        year_end = DateTime(yr + 1, 1, 1)
        hours_in_year = Int(Dates.value(year_end - year_start) ÷ 3_600_000)
        year_seen = falses(length(station_ids), hours_in_year)
        annual_expected += length(station_ids) * hours_in_year
        open(path, "r") do io
            eof(io) && throw(ArgumentError("empty ERA5 annual file: $path"))
            header = Symbol.(_csv_field.(split(chomp(readline(io)), ','; keepempty=true)))
            positions = Dict(name => findfirst(==(name), header) for name in needed)
            missing_columns = [name for (name, position) in positions if position === nothing]
            isempty(missing_columns) || throw(ArgumentError(
                "ERA5 file is missing columns: $(join(missing_columns, ", "))",
            ))
            for line in eachline(io)
                fields = split(chomp(line), ','; keepempty=true)
                sid = _csv_field(fields[positions[:station_id]])
                if !haskey(id_index, sid)
                    annual_unknown_stations += 1
                    continue
                end
                utc = _parse_datetime(_csv_field(fields[positions[:time_utc]]))
                slot = Int(Dates.value(utc - year_start) ÷ 3_600_000) + 1
                exact_hour = year_start <= utc < year_end &&
                    utc == year_start + Hour(slot - 1) && 1 <= slot <= hours_in_year
                if !exact_hour
                    annual_invalid_times += 1
                    continue
                end
                i = id_index[sid]
                if year_seen[i, slot]
                    annual_duplicates += 1
                else
                    year_seen[i, slot] = true
                    annual_matched += 1
                end
                bjt = _parse_datetime(_csv_field(fields[positions[:time_bjt]]))
                bjt == utc + Hour(8) || (offset_mismatch += 1)
                feature_time = align_feature_time(utc; offset_hours=offset_hours)
                is_target = haskey(time_index, feature_time)
                j = is_target ? time_index[feature_time] : 0
                target_duplicate = is_target && seen[i, j]
                if target_duplicate
                    duplicates += 1
                elseif is_target
                    seen[i, j] = true
                    matched += 1
                end
                for variable in vars
                    value = _finite_float(_csv_field(fields[positions[variable]]))
                    !isfinite(value) && (annual_nonfinite += 1)
                    is_target && !target_duplicate && (matrices[variable][i, j] = value)
                end
            end
        end
    end
    missing_cells = count(!, seen)
    nonfinite_cells = sum(count(x -> !isfinite(x), matrix) for matrix in values(matrices))
    qc = DataFrame(
        station_count=[length(station_ids)], target_time_count=[length(feature_times)],
        expected_station_hours=[length(station_ids) * length(feature_times)],
        matched_station_hours=[matched], duplicate_keys=[duplicates],
        missing_station_hours=[missing_cells], nonfinite_values=[nonfinite_cells],
        annual_expected_station_hours=[annual_expected], annual_matched_station_hours=[annual_matched],
        annual_duplicate_keys=[annual_duplicates],
        annual_missing_station_hours=[annual_expected - annual_matched],
        annual_nonfinite_values=[annual_nonfinite], annual_invalid_times=[annual_invalid_times],
        annual_unknown_station_rows=[annual_unknown_stations],
        time_bjt_offset_mismatches=[offset_mismatch], feature_time_offset_hours=[offset_hours],
        annual_complete=[annual_duplicates == 0 && annual_expected == annual_matched &&
            annual_nonfinite == 0 && annual_invalid_times == 0 && annual_unknown_stations == 0],
        complete=[duplicates == 0 && missing_cells == 0 && nonfinite_cells == 0 &&
            (!validate_annual_complete || (annual_duplicates == 0 && annual_expected == annual_matched &&
                annual_nonfinite == 0 && annual_invalid_times == 0 && annual_unknown_stations == 0)) &&
            offset_mismatch == 0],
    )
    if strict && !qc.complete[1]
        error("ERA5 station-hour quality control failed: $(NamedTuple(qc[1, :]))")
    end
    return (; values=matrices, qc)
end

include("era5/era5_panel_screen.jl")
include("era5/era5_panel_variability.jl")
include("era5/era5_selection_run.jl")

end
