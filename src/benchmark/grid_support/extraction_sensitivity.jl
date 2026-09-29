# ---------------------------------------------------------------------------------------------
# What nearest-pixel extraction discards
# ---------------------------------------------------------------------------------------------

"""
The four 1-based grid neighbours of the fractional 0-based coordinate `(x, y)`, with their
bilinear weights, as `(indices, weights)` ordered (lo,lo), (hi,lo), (lo,hi), (hi,hi).
"""
function bilinear_neighbours(x::Real, y::Real)
    ix = floor(Int, x); iy = floor(Int, y)
    fx = x - ix; fy = y - iy
    indices = ((ix + 1, iy + 1), (ix + 2, iy + 1), (ix + 1, iy + 2), (ix + 2, iy + 2))
    weights = ((1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy)
    return indices, weights
end

"""
Bilinear estimate from four neighbour values, renormalised over the ones that passed QC.

`NaN` marks a rejected neighbour: it drops out and the remaining weights are rescaled, so a pixel
beside a rejected one still yields a value rather than a hole. All four rejected, or the accepted
ones carrying no weight, gives `NaN`.
"""
function bilinear_combine(values::NTuple{4,Float64}, weights::NTuple{4,Float64})
    total = 0.0
    accumulated = 0.0
    @inbounds for k in 1:4
        isnan(values[k]) && continue
        total += weights[k]
        accumulated += weights[k] * values[k]
    end
    return total > 0 ? accumulated / total : NaN
end

"""`true` when a raw FY4B value and its quality flag pass `extract_precipitation`'s QC."""
_fy4b_accepts(value, quality) =
    !ismissing(value) && !ismissing(quality) && value != Float32(65534.0) &&
    value >= 0 && value <= 30 && Int(quality) < 2

"""
One FY4B segment sampled at `lon`/`lat` both ways, as `(nearest, bilinear)`.

`nearest` reproduces `FY4BPreprocessing.extract_precipitation` exactly, including its QC and its
satellite-longitude handling, so the pair is a controlled comparison: the two differ only in the
interpolation, never in which file, hour or quality rule produced them.
"""
function sample_precipitation_modes(
    nc_file::AbstractString, lon::AbstractVector{<:Real}, lat::AbstractVector{<:Real},
    coord_cache::AbstractDict,
)
    ds = Dataset(nc_file)
    try
        precip = ds["Precipitation"]
        quality = haskey(ds, "DQF") ? ds["DQF"] : nothing
        sat_lon = FY4BPreprocessing.dataset_satellite_lon(ds, basename(nc_file))
        if !haskey(coord_cache, sat_lon)
            coord_cache[sat_lon] = [FY4BPreprocessing.latlon_to_xy(lat[i], lon[i], sat_lon)
                                    for i in eachindex(lon)]
        end
        coords = coord_cache[sat_lon]
        size_x = FY4BPreprocessing.GRID_SIZE
        inside(ix, iy) = 1 <= ix <= size_x && 1 <= iy <= size_x
        read_value(ix, iy) = begin
            value = precip[ix, iy]
            flag = quality === nothing ? 0 : quality[ix, iy]
            _fy4b_accepts(value, flag) ? Float64(value) : NaN
        end
        nearest = fill(NaN, length(lon))
        bilinear = fill(NaN, length(lon))
        for i in eachindex(lon)
            x, y = coords[i]
            (isnan(x) || isnan(y)) && continue
            ix = round(Int, x) + 1; iy = round(Int, y) + 1
            inside(ix, iy) && (nearest[i] = read_value(ix, iy))
            indices, weights = bilinear_neighbours(x, y)
            values = ntuple(k -> inside(indices[k]...) ? read_value(indices[k]...) : NaN, 4)
            bilinear[i] = bilinear_combine(values, weights)
        end
        return nearest, bilinear
    finally
        close(ds)
    end
end

"""
Hourly FY4B at `lon`/`lat` by both sampling modes, as `(times, Y_nearest, Y_bilinear)` with the
matrices `[station, time]`.

`complete_groups` is what `FY4BPreprocessing.build_hourly_qc` returns: strict-complete hours only,
four 15-minute segments each. The aggregation matches `FY4BPreprocessing.aggregate_hourly` - the
four rates average to the hourly accumulation, and a single missing segment voids the hour - and
`times` carries the same hour-ending Beijing labels, so the result lines up with the shipped table.
Both modes come out of one pass, so each file is read once.
"""
function fy4b_hourly_modes(
    complete_groups::AbstractDict, lon::AbstractVector{<:Real}, lat::AbstractVector{<:Real};
    progress_every::Int=200,
)
    hours = sort(collect(keys(complete_groups)))
    n = length(lon)
    Y_nearest = fill(NaN, n, length(hours))
    Y_bilinear = fill(NaN, n, length(hours))
    coord_cache = Dict{Float64,Vector{Tuple{Float64,Float64}}}()
    times = DateTime[]
    for (index, hour_start) in enumerate(hours)
        sums = (zeros(Float64, n), zeros(Float64, n))
        gaps = (falses(n), falses(n))
        for nc_file in complete_groups[hour_start]
            nearest, bilinear = sample_precipitation_modes(nc_file, lon, lat, coord_cache)
            for (values, total, gap) in ((nearest, sums[1], gaps[1]), (bilinear, sums[2], gaps[2]))
                for i in 1:n
                    isnan(values[i]) ? (gap[i] = true) : (total[i] += values[i])
                end
            end
        end
        for i in 1:n
            Y_nearest[i, index] = gaps[1][i] ? NaN : sums[1][i] / 4
            Y_bilinear[i, index] = gaps[2][i] ? NaN : sums[2][i] / 4
        end
        push!(times, FY4BPreprocessing.aligned_hour_end(hour_start))
        if progress_every > 0 && (index % progress_every == 0 || index == length(hours))
            println("  sampled $index / $(length(hours)) strict hours")
        end
    end
    return times, Y_nearest, Y_bilinear
end

"""
How much the nearest-pixel rule changed the value a gauge was given, against bilinear.

One row per gauge plus an `all_stations` row. `n_wet_flip` counts hours the two modes disagree
about whether it rained at all, which is what a detection score would see. When `Y_reference` is
given - the shipped column on the same hours and gauges - `reference_mismatch` counts cells where
this module's nearest sampling failed to reproduce it, and must be zero for the rest to mean
anything.
"""
function extraction_sensitivity_table(
    station_ids::AbstractVector{<:AbstractString},
    Y_nearest::AbstractMatrix{<:Real}, Y_bilinear::AbstractMatrix{<:Real};
    Y_reference::Union{Nothing,AbstractMatrix{<:Real}}=nothing, wet_mm::Real=WET_MM,
)
    size(Y_nearest) == size(Y_bilinear) ||
        throw(DimensionMismatch("the two sampling modes must have the same shape"))
    size(Y_nearest, 1) == length(station_ids) ||
        throw(DimensionMismatch("matrices must have one row per station"))
    Y_reference === nothing || size(Y_reference) == size(Y_nearest) ||
        throw(DimensionMismatch("Y_reference must match the sampled shape"))
    function summarise(scope, id, a, b, reference)
        both = .!isnan.(a) .& .!isnan.(b)
        n = count(both)
        delta = n > 0 ? (b[both] .- a[both]) : Float64[]
        wet_a = both .& (a .>= wet_mm)
        wet_b = both .& (b .>= wet_mm)
        either_wet = wet_a .| wet_b
        n_either_wet = count(either_wet)
        wet_delta = b[either_wet] .- a[either_wet]
        flips = count(wet_a .!= wet_b)
        return (
            scope=scope, station_id=id, n=n,
            coverage_nearest=count(.!isnan.(a)) / length(a),
            coverage_bilinear=count(.!isnan.(b)) / length(b),
            mean_nearest=n > 0 ? mean(a[both]) : NaN,
            mean_bilinear=n > 0 ? mean(b[both]) : NaN,
            mean_abs_delta=n > 0 ? mean(abs.(delta)) : NaN,
            rmse_delta=n > 0 ? sqrt(mean(delta .^ 2)) : NaN,
            max_abs_delta=n > 0 ? maximum(abs.(delta)) : NaN,
            r=n > 1 ? cor(a[both], b[both]) : NaN,
            n_wet_nearest=count(wet_a), n_wet_bilinear=count(wet_b),
            n_wet_either=n_either_wet,
            # Over wet hours, not over the record: 93% of hours are dry at both, and averaging the
            # change over those hides it. These are the numbers a wet-hour or detection score sees.
            rmse_delta_wet=n_either_wet > 0 ? sqrt(mean(wet_delta .^ 2)) : NaN,
            mean_abs_delta_wet=n_either_wet > 0 ? mean(abs.(wet_delta)) : NaN,
            n_wet_flip=flips, fraction_wet_flip=n > 0 ? flips / n : NaN,
            fraction_wet_flip_of_wet=n_either_wet > 0 ? flips / n_either_wet : NaN,
            reference_mismatch=reference === nothing ? -1 :
                count(index -> !isequal(a[index], reference[index]), eachindex(a)),
        )
    end
    reference_row(rows...) = Y_reference === nothing ? nothing : vec(Y_reference[rows...])
    rows = NamedTuple[summarise("all_stations", "", vec(Y_nearest), vec(Y_bilinear),
        Y_reference === nothing ? nothing : vec(Y_reference))]
    for i in eachindex(station_ids)
        push!(rows, summarise("station", station_ids[i], vec(Y_nearest[i, :]),
            vec(Y_bilinear[i, :]), reference_row(i, :)))
    end
    return DataFrame(rows)
end
