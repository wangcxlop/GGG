# ---------------------------------------------------------------------------------------------
# What two gauges inside one cell disagree by
# ---------------------------------------------------------------------------------------------

"""
Disagreement between two gauge series, over all jointly observed hours and over their wet hours.

A wet hour is one where *either* gauge reaches `wet_mm`: restricting to hours both call wet would
discard exactly the disagreements that matter.
"""
function pair_disagreement(
    a::AbstractVector{<:Real}, b::AbstractVector{<:Real}; wet_mm::Real=WET_MM,
)
    both = .!isnan.(a) .& .!isnan.(b)
    count(both) == 0 && return (; n=0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN,
        n_wet=0, RMSE_wet=NaN, MAE_wet=NaN, mean_mm=NaN)
    all_hours = metric_continuous(a, b; mask=both)
    wet = both .& ((a .>= wet_mm) .| (b .>= wet_mm))
    n_wet = count(wet)
    wet_hours = n_wet > 0 ? metric_continuous(a, b; mask=wet) : (; RMSE=NaN, MAE=NaN)
    return (; n=all_hours.n, RMSE=all_hours.RMSE, MAE=all_hours.MAE, Bias=all_hours.Bias,
        r=all_hours.r, n_wet, RMSE_wet=wet_hours.RMSE, MAE_wet=wet_hours.MAE,
        mean_mm=mean(a[both]))
end

"""
One row per gauge pair closer than `max_separation_km`: how far apart they are, whether each grid
gives them the same cell, and how much their own observations disagree.

Pairs that share a cell are the point-to-pixel floor; pairs at the same separation that do not are
the control, so the table carries both and `pair_disagreement_summary` bins them together.
"""
function gauge_pair_table(
    station_ids::AbstractVector{<:AbstractString}, lonlat::AbstractMatrix{<:Real},
    Y_obs::AbstractMatrix{<:Real},
    cells_by_grid::AbstractDict{<:AbstractString,<:AbstractVector};
    max_separation_km::Real=60.0, wet_mm::Real=WET_MM, min_hours::Int=1,
)
    n = length(station_ids)
    size(Y_obs, 1) == n || throw(DimensionMismatch("Y_obs must have one row per station"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("lonlat must be n x 2"))
    distances = haversine_distance_matrix(lonlat, lonlat)
    grids = sort(collect(keys(cells_by_grid)))
    rows = NamedTuple[]
    for i in 1:(n - 1), j in (i + 1):n
        separation = distances[i, j]
        separation <= max_separation_km || continue
        stats = pair_disagreement(view(Y_obs, i, :), view(Y_obs, j, :); wet_mm)
        stats.n >= min_hours || continue
        shares = (; (Symbol("shares_", grid) =>
            string(cells_by_grid[grid][i]) == string(cells_by_grid[grid][j]) for grid in grids)...)
        push!(rows, merge(
            (station_a=station_ids[i], station_b=station_ids[j], separation_km=separation),
            shares,
            (n=stats.n, RMSE=stats.RMSE, MAE=stats.MAE, Bias=stats.Bias, r=stats.r,
                n_wet=stats.n_wet, RMSE_wet=stats.RMSE_wet, MAE_wet=stats.MAE_wet,
                mean_mm=stats.mean_mm),
        ))
    end
    return DataFrame(rows)
end

"""Label of the `edges` bin holding `value`, as `lo_hi` km; the last bin is closed on the right."""
function _separation_bin(value::Real, edges::AbstractVector{<:Real})
    index = clamp(searchsortedlast(edges, value), 1, length(edges) - 1)
    return "$(edges[index])_$(edges[index + 1])"
end

"""
`pairs` summarised by (grid, whether the pair shares a cell, separation bin).

The comparison to read is one row against the row beside it: same distance bin, differing only in
whether one satellite value has to serve both gauges.
"""
function pair_disagreement_summary(
    pairs::DataFrame; edges::AbstractVector{<:Real}=SEPARATION_BIN_EDGES,
)
    grids = [String(name)[(length("shares_") + 1):end]
             for name in names(pairs) if startswith(String(name), "shares_")]
    bins = _separation_bin.(pairs.separation_km, Ref(edges))
    rows = NamedTuple[]
    for grid in grids
        column = pairs[!, Symbol("shares_", grid)]
        for shares in (true, false), bin in sort(unique(bins))
            selected = (column .== shares) .& (bins .== bin)
            any(selected) || continue
            inside = pairs[selected, :]
            wet = filter(isfinite, inside.RMSE_wet)
            finite_r = filter(isfinite, inside.r)
            push!(rows, (
                grid=grid, shares_cell=shares, separation_bin_km=bin, n_pairs=nrow(inside),
                median_separation_km=median(inside.separation_km),
                median_RMSE=median(inside.RMSE), mean_RMSE=mean(inside.RMSE),
                median_MAE=median(inside.MAE),
                median_r=isempty(finite_r) ? NaN : median(finite_r),
                median_RMSE_wet=isempty(wet) ? NaN : median(wet),
                median_n_hours=median(inside.n),
            ))
        end
    end
    return sort(DataFrame(rows), [:grid, order(:shares_cell; rev=true), :separation_bin_km])
end
