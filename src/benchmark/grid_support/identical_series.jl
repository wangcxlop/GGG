# ---------------------------------------------------------------------------------------------
# The same question without assuming a grid
# ---------------------------------------------------------------------------------------------

"""
`(indistinguishable, overlap)` for rows `i` and `j` of `Y`: whether they ever differ on an hour
where both are finite, and how many such hours there are.

Equality is exact, because two gauges reading one cell are handed the identical stored value, so
any difference at all means different cells. A pair with fewer than `min_overlap` shared hours is
not called indistinguishable.
"""
function _series_indistinguishable(Y::AbstractMatrix{<:Real}, i::Int, j::Int, min_overlap::Int)
    overlap = 0
    @inbounds for t in axes(Y, 2)
        a = Y[i, t]; b = Y[j, t]
        (isnan(a) || isnan(b)) && continue
        a == b || return (false, overlap)
        overlap += 1
    end
    return (overlap >= min_overlap, overlap)
end

"""
`(groups, overlaps)`: gauges whose series are indistinguishable, as vectors of row indices into
`Y`, and the thinnest overlap each member was joined on.

Only groups of two or more are returned. Indistinguishability is not transitive when coverage
differs between gauges, so groups are the connected components of the relation rather than its
equivalence classes; the overlap counts are what make a component held together by a thin link
visible.
"""
function identical_series_groups(Y::AbstractMatrix{<:Real}; min_overlap::Int=1)
    n = size(Y, 1)
    parent = collect(1:n)
    root(a) = parent[a] == a ? a : (parent[a] = root(parent[a]))
    overlaps = Dict{Int,Int}()
    for i in 1:(n - 1), j in (i + 1):n
        same, overlap = _series_indistinguishable(Y, i, j, min_overlap)
        same || continue
        ri, rj = root(i), root(j)
        ri == rj || (parent[ri] = rj)
        overlaps[i] = haskey(overlaps, i) ? min(overlaps[i], overlap) : overlap
        overlaps[j] = haskey(overlaps, j) ? min(overlaps[j], overlap) : overlap
    end
    components = Dict{Int,Vector{Int}}()
    for i in 1:n
        push!(get!(components, root(i), Int[]), i)
    end
    groups = [sort(members) for members in values(components) if length(members) > 1]
    return sort(groups; by=first), overlaps
end

"""
One row per group of gauges that `product`'s shipped table cannot tell apart.

`min_overlap_hours` is the thinnest pairwise overlap holding the group together.
"""
function identical_series_table(
    station_ids::AbstractVector{<:AbstractString}, Y::AbstractMatrix{<:Real},
    product::AbstractString; min_overlap::Int=1,
)
    size(Y, 1) == length(station_ids) ||
        throw(DimensionMismatch("Y must have one row per station"))
    groups, overlaps = identical_series_groups(Y; min_overlap)
    rows = NamedTuple[]
    for (index, members) in enumerate(groups)
        push!(rows, (
            product=product, group=index, n_stations=length(members),
            station_ids=join(station_ids[members], ";"),
            min_overlap_hours=minimum(get(overlaps, member, 0) for member in members),
        ))
    end
    return DataFrame(rows)
end


"""
Per-station cell keys taken from `identical_series_groups` rather than from any grid.

The assumption-free counterpart to a grid's cell key: gauges the product cannot separate get one
key, everything else gets its own. Use it wherever "shares a cell" has to be right even though
which grid the product was sampled on is unknown - for GPM and GSMaP, that is everywhere.
"""
function series_group_keys(
    station_ids::AbstractVector{<:AbstractString},
    series_groups::AbstractVector{<:AbstractVector{Int}},
)
    keys = ["solo_$(id)" for id in station_ids]
    for (index, members) in enumerate(series_groups), member in members
        keys[member] = "group_$index"
    end
    return keys
end

"""
Whether the assumed grid explains the gauges the shipped table cannot tell apart, as one row.

Sharing a cell forces identical values, so every group of cell-sharing gauges has to sit inside a
single identical-series group: `cell_groups_confirmed` below `n_cell_groups` means the grid this
module assumed is not the grid the product was sampled on. The converse does not follow - two
gauges in different cells that both read zero on every shared hour are indistinguishable too - so
`stations_identical_across_cells` is reported as context rather than treated as a failure.

`series_groups` lets a caller supply the grouping once and test several candidate grids against
it: the relation costs a pass over every gauge pair and every hour, while swapping the grid costs
nothing.
"""
function cell_grouping_agreement(
    station_ids::AbstractVector{<:AbstractString}, cells::AbstractVector,
    Y::AbstractMatrix{<:Real}; grid::AbstractString, product::AbstractString, min_overlap::Int=1,
    series_groups::Union{Nothing,AbstractVector{<:AbstractVector{Int}}}=nothing,
)
    length(cells) == length(station_ids) == size(Y, 1) ||
        throw(DimensionMismatch("cells, station_ids and Y must describe the same stations"))
    by_cell = Dict{String,Vector{Int}}()
    for i in eachindex(station_ids)
        push!(get!(by_cell, string(cells[i]), Int[]), i)
    end
    cell_groups = [members for members in values(by_cell) if length(members) > 1]
    series_groups = series_groups === nothing ? first(identical_series_groups(Y; min_overlap)) :
        series_groups
    series_of = Dict{Int,Int}()
    for (index, members) in enumerate(series_groups), member in members
        series_of[member] = index
    end
    confirmed = count(cell_groups) do members
        first_group = get(series_of, first(members), 0)
        first_group != 0 && all(get(series_of, member, 0) == first_group for member in members)
    end
    across = 0
    for members in series_groups
        reference = string(cells[first(members)])
        across += count(member -> string(cells[member]) != reference, members)
    end
    return DataFrame([(
        product=product, grid=grid,
        n_cell_groups=length(cell_groups),
        n_stations_in_cell_groups=sum(length, cell_groups; init=0),
        n_series_groups=length(series_groups),
        n_stations_in_series_groups=sum(length, series_groups; init=0),
        cell_groups_confirmed=confirmed,
        stations_identical_across_cells=across,
        grid_explains_series=confirmed == length(cell_groups),
    )])
end
