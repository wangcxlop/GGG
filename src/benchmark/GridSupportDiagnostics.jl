"""
What the benchmark's satellite inputs lost when the gridded products were reduced to gauge points.

The benchmark never opens a satellite grid: every product reaches it as a wide `time x station_id`
table, sampled nearest-pixel from the FY4B disk (`FY4BPreprocessing.extract_precipitation`) or as a
0.1 degree cell mean from Earth Engine. This module measures the cost of that reduction, and it
measures only - it fits nothing, changes no benchmark output, and is not part of any run.

Six questions, one group of functions each:

- `pixel_assignment_table` / `effective_cell_keys`: which grid cell each gauge actually reads, on
  the FY4B disk and on the 0.1 degree lat/lon grid. FY4B needs the second function because the
  archive holds segments from two sub-satellite longitudes, so one gauge has two pixels.
- `cell_membership_table` / `cell_multiplicity_summary`: how many gauges share a cell, and so how
  many receive a satellite value that carries no information about where inside the cell they are.
- `identical_series_groups` / `identical_series_table` / `cell_grouping_agreement`: the same
  question answered from the shipped tables alone, with no grid geometry assumed, and then used to
  test the geometry. Both assumptions needed it: the single-projection FY4B grouping and the
  0.1 degree grid the Earth Engine exports were assumed to sample are each contradicted by what
  the stored values can actually distinguish. Where the two disagree, the data is the answer.
- `gauge_pair_table` / `pair_disagreement_summary`: how much two gauges inside one cell disagree
  with each other. That is the floor on point-to-pixel error - a satellite cannot match both - and
  the ceiling on what any downscaling could recover. Pairs that do *not* share a cell, at the same
  separation, are the control.
- `pixel_footprint_km` / `footprint_table`: FY4B's 4 km is at nadir, and the study area is 23
  degrees of longitude and 32 degrees of latitude away from the sub-satellite point. This is the
  real ground footprint there.
- `domain_raster_values` / `terrain_representativeness_table`: how biased the gauge network is as a
  sample of the terrain it is used to validate over.

Conventions shared with the rest of the benchmark: matrices are `[station, time]` with `NaN` for
missing; a wet hour is `>= WET_MM`; hour labels are hour-ending Beijing time.
"""
module GridSupportDiagnostics

using DataFrames, Dates, Statistics
using NCDatasets
using MixedGWR: metric_continuous
using Main: FY4BPreprocessing
using Main.TraditionalInterpolation: haversine_distance_matrix

export WET_MM, COARSE_CELL_DEG, FY4B_NADIR_RESOLUTION_M, SEPARATION_BIN_EDGES
export fy4b_pixel_index, latlon_cell_index, latlon_cell_center, pixel_footprint_km
export pixel_assignment_table, effective_cell_keys
export cell_membership_table, cell_multiplicity_summary
export identical_series_groups, identical_series_table
export series_group_keys, cell_grouping_agreement
export pair_disagreement, gauge_pair_table, pair_disagreement_summary
export footprint_table, domain_raster_values, terrain_representativeness_table
export bilinear_neighbours, bilinear_combine, sample_precipitation_modes
export fy4b_hourly_modes, extraction_sensitivity_table

"""A wet hour, matching the benchmark's `rain_threshold` and `wet_threshold`."""
const WET_MM = 0.1

"""
The lat/lon cell size GPM and GSMaP were assumed to be sampled on: `scale: 11132` m in the Earth
Engine exports, which is 0.1 degree at the equator. An assumption - see `cell_grouping_agreement`.
"""
const COARSE_CELL_DEG = 0.1

"""FY4B's nominal resolution, which `FY4BPreprocessing.RESOLUTION` states is *at nadir*."""
const FY4B_NADIR_RESOLUTION_M = 4000.0

"""
Earth radius for converting a grid step to ground km.

`FY4BPreprocessing`'s value rather than `TraditionalInterpolation`'s 6378.388: the footprint is a
property of the geostationary projection, so it has to use the same figure that projection does.
The two differ by 0.25 km in 6378, far below anything reported here.
"""
const EARTH_RADIUS_KM = 6378.137

"""Gauge separations, in km, that `pair_disagreement_summary` bins pairs into."""
const SEPARATION_BIN_EDGES = [0.0, 5.0, 10.0, 20.0, 40.0, 60.0]


# ---------------------------------------------------------------------------------------------
# Which cell a gauge reads
# ---------------------------------------------------------------------------------------------

"""
The 1-based FY4B disk pixel `FY4BPreprocessing.extract_precipitation` reads for a point, as
`(ix, iy, x, y)` where `(x, y)` is the fractional 0-based position it rounded.

`(0, 0)` marks the `NaN` that `latlon_to_xy` returns, which is rarer than it looks: the limb test
in `latlon_to_scan_angles` measures the angle *at the satellite*, and no point on Earth subtends
more than about 8.7 degrees from geostationary orbit, so it never fires and the far hemisphere
folds onto valid pixel indices instead of being rejected. What actually guards the disk is the
bounds check in `extract_precipitation`, repeated here wherever this module samples. Every gauge
in this study sits well inside the disk, so nothing below depends on the difference.
"""
function fy4b_pixel_index(lon::Real, lat::Real; sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON)
    x, y = FY4BPreprocessing.latlon_to_xy(lat, lon, sat_lon)
    (isnan(x) || isnan(y)) && return (0, 0, x, y)
    return (round(Int, x) + 1, round(Int, y) + 1, x, y)
end

"""
The index of the `step_deg` lat/lon cell containing a point, as `(ilon, ilat)`.

The grid is anchored at (-180, -90), which is how IMERG and GSMaP are gridded and therefore what
Earth Engine samples at `scale: 11132`. That anchoring is an assumption about the export, not
something this repository can read off the data - `identical_series_table` is what tests it.
"""
latlon_cell_index(lon::Real, lat::Real; step_deg::Real=COARSE_CELL_DEG) =
    (floor(Int, (Float64(lon) + 180) / step_deg), floor(Int, (Float64(lat) + 90) / step_deg))

"""Centre of the `step_deg` cell `(ilon, ilat)`, as `(lon, lat)`."""
latlon_cell_center(ilon::Integer, ilat::Integer; step_deg::Real=COARSE_CELL_DEG) =
    ((ilon + 0.5) * step_deg - 180, (ilat + 0.5) * step_deg - 90)


# ---------------------------------------------------------------------------------------------
# FY4B ground footprint
# ---------------------------------------------------------------------------------------------

"""Local east/north offset in km of `(lon, lat)` from `(lon0, lat0)`, for offsets of a few km."""
_local_km(lon, lat, lon0, lat0) = (
    (lon - lon0) * cosd(lat0) * (pi / 180) * EARTH_RADIUS_KM,
    (lat - lat0) * (pi / 180) * EARTH_RADIUS_KM,
)

"""
Ground size of one FY4B pixel at `(lon, lat)`, as
`(; along_scan_km, cross_scan_km, area_km2, nadir_ratio)`.

`RESOLUTION = 4000.0` in `FY4BPreprocessing` is the pixel size at the sub-satellite point. Away
from it the same scan-angle step subtends more ground, by a factor that grows with viewing angle.
This differentiates `latlon_to_xy` about the point to get the grid Jacobian, inverts it, and
measures the ground displacement of a one-pixel step along each grid axis. `nadir_ratio` is the
area relative to a 4 km square.

`h` is the finite-difference step in degrees: large enough that the scan angles differ in double
precision, small enough that the projection is linear across it.
"""
function pixel_footprint_km(
    lon::Real, lat::Real; sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON, h::Real=1e-3,
)
    blank = (; along_scan_km=NaN, cross_scan_km=NaN, area_km2=NaN, nadir_ratio=NaN)
    x_e, y_e = FY4BPreprocessing.latlon_to_xy(lat, lon + h, sat_lon)
    x_w, y_w = FY4BPreprocessing.latlon_to_xy(lat, lon - h, sat_lon)
    x_n, y_n = FY4BPreprocessing.latlon_to_xy(lat + h, lon, sat_lon)
    x_s, y_s = FY4BPreprocessing.latlon_to_xy(lat - h, lon, sat_lon)
    any(isnan, (x_e, y_e, x_w, y_w, x_n, y_n, x_s, y_s)) && return blank
    # J maps (dlon, dlat) to (dx, dy) in pixels.
    j11 = (x_e - x_w) / (2h); j12 = (x_n - x_s) / (2h)
    j21 = (y_e - y_w) / (2h); j22 = (y_n - y_s) / (2h)
    determinant = j11 * j22 - j12 * j21
    determinant == 0 && return blank
    # The inverse maps a one-pixel step back to (dlon, dlat).
    dlon_dx, dlat_dx = j22 / determinant, -j21 / determinant
    dlon_dy, dlat_dy = -j12 / determinant, j11 / determinant
    vx = _local_km(lon + dlon_dx, lat + dlat_dx, lon, lat)
    vy = _local_km(lon + dlon_dy, lat + dlat_dy, lon, lat)
    area = abs(vx[1] * vy[2] - vx[2] * vy[1])
    nadir = FY4B_NADIR_RESOLUTION_M / 1000
    return (; along_scan_km=hypot(vx...), cross_scan_km=hypot(vy...),
        area_km2=area, nadir_ratio=area / nadir^2)
end


# ---------------------------------------------------------------------------------------------
# Tables: assignment and multiplicity
# ---------------------------------------------------------------------------------------------

"""
One row per gauge: the FY4B pixel and the `step_deg` lat/lon cell it reads, how far it sits from
the FY4B pixel centre, and how large that pixel is on the ground.

`fy4b_offset_km` carries the fractional pixel offset through the same Jacobian
`pixel_footprint_km` uses, so it is the ground distance from the gauge to the centre of the value
it is handed.
"""
function pixel_assignment_table(
    station_ids::AbstractVector{<:AbstractString},
    lon::AbstractVector{<:Real}, lat::AbstractVector{<:Real};
    sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON, step_deg::Real=COARSE_CELL_DEG,
)
    length(station_ids) == length(lon) == length(lat) ||
        throw(DimensionMismatch("station_ids, lon and lat must have the same length"))
    rows = NamedTuple[]
    for i in eachindex(station_ids)
        ix, iy, x, y = fy4b_pixel_index(lon[i], lat[i]; sat_lon)
        footprint = pixel_footprint_km(lon[i], lat[i]; sat_lon)
        dx = isnan(x) ? NaN : x - round(x)
        dy = isnan(y) ? NaN : y - round(y)
        ilon, ilat = latlon_cell_index(lon[i], lat[i]; step_deg)
        clon, clat = latlon_cell_center(ilon, ilat; step_deg)
        push!(rows, (
            station_id=station_ids[i], lon=Float64(lon[i]), lat=Float64(lat[i]),
            sat_lon=Float64(sat_lon),
            fy4b_x=x, fy4b_y=y, fy4b_ix=ix, fy4b_iy=iy, fy4b_pixel="$(ix)_$(iy)",
            fy4b_offset_pixels=isnan(dx) ? NaN : hypot(dx, dy),
            fy4b_offset_km=isnan(dx) ? NaN :
                hypot(dx * footprint.along_scan_km, dy * footprint.cross_scan_km),
            fy4b_along_scan_km=footprint.along_scan_km,
            fy4b_cross_scan_km=footprint.cross_scan_km,
            fy4b_footprint_km2=footprint.area_km2,
            coarse_cell="$(ilon)_$(ilat)", coarse_cell_lon=clon, coarse_cell_lat=clat,
        ))
    end
    return DataFrame(rows)
end

"""
The cell keys a product's stored values actually follow when it was sampled through more than one
grid, as the intersection of the per-grid keys.

The FY4B archive carries segments from two sub-satellite longitudes - 105.0E and 133.0E, both
present for this record - and `extract_precipitation` reads each file's own subpoint. A gauge's
pixel therefore differs between the two projections, and two gauges are handed the same value
only where they share a pixel in *both*. The single-projection grouping overstates how much the
shipped column merges; this is the partition that governs it.
"""
function effective_cell_keys(keys_by_grid::AbstractVector{<:AbstractVector})
    isempty(keys_by_grid) && throw(ArgumentError("at least one grid is required"))
    allequal(length.(keys_by_grid)) ||
        throw(DimensionMismatch("every grid must supply a key for every station"))
    return [join((string(keys[i]) for keys in keys_by_grid), "|")
            for i in eachindex(first(keys_by_grid))]
end

"""
One row per (grid, occupied cell): which gauges fall in it.

`cells_by_grid` maps a grid name to that grid's cell key for each gauge, in `station_ids` order.
"""
function cell_membership_table(
    station_ids::AbstractVector{<:AbstractString},
    cells_by_grid::AbstractDict{<:AbstractString,<:AbstractVector},
)
    rows = NamedTuple[]
    for grid in sort(collect(keys(cells_by_grid)))
        cells = cells_by_grid[grid]
        length(cells) == length(station_ids) || throw(DimensionMismatch(
            "grid $grid has $(length(cells)) cells for $(length(station_ids)) stations"))
        members = Dict{String,Vector{String}}()
        for i in eachindex(station_ids)
            push!(get!(members, string(cells[i]), String[]), station_ids[i])
        end
        for cell in sort(collect(keys(members)))
            ids = sort(members[cell])
            push!(rows, (grid=grid, cell=cell, n_stations=length(ids), station_ids=join(ids, ";")))
        end
    end
    return DataFrame(rows)
end

"""
One row per grid: how much of the gauge network it resolves.

`n_stations_sharing` counts gauges that are *not* alone in their cell - the gauges whose satellite
value cannot distinguish them from a neighbour.
"""
function cell_multiplicity_summary(membership::DataFrame)
    rows = NamedTuple[]
    for grid in sort(unique(membership.grid))
        cells = membership[membership.grid .== grid, :]
        counts = cells.n_stations
        sharing = sum(counts[counts .> 1])
        total = sum(counts)
        push!(rows, (
            grid=grid, n_stations=total, n_occupied_cells=nrow(cells),
            max_stations_per_cell=maximum(counts),
            n_cells_with_multiple=count(>(1), counts),
            n_stations_sharing=sharing,
            fraction_stations_sharing=sharing / total,
            mean_stations_per_occupied_cell=total / nrow(cells),
        ))
    end
    return DataFrame(rows)
end


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


# ---------------------------------------------------------------------------------------------
# Footprint, and how much of the domain a gauge speaks for
# ---------------------------------------------------------------------------------------------

"""
FY4B's real ground footprint over the study area, plus how much of it the gauge network samples.

`probe_points` are labelled `(name, lon, lat)` reference locations - the study box corners and
centre - reported beside the gauge mean so the footprint is anchored to the domain rather than to
where gauges happen to be. The `coverage` row counts distinct FY4B pixels covering the box and,
in `nadir_ratio`, the fraction of them holding a gauge: the share of the domain a gauge-only score
can speak for. `samples` points per axis must be fine enough that no pixel is stepped over.
"""
function footprint_table(
    assignment::DataFrame, probe_points::AbstractVector{<:Tuple{<:AbstractString,<:Real,<:Real}},
    bounds::NamedTuple; sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON, samples::Int=400,
)
    rows = NamedTuple[]
    for (name, lon, lat) in probe_points
        footprint = pixel_footprint_km(lon, lat; sat_lon)
        push!(rows, (scope="probe", name=name, lon=Float64(lon), lat=Float64(lat),
            along_scan_km=footprint.along_scan_km, cross_scan_km=footprint.cross_scan_km,
            footprint_km2=footprint.area_km2, nadir_ratio=footprint.nadir_ratio))
    end
    nadir_area = (FY4B_NADIR_RESOLUTION_M / 1000)^2
    push!(rows, (scope="gauges", name="mean", lon=NaN, lat=NaN,
        along_scan_km=mean(assignment.fy4b_along_scan_km),
        cross_scan_km=mean(assignment.fy4b_cross_scan_km),
        footprint_km2=mean(assignment.fy4b_footprint_km2),
        nadir_ratio=mean(assignment.fy4b_footprint_km2) / nadir_area))
    covering = Set{Tuple{Int,Int}}()
    for lon in range(bounds.west, bounds.east; length=samples),
        lat in range(bounds.south, bounds.north; length=samples)
        ix, iy, _, _ = fy4b_pixel_index(lon, lat; sat_lon)
        ix == 0 && continue
        push!(covering, (ix, iy))
    end
    occupied = length(unique(assignment.fy4b_pixel))
    push!(rows, (scope="coverage", name="fy4b_pixels_in_study_box", lon=NaN, lat=NaN,
        along_scan_km=NaN, cross_scan_km=NaN, footprint_km2=Float64(length(covering)),
        nadir_ratio=occupied / length(covering)))
    return DataFrame(rows)
end

"""
Finite cell values of `raster_path` resampled onto the study grid, as a flat vector.

Shells out to `gdalwarp` then `gdal_translate -of XYZ`, the route `TerrainFeatures` and
`LandformClassification` already take; nothing is written outside a temporary directory.
"""
function domain_raster_values(
    raster_path::AbstractString, bounds::NamedTuple, step_deg::Real;
    resampling::AbstractString="average", nodata::Real=-9999.0,
)
    isfile(raster_path) || error("Missing raster: $raster_path")
    gdalwarp = Main.TerrainFeatures.find_tool("gdalwarp")
    gdal_translate = Main.TerrainFeatures.find_tool("gdal_translate")
    return mktempdir() do workdir
        tif = joinpath(workdir, "domain.tif")
        xyz = joinpath(workdir, "domain.xyz")
        run(pipeline(Cmd(Cmd([
            gdalwarp, "-q", "-overwrite", "-t_srs", "EPSG:4326",
            "-te", string(bounds.west), string(bounds.south),
            string(bounds.east), string(bounds.north),
            "-tr", string(step_deg), string(step_deg), "-r", resampling,
            "-dstnodata", string(Int(nodata)), raster_path, tif,
        ]); dir=workdir); stderr=devnull))
        run(pipeline(Cmd(Cmd([gdal_translate, "-q", "-of", "XYZ", tif, xyz]);
            dir=workdir); stderr=devnull))
        grid = Main.LandformClassification.parse_xyz_grid(eachline(xyz); nodata)
        return filter(isfinite, vec(grid.grid))
    end
end

"""
How the gauge network's terrain compares with the domain it is used to validate over.

Deciles and moments of both distributions, plus the share of the domain lying outside the range
any gauge occupies - the part of the study area no gauge-based score can speak for.
"""
function terrain_representativeness_table(
    gauge_values::AbstractVector{<:Real}, domain_values::AbstractVector{<:Real},
    variable::AbstractString; probs::AbstractVector{<:Real}=collect(0.0:0.1:1.0),
)
    gauges = filter(isfinite, Float64.(gauge_values))
    domain = filter(isfinite, Float64.(domain_values))
    (isempty(gauges) || isempty(domain)) && error("$variable has no finite values to compare")
    rows = NamedTuple[]
    for p in probs
        push!(rows, (variable=variable, statistic="p$(round(Int, 100p))",
            gauges=quantile(gauges, p), domain=quantile(domain, p)))
    end
    push!(rows, (variable=variable, statistic="mean", gauges=mean(gauges), domain=mean(domain)))
    push!(rows, (variable=variable, statistic="sd", gauges=std(gauges), domain=std(domain)))
    push!(rows, (variable=variable, statistic="n",
        gauges=Float64(length(gauges)), domain=Float64(length(domain))))
    lo, hi = extrema(gauges)
    push!(rows, (variable=variable, statistic="domain_fraction_below_lowest_gauge",
        gauges=NaN, domain=count(<(lo), domain) / length(domain)))
    push!(rows, (variable=variable, statistic="domain_fraction_above_highest_gauge",
        gauges=NaN, domain=count(>(hi), domain) / length(domain)))
    push!(rows, (variable=variable, statistic="domain_fraction_outside_gauge_range",
        gauges=NaN, domain=count(value -> value < lo || value > hi, domain) / length(domain)))
    return DataFrame(rows)
end


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

end # module
