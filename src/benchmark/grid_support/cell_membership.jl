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
