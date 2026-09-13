"""
Landform regions for the gauge network, from DEM local relief.

The classifier is relief amplitude - the elevation range (max - min) inside a square window - with
the breaks of the Chinese basic landform scheme: plain < 30 m, hills 30-200 m, and mountains split
at 500 and 1000 m into small-, medium- and large-relief mountains (the same breaks as "A
Geomorphological Regionalization using the Upscaled DEM: the Beijing-Tianjin-Hebei Area, China Case
Study", Sci. Rep. 2020, doi:10.1038/s41598-020-66993-9, which also splits 30-75 m off as tableland).
Altitude (low < 1000 m) is recorded beside it, and window-mean slope is reported as a coherence check;
neither reassigns a gauge.

Relief grows with the window, so the window is not picked by hand. `optimal_relief_window` applies
the mean change-point method to the area-wide mean relief over a fixed ladder of window sizes, and
that choice is made before any gauge is classified.

Why this network needs the sub-classes: sampled from the 30 m Copernicus DEM, the gauges hold at most
one "plain" (a flat reservoir surface) and 6-63 "hills" depending on the window. A plain/hills/mountains
split would leave one populated class, so the evaluation regions are the relief groups
`RELIEF_REGIONS`, with hills pooled into the < 500 m group.

GMT writes a `gmt.history` file into its working directory even with `--GMT_HISTORY=false`, so every
GMT call here runs inside a temporary directory rather than the caller's.
"""
module LandformClassification

using CSV, DataFrames, Statistics
using Main.TableIO: write_csv_atomic
using Main.TerrainFeatures: find_tool

export RELIEF_CLASSES, RELIEF_REGIONS, relief_class, relief_region, altitude_class
export mean_change_point, optimal_relief_window, block_relief, mean_block_relief, parse_xyz_grid
export prepare_station_landform

const RELIEF_BREAKS_M = [30.0, 200.0, 500.0, 1000.0]
const RELIEF_CLASSES = ["plain", "hills", "small_relief_mountain", "medium_relief_mountain", "large_relief_mountain"]
const RELIEF_REGIONS = ["relief_lt500", "relief_500_1000", "relief_ge1000"]
const NODATA = -9999.0
const UTM_PROJ = "+proj=utm +zone=49 +datum=WGS84 +units=m"

function _check_relief(relief_m::Real)
    isfinite(relief_m) && relief_m >= 0 ||
        throw(ArgumentError("relief must be finite and non-negative, got $relief_m"))
    return Float64(relief_m)
end

"""Basic landform class of a relief amplitude in metres; each break belongs to the class above it."""
relief_class(relief_m::Real) = RELIEF_CLASSES[searchsortedlast(RELIEF_BREAKS_M, _check_relief(relief_m)) + 1]

"""Evaluation region of a relief amplitude: hills and small-relief mountains pool below 500 m."""
function relief_region(relief_m::Real)
    relief = _check_relief(relief_m)
    return relief < 500 ? RELIEF_REGIONS[1] : relief < 1000 ? RELIEF_REGIONS[2] : RELIEF_REGIONS[3]
end

"""Altitude class of the same scheme: low < 1000 m, middle < 3500 m, high otherwise."""
function altitude_class(elevation_m::Real)
    isfinite(elevation_m) || throw(ArgumentError("elevation must be finite"))
    return elevation_m < 1000 ? "low_altitude" : elevation_m < 3500 ? "middle_altitude" : "high_altitude"
end

"""
Mean change-point of a sequence.

Splitting `x` before index `k` leaves the segments `x[1:k-1]` and `x[k:end]`. With `S` the sum of
squared deviations of `x` from its mean and `S_k` the sum of the two segments' squared deviations from
their own means, the change point is the `k` that maximises `S - S_k`; ties go to the earlier `k`.
Returns `(; index, gain)` with `gain[k] = S - S_k` (`NaN` at `k = 1`, which is not a split).
"""
function mean_change_point(x::AbstractVector{<:Real})
    n = length(x)
    n >= 3 || throw(ArgumentError("need at least three values"))
    all(isfinite, x) || throw(ArgumentError("values must be finite"))
    values = Float64.(x)
    sum_squares(v) = sum(abs2, v .- mean(v))
    total = sum_squares(values)
    gain = fill(NaN, n)
    for k in 2:n
        gain[k] = total - (sum_squares(@view values[1:k-1]) + sum_squares(@view values[k:n]))
    end
    return (; index=argmax(k -> gain[k], 2:n), gain)
end

"""
The relief window chosen by the mean change-point method.

`mean_relief_m[i]` is the area-wide mean relief for square windows of side `side_km[i]`. The analysed
sequence is the log of relief per unit area, `X_i = ln(RA_i / A_i)` with `A_i = side_km[i]^2`: it
falls steeply while the window is still growing into the terrain and flattens once it spans whole
hillslopes, and the change point marks that transition. The unit of `A` only shifts `X` and does not
move the change point.

Before quoting the result, check this sequence definition against the method's source (Liu Xinhua et
al. 2001, after Tu & Liu 1991); the implementation follows the usual description of it.
"""
function optimal_relief_window(side_km::AbstractVector{<:Real}, mean_relief_m::AbstractVector{<:Real})
    length(side_km) == length(mean_relief_m) ||
        throw(DimensionMismatch("side_km and mean_relief_m must have the same length"))
    all(>(0), side_km) && all(i -> side_km[i] < side_km[i + 1], 1:length(side_km) - 1) ||
        throw(ArgumentError("window sides must be positive and strictly increasing"))
    all(value -> isfinite(value) && value > 0, mean_relief_m) ||
        throw(ArgumentError("mean relief must be finite and positive"))
    log_relief_per_area = log.(Float64.(mean_relief_m) ./ Float64.(side_km) .^ 2)
    change = mean_change_point(log_relief_per_area)
    return (; change.index, side_km=Float64(side_km[change.index]), log_relief_per_area, change.gain)
end

"""
Relief of non-overlapping `k`x`k` blocks of a max/min pair of grids.

`max_grid` and `min_grid` hold each cell's maximum and minimum elevation, so the relief of a block is
the largest maximum minus the smallest minimum - exact, not an approximation, when the cells are whole
blocks of the finer DEM. Trailing rows and columns that do not fill a block are dropped, and a block
with any non-finite cell is `NaN`.
"""
function block_relief(max_grid::AbstractMatrix{<:Real}, min_grid::AbstractMatrix{<:Real}, k::Integer)
    size(max_grid) == size(min_grid) || throw(DimensionMismatch("max_grid and min_grid must have the same size"))
    k >= 1 || throw(ArgumentError("block size must be at least 1"))
    n_rows, n_cols = size(max_grid) .÷ k
    relief = fill(NaN, n_rows, n_cols)
    for block_col in 1:n_cols, block_row in 1:n_rows
        highest, lowest, complete = -Inf, Inf, true
        for j in (block_col - 1) * k + 1:block_col * k, i in (block_row - 1) * k + 1:block_row * k
            hi, lo = max_grid[i, j], min_grid[i, j]
            isfinite(hi) && isfinite(lo) || (complete = false; break)
            highest = max(highest, hi)
            lowest = min(lowest, lo)
        end
        complete && (relief[block_row, block_col] = highest - lowest)
    end
    return relief
end

"""Mean block relief and the number of complete blocks for each block size in `ks`."""
function mean_block_relief(max_grid::AbstractMatrix{<:Real}, min_grid::AbstractMatrix{<:Real}, ks)
    mean_relief = Float64[]
    n_blocks = Int[]
    for k in ks
        values = filter(isfinite, vec(block_relief(max_grid, min_grid, k)))
        push!(mean_relief, isempty(values) ? NaN : mean(values))
        push!(n_blocks, length(values))
    end
    return (; mean_relief, n_blocks)
end

"""
A regular grid from GDAL `XYZ` lines (`x y z`, one cell centre per line).

Returns `(; x, y, grid)` with `x` ascending, `y` descending (north up) and `grid[row, col]` the value
at `(x[col], y[row])`; the nodata value becomes `NaN`. A grid with missing cells is an error.
"""
function parse_xyz_grid(lines; nodata::Real=NODATA)
    xs, ys, zs = Float64[], Float64[], Float64[]
    for line in lines
        fields = split(line)
        isempty(fields) && continue
        length(fields) == 3 || throw(ArgumentError("expected `x y z`, got: $line"))
        push!(xs, parse(Float64, fields[1]))
        push!(ys, parse(Float64, fields[2]))
        push!(zs, parse(Float64, fields[3]))
    end
    x = sort(unique(xs))
    y = sort(unique(ys); rev=true)
    length(x) * length(y) == length(zs) ||
        throw(ArgumentError("XYZ lines do not form a complete $(length(y))x$(length(x)) grid"))
    column = Dict(value => j for (j, value) in enumerate(x))
    row = Dict(value => i for (i, value) in enumerate(y))
    grid = fill(NaN, length(y), length(x))
    for t in eachindex(zs)
        grid[row[ys[t]], column[xs[t]]] = zs[t] == nodata ? NaN : zs[t]
    end
    return (; x, y, grid)
end

"""
`value` with three decimals and never an exponent - Julia prints `3474245.66` as `3.47424566e6`.

Equal to `@sprintf("%.3f", value)` without depending on Printf, which the package does not list.
"""
function _fixed3(value::Real)
    # In BigFloat the product is exact, so a binary value just above a decimal tie rounds up, as Printf does.
    millis = Int(round(BigInt, big(Float64(value)) * 1000))
    sign = millis < 0 ? "-" : ""
    whole, fraction = divrem(abs(millis), 1000)
    return string(sign, whole, ".", lpad(fraction, 3, '0'))
end

"""Run `command` in `workdir`, returning stdout lines; stderr (GDAL/PROJ chatter) is discarded."""
_output_lines(command::Cmd, workdir::AbstractString) =
    readlines(pipeline(Cmd(command; dir=workdir); stderr=devnull))

"""Project `(a, b)` pairs between lon/lat and UTM 49N with `gmt mapproject` (`inverse` = UTM to lon/lat)."""
function _project(gmt::AbstractString, points::AbstractMatrix{<:Real}, workdir::AbstractString; inverse::Bool=false)
    path = joinpath(workdir, "points.txt")
    open(path, "w") do io
        for row in eachrow(points)
            # UTM northings print as `3.47e6` from Julia; whole metres are exact enough for 30 m cells.
            inverse ? println(io, round(Int, row[1]), ' ', round(Int, row[2])) : println(io, Float64(row[1]), ' ', Float64(row[2]))
        end
    end
    arguments = [gmt, "mapproject", path, "-J$UTM_PROJ", "-F"]
    inverse && push!(arguments, "-I")
    lines = _output_lines(Cmd(arguments), workdir)
    length(lines) == size(points, 1) || error("projected point count mismatch")
    return [parse(Float64, split(line)[k]) for line in lines, k in 1:2]
end

"""
Elevation range, and optionally mean slope, of a square window of side `side_m` centred on each point.

Uses `gmt grdinfo -R<window> -C`, which snaps the window outward to whole grid cells; the realised
cell count is returned. A window reaching nodata is an error rather than a smaller sample.
"""
function _window_stats(gmt, dem_path, slope_path, xy::AbstractMatrix{<:Real}, side_m::Real, workdir)
    rows = NamedTuple[]
    for (x, y) in eachrow(xy)
        half = side_m / 2
        region = "-R" * join((_fixed3(v) for v in (x - half, x + half, y - half, y + half)), "/")
        fields = split(only(_output_lines(Cmd([gmt, "grdinfo", dem_path, region, "-C"]), workdir)), '\t')
        z_min, z_max = parse(Float64, fields[6]), parse(Float64, fields[7])
        n_cells = parse(Int, fields[10]) * parse(Int, fields[11])
        z_min > NODATA / 2 || error("relief window at ($x, $y) reaches nodata")
        mean_slope = NaN
        if slope_path !== nothing
            slope_fields = split(only(_output_lines(Cmd([gmt, "grdinfo", slope_path, region, "-C", "-L2"]), workdir)), '\t')
            mean_slope = parse(Float64, slope_fields[12])
        end
        push!(rows, (; relief_m=z_max - z_min, mean_slope_deg=mean_slope, n_cells))
    end
    return rows
end

"""
Classify every gauge by the DEM relief around it, and write the diagnostics that justify the window.

1. `gdalwarp -r max` / `-r min` aggregate the 30 m DEM to `cell_m` cells over the gauge network's UTM
   bounding box. Cells are whole blocks of the source grid (`-tap`), so block relief is exact.
2. Area-wide mean relief over windows of `cell_m * ks` -> `optimal_relief_window` -> window `W`.
3. Gauge-centred relief at `W`, and at `W/2` and `2W` for sensitivity, from the 30 m DEM itself; mean
   slope at `W` from `slope_path`.
4. A region with fewer than `min_region_gauges` gauges is flagged `insufficient_gauges`.

Writes `output_path` (one row per gauge) and, under `diagnostics_dir`, `window_change_point.csv`,
`class_counts.csv`, `window_sensitivity.csv` and `relief_blocks.csv` (block relief at `W` with cell
centres in lon/lat, for mapping). Returns the gauge table and the chosen window.
"""
function prepare_station_landform(;
    dem_path::AbstractString, slope_path::AbstractString, terrain_path::AbstractString,
    output_path::AbstractString, diagnostics_dir::AbstractString,
    cell_m::Int=300, ks=1:30, min_region_gauges::Int=20,
)
    gmt, gdalwarp, gdal_translate = find_tool("gmt"), find_tool("gdalwarp"), find_tool("gdal_translate")
    terrain = CSV.read(terrain_path, DataFrame; types=Dict(:station_id => String))
    allunique(terrain.station_id) || error("duplicate station IDs in $terrain_path")
    lonlat = Float64.([terrain.lon terrain.lat])

    return mktempdir() do workdir
        xy = _project(gmt, lonlat, workdir)
        box = (
            floor(minimum(xy[:, 1]) / cell_m) * cell_m, floor(minimum(xy[:, 2]) / cell_m) * cell_m,
            ceil(maximum(xy[:, 1]) / cell_m) * cell_m, ceil(maximum(xy[:, 2]) / cell_m) * cell_m,
        )
        grids = Dict{String,Any}()
        for statistic in ("max", "min")
            tif = joinpath(workdir, "dem_$statistic.tif")
            xyz = joinpath(workdir, "dem_$statistic.xyz")
            println("  gdalwarp -r $statistic to $(cell_m) m ...")
            run(pipeline(Cmd(Cmd([
                gdalwarp, "-q", "-overwrite", "-r", statistic, "-tr", string(cell_m), string(cell_m), "-tap",
                "-te", string.(box)..., "-srcnodata", string(Int(NODATA)), "-dstnodata", string(Int(NODATA)),
                dem_path, tif,
            ]); dir=workdir); stderr=devnull))
            run(pipeline(Cmd(Cmd([gdal_translate, "-q", "-of", "XYZ", tif, xyz]); dir=workdir); stderr=devnull))
            grids[statistic] = parse_xyz_grid(eachline(xyz))
        end
        grids["max"].x == grids["min"].x && grids["max"].y == grids["min"].y ||
            error("max and min grids are not aligned")

        side_km = collect(ks) .* cell_m ./ 1000   # (k * 300) / 1000 is 2.7; k * (300 / 1000) is 2.6999999999999997
        area = mean_block_relief(grids["max"].grid, grids["min"].grid, ks)
        window = optimal_relief_window(side_km, area.mean_relief)
        side_m = window.side_km * 1000
        println("  mean change-point window: $(window.side_km) km")
        change_table = DataFrame(
            side_km=side_km, mean_relief_m=area.mean_relief, n_blocks=area.n_blocks,
            log_relief_per_area=window.log_relief_per_area, change_point_gain=window.gain,
            selected=eachindex(side_km) .== window.index,
        )

        windows = [("half", side_m / 2), ("selected", side_m), ("double", 2 * side_m)]
        stats = Dict{String,Vector{NamedTuple}}()
        for (label, side) in windows
            println("  gauge-centred relief, $(label) window $(side / 1000) km ...")
            stats[label] = _window_stats(gmt, dem_path, label == "selected" ? slope_path : nothing, xy, side, workdir)
        end

        relief(label) = [row.relief_m for row in stats[label]]
        stations = DataFrame(
            station_id=terrain.station_id, lon=lonlat[:, 1], lat=lonlat[:, 2],
            elevation_m=Float64.(terrain.elevation_m),
            window_side_km=fill(window.side_km, nrow(terrain)),
            relief_m=relief("selected"), window_cells=[row.n_cells for row in stats["selected"]],
            window_mean_slope_deg=[row.mean_slope_deg for row in stats["selected"]],
        )
        stations.relief_class = relief_class.(stations.relief_m)
        stations.relief_region = relief_region.(stations.relief_m)
        stations.altitude_class = altitude_class.(stations.elevation_m)
        stations.relief_m_half_window = relief("half")
        stations.relief_region_half_window = relief_region.(stations.relief_m_half_window)
        stations.relief_m_double_window = relief("double")
        stations.relief_region_double_window = relief_region.(stations.relief_m_double_window)
        stations.class_stable = (stations.relief_region .== stations.relief_region_half_window) .&
            (stations.relief_region .== stations.relief_region_double_window)
        region_count = Dict(region => count(==(region), stations.relief_region) for region in RELIEF_REGIONS)
        stations.region_gauges = [region_count[region] for region in stations.relief_region]
        stations.insufficient_gauges = stations.region_gauges .< min_region_gauges

        counts = DataFrame(scheme=String[], class=String[], n_gauges=Int[], mean_relief_m=Float64[],
            mean_window_slope_deg=Float64[], mean_elevation_m=Float64[], insufficient_gauges=Bool[])
        for (scheme, column, classes) in (("relief_class", :relief_class, RELIEF_CLASSES), ("relief_region", :relief_region, RELIEF_REGIONS))
            for class in classes
                members = stations[stations[!, column] .== class, :]
                n = nrow(members)
                average(v) = n == 0 ? NaN : mean(v)
                push!(counts, (scheme, class, n, average(members.relief_m), average(members.window_mean_slope_deg),
                    average(members.elevation_m), n < min_region_gauges))
            end
        end

        sensitivity = DataFrame(window=String[], side_km=Float64[], class=String[], n_gauges=Int[])
        for (label, side) in windows
            regions = relief_region.(relief(label))
            classes = relief_class.(relief(label))
            for class in RELIEF_CLASSES
                push!(sensitivity, (label, side / 1000, class, count(==(class), classes)))
            end
            for region in RELIEF_REGIONS
                push!(sensitivity, (label, side / 1000, region, count(==(region), regions)))
            end
        end

        k = window.index
        blocks = block_relief(grids["max"].grid, grids["min"].grid, ks[k])
        block_side = ks[k] * cell_m
        centres = [(grids["max"].x[1] - cell_m / 2 + (j - 0.5) * block_side, grids["max"].y[1] + cell_m / 2 - (i - 0.5) * block_side)
            for i in axes(blocks, 1), j in axes(blocks, 2) if isfinite(blocks[i, j])]
        block_values = [blocks[i, j] for i in axes(blocks, 1), j in axes(blocks, 2) if isfinite(blocks[i, j])]
        block_lonlat = _project(gmt, [first.(centres) last.(centres)], workdir; inverse=true)
        block_table = DataFrame(lon=block_lonlat[:, 1], lat=block_lonlat[:, 2], x_utm=first.(centres), y_utm=last.(centres),
            side_km=fill(window.side_km, length(block_values)), relief_m=block_values)
        block_table.relief_class = relief_class.(block_table.relief_m)
        block_table.relief_region = relief_region.(block_table.relief_m)

        write_csv_atomic(output_path, stations)
        write_csv_atomic(joinpath(diagnostics_dir, "window_change_point.csv"), change_table)
        write_csv_atomic(joinpath(diagnostics_dir, "class_counts.csv"), counts)
        write_csv_atomic(joinpath(diagnostics_dir, "window_sensitivity.csv"), sensitivity)
        write_csv_atomic(joinpath(diagnostics_dir, "relief_blocks.csv"), block_table)
        (; stations, window_side_km=window.side_km, counts, change_table)
    end
end

end # module
