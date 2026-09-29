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
