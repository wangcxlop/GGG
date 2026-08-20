module TerrainFeatures

using CSV
using DataFrames

export aspect_components, extract_station_terrain

const STUDY_BOUNDS = (west=109.4, east=111.6, south=31.2, north=33.4)
const TARGET_CRS = "EPSG:32649"
const NODATA = -9999.0

function aspect_components(aspect_deg::Real, slope_deg::Real; flat_threshold::Real=0.1)
    aspect = Float64(aspect_deg)
    slope = Float64(slope_deg)
    if !isfinite(aspect) || !isfinite(slope) || slope < flat_threshold
        return (0.0, 0.0)
    end
    return (sind(aspect), cosd(aspect))
end

function find_tool(name::AbstractString)
    path = Sys.which(name)
    path === nothing && error("Missing GDAL command: $name")
    return path
end

function expected_tiles(raw_dir::AbstractString)
    paths = [
        joinpath(
            raw_dir,
            "Copernicus_DSM_COG_10_N$(lpad(lat, 2, '0'))_00_" *
            "E$(lpad(lon, 3, '0'))_00_DEM.tif",
        )
        for lat in 31:33 for lon in 109:111
    ]
    missing = filter(path -> !isfile(path), paths)
    isempty(missing) || error("Missing $(length(missing)) study-area DEM tiles")
    return paths
end

function run_dem_derivatives(raw_dir::AbstractString, output_dir::AbstractString)
    mkpath(output_dir)
    vrt = joinpath(output_dir, "copernicus_glo30_model_mosaic.vrt")
    projected = joinpath(output_dir, "copernicus_glo30_utm49n_30m.tif")
    slope = joinpath(output_dir, "copernicus_glo30_slope_deg.tif")
    aspect = joinpath(output_dir, "copernicus_glo30_aspect_deg.tif")

    expected_tiles(raw_dir)
    if all(path -> isfile(path) && filesize(path) > 0, (projected, slope, aspect))
        return (; projected, slope, aspect)
    end

    gdalbuildvrt = find_tool("gdalbuildvrt")
    gdalwarp = find_tool("gdalwarp")
    gdaldem = find_tool("gdaldem")
    tiles = expected_tiles(raw_dir)

    run(Cmd(vcat([gdalbuildvrt, "-overwrite", vrt], tiles)))
    run(Cmd([
        gdalwarp, "-overwrite",
        "-t_srs", TARGET_CRS,
        "-te_srs", "EPSG:4326",
        "-te", string(STUDY_BOUNDS.west), string(STUDY_BOUNDS.south),
        string(STUDY_BOUNDS.east), string(STUDY_BOUNDS.north),
        "-tr", "30", "30", "-tap", "-r", "bilinear",
        "-dstnodata", string(Int(NODATA)),
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3",
        vrt, projected,
    ]))
    run(Cmd([
        gdaldem, "slope", projected, slope,
        "-compute_edges", "-of", "GTiff",
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3",
    ]))
    run(Cmd([
        gdaldem, "aspect", projected, aspect,
        "-compute_edges", "-zero_for_flat", "-of", "GTiff",
        "-co", "TILED=YES", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3",
    ]))
    return (; projected, slope, aspect)
end

function command_lines(command::Cmd)
    return mktemp() do _, io
        run(pipeline(command; stdout=io))
        flush(io)
        seekstart(io)
        readlines(io)
    end
end

function projected_station_points(stations::DataFrame)
    gmt = find_tool("gmt")
    return mktemp() do path, io
        for row in eachrow(stations)
            println(io, Float64(row.lon), ' ', Float64(row.lat))
        end
        close(io)
        command = Cmd([
            gmt, "mapproject", path,
            "-J+proj=utm +zone=49 +datum=WGS84 +units=m", "-F",
        ])
        lines = command_lines(command)
        length(lines) == nrow(stations) || error("Projected station count mismatch")
        points = Tuple{Float64, Float64}[]
        sizehint!(points, length(lines))
        for line in lines
            fields = split(strip(line))
            length(fields) >= 2 || error("Invalid projected coordinate: $line")
            x = parse(Float64, fields[1])
            y = parse(Float64, fields[2])
            push!(points, (x, y))
        end
        points
    end
end

function sample_raster_batch(
    path::AbstractString,
    points::Vector{Tuple{Float64, Float64}};
    resampling::String,
)
    interpolation = resampling == "bilinear" ? "-nl" :
                    resampling == "nearest" ? "-nn" :
                    error("Unsupported GMT resampling: $resampling")
    gmt = find_tool("gmt")
    return mktemp() do point_path, io
        for (x, y) in points
            println(io, x, ' ', y)
        end
        close(io)
        lines = command_lines(Cmd([
            gmt, "grdtrack", point_path, "-G$path", interpolation,
        ]))
        length(lines) == length(points) || error("Raster sample count mismatch")
        values = Float64[]
        sizehint!(values, length(lines))
        for line in lines
            fields = split(strip(line))
            isempty(fields) && error("Empty raster sample from $(basename(path))")
            value = parse(Float64, fields[end])
            isfinite(value) && value != NODATA || error("Invalid sample from $(basename(path))")
            push!(values, value)
        end
        values
    end
end

function write_csv_atomic(path::AbstractString, table)
    mkpath(dirname(path))
    temporary = string(path, ".tmp-", getpid())
    try
        CSV.write(temporary, table)
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

function extract_station_terrain(
    station_meta_path::AbstractString,
    raw_dir::AbstractString,
    output_dir::AbstractString,
    output_path::AbstractString;
    expected_station_count::Int=237,
)
    stations = CSV.read(station_meta_path, DataFrame; types=Dict(:station_id => String))
    required = [:station_id, :lon, :lat]
    missing_columns = setdiff(required, propertynames(stations))
    isempty(missing_columns) || error("Station metadata missing: $(join(missing_columns, ", "))")
    nrow(stations) == expected_station_count ||
        error("Expected $expected_station_count stations, found $(nrow(stations))")
    allunique(stations.station_id) || error("Duplicate station IDs")

    rasters = run_dem_derivatives(raw_dir, output_dir)
    points = projected_station_points(stations)
    elevation = sample_raster_batch(rasters.projected, points; resampling="bilinear")
    slope_deg = sample_raster_batch(rasters.slope, points; resampling="bilinear")
    aspect_deg = sample_raster_batch(rasters.aspect, points; resampling="nearest")
    aspect_sin = Float64[]
    aspect_cos = Float64[]
    sizehint!(aspect_sin, nrow(stations))
    sizehint!(aspect_cos, nrow(stations))
    for (slope, aspect) in zip(slope_deg, aspect_deg)
        sin_aspect, cos_aspect = aspect_components(aspect, slope)
        push!(aspect_sin, sin_aspect)
        push!(aspect_cos, cos_aspect)
    end

    result = DataFrame(
        station_id=stations.station_id,
        lon=Float64.(stations.lon),
        lat=Float64.(stations.lat),
        elevation_m=elevation,
        slope_deg=slope_deg,
        aspect_deg=aspect_deg,
        aspect_sin=aspect_sin,
        aspect_cos=aspect_cos,
    )
    all(isfinite, Matrix(result[:, Not(:station_id)])) || error("Non-finite terrain feature")
    all(value -> 0 <= value <= 90, result.slope_deg) || error("Slope outside 0-90 degrees")
    all(value -> 0 <= value <= 360, result.aspect_deg) || error("Aspect outside 0-360 degrees")
    write_csv_atomic(output_path, result)
    return (; table=result, rasters...)
end

end
