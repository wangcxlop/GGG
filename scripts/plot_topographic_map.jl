#!/usr/bin/env julia

using CSV
using DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW_DIR = joinpath(ROOT, "data", "raw", "copernicus_glo30")
const PROCESSED_DIR = joinpath(ROOT, "data", "processed", "dem")
const OUTPUT_DIR = joinpath(ROOT, "output", "terrain_map")
const STATION_FILE = joinpath(ROOT, "data", "processed", "study_area", "station_meta.csv")

const WEST, EAST = 109.4, 111.6
const SOUTH, NORTH = 31.2, 33.4
const DISPLAY_WIDTH = 3000

const VRT_FILE = joinpath(PROCESSED_DIR, "copernicus_glo30_mosaic.vrt")
const DISPLAY_DEM = joinpath(PROCESSED_DIR, "copernicus_glo30_study_area_map.tif")
const CPT_FILE = joinpath(OUTPUT_DIR, "elevation_green_yellow_red.cpt")
const MAP_NAME = "study_area_topographic_map_with_stations"
const MAP_FILE = joinpath(OUTPUT_DIR, "$MAP_NAME.png")

function find_gmt_bin()
    candidates = String[]
    haskey(ENV, "GMT_BIN") && push!(candidates, ENV["GMT_BIN"])
    Sys.iswindows() && push!(candidates, raw"C:\programs\gmt6\bin")

    gmt = Sys.which("gmt")
    gmt !== nothing && push!(candidates, dirname(gmt))

    for dir in candidates
        exe = joinpath(dir, Sys.iswindows() ? "gmt.exe" : "gmt")
        isfile(exe) && return normpath(dir)
    end

    error(
        "未找到 GMT。请先运行 `julia --project=. -e \"using Pkg; Pkg.build(\\\"GMT\\\")\"`，" *
        "或通过环境变量 GMT_BIN 指定 GMT 的 bin 目录。",
    )
end

function tool(bin::AbstractString, name::AbstractString)
    exe = joinpath(bin, Sys.iswindows() ? "$name.exe" : name)
    isfile(exe) || error("缺少 GMT/GDAL 工具：$exe")
    return exe
end

function expected_tiles()
    return [
        "Copernicus_DSM_COG_10_N$(lpad(lat, 2, '0'))_00_" *
        "E$(lpad(lon, 3, '0'))_00_DEM.tif"
        for lat in 31:33 for lon in 109:111
    ]
end

function validate_inputs()
    expected = expected_tiles()
    missing = filter(name -> !isfile(joinpath(RAW_DIR, name)), expected)
    isempty(missing) || error("缺少 $(length(missing)) 个 DEM 瓦片：$(join(missing, ", "))")

    parts = filter(name -> endswith(name, ".part"), readdir(RAW_DIR))
    isempty(parts) || error("原始数据目录仍有未完成下载：$(join(parts, ", "))")
    return joinpath.(RAW_DIR, expected)
end

function read_stations()
    isfile(STATION_FILE) || error("缺少站点文件：$STATION_FILE")
    stations = CSV.read(STATION_FILE, DataFrame)
    required = ["station_id", "lon", "lat"]
    missing_columns = setdiff(required, names(stations))
    isempty(missing_columns) || error("站点文件缺少字段：$(join(missing_columns, ", "))")

    any(ismissing, stations.lon) && error("站点经度存在缺失值")
    any(ismissing, stations.lat) && error("站点纬度存在缺失值")
    all(isfinite, stations.lon) || error("站点经度存在非有限值")
    all(isfinite, stations.lat) || error("站点纬度存在非有限值")

    inside =
        (WEST .<= stations.lon .<= EAST) .&
        (SOUTH .<= stations.lat .<= NORTH)
    selected = stations[inside, :]
    isempty(selected) && error("研究区范围内没有站点")
    return selected
end

function build_display_dem(bin::AbstractString, tiles::Vector{String})
    gdalbuildvrt = tool(bin, "gdalbuildvrt")
    gdal_translate = tool(bin, "gdal_translate")

    mkpath(PROCESSED_DIR)
    println("构建9瓦片虚拟镶嵌：$(VRT_FILE)")
    run(Cmd(vcat([gdalbuildvrt, "-overwrite", VRT_FILE], tiles)))

    println("裁剪研究区并生成绘图分辨率栅格：$(DISPLAY_DEM)")
    args = [
        gdal_translate,
        "-projwin", string(WEST), string(NORTH), string(EAST), string(SOUTH),
        "-outsize", string(DISPLAY_WIDTH), "0",
        "-r", "average",
        "-of", "GTiff",
        "-co", "TILED=YES",
        "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=3",
        "-co", "BIGTIFF=IF_SAFER",
        VRT_FILE,
        DISPLAY_DEM,
    ]
    run(Cmd(args))
    return DISPLAY_DEM
end

function elevation_range(bin::AbstractString, dem::AbstractString)
    gmt = tool(bin, "gmt")
    line = mktemp() do _, io
        run(pipeline(Cmd([gmt, "grdinfo", "-C", dem]); stdout=io))
        flush(io)
        seekstart(io)
        read(io, String)
    end
    fields = split(strip(line))
    length(fields) >= 7 || error("无法解析 grdinfo 输出：$line")
    return parse(Float64, fields[6]), parse(Float64, fields[7])
end

function write_elevation_cpt(path::AbstractString, zmin::Real, zmax::Real)
    low = min(0.0, floor(zmin / 100) * 100)
    high = max(2000.0, ceil(zmax / 500) * 500)

    # 常规地形色带：低地绿，中海拔黄/橙，高海拔红。
    colors = [
        (low,  "28/120/64"),
        (200.0, "101/190/90"),
        (500.0, "198/219/107"),
        (1000.0, "246/198/85"),
        (1500.0, "229/134/59"),
    ]
    high > 2000 && push!(colors, (2000.0, "202/70/52"))
    push!(colors, (high, "128/0/38"))

    open(path, "w") do io
        for i in 1:(length(colors) - 1)
            z1, c1 = colors[i]
            z2, c2 = colors[i + 1]
            println(io, "$(z1) $c1 $(z2) $c2")
        end
        println(io, "B 28/120/64")
        println(io, "F 128/0/38")
        println(io, "N 255/255/255")
    end
    return low, high
end

function write_station_table(path::AbstractString, stations::DataFrame)
    open(path, "w") do io
        for row in eachrow(stations)
            println(io, "$(row.lon)\t$(row.lat)")
        end
    end
    return path
end

function plot_map(
    bin::AbstractString,
    dem::AbstractString,
    cpt::AbstractString,
    stations::DataFrame,
)
    gmt = tool(bin, "gmt")
    region = "$(WEST)/$(EAST)/$(SOUTH)/$(NORTH)"

    mkpath(OUTPUT_DIR)
    mktempdir() do gmt_userdir
        station_table = write_station_table(joinpath(gmt_userdir, "stations.tsv"), stations)
        cd(OUTPUT_DIR) do
            withenv(
                "GMT_USERDIR" => gmt_userdir,
                "GMT_SESSION_NAME" => string(getpid()),
            ) do
                run(Cmd([gmt, "begin", MAP_NAME, "png"]))
                completed = false
                try
                    run(Cmd([
                        gmt, "set",
                        "MAP_FRAME_TYPE", "plain",
                        "MAP_FRAME_PEN", "1p,40/40/40",
                        "FONT_TITLE", "16p,Helvetica-Bold,35/35/35",
                        "FONT_LABEL", "11p,Helvetica,45/45/45",
                        "FONT_ANNOT_PRIMARY", "9p,Helvetica,55/55/55",
                        "MAP_TICK_PEN_PRIMARY", "0.8p,45/45/45",
                    ]))
                    run(Cmd([
                        gmt, "grdimage", dem,
                        "-R$region",
                        "-JM22c",
                        "-C$cpt",
                        "-E300",
                        "-Bxa0.5f0.25+lLongitude (deg E)",
                        "-Bya0.5f0.25+lLatitude (deg N)",
                        "-BWSen+tStudy Area Topography and Stations",
                    ]))
                    # 白色小光晕保证站点在绿色、黄色和红色地形背景上均清晰可见。
                    run(Cmd([
                        gmt, "plot", station_table,
                        "-R$region", "-JM22c", "-Sc0.12c", "-Gwhite", "-W0p",
                    ]))
                    run(Cmd([
                        gmt, "plot", station_table,
                        "-R$region", "-JM22c", "-Sc0.07c", "-G20/20/20", "-W0p",
                    ]))
                    run(Cmd([
                        gmt, "colorbar",
                        "-C$cpt",
                        "-DJBC+w14c/0.45c+o0c/0.5c+h",
                        "-Bxa500f250+lElevation (m)",
                    ]))
                    run(Cmd([gmt, "end"]))
                    completed = true
                finally
                    completed || run(Cmd([gmt, "end"]))
                end
            end
        end
    end

    isfile(MAP_FILE) || error("GMT 未生成预期图片：$MAP_FILE")
    return MAP_FILE
end

function main()
    bin = find_gmt_bin()
    tiles = validate_inputs()
    stations = read_stations()
    dem = build_display_dem(bin, tiles)
    zmin, zmax = elevation_range(bin, dem)

    mkpath(OUTPUT_DIR)
    color_min, color_max = write_elevation_cpt(CPT_FILE, zmin, zmax)
    map = plot_map(bin, dem, CPT_FILE, stations)

    println("研究区高程范围：$(round(zmin; digits=1))–$(round(zmax; digits=1)) m")
    println("色带范围：$(color_min)–$(color_max) m")
    println("已叠加站点：$(nrow(stations)) 个")
    println("地形图已生成：$map")
    return map
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
