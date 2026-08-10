#!/usr/bin/env julia
module FY4BPreprocessing

using NCDatasets, CSV, DataFrames, Dates

export aggregate_fy4b_hourly, build_hourly_qc, find_nc_files, parse_nc_filename

# =========================
# FY4B Aggregation: 15min NetCDF -> Hourly CSV
# =========================

const ROOT = normpath(joinpath(@__DIR__, ".."))
const FY4B_DATA_DIR = joinpath(ROOT, "data", "FY4B")
const STATION_META_FILE = joinpath(ROOT, "data", "hubei_station_meta.csv")
const OUTPUT_FILE = joinpath(ROOT, "data", "processed", "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv")
const QC_FILE = joinpath(ROOT, "output", "input_audit", "fy4b_hourly_qc_2022_2025_JunSep_strict_navcorrected.csv")

function write_csv_atomic(path, df)
    mkpath(dirname(path))
    temp_path = string(path, ".tmp-", getpid())
    CSV.write(temp_path, df)
    mv(temp_path, path; force=true)
end

# Satellite parameters for FY4B (from NC files)
const DEFAULT_SAT_LON = 133.0   # Fallback satellite longitude (degrees)
const SAT_HEIGHT = 42164.0e3    # Satellite height from Earth center (m)
const EARTH_RADIUS = 6378.137e3 # Earth equatorial radius (m)
const GRID_SIZE = 2748          # Grid dimensions
const RESOLUTION = 4000.0       # 4km at nadir (m)
const STATION_TIME_OFFSET = Hour(8) # FY4B filenames are UTC; station timestamps are Beijing time

# Computed constants
const SAT_DISTANCE = SAT_HEIGHT - EARTH_RADIUS  # Distance from surface to satellite
const SCALE_FACTOR = (GRID_SIZE - 1) / 2.0      # For scaling scan angles to grid indices

"""
    latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON)

Convert latitude/longitude to scan angles (x, y in radians) for geostationary projection.
Uses standard geostationary projection formula (similar to GOES-R ABI).
Returns: (scan_x, scan_y) in radians, or (NaN, NaN) if not visible.
"""
function latlon_to_scan_angles(lat, lon, sat_lon=DEFAULT_SAT_LON)
    # Convert to radians
    φ = deg2rad(lat)
    λ = deg2rad(lon - sat_lon)  # Relative longitude

    # WGS84 parameters
    re = EARTH_RADIUS           # Equatorial radius
    rp = 6356.75231414e3        # Polar radius (m)
    f = (re - rp) / re          # Flattening
    e2 = f * (2 - f)            # Eccentricity squared

    # Geocentric latitude (accounting for Earth flattening)
    φ_c = atan((1 - e2) * tan(φ))

    # Distance from Earth center to surface point
    rc = rp / sqrt(1 - e2 * cos(φ_c)^2)

    # Cartesian coordinates of surface point (Earth-centered, x-axis through 0°E)
    xs = rc * cos(φ_c) * cos(λ)
    ys = rc * cos(φ_c) * sin(λ)
    zs = rc * sin(φ_c)

    # Visibility check: point must be on Earth side facing the satellite
    # The satellite is at (SAT_HEIGHT, 0, 0) in this coordinate system
    # Point is visible if it's not behind the Earth's limb
    dot_product = xs * (xs - SAT_HEIGHT) + ys^2 + zs^2
    if dot_product > 0
        # Check if behind Earth
        # Actually, we need to check if the line from satellite to point intersects Earth
        # Simplified: check the angle from satellite
        d2 = (SAT_HEIGHT - xs)^2 + ys^2 + zs^2
        cos_angle = (SAT_HEIGHT - xs) / sqrt(d2)
        if cos_angle < 0.156  # cos(81°) - maximum view angle
            return (NaN, NaN)
        end
    end

    # Scan angles from satellite perspective
    # x (E-W, west negative): angle in the equatorial plane
    scan_x = atan(ys / (SAT_HEIGHT - xs))

    # y (N-S, north positive): angle perpendicular to equatorial plane
    scan_y = atan(zs / sqrt((SAT_HEIGHT - xs)^2 + ys^2))

    return (scan_x, scan_y)
end

"""
    scan_angles_to_grid(scan_x, scan_y)

Convert scan angles (radians) to FY4B grid indices (0-based).
The grid is centered at the satellite sub-point (nadir).
"""
function scan_angles_to_grid(scan_x, scan_y)
    if isnan(scan_x) || isnan(scan_y)
        return (NaN, NaN)
    end

    # Convert scan angle to linear distance at nadir
    # At nadir: angle = distance / (H - R)
    # scan_angle = distance / SAT_DISTANCE
    # distance = scan_angle * SAT_DISTANCE

    # Grid spacing in radians (at nadir, RESOLUTION meters = scan_angle * SAT_DISTANCE)
    grid_spacing_rad = RESOLUTION / SAT_DISTANCE

    # Grid center
    center = (GRID_SIZE - 1) / 2.0

    # Convert scan angles to grid indices
    # Negative scan_x (west of nadir) -> smaller x index
    # Positive scan_y (north of nadir) -> smaller y index (image coordinates)
    x_idx = center + scan_x / grid_spacing_rad
    y_idx = center - scan_y / grid_spacing_rad

    return (x_idx, y_idx)
end

"""
    parse_satellite_lon(filename)

Parse FY4B satellite subpoint longitude from the filename.
Pattern segment: N_DISK_1330E or N_DISK_1050E.
Returns degrees east as Float64.
"""
function parse_satellite_lon(filename)
    m = match(r"N_DISK_(\d{4})E", filename)
    if m === nothing
        @warn "Could not parse satellite longitude from filename; using fallback" filename fallback=DEFAULT_SAT_LON
        return DEFAULT_SAT_LON
    end
    return parse(Int, m.captures[1]) / 10.0
end

"""
    parse_nc_filename(filename)

Parse FY4B NetCDF filename to extract start/end datetime and satellite longitude.
Pattern: FY4B-_AGRI--_N_DISK_1330E_L2-_QPE-_MULT_NOM_YYYYMMDDhhmmss_YYYYMMDDhhmmss_4000M_V0001.NC
Returns: (start_datetime, end_datetime, satellite_longitude)
"""
function parse_nc_filename(filename)
    m = match(r"NOM_(\d{14})_(\d{14})_", filename)
    if m === nothing
        return nothing
    end
    start_str, end_str = m.captures
    start_dt = DateTime(start_str, dateformat"yyyymmddHHMMSS")
    end_dt = DateTime(end_str, dateformat"yyyymmddHHMMSS")
    return (start_dt, end_dt, parse_satellite_lon(filename))
end

"""
    is_15min_file(filename)

Check if file represents a 15-minute aggregation window.
True 15min files have a time span of exactly 15 minutes (0:14:59).
"""
function is_15min_file(filename)
    dt_info = parse_nc_filename(filename)
    if dt_info === nothing
        return false
    end
    start_dt, end_dt, _ = dt_info
    return second(start_dt) == 0 &&
           minute(start_dt) in (0, 15, 30, 45) &&
           end_dt == start_dt + Minute(15) - Second(1)
end

function dataset_satellite_lon(ds, filename)
    if haskey(ds, "nominal_satellite_subpoint_lon")
        value = ds["nominal_satellite_subpoint_lon"][1]
        if !ismissing(value) && isfinite(value)
            return Float64(value)
        end
    end
    return parse_satellite_lon(filename)
end

aligned_hour_end(hour_start::DateTime) = hour_start + Hour(1) + STATION_TIME_OFFSET

"""
    latlon_to_xy(lat, lon, sat_lon=DEFAULT_SAT_LON)

Convert latitude/longitude to FY4B fixed grid x/y coordinates.
Uses geostationary projection formulas.
Returns: (x, y) in grid coordinates (0-based index)
"""
function latlon_to_xy(lat, lon, sat_lon=DEFAULT_SAT_LON)
    scan_x, scan_y = latlon_to_scan_angles(lat, lon, sat_lon)
    return scan_angles_to_grid(scan_x, scan_y)
end

"""
    station_grid_coords(stations_df, sat_lon, coord_cache)

Compute and cache station grid coordinates for a satellite subpoint longitude.
The default aggregation keeps the 133.0E projection used by the source metadata
and the published validation run. Filename-based 105.0E/133.0E re-extraction is
handled as a separate audit candidate in aggregate_fy4b_python.py.
"""
function station_grid_coords(stations_df, sat_lon, coord_cache)
    if !haskey(coord_cache, sat_lon)
        x_coords = Float64[]
        y_coords = Float64[]
        sizehint!(x_coords, nrow(stations_df))
        sizehint!(y_coords, nrow(stations_df))

        for row in eachrow(stations_df)
            x, y = latlon_to_xy(row.lat, row.lon, sat_lon)
            push!(x_coords, x)
            push!(y_coords, y)
        end

        coord_cache[sat_lon] = (x_coords, y_coords)
    end

    return coord_cache[sat_lon]
end

"""
    load_stations(meta_file)

Load station metadata.
Returns DataFrame with station_id, lon, lat.
"""
function load_stations(meta_file)
    df = CSV.read(meta_file, DataFrame)
    # Ensure correct column names
    rename!(df, lowercase.(names(df)))
    return df
end

"""
    extract_precipitation(nc_file, stations_df, coord_cache)

Extract precipitation values from NetCDF file at station locations.
Returns: Vector of precipitation values (mm) matching stations_df order.
"""
function extract_precipitation(nc_file, stations_df, coord_cache)
    ds = Dataset(nc_file)
    try
        precip = ds["Precipitation"]
        fill_val = Float32(65534.0)  # From NC attributes
        quality = haskey(ds, "DQF") ? ds["DQF"] : nothing

        sat_lon = dataset_satellite_lon(ds, basename(nc_file))
        x_coords, y_coords = station_grid_coords(stations_df, sat_lon, coord_cache)

        values = Float64[]
        sizehint!(values, nrow(stations_df))

        for i in 1:nrow(stations_df)
            x_idx = round(Int, x_coords[i]) + 1  # Convert 0-based to 1-based
            y_idx = round(Int, y_coords[i]) + 1

            # Bounds checking
            if x_idx < 1 || x_idx > GRID_SIZE || y_idx < 1 || y_idx > GRID_SIZE
                push!(values, NaN)
                continue
            end

            # NCDatasets exposes the NetCDF y/x dimensions in Julia indexing
            # order as x/y. Navigation audits against gauges confirm [x, y].
            val = precip[x_idx, y_idx]
            quality_value = quality === nothing ? 0 : quality[x_idx, y_idx]

            # Check for fill value or invalid
            if ismissing(val) || ismissing(quality_value) ||
               val == fill_val || val < 0 || val > 30 || Int(quality_value) >= 2
                push!(values, NaN)
            else
                push!(values, Float64(val))
            end
        end

        return values
    finally
        close(ds)
    end
end

"""
    aggregate_hourly(file_groups::Dict)

Aggregate 15min precipitation files to hourly sums.
file_groups: Dict mapping hour_start_datetime -> list of 4 NC file paths
Returns: DataFrame with time, station_id, fy4b_mm_h
"""
function aggregate_hourly(file_groups, stations_df)
    results = DataFrame(time=DateTime[], station_id=String[], fy4b_mm_h=Float64[])
    coord_cache = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()

    hours = sort(collect(keys(file_groups)))

    for (hour_index, hour_start) in enumerate(hours)
        files = file_groups[hour_start]

        starts = [first(parse_nc_filename(basename(file))) for file in files]
        minutes = sort(minute.(starts))
        @assert length(files) == 4 && minutes == [0, 15, 30, 45] "Non-complete hour reached aggregation: $hour_start"

        # The 15-minute QPE fields are rainfall rates. Their four-value mean is
        # numerically the one-hour accumulation in mm (sum(rate * 0.25 h)).
        hourly_sum = zeros(Float64, nrow(stations_df))
        has_missing = falses(nrow(stations_df))

        for nc_file in files
            vals = try
                extract_precipitation(nc_file, stations_df, coord_cache)
            catch error_value
                error("Failed to extract $(nc_file): $(sprint(showerror, error_value))")
            end
            for i in 1:length(vals)
                if isnan(vals[i])
                    has_missing[i] = true
                else
                    hourly_sum[i] += vals[i]
                end
            end
        end

        # Mark as NaN if any 15min value was missing
        for i in 1:length(hourly_sum)
            if has_missing[i]
                hourly_sum[i] = NaN
            else
                hourly_sum[i] /= 4.0
            end
        end

        # Add to results
        hour_end = aligned_hour_end(hour_start)
        time_str = Dates.format(hour_end, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"

        for (i, row) in enumerate(eachrow(stations_df))
            push!(results, (DateTime(replace(time_str, "Z" => "")), string(row.station_id), hourly_sum[i]))
        end
        if hour_index % 500 == 0 || hour_index == length(hours)
            println("  Aggregated $hour_index / $(length(hours)) strict hours")
        end
    end

    return results
end

"""
    find_nc_files(data_dir)

Find all NC files and group 15min files by hour.
Returns: Dict mapping hour_start -> list of NC file paths
"""
function find_nc_files(data_dir; years=2022:2025, months=6:9)
    file_groups = Dict{DateTime, Vector{String}}()
    candidates = Dict{Tuple{DateTime, Float64}, Vector{String}}()

    # Downloads are stored below per-year batch directories, not year/month.
    for (dir, _, files) in walkdir(data_dir), file in files
        endswith(uppercase(file), ".NC") || continue
        is_15min_file(file) || continue
        start_dt, _, sat_lon = parse_nc_filename(file)
        year(start_dt) in years || continue
        month(start_dt) in months || continue
        push!(get!(candidates, (start_dt, sat_lon), String[]), joinpath(dir, file))
    end

    duplicate_file_count = 0
    for ((start_dt, _), paths) in candidates
        sort!(paths)
        duplicate_file_count += length(paths) - 1
        hour_start = DateTime(year(start_dt), month(start_dt), day(start_dt), hour(start_dt))
        push!(get!(file_groups, hour_start, String[]), first(paths))
    end

    # Sort files within each hour group
    for (hour, files) in file_groups
        sort!(files)
    end
    duplicate_file_count > 0 && @info "Deduplicated repeated FY4B 15-minute files" duplicate_file_count

    return file_groups
end

"""Open every selected NetCDF and return the paths that fail structural validation."""
function find_unreadable_files(file_groups)
    files = sort(unique(vcat(values(file_groups)...)))
    unreadable = Dict{String, String}()
    for (index, file) in enumerate(files)
        try
            ds = Dataset(file)
            try
                haskey(ds, "Precipitation") || error("missing Precipitation variable")
                precipitation_size = size(ds["Precipitation"])
                precipitation_size == (GRID_SIZE, GRID_SIZE) ||
                    error("unexpected Precipitation size $precipitation_size")
            finally
                close(ds)
            end
        catch error_value
            unreadable[file] = sprint(showerror, error_value)
        end
        if index % 5_000 == 0 || index == length(files)
            println("  Validated $index / $(length(files)) NetCDF files")
        end
    end
    return unreadable
end

"""Return hour starts for exact native one-hour FY4B products in the raw download tree."""
function find_native_hourly_times(data_dir)
    native_hours = Set{DateTime}()
    for (dir, _, files) in walkdir(data_dir), file in files
        endswith(file, ".NC") || continue
        dt_info = parse_nc_filename(file)
        dt_info === nothing && continue
        start_dt, end_dt, _ = dt_info
        if minute(start_dt) == 0 && second(start_dt) == 0 &&
           end_dt == start_dt + Hour(1) - Second(1)
            push!(native_hours, start_dt)
        end
    end
    return native_hours
end

"""Build one auditable row for every expected hour in June--September 2022--2025."""
function build_hourly_qc(
    file_groups, native_hours;
    unreadable_files=Dict{String, String}(), years=2022:2025,
    start_month::Int=6, stop_month::Int=10,
)
    rows = NamedTuple[]
    complete_groups = Dict{DateTime, Vector{String}}()
    expected_minutes = [0, 15, 30, 45]

    for year_value in years
        hour_start = DateTime(year_value, start_month, 1)
        stop = DateTime(year_value, stop_month, 1)
        while hour_start < stop
            files = sort(get(file_groups, hour_start, String[]))
            readable_files = [file for file in files if !haskey(unreadable_files, file)]
            unreadable_hour_files = [file for file in files if haskey(unreadable_files, file)]
            starts = DateTime[]
            all_starts = DateTime[]
            for file in readable_files
                info = parse_nc_filename(basename(file))
                info === nothing || push!(starts, info[1])
            end
            for file in files
                info = parse_nc_filename(basename(file))
                info === nothing || push!(all_starts, info[1])
            end
            present_minutes = sort(unique(minute.(starts)))
            missing_minutes = setdiff(expected_minutes, present_minutes)
            duplicate_count = length(all_starts) - length(unique(minute.(all_starts)))
            complete = length(readable_files) == 4 && present_minutes == expected_minutes && duplicate_count == 0
            if complete
                complete_groups[hour_start] = readable_files
            end

            present_text = join([lpad(string(value), 2, '0') for value in present_minutes], ";")
            missing_text = join([lpad(string(value), 2, '0') for value in missing_minutes], ";")
            source_text = join(basename.(files), ";")
            push!(rows, (
                hour_start = Dates.format(hour_start, dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
                hour_end = Dates.format(hour_start + Hour(1), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
                status = complete ? "complete" : "incomplete",
                file_count = length(files),
                segment_count = length(readable_files),
                present_minutes = present_text,
                missing_minutes = missing_text,
                missing_count = length(missing_minutes),
                duplicate_count = duplicate_count,
                unreadable_count = length(unreadable_hour_files),
                unreadable_files = join(basename.(unreadable_hour_files), ";"),
                native_hourly_available = hour_start in native_hours,
                source_files = source_text,
            ))
            hour_start += Hour(1)
        end
    end

    return DataFrame(rows), complete_groups
end

"""
    pivot_to_wide(long_df)

Convert long format (time, station_id, fy4b_mm_h) to wide format
(time, station_1, station_2, ...)
"""
function pivot_to_wide(long_df)
    # Ensure station_id is string
    long_df.station_id = string.(long_df.station_id)

    # Pivot to wide format
    wide = unstack(long_df, :time, :station_id, :fy4b_mm_h)

    # Sort by time
    sort!(wide, :time)

    return wide
end

"""
    aggregate_fy4b_hourly(; kwargs...)

Main pipeline to aggregate FY4B data.
"""
function aggregate_fy4b_hourly(;
    data_dir::AbstractString=FY4B_DATA_DIR,
    station_meta_file::AbstractString=STATION_META_FILE,
    output_file::AbstractString=OUTPUT_FILE,
    qc_file::AbstractString=QC_FILE,
    years=2022:2025,
    start_month::Int=6,
    stop_month::Int=10,
)
    println("="^60)
    println("FY4B 15min -> Hourly Aggregation Pipeline")
    println("="^60)

    # Step 1: Load station metadata
    println("\n[1/5] Loading station metadata...")
    stations_df = load_stations(station_meta_file)
    all(name -> name in names(stations_df), ["station_id", "lon", "lat"]) ||
        error("Station metadata must contain station_id, lon, and lat columns")
    nrow(unique(stations_df, :station_id)) == nrow(stations_df) ||
        error("Station metadata contains duplicate station_id values")
    println("  Loaded $(nrow(stations_df)) stations")

    # Step 2: Find and group NC files
    println("\n[2/5] Finding 15min NC files...")
    selected_months = start_month:(stop_month - 1)
    file_groups = find_nc_files(data_dir; years, months=selected_months)
    println("  Validating selected NetCDF structures...")
    unreadable_files = find_unreadable_files(file_groups)
    for (file, message) in sort(collect(unreadable_files); by=first)
        @warn "Unreadable FY4B segment; treating as missing" file message
    end
    native_hours = find_native_hourly_times(data_dir)
    qc_df, complete_groups = build_hourly_qc(
        file_groups, native_hours;
        unreadable_files, years, start_month, stop_month,
    )
    total_hours = length(file_groups)
    total_files = sum(length.(values(file_groups)))
    complete_count = count(==("complete"), qc_df.status)
    incomplete_count = count(==("incomplete"), qc_df.status)
    missing_segment_count = sum(qc_df.missing_count)
    native_fallback_count = count((qc_df.status .== "incomplete") .& qc_df.native_hourly_available)
    println("  Found $total_hours represented hours with $total_files total 15min files")
    println("  Strict complete hours: $complete_count; incomplete hours: $incomplete_count")
    println("  Missing 15min segments: $missing_segment_count")
    println("  Incomplete hours with native hourly product (QC only): $native_fallback_count")
    @assert complete_count + incomplete_count == nrow(qc_df)
    @assert complete_count == length(complete_groups)
    write_csv_atomic(qc_file, qc_df)
    println("  QC saved: $qc_file")
    sat_lons = sort(unique(parse_satellite_lon(basename(file)) for files in values(file_groups) for file in files))
    sat_lon_text = join(string.(sat_lons), ", ")
    println("  Satellite subpoints in selected files: $sat_lon_text E")

    # Step 3: Aggregate to hourly
    println("\n[3/5] Aggregating precipitation to hourly...")
    long_df = aggregate_hourly(complete_groups, stations_df)
    println("  Aggregated $(nrow(long_df)) station-hour values")

    # Step 4: Convert to wide format
    println("\n[4/5] Converting to wide format...")
    wide_df = pivot_to_wide(long_df)
    println("  Wide format: $(size(wide_df, 1)) time rows × $(size(wide_df, 2)-1) station columns")

    # Step 5: Write output
    println("\n[5/5] Writing output...")
    # Format time column as ISO with Z
    wide_df.time = Dates.format.(wide_df.time, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z"
    write_csv_atomic(output_file, wide_df)
    println("  Saved: $output_file")

    # Summary statistics
    println("\n" * "="^60)
    println("Summary:")
    println("  Output dimensions: $(size(wide_df))")
    println("  Time range: $(wide_df.time[1]) to $(wide_df.time[end])")
    println("  Strict complete rows: $complete_count")
    println("  Station columns: $(size(wide_df, 2)-1)")
    println("="^60)
    return (; output_file, qc_file, complete_count, incomplete_count, unreadable_files)
end

end # module
