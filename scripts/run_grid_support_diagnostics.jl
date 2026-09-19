#!/usr/bin/env julia

# What the benchmark's satellite inputs lost when the gridded products were reduced to gauge
# points.
#
#   julia --project=. scripts/run_grid_support_diagnostics.jl              # D7 over 2022-06
#   julia --project=. scripts/run_grid_support_diagnostics.jl --window full   # D7 over the record
#   julia --project=. scripts/run_grid_support_diagnostics.jl --window none   # skip D7
#
# Reads the same tables the benchmark reads, through the same loader, so the station set and the
# hour grid are the published run's: 237 gauges on 13,471 common hours. Writes CSVs to
# output/grid_support_diagnostics/ and touches nothing else. Nothing under data/ is written.
#
# This measures; it decides nothing. `--window` only bounds D7, the nearest-versus-bilinear
# comparison, which is the one step that re-reads the FY4B NetCDF archive: one month is ~2,900
# files, the full record is ~160,000.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("TableIO", "HeavyRainEvents", "GridSupportDiagnostics")
using Main.TableIO: write_csv_atomic
using Main.HeavyRainEvents: align_to_reference
using Main.GridSupportDiagnostics
const GSD = Main.GridSupportDiagnostics
const FY4B = Main.FY4BPreprocessing

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const COVARIATES = joinpath(ROOT, "data", "processed", "covariates")
const DEM_DIR = joinpath(ROOT, "data", "processed", "dem")
const FY4B_DIR = joinpath(ROOT, "data", "FY4B")
const OUTDIR = joinpath(ROOT, "output", "grid_support_diagnostics")

# The benchmark's own window and hour count (scripts/run_interpolation_benchmark.jl:140-146). The
# expected count is kept so this fails loudly rather than quietly describing a different dataset
# from the one the published metrics came from.
const ANALYSIS_START = DateTime(2022, 1, 1, 9)
const ANALYSIS_END = DateTime(2025, 1, 1, 8)
const EXPECTED_COMMON_HOURS = 13471
const EXPECTED_STATIONS = 237

const PRODUCT_FILES = Dict(
    "FY4B" => "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
    "GPM" => "hubei_gpm_hourly_2022_2024_full_aligned.csv",
    "GSMaP" => "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
)
# The study box, from `StudyArea.STUDY_BOUNDS`; the footprint is probed at its corners and centre.
const BOUNDS = (west=109.4, east=111.6, south=31.2, north=33.4)
# The 1 km lat/lon step the domain rasters are aggregated to: 1/120 degree, the grid
# `data/dem_ShiYan_1km.tif` already sits on and a natural "fine" scale for this area.
const DOMAIN_STEP_DEG = 1 / 120
const D7_WINDOWS = Dict(
    "smoke" => (years=2022:2022, months=6:6),
    "full" => (years=2022:2024, months=1:12),
)
# Candidate lat/lon grids for GPM and GSMaP. The exports ask Earth Engine for `scale: 11132` m,
# which is 0.1 degree at the equator, but `reduceRegions` reprojects rather than reading the
# product's native cells, so which grid the stored values actually follow is an open question.
# `cell_grouping_agreement` answers it against the series the data can distinguish.
const COARSE_GRID_CANDIDATES = [
    ("0p1deg", 0.1), ("0p05deg", 0.05), ("0p2deg", 0.2), ("0p25deg", 0.25),
]

"""Parse `--window smoke|full|none`; absent means `smoke`, the bounded default."""
function parse_window(args)
    index = findfirst(startswith("--window"), args)
    index === nothing && return "smoke"
    text = args[index] == "--window" ? get(args, index + 1, "smoke") :
        split(args[index], "=", limit=2)[2]
    text in ("smoke", "full", "none") ||
        throw(ArgumentError("--window expects smoke, full or none, got $text"))
    return text
end

"""
The benchmark's own view of its inputs: the common station set, the common hour grid, and each
product on it.

Goes through `load_global_common_product_data` rather than reading the tables directly, so the
gauges and hours described below are the ones the published run scored, not a near-miss.
"""
function load_benchmark_inputs()
    cfg = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(product => joinpath(STUDY_DATA, file)
                       for (product, file) in PRODUCT_FILES),
        outdir=OUTDIR,
        analysis_start=ANALYSIS_START, analysis_end=ANALYSIS_END,
        expected_common_time_count=EXPECTED_COMMON_HOURS,
    )
    mkpath(OUTDIR)
    products, common_ids, product_data = load_global_common_product_data(cfg)
    length(common_ids) == EXPECTED_STATIONS || error(
        "Expected $EXPECTED_STATIONS common stations, got $(length(common_ids))")
    meta = load_station_meta(cfg.station_meta_path)
    lonlat = build_X_lonlat(meta, common_ids)
    times = product_data[first(products)].times
    return (; products, ids=common_ids, lonlat, times,
        Y_obs=product_data[first(products)].Y_obs,
        Y_sat=Dict(product => product_data[product].Y_sat for product in products))
end

"""
The sub-satellite longitudes present in the FY4B archive, ascending.

Read from the filenames rather than assumed: `extract_precipitation` takes each file's own
subpoint, so every distinct one is a separate projection the shipped column was built through.
"""
function archive_satellite_longitudes(data_dir::AbstractString)
    found = Set{Float64}()
    for (_, _, files) in walkdir(data_dir), file in files
        endswith(uppercase(file), ".NC") || continue
        push!(found, FY4B.parse_satellite_lon(file))
    end
    isempty(found) && error("No FY4B NetCDF files found under $data_dir")
    return sort(collect(found))
end

"""FY4B sampled both ways over `window`, against the shipped column on the same hours."""
function extraction_sensitivity(window::AbstractString, ids, lonlat, sat_times, sat_Y)
    years, months = D7_WINDOWS[window]
    println("[D7] re-sampling FY4B over $window ($(first(years))-$(last(years)), months $months)")
    file_groups = FY4B.find_nc_files(FY4B_DIR; years, months)
    isempty(file_groups) && error("No FY4B files found under $FY4B_DIR for $window")
    unreadable_files = FY4B.find_unreadable_files(file_groups)
    native_hours = FY4B.find_native_hourly_times(FY4B_DIR)
    _, complete_groups = FY4B.build_hourly_qc(
        file_groups, native_hours; unreadable_files, years, months)
    println("  $(length(complete_groups)) strict-complete hours of $(length(file_groups))")
    times, Y_nearest, Y_bilinear = GSD.fy4b_hourly_modes(
        complete_groups, lonlat[:, 1], lonlat[:, 2])
    reference = align_to_reference(times, ids, sat_times, ids, sat_Y)
    return GSD.extraction_sensitivity_table(ids, Y_nearest, Y_bilinear; Y_reference=reference)
end

"""
One flat table of the numbers this exists to produce.

Only the groupings that survived D3 are promoted here: `identical_<product>` for how much each
product merges, `fy4b_pixel_effective` for the geometry, and the best-scoring candidate grid for
each product so the shortfall is on the record rather than buried in a CSV.
"""
function summary_table(
    multiplicity, agreement, identical, disagreement, footprint, terrain, sensitivity,
)
    rows = NamedTuple[]
    add(question, value, unit) = push!(rows, (; question, value=Float64(value), unit))
    for product in unique(agreement.product)
        groups = identical[identical.product .== product, :]
        merged = sum(groups.n_stations; init=0)
        add("$product: gauges it cannot tell apart from a neighbour", merged, "stations")
        add("$product: fraction of gauges it cannot tell apart", merged / EXPECTED_STATIONS,
            "fraction")
        add("$product: largest set of gauges given one value",
            isempty(groups.n_stations) ? 1 : maximum(groups.n_stations), "stations")
        # A grid that accounts for every cell group wins outright; among those, the finest one.
        # Otherwise report the closest miss, which is the honest statement about that product.
        candidates = agreement[agreement.product .== product, :]
        explaining = candidates[candidates.grid_explains_series, :]
        row = nrow(explaining) > 0 ? explaining[argmax(explaining.n_cell_groups), :] :
            candidates[argmax(candidates.cell_groups_confirmed), :]
        add("$product: best candidate grid ($(row.grid)) explains", row.cell_groups_confirmed,
            "of $(row.n_cell_groups) cell groups")
        add("$product: that grid accounts for the whole grouping",
            row.grid_explains_series ? 1 : 0,
            "1 = yes")
    end
    for row in eachrow(multiplicity)
        row.grid == "fy4b_pixel_effective" || continue
        add("gauges sharing an FY4B pixel in every projection", row.n_stations_sharing, "stations")
        add("most gauges in one FY4B pixel", row.max_stations_per_cell, "stations")
    end
    for row in eachrow(disagreement)
        row.separation_bin_km == "0.0_5.0" || continue
        add("gauge-vs-gauge RMSE under 5 km, $(row.grid) " *
            (row.shares_cell ? "same cell" : "different cells"), row.median_RMSE, "mm/h")
    end
    for row in eachrow(footprint)
        if row.scope == "gauges"
            add("FY4B footprint at the gauges from $(row.sat_lon)E", row.footprint_km2, "km2")
            add("FY4B footprint vs its 4 km nadir spec from $(row.sat_lon)E", row.nadir_ratio, "x")
        elseif row.name == "fy4b_pixels_in_study_box"
            add("FY4B pixels covering the study box from $(row.sat_lon)E",
                row.footprint_km2, "pixels")
            add("fraction of those pixels holding a gauge from $(row.sat_lon)E",
                row.nadir_ratio, "fraction")
        end
    end
    for variable in unique(terrain.variable)
        rows_v = terrain[terrain.variable .== variable, :]
        value(statistic) = only(rows_v[rows_v.statistic .== statistic, :])
        # The range fraction understates the bias badly - the gauges sit inside the terrain's
        # range but low in it - so the central tendency and the upper tail go in too.
        add("$variable: gauge mean", value("mean").gauges, "units")
        add("$variable: domain mean", value("mean").domain, "units")
        add("$variable: gauge median vs domain median", value("p50").gauges - value("p50").domain,
            "units")
        add("$variable: highest gauge vs highest domain cell",
            value("p100").gauges - value("p100").domain, "units")
        add("$variable: range of the domain no gauge occupies",
            value("domain_fraction_outside_gauge_range").domain, "fraction")
    end
    if sensitivity !== nothing
        row = first(eachrow(sensitivity))
        add("nearest-vs-bilinear mean absolute difference", row.mean_abs_delta, "mm/h")
        add("nearest-vs-bilinear RMSE", row.rmse_delta, "mm/h")
        add("nearest-vs-bilinear RMSE over wet hours", row.rmse_delta_wet, "mm/h")
        add("nearest-vs-bilinear largest single-hour difference", row.max_abs_delta, "mm/h")
        add("wet hours where the two disagree about rain at all",
            row.fraction_wet_flip_of_wet, "fraction")
        add("cells where this run failed to reproduce the shipped column",
            row.reference_mismatch, "cells")
    end
    return DataFrame(rows)
end

function main(args=ARGS)
    window = parse_window(args)
    println("Grid support diagnostics: D7 window=$window, output=$OUTDIR")
    input = load_benchmark_inputs()
    println("  $(length(input.ids)) gauges on $(length(input.times)) common hours")

    # D1: the cell every gauge actually reads, once per projection the archive was built through.
    subpoints = archive_satellite_longitudes(FY4B_DIR)
    println("  FY4B sub-satellite longitudes in the archive: $(join(subpoints, ", "))E")
    by_subpoint = Dict(sat_lon => GSD.pixel_assignment_table(
        input.ids, input.lonlat[:, 1], input.lonlat[:, 2]; sat_lon) for sat_lon in subpoints)
    assignment = reduce(vcat, [by_subpoint[sat_lon] for sat_lon in subpoints])
    write_csv_atomic(joinpath(OUTDIR, "pixel_assignment.csv"), assignment)

    # D2: how many gauges a cell cannot tell apart. A gauge has one FY4B pixel per projection, and
    # two gauges are handed the same value only where they share a pixel in every one of them, so
    # `fy4b_pixel_effective` is the grouping that governs. The coarse grids are candidates that D3
    # then tests, because what Earth Engine sampled GPM and GSMaP on is not recorded anywhere.
    fy4b_label(sat_lon) = "fy4b_pixel_$(replace(string(sat_lon), "." => "p"))E"
    grids = Dict{String,Vector{String}}(
        fy4b_label(sat_lon) => by_subpoint[sat_lon].fy4b_pixel for sat_lon in subpoints)
    grids["fy4b_pixel_effective"] = GSD.effective_cell_keys(
        [by_subpoint[sat_lon].fy4b_pixel for sat_lon in subpoints])
    for (label, step) in COARSE_GRID_CANDIDATES
        grids["coarse_$label"] = [string(GSD.latlon_cell_index(
            input.lonlat[i, 1], input.lonlat[i, 2]; step_deg=step))
            for i in eachindex(input.ids)]
    end
    membership = GSD.cell_membership_table(input.ids, grids)
    multiplicity = GSD.cell_multiplicity_summary(membership)
    write_csv_atomic(joinpath(OUTDIR, "cell_membership.csv"), membership)
    write_csv_atomic(joinpath(OUTDIR, "cell_multiplicity.csv"), multiplicity)
    println("  cell multiplicity:")
    println(multiplicity)

    # D3: the same question from the shipped tables alone, and which candidate grid it endorses.
    # The series relation is computed once per product and reused across candidates.
    series = Dict(product => first(GSD.identical_series_groups(input.Y_sat[product]))
                  for product in input.products)
    identical = reduce(vcat, [GSD.identical_series_table(input.ids, input.Y_sat[product], product)
                              for product in input.products])
    candidates(product) = product == "FY4B" ?
        [key for key in keys(grids) if startswith(key, "fy4b_pixel")] :
        [key for key in keys(grids) if startswith(key, "coarse_")]
    agreement = reduce(vcat, [GSD.cell_grouping_agreement(
        input.ids, grids[grid], input.Y_sat[product];
        grid, product, series_groups=series[product])
        for product in input.products for grid in sort(candidates(product))])
    write_csv_atomic(joinpath(OUTDIR, "identical_series.csv"), identical)
    write_csv_atomic(joinpath(OUTDIR, "cell_grouping_agreement.csv"), agreement)
    for product in input.products
        rows = agreement[agreement.product .== product, :]
        any(rows.grid_explains_series) || @warn(
            "No candidate grid accounts for every group of gauges this product cannot separate; " *
            "read identical_series.csv as the answer and treat the cell counts as candidates",
            product, best=maximum(rows.cell_groups_confirmed))
    end
    println("  grid/series agreement:")
    println(agreement)

    # D4: what two gauges inside one cell disagree by, with non-sharing pairs as the control.
    # Flagged by what each product actually merges, not by a candidate grid: for FY4B that is the
    # effective pixel, and for the Earth Engine products only the series grouping is trustworthy.
    share_flags = Dict{String,Vector{String}}(
        "fy4b_pixel" => grids["fy4b_pixel_effective"])
    for product in input.products
        share_flags["identical_$product"] = GSD.series_group_keys(input.ids, series[product])
    end
    pairs = GSD.gauge_pair_table(input.ids, input.lonlat, input.Y_obs, share_flags)
    disagreement = GSD.pair_disagreement_summary(pairs)
    write_csv_atomic(joinpath(OUTDIR, "gauge_pair_disagreement.csv"), pairs)
    write_csv_atomic(joinpath(OUTDIR, "subcell_disagreement_summary.csv"), disagreement)
    println("  sub-cell gauge disagreement:")
    println(disagreement)

    # D5: FY4B's real footprint here, and how much of the domain holds a gauge at all. Once per
    # projection: from 105E the study area is nearly under the satellite, from 133E it is 23
    # degrees off, and the shipped column is built from both.
    probes(sat_lon) = [
        ("study_box_center", (BOUNDS.west + BOUNDS.east) / 2, (BOUNDS.south + BOUNDS.north) / 2),
        ("study_box_nw", BOUNDS.west, BOUNDS.north), ("study_box_ne", BOUNDS.east, BOUNDS.north),
        ("study_box_sw", BOUNDS.west, BOUNDS.south), ("study_box_se", BOUNDS.east, BOUNDS.south),
        ("sub_satellite_point", sat_lon, 0.0)]
    footprint = reduce(vcat, [insertcols!(
        GSD.footprint_table(by_subpoint[sat_lon], probes(sat_lon), BOUNDS; sat_lon),
        1, :sat_lon => sat_lon) for sat_lon in subpoints])
    write_csv_atomic(joinpath(OUTDIR, "fy4b_footprint.csv"), footprint)
    println("  FY4B footprint:")
    println(footprint)

    # D6: how biased the gauge network is as a sample of the terrain it validates over. Gauge and
    # domain values come from the same Copernicus rasters, so the comparison is of one dataset
    # with itself at two supports rather than of two sources.
    station_terrain = CSV.read(joinpath(COVARIATES, "station_terrain.csv"), DataFrame;
        types=Dict(:station_id => String))
    terrain = reduce(vcat, [GSD.terrain_representativeness_table(
        station_terrain[!, column],
        GSD.domain_raster_values(joinpath(DEM_DIR, raster), BOUNDS, DOMAIN_STEP_DEG),
        String(column))
        for (column, raster) in ((:elevation_m, "copernicus_glo30_utm49n_30m.tif"),
            (:slope_deg, "copernicus_glo30_slope_deg.tif"))])
    write_csv_atomic(joinpath(OUTDIR, "gauge_terrain_representativeness.csv"), terrain)

    # D7: what nearest-pixel extraction discarded, against bilinear on the same files and hours.
    sensitivity = window == "none" ? nothing :
        extraction_sensitivity(window, input.ids, input.lonlat, input.times, input.Y_sat["FY4B"])
    sensitivity_path = joinpath(OUTDIR, "extraction_sensitivity.csv")
    # A skipped D7 must not leave the previous run's table beside a summary that no longer
    # mentions it: the directory has to describe one run, not two.
    if sensitivity === nothing && isfile(sensitivity_path)
        rm(sensitivity_path)
        println("  removed a previous run's extraction_sensitivity.csv (D7 was skipped)")
    end
    if sensitivity !== nothing
        write_csv_atomic(sensitivity_path, sensitivity)
        mismatch = first(sensitivity).reference_mismatch
        mismatch == 0 || @warn(
            "Nearest sampling did not reproduce the shipped FY4B column; the bilinear " *
            "comparison is not controlled", mismatch)
        println("  extraction sensitivity (all stations):")
        println(sensitivity[1:1, :])
    end

    # What produced these tables, so a reader a month from now does not have to guess which
    # window D7 ran over or whether the inputs were still the published run's.
    write_csv_atomic(joinpath(OUTDIR, "run_scope.csv"), DataFrame([
        (item="generated_utc", value=string(now(UTC))),
        (item="script", value="scripts/run_grid_support_diagnostics.jl"),
        (item="d7_window", value=window),
        (item="analysis_start", value=string(ANALYSIS_START)),
        (item="analysis_end", value=string(ANALYSIS_END)),
        (item="common_hours", value=string(length(input.times))),
        (item="common_stations", value=string(length(input.ids))),
        (item="products", value=join(input.products, ";")),
        (item="fy4b_subpoints_E", value=join(subpoints, ";")),
        (item="domain_grid_step_deg", value=string(DOMAIN_STEP_DEG)),
        (item="measures", value="satellite inputs only; no benchmark output is read or written"),
    ]))

    summary = summary_table(multiplicity, agreement, identical, disagreement, footprint,
        terrain, sensitivity)
    write_csv_atomic(joinpath(OUTDIR, "summary.csv"), summary)
    println("\nSummary:")
    show(summary; allrows=true, allcols=true, truncate=0)
    println("\n\nWritten to $OUTDIR")
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
