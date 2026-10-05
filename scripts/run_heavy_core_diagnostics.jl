#!/usr/bin/env julia

# Why do GPM and GSMaP underestimate the heavy-rain cores of the selected events?
#
#   julia --project=. scripts/run_heavy_core_diagnostics.jl
#
# Reads the six events chosen by `run_heavy_rain_event_evaluation.jl --gpm-gsmap`, takes the core
# gauges of each (>= 50 mm on a localized day, >= 100 mm on a widespread one), and splits each core
# gauge's deficit into representativeness, displacement and magnitude (`HeavyCoreDiagnostics`).
# Three products: IMERG Final and GSMaP as used in the paper, and IMERG Uncal (no monthly gauge
# calibration), which is also on disk as the native 0.1 deg half-hourly grid - the only product
# here whose neighbourhood can be searched off the gauge pixels. Writes CSVs to
# output/heavy_core_diagnostics/.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics, NCDatasets

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("TableIO", "HeavyRainEvents", "NoRainEvaluation", "HeavyCoreDiagnostics")
using Main.TableIO: write_csv_atomic
using Main.HeavyRainEvents: met_day, align_to_reference, phase_spatial_metrics, small_scale_variance_share
using Main.NoRainEvaluation: long_to_wide_hourly
using Main.HeavyCoreDiagnostics

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const EVENTS_FILE = joinpath(ROOT, "output", "heavy_rain_events_gpm_gsmap", "selected_events.csv")
const UNCAL_STATION_DIR = joinpath(ROOT, "output", "gpm_imerg_uncal_stations_2022_2024")
const UNCAL_GRID_DIR = joinpath(ROOT, "data", "raw", "gpm_imerg_v07_uncal_hubei_2022_2024")
const OUTDIR = joinpath(ROOT, "output", "heavy_core_diagnostics")
const ANALYSIS_START = DateTime(2022, 1, 1, 9)
const ANALYSIS_END = DateTime(2025, 1, 1, 8)
const PRODUCT_FILES = [
    "GPM" => "hubei_gpm_hourly_2022_2024_full_aligned.csv",
    "GSMaP" => "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
]
const PRODUCTS = ["GPM", "GSMaP", "GPM_uncal"]
const CORE_MM = Dict("localized" => 50.0, "widespread" => 100.0)
const AREAL_RADII_KM = (5.0, 10.0, 15.0)
const NEARBY_RADII_KM = (10.0, 25.0, 50.0)
# The attribution uses a 10 km disc as the areal truth (~3 pixels, so it errs towards crediting
# pixel smoothing) and a 25 km search for displacement.
const AREAL_KM = 10.0
const NEARBY_KM = 25.0
const PAD_H = 3
# Three 8 h phases of the 08-08 day: 09-16, 17-00, 01-08 BJT (hour-ending labels).
const N_PHASES = 3

function load_inputs()
    obs_times, obs_ids, Y_obs = read_hourly_wide(joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"))
    obs_times, Y_obs = subset_time_window(obs_times, Y_obs, ANALYSIS_START, ANALYSIS_END)
    products = Dict{String,Matrix{Float64}}()
    for (product, file) in PRODUCT_FILES
        times, ids, Y = read_hourly_wide(joinpath(STUDY_DATA, file))
        times, Y = subset_time_window(times, Y, ANALYSIS_START, ANALYSIS_END)
        products[product] = align_to_reference(obs_times, obs_ids, times, ids, Y)
    end
    lonlat = build_X_lonlat(load_station_meta(joinpath(STUDY_DATA, "station_meta.csv")), obs_ids)
    return obs_times, obs_ids, Y_obs, products, lonlat
end

"""IMERG Uncal at the gauges for the months the events touch (UTC hour start + 9 h = BJT hour ending)."""
function load_uncalibrated_gpm(obs_times, obs_ids, days)
    months = unique(Dates.format(d + Day(k), "yyyymm") for d in days for k in -1:1)
    time, station, value = DateTime[], String[], Float64[]
    for month in sort(months)
        path = joinpath(UNCAL_STATION_DIR, "gpm_hubei_hourly_long_$(month).csv")
        isfile(path) || continue
        table = CSV.read(path, DataFrame; types=Dict(:station_id => String, :time => String))
        append!(time, DateTime.(first.(table.time, 19)))
        append!(station, table.station_id)
        append!(value, [ismissing(v) ? NaN : Float64(v) for v in table.gpm_mm_h])
    end
    return long_to_wide_hourly(time, station, value, obs_times, obs_ids; offset_hours=9).Y
end

"""
The native IMERG Uncal 08-08 BJT total for `day` (= UTC day `day`), `(lon, lat, total[lat, lon])`.

Times come from the file names: the files' `time` variable declares a Julian calendar, which
decodes 13 days late.
"""
function uncal_grid_total(day::Date)
    dir = joinpath(UNCAL_GRID_DIR, string(year(day)), lpad(dayofyear(day), 3, '0'))
    files = sort(filter(f -> occursin(Dates.format(day, "yyyymmdd") * "-S", f) && endswith(f, ".nc4"), readdir(dir)))
    length(files) == 48 || error("expected 48 half-hourly IMERG files for $day in $dir, found $(length(files))")
    lon = lat = total = nothing
    for file in files
        NCDataset(joinpath(dir, file)) do ds
            dimnames(ds["precipitationUncal"]) == ("lat", "lon", "time") || error("unexpected dimensions in $file")
            rate = Float64.(coalesce.(ds["precipitationUncal"][:, :, 1], NaN))
            if total === nothing
                lon, lat, total = Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), zeros(size(rate))
            end
            total .+= 0.5 .* rate  # mm/h over a half hour
        end
    end
    return lon, lat, total
end

nearest_cell(lon, lat, point) = (argmin(abs.(lat .- point[2])), argmin(abs.(lon .- point[1])))

pick(rows, radius) = only(filter(r -> r.radius_km == radius, rows))

function main()
    println("Reading study-area hourly tables ...")
    obs_times, obs_ids, Y_obs, products, lonlat = load_inputs()
    events = CSV.read(EVENTS_FILE, DataFrame)
    events = events[events.selected, :]
    products["GPM_uncal"] = load_uncalibrated_gpm(obs_times, obs_ids, events.day)

    point_of(i) = (lonlat[i, 1], lonlat[i, 2])
    core_rows, areal_rows, nearby_rows, hourly_rows, attribution_rows, check_rows = (NamedTuple[] for _ in 1:6)
    phase_tables = DataFrame[]
    for event in eachrow(events)
        day = event.day
        hour_idx = findall(t -> met_day(t) == day, obs_times)
        padded_idx = (first(hour_idx) - PAD_H):(last(hour_idx) + PAD_H)
        total(Y) = [all(isfinite, @view Y[i, hour_idx]) ? sum(@view Y[i, hour_idx]) : NaN for i in axes(Y, 1)]
        obs = total(Y_obs)
        sat = Dict(p => total(products[p]) for p in PRODUCTS)

        # Phase decomposition: does a poor day-total pattern hide well-captured phases?
        phases = phase_spatial_metrics(Y_obs, Dict(p => products[p] for p in ("GPM", "GSMaP")), lonlat, hour_idx; n_phases=N_PHASES)
        insertcols!(phases, 1, :event_day => day, :event_type => event.event_type,
            :obs_small_scale_share => small_scale_variance_share(lonlat, obs))
        phases.first_hour = obs_times[phases.first_hour]
        push!(phase_tables, phases)
        core = sort(findall(i -> isfinite(obs[i]) && obs[i] >= CORE_MM[event.event_type]
            && all(p -> isfinite(sat[p][i]), PRODUCTS), eachindex(obs)); by=i -> -obs[i])
        println("$(day) $(event.event_type): $(length(core)) core gauges (>= $(CORE_MM[event.event_type]) mm)")

        grid_lon, grid_lat, grid_total = uncal_grid_total(day)
        grid_points = [grid_lon[c] for r in eachindex(grid_lat), c in eachindex(grid_lon)]
        grid_lonlat = hcat(vec(grid_points), vec([grid_lat[r] for r in eachindex(grid_lat), c in eachindex(grid_lon)]))
        grid_values = vec(grid_total)

        # Sanity check: the grid at each gauge's nearest cell against the gauge-sampled Uncal series.
        # A gauge exactly on a cell edge (lon or lat a multiple of 0.1 deg) is a tie, which the two
        # samplings may break differently; only interior gauges have to match.
        for i in eachindex(obs_ids)
            isfinite(sat["GPM_uncal"][i]) || continue
            cell = nearest_cell(grid_lon, grid_lat, point_of(i))
            push!(check_rows, (; event_day=day, station_id=obs_ids[i], station_series_mm=sat["GPM_uncal"][i],
                on_cell_edge=any(x -> abs(10x - round(10x)) < 1e-6, point_of(i)),
                grid_cell_mm=grid_total[cell...], diff_mm=grid_total[cell...] - sat["GPM_uncal"][i]))
        end

        for i in core
            point = (lonlat[i, 1], lonlat[i, 2])
            consistency = neighbour_consistency(lonlat, obs, i)
            peak_hour = obs_times[hour_idx[argmax(Y_obs[i, hour_idx])]]
            push!(core_rows, merge((; event_day=day, event_type=event.event_type, station_id=obs_ids[i],
                lon=point[1], lat=point[2], obs_mm=obs[i], obs_peak_mm_h=maximum(Y_obs[i, hour_idx]),
                obs_peak_hour_bjt=peak_hour), consistency,
                NamedTuple{Tuple(Symbol.(PRODUCTS, "_mm"))}(Tuple(sat[p][i] for p in PRODUCTS))))
            areal = upscaled_gauge_truth(lonlat, obs, i; radii_km=AREAL_RADII_KM)
            for a in areal
                push!(areal_rows, merge((; event_day=day, station_id=obs_ids[i], obs_mm=obs[i]), a))
            end
            for p in PRODUCTS
                sources = [("gauge_pixels", lonlat, sat[p])]
                p == "GPM_uncal" && push!(sources, ("native_grid", grid_lonlat, grid_values))
                nearby_by_source = Dict{String,Float64}()
                for (source, points, values) in sources
                    rows = neighbourhood_max(points, values, point; radii_km=NEARBY_RADII_KM)
                    for r in rows
                        push!(nearby_rows, merge((; event_day=day, station_id=obs_ids[i], product=p, source), r))
                    end
                    nearby_by_source[source] = pick(rows, NEARBY_KM).max
                end
                # The native grid is the honest search where it exists; gauge pixels are all the
                # aligned products offer.
                nearby = get(nearby_by_source, "native_grid", nearby_by_source["gauge_pixels"])
                parts = deficit_attribution(obs[i], sat[p][i], pick(areal, AREAL_KM).mean, nearby)
                push!(attribution_rows, merge((; event_day=day, event_type=event.event_type, station_id=obs_ids[i],
                    product=p, obs_mm=obs[i], sat_mm=sat[p][i], areal_mm=pick(areal, AREAL_KM).mean,
                    nearby_mm=nearby, nearby_source=haskey(nearby_by_source, "native_grid") ? "native_grid" : "gauge_pixels"),
                    parts))
                hourly = core_hourly_diagnostics(Y_obs[i, padded_idx], products[p][i, padded_idx]; pad=PAD_H)
                push!(hourly_rows, merge((; event_day=day, event_type=event.event_type, station_id=obs_ids[i], product=p), hourly))
            end
        end
    end

    core = DataFrame(core_rows)
    attribution = DataFrame(attribution_rows)
    hourly = DataFrame(hourly_rows)
    checks = DataFrame(check_rows)
    summary = combine(groupby(innerjoin(attribution, select(hourly, :event_day, :station_id, :product,
            :peak_ratio, :top_share_obs, :top_deficit_share, :peak_lag_h, :best_lag_h, :sat_total_padded),
            on=[:event_day, :station_id, :product]), [:event_day, :event_type, :product]; sort=true),
        nrow => :n_core,
        :obs_mm => mean => :obs_mm, :sat_mm => mean => :sat_mm, :areal_mm => mean => :areal_mm,
        :nearby_mm => mean => :nearby_mm,
        [:representativeness, :deficit] => ((x, d) -> sum(x) / sum(d)) => :representativeness_share,
        [:displacement, :deficit] => ((x, d) -> sum(x) / sum(d)) => :displacement_share,
        [:magnitude, :deficit] => ((x, d) -> sum(x) / sum(d)) => :magnitude_share,
        :peak_ratio => median => :median_peak_ratio,
        :top_share_obs => median => :median_top3_share_of_obs,
        :top_deficit_share => median => :median_top3_share_of_deficit,
        :peak_lag_h => median => :median_peak_lag_h,
        :best_lag_h => median => :median_best_lag_h,
        :sat_total_padded => mean => :sat_padded_mm,
    )
    uncal = select(filter(:product => ==("GPM_uncal"), summary), :event_day, :sat_mm => :uncal_mm)
    summary = leftjoin(summary, uncal; on=:event_day)
    summary.final_over_uncal = ifelse.(summary.product .== "GPM", summary.sat_mm ./ summary.uncal_mm, NaN)
    select!(summary, Not(:uncal_mm))
    sort!(summary, [:event_day, :product])

    write_csv_atomic(joinpath(OUTDIR, "core_gauges.csv"), core)
    write_csv_atomic(joinpath(OUTDIR, "upscaling.csv"), DataFrame(areal_rows))
    write_csv_atomic(joinpath(OUTDIR, "neighbourhood_max.csv"), DataFrame(nearby_rows))
    write_csv_atomic(joinpath(OUTDIR, "core_hourly.csv"), hourly)
    write_csv_atomic(joinpath(OUTDIR, "deficit_attribution.csv"), attribution)
    write_csv_atomic(joinpath(OUTDIR, "grid_station_check.csv"), checks)
    phase_table = vcat(phase_tables...)
    write_csv_atomic(joinpath(OUTDIR, "phase_decomposition.csv"), phase_table)
    write_csv_atomic(joinpath(OUTDIR, "summary.csv"), summary)

    interior = checks[.!checks.on_cell_edge, :]
    println("\nNative grid vs gauge-sampled IMERG Uncal, daily totals: max |diff| = ",
        round(maximum(abs, interior.diff_mm); digits=4), " mm over $(nrow(interior)) interior gauge-days ",
        "($(count(checks.on_cell_edge)) gauge-days on a cell edge not compared)")
    println("\nCore-deficit attribution (areal truth = $(AREAL_KM) km disc, displacement search = $(NEARBY_KM) km):")
    show(summary; allrows=true, allcols=true)
    println("

Spatial r per 8 h phase and over the day (GPM, GSMaP):")
    show(select(phase_table, :event_day, :phase, :product, :obs_share, :r, :Bias, :obs_r_lat, :est_r_lat,
        :obs_r_lon, :est_r_lon, :obs_small_scale_share); allrows=true, allcols=true)
    println("\n\nWrote $(OUTDIR)")
    return OUTDIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
