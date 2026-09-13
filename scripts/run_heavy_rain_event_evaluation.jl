#!/usr/bin/env julia

# Heavy-rainfall event evaluation: do FY4B, GPM and GSMaP reproduce the spatial pattern of
# gauge-observed rain on representative heavy-rain days?
#
#   julia --project=. scripts/run_heavy_rain_event_evaluation.jl
#
# Screens 2022-2024 for days with more than three gauges above 50 mm (08-08 BJT), selects four
# widespread and two localized events, and scores each product against the gauges over a common
# hour window. Writes CSVs to output/heavy_rain_events/; the figures are drawn from those by
# `py -3.13 scripts/plot_heavy_rain_events.py`.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("TableIO", "StudyArea", "HeavyRainEvents")
using Main.TableIO: write_csv_atomic
using Main.StudyArea: STUDY_BOUNDS
using Main.HeavyRainEvents

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const OUTDIR = joinpath(ROOT, "output", "heavy_rain_events")
# The full benchmark window: the first and last hour-ending labels of the 08-08 BJT days.
const ANALYSIS_START = DateTime(2022, 1, 1, 9)
const ANALYSIS_END = DateTime(2025, 1, 1, 8)
const PRODUCT_FILES = [
    "FY4B" => "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
    "GPM" => "hubei_gpm_hourly_2022_2024_full_aligned.csv",
    "GSMaP" => "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
]
const PRODUCTS = first.(PRODUCT_FILES)
# Products with a complete 24 h on the selected days, scored a second time over the full day to
# check that the common window FY4B's gaps impose does not drive the conclusions.
const FULL_DAY_PRODUCTS = ["GPM", "GSMaP"]
const IDW_STEP_DEG = 0.025
const IDW_POWER = 2.0
const ID_COLUMNS = [:event_day, :event_type, :product, :window, :hours]

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

"""Mean and median of every metric across events, plus the metrics of all event-station pairs pooled."""
function summarize_metrics(metrics::DataFrame, stations::DataFrame, lonlat_by_id::Dict{String,Tuple{Float64,Float64}})
    metric_columns = setdiff(propertynames(metrics), ID_COLUMNS)
    rows = DataFrame[]
    for group in groupby(metrics, [:product, :window]; sort=true)
        for (statistic, reduce) in (("mean_over_events", mean), ("median_over_events", median))
            row = DataFrame(product=group.product[1], window=group.window[1], statistic=statistic, n_events=nrow(group))
            for column in metric_columns
                values = filter(isfinite, Float64.(group[!, column]))
                row[!, column] = [isempty(values) ? NaN : reduce(values)]
            end
            push!(rows, row)
        end
        pairs = stations[(stations.product .== group.product[1]) .& (stations.window .== group.window[1]), :]
        lonlat = [lonlat_by_id[id][k] for id in pairs.station_id, k in 1:2]
        pooled = spatial_metrics(pairs.obs_mm, pairs.sat_mm, lonlat)
        row = DataFrame(product=group.product[1], window=group.window[1], statistic="pooled_pairs", n_events=nrow(group))
        for column in metric_columns
            row[!, column] = [Float64(pooled[column])]
        end
        # Pooling repeats each gauge once per event, so a single rain centre or peak is meaningless.
        row.centroid_shift_km .= NaN
        row.peak_shift_km .= NaN
        push!(rows, row)
    end
    return vcat(rows...)
end

function main()
    println("Reading study-area hourly tables ...")
    obs_times, obs_ids, Y_obs, products, lonlat = load_inputs()
    println("  $(length(obs_ids)) gauges x $(length(obs_times)) hours, products: $(join(PRODUCTS, ", "))")

    days, obs_daily = daily_totals(obs_times, Y_obs)
    candidates = tag_rain_processes!(screen_heavy_rain_days(days, obs_daily, obs_ids))
    coverage = Dict(day => event_hours(obs_times, Y_obs, products, day) for day in candidates.day)
    for product in PRODUCTS
        candidates[!, Symbol(lowercase(product), "_hours")] = [coverage[day].product_hours[product] for day in candidates.day]
    end
    candidates.common_hours = [coverage[day].common_hours for day in candidates.day]
    candidates.excluded_rain_share = [coverage[day].excluded_rain_share for day in candidates.day]
    candidates.peak_excluded_rain_share = [coverage[day].peak_excluded_rain_share for day in candidates.day]
    candidates = select_representative_events(candidates)
    selected = candidates[candidates.selected, :]
    println("  $(nrow(candidates)) heavy-rain days, $(count(candidates.eligible)) eligible, selected:")
    for event in eachrow(selected)
        println("    $(event.day)  $(rpad(event.event_type, 10))  gauges>50mm=$(event.n_heavy)  " *
            "max=$(event.max_mm) mm  common hours=$(event.common_hours)  " *
            "excluded rain=$(round(100 * event.excluded_rain_share; digits=1))% " *
            "(peak gauge $(round(100 * event.peak_excluded_rain_share; digits=1))%)")
    end

    day_column = Dict(day => k for (k, day) in enumerate(days))
    metric_rows = NamedTuple[]
    station_rows = NamedTuple[]
    surfaces = DataFrame[]
    bounds = (STUDY_BOUNDS.west, STUDY_BOUNDS.east, STUDY_BOUNDS.south, STUDY_BOUNDS.north)
    for event in eachrow(selected)
        hours = coverage[event.day]
        windows = [
            ("common_hours", hours.hour_idx, PRODUCTS),
            ("full_day", hours.day_hour_idx, FULL_DAY_PRODUCTS),
        ]
        for (window, hour_idx, window_products) in windows
            totals = event_station_totals(Y_obs, Dict(p => products[p] for p in window_products), hour_idx)
            for product in window_products
                scores = spatial_metrics(totals.obs, totals.sat[product], lonlat)
                push!(metric_rows, merge(
                    (; event_day=event.day, event_type=event.event_type, product, window, hours=length(hour_idx)),
                    scores,
                ))
                for i in findall(totals.keep)
                    push!(station_rows, (;
                        event_day=event.day, event_type=event.event_type, window, product,
                        station_id=obs_ids[i], lon=lonlat[i, 1], lat=lonlat[i, 2],
                        obs_24h_mm=obs_daily[i, day_column[event.day]],
                        obs_mm=totals.obs[i], sat_mm=totals.sat[product][i],
                        diff_mm=totals.sat[product][i] - totals.obs[i],
                    ))
                end
            end
            window == "common_hours" || continue
            for (source, values) in vcat(["Gauge" => totals.obs], [p => totals.sat[p] for p in PRODUCTS])
                surface = idw_surface(lonlat[totals.keep, :], values[totals.keep]; bounds, step_deg=IDW_STEP_DEG, power=IDW_POWER)
                insertcols!(surface, 1, :event_day => event.day, :source => source)
                push!(surfaces, surface)
            end
        end
    end
    metrics = DataFrame(metric_rows)
    stations = DataFrame(station_rows)
    lonlat_by_id = Dict(id => (lonlat[i, 1], lonlat[i, 2]) for (i, id) in enumerate(obs_ids))
    summary = summarize_metrics(metrics, stations, lonlat_by_id)

    write_csv_atomic(joinpath(OUTDIR, "candidate_days.csv"), candidates)
    write_csv_atomic(joinpath(OUTDIR, "selected_events.csv"), selected)
    write_csv_atomic(joinpath(OUTDIR, "station_event_totals.csv"), stations)
    write_csv_atomic(joinpath(OUTDIR, "event_metrics.csv"), metrics)
    write_csv_atomic(joinpath(OUTDIR, "event_metrics_summary.csv"), summary)
    write_csv_atomic(joinpath(OUTDIR, "idw_surfaces.csv"), vcat(surfaces...))
    write_csv_atomic(joinpath(OUTDIR, "idw_settings.csv"), DataFrame(
        method="inverse distance weighting", power=IDW_POWER, neighbors="all gauges with a finite value",
        step_deg=IDW_STEP_DEG, west=bounds[1], east=bounds[2], south=bounds[3], north=bounds[4],
        window="common_hours",
        note="Visualization of satellite values sampled at gauge locations, not native satellite " *
            "fields. No metric is computed from these surfaces.",
    ))

    println("\nSpatial agreement over the common window (gauge-pixel pairs):")
    show(select(metrics[metrics.window .== "common_hours", :],
        :event_day, :event_type, :product, :hours, :n, :r, :rho, :RMSE, :Bias, :CSI_50, :cv_ratio, :centroid_shift_km);
        allrows=true, allcols=true)
    println("\n\nWrote $(OUTDIR)")
    return OUTDIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
