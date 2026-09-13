#!/usr/bin/env julia

# Temporal evaluation: do FY4B, GPM and GSMaP capture how rainfall evolves hour by hour, and does that
# depend on rain intensity, season and landform region?
#
#   julia --project=. scripts/prepare_station_landform.jl        # once: the landform regions
#   julia --project=. scripts/run_satellite_temporal_evaluation.jl
#
# Scores 2022-2024 (08-08 BJT days) in two samples: `all_products`, the cells where the gauge and all
# three products are finite (FY4B's coverage), and `gpm_gsmap_full`, the full record for GPM and GSMaP.
# Writes CSVs to output/satellite_temporal_evaluation/; the figures are drawn from those by
# `py -3.13 scripts/plot_satellite_temporal_evaluation.py`.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("TableIO", "HeavyRainEvents", "LandformClassification", "SatelliteTemporalEvaluation")
using Main.TableIO: write_csv_atomic
using Main.HeavyRainEvents: align_to_reference
using Main.LandformClassification: RELIEF_REGIONS
using Main.SatelliteTemporalEvaluation

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const LANDFORM_FILE = joinpath(ROOT, "data", "processed", "covariates", "station_landform.csv")
const OUTDIR = joinpath(ROOT, "output", "satellite_temporal_evaluation")
# The first and last hour-ending labels of the 08-08 BJT days 2022-01-01 .. 2024-12-31.
const ANALYSIS_START = DateTime(2022, 1, 1, 9)
const ANALYSIS_END = DateTime(2025, 1, 1, 8)
const PRODUCT_FILES = [
    "FY4B" => "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
    "GPM" => "hubei_gpm_hourly_2022_2024_full_aligned.csv",
    "GSMaP" => "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
]
const SAMPLES = ["all_products" => ["FY4B", "GPM", "GSMaP"], "gpm_gsmap_full" => ["GPM", "GSMaP"]]
const BOOTSTRAP_REPS = 1000
const SEED = 20260913
const PRIMARY_DRY_GAP_H = 3        # 0.5 mm buckets tip only every few hours in drizzle
const SENSITIVITY_DRY_GAP_H = 1
const EVENT_PAD_H = 3
const EVENT_MAX_LAG_H = 3
const SERIES_MAX_LAG_H = 6
const MIN_STATION_FRACTION = 0.9
const MIN_STATION_HOURS = 50
# An independent count of gauge hours >= 20 mm/h over the analysis window, taken from the raw CSV.
const GAUGE_SEVERE_HOURS_UPPER_BOUND = 1626

function load_inputs()
    obs_times, obs_ids, Y_obs = read_hourly_wide(joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"))
    obs_times, Y_obs = subset_time_window(obs_times, Y_obs, ANALYSIS_START, ANALYSIS_END)
    products = Dict{String,Matrix{Float64}}()
    for (product, file) in PRODUCT_FILES
        times, ids, Y = read_hourly_wide(joinpath(STUDY_DATA, file))
        times, Y = subset_time_window(times, Y, ANALYSIS_START, ANALYSIS_END)
        products[product] = align_to_reference(obs_times, obs_ids, times, ids, Y)
    end
    landform = CSV.read(LANDFORM_FILE, DataFrame; types=Dict(:station_id => String))
    row_of = Dict(id => k for (k, id) in enumerate(landform.station_id))
    absent = filter(id -> !haskey(row_of, id), obs_ids)
    isempty(absent) || error("$(length(absent)) gauges have no landform row, e.g. $(first(absent))")
    landform = landform[[row_of[id] for id in obs_ids], :]
    return obs_times, obs_ids, Y_obs, products, landform
end

"""Tag every region-level row with whether its region has too few gauges to interpret."""
function flag_regions!(table::DataFrame, insufficient::Dict{String,Bool})
    table.insufficient_gauges = [get(insufficient, region, false) for region in table.region]
    return table
end

function check_invariants(coverage, metrics, confusion, event_summaries)
    pooled(t, sample) = t[(t.sample .== sample) .& (t.season .== "all") .& (t.region .== "all"), :]
    for (sample, members) in SAMPLES, product in members
        rows = pooled(metrics, sample)
        rows = rows[rows.product .== product, :]
        wet_cells = only(pooled(coverage, sample).n_gauge_wet)
        only(rows[rows.gauge_class .== "all_wet", :n]) == wet_cells ||
            error("$sample/$product: intensity table n differs from the sample's gauge-wet cells")
        sum(rows[rows.gauge_class .!= "all_wet", :n]) == wet_cells ||
            error("$sample/$product: the five classes do not partition the gauge-wet cells")
    end
    for group in groupby(confusion, [:sample, :product, :season, :region, :gauge_class])
        n = only(metrics[(metrics.sample .== group.sample[1]) .& (metrics.product .== group.product[1]) .&
            (metrics.season .== group.season[1]) .& (metrics.region .== group.region[1]) .&
            (metrics.gauge_class .== group.gauge_class[1]), :n])
        sum(group.n) == n || error("confusion row $(Tuple(group[1, 1:5])) does not sum to the class n")
    end
    severe = pooled(metrics, "gpm_gsmap_full")
    severe = only(severe[(severe.product .== "GPM") .& (severe.gauge_class .== "severe_rainstorm"), :n])
    severe <= GAUGE_SEVERE_HOURS_UPPER_BOUND || error("more severe-rainstorm gauge hours than the raw file holds")
    primary = event_summaries[(event_summaries.min_dry_gap_h .== PRIMARY_DRY_GAP_H) .& (event_summaries.season .== "all") .&
        (event_summaries.region .== "all") .& (event_summaries.event_class .== "all") .& (event_summaries.product .== "GPM"), :]
    only(primary[primary.sample .== "all_products", :n_events]) <= only(primary[primary.sample .== "gpm_gsmap_full", :n_events]) ||
        error("the three-product sample scored more GPM events than the GPM/GSMaP sample")
    println("  invariants hold (class partition, confusion sums, severe-hour bound, event retention)")
end

"""Print one scored heavy event with its hourly values, so the event scores can be checked by hand."""
function show_example_event(obs_times, obs_ids, Y_obs, products, events, scores)
    candidates = scores[(scores.product .== "GPM") .& (scores.event_class .== "rainstorm") .& scores.detected, :]
    isempty(candidates) && return
    score = candidates[1:1, :]
    event = only(events[events.event_id .== only(score.event_id), :])
    row = findfirst(==(event.station_id), obs_ids)
    columns = findfirst(==(event.window_start), obs_times):findfirst(==(event.window_stop), obs_times)
    println("\nHand-check: event $(event.event_id) at $(event.station_id), $(event.start_time) .. $(event.stop_time)")
    for j in columns
        println("  $(obs_times[j])  gauge $(rpad(Y_obs[row, j], 6))  GPM $(round(products["GPM"][row, j]; digits=2))")
    end
    println("  scores: ", select(score, :peak_error_h, :centroid_error_h, :onset_error_h, :end_error_h,
        :volume_rel_bias, :peak_ratio, :r_event, :best_lag_h)[1, :])
end

function main()
    println("Reading hourly tables and landform regions ...")
    obs_times, obs_ids, Y_obs, products, landform = load_inputs()
    ctx = evaluation_context(obs_times, landform.relief_region; regions=RELIEF_REGIONS)
    insufficient = Dict(String(r.relief_region) => r.insufficient_gauges for r in eachrow(landform))
    println("  $(length(obs_ids)) gauges x $(length(obs_times)) hours; regions: ",
        join(["$r ($(count(==(r), landform.relief_region)))" for r in RELIEF_REGIONS], ", "))
    masks = Dict(sample => sample_mask(Y_obs, (products[p] for p in members)...) for (sample, members) in SAMPLES)

    println("Coverage, intensity classes and false alarms ...")
    coverage = vcat([sample_coverage_table(ctx, Y_obs, masks[s]; sample=s) for (s, _) in SAMPLES]...)
    metrics, confusion, false_alarms = DataFrame[], DataFrame[], DataFrame[]
    for (sample, members) in SAMPLES, product in members
        m, c = intensity_class_table(ctx, Y_obs, products[product], masks[sample]; sample, product, reps=BOOTSTRAP_REPS, seed=SEED)
        push!(metrics, m)
        push!(confusion, c)
        push!(false_alarms, false_alarm_table(ctx, Y_obs, products[product], masks[sample]; sample, product))
    end
    metrics, confusion, false_alarms = vcat(metrics...), vcat(confusion...), vcat(false_alarms...)

    println("Station rain events ...")
    event_summaries = DataFrame[]
    primary_events = nothing
    for gap in (PRIMARY_DRY_GAP_H, SENSITIVITY_DRY_GAP_H)
        tables = event_tables(ctx, Y_obs, products, SAMPLES, obs_ids; min_dry_gap=gap, pad=EVENT_PAD_H, max_lag=EVENT_MAX_LAG_H)
        summary = event_summary(ctx, tables.events, tables.scores, SAMPLES;
            reps=BOOTSTRAP_REPS, seed=SEED, bootstrap=gap == PRIMARY_DRY_GAP_H)
        insertcols!(summary, 1, :min_dry_gap_h => gap)
        push!(event_summaries, summary)
        gap == PRIMARY_DRY_GAP_H && (primary_events = tables)
        println("  dry gap $(gap) h: $(nrow(tables.events)) gauge events, $(count(tables.events.gauge_complete)) complete")
    end
    event_summaries = vcat(event_summaries...)

    println("Diurnal cycles, station correlations and regional series ...")
    diurnal = vcat([diurnal_cycle_table(ctx, ["Gauge" => Y_obs; [p => products[p] for p in members]], masks[sample]; sample)
        for (sample, members) in SAMPLES]...)
    diurnal_scores = diurnal_summary(diurnal, first.(PRODUCT_FILES))
    stations = vcat([station_correlation_table(ctx, Y_obs, products[product], masks[sample]; sample, product,
        station_ids=obs_ids, min_n=MIN_STATION_HOURS) for (sample, members) in SAMPLES for product in members]...)
    station_summary = station_correlation_summary(stations)
    regional = vcat([regional_series_table(ctx, Y_obs, Dict(p => products[p] for p in members), masks[sample];
        sample, max_lag=SERIES_MAX_LAG_H, min_station_fraction=MIN_STATION_FRACTION) for (sample, members) in SAMPLES]...)

    check_invariants(coverage, metrics, confusion, event_summaries)
    show_example_event(obs_times, obs_ids, Y_obs, products, primary_events.events, primary_events.scores)

    for table in (coverage, metrics, confusion, false_alarms, primary_events.events, primary_events.scores,
            event_summaries, diurnal, diurnal_scores, stations, station_summary, regional)
        flag_regions!(table, insufficient)
    end
    outputs = [
        "sample_coverage.csv" => coverage, "intensity_class_metrics.csv" => metrics,
        "intensity_confusion.csv" => confusion, "false_alarm_metrics.csv" => false_alarms,
        "station_events.csv" => primary_events.events, "event_scores.csv" => primary_events.scores,
        "event_timing_summary.csv" => event_summaries, "diurnal_cycle.csv" => diurnal,
        "diurnal_summary.csv" => diurnal_scores, "station_hourly_correlation.csv" => stations,
        "station_correlation_summary.csv" => station_summary, "regional_series_metrics.csv" => regional,
    ]
    for (name, table) in outputs
        write_csv_atomic(joinpath(OUTDIR, name), table)
    end
    write_csv_atomic(joinpath(OUTDIR, "run_settings.csv"), DataFrame(
        setting=["analysis_start", "analysis_end", "samples", "intensity_bounds_mm_h", "wet_threshold_mm_h",
            "gauge_resolution_mm", "seasons", "regions", "relief_window_km", "min_region_gauges_flag",
            "bootstrap_reps", "bootstrap_block", "seed", "event_dry_gap_h", "event_dry_gap_sensitivity_h",
            "event_pad_h", "event_max_lag_h", "series_max_lag_h", "regional_min_station_fraction",
            "station_correlation_min_hours"],
        value=[string(ANALYSIS_START), string(ANALYSIS_END),
            join(["$s=" * join(m, "+") for (s, m) in SAMPLES], "; "), join(INTENSITY_BOUNDS, ","), string(WET_MM),
            string(GAUGE_RESOLUTION_MM), join(SEASONS, ","), join(RELIEF_REGIONS, ","), string(first(landform.window_side_km)),
            "20 gauges", string(BOOTSTRAP_REPS), "08-08 BJT met-day", string(SEED), string(PRIMARY_DRY_GAP_H),
            string(SENSITIVITY_DRY_GAP_H), string(EVENT_PAD_H), string(EVENT_MAX_LAG_H), string(SERIES_MAX_LAG_H),
            string(MIN_STATION_FRACTION), string(MIN_STATION_HOURS)],
    ))

    pooled(t) = t[(t.season .== "all") .& (t.region .== "all"), :]
    println("\nIntensity classes, all seasons and regions:")
    show(select(pooled(metrics), :sample, :product, :gauge_class, :n, :RB_pct, :POD_rain, :class_hit, :under_class, :over_class);
        allrows=true, allcols=true)
    println("\n\nFalse alarms, all seasons and regions:")
    show(select(pooled(false_alarms), :sample, :product, :false_alarm_rate, :false_alarm_rate_res, :FAR, :dry_volume_share);
        allrows=true, allcols=true)
    primary = event_summaries[(event_summaries.min_dry_gap_h .== PRIMARY_DRY_GAP_H), :]
    println("\n\nEvents, all seasons and regions (dry gap $(PRIMARY_DRY_GAP_H) h):")
    show(select(pooled(primary), :sample, :product, :event_class, :n_events, :retained_share, :POD_event,
        :median_peak_error_h, :peak_within_1h, :pooled_volume_rel_bias, :median_r_event); allrows=true, allcols=true)
    println("\n\nDiurnal cycle, all regions:")
    show(select(diurnal_scores[diurnal_scores.region .== "all", :], :sample, :product, :season, :r_amount,
        :peak_hour_gauge, :peak_hour_diff_h, :phase_diff_h, :rel_amplitude_ratio); allrows=true, allcols=true)
    println("\n\nWrote $(OUTDIR)")
    return OUTDIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
