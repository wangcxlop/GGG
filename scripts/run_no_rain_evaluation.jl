#!/usr/bin/env julia

# No-rain evaluation, part 1: what FY4B, GPM and GSMaP report on the station-hours the gauges call
# dry, whether that rain is really false, and how dry spells and dry days compare.
#
#   julia --project=. scripts/run_no_rain_evaluation.jl
#
# Same gauges, products, window and samples as scripts/run_satellite_temporal_evaluation.jl, whose
# gauge-dry numbers this reproduces before writing anything. It adds two inputs: ERA5-Land 2 m
# temperature and dew point at the gauges, and the uncalibrated IMERG series
# (output/gpm_imerg_uncal_stations_2022_2024/). CSVs go to output/no_rain_evaluation/. Part 2, the
# benchmark's methods on dry hours, is scripts/run_no_rain_fusion_evaluation.jl; both are drawn by
# `py -3.13 scripts/plot_no_rain_evaluation.py`.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("TableIO", "HeavyRainEvents", "LandformClassification", "SatelliteTemporalEvaluation",
    "ERA5VariableSelection", "NoRainEvaluation")
using Main.TableIO: write_csv_atomic
using Main.HeavyRainEvents: align_to_reference
using Main.LandformClassification: RELIEF_REGIONS
using Main.SatelliteTemporalEvaluation: evaluation_context, sample_mask, false_alarm_table, SEASONS, WET_MM
using Main.ERA5VariableSelection: load_era5_panel
using Main.NoRainEvaluation

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const LANDFORM_FILE = joinpath(ROOT, "data", "processed", "covariates", "station_landform.csv")
const ERA5_FILES = Dict(year => joinpath(ROOT, "data", "processed", "covariates", "era5_land",
    "era5_land_station_hourly_utc_$(year).csv") for year in 2022:2024)
const UNCAL_DIR = joinpath(ROOT, "output", "gpm_imerg_uncal_stations_2022_2024")
const TEMPORAL_FALSE_ALARMS = joinpath(ROOT, "output", "satellite_temporal_evaluation", "false_alarm_metrics.csv")
const OUTDIR = joinpath(ROOT, "output", "no_rain_evaluation")
const ANALYSIS_START = DateTime(2022, 1, 1, 9)
const ANALYSIS_END = DateTime(2025, 1, 1, 8)
const PRODUCT_FILES = [
    "FY4B" => "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
    "GPM" => "hubei_gpm_hourly_2022_2024_full_aligned.csv",
    "GSMaP" => "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
]
# The first two are the temporal report's samples; the third pairs calibrated and uncalibrated IMERG
# on the cells both have, so adding the uncalibrated series never changes the other two.
const SAMPLES = ["all_products" => ["FY4B", "GPM", "GSMaP"], "gpm_gsmap_full" => ["GPM", "GSMaP"],
    "gpm_cal_uncal" => ["GPM", "GPM_uncal"]]
const TEMPORAL_SAMPLES = ("all_products", "gpm_gsmap_full")
# Hourly spells need a contiguous record; FY4B never has the 23:00 label, so only these two qualify.
const SPELL_SAMPLES = ("gpm_gsmap_full", "gpm_cal_uncal")
# A day counts when this many of its 24 hours are scored for both sides (FY4B lacks 23:00 every day).
const MIN_DAY_HOURS = Dict("all_products" => 20, "gpm_gsmap_full" => 24, "gpm_cal_uncal" => 24)
const DAILY_THRESHOLDS = [0.1, 1.0]
const BOOTSTRAP_REPS = 1000
const SEED = 20260919
const NEIGHBOUR_RADIUS_KM = 15.0
const NETWORK_MIN_REPORTING = 0.9
const SIDE_WINDOW_H = 3
const DPD_EDGES = [-Inf, 1.0, 3.0, 6.0, 10.0, Inf]
const DPD_LEVELS = ["< 1", "1-3", "3-6", "6-10", ">= 10"]
const T2M_EDGES = [-Inf, 0.0, 5.0, 15.0, 25.0, Inf]
const T2M_LEVELS = ["< 0", "0-5", "5-15", "15-25", ">= 25"]
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
# Silent-gauge screen: a run of >= 7 complete days < 1 mm at a gauge is flagged when the median of its
# 3 nearest reporting gauges reaches 5 mm on >= 3 of those days. 2 and 5 days are reported as sensitivity.
const SILENT_RULE = (; dry_mm=1.0, wet_mm=5.0, min_days=7, min_wet_days=3, k=3)
const SILENT_SENSITIVITY = [2, 3, 5]
const SCREEN_LEVELS = ["gauge working", "silent-gauge spell"]

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
    lonlat = build_X_lonlat(load_station_meta(joinpath(STUDY_DATA, "station_meta.csv")), obs_ids)
    return obs_times, obs_ids, Y_obs, products, landform, lonlat
end

"""The uncalibrated IMERG monthly long files on the gauge grid (UTC hour start + 9 h = BJT hour ending)."""
function load_uncalibrated_gpm(obs_times, obs_ids)
    files = sort(filter(f -> startswith(f, "gpm_hubei_hourly_long_") && endswith(f, ".csv"), readdir(UNCAL_DIR)))
    length(files) == 36 || error("expected 36 monthly uncalibrated IMERG files in $UNCAL_DIR, found $(length(files))")
    time, station, value = DateTime[], String[], Float64[]
    for file in files
        table = CSV.read(joinpath(UNCAL_DIR, file), DataFrame; types=Dict(:station_id => String, :time => String))
        append!(time, DateTime.(first.(table.time, 19)))
        append!(station, table.station_id)
        append!(value, [ismissing(v) ? NaN : Float64(v) for v in table.gpm_mm_h])
    end
    wide = long_to_wide_hourly(time, station, value, obs_times, obs_ids; offset_hours=9)
    wide.matched == length(obs_ids) * length(obs_times) ||
        error("uncalibrated IMERG fills $(wide.matched) of $(length(obs_ids) * length(obs_times)) gauge cells")
    return wide.Y
end

"""Every stratifier as `name => (labels, levels)`; `labels` is a matrix or a `(i, j) -> label` function."""
function stratifiers(ctx, obs_ids, Y_obs, lonlat, era5, silent)
    n_station = length(obs_ids)
    month_of_column = month.(ctx.days)[ctx.column_day]
    before, after = hours_to_rain(Y_obs)
    neighbour_share, _ = neighbour_wet_share(Y_obs, lonlat; radius_km=NEIGHBOUR_RADIUS_KM)
    network_share = network_wet_share(Y_obs; min_reporting=NETWORK_MIN_REPORTING)
    t2m, dpd = era5[:t2m_c], era5[:t2m_c] .- era5[:d2m_c]
    dpd_labels = bin_labels(dpd, DPD_EDGES)
    season_dpd = [dpd_labels[i, j] == 0 ? 0 : (ctx.column_season[j] - 1) * length(DPD_LEVELS) + dpd_labels[i, j]
        for i in 1:n_station, j in eachindex(ctx.times)]
    return [
        "all" => ((i, j) -> 1, ["all"]),
        "season" => ((i, j) -> ctx.column_season[j], SEASONS),
        "month" => ((i, j) -> month_of_column[j], MONTHS),
        "hour" => ((i, j) -> ctx.column_hour[j] + 1, [lpad(h, 2, '0') for h in 0:23]),
        "season_hour" => ((i, j) -> (ctx.column_season[j] - 1) * 24 + ctx.column_hour[j] + 1,
            ["$(s)|$(lpad(h, 2, '0'))" for s in SEASONS for h in 0:23]),
        "region" => ((i, j) -> ctx.station_region_index[i], RELIEF_REGIONS),
        "station" => ((i, j) -> i, obs_ids),
        "rain_proximity" => (proximity_labels(before, after), PROXIMITY_LEVELS),
        "rain_side" => (side_labels(before, after; window=SIDE_WINDOW_H), SIDE_LEVELS),
        "neighbours" => (neighbour_labels(neighbour_share, Y_obs), NEIGHBOUR_LEVELS),
        "network" => (network_labels(network_share, Y_obs), NETWORK_LEVELS),
        "dewpoint_depression" => (dpd_labels, DPD_LEVELS),
        "t2m" => (bin_labels(t2m, T2M_EDGES), T2M_LEVELS),
        "season_dewpoint" => (season_dpd, ["$(s)|$(b)" for s in SEASONS for b in DPD_LEVELS]),
        "gauge_screen" => ((i, j) -> silent[i, ctx.column_day[j]] ? 2 : 1, SCREEN_LEVELS),
        "season_screened" => ((i, j) -> silent[i, ctx.column_day[j]] ? 0 : ctx.column_season[j], SEASONS),
    ], network_share
end

"""How many of the *other* products in the sample report >= 0.1 on each gauge-dry cell (+1; 0 = gauge wet)."""
function other_products_labels(Y_obs, products, members, product)
    others = [products[p] for p in members if p != product]
    return [Y_obs[k] < WET_MM ? 1 + count(Y -> Y[k] >= WET_MM, others) : 0 for k in CartesianIndices(Y_obs)]
end

"""Reproduce the temporal report's gauge-dry table from `occurrence_table` rows, then check it against the CSV."""
function check_against_temporal(ctx, Y_obs, products, masks, occurrence)
    reference = CSV.read(TEMPORAL_FALSE_ALARMS, DataFrame)
    for (sample, members) in SAMPLES
        sample in TEMPORAL_SAMPLES || continue
        for product in members
            table = false_alarm_table(ctx, Y_obs, products[product], masks[sample]; sample, product)
            for season in vcat(["all"], SEASONS)
                expected = only(table[(table.season .== season) .& (table.region .== "all"), :])
                published = only(reference[(reference.sample .== sample) .& (reference.product .== product) .&
                    (reference.season .== season) .& (reference.region .== "all"), :])
                stratifier, level = season == "all" ? ("all", "all") : ("season", season)
                rows = occurrence[(occurrence.sample .== sample) .& (occurrence.product .== product) .&
                    (occurrence.stratifier .== stratifier) .& (occurrence.level .== level), :]
                at01 = only(rows[rows.threshold .== 0.1, :])
                at05 = only(rows[rows.threshold .== 0.5, :])
                for (ours, theirs) in ((at01.n_dry, expected.n_gauge_dry), (at01.POFD, expected.false_alarm_rate),
                        (at05.POFD, expected.false_alarm_rate_res), (at01.FAR, expected.FAR),
                        (at05.FAR, expected.FAR_res), (at01.dry_volume_share, expected.dry_volume_share))
                    isapprox(ours, theirs; rtol=1e-10) || error("$sample/$product/$season: $ours vs false_alarm_table $theirs")
                end
                for column in (:n_gauge_dry, :false_alarm_rate, :false_alarm_rate_res, :FAR, :FAR_res, :dry_volume_share)
                    isapprox(expected[column], published[column]; rtol=1e-9) ||
                        error("$sample/$product/$season/$column differs from the published temporal report")
                end
            end
        end
    end
    println("  gate: every gauge-dry count and rate matches false_alarm_metrics.csv")
end

"""Calibrated against uncalibrated IMERG on their shared cells: the clock check and the zero pattern."""
function calibration_check(cal, uncal, mask)
    lags = lagged_correlation(cal, uncal, mask; lags=-3:3)
    n = both_zero = cal_zero_only = uncal_zero_only = both_wet = cal_wet_only = uncal_wet_only = 0
    volume_cal = volume_uncal = 0.0
    for k in eachindex(cal)
        mask[k] || continue
        a, b = cal[k], uncal[k]
        n += 1
        both_zero += a == 0 && b == 0
        cal_zero_only += a == 0 && b > 0
        uncal_zero_only += a > 0 && b == 0
        both_wet += a >= 0.1 && b >= 0.1
        cal_wet_only += a >= 0.1 && b < 0.1
        uncal_wet_only += a < 0.1 && b >= 0.1
        volume_cal += a
        volume_uncal += b
    end
    pattern = DataFrame(n=[n], both_zero=[both_zero / n], zero_in_calibrated_only=[cal_zero_only / n],
        zero_in_uncalibrated_only=[uncal_zero_only / n], both_wet=[both_wet / n], wet_in_calibrated_only=[cal_wet_only / n],
        wet_in_uncalibrated_only=[uncal_wet_only / n], volume_ratio_calibrated_to_uncalibrated=[volume_cal / volume_uncal])
    return lags, pattern
end

"""Every screened dry run with its dates and what the products recorded over it."""
function silent_spell_table(screen, ctx, obs_ids, products, masks)
    table = copy(screen.spells)
    row_of = Dict(id => i for (i, id) in enumerate(obs_ids))
    table.first_date = ctx.days[table.first_day]
    table.last_date = ctx.days[table.last_day]
    table.season = [SEASONS[ctx.day_season[d]] for d in table.first_day]
    for product in ("GPM", "GSMaP")
        D = daily_sums(ctx.column_day, products[product], masks["gpm_gsmap_full"]; min_hours=24)
        table[!, Symbol(product, "_total")] = [sum(filter(isfinite, D[row_of[r.station_id], r.first_day:r.last_day]); init=0.0)
            for r in eachrow(table)]
    end
    return select!(table, :station_id, :first_date, :last_date, :season, :n_days, :gauge_total, :neighbour_total,
        :neighbour_wet_days, :GPM_total, :GSMaP_total, :flagged, :first_day, :last_day)
end

"""How the screen's reach depends on the number of contradicted days it asks for."""
function silent_sensitivity(D_full, lonlat, obs_ids, ctx, Y_obs, mask)
    rows = NamedTuple[]
    dry = mask .& (Y_obs .< WET_MM)
    n_dry = count(dry)
    for min_wet_days in SILENT_SENSITIVITY
        screen = silent_gauge_spells(D_full, lonlat, obs_ids; merge(SILENT_RULE, (; min_wet_days))...)
        in_spell = count(k -> dry[k] && screen.flagged[k[1], ctx.column_day[k[2]]], CartesianIndices(dry))
        push!(rows, (; min_wet_days, n_spells=count(screen.spells.flagged), n_stations=length(unique(
            screen.spells.station_id[screen.spells.flagged])), station_days=count(screen.flagged),
            dry_hours_in_spells=in_spell, dry_hour_share=in_spell / n_dry, p_dry_given_wet=screen.p_dry_given_wet))
    end
    return DataFrame(rows)
end

function main()
    println("Reading gauges, products, landform regions and station coordinates ...")
    obs_times, obs_ids, Y_obs, products, landform, lonlat = load_inputs()
    ctx = evaluation_context(obs_times, landform.relief_region; regions=RELIEF_REGIONS)
    println("Reading the uncalibrated IMERG series ...")
    products["GPM_uncal"] = load_uncalibrated_gpm(obs_times, obs_ids)
    println("Reading ERA5-Land 2 m temperature and dew point ...")
    era5 = load_era5_panel(ERA5_FILES, obs_ids, obs_times; variables=[:t2m_c, :d2m_c], offset_hours=9).values
    masks = Dict(sample => sample_mask(Y_obs, (products[p] for p in members)...) for (sample, members) in SAMPLES)
    println("  $(length(obs_ids)) gauges x $(length(obs_times)) hours; cells per sample: ",
        join(["$s $(count(masks[s]))" for (s, _) in SAMPLES], ", "))

    lags, zero_pattern = calibration_check(products["GPM"], products["GPM_uncal"], masks["gpm_cal_uncal"])
    best_lag = lags.lag[argmax(lags.r)]
    println("  calibrated vs uncalibrated IMERG: r by lag ", join(["$(l)h $(round(r; digits=3))" for (l, r) in zip(lags.lag, lags.r)], ", "))
    best_lag == 0 || error("uncalibrated IMERG peaks at lag $(best_lag) h against the calibrated series; its clock is wrong")

    println("Silent-gauge screen ...")
    D_full = daily_sums(ctx.column_day, Y_obs, isfinite.(Y_obs); min_hours=24)
    screen = silent_gauge_spells(D_full, lonlat, obs_ids; SILENT_RULE...)
    silent_table = silent_spell_table(screen, ctx, obs_ids, products, masks)
    sensitivity = silent_sensitivity(D_full, lonlat, obs_ids, ctx, Y_obs, masks["gpm_gsmap_full"])
    println("  $(count(screen.spells.flagged)) flagged spells, $(count(screen.flagged)) station-days; ",
        "P(gauge < 1 mm | neighbours >= 5 mm) = $(round(screen.p_dry_given_wet; digits=3))")

    println("Occurrence and dry-hour tables ...")
    strata, network_share = stratifiers(ctx, obs_ids, Y_obs, lonlat, era5, screen.flagged)
    occurrence, exceedance, extent = DataFrame[], DataFrame[], NamedTuple[]
    for (s_index, (sample, members)) in enumerate(SAMPLES), (p_index, product) in enumerate(members)
        mask, Y = masks[sample], products[product]
        for (k, (name, (labels, levels))) in enumerate(strata)
            reps = name == "station" ? 0 : BOOTSTRAP_REPS
            push!(occurrence, occurrence_table(labels, levels, Y_obs, Y, mask, ctx.column_day; reps,
                seed=SEED + 10_000 * s_index + 1_000 * p_index + k, key=(; sample, product, stratifier=name)))
        end
        if sample == "all_products"
            labels = other_products_labels(Y_obs, products, members, product)
            push!(occurrence, occurrence_table(labels, ["0 others wet", "1 other wet", "2 others wet"], Y_obs, Y, mask,
                ctx.column_day; reps=BOOTSTRAP_REPS, seed=SEED + 10_000 * s_index + 1_000 * p_index + 999,
                key=(; sample, product, stratifier="other_products")))
        end
        for (name, (labels, levels)) in strata[1:2]
            push!(exceedance, occurrence_table(labels, levels, Y_obs, Y, mask, ctx.column_day;
                thresholds=EXCEEDANCE_THRESHOLDS, reps=0, key=(; sample, product, stratifier=name)))
        end
        push!(extent, phantom_extent(network_share, Y, mask; key=(; sample, product)))
        println("  $sample / $product")
    end
    occurrence, exceedance = vcat(occurrence...), vcat(exceedance...)
    check_against_temporal(ctx, Y_obs, products, masks, occurrence)

    println("Neighbouring gauges as an estimate ...")
    baseline = vcat([neighbour_gauge_baseline(Y_obs, lonlat, masks[sample], ctx.column_day; reps=BOOTSTRAP_REPS,
        seed=SEED + 50_000 + k, key=(; sample)) for (k, sample) in enumerate(TEMPORAL_SAMPLES)]...)

    println("Dry spells and dry days ...")
    hourly_summary, hourly_hist = DataFrame[], DataFrame[]
    for (sample, members) in SAMPLES
        sample in SPELL_SAMPLES || continue
        sources = vcat(["Gauge" => (Y_obs, [0.1])], [p => (products[p], DRY_THRESHOLDS) for p in members])
        for (source, (Y, thresholds)) in sources, threshold in thresholds
            spells = dry_spell_table(Y, masks[sample]; threshold, station_ids=obs_ids, key=(; sample, source, threshold))
            push!(hourly_summary, spells.summary)
            push!(hourly_hist, spells.histogram)
        end
    end
    daily, daily_summary, daily_hist = DataFrame[], DataFrame[], DataFrame[]
    day_season = ctx.day_season
    for (s_index, (sample, members)) in enumerate(SAMPLES)
        min_hours = MIN_DAY_HOURS[sample]
        D_obs = daily_sums(ctx.column_day, Y_obs, masks[sample]; min_hours)
        day_axis = collect(1:size(D_obs, 2))
        for (p_index, product) in enumerate(members)
            D = daily_sums(ctx.column_day, products[product], masks[sample]; min_hours)
            valid = isfinite.(D_obs) .& isfinite.(D)
            for (t_index, threshold) in enumerate(DAILY_THRESHOLDS)
                for (name, labels, levels) in (("all", (i, d) -> 1, ["all"]), ("season", (i, d) -> day_season[d], SEASONS))
                    push!(daily, occurrence_table(labels, levels, D_obs, D, valid, day_axis; thresholds=[threshold],
                        wet=threshold, reps=BOOTSTRAP_REPS, seed=SEED + 70_000 + 1_000 * s_index + 100 * p_index + 10 * t_index,
                        key=(; sample, product, min_day_hours=min_hours, stratifier=name)))
                end
            end
            for (source, Y) in (("Gauge", D_obs), (product, D)), threshold in DAILY_THRESHOLDS
                spells = dry_spell_table(Y, valid; threshold, station_ids=obs_ids,
                    key=(; sample, pair=product, source, threshold))
                push!(daily_summary, spells.summary)
                push!(daily_hist, spells.histogram)
            end
        end
    end

    stations = DataFrame(station_id=obs_ids, lon=lonlat[:, 1], lat=lonlat[:, 2], elevation_m=landform.elevation_m,
        relief_m=landform.relief_m, relief_region=landform.relief_region)
    network_hours = DataFrame(n_hours=[length(network_share)], n_hours_reporting=[count(isfinite, network_share)],
        n_network_dry=[count(==(0), network_share)])

    mkpath(OUTDIR)
    outputs = [
        "dry_occurrence.csv" => occurrence, "dry_exceedance.csv" => exceedance,
        "phantom_extent.csv" => DataFrame(extent), "network_hours.csv" => network_hours,
        "neighbour_gauge_baseline.csv" => baseline,
        "dry_spells_hourly_summary.csv" => vcat(hourly_summary...), "dry_spells_hourly_histogram.csv" => vcat(hourly_hist...),
        "dry_days.csv" => vcat(daily...), "dry_spells_daily_summary.csv" => vcat(daily_summary...),
        "dry_spells_daily_histogram.csv" => vcat(daily_hist...),
        "gpm_calibration_lags.csv" => lags, "gpm_calibration_zero_pattern.csv" => zero_pattern,
        "silent_gauge_spells.csv" => silent_table, "silent_gauge_sensitivity.csv" => sensitivity,
        "stations.csv" => stations,
    ]
    for (name, table) in outputs
        write_csv_atomic(joinpath(OUTDIR, name), table)
    end
    write_csv_atomic(joinpath(OUTDIR, "run_settings.csv"), DataFrame(
        setting=["analysis_start", "analysis_end", "samples", "wet_threshold_mm_h", "estimate_thresholds_mm_h",
            "exceedance_thresholds_mm_h", "bootstrap_reps", "bootstrap_block", "seed", "neighbour_radius_km",
            "network_min_reporting", "rain_side_window_h", "proximity_edges_h", "dewpoint_depression_edges_c",
            "t2m_edges_c", "daily_thresholds_mm", "min_day_hours", "era5_offset_hours", "uncalibrated_imerg_offset_hours",
            "silent_gauge_rule"],
        value=[string(ANALYSIS_START), string(ANALYSIS_END), join(["$s=" * join(m, "+") for (s, m) in SAMPLES], "; "),
            "0.1", join(DRY_THRESHOLDS, ","), join(EXCEEDANCE_THRESHOLDS, ","), string(BOOTSTRAP_REPS), "08-08 BJT met-day",
            string(SEED), string(NEIGHBOUR_RADIUS_KM), string(NETWORK_MIN_REPORTING), string(SIDE_WINDOW_H),
            join(PROXIMITY_EDGES, ","), join(DPD_EDGES, ","), join(T2M_EDGES, ","), join(DAILY_THRESHOLDS, ","),
            join(["$k=$v" for (k, v) in MIN_DAY_HOURS], "; "), "9", "9",
            join(["$k=$v" for (k, v) in pairs(SILENT_RULE)], "; ")],
    ))

    pooled = occurrence[occurrence.stratifier .== "all", :]
    println("\nGauge-dry hours, all seasons and gauges:")
    show(select(pooled, :sample, :product, :threshold, :n_dry, :dry_share, :POFD, :FAR, :NPV, :freq_bias, :HSS,
        :mean_dry, :spurious_mm_per_year, :dry_volume_share, :zero_share_dry); allrows=true, allcols=true)
    println("\n\nWrote $(OUTDIR)")
    return OUTDIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
