#!/usr/bin/env julia

# No-rain evaluation, part 2: how the interpolation benchmark's methods behave on the gauge-dry hours.
#
#   julia --project=. scripts/run_no_rain_fusion_evaluation.jl
#   julia --project=. scripts/run_no_rain_fusion_evaluation.jl <benchmark_output_dir>
#
# Reads the held-out predictions a finished `run_interpolation_benchmark.jl full --nested-covariates
# --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor` run saved (`oof_<method>.csv`
# for every anchor under both CV schemes) and re-scores them on dry hours; it fits nothing. Before
# writing, every method's dry-hour n, RMSE and bias are checked against the run's own
# `metrics_pooled.csv`. CSVs go to output/no_rain_evaluation/fusion/; part 1 is
# scripts/run_no_rain_evaluation.jl, and `py -3.13 scripts/plot_no_rain_evaluation.py` draws both.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("BenchmarkDiagnostics", "TableIO", "NoRainEvaluation")
using Main.BenchmarkDiagnostics: read_wide_matrix, read_mask_matrix, read_fold_map
using Main.TableIO: write_csv_atomic
using Main.HeavyRainEvents: met_day
using Main.SatelliteTemporalEvaluation: season_of, SEASONS, WET_MM
using Main.NoRainEvaluation

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(
    ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_" *
        "blendagrenv_fusedanchor",
)
const OUTDIR = joinpath(ROOT, "output", "no_rain_evaluation", "fusion")
const SCHEMES = ["balanced_spatial", "random"]
const PRIMARY_SCHEME = "balanced_spatial"
const ANCHORS = ["FY4B", "GPM", "GSMaP", "MERGED_MEAN", "MERGED_OLS"]
const REFERENCE = "adw"
# The methods broken down further; the summary and the paired comparison cover every saved method.
const PANEL = ["raw", "adw", "idw", "tps", "gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
const DISTANCE_EDGES = [0.0, 20.0, 50.0, 100.0, Inf]
const DISTANCE_LEVELS = ["0-20 km", "20-50 km", "50-100 km", ">= 100 km"]
const QUADRANTS = ["gauge dry, anchor dry", "gauge dry, anchor wet", "gauge wet, anchor dry", "gauge wet, anchor wet"]
const ZERO_TAUS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
const DAILY_THRESHOLDS = [0.1, 1.0]
const MIN_DAY_HOURS = 20          # the grid never has 23:00, so no day is complete
const BOOTSTRAP_REPS = 1000
const SEED = 20260919
const NEIGHBOUR_RADIUS_KM = 15.0
const SIDE_WINDOW_H = 3
# The same silent-gauge screen as scripts/run_no_rain_evaluation.jl, on the full contiguous record.
const SILENT_RULE = (; dry_mm=1.0, wet_mm=5.0, min_days=7, min_wet_days=3, k=3)
const SCREEN_LEVELS = ["gauge working", "silent-gauge spell"]

"""The `MGERConfig` the full benchmark run used, so the common station/time grid matches."""
function study_config(outdir::AbstractString)
    return MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(STUDY_DATA, "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv"),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2024_full_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2024_full_aligned.csv"),
        ),
        outdir=outdir,
        rain_threshold=0.1,
        analysis_start=DateTime(2022, 1, 1, 9),
        analysis_end=DateTime(2025, 1, 1, 8),
        expected_common_time_count=13471,
    )
end

"""
Gauge-side stratifiers on the benchmark grid. Distances to rain and the silent-gauge screen are computed on
the full contiguous gauge record and then picked out at the grid's hours: on the grid itself a column step can
span many hours, and no day is complete.
"""
function gauge_strata(mger, ids, grid_times, lonlat)
    times, full_ids, Y = read_hourly_wide(mger.obs_hourly_wide_path)
    times, Y = subset_time_window(times, Y, mger.analysis_start, mger.analysis_end)
    all(==(Hour(1)), diff(times)) || error("the gauge record is not contiguous")
    row_of = Dict(id => k for (k, id) in enumerate(full_ids))
    column_of = Dict(t => j for (j, t) in enumerate(times))
    rows, columns = [row_of[id] for id in ids], [column_of[t] for t in grid_times]
    Y = Y[rows, :]
    before, after = hours_to_rain(Y)
    share, _ = neighbour_wet_share(Y, lonlat; radius_km=NEIGHBOUR_RADIUS_KM)
    days = met_day.(times)
    day_number = Dict(d => k for (k, d) in enumerate(unique(days)))
    full_day = [day_number[d] for d in days]
    screen = silent_gauge_spells(daily_sums(full_day, Y, isfinite.(Y); min_hours=24), lonlat, ids; SILENT_RULE...)
    pick(M) = M[:, columns]
    Y_grid = pick(Y)
    silent = BitMatrix(screen.flagged[:, full_day[columns]])
    return [
        "rain_proximity" => (proximity_labels(pick(before), pick(after)), PROXIMITY_LEVELS),
        "rain_side" => (side_labels(pick(before), pick(after); window=SIDE_WINDOW_H), SIDE_LEVELS),
        "neighbours" => (neighbour_labels(pick(share), Y_grid), NEIGHBOUR_LEVELS),
        "gauge_screen" => ((i, j) -> silent[i, j] ? 2 : 1, SCREEN_LEVELS),
    ], Y_grid, silent
end

"""
The evaluation mask and the saved method names of one `(scheme, anchor)` directory, checked to sit on the
grid's hours. Predictions are read one method at a time by the caller: all twenty at once would hold
about 0.5 GB per anchor for no benefit.
"""
function open_anchor(scheme_dir, anchor, ids, grid_times)
    dir = joinpath(scheme_dir, lowercase(anchor))
    stamps = CSV.read(joinpath(dir, "common_evaluation_mask.csv"), DataFrame; select=[:time], types=Dict(:time => String))
    DateTime.(first.(stamps.time, 19)) == grid_times || error("$dir: the saved hours are not the benchmark grid")
    mask = read_mask_matrix(joinpath(dir, "common_evaluation_mask.csv"), ids)
    methods = sort([f[5:end-4] for f in readdir(dir) if startswith(f, "oof_") && endswith(f, ".csv")])
    return dir, mask, methods
end

"""Fail loudly unless every dry-hour n, RMSE and bias equal the run's own pooled `no_rain` rows."""
function check_against_run(summary::DataFrame, pooled::DataFrame)
    index = Dict((row.scheme, row.product, row.method) => row for row in eachrow(pooled)
        if row.group == "rain_intensity" && row.level == "no_rain")
    checked = 0
    for row in eachrow(summary[summary.threshold .== WET_MM, :])
        reference = get(index, (row.scheme, row.anchor, row.method), nothing)
        reference === nothing && error("no pooled no_rain row for $(row.scheme)/$(row.anchor)/$(row.method)")
        row.n_dry == reference.n || error("$(row.scheme)/$(row.anchor)/$(row.method): n $(row.n_dry) vs $(reference.n)")
        isapprox(row.RMSE_dry, reference.RMSE; rtol=1e-9) ||
            error("$(row.scheme)/$(row.anchor)/$(row.method): dry RMSE $(row.RMSE_dry) vs $(reference.RMSE)")
        isapprox(row.mean_dry - row.mean_obs_dry, reference.Bias; rtol=1e-8, atol=1e-12) ||
            error("$(row.scheme)/$(row.anchor)/$(row.method): dry bias $(row.mean_dry - row.mean_obs_dry) vs $(reference.Bias)")
        checked += 1
    end
    println("  gate: $(checked) dry-hour rows match metrics_pooled.csv (n, RMSE, bias)")
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    mkpath(OUTDIR)

    println("Rebuilding the gauges on the benchmark grid ...")
    mger = study_config(OUTDIR)
    _, ids, product_data = load_global_common_product_data(mger)
    grid = product_data["GPM"]
    Y_obs, grid_times = grid.Y_obs, grid.times
    lonlat = build_X_lonlat(load_station_meta(mger.station_meta_path), ids)
    days = met_day.(grid_times)
    unique_days = sort(unique(days))
    day_number = Dict(d => k for (k, d) in enumerate(unique_days))
    column_day = [day_number[d] for d in days]
    column_season = [findfirst(==(season_of(d)), SEASONS) for d in days]
    column_hour = hour.(grid_times)
    day_season = [findfirst(==(season_of(d)), SEASONS) for d in unique_days]
    gauge_side, Y_picked, silent = gauge_strata(mger, ids, grid_times, lonlat)
    isequal(Y_picked, Y_obs) || error("the gauge record at the grid's hours differs from the benchmark's gauges")
    println("  $(length(ids)) gauges x $(length(grid_times)) hours on $(length(unique_days)) met-days; ",
        "$(count(silent)) station-hours inside silent-gauge spells")

    pooled = filter(row -> ismissing(row.fold), CSV.read(joinpath(run_dir, "metrics_pooled.csv"), DataFrame;
        types=Dict(:level => String)))
    summary, paired, strata, daily, zeroing = DataFrame[], DataFrame[], DataFrame[], DataFrame[], DataFrame[]
    for (s_index, scheme) in enumerate(SCHEMES)
        scheme_dir = joinpath(run_dir, scheme)
        fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))
        distance = nearest_train_km([fold_map[id] for id in ids], lonlat)
        distance_level = [searchsortedlast(DISTANCE_EDGES, d) for d in distance]
        for (a_index, anchor) in enumerate(ANCHORS)
            dir, mask, methods = open_anchor(scheme_dir, anchor, ids, grid_times)
            read_method(method) = read_wide_matrix(joinpath(dir, "oof_$(method).csv"), ids)
            raw = read_method("raw")
            reference = read_method(REFERENCE)
            anchor_wet = isfinite.(raw) .& (raw .>= WET_MM)
            quadrant = [(Y_obs[i, j] >= WET_MM ? 2 : 0) + (anchor_wet[i, j] ? 2 : 1) for i in axes(raw, 1), j in axes(raw, 2)]
            seed = SEED + 100_000 * s_index + 1_000 * a_index
            for (m_index, method) in enumerate(methods)
                Y = method == "raw" ? raw : method == REFERENCE ? reference : read_method(method)
                scored = mask .& isfinite.(Y)
                key = (; scheme, anchor, method)
                push!(summary, occurrence_table((i, j) -> 1, ["all"], Y_obs, Y, scored, column_day;
                    reps=BOOTSTRAP_REPS, seed=seed + m_index, key=merge(key, (; stratifier="all"))))
                if method != REFERENCE
                    both = scored .& isfinite.(reference)
                    push!(paired, paired_error_table((i, j) -> Y_obs[i, j] >= WET_MM ? 2 : 1, ["gauge dry", "gauge wet"],
                        Y_obs, Y, reference, both, column_day; reps=BOOTSTRAP_REPS, seed=seed + 100 + m_index,
                        key=merge(key, (; stratifier="gauge_state"))))
                    push!(paired, paired_error_table(quadrant, QUADRANTS, Y_obs, Y, reference, both, column_day;
                        reps=BOOTSTRAP_REPS, seed=seed + 200 + m_index, key=merge(key, (; stratifier="quadrant"))))
                    push!(paired, paired_error_table((i, j) -> Y_obs[i, j] >= WET_MM ? 3 : silent[i, j] ? 2 : 1,
                        ["gauge dry", "gauge dry, silent-gauge spell", "gauge wet"], Y_obs, Y, reference, both, column_day;
                        reps=BOOTSTRAP_REPS, seed=seed + 400 + m_index, key=merge(key, (; stratifier="screened_state"))))
                end
                if scheme == PRIMARY_SCHEME && method in PANEL
                    p_index = findfirst(==(method), PANEL)
                    labelled = vcat([
                        "anchor_state" => ((i, j) -> anchor_wet[i, j] ? 2 : 1, ["anchor dry", "anchor wet"]),
                        "distance" => ((i, j) -> distance_level[i], DISTANCE_LEVELS),
                        "season" => ((i, j) -> column_season[j], SEASONS),
                        "hour" => ((i, j) -> column_hour[j] + 1, [lpad(h, 2, '0') for h in 0:23]),
                    ], gauge_side)
                    for (k, (name, (labels, levels))) in enumerate(labelled)
                        push!(strata, occurrence_table(labels, levels, Y_obs, Y, scored, column_day;
                            reps=BOOTSTRAP_REPS, seed=seed + 300 + 20 * p_index + k,
                            key=merge(key, (; stratifier=name))))
                    end
                    push!(zeroing, zero_threshold_sensitivity(Y_obs, Y, scored; taus=ZERO_TAUS, key))
                    D_obs = daily_sums(column_day, Y_obs, scored; min_hours=MIN_DAY_HOURS)
                    D = daily_sums(column_day, Y, scored; min_hours=MIN_DAY_HOURS)
                    valid = isfinite.(D_obs) .& isfinite.(D)
                    for (t_index, threshold) in enumerate(DAILY_THRESHOLDS)
                        for (name, labels, levels) in (("all", (i, d) -> 1, ["all"]),
                                ("season", (i, d) -> day_season[d], SEASONS))
                            push!(daily, occurrence_table(labels, levels, D_obs, D, valid, collect(1:size(D, 2));
                                thresholds=[threshold], wet=threshold, reps=BOOTSTRAP_REPS,
                                seed=seed + 700 + 10 * p_index + t_index, key=merge(key, (; stratifier=name))))
                        end
                    end
                end
            end
            println("  $scheme / $anchor: $(length(methods)) methods")
            GC.gc()
        end
    end
    summary = vcat(summary...)
    check_against_run(summary, pooled)

    scope = CSV.read(joinpath(run_dir, "benchmark_scope.csv"), DataFrame)
    provenance = vcat(DataFrame(key=["run_directory"], value=[basename(run_dir)]),
        filter(:key => in(("git_commit", "git_branch", "git_dirty")), scope))
    outputs = [
        "fusion_dry_summary.csv" => summary, "fusion_dry_paired.csv" => vcat(paired...),
        "fusion_dry_strata.csv" => vcat(strata...), "fusion_dry_days.csv" => vcat(daily...),
        "fusion_zero_threshold.csv" => vcat(zeroing...), "source_provenance.csv" => provenance,
    ]
    for (name, table) in outputs
        write_csv_atomic(joinpath(OUTDIR, name), table)
    end

    primary = summary[(summary.scheme .== PRIMARY_SCHEME) .& (summary.threshold .== WET_MM) .& in.(summary.method, Ref(PANEL)), :]
    println("\nGauge-dry hours, balanced spatial CV:")
    show(select(primary, :anchor, :method, :n_dry, :mean_dry, :RMSE_dry, :POFD, :zero_share_dry, :dry_sse_share,
        :spurious_mm_per_year); allrows=true, allcols=true)
    println("\n\nWrote $(OUTDIR)")
    return OUTDIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
