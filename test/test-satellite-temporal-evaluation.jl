using Test, DataFrames, Dates, Random, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("SatelliteTemporalEvaluation")
const STE = Main.SatelliteTemporalEvaluation

"""`n` hour-ending labels starting at the first hour of the 08-08 BJT day `day`."""
hourly_axis(day::Date, n::Int) = [DateTime(day) + Hour(9 + h) for h in 0:(n - 1)]

@testset "intensity_class is lower-bound inclusive" begin
    @test STE.intensity_class.([0.0, 0.09, 0.1, 1.99, 2.0, 3.99, 4.0, 8.0, 19.99, 20.0, 150.0]) ==
        [0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 5]
    @test STE.intensity_class(NaN) == -1
    @test STE.season_of(Date(2023, 2, 28)) == "DJF"
    @test STE.season_of(Date(2023, 3, 1)) == "MAM"
    @test STE.season_of(Date(2023, 12, 1)) == "DJF"
end

@testset "evaluation_context assigns met-days, seasons and day offsets" begin
    times = hourly_axis(Date(2023, 2, 28), 48)
    ctx = STE.evaluation_context(times, ["b", "a", "b"]; regions=["a", "b"])
    @test ctx.days == [Date(2023, 2, 28), Date(2023, 3, 1)]
    @test ctx.column_day == vcat(fill(1, 24), fill(2, 24))
    @test ctx.column_season == vcat(fill(4, 24), fill(1, 24))      # DJF, then MAM
    @test ctx.column_offset[1] == 0 && ctx.column_offset[24] == 23
    @test ctx.column_hour[1] == 9
    @test ctx.station_region_index == [2, 1, 2]
    @test_throws ArgumentError STE.evaluation_context(times[[1, 3, 4]], ["a"])
    @test_throws ArgumentError STE.evaluation_context(times, ["c"]; regions=["a", "b"])
end

@testset "station_rain_events splits on the dry gap and rejects unconfirmed extents" begin
    y = [0.0, 0, 0, 1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0]
    events = STE.station_rain_events(y; min_dry_gap=3)
    @test [(e.start, e.stop) for e in events] == [(4, 6), (10, 10)]
    @test all(e -> e.complete, events)
    @test events[1].total_mm == 2.0 && events[1].peak_index == 4    # the first of two equal maxima
    @test length(STE.station_rain_events(y; min_dry_gap=1)) == 3
    y[8] = NaN
    events = STE.station_rain_events(y; min_dry_gap=3)
    @test [(e.start, e.stop) for e in events] == [(4, 6), (10, 10)]
    @test !events[1].complete && !events[2].complete
    @test !STE.station_rain_events([1.0, 0, 0, 0, 0])[1].complete    # starts at the edge of the axis
    @test isempty(STE.station_rain_events(zeros(5)))
end

@testset "score_event and best_lag: positive means the satellite is late" begin
    obs = [0.0, 0, 0, 1, 4, 1, 0, 0, 0]
    sat = [0.0, 0, 0, 0, 0, 1, 4, 1, 0]
    s = STE.score_event(obs, sat, 4, 6)
    @test s.detected
    @test s.peak_error_h == 2 && s.onset_error_h == 2 && s.end_error_h == 2
    @test s.centroid_error_h ≈ 2
    @test s.best_lag_h == 2
    @test STE.best_lag(sat, obs; max_lag=3) == -2
    @test s.volume_rel_bias == 0 && s.peak_ratio == 1
    @test s.duration_obs_h == 3 && s.duration_sat_h == 3
    drizzle = STE.score_event(obs, fill(0.3, 9), 4, 6)
    @test drizzle.detected && !drizzle.detected_res
    @test isnan(drizzle.onset_error_h_res)
    missed = STE.score_event(obs, zeros(9), 4, 6)
    @test !missed.detected && isnan(missed.peak_error_h) && isnan(missed.r_event) && isnan(missed.best_lag_h)
    @test_throws ArgumentError STE.score_event(obs, [NaN; sat[2:end]], 4, 6)
end

@testset "day_block_bootstrap is seeded and resamples whole days" begin
    day = [1, 1, 2, 2, 3]
    sums = [1.0 2.0; 1.0 2.0; 1.0 2.0; 1.0 2.0; 1.0 2.0]
    ci = STE.day_block_bootstrap(day, sums, t -> (t[2] / t[1],); reps=200, rng=MersenneTwister(1))
    @test ci.lo == [2.0] && ci.hi == [2.0]      # a ratio that is 2 on every day cannot move
    varying = [1.0 1.0; 1.0 1.0; 1.0 3.0; 1.0 3.0; 1.0 10.0]
    a = STE.day_block_bootstrap(day, varying, t -> (t[2] / t[1],); reps=200, rng=MersenneTwister(7))
    b = STE.day_block_bootstrap(day, varying, t -> (t[2] / t[1],); reps=200, rng=MersenneTwister(7))
    @test a == b
    @test a.lo[1] < mean(varying[:, 2]) < a.hi[1]
    @test STE.day_block_bootstrap(Int[], zeros(0, 2), t -> (t[1],); reps=10, rng=MersenneTwister(1)) === nothing
end

@testset "intensity_class_table scores only gauge-wet cells and agrees with its confusion counts" begin
    times = hourly_axis(Date(2023, 7, 1), 48)
    ctx = STE.evaluation_context(times, ["a", "b"])
    Y_obs = zeros(2, 48)
    Y_sat = zeros(2, 48)
    Y_obs[1, 1:6] = [0.5, 1.5, 2.0, 5.0, 10.0, 25.0]
    Y_sat[1, 1:6] = [0.0, 1.0, 3.0, 2.0, 12.0, 30.0]
    Y_obs[2, 30] = 3.0
    Y_sat[2, 30] = NaN                                  # outside the mask
    Y_sat[2, 40] = 5.0                                  # a false alarm: not in this table
    mask = STE.sample_mask(Y_obs, Y_sat)
    metrics, confusion = STE.intensity_class_table(ctx, Y_obs, Y_sat, mask; sample="s", product="p", reps=50)
    all_wet = only(metrics[(metrics.season .== "all") .& (metrics.region .== "all") .& (metrics.gauge_class .== "all_wet"), :])
    @test all_wet.n == 6
    @test all_wet.Bias ≈ mean(Y_sat[1, 1:6] .- Y_obs[1, 1:6])
    @test all_wet.POD_rain ≈ 5 / 6
    @test all_wet.class_hit ≈ 4 / 6                     # 1.5->1.0, 10->12, 25->30 and 2->3 keep their class
    @test all_wet.under_class ≈ 2 / 6 && all_wet.over_class == 0
    heavy = only(metrics[(metrics.season .== "all") .& (metrics.region .== "all") .& (metrics.gauge_class .== "heavy"), :])
    @test heavy.n == 1 && heavy.under_class == 1
    @test only(metrics[(metrics.season .== "JJA") .& (metrics.region .== "b") .& (metrics.gauge_class .== "all_wet"), :n]) == 0
    pooled = confusion[(confusion.season .== "all") .& (confusion.region .== "all"), :]
    for group in groupby(pooled, :gauge_class)
        @test sum(group.n) == only(metrics[(metrics.season .== "all") .& (metrics.region .== "all") .&
            (metrics.gauge_class .== group.gauge_class[1]), :n])
    end
    again, _ = STE.intensity_class_table(ctx, Y_obs, Y_sat, mask; sample="s", product="p", reps=50)
    @test isequal(metrics, again)
end

@testset "false_alarm_table splits gauge-dry satellite rain at the gauge resolution" begin
    times = hourly_axis(Date(2023, 7, 1), 24)
    ctx = STE.evaluation_context(times, ["a"])
    Y_obs = zeros(1, 24)
    Y_sat = zeros(1, 24)
    Y_sat[1, 1:4] = [0.3, 1.0, 5.0, 0.05]
    Y_obs[1, 10] = 2.0
    Y_sat[1, 10] = 2.0
    table = STE.false_alarm_table(ctx, Y_obs, Y_sat, STE.sample_mask(Y_obs, Y_sat); sample="s", product="p")
    row = only(table[(table.season .== "all") .& (table.region .== "all"), :])
    @test row.n_gauge_dry == 23
    @test row.false_alarm_rate ≈ 3 / 23
    @test row.false_alarm_rate_res ≈ 2 / 23
    @test row.dry_share_below_gauge_resolution ≈ 1 / 23
    @test row.dry_share_heavy ≈ 1 / 23
    @test row.FAR ≈ 3 / 4
    @test row.dry_volume_share ≈ 6.35 / 8.35
    @test row.gauge_dry_given_moderate ≈ 0.0
end

@testset "event_tables scores a sample only where all its products are finite" begin
    times = hourly_axis(Date(2023, 7, 1), 24)
    ctx = STE.evaluation_context(times, ["a"])
    Y_obs = zeros(1, 24)
    Y_obs[1, 8:10] = [1.0, 3.0, 1.0]
    Y_obs[1, 18] = 0.5
    good = copy(Y_obs)
    gappy = copy(Y_obs)
    gappy[1, 21] = NaN                                  # inside the second event's window only
    products = Dict("A" => good, "B" => gappy)
    samples = ["both" => ["A", "B"], "a_only" => ["A"]]
    tables = STE.event_tables(ctx, Y_obs, products, samples, ["s1"])
    @test nrow(tables.events) == 2
    @test tables.events.in_both == [true, false]
    @test tables.events.in_a_only == [true, true]
    @test tables.events.event_class == ["moderate", "light"]
    @test nrow(tables.scores) == 3                      # A twice, B once
    @test all(tables.scores.peak_error_h .== 0)
    # The first window stops before the second event's wet hour, whatever the pad.
    wide = STE.event_tables(ctx, Y_obs, products, samples, ["s1"]; pad=10)
    @test wide.events.window_stop[1] == times[17]
    summary = STE.event_summary(ctx, tables.events, tables.scores, samples; reps=20)
    both = only(summary[(summary.sample .== "both") .& (summary.product .== "B") .& (summary.season .== "all") .&
        (summary.region .== "all") .& (summary.event_class .== "all"), :])
    @test both.n_gauge_events == 2 && both.n_events == 1 && both.retained_share == 0.5
    @test both.POD_event == 1 && both.peak_within_1h == 1
end

@testset "diurnal cycle: harmonic phase and circular hour differences" begin
    hours = 0:22                                        # one hour missing, as FY4B's 23:00
    fit = STE.harmonic_phase(collect(hours), 2 .+ cos.(2π .* (hours .- 17) ./ 24))
    @test fit.phase_hour ≈ 17 atol=1e-9
    @test fit.amplitude ≈ 1 && fit.mean ≈ 2
    @test STE.circular_hour_difference(1, 23) == 2
    @test STE.circular_hour_difference(23, 1) == -2

    times = hourly_axis(Date(2023, 7, 1), 24 * 10)
    ctx = STE.evaluation_context(times, fill("a", 12))
    Y = [max(0.0, cos(2π * (hour(t) - 16) / 24)) * (1 + s / 12) for s in 1:12, t in times]
    shifted = [max(0.0, cos(2π * (hour(t) - 19) / 24)) * (1 + s / 12) for s in 1:12, t in times]
    mask = STE.sample_mask(Y, shifted)
    cycle = STE.diurnal_cycle_table(ctx, ["Gauge" => Y, "Late" => shifted], mask; sample="s")
    @test nrow(cycle) == 2 * 5 * 2 * 24
    summary = STE.diurnal_summary(cycle, ["Late"])
    row = only(summary[(summary.season .== "all") .& (summary.region .== "all"), :])
    @test row.phase_diff_h ≈ 3 atol=0.3
    @test row.peak_hour_diff_h == 3
    @test row.mean_amount_ratio ≈ 1 atol=1e-6
end

@testset "regional series: station coverage, block sums and correlation" begin
    times = hourly_axis(Date(2023, 7, 1), 48)
    ctx = STE.evaluation_context(times, ["a", "a", "b"]; regions=["a", "b"])
    Y = repeat(collect(1.0:48)', 3)
    mask = trues(3, 48)
    mask[1, 5] = false
    @test STE.regional_mean_series(ctx, Y, mask, 1; min_station_fraction=0.9)[5] |> isnan
    @test STE.regional_mean_series(ctx, Y, mask, 1; min_station_fraction=0.5)[5] == 5.0
    @test STE.regional_mean_series(ctx, Y, mask, 0)[6] == 6.0
    three = STE.aggregate_series(ctx, [NaN; collect(2.0:48)], 3)
    @test length(three.values) == 16
    @test isnan(three.values[1]) && three.values[2] == 4 + 5 + 6
    daily = STE.aggregate_series(ctx, collect(1.0:48), 24)
    @test daily.values == [sum(1:24), sum(25:48)] && daily.day == [1, 2]
    x = [rand(MersenneTwister(3), 40); zeros(8)]
    late = [zeros(2); x[1:end-2]]
    scores = STE.series_metrics(x, late; max_lag=4)
    @test scores.best_lag_h == 2
    @test STE.series_metrics(x[1:5], x[1:5]).n == 5 && isnan(STE.series_metrics(x[1:5], x[1:5]).r)
end
