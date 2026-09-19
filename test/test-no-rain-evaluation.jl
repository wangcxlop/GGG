using Test, DataFrames, Dates, Random, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("NoRainEvaluation")
const NRE = Main.NoRainEvaluation
const STE_NR = Main.SatelliteTemporalEvaluation

@testset "hours_to_rain marks missing hours as an unknown distance" begin
    y = [0.0, 0.5, 0.0, 0.0, NaN, 0.0, 1.0, 0.0]
    before, after = NRE.hours_to_rain(y)
    @test isequal(before, [Inf, 0, 1, 2, NaN, -1, 0, 1])
    @test isequal(after, [1, 0, -2, -1, NaN, 1, 0, Inf])
    B, A = NRE.hours_to_rain(permutedims(y))
    @test isequal(vec(B), before) && isequal(vec(A), after)
    # rain 2 h back, missing 1 h ahead: the distance lies in [1, 2], which spans two classes
    labels = NRE.proximity_labels(B, A; edges=[1, 3, 6, 24])
    @test vec(labels) == [1, 0, 1, 0, 0, 1, 0, 1]
    far = zeros(1, 40)
    far[1, 1] = 1.0
    @test NRE.proximity_labels(NRE.hours_to_rain(far)...)[1, 30] == 5   # > 24 h, no rain ahead
end

@testset "side_labels separates event gaps, onsets and ends" begin
    y = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    labels = vec(NRE.side_labels(NRE.hours_to_rain(permutedims(y))...; window=3))
    @test labels[2] == 1      # rain on both sides
    @test labels[4] == 3      # rain 1 h behind, 8 h ahead: after rain
    @test labels[6] == 3      # rain exactly 3 h behind is still inside the window
    @test labels[7] == 4      # 4 h behind, 5 h ahead: away from rain
    @test labels[9] == 2      # 6 h behind, 3 h ahead: before rain
    @test labels[1] == 0 && labels[12] == 0
    y2 = [0.0, 0.0, NaN, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test vec(NRE.side_labels(NRE.hours_to_rain(permutedims(y2))...; window=3))[4] == 0   # missing 1 h back
end

@testset "neighbour and network labels describe the gauges around a dry hour" begin
    lonlat = [110.0 32.0; 110.05 32.0; 110.1 32.0; 111.5 33.0]   # three ~5 km apart, one far away
    Y = [0.0 0.0 0.0;
         0.5 0.0 0.0;
         0.0 0.0 NaN;
         0.0 0.0 0.0]
    share, neighbours = NRE.neighbour_wet_share(Y, lonlat; radius_km=15.0)
    @test sort(neighbours[1]) == [2, 3] && isempty(neighbours[4])
    @test share[1, 1] ≈ 0.5 && share[3, 1] ≈ 0.5 && share[1, 2] == 0.0
    @test isnan(share[4, 1])
    @test share[1, 3] == 0.0        # station 3 missing: only station 2 reports
    labels = NRE.neighbour_labels(share, Y)
    @test labels[1, 1] == 3 && labels[1, 2] == 1 && labels[2, 1] == 0 && labels[4, 1] == 0
    network = NRE.network_wet_share(Y; min_reporting=0.9)
    @test network[1] ≈ 0.25 && network[2] == 0.0 && isnan(network[3])   # hour 3: 3 of 4 report
    nl = NRE.network_labels(network, Y)
    @test nl[1, 1] == 3 && nl[1, 2] == 1 && nl[1, 3] == 0 && nl[2, 1] == 0
    @test NRE.bin_labels([-1.0 0.0 2.9; 3.0 NaN 50.0], [-Inf, 0, 3, Inf]) == [1 2 2; 3 0 3]
    est = [0.0 0.3 0.2; 0.0 0.0 0.2; 0.0 0.0 0.2; 0.5 0.0 0.2]
    extent = NRE.phantom_extent(network, est, trues(size(est)); cutoffs=[0.0, 0.5])
    @test extent.n_hours == 1 && extent.mean_wet_share ≈ 0.25        # only hour 2 is network-dry
    @test extent.share_hours_gt_0 == 1.0 && extent.share_hours_gt_50 == 0.0
end

@testset "occurrence_table reproduces false_alarm_table and the contingency scores" begin
    times = [DateTime(2023, 7, 1, 9) + Hour(h) for h in 0:47]
    ctx = STE_NR.evaluation_context(times, ["a", "a"])
    rng = MersenneTwister(3)
    Y_obs = [rand(rng) < 0.2 ? 0.5 * rand(rng, 1:6) : 0.0 for _ in 1:2, _ in 1:48]
    Y_sat = [rand(rng) < 0.3 ? 3 * rand(rng) : 0.0 for _ in 1:2, _ in 1:48]
    Y_sat[1, 5] = NaN
    mask = STE_NR.sample_mask(Y_obs, Y_sat)
    reference = only(filter(r -> r.season == "all" && r.region == "all",
        STE_NR.false_alarm_table(ctx, Y_obs, Y_sat, mask; sample="s", product="p")))
    table = NRE.occurrence_table((i, j) -> 1, ["all"], Y_obs, Y_sat, mask, ctx.column_day; reps=200)
    at01, at05 = eachrow(table)
    @test at01.n_dry == reference.n_gauge_dry
    @test at01.POFD ≈ reference.false_alarm_rate && at05.POFD ≈ reference.false_alarm_rate_res
    @test at01.FAR ≈ reference.FAR && at05.FAR ≈ reference.FAR_res
    @test at01.dry_volume_share ≈ reference.dry_volume_share
    @test at01.mean_dry_lo <= at01.mean_dry <= at01.mean_dry_hi

    # A hand-counted 2x2: 2 hits, 1 miss, 1 false alarm, 4 correct negatives.
    obs = [1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0]
    est = [2.0 0.3 0.0 0.2 0.0 0.0 0.05 0.0]
    row = first(eachrow(NRE.occurrence_table(ones(Int, 1, 8), ["all"], obs, est, trues(1, 8), 1:8;
        thresholds=[0.1], reps=0, key=(; product="p"))))
    a, b, c, d = 2, 1, 1, 4
    @test row.product == "p" && row.n_dry == 5 && row.n_wet == 3
    @test row.POFD ≈ b / (b + d) && row.specificity ≈ d / (b + d) && row.FAR ≈ b / (a + b)
    @test row.NPV ≈ d / (c + d) && row.POD ≈ a / (a + c) && row.freq_bias ≈ (a + b) / (a + c)
    @test row.HSS ≈ 2 * (a * d - b * c) / ((a + c) * (c + d) + (a + b) * (b + d))
    @test row.mean_dry ≈ 0.25 / 5 && row.zero_share_dry ≈ 3 / 5
    @test row.mean_obs_dry == 0.0
    @test row.RMSE_dry ≈ sqrt((0.2^2 + 0.05^2) / 5)
    @test row.spurious_mm_per_year ≈ 0.25 / 8 * NRE.HOURS_PER_YEAR
    @test isnan(row.POFD_lo)
    # A level with only dry cells leaves the wet-side scores undefined.
    dry_only = first(eachrow(NRE.occurrence_table([0 0 0 1 1 1 1 1], ["dry"], obs, est, trues(1, 8), 1:8;
        thresholds=[0.1], reps=0)))
    @test dry_only.n_wet == 0 && dry_only.POFD ≈ 1 / 5 && isnan(dry_only.FAR) && isnan(dry_only.HSS)
end

@testset "paired_error_table splits the squared-error gap across levels" begin
    obs = [0.0 0.0 1.0 1.0]
    method = [0.4 0.0 1.0 0.5]
    ref = [0.2 0.0 0.0 1.0]
    table = NRE.paired_error_table([1 1 2 2], ["dry", "wet"], obs, method, ref, trues(1, 4), 1:4; reps=0)
    @test table.RMSE_method[1] ≈ sqrt(0.16 / 2) && table.RMSE_ref[1] ≈ sqrt(0.04 / 2)
    @test table.improvement[1] ≈ 1 - sqrt(0.16 / 0.04)
    @test table.sse_gap ≈ [0.12, 0.25 - 1.0]
    @test sum(table.sse_gap_share) ≈ 1.0
    @test table.mean_diff[1] ≈ 0.1
end

@testset "neighbour_gauge_baseline scores one gauge as the other's estimate" begin
    lonlat = [110.0 32.0; 110.03 32.0; 110.5 32.0]      # 1-2 ~2.8 km apart, 3 ~45 km away
    Y = [0.0 0.0 0.5 0.0;
         0.5 0.0 0.5 0.0;
         0.0 0.0 0.0 0.0]
    table = NRE.neighbour_gauge_baseline(Y, lonlat, trues(size(Y)), 1:4; edges_km=[0.0, 5.0, 30.0], reps=0)
    near = only(filter(r -> r.level == "0-5 km", table))
    # Pair (1 <- 2): station 1 dry on hours 1, 2, 4; gauge 2 wet on hour 1. Pair (2 <- 1): station 2 dry on 2, 4.
    @test near.n_pairs == 2 && near.n_dry == 5 && near.POFD ≈ 1 / 5
    @test only(filter(r -> r.level == "5-30 km", table)).n_pairs == 0
    nearest = only(filter(r -> r.level == "nearest gauge", table))
    @test nearest.n_pairs == 3 && nearest.median_km < 5
end

@testset "run_lengths censors runs that touch a gap" begin
    flags = Bool[1, 1, 0, 1, 1, 1, 0, 1, 1, 0]
    valid = Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    runs = NRE.run_lengths(flags, valid)
    @test runs.lengths == [2, 3, 2] && runs.censored == [true, false, false]
    valid[6] = false
    runs = NRE.run_lengths(flags, valid)
    @test runs.lengths == [2, 2, 2] && runs.censored == [true, true, false]
    breaks = falses(10)
    breaks[9] = true
    runs = NRE.run_lengths(flags, trues(10); breaks)
    @test runs.lengths == [2, 3, 1, 1] && runs.censored == [true, false, true, true]
end

@testset "dry_spell_table and daily_sums" begin
    Y = [1.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0;
         0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    spells = NRE.dry_spell_table(Y, trues(size(Y)); threshold=0.1, station_ids=["a", "b"], bin_edges=[1, 3])
    pooled = only(filter(r -> r.level == "all", spells.summary))
    @test pooled.n_spells == 2 && pooled.n_censored == 1 && pooled.mean_length == 2.5
    @test pooled.dry_share ≈ 13 / 16
    @test spells.histogram.n_spells == [1, 1] && spells.histogram.cells_in_spells == [2, 3]
    day = [1, 1, 1, 2, 2, 2]
    Z = [1.0 2.0 NaN 1.0 1.0 1.0]
    valid = trues(1, 6)
    @test isequal(NRE.daily_sums(day, Z, valid; min_hours=3), [NaN 3.0])
    @test NRE.daily_sums(day, Z, valid; min_hours=2) == [3.0 3.0]
end

@testset "silent_gauge_spells flags zeros its neighbours contradict" begin
    lonlat = [110.0 32.0; 110.02 32.0; 110.04 32.0; 110.06 32.0]
    D = zeros(4, 12)
    D[2:4, [2, 5, 8]] .= 10.0          # the neighbours are wet on three days
    D[1, 10] = 3.0                     # station 1 is dry for days 1-9, then records rain
    D[2, 11] = NaN
    screen = NRE.silent_gauge_spells(D, lonlat, ["a", "b", "c", "d"]; min_days=7, min_wet_days=3, k=3)
    a = only(filter(r -> r.station_id == "a", screen.spells))
    @test a.first_day == 1 && a.last_day == 9 && a.neighbour_wet_days == 3 && a.flagged
    @test all(screen.flagged[1, 1:9]) && !any(screen.flagged[1, 10:12]) && !any(screen.flagged[2:4, :])
    @test screen.p_dry_given_wet ≈ 3 / 12           # 12 station-days with a wet neighbour median; only a's are dry
    strict = NRE.silent_gauge_spells(D, lonlat, ["a", "b", "c", "d"]; min_days=7, min_wet_days=4)
    @test !any(strict.flagged)
end

@testset "long_to_wide_hourly and lagged_correlation line up a long product" begin
    ref_times = [DateTime(2022, 1, 1, 9) + Hour(h) for h in 0:3]
    utc = [DateTime(2022, 1, 1, 0), DateTime(2022, 1, 1, 1), DateTime(2022, 1, 1, 0), DateTime(2022, 1, 5)]
    wide = NRE.long_to_wide_hourly(utc, ["a", "a", "b", "a"], [1.0, 2.0, 3.0, 9.0], ref_times, ["b", "a"])
    @test wide.matched == 3
    @test isequal(wide.Y, [3.0 NaN NaN NaN; 1.0 2.0 NaN NaN])
    @test_throws ArgumentError NRE.long_to_wide_hourly(utc[[1, 1]], ["a", "a"], [1.0, 2.0], ref_times, ["a"])
    rng = MersenneTwister(1)
    A = rand(rng, 3, 50)
    B = hcat(zeros(3, 1), A[:, 1:(end - 1)])         # B lags A by one hour
    lags = NRE.lagged_correlation(A, B, trues(size(A)); lags=-2:2)
    @test lags.lag[argmax(lags.r)] == 1 && maximum(lags.r) ≈ 1.0
end

@testset "zero_threshold_sensitivity" begin
    obs = [0.0 0.0 1.0 0.5]
    est = [0.05 0.3 1.2 0.2]
    table = NRE.zero_threshold_sensitivity(obs, est, trues(1, 4); taus=[0.0, 0.1, 0.4])
    @test table.mean_dry ≈ [0.175, 0.15, 0.0]
    @test table.POFD ≈ [0.5, 0.5, 0.0] && table.POD ≈ [1.0, 1.0, 0.5]
    @test table.RMSE[1] ≈ sqrt((0.05^2 + 0.3^2 + 0.2^2 + 0.3^2) / 4)
end
