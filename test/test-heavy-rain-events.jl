using Test, DataFrames, Dates, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("HeavyRainEvents")
const HRE = Main.HeavyRainEvents

"""The hour-ending labels of `n_days` consecutive 08-08 BJT days starting at `day`."""
event_day_times(day::Date, n_days::Int=1) = [DateTime(day) + Hour(9 + h) for h in 0:(24 * n_days - 1)]

@testset "met_day assigns hour-ending labels to the 08-08 BJT day" begin
    day = Date(2023, 8, 26)
    @test HRE.met_day(DateTime(2023, 8, 26, 9)) == day
    @test HRE.met_day(DateTime(2023, 8, 27, 8)) == day
    @test HRE.met_day(DateTime(2023, 8, 26, 8)) == day - Day(1)
    @test HRE.met_day(DateTime(2023, 8, 26, 21); day_start_hour=20) == day
    @test HRE.met_day(DateTime(2023, 8, 26, 20); day_start_hour=20) == day - Day(1)
end

@testset "align_to_reference reorders stations and leaves absent hours NaN" begin
    ref_times = event_day_times(Date(2023, 8, 26))[1:3]
    times = [ref_times[3], ref_times[1]]
    Y = [1.0 2.0; 3.0 4.0; 5.0 6.0]   # stations b, a, c
    aligned = HRE.align_to_reference(ref_times, ["a", "b"], times, ["b", "a", "c"], Y)
    @test aligned[1, 1] == 4.0 && aligned[1, 3] == 3.0
    @test aligned[2, 1] == 2.0 && aligned[2, 3] == 1.0
    @test all(isnan, aligned[:, 2])
    @test_throws ArgumentError HRE.align_to_reference(ref_times, ["a", "z"], times, ["b", "a", "c"], Y)
end

@testset "daily_totals needs all 24 hours of a station-day" begin
    times = vcat(event_day_times(Date(2023, 8, 26), 2), [DateTime(2023, 8, 28, 9)])
    Y = ones(3, length(times))
    Y[2, 5] = NaN
    days, totals = HRE.daily_totals(times, Y)
    @test days == [Date(2023, 8, 26), Date(2023, 8, 27)]
    @test totals[1, :] == [24.0, 24.0]
    @test isnan(totals[2, 1]) && totals[2, 2] == 24.0
end

@testset "screen_heavy_rain_days needs more than three gauges strictly above 50 mm" begin
    days = [Date(2023, 7, 1), Date(2023, 7, 2), Date(2023, 7, 3), Date(2023, 7, 5)]
    ids = ["s$i" for i in 1:6]
    daily = zeros(6, 4)
    daily[1:3, 1] .= 60.0                        # three gauges: rejected
    daily[1:4, 2] .= 60.0; daily[5, 2] = 30.0    # four gauges: accepted
    daily[1:3, 3] .= 60.0; daily[4, 3] = 50.0    # exactly 50 mm does not count
    daily[1:4, 4] .= [51.0, 70.0, 55.0, 90.0]; daily[6, 4] = NaN
    table = HRE.screen_heavy_rain_days(days, daily, ids)
    @test table.day == [Date(2023, 7, 2), Date(2023, 7, 5)]
    @test table.n_heavy == [4, 4]
    @test table.n_valid == [6, 5]
    @test table.n_wide == [5, 4]
    @test table.frac_wide ≈ [5 / 6, 4 / 5]
    @test table.max_mm == [60.0, 90.0]
    @test table.peak_station_id[2] == "s4"
    @test nrow(HRE.screen_heavy_rain_days(days, zeros(6, 4), ids)) == 0
end

@testset "tag_rain_processes! joins consecutive days only" begin
    table = DataFrame(day=Date(2023, 7, 2) .+ Day.([0, 1, 3, 4, 5]))
    HRE.tag_rain_processes!(table)
    @test table.process_id == [1, 1, 2, 2, 2]
    @test_throws ArgumentError HRE.tag_rain_processes!(DataFrame(day=[Date(2023, 7, 2), Date(2023, 7, 1)]))
end

@testset "event_hours drops an hour for every source when one product lacks it" begin
    day = Date(2023, 8, 26)
    n = 10
    Y_obs = zeros(n, 24); Y_obs[:, 1] .= 3.0; Y_obs[:, 15] .= 1.0
    A = zeros(n, 24); A[:, 15] .= NaN               # product A misses hour 15 everywhere
    B = zeros(n, 24); B[1, 2] = NaN; B[1:2, 3] .= NaN # 90% is still available, 80% is not
    # A second day on the time axis must not leak into the first.
    times = event_day_times(day, 2)
    pad(M) = hcat(M, fill(99.0, n, 24))
    coverage = HRE.event_hours(times, pad(Y_obs), Dict("A" => pad(A), "B" => pad(B)), day)
    @test coverage.day_hour_idx == 1:24
    @test coverage.hour_idx == setdiff(1:24, [3, 15])
    @test coverage.common_hours == 22
    @test coverage.product_hours == Dict("A" => 23, "B" => 23)
    @test coverage.excluded_rain_share ≈ 10 / 40
    @test coverage.peak_excluded_rain_share ≈ 1 / 4   # tied totals: the first gauge is the peak

    # A peak gauge whose storm splits across a kept and a dropped hour, beside a gauge that is
    # larger in the excluded hour alone but incomplete, so it cannot be the peak.
    heavy_obs = copy(Y_obs); heavy_obs[1, 1] = 60.0; heavy_obs[1, 15] = 30.0
    heavy_obs[2, 15] = 100.0; heavy_obs[2, 7] = NaN
    heavy = HRE.event_hours(times, pad(heavy_obs), Dict("A" => pad(A), "B" => pad(B)), day)
    @test heavy.excluded_rain_share ≈ (30 + 100 + 8) / (60 + 30 + 100 + 3 + 8 * 4)
    @test heavy.peak_excluded_rain_share ≈ 30 / 90

    no_rain = HRE.event_hours(times, zeros(n, 48), Dict("A" => pad(A)), day)
    @test isnan(no_rain.peak_excluded_rain_share) && isnan(no_rain.excluded_rain_share)
    @test_throws ArgumentError HRE.event_hours(times[1:30], pad(Y_obs)[:, 1:30], Dict{String,Matrix{Float64}}(), day + Day(1))
end

@testset "event_station_totals scores every product on one station set" begin
    obs = [1.0 2.0 3.0; 1.0 NaN 1.0; 2.0 2.0 2.0]
    A = [0.5 0.5 0.5; 1.0 1.0 1.0; 1.0 NaN 1.0]
    B = ones(3, 3)
    totals = HRE.event_station_totals(obs, Dict("A" => A, "B" => B), [1, 3])
    @test totals.keep == [true, true, true]
    @test totals.obs == [4.0, 2.0, 4.0]
    @test totals.sat["A"] == [1.0, 2.0, 2.0]
    full = HRE.event_station_totals(obs, Dict("A" => A, "B" => B), 1:3)
    @test full.keep == [true, false, false]
    @test isnan(full.obs[2]) && isnan(full.sat["B"][3])
    @test full.sat["B"][1] == 3.0
end

@testset "select_representative_events splits widespread from localized, one day per process" begin
    table = DataFrame(
        day=Date(2023, 7, 1) .+ Day.([0, 1, 10, 20, 30, 40, 50, 60]),
        n_heavy=[100, 90, 80, 70, 5, 6, 4, 5],
        frac_wide=[0.9, 0.9, 0.8, 0.7, 0.1, 0.2, 0.2, 0.1],
        max_mm=[120.0, 130.0, 110.0, 100.0, 140.0, 150.0, 160.0, 130.0],
        process_id=[1, 1, 2, 3, 4, 5, 6, 7],
        common_hours=[23, 23, 23, 10, 22, 22, 22, 22],
        excluded_rain_share=[0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.3, 0.1],
        # Row 5 passes network-wide but its peak gauge loses half its rain.
        peak_excluded_rain_share=[0.05, 0.05, 0.05, 0.05, 0.5, 0.1, 0.1, 0.1],
    )
    selection = HRE.select_representative_events(table; n_widespread=2, n_localized=2)
    @test selection.eligible == [true, true, true, false, false, true, false, true]
    @test selection.selected == [true, false, true, false, false, true, false, true]
    @test selection.event_type[selection.selected] == ["widespread", "widespread", "localized", "localized"]
    @test !(:selected in propertynames(table))
    @test_logs (:warn, r"localized") HRE.select_representative_events(table; n_widespread=2, n_localized=3)
    @test_throws ArgumentError HRE.select_representative_events(select(table, Not(:process_id)))
end

@testset "spatial_metrics" begin
    lonlat = [110.0 32.0; 110.5 32.0; 111.0 32.5; 110.2 31.5; 110.8 33.0]
    obs = [10.0, 60.0, 30.0, 0.0, 80.0]

    perfect = HRE.spatial_metrics(obs, copy(obs), lonlat)
    @test perfect.n == 5
    @test perfect.r ≈ 1 && perfect.rho ≈ 1 && perfect.KGE ≈ 1
    @test perfect.RMSE == 0
    @test perfect.CRMSE ≈ 0 atol=1e-12
    @test perfect.sd_ratio ≈ 1 && perfect.cv_ratio ≈ 1
    @test perfect.CSI_25 == 1 && perfect.CSI_50 == 1
    @test perfect.centroid_shift_km ≈ 0 atol=1e-9
    @test perfect.peak_shift_km == 0

    shifted = HRE.spatial_metrics(obs, obs .+ 5, lonlat)
    @test shifted.r ≈ 1
    @test shifted.Bias ≈ 5
    @test shifted.CRMSE ≈ 0 atol=1e-9
    @test shifted.RB_pct ≈ 100 * 5 / 36

    # Taylor-diagram identity, which is why sd_ratio uses the population std.
    noisy = [20.0, 40.0, 35.0, 5.0, 60.0]
    m = HRE.spatial_metrics(obs, noisy, lonlat)
    sd_obs = std(obs; corrected=false)
    @test m.CRMSE^2 ≈ sd_obs^2 * (1 + m.sd_ratio^2 - 2 * m.sd_ratio * m.r)

    moved = HRE.spatial_metrics(obs, [10.0, 90.0, 30.0, 0.0, 80.0], lonlat)
    @test moved.peak_shift_km ≈ Main.TraditionalInterpolation.haversine_distance_matrix(lonlat[5:5, :], lonlat[2:2, :])[1, 1]

    @test HRE._tied_ranks([2.0, 1.0, 2.0, 3.0]) == [2.5, 1.0, 2.5, 4.0]
    tied = HRE.spatial_metrics([1.0, 2.0, 2.0, 3.0], [1.0, 1.0, 2.0, 3.0], lonlat[1:4, :])
    @test tied.rho ≈ 5 / 6

    dropped = HRE.spatial_metrics([obs; NaN], [obs; 1.0], vcat(lonlat, [110.0 33.0]))
    @test dropped.n == 5
    @test haskey(HRE.spatial_metrics(obs, obs, lonlat; thresholds=(2.5,)), :CSI_2p5)
    @test_throws ArgumentError HRE.spatial_metrics([1.0, NaN], [NaN, 1.0], lonlat[1:2, :])
end

@testset "idw_surface interpolates finite gauges onto the grid" begin
    lonlat = [110.0 32.0; 111.0 33.0; 110.5 32.5]
    surface = HRE.idw_surface(lonlat[1:2, :], [10.0, 20.0]; bounds=(110.0, 111.0, 32.0, 33.0), step_deg=0.5)
    @test nrow(surface) == 9
    @test only(surface.value[(surface.lon .== 110.0) .& (surface.lat .== 32.0)]) ≈ 10.0
    @test all(10.0 .<= surface.value .<= 20.0)
    with_gap = HRE.idw_surface(lonlat, [10.0, 20.0, NaN]; bounds=(110.0, 111.0, 32.0, 33.0), step_deg=0.5)
    @test with_gap.value ≈ surface.value
end

@testset "phase_spatial_metrics: phases with opposite gradients cancel in the day total" begin
    n = 60
    lat = collect(range(31.0, 33.0; length=n))
    lonlat = hcat(fill(110.0, n), lat)
    noise1, noise2 = 2.0 .* sin.(2.3 .* (1:n)), 2.0 .* cos.(1.7 .* (1:n))   # small-scale structure
    north, south = 10 .* (lat .- 31), 10 .* (33 .- lat)
    Y_obs = hcat(north .+ noise1, south .+ noise2)            # day total = 20 + noise: no gradient
    sat = hcat(2.0 .* north, 0.5 .* south)                     # right pattern, wrong amplitude per phase
    table = HRE.phase_spatial_metrics(Y_obs, Dict("SAT" => sat), lonlat, [1, 2]; n_phases=2)
    @test table.phase == ["1", "2", "all"]
    @test all(table.r[1:2] .> 0.9)
    @test abs(table.r[3]) < 0.3
    @test isapprox(table.obs_share[1:2], [0.5, 0.5]; atol=0.02) && table.obs_share[3] == 1.0
    @test table.obs_r_lat[1] > 0.9 && table.obs_r_lat[2] < -0.9
    @test table.est_r_lat[3] ≈ 1.0                             # the product keeps a gradient the gauges lack
    @test_throws ArgumentError HRE.phase_spatial_metrics(Y_obs, Dict("SAT" => sat), lonlat, [1, 2]; n_phases=3)
end

@testset "small_scale_variance_share separates a smooth field from a noisy one" begin
    lonlat = hcat(fill(110.0, 40), collect(range(31.0, 33.0; length=40)))   # ~5.7 km spacing
    smooth = HRE.small_scale_variance_share(lonlat, lonlat[:, 2]; max_km=10.0)
    noisy = HRE.small_scale_variance_share(lonlat, [isodd(i) ? 1.0 : 0.0 for i in 1:40]; max_km=10.0)
    @test smooth < 0.01
    @test noisy > 1.0
    @test isnan(HRE.small_scale_variance_share(lonlat, lonlat[:, 2]; max_km=1.0))
end
