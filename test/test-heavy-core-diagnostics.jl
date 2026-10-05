using Test, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("HeavyCoreDiagnostics")
const HCD = Main.HeavyCoreDiagnostics

"""A 0.1 deg grid around (110, 32) flattened to `(lonlat, values)`, with a Gaussian core at `peak`."""
function synthetic_field(peak; amplitude=100.0, sigma_deg=0.03)
    lons, lats = 109.5:0.1:110.5, 31.5:0.1:32.5
    lonlat = [p[k] for p in [(x, y) for x in lons for y in lats], k in 1:2]
    values = [amplitude * exp(-((x - peak[1])^2 + (y - peak[2])^2) / (2sigma_deg^2)) for (x, y) in eachrow(lonlat)]
    return lonlat, values
end

@testset "neighbourhood_max finds a displaced core and its offset" begin
    lonlat, values = synthetic_field((110.2, 32.0))
    near, far = HCD.neighbourhood_max(lonlat, values, (110.0, 32.0); radii_km=(10.0, 25.0))
    @test far.max ≈ 100.0
    @test (far.max_lon, far.max_lat) == (110.2, 32.0)
    @test isapprox(far.offset_km, 18.9; atol=0.2)   # 0.2 deg of longitude at 32 N
    @test near.max < 0.1 * far.max                   # the core lies outside the 10 km disc
    @test HCD.neighbourhood_max(lonlat, values, (0.0, 0.0); radii_km=(10.0,))[1].n == 0
end

@testset "upscaled_gauge_truth averages the gauges inside each radius" begin
    lonlat = [110.0 32.0; 110.03 32.0; 110.2 32.0; 110.0 32.03]
    totals = [100.0, 60.0, 10.0, NaN]
    within5, within25 = HCD.upscaled_gauge_truth(lonlat, totals, 1; radii_km=(5.0, 25.0))
    @test (within5.n, within5.mean) == (2, 80.0)   # the NaN gauge is skipped
    @test within25.n == 3 && within25.mean ≈ 170 / 3
end

@testset "deficit_attribution splits the deficit into parts that sum to it" begin
    # A uniform field: the areal truth equals the point value, so representativeness is zero.
    uniform = HCD.deficit_attribution(80.0, 20.0, 80.0, 20.0)
    @test uniform.representativeness_share == 0.0
    @test uniform.magnitude_share == 1.0
    parts = HCD.deficit_attribution(120.0, 20.0, 80.0, 50.0)
    @test parts.representativeness + parts.displacement + parts.magnitude ≈ parts.deficit
    @test (parts.representativeness, parts.displacement, parts.magnitude) == (40.0, 30.0, 30.0)
    # A nearby maximum above the areal truth cannot explain more than the pixel-scale deficit.
    @test HCD.deficit_attribution(120.0, 20.0, 80.0, 200.0).displacement == 60.0
    @test isnan(HCD.deficit_attribution(20.0, 30.0, 25.0, 30.0).magnitude_share)
end

@testset "core_hourly_diagnostics separates a missed burst from a timing shift" begin
    obs = [0, 0, 0, 1, 2, 40, 2, 1, 0, 0, 0, 0.0]          # pad 3 either side of a 6 h event
    damped = [0, 0, 0, 1, 2, 8, 2, 1, 0, 0, 0, 0.0]
    burst = HCD.core_hourly_diagnostics(obs, damped; pad=3)
    @test burst.obs_total == 46.0 && burst.sat_total == 14.0
    @test burst.peak_ratio == 0.2 && burst.peak_lag_h == 0
    @test burst.top_deficit_share == 1.0                     # all of the deficit is in the burst
    late = HCD.core_hourly_diagnostics(obs, circshift(obs, 2); pad=3)
    @test late.best_lag_h == 2 && late.best_lag_r ≈ 1.0
    @test late.peak_lag_h == 2
    @test late.sat_total_padded == 46.0                      # the shifted rain is still in the padded window
end

@testset "neighbour_consistency compares a gauge with its nearest neighbours" begin
    lonlat = [110.0 32.0; 110.02 32.0; 110.0 32.02; 110.5 32.5]
    result = HCD.neighbour_consistency(lonlat, [120.0, 100.0, 80.0, 0.0], 1; k=2)
    @test result.k == 2 && result.neighbour_mean == 90.0 && result.ratio ≈ 4 / 3
end
