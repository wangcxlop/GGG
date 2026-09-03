using Test, DataFrames, Dates, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("BenchmarkDiagnostics")
const BD = Main.BenchmarkDiagnostics

"""
A fixture whose NaN pattern is known by construction.

4 stations x 6 hours, one fold per pair of stations. Every cell is evaluable (gauge and
satellite finite) except one, which is knocked out on the satellite side so the evaluable
denominator differs from the raw cell count. `mgwr` then drops:

- hour 2 entirely for fold 1 (both its stations) — a whole-hour dropout,
- hour 4 for a single station of fold 1 — a scattered dropout,
- nothing at all for fold 2.
"""
function dropout_fixture()
    ids = ["s1", "s2", "s3", "s4"]
    fold_map = Dict("s1" => 1, "s2" => 1, "s3" => 2, "s4" => 2)
    y_obs = fill(1.0, 4, 6)
    y_sat = fill(0.5, 4, 6)
    # Station 1, hour 6 is not evaluable: the satellite is missing there.
    y_sat[1, 6] = NaN
    mgwr = fill(0.7, 4, 6)
    mgwr[1:2, 2] .= NaN      # whole hour for fold 1
    mgwr[1, 4] = NaN         # scattered
    clean = fill(0.7, 4, 6)
    predictions = Dict{String,Matrix{Float64}}("mgwr" => mgwr, "clean" => clean)
    return (; ids, fold_map, y_obs, y_sat, predictions)
end

@testset "dropout_table splits whole-hour from scattered failures" begin
    fixture = dropout_fixture()
    table = BD.dropout_table(
        fixture.y_obs, fixture.y_sat, fixture.predictions, fixture.ids, fixture.fold_map;
        scheme="balanced_spatial", product="GSMaP", methods=["mgwr", "clean"],
    )
    @test nrow(table) == 4

    fold1 = only(filter(row -> row.method == "mgwr" && row.fold == 1, table))
    # Fold 1 has 2 stations x 6 hours = 12 cells, minus the one non-evaluable satellite cell.
    @test fold1.n_evaluable == 11
    @test fold1.n_dropped == 3
    @test fold1.cells_in_all_dropped_hours == 2
    @test fold1.cells_scattered == 1
    @test fold1.hours_evaluable == 6
    @test fold1.hours_all_dropped == 1
    @test fold1.own_coverage ≈ 1 - 3 / 11

    fold2 = only(filter(row -> row.method == "mgwr" && row.fold == 2, table))
    @test fold2.n_dropped == 0
    @test fold2.own_coverage == 1.0

    for row in eachrow(filter(:method => ==("clean"), table))
        @test row.n_dropped == 0
        @test row.hours_all_dropped == 0
        @test row.own_coverage == 1.0
    end
end

@testset "dropout_table counts an hour as fully dropped only over evaluable stations" begin
    fixture = dropout_fixture()
    # Knock the *evaluable* station of fold 1 out at hour 6, leaving the non-evaluable one.
    # Hour 6 then has one evaluable station and it is NaN, so the hour is fully dropped even
    # though the other station's cell is untouched.
    fixture.predictions["mgwr"][2, 6] = NaN
    table = BD.dropout_table(
        fixture.y_obs, fixture.y_sat, fixture.predictions, fixture.ids, fixture.fold_map;
        scheme="balanced_spatial", product="GSMaP", methods=["mgwr"],
    )
    fold1 = only(filter(row -> row.fold == 1, table))
    @test fold1.hours_all_dropped == 2
    @test fold1.cells_in_all_dropped_hours == 3
    @test fold1.cells_scattered == 1
end

@testset "rebuild_common_mask mirrors the benchmark's rule" begin
    y_obs = fill(1.0, 3, 4)
    y_obs[1, 1] = NaN
    a = fill(0.0, 3, 4)
    b = fill(0.0, 3, 4)
    b[2, 2] = NaN
    predictions = Dict{String,Matrix{Float64}}("a" => a, "b" => b, "c" => fill(NaN, 3, 4))
    mask = BD.rebuild_common_mask(y_obs, predictions; methods=["a", "b"])
    @test count(mask) == 10
    @test !mask[1, 1] && !mask[2, 2]
    # A method outside the defining set cannot shrink the mask, however broken it is.
    @test BD.rebuild_common_mask(y_obs, predictions; methods=["a"]) == (.!isnan.(y_obs))
    @test_throws ArgumentError BD.rebuild_common_mask(y_obs, predictions; methods=["missing"])
end

@testset "run_comparison_table pairs the two runs on cells both evaluated" begin
    y_obs = fill(4.0, 3, 4)
    # The "after" run has the wider mask: one extra cell the "before" run had excluded.
    before_mask = trues(3, 4); before_mask[1, 1] = false
    after_mask = trues(3, 4)
    # `steady` predicts identically in both runs, so its paired delta must be exactly zero even
    # though its own-mask RMSE differs — that difference is the mask moving, not the method.
    steady_before = fill(4.0, 3, 4); steady_before[1, 1] = 1.0
    steady_after = copy(steady_before)
    # `improved` halves its error between runs.
    improved_before = fill(6.0, 3, 4)
    improved_after = fill(5.0, 3, 4)
    # `patchy` is NaN at one cell in the after run; the pair must drop it on both sides.
    patchy_before = fill(4.0, 3, 4); patchy_before[3, 3] = 10.0
    patchy_after = fill(4.0, 3, 4); patchy_after[3, 3] = NaN

    before = Dict{String,Matrix{Float64}}(
        "steady" => steady_before, "improved" => improved_before,
        "patchy" => patchy_before, "gone" => fill(4.0, 3, 4),
    )
    after = Dict{String,Matrix{Float64}}(
        "steady" => steady_after, "improved" => improved_after,
        "patchy" => patchy_after, "fresh" => fill(4.0, 3, 4),
    )
    table = BD.run_comparison_table(
        y_obs, before, before_mask, after, after_mask;
        scheme="balanced_spatial", product="GPM",
    )
    @test nrow(table) == 5
    @test all(table.mask_cells_shared .== 11)

    steady = only(filter(:method => ==("steady"), table))
    @test steady.delta_paired == 0.0
    @test steady.n_paired == 11
    # The own-mask numbers do differ, which is exactly the confound the pairing removes.
    @test steady.RMSE_before != steady.RMSE_after

    improved = only(filter(:method => ==("improved"), table))
    @test improved.RMSE_paired_before ≈ 2.0
    @test improved.RMSE_paired_after ≈ 1.0
    @test improved.delta_paired ≈ -1.0
    @test improved.relative_paired ≈ -0.5

    patchy = only(filter(:method => ==("patchy"), table))
    @test patchy.n_paired == 10          # the after-run NaN drops that cell from both sides
    @test patchy.RMSE_paired_before == 0.0

    @test only(filter(:method => ==("gone"), table)).present_in == "before"
    @test only(filter(:method => ==("fresh"), table)).present_in == "after"
    @test isnan(only(filter(:method => ==("gone"), table)).delta_paired)
end

@testset "mask_cost_table reports the cells one method's failures remove" begin
    y_obs = fill(2.0, 3, 4)
    good = fill(2.0, 3, 4)
    broken = fill(2.0, 3, 4)
    broken[1, 1] = NaN
    broken[2, 3] = NaN
    predictions = Dict{String,Matrix{Float64}}("good" => good, "mgwr" => broken)
    table = BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded="mgwr", methods=["good", "mgwr"],
    )
    @test nrow(table) == 2
    @test all(table.mask_cells_full .== 10)
    @test all(table.mask_cells_reduced .== 12)
    @test all(table.cells_recovered .== 2)

    good_row = only(filter(:method => ==("good"), table))
    @test good_row.n_full == 10
    @test good_row.n_reduced == 12          # scored on every cell once mgwr stops defining the mask
    mgwr_row = only(filter(:method => ==("mgwr"), table))
    @test mgwr_row.n_reduced == 10          # its own NaNs are still skipped by the scorer

    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded="good", methods=["mgwr"],
    )
end

@testset "mask_cost_table sees correlated failure only when excluded as a group" begin
    # The failure mode on the real run: the joint models share their predictor matrices and
    # their `valid` guard, so a station-hour with a missing covariate is dropped by all three at
    # once. Excluding any one of them recovers nothing, because the other two still mask the
    # cell — which is why the one-at-a-time table reported `cells_recovered = 0` and read as
    # "the shared mask is harmless".
    y_obs = fill(2.0, 3, 4)
    good = fill(2.0, 3, 4)
    joint = ["residual_gwr", "mixed_gwr", "mgwr"]
    predictions = Dict{String,Matrix{Float64}}("good" => good)
    for name in joint
        shared_failure = fill(2.0, 3, 4)
        shared_failure[1, 1] = NaN
        shared_failure[2, 3] = NaN
        predictions[name] = shared_failure
    end
    methods = vcat("good", joint)

    for name in joint
        one_at_a_time = BD.mask_cost_table(
            y_obs, predictions; scheme="balanced_spatial", product="GPM",
            excluded=name, methods,
        )
        @test all(one_at_a_time.cells_recovered .== 0)
    end
    as_a_group = BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded=joint, methods,
    )
    @test all(as_a_group.cells_recovered .== 2)
    @test all(as_a_group.mask_cells_full .== 10)
    @test all(as_a_group.mask_cells_reduced .== 12)
    # One label per row, so the CSV keeps its shape and a single-name call is unchanged.
    @test all(as_a_group.excluded_method .== "residual_gwr,mixed_gwr,mgwr")
    @test all(BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded="mgwr", methods,
    ).excluded_method .== "mgwr")
    # `JOINT_MASK_METHODS` is the set the benchmark actually needs excluded together.
    @test BD.JOINT_MASK_METHODS == joint

    # A group must still name only mask-defining methods, name none of them twice, and leave
    # something behind to compare against.
    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded=["mgwr", "absent"], methods,
    )
    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded=["mgwr", "mgwr"], methods,
    )
    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded=methods, methods,
    )
    @test_throws ArgumentError BD.mask_cost_table(
        y_obs, predictions; scheme="balanced_spatial", product="GPM",
        excluded=String[], methods,
    )
end

@testset "satellite_quadrant_table splits the gap by gauge and satellite state" begin
    # One cell per quadrant, plus a fifth that would land in `dry_wet` but is masked out.
    y_obs = [0.0 0.0 5.0 5.0 0.0]
    y_sat = [0.0 2.0 0.0 3.0 2.0]
    adw = [0.0 0.0 4.0 4.0 0.0]
    mgwr = [1.0 2.0 5.0 5.0 9.0]
    mask = BitMatrix([true true true true false])
    predictions = Dict{String,Matrix{Float64}}("raw" => y_sat, "adw" => adw, "mgwr" => mgwr)

    table = BD.satellite_quadrant_table(
        y_obs, predictions, mask; scheme="balanced_spatial", product="GSMaP",
    )
    rows = Dict(row.quadrant => row for row in eachrow(filter(r -> r.method == "mgwr", table)))
    @test sort(collect(keys(rows))) == ["dry_dry", "dry_wet", "wet_dry", "wet_wet"]

    # The masked-out cell would have been `dry_wet`; it is neither counted nor scored.
    @test all(row.n == 1 for row in values(rows))
    @test all(row.sample_share == 0.25 for row in values(rows))

    @test rows["dry_dry"].MSE == 1.0 && rows["dry_dry"].reference_MSE == 0.0
    @test rows["dry_wet"].MSE == 4.0 && rows["dry_wet"].reference_MSE == 0.0
    @test rows["wet_dry"].MSE == 0.0 && rows["wet_dry"].reference_MSE == 1.0
    @test rows["wet_wet"].MSE == 0.0 && rows["wet_wet"].reference_MSE == 1.0

    @test rows["dry_dry"].mean_satellite == 0.0
    @test rows["dry_wet"].mean_satellite == 2.0
    @test rows["wet_wet"].mean_satellite == 3.0
    # The satellite's own error in each quadrant, i.e. how much the residual anchor is off by.
    @test rows["dry_wet"].satellite_MSE == 4.0
    @test rows["wet_dry"].satellite_MSE == 25.0

    # The whole point of the table: the contributions decompose the pooled MSE gap exactly.
    total_gap = 5.0 / 4 - 2.0 / 4
    @test sum(row.mse_gap_contribution for row in values(rows)) ≈ total_gap
    @test rows["dry_wet"].gap_share ≈ 1.0 / total_gap
end

@testset "satellite_quadrant_table counts a cell at the threshold as wet" begin
    y_obs = [0.0 0.0]
    y_sat = [0.1 0.09999]
    predictions = Dict{String,Matrix{Float64}}(
        "raw" => y_sat, "adw" => [0.0 0.0], "mgwr" => [1.0 1.0],
    )
    table = BD.satellite_quadrant_table(
        y_obs, predictions, trues(1, 2); scheme="random", product="GPM", threshold=0.1,
    )
    rows = Dict(row.quadrant => row for row in eachrow(filter(r -> r.method == "mgwr", table)))
    @test rows["dry_wet"].n == 1
    @test rows["dry_dry"].n == 1
end

@testset "read_run_grid and load_gauge_matrix align a gauge file onto a run's grid" begin
    dir = mktempdir()
    # A run artefact: no `Z` on its timestamps, two stations, three hours.
    run_path = joinpath(dir, "oof_raw.csv")
    write(run_path, """
    time,s1,s2
    2022-06-01T09:00:00,0.0,1.0
    2022-06-01T10:00:00,2.0,3.0
    2022-06-01T11:00:00,4.0,5.0
    """)
    ids, times = BD.read_run_grid(run_path)
    @test ids == ["s1", "s2"]
    @test times == [DateTime(2022, 6, 1, 9), DateTime(2022, 6, 1, 10), DateTime(2022, 6, 1, 11)]

    # A gauge file: `Z`-suffixed, out of order, wider than the run, and missing the run's last
    # hour. Column order differs from the run's, so ordering by `ids` is what is being checked.
    gauge_path = joinpath(dir, "obs.csv")
    write(gauge_path, """
    time,s2,s1,s3
    2022-06-01T10:00:00Z,30.0,20.0,0.0
    2022-06-01T08:00:00Z,99.0,99.0,0.0
    2022-06-01T09:00:00Z,10.0,,0.0
    """)
    y_obs, unmatched = BD.load_gauge_matrix(gauge_path, ids, times)
    @test size(y_obs) == (2, 3)
    @test unmatched == 1                      # 11:00 is absent from the gauge file
    @test y_obs[1, 1] === NaN                 # s1 at 09:00 is empty in the gauge file
    @test y_obs[2, 1] == 10.0
    @test y_obs[:, 2] == [20.0, 30.0]
    @test all(isnan, y_obs[:, 3])             # the unmatched hour stays NaN, never zero
end

"""
Six stations on a line, far enough apart to be resolvable and close enough that a boxcar of 500 km
covers all of them, so the weighted fit is an ordinary one over the other five.
"""
line_stations() = Float64[0.0 30.0; 0.2 30.0; 0.4 30.0; 0.6 30.0; 0.8 30.0; 1.0 30.0]

@testset "station_weight_matrix excludes each station from its own neighbourhood" begin
    lonlat = line_stations()
    weights = BD.station_weight_matrix(lonlat; kernel=4, bw=500.0, adaptive=false)
    @test size(weights) == (6, 6)
    @test all(weights[i, i] == 0.0 for i in 1:6)
    @test all(weights[i, j] == 1.0 for i in 1:6, j in 1:6 if i != j)

    # The diagonal is Inf *before* weighting, so an adaptive bandwidth cannot spend one of its
    # slots on the target itself.
    adaptive = BD.station_weight_matrix(lonlat; kernel=2, bw=3.0, adaptive=true)
    @test all(adaptive[i, i] == 0.0 for i in 1:6)
    @test all(count(>(0.0), adaptive[:, j]) <= 3 for j in 1:6)
    @test_throws ArgumentError BD.station_weight_matrix(
        lonlat; kernel=2, bw=6.0, adaptive=true,
    )
end

@testset "local_satellite_slope recovers a planted local slope" begin
    lonlat = line_stations()
    weights = BD.station_weight_matrix(lonlat; kernel=4, bw=500.0, adaptive=false)
    y_sat = Float64[0.0 1.0; 2.0 3.0; 4.0 5.0; 6.0 7.0; 8.0 9.0; 10.0 11.0]
    # r = y_obs - y_sat is exactly 3 - 0.5 * y_sat, so every neighbourhood must return -0.5.
    y_obs = (3.0 .- 0.5 .* y_sat) .+ y_sat
    fitted = BD.local_satellite_slope(y_obs, y_sat, weights)
    @test size(fitted.slope) == (6, 2)
    @test all(isapprox.(fitted.slope, -0.5; atol=1e-10))
    @test fitted.degenerate_cells == 0
    # `satellite_mean` is the neighbourhood's weighted mean, which excludes the target.
    @test fitted.satellite_mean[1, 1] ≈ mean(y_sat[2:6, 1])

    # A constant satellite field has no variance to regress on: the slope falls back to 0, the
    # forced anchor, and the cell is counted rather than silently reported as a fit.
    flat = BD.local_satellite_slope(fill(5.0, 6, 1), fill(2.0, 6, 1), weights)
    @test flat.degenerate_cells == 6
    @test all(iszero, flat.slope)

    # Blocking must not change the answer.
    blocked = BD.local_satellite_slope(y_obs, y_sat, weights; block=1)
    @test blocked.slope == fitted.slope
end

@testset "local_anchor_prediction is exact at a zero slope and floors at zero" begin
    lonlat = line_stations()
    weights = BD.station_weight_matrix(lonlat; kernel=4, bw=500.0, adaptive=false)
    y_sat = Float64[0.0 1.5; 2.0 3.0; 4.0 5.0; 6.0 7.0; 8.0 9.0; 10.0 11.0]
    y_obs = (3.0 .- 0.5 .* y_sat) .+ y_sat
    fitted = BD.local_satellite_slope(y_obs, y_sat, weights)
    prediction = Float64[0.0 1.5; 2.0 0.0; 4.0 5.0; 6.0 7.0; 0.0 9.0; 10.0 11.0]

    # The forced anchor is a zero slope; it must return the stored prediction bit for bit, which
    # is what the replay script asserts before reporting anything.
    identity = BD.local_anchor_prediction(
        y_sat, prediction, zeros(6, 2), fitted.satellite_mean; shrink=0.8,
    )
    @test identity.prediction == prediction
    @test identity.clipped_cells == 0

    # The increment is `shrink * slope * (y_sat - neighbourhood mean)`, so it vanishes exactly
    # where the target's satellite value is what its neighbours would have predicted.
    slope = fill(-0.5, 6, 2)
    replayed = BD.local_anchor_prediction(
        y_sat, prediction, slope, fitted.satellite_mean; shrink=1.0,
    )
    expected = max.(
        prediction .+ (-0.5) .* (y_sat .- fitted.satellite_mean), 0.0,
    )
    @test replayed.prediction ≈ expected

    # Flooring stops the effective anchor going below zero: the satellite may be removed, never
    # inverted, so the increment is never more negative than -y_sat.
    steep = fill(-5.0, 6, 2)
    floored = BD.local_anchor_prediction(
        y_sat, prediction, steep, fitted.satellite_mean; shrink=1.0, floor_at_zero=true,
    )
    unfloored = BD.local_anchor_prediction(
        y_sat, prediction, steep, fitted.satellite_mean; shrink=1.0,
    )
    @test all(floored.prediction .>= unfloored.prediction)
    @test all(floored.prediction .>= max.(prediction .- y_sat, 0.0) .- 1e-9)

    # NaN in either input propagates rather than being treated as zero.
    holed = copy(prediction)
    holed[1, 1] = NaN
    @test isnan(BD.local_anchor_prediction(
        y_sat, holed, slope, fitted.satellite_mean,
    ).prediction[1, 1])
end

@testset "false_alarm_coherence_table measures lift against the same-hour null" begin
    lonlat = line_stations()
    weights = BD.station_weight_matrix(lonlat; kernel=4, bw=500.0, adaptive=false)
    # Two hours, each with three gauge-dry/satellite-wet cells out of six. The boxcar covers every
    # other station, so the neighbourhood share and the same-hour share are the same quantity by
    # construction and the lift must be 1 - the null the table exists to compare against.
    y_sat = fill(0.0, 6, 2)
    y_obs = fill(0.0, 6, 2)
    y_sat[1:3, :] .= 5.0
    mask = trues(6, 2)
    fitted = BD.local_satellite_slope(y_obs, y_sat, weights)
    table = BD.false_alarm_coherence_table(
        y_obs, y_sat, mask, weights, fitted.slope, fitted.satellite_mean;
        scheme="balanced_spatial", product="fy4b", method="mgwr", kernel=4, bw=500.0,
        adaptive=false,
    )
    wet = only(filter(row -> row.quadrant == "dry_wet", table))
    @test wet.n == 6
    @test wet.neighbour_share ≈ wet.hour_share
    @test wet.lift ≈ 1.0
    dry = only(filter(row -> row.quadrant == "dry_dry", table))
    @test dry.n == 6
    @test dry.lift ≈ 1.0
    @test nrow(table) == 2                    # no wet gauge anywhere, so no wet_* quadrants
end
