using Test, DataFrames

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
