#!/usr/bin/env julia

using Test
using CSV, DataFrames, Dates
using Random

const BENCHMARK_TEST_ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
include(joinpath(BENCHMARK_TEST_ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")

function _write_benchmark_wide(path, times, ids, values)
    df = DataFrame(time=times)
    for (index, id) in enumerate(ids)
        df[!, Symbol(id)] = values[index, :]
    end
    CSV.write(path, df)
end

@testset "Traditional interpolation" begin
    lonlat = [
        110.0 30.0
        111.0 30.0
        110.0 31.0
        111.0 31.0
        110.5 30.5
    ]
    constant_values = fill(3.5, 5, 3)
    @test idw_predict(lonlat, constant_values, lonlat[1:2, :]) ≈ fill(3.5, 2, 3)
    @test adw_predict(lonlat, constant_values, lonlat[1:2, :]) ≈ fill(3.5, 2, 3)

    exact_values = reshape(collect(1.0:5.0), 5, 1)
    @test idw_predict(lonlat, exact_values, lonlat[3:3, :])[1] == 3.0
    @test adw_predict(lonlat, exact_values, lonlat[3:3, :])[1] == 3.0
    @test isfinite(idw_predict(lonlat, exact_values, lonlat; exclude_self=true)[1])
    asymmetric_target = reshape([110.2, 30.4], 1, 2)
    @test adw_predict(lonlat, exact_values, asymmetric_target)[1] !=
        idw_predict(lonlat, exact_values, asymmetric_target)[1]
    @test isfinite(idw_predict(lonlat, exact_values, asymmetric_target; neighbors=99)[1])

    center = (mean(lonlat[:, 1]), mean(lonlat[:, 2]))
    xy = local_km_coordinates(lonlat; center=center)
    plane = reshape(2.0 .+ 0.03 .* xy[:, 1] .- 0.02 .* xy[:, 2], :, 1)
    tps_fit = tps_predict(lonlat, plane, lonlat; smooth=0.01)
    @test maximum(abs.(tps_fit .- plane)) < 1e-7
    @test all(isfinite, tps_loo_predict(lonlat, plane; smooth=0.01))

    missing_values = copy(constant_values)
    missing_values[1:3, 2] .= NaN
    missing_prediction = tps_predict(lonlat, missing_values, lonlat[1:1, :]; smooth=0.01)
    @test isfinite(missing_prediction[1, 1])
    @test isnan(missing_prediction[1, 2])
    collinear = [110.0 30.0; 110.1 30.0; 110.2 30.0; 110.3 30.0]
    @test all(isnan, tps_predict(collinear, fill(1.0, 4, 1), collinear; smooth=0.01))

    selected_idw = (; power=2.0, neighbors=0)
    target_lonlat = reshape([110.25, 30.25], 1, 2)
    direct_a = predict_selected(
        selected_idw, "idw", "direct", lonlat, target_lonlat,
        constant_values, zeros(size(constant_values)), zeros(1, 3),
    )
    direct_b = predict_selected(
        selected_idw, "idw", "direct", lonlat, target_lonlat,
        constant_values, fill(100.0, size(constant_values)), fill(200.0, 1, 3),
    )
    @test direct_a == direct_b
    @test_throws ArgumentError predict_selected(
        selected_idw, "idw", "residual", lonlat, target_lonlat,
        constant_values, zeros(size(constant_values)), zeros(1, 3),
    )
end

@testset "mixed_gwr and mgwr sweep every configured kernel" begin
    # `gwr`/`residual_gwr` have always swept every configured kernel; `mixed_gwr`/`mgwr` used to
    # be hardcoded to BISQUARE regardless of `cfg.mger.kernels`. This is the parity guard: with
    # more than one kernel configured, both methods' scans must actually vary kernel, not just
    # report it.
    rng = MersenneTwister(2027)
    n = 30
    train_lonlat = hcat(110.0 .+ rand(rng, n), 30.0 .+ rand(rng, n))
    n_time = 6
    y_sat = 5.0 .+ 0.3 .* randn(rng, n, n_time)
    y_obs = y_sat .+ 0.2 .* randn(rng, n, n_time)

    mger = MGERConfig(
        station_meta_path="unused.csv", obs_hourly_wide_path="unused.csv",
        sat_paths=Dict("FY4B" => "unused.csv"), outdir="unused",
        kernels=[GAUSSIAN, BISQUARE], bw_adaptive=[8.0, 16.0], bw_fixed_km=[20.0, 50.0],
    )
    cfg = InterpolationBenchmarkConfig(mger=mger, mgwr_max_tuning_iterations=3)

    mixed_rows = NamedTuple[]
    mixed_selected = select_interpolation_parameter!(
        mixed_rows, cfg, "mixed_gwr", "residual", :balanced_spatial, "TEST", 1,
        train_lonlat, y_obs, y_sat,
    )
    @test length(unique(row.kernel for row in mixed_rows)) > 1
    @test length(unique(row.adaptive for row in mixed_rows)) > 1
    @test mixed_selected.kernel in cfg.mger.kernels
    @test mixed_selected.adaptive isa Bool
    @test count(row.selected for row in mixed_rows) == 1

    mgwr_rows = NamedTuple[]
    mgwr_selected = select_interpolation_parameter!(
        mgwr_rows, cfg, "mgwr", "residual", :balanced_spatial, "TEST", 1,
        train_lonlat, y_obs, y_sat,
    )
    @test length(unique(row.kernel for row in mgwr_rows)) > 1
    @test length(unique(row.adaptive for row in mgwr_rows)) > 1
    @test mgwr_selected.kernel in cfg.mger.kernels
    @test mgwr_selected.adaptive isa Bool
    @test count(row.selected for row in mgwr_rows) == length(mgwr_selected.bandwidths)
end

@testset "Balanced spatial folds" begin
    ids = string.(1:23)
    lonlat = hcat(
        [110.0 + 0.1 * mod(i - 1, 6) for i in 1:23],
        [30.0 + 0.1 * div(i - 1, 6) for i in 1:23],
    )
    folds1 = split_stations_balanced_spatial_kfold(ids, lonlat; k=5, seed=42)
    folds2 = split_stations_balanced_spatial_kfold(ids, lonlat; k=5, seed=42)
    @test folds1 == folds2
    @test sort(reduce(vcat, folds1)) == sort(ids)
    @test maximum(length.(folds1)) - minimum(length.(folds1)) <= 1

    # Both initialisations must stay valid partitions with balanced fold sizes.
    for init in (:kmeanspp, :farthest)
        folds = split_stations_balanced_spatial_kfold(ids, lonlat; k=5, seed=7, center_init=init)
        @test sort(reduce(vcat, folds)) == sort(ids)
        @test maximum(length.(folds)) - minimum(length.(folds)) <= 1
    end
    @test split_stations_balanced_spatial_kfold(ids, lonlat; k=5, seed=7, center_init=:farthest) ==
        split_stations_balanced_spatial_kfold(ids, lonlat; k=5, seed=7, center_init=:farthest)
end

@testset "Different seeds give different spatial partitions" begin
    # Repeated cross-validation only measures split sensitivity if seeds actually move stations
    # between folds. The farthest-point initialisation randomises just the first center, which on
    # the real 237-station network collapsed 20 seeds into 10 partitions differing by 1-5 stations.
    rng = MersenneTwister(2026)
    ids = string.(1:120)
    lonlat = hcat(109.5 .+ 2.5 .* rand(rng, 120), 29.5 .+ 2.5 .* rand(rng, 120))
    signature(folds) = sort([sort(f) for f in folds])
    seeds = [20260627 + 1000 * i for i in 0:9]
    partitions = [split_stations_balanced_spatial_kfold(
        ids, lonlat; k=5, seed=s, center_init=:kmeanspp) for s in seeds]
    @test length(unique(signature.(partitions))) >= 8
    # Every partition stays a valid, capacity-balanced split.
    @test all(p -> sort(reduce(vcat, p)) == sort(ids), partitions)
    @test all(p -> maximum(length.(p)) - minimum(length.(p)) <= 1, partitions)
end

@testset "Hilbert curve index" begin
    # The whole seed-free scheme rests on this ordering being a genuine Hilbert curve: a bijection
    # onto 0:n^2-1 whose consecutive cells are grid neighbours. If it degenerated into a row-major
    # or Morton order the fold blocks would stop being compact and nothing downstream would say so.
    for order in 2:5
        n = 1 << order
        cells = [(x, y) for x in 0:(n - 1), y in 0:(n - 1)]
        located = Dict(_hilbert_index(x, y, order) => (x, y) for (x, y) in cells)
        @test length(located) == n * n
        @test extrema(keys(located)) == (0, n * n - 1)
        @test all(0:(n * n - 2)) do d
            (x1, y1), (x2, y2) = located[d], located[d + 1]
            abs(x1 - x2) + abs(y1 - y2) == 1
        end
    end
    # A degenerate axis must not divide by zero.
    @test _hilbert_axis([3.0, 3.0, 3.0]) == [0, 0, 0]
end

@testset "Hilbert folds are seed-free and rotation-indexed" begin
    rng = MersenneTwister(2026)
    ids = string.(1:120)
    lonlat = hcat(109.5 .+ 2.5 .* rand(rng, 120), 29.5 .+ 2.5 .* rand(rng, 120))
    signature(folds) = sort([sort(f) for f in folds])

    # The point of the scheme: the seed argument is inert, so no reported partition is an RNG draw.
    @test split_stations_balanced_spatial_kfold(
            ids, lonlat; k=5, seed=1, center_init=:hilbert, rotation=3) ==
        split_stations_balanced_spatial_kfold(
            ids, lonlat; k=5, seed=999_999, center_init=:hilbert, rotation=3)
    # Rotation, not the seed, is what moves stations between folds.
    @test signature(split_stations_balanced_spatial_kfold(
            ids, lonlat; k=5, center_init=:hilbert, rotation=0)) !=
        signature(split_stations_balanced_spatial_kfold(
            ids, lonlat; k=5, center_init=:hilbert, rotation=1))

    rotations = [split_stations_balanced_spatial_kfold(
        ids, lonlat; k=5, center_init=:hilbert, rotation=r) for r in 0:19]
    # Repeated cross-validation still measures split sensitivity: the rotations must not collapse
    # the way `:farthest` did. Capacity-constrained Lloyd pulls different initialisations toward the
    # same optima, so the two schemes are close rather than identical. Measured distinct partitions,
    # hilbert rotations against `:kmeanspp` seeds:
    #
    #   this uniform 120-station cloud   7/10 and 16/20, against 9/10 and 13/20
    #   the real clustered 237-station network   8/10 and 14/20, against 8/10 and 15/20
    #
    # Comparable in both directions, which is the bar: the seed goes away without costing the
    # partition diversity repeated CV needs. `scripts/verify_fold_rotations.jl` regates this on the
    # real network, where the comparison that matters is made against the incumbent directly.
    @test length(unique(signature.(rotations[1:10]))) >= 7
    @test length(unique(signature.(rotations))) >= 15
    @test all(p -> sort(reduce(vcat, p)) == sort(ids), rotations)
    @test all(p -> maximum(length.(p)) - minimum(length.(p)) <= 1, rotations)
    # Every rotation is reproducible on its own.
    @test all(r -> split_stations_balanced_spatial_kfold(
        ids, lonlat; k=5, center_init=:hilbert, rotation=r) == rotations[r + 1], 0:19)
end

"""Path-only MGERConfig for tests that exercise configuration validation, not data loading."""
_dummy_mger() = MGERConfig(
    station_meta_path="unused.csv", obs_hourly_wide_path="unused.csv",
    sat_paths=Dict("FY4B" => "unused.csv"), outdir="unused",
)

@testset "Repeated cross-validation seeds" begin
    base = InterpolationBenchmarkConfig(mger=_dummy_mger())
    @test benchmark_seeds(base) == [base.seed]
    repeated = InterpolationBenchmarkConfig(
        mger=_dummy_mger(), seeds=[11, 22, 33],
    )
    @test benchmark_seeds(repeated) == [11, 22, 33]
    @test_throws ArgumentError _validate_benchmark_config(
        InterpolationBenchmarkConfig(mger=_dummy_mger(), seeds=[5, 5]), 50,
    )
    @test_throws ArgumentError _validate_benchmark_config(
        InterpolationBenchmarkConfig(
            mger=_dummy_mger(), fold_center_init=:nope,
        ), 50,
    )
    # The seed-free initialisation is the default and must validate; the seeded ones stay reachable.
    for init in (:hilbert, :kmeanspp, :farthest)
        @test _validate_benchmark_config(
            InterpolationBenchmarkConfig(mger=_dummy_mger(), fold_center_init=init), 50,
        ) === nothing
    end
    @test InterpolationBenchmarkConfig(mger=_dummy_mger()).fold_center_init === :hilbert
    # A single repeat keeps the historical output layout; repeats nest under repeat_<i>.
    @test _repeat_dir("/out", 1, 1) == "/out"
    @test _repeat_dir("/out", 2, 5) == joinpath("/out", "repeat_02")
end

@testset "Repeat summaries report across-partition spread" begin
    rows = NamedTuple[]
    for (repeat_value, seed_value) in [(1, 11), (2, 22)], method in ["adw", "residual_gwr"]
        # adw wins the first partition, residual_gwr the second: an unstable ranking.
        rmse = method == "adw" ? (repeat_value == 1 ? 1.0 : 2.0) : (repeat_value == 1 ? 2.0 : 1.0)
        push!(rows, _empty_metric_row(;
            scheme="balanced_spatial", product="FY4B", method, fold=missing,
            repeat=repeat_value, seed=seed_value, group="overall", level="all",
            n=100, coverage=1.0, RMSE=rmse, MAE=rmse,
        ))
    end
    metrics = DataFrame(rows)

    summary = summarize_repeats(metrics)
    @test nrow(summary) == 2
    @test all(summary.n_repeats .== 2)
    adw = only(filter(:method => ==("adw"), summary))
    @test adw.RMSE_mean ≈ 1.5
    @test adw.RMSE_min ≈ 1.0
    @test adw.RMSE_max ≈ 2.0
    @test adw.RMSE_std > 0

    stability = method_rank_stability(metrics)
    @test nrow(stability) == 2
    @test all(stability.n_repeats .== 2)
    # Each method wins exactly once, so neither ranking survives the second partition.
    @test all(stability.win_count .== 1)
    @test all(stability.mean_rank .== 1.5)
end

"""One pooled `metrics` row for `assess_gwr_claim` test fixtures."""
_claim_metric_row(; product, method, group, level, RMSE, CSI=0.5, FAR=0.1, repeat=1) =
    _empty_metric_row(;
        scheme="balanced_spatial", product, method, fold=missing, repeat, seed=10 + repeat,
        group, level, n=100, coverage=1.0, RMSE, MAE=RMSE, CSI, FAR,
    )

"""One `paired_bootstrap_rows`-shaped row for `assess_gwr_claim` test fixtures."""
_claim_bootstrap_row(; product, baseline, stratum, ci_low, pvalue_holm, repeat=1) = (;
    scheme="balanced_spatial", product, stratum, baseline, repeat, seed=10 + repeat,
    n=100, n_day=10, reps=200, RMSE_baseline=1.0, RMSE_gwr=0.9, delta_RMSE=0.1,
    relative_improvement=0.1, ci_low, ci_high=0.2, pvalue=0.01, pvalue_holm,
)

@testset "assess_gwr_claim requires beating every baseline, not just one" begin
    products = ["TEST"]
    rows = NamedTuple[]
    # overall: adw has the lowest RMSE (the historical argmin pick); GWR beats all three.
    push!(rows, _claim_metric_row(product="TEST", method="residual_gwr", group="overall", level="all", RMSE=0.9))
    push!(rows, _claim_metric_row(product="TEST", method="idw", group="overall", level="all", RMSE=1.5))
    push!(rows, _claim_metric_row(product="TEST", method="adw", group="overall", level="all", RMSE=1.0))
    push!(rows, _claim_metric_row(product="TEST", method="tps", group="overall", level="all", RMSE=1.3))
    # heavy, moderate, year, event: GWR comfortably beats/matches all three baselines.
    for (method, RMSE) in (("residual_gwr", 0.8), ("idw", 1.5), ("adw", 1.2), ("tps", 1.4))
        push!(rows, _claim_metric_row(product="TEST", method=method, group="rain_intensity", level="heavy", RMSE=RMSE))
    end
    for (method, RMSE) in (("residual_gwr", 1.0), ("idw", 1.02), ("adw", 1.01), ("tps", 1.03))
        push!(rows, _claim_metric_row(product="TEST", method=method, group="rain_intensity", level="moderate", RMSE=RMSE))
    end
    for (method, RMSE) in (("residual_gwr", 0.85), ("idw", 1.4), ("adw", 1.1), ("tps", 1.3))
        push!(rows, _claim_metric_row(product="TEST", method=method, group="year", level="2022", RMSE=RMSE))
    end
    for (method, CSI, FAR) in (("residual_gwr", 0.6, 0.10), ("idw", 0.55, 0.12), ("adw", 0.58, 0.11), ("tps", 0.56, 0.13))
        push!(rows, _claim_metric_row(product="TEST", method=method, group="event_threshold", level="0.1", RMSE=1.0, CSI=CSI, FAR=FAR))
    end
    metrics = DataFrame(rows)

    # idw and adw (the argmin/toughest baseline by RMSE) are significant; tps is not — a
    # baseline the old argmin selection never even looked at, since it only tested adw.
    bootstrap_mixed = DataFrame([
        _claim_bootstrap_row(product="TEST", baseline="idw", stratum="overall", ci_low=0.3, pvalue_holm=0.01),
        _claim_bootstrap_row(product="TEST", baseline="adw", stratum="overall", ci_low=0.05, pvalue_holm=0.02),
        _claim_bootstrap_row(product="TEST", baseline="tps", stratum="overall", ci_low=-0.05, pvalue_holm=0.2),
    ])
    claim_mixed = assess_gwr_claim(metrics, bootstrap_mixed, products)
    row = only(eachrow(claim_mixed))
    @test row.best_traditional == "adw" # still the argmin baseline, but now descriptive-only
    @test row.overall_win_count == 2
    # The old logic tested only the argmin (adw) baseline, which IS significant here, and would
    # have reported this product as significant. Requiring all three catches the tps failure.
    @test row.paired_significant == false
    @test row.product_supported == false

    # All three baselines significant: the claim can now be supported.
    bootstrap_all = DataFrame([
        _claim_bootstrap_row(product="TEST", baseline="idw", stratum="overall", ci_low=0.3, pvalue_holm=0.01),
        _claim_bootstrap_row(product="TEST", baseline="adw", stratum="overall", ci_low=0.05, pvalue_holm=0.02),
        _claim_bootstrap_row(product="TEST", baseline="tps", stratum="overall", ci_low=0.1, pvalue_holm=0.03),
    ])
    claim_all = assess_gwr_claim(metrics, bootstrap_all, products)
    row_all = only(eachrow(claim_all))
    @test row_all.overall_win_count == 3
    @test row_all.paired_significant == true
    @test row_all.heavy_win_count == 3
    @test row_all.moderate_win_count == 3
    @test row_all.year_win_count == 1
    @test row_all.event_win_count == 3
    @test row_all.product_supported == true

    # A missing baseline comparison (e.g. its bootstrap deltas were empty) must fail closed,
    # not be silently treated as passing.
    bootstrap_incomplete = DataFrame([
        _claim_bootstrap_row(product="TEST", baseline="idw", stratum="overall", ci_low=0.3, pvalue_holm=0.01),
        _claim_bootstrap_row(product="TEST", baseline="adw", stratum="overall", ci_low=0.05, pvalue_holm=0.02),
    ])
    claim_incomplete = assess_gwr_claim(metrics, bootstrap_incomplete, products)
    @test only(eachrow(claim_incomplete)).paired_significant == false

    # Partial win in a point-estimate stratum: GWR beats two of three heavy-rain baselines by
    # less than the 5% threshold for the third.
    partial_rows = copy(rows)
    partial_idx = findfirst(r -> r.method == "tps" && r.group == "rain_intensity" && r.level == "heavy", partial_rows)
    partial_rows[partial_idx] = _claim_metric_row(
        product="TEST", method="tps", group="rain_intensity", level="heavy", RMSE=0.82,
    )
    claim_partial_heavy = assess_gwr_claim(DataFrame(partial_rows), bootstrap_all, products)
    row_partial = only(eachrow(claim_partial_heavy))
    @test row_partial.heavy_win_count == 2
    @test row_partial.product_supported == false
end

@testset "DEM screening uses training stations only" begin
    rng = MersenneTwister(91)
    n = 24
    train_idx = 1:20
    validation_idx = 21:24
    ids = string.(3001:(3000 + n))
    lonlat = hcat(
        collect(range(109.5, 112.0; length=n)),
        30.0 .+ rand(rng, n),
    )
    elevation = collect(range(100.0, 1800.0; length=n))
    aspect = rand(rng, n) .* 360
    terrain = DataFrame(
        station_id=ids, elevation_m=elevation, slope_deg=5.0 .+ rand(rng, n) .* 20,
        aspect_deg=aspect, aspect_sin=sind.(aspect), aspect_cos=cosd.(aspect),
    )
    times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 11))
    y_obs = fill(5.0, n, length(times))
    residual = 0.002 .* elevation
    y_sat = y_obs .- reshape(residual, :, 1)
    dem = DEMTerrainExperiment.DEMExperimentConfig(
        outdir="unused", min_wet_hours=1, bandwidth_candidates=[8, 12],
        screen_permutations=19, spatial_permutations=19, q_threshold=0.2,
        max_iterations=50,
    )
    first = screen_dem_subset(
        terrain[train_idx, :], lonlat[train_idx, :], times,
        y_obs[train_idx, :], y_sat[train_idx, :], dem;
        scheme="test", product="FY4B", fold=1, phase="cv", seed=17,
    )
    changed_obs = copy(y_obs)
    changed_sat = copy(y_sat)
    changed_obs[validation_idx, :] .= 10_000
    changed_sat[validation_idx, :] .= -10_000
    second = screen_dem_subset(
        terrain[train_idx, :], lonlat[train_idx, :], times,
        changed_obs[train_idx, :], changed_sat[train_idx, :], dem;
        scheme="test", product="FY4B", fold=1, phase="cv", seed=17,
    )
    @test isequal(first.screen, second.screen)
    @test isequal(first.roles, second.roles)
    @test first.response == second.response
end

"""Small in-memory joint-covariate fixture: population-wide matrices/tables plus a `train_idx`
subset, mirroring what `load_joint_benchmark_inputs` would hand the fold loop, without going
through file I/O (`ERA5VariableSelection.load_era5_panel` enforces full-year file completeness,
which is disproportionate to fixture here)."""
function synthetic_joint_benchmark_fixture()
    rng = MersenneTwister(2026)
    n, nt = 30, 16
    ids = string.(4001:(4000 + n))
    times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, nt - 1))
    lonlat = hcat(collect(range(109.5, 112.0; length=n)), 30.0 .+ rand(rng, n))
    terrain = DataFrame(
        station_id=ids, elevation_m=collect(range(100.0, 1800.0; length=n)),
        slope_deg=5.0 .+ rand(rng, n) .* 20, aspect_sin=randn(rng, n), aspect_cos=randn(rng, n),
    )
    era5 = Dict(variable => randn(rng, n, nt)
        for variable in JointVariableSelection.ERA5VariableSelection.ERA5_VARIABLES)
    y_obs = 2.0 .+ abs.(randn(rng, n, nt))
    y_sat = y_obs .- 0.3 .* randn(rng, n, nt)
    cfg = JointVariableSelection.JointSelectionConfig(
        outdir="unused", min_wet_hours=1, min_stations_per_time=3,
        independent_permutations=9, spatial_permutations=9, bandwidth_candidates=[6, 8, 10],
    )
    joint_cfg = JointCovariateModels.JointCovariateBenchmarkConfig(
        spec_path=nothing, terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[6, 8, 10],
    )
    return (; ids, times, lonlat, terrain, era5, y_obs, y_sat, cfg, joint_cfg)
end

@testset "Nested joint covariate selection: screening uses training stations only" begin
    f = synthetic_joint_benchmark_fixture()
    train_idx = collect(1:22)
    val_idx = collect(23:30)
    first = _screen_joint_subset(
        "FY4B", f.y_obs, f.y_sat, f.terrain, f.era5, nothing, train_idx, f.ids, f.times,
        f.lonlat, f.cfg; scheme="test", fold=1, repeat=1, seed=17, run_seed=4, dem_seed=81,
    )
    changed_obs, changed_sat = copy(f.y_obs), copy(f.y_sat)
    changed_obs[val_idx, :] .= 10_000
    changed_sat[val_idx, :] .= -10_000
    second = _screen_joint_subset(
        "FY4B", changed_obs, changed_sat, f.terrain, f.era5, nothing, train_idx, f.ids, f.times,
        f.lonlat, f.cfg; scheme="test", fold=1, repeat=1, seed=17, run_seed=4, dem_seed=81,
    )
    @test first.role_map == second.role_map
    @test isequal(first.candidates, second.candidates)

    # The tagged tables carry scheme/product/fold/repeat/seed for downstream CSV writing.
    @test all(==("test"), first.roles.scheme)
    @test all(==("FY4B"), first.roles.product)
    @test all(==(1), first.roles.fold)
    @test all(==(1), first.roles.repeat)
    @test all(==(17), first.roles.seed)
    @test Set(first.roles.variable_group) == Set(JointCovariateModels.JOINT_GROUP_ORDER)

    # The returned role_map is a valid drop-in for `build_joint_fold_context`'s `roles` argument.
    context = build_joint_fold_context(
        "FY4B", first.role_map, train_idx, val_idx, f.lonlat, f.y_obs, f.y_sat,
        f.terrain, f.era5, nothing, f.joint_cfg,
    )
    @test Set(context.variables) == Set(keys(first.role_map))
end

@testset "Joint role stability and output writing" begin
    f = synthetic_joint_benchmark_fixture()
    store = _empty_joint_store()
    for (fold, train_idx) in ((1, collect(1:22)), (2, collect(9:30)))
        selection = _screen_joint_subset(
            "FY4B", f.y_obs, f.y_sat, f.terrain, f.era5, nothing, train_idx, f.ids, f.times,
            f.lonlat, f.cfg; scheme="balanced_spatial", fold, repeat=1, seed=7,
            run_seed=fold, dem_seed=100 + fold,
        )
        _store_joint_selection!(store, selection)
    end
    roles = vcat(store[:roles]...; cols=:union)
    @test nrow(roles) == 2 * length(JointCovariateModels.JOINT_GROUP_ORDER)
    stability = _joint_role_stability(roles)
    @test all(==(2), stability.fold_count)
    @test nrow(stability) == length(JointCovariateModels.JOINT_GROUP_ORDER)

    mktempdir() do dir
        _write_joint_selection_outputs(dir, store)
        for filename in (
            "joint_fold_candidates.csv", "joint_fold_vif.csv",
            "joint_fold_spatial_variability.csv", "joint_fold_roles.csv",
            "joint_role_stability.csv",
        )
            @test isfile(joinpath(dir, filename))
        end
    end
end

@testset "Residual shrinkage" begin
    f = synthetic_joint_benchmark_fixture()
    train_idx = collect(1:22)
    val_idx = collect(23:30)
    roles = Dict("elevation" => "local")
    context = JointCovariateModels.build_joint_fold_context(
        "FY4B", roles, train_idx, val_idx, f.lonlat, f.y_obs, f.y_sat,
        f.terrain, f.era5, nothing, f.joint_cfg,
    )
    residuals = f.y_obs[train_idx, :] .- f.y_sat[train_idx, :]
    y_obs_train = f.y_obs[train_idx, :]
    y_sat_train = f.y_sat[train_idx, :]
    times = collect(1:size(residuals, 2))

    bisquare = _kernel_function(BISQUARE)
    unshrunk = _joint_candidate_metrics(
        context, residuals, y_obs_train, y_sat_train, "residual_gwr",
        [8.0], bisquare, times, nothing, nothing, [1.0],
    )
    @test unshrunk.shrink == 1.0

    ladder = _joint_candidate_metrics(
        context, residuals, y_obs_train, y_sat_train, "residual_gwr",
        [8.0], bisquare, times, nothing, nothing, [0.25, 0.5, 0.75, 1.0],
    )
    # The ladder always contains 1.0, so it can never score worse than the unshrunk fit.
    @test ladder.RMSE <= unshrunk.RMSE + 1e-12
    @test ladder.shrink in (0.25, 0.5, 0.75, 1.0)
    # Shrinking rescales the correction; it cannot change which cells are predictable.
    @test ladder.n == unshrunk.n
    @test ladder.coverage == unshrunk.coverage

    # Inflate the residual the model is asked to reproduce so the fitted correction comes out
    # roughly twice as large as the error it should remove. The selector must then pull it back.
    inflated_obs = y_sat_train .+ 2.0 .* (y_obs_train .- y_sat_train)
    inflated = _joint_candidate_metrics(
        context, 2.0 .* residuals, inflated_obs, y_sat_train, "residual_gwr",
        [8.0], bisquare, times, nothing, nothing, [0.25, 0.5, 0.75, 1.0],
    )
    over_corrected = _joint_candidate_metrics(
        context, 2.0 .* residuals, y_obs_train, y_sat_train, "residual_gwr",
        [8.0], bisquare, times, nothing, nothing, [0.25, 0.5, 0.75, 1.0],
    )
    @test inflated.shrink >= over_corrected.shrink
    @test over_corrected.shrink < 1.0

    @test_throws ArgumentError _joint_candidate_metrics(
        context, residuals, y_obs_train, y_sat_train, "residual_gwr",
        [8.0], bisquare, times, nothing, nothing, Float64[],
    )
end

@testset "DEM and joint mixed_gwr/mgwr sweep every kernel and bandwidth family" begin
    # Mirrors the baseline-path parity test above: `select_dem_parameter!`/
    # `select_joint_parameter!` used to fix `mixed_gwr`/`mgwr` at kernel=BISQUARE and, even after
    # the kernel sweep was added, only ever searched the adaptive (neighbor-count) bandwidth
    # family. With more than one kernel and a non-empty fixed-km grid configured, both dimensions
    # must actually vary, matching the baseline path's `gwr`.
    mger = MGERConfig(
        station_meta_path="unused.csv", obs_hourly_wide_path="unused.csv",
        sat_paths=Dict("FY4B" => "unused.csv"), outdir="unused",
        kernels=[GAUSSIAN, BISQUARE], bw_adaptive=[8.0, 16.0], bw_fixed_km=[20.0, 50.0],
    )

    @testset "DEM path" begin
        rng = MersenneTwister(2029)
        n = 40
        lonlat = hcat(collect(range(109.5, 112.0; length=n)), 30.0 .+ rand(rng, n))
        elevation = collect(range(100.0, 1800.0; length=n))
        aspect = rand(rng, n) .* 360
        terrain = DataFrame(
            station_id=string.(1:n), elevation_m=elevation, slope_deg=5.0 .+ rand(rng, n) .* 20,
            aspect_deg=aspect, aspect_sin=sind.(aspect), aspect_cos=cosd.(aspect),
        )
        train_idx, val_idx = 1:32, 33:40
        role_map = Dict("elevation" => "local")
        response = 0.002 .* elevation[train_idx] .+ 0.01 .* randn(rng, length(train_idx))
        selection = (; role_map, response, selection_status="selected")
        dem = DEMTerrainExperiment.DEMExperimentConfig(
            outdir="unused", bandwidth_candidates=[16, 24],
        )
        dem_context = build_dem_fold_context(
            selection, terrain[train_idx, :], terrain[val_idx, :],
            lonlat[train_idx, :], lonlat[val_idx, :], dem,
        )
        cfg = InterpolationBenchmarkConfig(mger=mger, terrain_path="unused", dem=dem)

        for method in ("mixed_gwr", "mgwr")
            scan_rows = NamedTuple[]
            selected = select_dem_parameter!(
                scan_rows, cfg, method, "residual", :balanced_spatial, "TEST", 1, dem_context,
            )
            @test length(unique(row.kernel for row in scan_rows)) > 1
            @test length(unique(row.adaptive for row in scan_rows)) > 1
            @test selected.kernel in cfg.mger.kernels
            @test selected.adaptive isa Bool
        end
    end

    @testset "Joint-covariate path" begin
        f = synthetic_joint_benchmark_fixture()
        train_idx = collect(1:22)
        val_idx = collect(23:30)
        roles = Dict("elevation" => "local")
        joint_context = JointCovariateModels.build_joint_fold_context(
            "FY4B", roles, train_idx, val_idx, f.lonlat, f.y_obs, f.y_sat,
            f.terrain, f.era5, nothing, f.joint_cfg,
        )
        cfg = InterpolationBenchmarkConfig(
            mger=mger, joint_covariates=f.joint_cfg, mgwr_max_tuning_iterations=3,
        )
        y_obs_train = f.y_obs[train_idx, :]
        y_sat_train = f.y_sat[train_idx, :]

        for method in ("mixed_gwr", "mgwr")
            scan_rows = NamedTuple[]
            selected = select_joint_parameter!(
                scan_rows, cfg, method, "residual", :balanced_spatial, "FY4B", 1,
                joint_context, y_obs_train, y_sat_train,
            )
            @test length(unique(row.kernel for row in scan_rows)) > 1
            @test length(unique(row.adaptive for row in scan_rows)) > 1
            @test selected.kernel in cfg.mger.kernels
            @test selected.adaptive isa Bool
        end
    end
end

@testset "Explicit global bandwidth candidate" begin
    # The whole design rests on `bw = Inf` meaning "unweighted OLS over the training stations"
    # for every kernel, so that global can ride the existing bandwidth plumbing instead of
    # needing its own code path. Pin that rather than leaving it implicit in `kernel.jl`.
    @testset "bw = Inf gives every kernel a flat weight of one" begin
        distances = [0.0 5.0; 12.5 0.0; 40.0 120.0]
        for kernel in (GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR)
            @test gw_weight(distances, Inf; kernel=kernel, adaptive=false) == ones(3, 2)
        end
    end

    @testset "global still excludes the held-out station under exclude_self" begin
        # `_gwr_predict` drops self by setting its distance to Inf, which is enough to keep it
        # out of the adaptive neighbour ranking but NOT enough to zero its weight at bw = Inf:
        # `kernel(Inf, Inf)` is NaN for four kernels and 1.0 for boxcar. Boxcar is the one that
        # fails silently - it would hand each station its own value back.
        rng = MersenneTwister(4711)
        n = 12
        lonlat = hcat(110.0 .+ rand(rng, n), 30.0 .+ rand(rng, n))
        values = reshape(collect(1.0:n), n, 1)
        center = (mean(lonlat[:, 1]), mean(lonlat[:, 2]))
        X = build_X_intercept_centered(lonlat; center=center)
        # At bw = Inf every weight is one, so leave-one-out GWR is exactly OLS refitted on the
        # other n-1 stations. That is the invariant a leaked self-weight breaks.
        expected = [only(X[i:i, :] * (X[setdiff(1:n, i), :] \ values[setdiff(1:n, i), 1]))
                    for i in 1:n]
        for kernel in (GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR)
            loo = _gwr_predict(
                lonlat, values, lonlat; kernel=kernel, adaptive=false, bw=Inf, exclude_self=true,
            )
            @test all(isfinite, loo)
            @test vec(loo) ≈ expected
        end
    end

    mger = MGERConfig(
        station_meta_path="unused.csv", obs_hourly_wide_path="unused.csv",
        sat_paths=Dict("FY4B" => "unused.csv"), outdir="unused",
        kernels=[GAUSSIAN, BISQUARE], bw_adaptive=[8.0, 16.0], bw_fixed_km=[20.0, 50.0],
    )

    @testset "_bandwidth_families appends global only when enabled" begin
        with_global = InterpolationBenchmarkConfig(mger=mger)
        without_global = InterpolationBenchmarkConfig(mger=mger, bw_include_global=false)
        @test with_global.bw_include_global

        (adaptive_on, fixed_on) = _bandwidth_families(with_global, mger.bw_adaptive)
        @test adaptive_on == (true, [8.0, 16.0])
        @test fixed_on == (false, [20.0, 50.0, Inf])

        (_, fixed_off) = _bandwidth_families(without_global, mger.bw_adaptive)
        @test fixed_off == (false, [20.0, 50.0])

        # The adaptive family is passed through untouched: the joint path has already filtered
        # its own grid by the time it calls this.
        (adaptive_filtered, _) = _bandwidth_families(with_global, [12])
        @test adaptive_filtered == (true, [12.0])
    end

    @testset "adaptive candidates are capped by the inner selection split" begin
        # 30 training stations split into three inner groups of ten: a candidate is only ever
        # fitted on 20, so anything above that would fall through gw_weight's `dn > 1` branch
        # to a near-global fit during scoring and then refit as a genuine 24-neighbour model.
        groups = [collect(1:10), collect(11:20), collect(21:30)]
        @test _selection_train_minimum(30, groups) == 20
        @test _usable_adaptive_candidates([8.0, 16.0, 20.0, 24.0], 30, groups) ==
            [8.0, 16.0, 20.0]
        # `:loocv` geometry holds out a single station instead of a whole group.
        @test _selection_train_minimum(30, nothing) == 29
        @test _usable_adaptive_candidates([8.0, 24.0], 30, nothing) == [8.0, 24.0]
    end

    @testset "the joint scan actually offers global" begin
        f = synthetic_joint_benchmark_fixture()
        train_idx = collect(1:22)
        val_idx = collect(23:30)
        joint_context = JointCovariateModels.build_joint_fold_context(
            "FY4B", Dict("elevation" => "local"), train_idx, val_idx,
            f.lonlat, f.y_obs, f.y_sat, f.terrain, f.era5, nothing, f.joint_cfg,
        )
        cfg = InterpolationBenchmarkConfig(
            mger=mger, joint_covariates=f.joint_cfg, mgwr_max_tuning_iterations=3,
        )
        scan_rows = NamedTuple[]
        select_joint_parameter!(
            scan_rows, cfg, "gwr", "residual", :balanced_spatial, "FY4B", 1,
            joint_context, f.y_obs[train_idx, :], f.y_sat[train_idx, :],
        )
        global_rows = filter(row -> row.bw == Inf, scan_rows)
        @test !isempty(global_rows)
        @test all(row -> row.adaptive === false, global_rows)
        @test any(row -> row.status == "success", global_rows)
    end
end

@testset "Joint covariate config requires exactly one selection mode" begin
    mger = _dummy_mger()
    dummy_selection = JointVariableSelection.JointSelectionConfig(outdir="unused")
    both = InterpolationBenchmarkConfig(
        mger=mger,
        joint_covariates=JointCovariateModels.JointCovariateBenchmarkConfig(
            spec_path="unused.csv", terrain_path="unused.csv", era5_paths=Dict{Int,String}(),
        ),
        joint_selection=dummy_selection,
    )
    err = try
        _validate_benchmark_config(both, 50)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("exactly one of spec_path", err.msg)

    neither = InterpolationBenchmarkConfig(
        mger=mger,
        joint_covariates=JointCovariateModels.JointCovariateBenchmarkConfig(
            spec_path=nothing, terrain_path="unused.csv", era5_paths=Dict{Int,String}(),
        ),
    )
    err2 = try
        _validate_benchmark_config(neither, 50)
        nothing
    catch e
        e
    end
    @test err2 isa ArgumentError
    @test occursin("exactly one of spec_path", err2.msg)
end

@testset "The claim is assessed on every partition, not just the first" begin
    products = ["TEST"]
    # Partition 1: GWR beats all three baselines on every gate. Partition 2: it is worse overall
    # and loses the heavy-rain gate. A run pinned to repeat 1 would call this supported outright.
    rows = NamedTuple[]
    for repeat_value in 1:2
        gwr_overall = repeat_value == 1 ? 0.9 : 1.6
        gwr_heavy = repeat_value == 1 ? 0.8 : 1.3
        push!(rows, _claim_metric_row(product="TEST", method="residual_gwr",
            group="overall", level="all", RMSE=gwr_overall, repeat=repeat_value))
        for (method, RMSE) in (("idw", 1.5), ("adw", 1.0), ("tps", 1.3))
            push!(rows, _claim_metric_row(product="TEST", method=method,
                group="overall", level="all", RMSE=RMSE, repeat=repeat_value))
        end
        for (method, RMSE) in (
            ("residual_gwr", gwr_heavy), ("idw", 1.5), ("adw", 1.2), ("tps", 1.4),
        )
            push!(rows, _claim_metric_row(product="TEST", method=method,
                group="rain_intensity", level="heavy", RMSE=RMSE, repeat=repeat_value))
        end
        for (method, RMSE) in (("residual_gwr", 1.0), ("idw", 1.02), ("adw", 1.01), ("tps", 1.03))
            push!(rows, _claim_metric_row(product="TEST", method=method,
                group="rain_intensity", level="moderate", RMSE=RMSE, repeat=repeat_value))
        end
        for (method, RMSE) in (("residual_gwr", 0.85), ("idw", 1.4), ("adw", 1.1), ("tps", 1.3))
            push!(rows, _claim_metric_row(product="TEST", method=method,
                group="year", level="2022", RMSE=RMSE, repeat=repeat_value))
        end
        for (method, CSI, FAR) in (
            ("residual_gwr", 0.6, 0.10), ("idw", 0.55, 0.12),
            ("adw", 0.58, 0.11), ("tps", 0.56, 0.13),
        )
            push!(rows, _claim_metric_row(product="TEST", method=method,
                group="event_threshold", level="0.1", RMSE=1.0, CSI=CSI, FAR=FAR,
                repeat=repeat_value))
        end
    end
    metrics = DataFrame(rows)
    bootstrap = DataFrame([
        _claim_bootstrap_row(product="TEST", baseline=baseline, stratum="overall",
            ci_low=0.3, pvalue_holm=0.01, repeat=repeat_value)
        for repeat_value in 1:2 for baseline in ("idw", "adw", "tps")
    ])

    claim = assess_gwr_claim(metrics, bootstrap, products)
    @test nrow(claim) == 2
    @test sort(claim.repeat) == [1, 2]
    first_partition = only(filter(:repeat => ==(1), claim))
    second_partition = only(filter(:repeat => ==(2), claim))
    @test first_partition.direction_improved
    @test first_partition.product_supported
    # The second partition disagrees, which is exactly what pinning to repeat 1 used to hide.
    @test !second_partition.direction_improved
    @test !second_partition.product_supported

    agreement = only(claim_agreement(claim))
    @test agreement.n_repeats == 2
    @test agreement.improved_repeats == 1
    @test agreement.supported_repeats == 1
    @test !agreement.unanimous
    @test agreement.min_overall_relative_improvement <
        agreement.max_overall_relative_improvement
end

@testset "Fixed full-data covariate selection is refused unless declared exploratory" begin
    mger = _dummy_mger()
    leaky = InterpolationBenchmarkConfig(
        mger=mger,
        joint_covariates=JointCovariateModels.JointCovariateBenchmarkConfig(
            spec_path="unused.csv", terrain_path="unused.csv", era5_paths=Dict{Int,String}(),
        ),
    )
    err = try
        _validate_benchmark_config(leaky, 50)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("held-out fold", err.msg)

    # The escape hatch still works - it is an admission, not a fix - and the failure it then hits
    # is the ordinary missing-file one, which proves the leak check ran first and passed.
    declared = InterpolationBenchmarkConfig(
        mger=mger, exploratory_only=true,
        joint_covariates=JointCovariateModels.JointCovariateBenchmarkConfig(
            spec_path="unused.csv", terrain_path="unused.csv", era5_paths=Dict{Int,String}(),
        ),
    )
    err2 = try
        _validate_benchmark_config(declared, 50)
        nothing
    catch e
        e
    end
    @test err2 isa ArgumentError
    @test occursin("joint specification does not exist", err2.msg)
end

@testset "Fold spread summarises metrics_folds, not the pooled row" begin
    rows = NamedTuple[]
    for fold_value in 1:3
        push!(rows, _empty_metric_row(;
            scheme="balanced_spatial", product="FY4B", method="adw", fold=fold_value,
            repeat=1, seed=11, group="overall", level="all",
            n=100, coverage=1.0, RMSE=Float64(fold_value), MAE=Float64(fold_value),
        ))
    end
    # A pooled row sits in the same table and must be ignored, or the spread would be computed
    # over an estimate of itself.
    push!(rows, _empty_metric_row(;
        scheme="balanced_spatial", product="FY4B", method="adw", fold=missing,
        repeat=1, seed=11, group="overall", level="all",
        n=300, coverage=1.0, RMSE=99.0, MAE=99.0,
    ))
    summary = summarize_fold_spread(DataFrame(rows))
    @test nrow(summary) == 1
    row = only(summary)
    @test row.n_folds == 3
    @test row.total_n == 300
    @test row.RMSE_mean ≈ 2.0
    @test row.RMSE_min ≈ 1.0
    @test row.RMSE_max ≈ 3.0
    @test row.RMSE_std ≈ 1.0

    # One fold has no spread to report, and must say NaN rather than 0.
    single = summarize_fold_spread(DataFrame([_empty_metric_row(;
        scheme="balanced_spatial", product="FY4B", method="adw", fold=1,
        repeat=1, seed=11, group="overall", level="all",
        n=100, coverage=1.0, RMSE=1.0, MAE=1.0,
    )]))
    @test isnan(only(single).RMSE_std)
end

@testset "auto compares contenders on one shared mask" begin
    y_obs = [1.0 2.0 3.0 4.0; 1.0 2.0 3.0 4.0]
    y_sat = fill(1.0, 2, 4)
    # `patchy` declines the two hard columns and is mediocre on the two it keeps. `steady`
    # predicts everything: much better than `patchy` on the easy columns, hopeless on the hard
    # ones. This is the shape of the bias the shared mask exists to remove — a method that fails
    # where the problem is hard is never charged for it, and its scan-row RMSE flatters it.
    patchy = y_obs .+ 0.5
    patchy[:, 3:4] .= NaN
    steady = copy(y_obs)
    steady[:, 1:2] .+= 0.1
    steady[:, 3:4] .+= 10.0
    # Scored on its own cells — which is what each method's `selected` scan row carries — the
    # method that gave up looks like the better one.
    @test _candidate_metrics(y_obs, y_sat, patchy).RMSE <
        _candidate_metrics(y_obs, y_sat, steady).RMSE
    # Scored on the shared mask, where neither is charged for the columns `patchy` skipped, it is
    # not: `steady` is genuinely better on the cells both of them predicted.
    choice = select_auto_method(
        [(; method="patchy", prediction=patchy), (; method="steady", prediction=steady)],
        y_obs, y_sat,
    )
    @test choice.chosen == "steady"
    @test choice.runner_up == "patchy"
    @test choice.n == 4
    @test choice.shared_mask_coverage ≈ 0.5
    @test choice.chosen_rmse ≈ 0.1

    # An undeduplicated contender list makes `auto` look like a wider choice than it was:
    # `residual_gwr` and `mixed_gwr` are the same model whenever the fold's role map has no
    # "global" role, so one of them is the other's twin at a margin of zero.
    twin = select_auto_method(
        [(; method="residual_gwr", prediction=steady),
         (; method="mixed_gwr", prediction=copy(steady)),
         (; method="patchy", prediction=patchy)],
        y_obs, y_sat,
    )
    @test twin.n_contenders == 2                       # not 3
    @test twin.collapsed == "mixed_gwr=residual_gwr"   # first occurrence keeps the name
    @test twin.chosen == "residual_gwr"
    @test twin.runner_up == "patchy"                   # a real alternative, not the twin
    # Collapsing must not change what wins or what it scored.
    plain = select_auto_method(
        [(; method="residual_gwr", prediction=steady), (; method="patchy", prediction=patchy)],
        y_obs, y_sat,
    )
    @test twin.chosen == plain.chosen
    @test twin.chosen_rmse == plain.chosen_rmse
    @test twin.n == plain.n
    @test isempty(plain.collapsed)
    # `isequal`, not `==`: these matrices carry NaN, and two contenders that both gave up on the
    # same cells are still the same contender.
    @test select_auto_method(
        [(; method="a", prediction=patchy), (; method="b", prediction=copy(patchy))],
        y_obs, y_sat,
    ).n_contenders == 1
    # Predictions that differ anywhere are two contenders, however close.
    nudged = copy(steady); nudged[1, 1] = nextfloat(nudged[1, 1])
    @test select_auto_method(
        [(; method="a", prediction=steady), (; method="b", prediction=nudged)], y_obs, y_sat,
    ).n_contenders == 2

    # Nothing to score on: the caller must be told, not handed a default.
    @test select_auto_method(
        [(; method="empty", prediction=fill(NaN, 2, 4))], y_obs, y_sat,
    ) === nothing
    @test select_auto_method(NamedTuple[], y_obs, y_sat) === nothing
end

@testset "Interpolation benchmark fixture" begin
    mktempdir() do temp_dir
        ids = string.(2001:2016)
        lon = [110.0 + 0.12 * mod(i - 1, 4) for i in 1:16]
        lat = [30.0 + 0.11 * div(i - 1, 4) for i in 1:16]
        lonlat = hcat(lon, lat)
        times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 7))
        station_path = joinpath(temp_dir, "stations.csv")
        CSV.write(station_path, DataFrame(station_id=ids, lon=lon, lat=lat))
        aspect = collect(range(1.0, 359.0; length=length(ids)))
        terrain_path = joinpath(temp_dir, "station_terrain.csv")
        CSV.write(terrain_path, DataFrame(
            station_id=ids,
            elevation_m=200.0 .+ 30.0 .* collect(1:length(ids)),
            slope_deg=2.0 .+ mod.(collect(1:length(ids)), 7),
            aspect_deg=aspect, aspect_sin=sind.(aspect), aspect_cos=cosd.(aspect),
        ))

        obs = Matrix{Float64}(undef, 16, length(times))
        for station in 1:16, time in eachindex(times)
            obs[station, time] = max(
                0.0, 0.5 + 0.08 * station + 0.3 * sin(time / 2) + 0.04 * mod(station * time, 3),
            )
        end
        obs_path = joinpath(temp_dir, "obs.csv")
        _write_benchmark_wide(obs_path, times, ids, obs)
        sat_paths = Dict{String,String}()
        for (product, bias) in (("FY4B", 0.3), ("GPM", -0.2), ("GSMaP", 0.1))
            sat = max.(obs .+ bias .+ 0.03 .* reshape(1:16, :, 1), 0.0)
            path = joinpath(temp_dir, "$(lowercase(product)).csv")
            _write_benchmark_wide(path, times, ids, sat)
            sat_paths[product] = path
        end

        outdir = joinpath(temp_dir, "benchmark")
        mger = MGERConfig(
            station_meta_path=station_path,
            obs_hourly_wide_path=obs_path,
            sat_paths=sat_paths,
            outdir=outdir,
            kernels=[GAUSSIAN],
            bw_adaptive=[5.0],
            bw_fixed_km=[100.0],
            expected_common_time_count=length(times),
        )
        cfg = InterpolationBenchmarkConfig(
            mger=mger,
            terrain_path=terrain_path,
            dem=DEMTerrainExperiment.DEMExperimentConfig(
                outdir=outdir, min_wet_hours=1, k=2,
                bandwidth_candidates=[5], screen_permutations=9,
                spatial_permutations=9, q_threshold=0.05, max_iterations=200,
            ),
            k=2,
            seed=7,
            cv_schemes=[:balanced_spatial],
            idw_powers=[2.0],
            neighbor_candidates=Union{Nothing,Int}[4],
            tps_smooth_candidates=[0.01],
            min_tuning_coverage=0.8,
            bootstrap_reps=20,
        )
        result = run_interpolation_benchmark(cfg)
        @test !isempty(result.metrics)
        @test nrow(result.status) == 2 * 3 * 8
        @test Set(result.status.method) == Set([
            "idw", "adw", "tps", "gwr",
            "residual_gwr", "mixed_gwr", "mgwr", "auto",
        ])
        fitted_status = filter(:method => !=("auto"), result.status)
        @test all(fitted_status.status .== "success")
        # This fixture is the legacy DEM path with 8 training stations per fold, so `auto` has
        # neither a `dem_context` nor an inner split to choose on. It must say so rather than
        # failing: "skipped" is what keeps a real breakage visible elsewhere.
        auto_status = filter(:method => ==("auto"), result.status)
        @test all(auto_status.status .== "skipped")
        @test all(occursin("legacy DEM path", row.error) for row in eachrow(auto_status))
        @test isempty(result.auto_selection)
        @test Set(result.metrics.method) == Set(BENCHMARK_METHODS)
        @test Set((row.mode, row.method) for row in eachrow(result.scans)) == Set(BENCHMARK_RUNS)
        @test Set(result.scans.group[result.scans.method .== "mgwr"]) ==
            Set(["intercept", "longitude", "latitude"])
        @test all(group -> count(group.selected) == 1,
            groupby(result.scans, [:repeat, :scheme, :product, :fold, :mode, :method, :group]))
        @test all(method -> method in TRADITIONAL_METHODS, result.claim.best_traditional)
        @test isfile(joinpath(outdir, "metrics_stratified.csv"))
        @test isfile(joinpath(outdir, "parameter_scan.csv"))
        @test isfile(joinpath(outdir, "paired_comparisons.csv"))
        # A single repeat keeps the historical layout and stamps every row with repeat 1.
        @test all(result.metrics.repeat .== 1)
        @test all(result.scans.repeat .== 1)
        @test all(result.status.repeat .== 1)
        @test all(result.metrics.seed .== cfg.seed)
        @test isdir(joinpath(outdir, "balanced_spatial"))
        @test !isdir(joinpath(outdir, "repeat_01"))
        @test all(result.repeat_summary.n_repeats .== 1)
        for filename in (
            "dem_fold_quality_control.csv", "dem_fold_correlation.csv",
            "dem_fold_aspect_joint_test.csv", "dem_fold_vif.csv",
            "dem_fold_monthly_correlation.csv", "dem_fold_spatial_variability.csv",
            "dem_fold_roles.csv", "dem_role_stability.csv",
            "dem_final_full_data_spec.csv", "dem_bandwidths.csv", "model_status.csv",
        )
            @test isfile(joinpath(outdir, filename))
        end
        roles = CSV.read(joinpath(outdir, "dem_fold_roles.csv"), DataFrame)
        @test nrow(roles) == 3 * 3 + 2 * 3 * 3
        @test Set(roles.variable_group) == Set(["elevation", "slope", "aspect"])
        final_spec = CSV.read(joinpath(outdir, "dem_final_full_data_spec.csv"), DataFrame)
        @test nrow(final_spec) == 3 * 3
        stability = CSV.read(joinpath(outdir, "dem_role_stability.csv"), DataFrame)
        @test nrow(stability) == 3 * 3
        @test all(stability.fold_count .== 2)
        dem_bandwidths = CSV.read(joinpath(outdir, "dem_bandwidths.csv"), DataFrame)
        @test Set(dem_bandwidths.method) ==
            Set(["residual_gwr", "mixed_gwr", "mgwr"])
        residual_status = filter(:method => in(["residual_gwr", "mixed_gwr", "mgwr"]), result.status)
        @test all(residual_status.dem_variable_count .== 0)
        @test all(residual_status.dem_selection_status .== "no_dem_selected")
        @test all(group -> length(unique(group.dem_variables)) == 1,
            groupby(residual_status, [:scheme, :product, :fold]))
        for product in ("fy4b", "gpm", "gsmap"), method in BENCHMARK_METHODS
            path = joinpath(outdir, "balanced_spatial", product, "oof_$(method).csv")
            @test isfile(path)
            output = CSV.read(path, DataFrame)
            values = Matrix{Float64}(output[:, Not(:time)])
            @test all(value -> isnan(value) || value >= 0, values)
        end
        for product in ("fy4b", "gpm", "gsmap"), method in
            ("residual_idw", "residual_adw", "residual_tps")
            @test !isfile(joinpath(outdir, "balanced_spatial", product, "oof_$(method).csv"))
        end
        scope = CSV.read(joinpath(outdir, "benchmark_scope.csv"), DataFrame)
        common_mask_scope = only(scope.value[scope.key .== "common_evaluation_mask"])
        # The scope row now states what the rule produced, not only what the rule was. The mask
        # is an intersection over MASK_METHODS, so one masked method's coverage moving re-scores
        # every method including the traditional baselines; without a number here, nothing in a
        # run's own output changes when that happens.
        @test startswith(common_mask_scope,
            "true across all eight methods; cells kept of obs+satellite evaluable: ")
        for product in ("fy4b", "gpm", "gsmap")
            written = CSV.read(
                joinpath(outdir, "balanced_spatial", product, "common_evaluation_mask.csv"),
                DataFrame,
            )
            kept = count(Matrix(written[:, Not(:time)]))
            # `product` here is the directory name, which is `lowercase(product)`; the scope
            # row carries the product's own casing (GSMaP, not GSMAP), so compare lowercased.
            @test occursin(
                "balanced_spatial/$(product) $(kept)/", lowercase(common_mask_scope),
            )
        end
        @test only(scope.value[scope.key .== "repeated_cv_partitions"]) == "1"

        # Two partitions: rows are stamped per repeat and the spread tables are populated.
        repeat_outdir = joinpath(temp_dir, "benchmark_repeats")
        repeated = InterpolationBenchmarkConfig(
            mger=MGERConfig(
                station_meta_path=station_path, obs_hourly_wide_path=obs_path,
                sat_paths=sat_paths, outdir=repeat_outdir, kernels=[GAUSSIAN],
                bw_adaptive=[5.0], bw_fixed_km=[100.0],
                expected_common_time_count=length(times),
            ),
            k=2, seed=7, seeds=[7, 8], cv_schemes=[:balanced_spatial],
            idw_powers=[2.0], neighbor_candidates=Union{Nothing,Int}[4],
            tps_smooth_candidates=[0.01], min_tuning_coverage=0.8, bootstrap_reps=0,
        )
        repeat_result = run_interpolation_benchmark(repeated)
        @test sort(unique(repeat_result.metrics.repeat)) == [1, 2]
        @test sort(unique(repeat_result.metrics.seed)) == [7, 8]
        @test sort(unique(repeat_result.scans.repeat)) == [1, 2]
        @test all(group -> count(group.selected) == 1, groupby(
            repeat_result.scans, [:repeat, :scheme, :product, :fold, :mode, :method, :group]))
        @test all(repeat_result.repeat_summary.n_repeats .== 2)
        @test all(repeat_result.rank_stability.n_repeats .== 2)
        @test isfile(joinpath(repeat_outdir, "metrics_repeat_summary.csv"))
        @test isfile(joinpath(repeat_outdir, "method_rank_stability.csv"))
        # Repeats nest under repeat_<i>; only the first writes the large OOF tables.
        @test isdir(joinpath(repeat_outdir, "repeat_01", "balanced_spatial"))
        @test isdir(joinpath(repeat_outdir, "repeat_02", "balanced_spatial"))
        @test isfile(joinpath(repeat_outdir, "repeat_01", "balanced_spatial", "fy4b", "oof_adw.csv"))
        @test !isfile(joinpath(repeat_outdir, "repeat_02", "balanced_spatial", "fy4b", "oof_adw.csv"))
        repeat_scope = CSV.read(joinpath(repeat_outdir, "benchmark_scope.csv"), DataFrame)
        @test only(repeat_scope.value[repeat_scope.key .== "repeated_cv_partitions"]) == "2"
        @test only(repeat_scope.value[repeat_scope.key .== "repeated_cv_seeds"]) == "7,8"
        @test only(repeat_scope.value[repeat_scope.key .== "tuning_time_weighting"]) ==
            only(scope.value[scope.key .== "tuning_time_weighting"])
        @test startswith(
            only(scope.value[scope.key .== "tuning_time_weighting"]), "uniform")
    end
end

@testset "Inner spatial selection split" begin
    n_station = 60
    ids = string.(5001:(5000 + n_station))
    lon = [110.0 + 0.15 * mod(i - 1, 10) for i in 1:n_station]
    lat = [30.0 + 0.15 * div(i - 1, 10) for i in 1:n_station]
    lonlat = hcat(lon, lat)
    cfg = InterpolationBenchmarkConfig(
        mger=MGERConfig(station_meta_path="a", obs_hourly_wide_path="b",
            sat_paths=Dict("FY4B" => "c"), outdir="d"),
        k=5,
    )

    @testset "partitions the training set exactly once" begin
        groups = selection_folds(cfg, :balanced_spatial, ids, lonlat, 1, cfg.seed)
        @test groups !== nothing
        @test length(groups) == cfg.k
        flat = reduce(vcat, groups)
        @test sort(flat) == collect(1:n_station)
        @test length(unique(flat)) == n_station
        # Same splitter as the outer partition, so the inner geometry matches by construction.
        @test all(group -> !isempty(group), groups)
        # The partition must not be the same one in every outer fold.
        other = selection_folds(cfg, :balanced_spatial, ids, lonlat, 2, cfg.seed)
        @test Set(Set.(groups)) != Set(Set.(other))
        # `:random` gets a random inner split — its reporting geometry already matches LOOCV.
        random_groups = selection_folds(cfg, :random, ids, lonlat, 1, cfg.seed)
        @test sort(reduce(vcat, random_groups)) == collect(1:n_station)
    end

    @testset "declines to split a fold that is too small" begin
        # 8 stations split five ways would leave too few to fit the widest design; the caller
        # falls back to leave-one-out and records it rather than failing or degrading silently.
        @test selection_folds(cfg, :balanced_spatial, ids[1:8], lonlat[1:8, :], 1, cfg.seed) ===
            nothing
        @test selection_folds(cfg, :balanced_spatial, ids[1:12], lonlat[1:12, :], 1, cfg.seed) ===
            nothing
        @test selection_folds(cfg, :balanced_spatial, ids[1:40], lonlat[1:40, :], 1, cfg.seed) !==
            nothing
    end

    @testset "_selection_oof assembles every training row" begin
        groups = selection_folds(cfg, :balanced_spatial, ids, lonlat, 1, cfg.seed)
        n_time = 4
        seen_train = Set{Int}()
        oof = _selection_oof(groups, n_station, n_time, (tr, va) -> begin
            # Inner train and inner val must be complementary, never overlapping.
            @test isempty(intersect(tr, va))
            @test sort(vcat(tr, va)) == collect(1:n_station)
            union!(seen_train, tr)
            fill(7.5, length(va), n_time)
        end)
        @test size(oof) == (n_station, n_time)
        @test all(oof .== 7.5)
        @test !any(isnan, oof)
        # Every station is used for fitting in at least one inner group.
        @test length(seen_train) == n_station
    end

    @testset "config validation" begin
        base = InterpolationBenchmarkConfig(
            mger=MGERConfig(station_meta_path="a", obs_hourly_wide_path="b",
                sat_paths=Dict("FY4B" => "c"), outdir="d"))
        @test base.tuning_geometry === :inner_spatial
        @test base.tuning_inner_k == 0
        for bad in (:loo, :spatial, :none)
            @test_throws ArgumentError _validate_benchmark_config(
                InterpolationBenchmarkConfig(mger=base.mger, tuning_geometry=bad), 60)
        end
        @test_throws ArgumentError _validate_benchmark_config(
            InterpolationBenchmarkConfig(mger=base.mger, tuning_inner_k=1), 60)
        # 0 means "use cfg.k", and any k_inner >= 2 is accepted.
        @test _validate_benchmark_config(
            InterpolationBenchmarkConfig(mger=base.mger, tuning_inner_k=3), 60) === nothing
    end
end

@testset "Benchmark runs end to end on the inner spatial split" begin
    mktempdir() do temp_dir
        # Deliberately larger than the main fixture: 60 stations with k=3 leave 40 training
        # stations and 27 inner-training stations, enough for the widest design, so this exercises
        # the shipped default rather than falling back to leave-one-out.
        n_station = 60
        ids = string.(6001:(6000 + n_station))
        lon = [110.0 + 0.15 * mod(i - 1, 10) for i in 1:n_station]
        lat = [30.0 + 0.15 * div(i - 1, 10) for i in 1:n_station]
        times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 7))
        station_path = joinpath(temp_dir, "stations.csv")
        CSV.write(station_path, DataFrame(station_id=ids, lon=lon, lat=lat))

        obs = Matrix{Float64}(undef, n_station, length(times))
        for station in 1:n_station, time in eachindex(times)
            obs[station, time] = max(
                0.0, 0.5 + 0.03 * station + 0.3 * sin(time / 2) + 0.04 * mod(station * time, 3))
        end
        obs_path = joinpath(temp_dir, "obs.csv")
        _write_benchmark_wide(obs_path, times, ids, obs)
        sat_paths = Dict{String,String}()
        for (product, bias) in (("FY4B", 0.3), ("GPM", -0.2))
            path = joinpath(temp_dir, "$(lowercase(product)).csv")
            _write_benchmark_wide(path, times, ids,
                max.(obs .+ bias .+ 0.02 .* reshape(1:n_station, :, 1), 0.0))
            sat_paths[product] = path
        end

        function run_with(geometry::Symbol)
            outdir = joinpath(temp_dir, "benchmark_$(geometry)")
            cfg = InterpolationBenchmarkConfig(
                mger=MGERConfig(
                    station_meta_path=station_path, obs_hourly_wide_path=obs_path,
                    sat_paths=sat_paths, outdir=outdir, kernels=[GAUSSIAN],
                    bw_adaptive=[8.0, 16.0], bw_fixed_km=[20.0, 100.0],
                    expected_common_time_count=length(times),
                ),
                k=3, seed=11, cv_schemes=[:balanced_spatial],
                idw_powers=[2.0], neighbor_candidates=Union{Nothing,Int}[8],
                tps_smooth_candidates=[0.01], min_tuning_coverage=0.8, bootstrap_reps=0,
                tuning_geometry=geometry,
            )
            return cfg, run_interpolation_benchmark(cfg)
        end

        cfg_inner, inner = run_with(:inner_spatial)
        @test !isempty(inner.metrics)
        # Every method that ran to completion still produces exactly one selected candidate per
        # cell. (`mgwr` legitimately fails on a fixture this small — see the parity check below,
        # which is what guards against the inner split breaking a method.)
        succeeded = Set((row.product, row.fold, row.method)
            for row in eachrow(inner.status) if row.status == "success")
        @test all(groupby(inner.scans, [:scheme, :product, :fold, :mode, :method, :group])) do group
            method = group.method[1] in ("mixed_gwr", "mgwr") ? group.method[1] :
                (group.mode[1] == "direct" ? group.method[1] : "residual_$(group.method[1])")
            (group.product[1], group.fold[1], method) in succeeded || return true
            return count(group.selected) == 1
        end
        scope = CSV.read(joinpath(cfg_inner.mger.outdir, "benchmark_scope.csv"), DataFrame)
        geometry_note = only(scope.value[scope.key .== "tuning_geometry"])
        @test startswith(geometry_note, "inner_spatial")
        # The fold was big enough, so nothing fell back — otherwise the scope row says so and
        # this fixture would not be testing the shipped default at all.
        @test !occursin("FELL BACK", geometry_note)

        # The legacy criterion still runs and reaches different conclusions.
        cfg_loocv, loocv = run_with(:loocv)
        legacy_scope = CSV.read(joinpath(cfg_loocv.mger.outdir, "benchmark_scope.csv"), DataFrame)
        @test startswith(
            only(legacy_scope.value[legacy_scope.key .== "tuning_geometry"]), "loocv")
        # Parity: switching the criterion must not make any method start failing. Whatever fails
        # under the inner split has to fail under leave-one-out too, or the split broke it.
        failures(result) = Set((row.product, row.fold, row.method)
            for row in eachrow(result.status) if row.status != "success")
        @test failures(inner) ⊆ failures(loocv)

        key = [:scheme, :product, :fold, :mode, :method, :group]
        picks(scans) = select(filter(:selected => identity, scans), key..., :bw, :power, :smooth)
        joined = innerjoin(picks(inner.scans), picks(loocv.scans), on=key, makeunique=true)
        @test nrow(joined) > 0
        # The two criteria are different estimands, so they must not be trivially identical.
        @test any(.!isequal.(joined.bw, joined.bw_1)) ||
            any(.!isequal.(joined.power, joined.power_1)) ||
            any(.!isequal.(joined.smooth, joined.smooth_1))
    end
end

@testset "Tuning time weighting" begin
    Random.seed!(20260822)
    n_station, n_time, n_wet = 20, 5000, 350
    wet_hours = randperm(n_time)[1:n_wet]
    is_wet = falses(n_time)
    is_wet[wet_hours] .= true

    y_obs = zeros(Float64, n_station, n_time)
    prediction = zeros(Float64, n_station, n_time)
    for time in 1:n_time
        # Errors differ by an order of magnitude between the classes, so an estimator that gets
        # the class mixture wrong cannot accidentally land on the right pooled answer.
        y_obs[:, time] .= is_wet[time] ? 8.0 : 0.0
        prediction[:, time] .= y_obs[1, time] + (is_wet[time] ? 4.0 : 0.5)
    end
    y_sat = zeros(Float64, n_station, n_time)

    @testset "stratified sample" begin
        indices, weights = _tuning_time_sample(y_obs, 336, :stratified)
        # The historical sample took the union of two overlapping halves, so it could collapse
        # below the requested size; disjoint strata cannot.
        @test length(indices) == 336
        @test length(unique(indices)) == 336
        @test issorted(indices)
        @test length(weights) == 336
        # Certainty stratum weighs 1, the subsampled remainder carries the inverse inclusion
        # probability, and the design totals the full hour count.
        certainty = div(336, 8)
        remainder = 336 - certainty
        @test sort(unique(round.(weights; digits=6))) ==
            sort(unique(round.([1.0, (n_time - certainty) / remainder]; digits=6)))
        @test count(≈(1.0), weights) == certainty
        @test sum(weights) ≈ n_time
        # The certainty stratum holds only wet hours; the weighted remainder is drawn from the
        # dry majority, which is where the reported metric actually lives.
        @test all(is_wet[indices[weights .≈ 1.0]])
        @test count(is_wet[indices[weights .> 1.0]]) < 20
    end

    @testset "uniform sample reproduces the legacy objective" begin
        indices, weights = _tuning_time_sample(y_obs, 336, :uniform)
        @test weights === nothing
        @test indices == _tuning_time_indices(y_obs, 336)
        # No subsetting at all when the full record already fits inside the budget.
        short, short_weights = _tuning_time_sample(y_obs[:, 1:100], 336, :stratified)
        @test short == collect(1:100)
        @test short_weights === nothing
    end

    @testset "weighted metric estimates the full-sample metric" begin
        full = _candidate_metrics(y_obs, y_sat, prediction; require_satellite=false)

        indices, weights = _tuning_time_sample(y_obs, 336, :stratified)
        fixed = _candidate_metrics(
            y_obs[:, indices], y_sat[:, indices], prediction[:, indices];
            require_satellite=false, time_weights=weights,
        )
        legacy_indices, _ = _tuning_time_sample(y_obs, 336, :uniform)
        legacy = _candidate_metrics(
            y_obs[:, legacy_indices], y_sat[:, legacy_indices], prediction[:, legacy_indices];
            require_satellite=false,
        )

        # The weighted subsample recovers the pooled metric; the legacy one scores the wet class
        # it oversampled and lands more than twice too high.
        @test isapprox(fixed.RMSE, full.RMSE; rtol=0.05)
        @test isapprox(fixed.MAE, full.MAE; rtol=0.10)
        @test legacy.RMSE > 2 * full.RMSE
        @test legacy.MAE > 2 * full.MAE

        # Coverage stays a raw cell count so `min_tuning_coverage` keeps its meaning, and
        # all-ones weights must not perturb the unweighted path.
        @test fixed.n == n_station * 336
        @test fixed.coverage == 1.0
        ones_weighted = _candidate_metrics(
            y_obs, y_sat, prediction;
            require_satellite=false, time_weights=ones(Float64, n_time),
        )
        @test ones_weighted.RMSE ≈ full.RMSE
        @test ones_weighted.MAE ≈ full.MAE
        @test ones_weighted.n == full.n
        @test ones_weighted.coverage == full.coverage
        @test_throws ArgumentError _candidate_metrics(
            y_obs, y_sat, prediction;
            require_satellite=false, time_weights=ones(Float64, n_time - 1),
        )
    end

    @testset "config validation" begin
        base = InterpolationBenchmarkConfig(
            mger=MGERConfig(
                station_meta_path="a", obs_hourly_wide_path="b",
                sat_paths=Dict("FY4B" => "c"), outdir="d",
            ),
        )
        # `:uniform` is the default on purpose — see the note in `_tuning_time_sample`.
        @test base.tuning_time_weighting === :uniform
        @test_throws ArgumentError _validate_benchmark_config(
            InterpolationBenchmarkConfig(
                mger=base.mger, tuning_time_weighting=:wet_only,
            ), 10,
        )
    end
end

@testset "No selection step sees the outer held-out fold" begin
    # The invariant the whole selection design rests on: corrupt the observations at one fold's
    # validation stations and every hyperparameter chosen for that fold must be bit-identical.
    #
    # Scoped to the target fold on purpose. Those stations are *training* data for the other two
    # folds, so their selections legitimately move - that is also what proves the corruption
    # reached the pipeline rather than being written to a file nobody read.
    #
    # This covers kernel, bandwidth family, bandwidth, IDW power/neighbours, TPS smoothing,
    # residual shrinkage and the `auto` model choice. Covariate and local/global role selection
    # have their own equivalent at "Nested joint covariate selection: screening uses training
    # stations only", which corrupts validation rows of population-wide matrices the same way.
    mktempdir() do temp_dir
        n_station, k, seed, target_fold = 60, 3, 11, 2
        # Pinned, not left to the default: `target_rows` below only names the target fold's
        # stations if the partition here is the one the pipeline builds, and a fixture that
        # silently corrupts another fold's *training* data would pass this test vacuously.
        center_init = :hilbert
        ids = string.(7001:(7000 + n_station))
        lon = [110.0 + 0.15 * mod(i - 1, 10) for i in 1:n_station]
        lat = [30.0 + 0.15 * div(i - 1, 10) for i in 1:n_station]
        lonlat = hcat(lon, lat)
        times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 7))
        station_path = joinpath(temp_dir, "stations.csv")
        CSV.write(station_path, DataFrame(station_id=ids, lon=lon, lat=lat))

        obs = Matrix{Float64}(undef, n_station, length(times))
        for station in 1:n_station, time in eachindex(times)
            obs[station, time] = max(
                0.0, 0.5 + 0.03 * station + 0.3 * sin(time / 2) + 0.04 * mod(station * time, 3))
        end
        sat = max.(obs .+ 0.3 .+ 0.02 .* reshape(1:n_station, :, 1), 0.0)

        folds = benchmark_folds(:balanced_spatial, ids, lonlat; k=k, seed=seed, center_init)
        position = Dict(id => index for (index, id) in enumerate(ids))
        target_rows = [position[id] for id in folds[target_fold]]
        @test !isempty(target_rows)

        function run_with(label::String, observations::Matrix{Float64})
            obs_path = joinpath(temp_dir, "obs_$(label).csv")
            sat_path = joinpath(temp_dir, "fy4b_$(label).csv")
            _write_benchmark_wide(obs_path, times, ids, observations)
            _write_benchmark_wide(sat_path, times, ids, sat)
            cfg = InterpolationBenchmarkConfig(
                mger=MGERConfig(
                    station_meta_path=station_path, obs_hourly_wide_path=obs_path,
                    sat_paths=Dict("FY4B" => sat_path),
                    outdir=joinpath(temp_dir, "run_$(label)"),
                    kernels=[GAUSSIAN, BISQUARE],
                    bw_adaptive=[8.0, 16.0], bw_fixed_km=[20.0, 100.0],
                    expected_common_time_count=length(times),
                ),
                k=k, seed=seed, fold_center_init=center_init, cv_schemes=[:balanced_spatial],
                idw_powers=[1.5, 2.0], neighbor_candidates=Union{Nothing,Int}[8, 16],
                tps_smooth_candidates=[0.01, 0.1], min_tuning_coverage=0.8, bootstrap_reps=0,
            )
            return run_interpolation_benchmark(cfg)
        end

        clean = run_with("clean", obs)
        corrupted_obs = copy(obs)
        corrupted_obs[target_rows, :] .= 9_999.0
        corrupted = run_with("corrupted", corrupted_obs)

        key = [:scheme, :product, :fold, :mode, :method, :group]
        tuned = [:kernel, :adaptive, :bw, :power, :smooth, :shrink]
        picks(scans) = select(
            filter(row -> row.selected && row.fold == target_fold, scans), key..., tuned...,
        )
        clean_picks = picks(clean.scans)
        @test nrow(clean_picks) > 0
        joined = innerjoin(clean_picks, picks(corrupted.scans), on=key, makeunique=true)
        @test nrow(joined) == nrow(clean_picks)
        for column in tuned
            @test all(isequal.(joined[!, column], joined[!, Symbol(column, :_1)]))
        end

        # `auto` is a selection step too, and its input is the inner split alone.
        auto_pick(result) = only(filter(
            row -> row.fold == target_fold, result.auto_selection,
        ).chosen)
        @test auto_pick(clean) == auto_pick(corrupted)

        # Sentinel: the corruption really did reach the run. The held-out fold's own score has to
        # move (its observations changed), and the other folds - for which those stations are
        # training data - are free to select differently.
        fold_rmse(result) = only(filter(row ->
            row.fold === target_fold && row.group == "overall" && row.level == "all" &&
            row.method == "adw", result.metrics,
        ).RMSE)
        @test !isapprox(fold_rmse(clean), fold_rmse(corrupted))
    end
end
