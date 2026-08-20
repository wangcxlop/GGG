#!/usr/bin/env julia

using Test
using CSV, DataFrames, Dates
using Random

const BENCHMARK_TEST_ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(BENCHMARK_TEST_ROOT, "src"))

using MixedGWR
include(joinpath(BENCHMARK_TEST_ROOT, "src", "InterpolationBenchmark.jl"))

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
_claim_metric_row(; product, method, group, level, RMSE, CSI=0.5, FAR=0.1) = _empty_metric_row(;
    scheme="balanced_spatial", product, method, fold=missing, repeat=1, seed=1,
    group, level, n=100, coverage=1.0, RMSE, MAE=RMSE, CSI, FAR,
)

"""One `paired_bootstrap_rows`-shaped row for `assess_gwr_claim` test fixtures."""
_claim_bootstrap_row(; product, baseline, stratum, ci_low, pvalue_holm) = (;
    scheme="balanced_spatial", product, stratum, baseline, repeat=1, seed=1,
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
        @test nrow(result.status) == 2 * 3 * 7
        @test all(result.status.status .== "success")
        @test Set(result.status.method) ==
            Set(["idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr"])
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
        @test common_mask_scope == "true across all eight methods"
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
    end
end
