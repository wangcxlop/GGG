#!/usr/bin/env julia

using Test
using CSV, DataFrames, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
include(joinpath(ROOT, "src", "MGERPipeline.jl"))

const TEST_KERNELS = [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR]
const TEST_KERNEL_NAMES = ["gaussian", "exponential", "bisquare", "tricube", "boxcar"]

function write_wide_fixture(path, times, station_ids, values)
    df = DataFrame(time=times)
    for (index, station_id) in enumerate(station_ids)
        df[!, Symbol(station_id)] = values[index, :]
    end
    CSV.write(path, df)
end

function make_fixture(root)
    station_ids = string.(1001:1020)
    n_station = length(station_ids)
    times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 11))
    n_time = length(times)

    lon = [110.0 + 0.18 * mod(i - 1, 5) for i in 1:n_station]
    lat = [30.0 + 0.16 * div(i - 1, 5) for i in 1:n_station]
    station_path = joinpath(root, "stations.csv")
    CSV.write(station_path, DataFrame(station_id=station_ids, lon=lon, lat=lat))

    obs = Matrix{Float64}(undef, n_station, n_time)
    for i in 1:n_station, t in 1:n_time
        obs[i, t] = 1.2 + 0.04 * i + 0.25 * sin(t / 2) + 0.12 * mod(i + t, 4)
    end
    obs_path = joinpath(root, "obs.csv")
    write_wide_fixture(obs_path, times, station_ids, obs)

    sat_paths = Dict{String, String}()
    for (product, bias, scale) in (("FY4B", 0.45, 0.08), ("GPM", -0.25, 0.06), ("GSMaP", 0.20, 0.05))
        sat = Matrix{Float64}(undef, n_station, n_time)
        for i in 1:n_station, t in 1:n_time
            sat[i, t] = max(obs[i, t] + bias + scale * cos(i * t / 3), 0.0)
        end
        path = joinpath(root, "$(lowercase(product)).csv")
        write_wide_fixture(path, times, station_ids, sat)
        sat_paths[product] = path
    end
    return (; station_ids, times, station_path, obs_path, sat_paths)
end

function fixture_config(fixture, outdir; kernels=copy(TEST_KERNELS))
    return MGERConfig(
        station_meta_path=fixture.station_path,
        obs_hourly_wide_path=fixture.obs_path,
        sat_paths=copy(fixture.sat_paths),
        outdir=outdir,
        kernels=kernels,
        bw_adaptive=[6.0, 9.0],
        bw_fixed_km=[60.0, 200.0],
        expected_common_time_count=length(fixture.times),
    )
end

@testset "Target-centred local linear residual GWR" begin
    train_lonlat = [
        109.9 30.0
        110.1 30.0
        110.0 29.9
        110.0 30.1
    ]
    target_lonlat = reshape([110.0, 30.0], 1, 2)
    weights = ones(Float64, 4, 1)

    constant_residual = fill(2.5, 4, 2)
    constant_residual[1, 2] = NaN
    constant_prediction = local_linear_residual_predict(
        train_lonlat, constant_residual, target_lonlat, weights;
        slope_ridge=1e-6,
    )
    @test constant_prediction[1, 1] ≈ 2.5 atol=1e-12
    @test constant_prediction[1, 2] ≈ 2.5 atol=1e-10

    east, north = target_centered_offsets_km(train_lonlat, 110.0, 30.0)
    planar_residual = reshape(3.0 .+ 0.2 .* east .- 0.1 .* north, :, 1)
    planar_prediction = local_linear_residual_predict(
        train_lonlat, planar_residual, target_lonlat, weights;
        slope_ridge=1e-8,
    )
    @test planar_prediction[1, 1] ≈ 3.0 atol=1e-10

    insufficient_weights = reshape([1.0, 1.0, 0.0, 0.0], 4, 1)
    insufficient = local_linear_residual_predict(
        train_lonlat, constant_residual[:, 1:1], target_lonlat, insufficient_weights;
        slope_ridge=1e-6,
    )
    @test isnan(insufficient[1, 1])

    repeated_lonlat = repeat(target_lonlat, 3, 1)
    degenerate = local_linear_residual_predict(
        repeated_lonlat, fill(1.0, 3, 1), target_lonlat, ones(3, 1);
        slope_ridge=1e-6,
    )
    @test isnan(degenerate[1, 1])

    distances = pairwise_haversine_km(train_lonlat)
    loocv_weights = gw_weight(
        make_loocv_dist(distances), 200.0; kernel=GAUSSIAN, adaptive=false,
    )
    @test all(loocv_weights[i, i] == 0.0 for i in axes(loocv_weights, 1))

    sat_train = fill(1.0, 4, 1)
    obs_train = fill(3.0, 4, 1)
    corrected_train, predicted_train_residual, _ = bias_correct_stgwr(
        train_lonlat, obs_train, sat_train, distances;
        kernel=GAUSSIAN, adaptive=false, bw=200.0,
        use_loocv=false, slope_ridge=1e-6,
    )
    @test all(predicted_train_residual .≈ 2.0)
    @test all(corrected_train .≈ 3.0)

    predicted_residual = local_linear_residual_predict(
        train_lonlat, fill(-2.0, 4, 1), target_lonlat, weights;
        slope_ridge=1e-6,
    )
    corrected = fill(0.5, 1, 1) .+ predicted_residual
    @test corrected[1, 1] ≈ -1.5 atol=1e-12
    negative = negative_output_stats(corrected)
    @test negative.negative_n == 1
    @test negative.negative_fraction == 1.0
    @test negative.min_corrected ≈ -1.5

    corrected_from_obs_minus_sat = fill(0.5, 1, 1) .+ local_linear_residual_predict(
        train_lonlat, fill(2.0, 4, 1), target_lonlat, weights;
        slope_ridge=1e-6,
    )
    @test corrected_from_obs_minus_sat[1, 1] ≈ 2.5 atol=1e-12
end

@testset "MGER scan coverage threshold" begin
    lonlat = [
        110.00 30.00
        110.01 30.00
        110.00 30.01
        110.01 30.01
        111.00 31.00
    ]
    observations = fill(2.0, 5, 2)
    satellite = fill(1.0, 5, 2)
    distances = pairwise_haversine_km(lonlat)
    common_args = (
        kernels=[BOXCAR],
        bw_adaptive=Float64[],
        bw_fixed_km=[5.0],
        slope_ridge_candidates=[1e-6],
        use_loocv=true,
    )

    scan, best = scan_params(
        lonlat, observations, satellite, distances;
        common_args...,
        min_scan_coverage=0.75,
    )
    @test nrow(scan) == 1
    @test best.coverage == 0.8
    @test_throws ArgumentError scan_params(
        lonlat, observations, satellite, distances;
        common_args...,
        min_scan_coverage=0.95,
    )
    @test_throws ArgumentError scan_params(
        lonlat, observations, satellite, distances;
        common_args...,
        min_scan_coverage=1.1,
    )
end

@testset "MGER five kernels" begin
    @test kernel_name.(TEST_KERNELS) == TEST_KERNEL_NAMES
    @test_throws ArgumentError kernel_name(-1)
    @test_throws ArgumentError kernel_name(5)

    distances = reshape([0.0, 2.0], 2, 1)
    gaussian = gw_weight(distances, 1.0; kernel=GAUSSIAN, adaptive=false)
    exponential = gw_weight(distances, 1.0; kernel=EXPONENTIAL, adaptive=false)
    @test gaussian[2] > 0.0
    @test exponential[2] > 0.0
    for kernel in (BISQUARE, TRICUBE, BOXCAR)
        weights = gw_weight(distances, 1.0; kernel=kernel, adaptive=false)
        @test weights[1] == 1.0
        @test weights[2] == 0.0
    end

    mktempdir() do temp_dir
        fixture = make_fixture(temp_dir)
        result_dir = joinpath(temp_dir, "five_kernels")
        cfg = fixture_config(fixture, result_dir)
        summary = run_multikernel_spatial_kfold_pipeline(
            cfg; k=2, seed=20260627, fold_scheme=:random,
        )

        @test nrow(summary) == 15
        @test sort(unique(summary.kernel)) == TEST_KERNELS
        @test sort(unique(summary.product)) == ["FY4B", "GPM", "GSMaP"]

        status = CSV.read(joinpath(result_dir, "kernel_run_status.csv"), DataFrame)
        @test nrow(status) == 5
        @test all(status.status .== "success")
        @test status.kernel_name == TEST_KERNEL_NAMES

        pooled = CSV.read(joinpath(result_dir, "summary_five_kernels_pooled.csv"), DataFrame)
        folds = CSV.read(joinpath(result_dir, "summary_five_kernels_folds.csv"), DataFrame)
        fold_stats = CSV.read(joinpath(result_dir, "summary_five_kernels_fold_stats.csv"), DataFrame)
        @test nrow(pooled) == 15
        @test nrow(folds) == 30
        @test nrow(fold_stats) == 15
        @test :slope_ridge_by_fold in propertynames(pooled)
        @test :slope_ridge in propertynames(folds)
        @test all(col in propertynames(pooled) for col in (
            :negative_n, :negative_fraction, :min_corrected,
        ))

        reference_split = DataFrame()
        for (kernel, name) in zip(TEST_KERNELS, TEST_KERNEL_NAMES)
            kernel_dir = joinpath(result_dir, name)
            split = CSV.read(joinpath(kernel_dir, "split_common_spatial5fold.csv"), DataFrame)
            if isempty(reference_split)
                reference_split = split
            else
                @test split == reference_split
            end

            for product in ("FY4B", "GPM", "GSMaP")
                corrected = CSV.read(
                    joinpath(kernel_dir, "corr_$(product)_spatial5fold_val.csv"), DataFrame,
                )
                @test nrow(corrected) == length(fixture.times)
                @test ncol(corrected) - 1 == length(fixture.station_ids)
                for fold in 1:2
                    fold_dir = joinpath(kernel_dir, "fold_$(fold)")
                    scan = CSV.read(joinpath(fold_dir, "scan_$(product)_spatialcv.csv"), DataFrame)
                    failures = CSV.read(
                        joinpath(fold_dir, "scan_$(product)_spatialcv_failures.csv"), DataFrame,
                    )
                    attempts = vcat(
                        select(scan, :kernel, :adaptive, :bw, :slope_ridge),
                        select(failures, :kernel, :adaptive, :bw, :slope_ridge),
                    )
                    @test nrow(attempts) == 16
                    @test all(attempts.kernel .== kernel)
                    @test Set(attempts.adaptive) == Set([false, true])
                    @test Set(attempts.slope_ridge) == Set(cfg.slope_ridge_candidates)
                end
            end
        end

        scope = CSV.read(
            joinpath(result_dir, "gaussian", "validation_scope.csv"), DataFrame,
        )
        scope_map = Dict(scope.key .=> scope.value)
        @test scope_map["residual_definition"] == "R = P_obs - P_sat; P_corr = P_sat + R_hat"
        @test scope_map["negative_precipitation_policy"] == "retain raw negative corrected values without clipping"
        @test parse(Float64, scope_map["min_scan_coverage"]) == cfg.min_scan_coverage
        # This run passes `fold_scheme=:random`, which interleaves held-out stations with training
        # ones, so the scope must not advertise spatial generalisation. The artefact filenames
        # still say "spatial5fold" for continuity with existing output trees; the recorded claim is
        # what has to be true.
        @test scope_map["fold_scheme"] == "random"
        @test occursin("NOT supported", scope_map["supported_claim"])
        @test !occursin("spatial bias correction", scope_map["supported_claim"])

        duplicate_cfg = fixture_config(fixture, joinpath(temp_dir, "duplicate"); kernels=[GAUSSIAN, GAUSSIAN])
        invalid_cfg = fixture_config(fixture, joinpath(temp_dir, "invalid"); kernels=[9])
        invalid_coverage_cfg = MGERConfig(
            station_meta_path=fixture.station_path,
            obs_hourly_wide_path=fixture.obs_path,
            sat_paths=copy(fixture.sat_paths),
            outdir=joinpath(temp_dir, "invalid_coverage"),
            kernels=[GAUSSIAN],
            min_scan_coverage=1.1,
        )
        @test_throws ArgumentError run_multikernel_spatial_kfold_pipeline(
            duplicate_cfg; k=2, seed=1,
        )
        @test_throws ArgumentError run_multikernel_spatial_kfold_pipeline(
            invalid_cfg; k=2, seed=1,
        )
        @test_throws ArgumentError run_multikernel_spatial_kfold_pipeline(
            invalid_coverage_cfg; k=2, seed=1,
        )
    end
end

@testset "_kernel_selection_stability" begin
    folds = DataFrame(
        product=["FY4B", "FY4B", "FY4B", "GPM", "GPM", "GPM"],
        fold=[1, 2, 3, 1, 2, 3],
        kernel=[GAUSSIAN, GAUSSIAN, BISQUARE, BOXCAR, BOXCAR, BOXCAR],
    )
    stability = _kernel_selection_stability(folds)
    @test sort(unique(stability.product)) == ["FY4B", "GPM"]
    @test all(combine(groupby(stability, :product), :win_count => sum => :total).total .== 3)

    fy4b = sort(filter(:product => ==("FY4B"), stability), :kernel)
    gaussian_row = only(filter(:kernel => ==(GAUSSIAN), fy4b))
    bisquare_row = only(filter(:kernel => ==(BISQUARE), fy4b))
    @test gaussian_row.win_count == 2
    @test bisquare_row.win_count == 1
    @test gaussian_row.selection_frequency ≈ 2 / 3
    # Sorted by win_count descending within product: the 2/3 majority kernel comes first.
    @test first(filter(:product => ==("FY4B"), stability)).kernel == GAUSSIAN

    gpm = filter(:product => ==("GPM"), stability)
    @test nrow(gpm) == 1
    @test only(gpm).win_count == 3
    @test only(gpm).selection_frequency == 1.0

    @test isempty(_kernel_selection_stability(DataFrame()))
end

@testset "Nested kernel selection closes the argmax-on-held-out-scores leak" begin
    mktempdir() do temp_dir
        fixture = make_fixture(temp_dir)
        result_dir = joinpath(temp_dir, "nested_kernels")
        cfg = fixture_config(fixture, result_dir)
        pooled = run_nested_kernel_spatial_kfold_pipeline(
            cfg; k=2, seed=20260627, fold_scheme=:random,
        )

        @test nrow(pooled) == 3
        @test sort(unique(pooled.product)) == ["FY4B", "GPM", "GSMaP"]

        stability = CSV.read(joinpath(result_dir, "kernel_selection_stability.csv"), DataFrame)
        @test sort(unique(stability.product)) == ["FY4B", "GPM", "GSMaP"]
        totals = combine(groupby(stability, :product), :win_count => sum => :total)
        @test all(==(2), totals.total) # k=2 folds

        folds = CSV.read(joinpath(result_dir, "summary_three_products_folds.csv"), DataFrame)
        @test nrow(folds) == 3 * 2
        # The per-fold winner recorded in the pooled summary must be one of the folds' actual
        # winners for that product (sanity link between the two output tables).
        for product in ("FY4B", "GPM", "GSMaP")
            product_folds = filter(:product => ==(product), folds)
            @test Set(product_folds.kernel) ⊆ Set(TEST_KERNELS)
        end

        # The core leak-closing property: within each training fold, the scan actually
        # considered every configured kernel (unlike the single-kernel-forced path, where
        # the existing "MGER five kernels" test above asserts `all(attempts.kernel .== kernel)`
        # for exactly this file).
        for product in ("FY4B", "GPM", "GSMaP"), fold in 1:2
            fold_dir = joinpath(result_dir, "fold_$(fold)")
            scan = CSV.read(joinpath(fold_dir, "scan_$(product)_spatialcv.csv"), DataFrame)
            failures = CSV.read(
                joinpath(fold_dir, "scan_$(product)_spatialcv_failures.csv"), DataFrame,
            )
            attempts = vcat(select(scan, :kernel), select(failures, :kernel))
            @test Set(attempts.kernel) == Set(TEST_KERNELS)
            @test nrow(attempts) == 16 * length(TEST_KERNELS)
        end

        single_kernel_cfg = fixture_config(
            fixture, joinpath(temp_dir, "single"); kernels=[GAUSSIAN],
        )
        @test_throws ArgumentError run_nested_kernel_spatial_kfold_pipeline(
            single_kernel_cfg; k=2, seed=1,
        )
    end
end
