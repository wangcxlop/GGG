#!/usr/bin/env julia

using Test
using CSV, DataFrames, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

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
                        select(scan, :kernel, :adaptive, :bw),
                        select(failures, :kernel, :adaptive, :bw),
                    )
                    @test nrow(attempts) == 4
                    @test all(attempts.kernel .== kernel)
                    @test Set(attempts.adaptive) == Set([false, true])
                end
            end
        end

        duplicate_cfg = fixture_config(fixture, joinpath(temp_dir, "duplicate"); kernels=[GAUSSIAN, GAUSSIAN])
        invalid_cfg = fixture_config(fixture, joinpath(temp_dir, "invalid"); kernels=[9])
        @test_throws ArgumentError run_multikernel_spatial_kfold_pipeline(
            duplicate_cfg; k=2, seed=1,
        )
        @test_throws ArgumentError run_multikernel_spatial_kfold_pipeline(
            invalid_cfg; k=2, seed=1,
        )
    end
end

