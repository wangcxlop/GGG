#!/usr/bin/env julia

using Test
using CSV, DataFrames, Dates

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
            bw_adaptive=[4.0],
            bw_fixed_km=[100.0],
            expected_common_time_count=length(times),
        )
        cfg = InterpolationBenchmarkConfig(
            mger=mger,
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
        @test all(result.status.status .== "success")
        @test Set(result.metrics.method) == Set(BENCHMARK_METHODS)
        @test isfile(joinpath(outdir, "metrics_stratified.csv"))
        @test isfile(joinpath(outdir, "parameter_scan.csv"))
        @test isfile(joinpath(outdir, "paired_comparisons.csv"))
        for product in ("fy4b", "gpm", "gsmap"), method in BENCHMARK_METHODS
            path = joinpath(outdir, "balanced_spatial", product, "oof_$(method).csv")
            @test isfile(path)
            output = CSV.read(path, DataFrame)
            values = Matrix{Float64}(output[:, Not(:time)])
            @test all(value -> isnan(value) || value >= 0, values)
        end
    end
end
