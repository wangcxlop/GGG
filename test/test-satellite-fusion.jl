#!/usr/bin/env julia

# The fused anchor and the banded blending weight.
#
# Kept out of `test-interpolation-benchmark.jl` only because that file is already 2000 lines; it
# uses the same fixtures and the same `load_pipeline` entry point, and runs after it in
# `runtests.jl` so the pipeline is already compiled.

using Test
using CSV, DataFrames, Dates

const FUSION_TEST_ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
include(joinpath(FUSION_TEST_ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")

if !isdefined(@__MODULE__, :_write_benchmark_wide)
    function _write_benchmark_wide(path, times, ids, values)
        df = DataFrame(time=times)
        for (index, id) in enumerate(ids)
            df[!, Symbol(id)] = values[index, :]
        end
        CSV.write(path, df)
    end
end

@testset "Satellite fusion" begin
    n_station, n_time = 12, 9
    sources = [
        [0.5 + 0.1 * station + 0.2 * time for station in 1:n_station, time in 1:n_time],
        [1.0 + 0.05 * station * time for station in 1:n_station, time in 1:n_time],
        [2.0 - 0.03 * station + 0.11 * time for station in 1:n_station, time in 1:n_time],
    ]
    truth = [0.25, 0.4, 0.3, 0.2]
    y_obs = truth[1] .+ truth[2] .* sources[1] .+ truth[3] .* sources[2] .+ truth[4] .* sources[3]
    all_rows = collect(1:n_station)

    @testset "recovers the coefficients it was generated from" begin
        coefficients, used = fit_satellite_fusion(y_obs, sources, all_rows)
        @test used == n_station * n_time
        @test coefficients ≈ truth
        @test apply_satellite_fusion(sources, coefficients) ≈ y_obs
    end

    @testset "a single-product fusion is that product" begin
        # The identity every reading of a fused anchor rests on, and the one
        # `scripts/verify_fused_anchor_bounds.jl` asserts before reporting anything.
        identity = apply_satellite_fusion(sources, [0.0, 1.0, 0.0, 0.0])
        @test identity == sources[1]
    end

    @testset "clips at zero and propagates a missing source" begin
        @test all(==(0.0), apply_satellite_fusion(sources, [-100.0, 0.0, 0.0, 0.0]))
        holed = [copy(source) for source in sources]
        holed[2][3, 4] = NaN
        fused = apply_satellite_fusion(holed, truth)
        @test isnan(fused[3, 4])
        @test count(isnan, fused) == 1
        # A NaN gauge value drops its cell from the fit rather than the whole station.
        gapped = copy(y_obs)
        gapped[1, 1] = NaN
        _, used = fit_satellite_fusion(gapped, sources, all_rows)
        @test used == n_station * n_time - 1
    end

    @testset "variants" begin
        equal, used, fell_back = fusion_coefficients(:mean, y_obs, sources, all_rows)
        @test equal == [0.0, 1 / 3, 1 / 3, 1 / 3]
        @test used == 0
        @test !fell_back
        @test apply_satellite_fusion(sources, equal) ≈ (sources[1] .+ sources[2] .+ sources[3]) ./ 3

        fitted, _, fitted_fell_back = fusion_coefficients(:ols, y_obs, sources, all_rows)
        @test fitted ≈ truth
        @test !fitted_fell_back

        # A singular design falls back rather than throwing, and says so.
        collinear = [sources[1], sources[1], sources[2]]
        _, _, singular_fell_back = fusion_coefficients(:ols, y_obs, collinear, all_rows)
        @test singular_fell_back
        @test_throws ArgumentError fusion_coefficients(:nope, y_obs, sources, all_rows)
    end

    @testset "the fit never reads a station outside station_rows" begin
        # The leak the whole design is arranged to avoid: the benchmark fits these coefficients on
        # a fold's training stations and applies them to its held-out ones, so a held-out gauge
        # must not reach the anchor its own prediction is built on.
        train_rows = collect(1:8)
        held_out = collect(9:n_station)
        clean, _ = fit_satellite_fusion(y_obs, sources, train_rows)
        corrupted_obs = copy(y_obs)
        corrupted_obs[held_out, :] .+= 37.0
        corrupted, _ = fit_satellite_fusion(corrupted_obs, sources, train_rows)
        @test clean == corrupted
        # Sentinel: the corruption is real, and a fit that *does* see those rows moves.
        everything, _ = fit_satellite_fusion(corrupted_obs, sources, collect(1:n_station))
        @test everything != clean
    end
end

@testset "Grouped satellite fusion" begin
    n_station, n_time = 5, 6
    times = [DateTime(2022, 6, 1) + Hour(h) for h in (0, 1, 2, 5, 6, 7)]  # a gap after hour 2
    source = [10.0 * station + time for station in 1:n_station, time in 1:n_time]
    neighbours = [filter(!=(station), [3, 1, 2, 5, 4]) for station in 1:n_station]

    @testset "lags respect grid gaps and neighbour means skip NaN" begin
        holed = copy(source)
        holed[3, 2] = NaN
        features = fusion_features([holed], times, neighbours; k=2)
        @test length(features) == 4
        base, lag, lead, neighbour = features
        @test isequal(base, holed)
        @test lag[1, 2] == holed[1, 1]
        @test lag[1, 4] == holed[1, 4]       # 05:00 has no 04:00: falls back to its own value
        @test lead[1, 3] == holed[1, 3]      # 02:00 has no 03:00
        @test lead[1, 4] == holed[1, 5]
        @test lag[3, 3] == holed[3, 3]       # the lagged hour is NaN: own value
        @test neighbour[1, 1] == (holed[3, 1] + holed[2, 1]) / 2
        @test neighbour[1, 2] == holed[2, 2]  # station 3 is NaN at hour 2 and is skipped
        @test isnan(neighbour[3, 2]) == false
    end

    n_station, n_time = 40, 200
    lookup = [
        [mod(7 * station + 3 * time, 11) / 3 for station in 1:n_station, time in 1:n_time],
        [mod(5 * station + 2 * time, 13) / 4 for station in 1:n_station, time in 1:n_time],
    ]
    groups = [station <= 20 ? 1 : 2 for station in 1:n_station, time in 1:n_time]
    truth = [0.5 1.0; 2.0 0.5; -0.3 1.5]
    y_obs = [truth[1, groups[s, t]] + truth[2, groups[s, t]] * lookup[1][s, t] +
             truth[3, groups[s, t]] * lookup[2][s, t] for s in 1:n_station, t in 1:n_time]

    @testset "recovers each group's coefficients" begin
        coefficients, used, fell_back = fit_grouped_fusion(
            y_obs, lookup, lookup, groups, 2, collect(1:n_station); min_cells=10)
        @test coefficients ≈ truth
        @test used == [20 * n_time, 20 * n_time]
        @test !any(fell_back)
        @test apply_grouped_fusion(lookup, lookup, coefficients, groups) ≈ max.(y_obs, 0.0)
    end

    @testset "a sparse group takes the pooled fit" begin
        coefficients, _, fell_back = fit_grouped_fusion(
            y_obs, lookup, lookup, groups, 3, collect(1:n_station); min_cells=10)
        @test fell_back == [false, false, true]
        pooled, _, _ = fit_grouped_fusion(
            y_obs, lookup, lookup, ones(Int, n_station, n_time), 1, collect(1:n_station))
        @test coefficients[:, 3] ≈ pooled[:, 1]
    end

    @testset "the fit never reads a station outside station_rows" begin
        train_rows = collect(1:30)
        clean, _, _ = fit_grouped_fusion(y_obs, lookup, lookup, groups, 2, train_rows)
        corrupted = copy(y_obs)
        corrupted[31:end, :] .+= 37.0
        @test fit_grouped_fusion(corrupted, lookup, lookup, groups, 2, train_rows)[1] == clean
    end

    @testset "a NaN source drops the cell from the fit and the anchor" begin
        holed = [copy(matrix) for matrix in lookup]
        holed[2][4, 7] = NaN
        _, used, _ = fit_grouped_fusion(y_obs, holed, lookup, groups, 2, collect(1:n_station))
        @test sum(used) == n_station * n_time - 1
        fused = apply_grouped_fusion(holed, lookup, truth, groups)
        @test isnan(fused[4, 7]) && count(isnan, fused) == 1
    end
end

@testset "Banded blending weights" begin
    threshold = 0.1
    y_sat = [0.0 0.5 3.0; 9.0 0.0 0.2]
    sources = [
        [0.0 0.5 3.0; 9.0 0.0 0.2],
        [0.0 0.0 4.0; 12.0 0.0 0.0],
        [0.0 0.0 0.0; 20.0 0.0 0.05],
    ]
    prediction = [1.0 2.0 3.0; 4.0 5.0 6.0]
    fallback = [0.0 0.0 0.0; 0.0 0.0 0.0]

    @testset "the constant axis is the shipped intervention set" begin
        band = blend_band_matrix(:constant, y_sat, sources, threshold)
        @test band == [0 1 1; 1 0 1]
        @test blend_band_count(:constant) == 1
    end

    @testset "agreement crosses with the envelope" begin
        band = blend_band_matrix(:agreement_envelope, y_sat, sources, threshold)
        # (1,1): nothing wet -> untouched. (1,2): one product wet, max 0.5 -> light -> band 1.
        # (1,3): two wet, max 4.0 -> moderate -> (2-1)*3+2 = 5.
        # (2,1): three wet, max 20 -> heavy -> (3-1)*3+3 = 9.
        # (2,2): nothing wet -> untouched. (2,3): two wet (0.2, 0.05 is below), max 0.2 -> band 1.
        @test band == [0 1 5; 9 0 1]
        @test blend_band_count(:agreement_envelope) == 9
    end

    @testset "a zero weight is the identity, and band 0 is never touched" begin
        band = blend_band_matrix(:agreement_envelope, y_sat, sources, threshold)
        @test banded_blend_prediction(band, prediction, fallback, zeros(9)) == prediction
        ones_everywhere = banded_blend_prediction(band, prediction, fallback, ones(9))
        @test ones_everywhere[1, 1] == prediction[1, 1]
        @test ones_everywhere[2, 2] == prediction[2, 2]
        @test ones_everywhere[1, 2] == fallback[1, 2]
    end

    @testset "a constant weight reproduces the scalar blend" begin
        band = blend_band_matrix(:constant, y_sat, sources, threshold)
        for lambda in (0.0, 0.3, 1.0)
            @test banded_blend_prediction(band, prediction, fallback, [lambda]) ==
                satellite_wet_blend_prediction(y_sat, prediction, fallback, lambda, threshold)
        end
    end

    @testset "NaN rules match the scalar path" begin
        band = blend_band_matrix(:constant, y_sat, sources, threshold)
        holed_fallback = copy(fallback)
        holed_fallback[1, 2] = NaN
        # A NaN fallback leaves the cell unblended rather than poisoning it, which is what keeps
        # the scored population from moving with the weights.
        @test banded_blend_prediction(band, prediction, holed_fallback, [1.0])[1, 2] ==
            prediction[1, 2]
        holed_prediction = copy(prediction)
        holed_prediction[2, 3] = NaN
        @test isnan(banded_blend_prediction(band, holed_prediction, fallback, [1.0])[2, 3])
    end

    @testset "each band's weight is minimised on its own cells" begin
        # Two bands, opposite answers: the fallback is exact in band 1 and the prediction is exact
        # in band 2, so the only weight vector that minimises both is [1, 0].
        band = [1 1 2 2]
        observations = [5.0 5.0 7.0 7.0]
        predicted = [0.0 0.0 7.0 7.0]
        fell_back_to = [5.0 5.0 0.0 0.0]
        anchor = [1.0 1.0 1.0 1.0]
        choice = select_blend_lambdas(
            observations, anchor, band, predicted, fell_back_to, collect(0.0:0.1:1.0), 2,
        )
        @test choice !== nothing
        @test choice.lambdas == [1.0, 0.0]
        @test choice.inner_RMSE ≈ 0.0 atol = 1e-12
        @test choice.unblended_RMSE > choice.inner_RMSE
        @test choice.wet_cells == 4

        # A band with no scorable cell keeps weight 0, which is "leave it alone".
        sparse_band = [1 1 0 0]
        sparse = select_blend_lambdas(
            observations, anchor, sparse_band, predicted, fell_back_to,
            collect(0.0:0.1:1.0), 2,
        )
        @test sparse.lambdas == [1.0, 0.0]
        @test sparse.wet_cells == 2

        # Nothing scorable at all is a skipped fold, not a default weight.
        @test select_blend_lambdas(
            fill(NaN, 1, 4), anchor, band, predicted, fell_back_to, collect(0.0:0.1:1.0), 2,
        ) === nothing
    end
end

@testset "Fused anchor products run end to end" begin
    mktempdir() do temp_dir
        n_station, k, seed = 24, 3, 5
        ids = string.(8001:(8000 + n_station))
        lon = [110.0 + 0.2 * mod(index - 1, 6) for index in 1:n_station]
        lat = [30.0 + 0.2 * div(index - 1, 6) for index in 1:n_station]
        CSV.write(joinpath(temp_dir, "stations.csv"),
            DataFrame(station_id=ids, lon=lon, lat=lat))
        times = collect(DateTime(2022, 6, 1):Hour(1):DateTime(2022, 6, 1, 5))

        obs = [max(0.0, 0.4 + 0.05 * station + 0.3 * sin(time / 2))
               for station in 1:n_station, time in eachindex(times)]
        # Each product carries its own structure as well as the shared signal. Three affine
        # functions of one field are perfectly collinear, which makes the fusion's normal
        # equations singular and sends every fold down the equal-weight fallback - a fixture that
        # would test the fallback and call it a fit.
        wobble(a, b) = [mod(a * station + b * time, 5) / 20
                        for station in 1:n_station, time in eachindex(times)]
        satellites = Dict(
            "FY4B" => max.(obs .* 0.6 .+ 0.4 .+ wobble(3, 1), 0.0),
            "GPM" => max.(obs .* 1.3 .- 0.1 .+ wobble(1, 2), 0.0),
            "GSMaP" => max.(obs .* 0.9 .+ 0.2 .+ wobble(2, 3), 0.0),
        )
        sat_paths = Dict{String,String}()
        for (product, values) in satellites
            path = joinpath(temp_dir, "$(lowercase(product)).csv")
            _write_benchmark_wide(path, times, ids, values)
            sat_paths[product] = path
        end
        obs_path = joinpath(temp_dir, "obs.csv")
        _write_benchmark_wide(obs_path, times, ids, obs)

        outdir = joinpath(temp_dir, "run")
        cfg = InterpolationBenchmarkConfig(
            mger=MGERConfig(
                station_meta_path=joinpath(temp_dir, "stations.csv"),
                obs_hourly_wide_path=obs_path, sat_paths=sat_paths, outdir=outdir,
                kernels=[GAUSSIAN], bw_adaptive=[8.0], bw_fixed_km=[40.0],
                expected_common_time_count=length(times),
            ),
            k=k, seed=seed, cv_schemes=[:balanced_spatial],
            idw_powers=[2.0], neighbor_candidates=Union{Nothing,Int}[8],
            tps_smooth_candidates=[1e-2], bootstrap_reps=0,
            fused_anchor_variants=[:ols, :mean, :ols_lagnbr],
        )
        result = run_interpolation_benchmark(cfg)

        @testset "every derived product is scored" begin
            products = unique(result.metrics.product)
            @test "MERGED_OLS" in products
            @test "MERGED_MEAN" in products
            @test "MERGED_OLS_LAGNBR" in products
            @test issubset(["FY4B", "GPM", "GSMaP"], products)
            @test isdir(joinpath(outdir, "balanced_spatial", "merged_ols"))
        end

        @testset "the coefficients are recorded per fold and per variant" begin
            path = joinpath(outdir, "fused_anchor_selection.csv")
            @test isfile(path)
            table = CSV.read(path, DataFrame)
            @test nrow(table) == 2 * k
            @test sort(unique(table.variant)) == ["mean", "ols"]
            @test issubset([:beta_fy4b, :beta_gpm, :beta_gsmap], propertynames(table))
            @test !any(table.fell_back)
            mean_rows = filter(row -> row.variant == "mean", table)
            @test all(row -> row.intercept == 0.0, eachrow(mean_rows))
            @test all(row -> row.beta_gpm ≈ 1 / 3, eachrow(mean_rows))
        end

        @testset "a fold's anchor is the fusion fitted without that fold" begin
            # The end-to-end form of the leakage test above: what lands in `oof_raw.csv` for a
            # held-out station has to be the anchor refitted here from the other folds alone.
            # The run's own recorded partition, not one rebuilt here: rebuilding it would test
            # this file's ability to reproduce `benchmark_folds`' defaults rather than the anchor.
            split = CSV.read(joinpath(outdir, "balanced_spatial", "split_common.csv"),
                DataFrame; types=Dict(:station_id => String))
            position = Dict(id => index for (index, id) in enumerate(ids))
            fold_of = Dict(row.station_id => Int(row.fold) for row in eachrow(split))
            source_order = sort(collect(keys(sat_paths)))
            sources = [satellites[product] for product in source_order]
            stored = CSV.read(
                joinpath(outdir, "balanced_spatial", "merged_ols", "oof_raw.csv"), DataFrame)
            for fold in 1:k
                validation = [position[id] for id in ids if fold_of[id] == fold]
                @test !isempty(validation)
                training = setdiff(1:n_station, validation)
                coefficients, _, _ = fusion_coefficients(:ols, obs, sources, training)
                expected = apply_satellite_fusion(sources, coefficients)
                for station in validation, time in eachindex(times)
                    @test stored[time, Symbol(ids[station])] ≈ expected[station, time]
                end
            end
        end

        @testset "the grouped anchor is refitted per fold from training stations" begin
            split = CSV.read(joinpath(outdir, "balanced_spatial", "split_common.csv"),
                DataFrame; types=Dict(:station_id => String))
            position = Dict(id => index for (index, id) in enumerate(ids))
            fold_of = Dict(row.station_id => Int(row.fold) for row in eachrow(split))
            sources = [satellites[product] for product in sort(collect(keys(sat_paths)))]
            lonlat = hcat(lon, lat)
            distance = haversine_distance_matrix(lonlat, lonlat)
            neighbours = [filter(!=(s), sortperm(distance[s, :])) for s in 1:n_station]
            features = fusion_features(sources, times, neighbours)
            groups = blend_band_matrix(:agreement_envelope, first(sources), sources, 0.1) .+ 1
            stored = CSV.read(
                joinpath(outdir, "balanced_spatial", "merged_ols_lagnbr", "oof_raw.csv"), DataFrame)
            for fold in 1:k
                validation = [position[id] for id in ids if fold_of[id] == fold]
                training = setdiff(1:n_station, validation)
                coefficients, _, _ = fit_grouped_fusion(obs, sources, features, groups, 10, training)
                expected = apply_grouped_fusion(sources, features, coefficients, groups)
                for station in validation, time in eachindex(times)
                    @test stored[time, Symbol(ids[station])] ≈ expected[station, time]
                end
            end
            table = CSV.read(joinpath(outdir, "fused_anchor_grouped_coefficients.csv"), DataFrame)
            @test nrow(table) == k * 10 * 13
            @test sort(unique(table.group)) == collect(0:9)
        end
    end
end
