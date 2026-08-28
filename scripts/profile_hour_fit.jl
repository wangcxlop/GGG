#!/usr/bin/env julia

"""
Cost of one hour-fit, and how the benchmark's dominant loop scales with threads.

`dynamic_covariate_predict` is where an interpolation-benchmark run spends most of its time: one
model fit per hour, per inner selection group, per tuning candidate. Each of those fits allocates
a couple of n×n matrices, and past a handful of threads that allocation contends worse than the
extra cores help — so the loop has a wall-clock optimum well below one thread per core. Where that
optimum sits is a property of the machine, not of the code, so it has to be measured.

Run it once per candidate thread count and read off the loop line:

    julia -t 4  --project=. scripts/profile_hour_fit.jl full
    julia -t 8  --project=. scripts/profile_hour_fit.jl full
    julia -t 24 --project=. scripts/profile_hour_fit.jl full

`smoke` uses the June-2022 grid and loads much faster; `full` is what the reported runs use and is
the one to tune against. Writes nothing.
"""

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using DataFrames, Dates, LinearAlgebra, Printf, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
using Main: DEMTerrainExperiment, JointCovariateModels, JointVariableSelection
using Main.JointCovariateModels

include(joinpath(ROOT, "scripts", "run_interpolation_benchmark.jl"))

"""Best of `repeats` timed runs of `f`, after one warm-up call to keep compilation out of it."""
function probe(label, repeats, f)
    f()
    best = nothing
    for _ in 1:repeats
        result = @timed f()
        (best === nothing || result.time < best.time) && (best = result)
    end
    @printf("%-46s %9.3f ms %8.2f MiB\n", label, 1000 * best.time, best.bytes / 2^20)
    return best
end

"""One fold's `JointFoldContext`, built the way `_run_benchmark_fold!` builds it."""
function build_probe_context(cfg)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    joint_inputs = load_joint_benchmark_inputs(
        something(cfg.joint_covariates), products, ids, product_data[first(products)].times)
    product = first(products)
    data = product_data[product]
    folds = benchmark_folds(:balanced_spatial, ids, lonlat;
        k=cfg.k, seed=cfg.seed, center_init=cfg.fold_center_init, rotation=0)
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    train_idx = [id_map[id] for id in reduce(vcat, (folds[i] for i in 2:cfg.k))]
    val_idx = [id_map[id] for id in folds[1]]
    ndvi = joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned
    selection = JointVariableSelection.select_joint_covariates(
        product, data.Y_obs, data.Y_sat, joint_inputs.terrain, joint_inputs.era5.values,
        ndvi, train_idx, ids, data.times, lonlat, something(cfg.joint_selection),
        cfg.seed, cfg.seed + 20260815)
    context = build_joint_fold_context(
        product, selection.role_map, train_idx, val_idx, lonlat, data.Y_obs, data.Y_sat,
        joint_inputs.terrain, joint_inputs.era5.values, ndvi, something(cfg.joint_covariates))
    return (; context, data, train_idx)
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    cfg = benchmark_config(mode; nested_covariates=true)
    (; context, data, train_idx) = build_probe_context(cfg)

    residuals = data.Y_obs[train_idx, :] .- data.Y_sat[train_idx, :]
    kernel = _kernel_function(BISQUARE)
    bandwidth = Float64(first(context.bandwidth_candidates))
    tuning = first(_tuning_time_sample(
        data.Y_obs[train_idx, :], cfg.tuning_max_times, cfg.tuning_time_weighting))

    # The hour with the most usable stations: the most expensive one, and so the right unit cost.
    time = argmax([count(isfinite, @view residuals[:, t]) for t in axes(residuals, 2)])
    design = JointCovariateModels._design_at(context, "mixed_gwr", time; target=false)
    y = Vector{Float64}(residuals[:, time])
    valid = isfinite.(y)
    for X in design.local_groups
        valid .&= vec(all(isfinite, X; dims=2))
    end
    Xlocal = Matrix{Float64}(design.local_groups[1][valid, :])
    train_lonlat = context.train_lonlat[valid, :]
    distances = DEMTerrainExperiment._haversine_matrix(train_lonlat, train_lonlat)
    response = reshape(y[valid], :, 1)
    adjusted = min(bandwidth, count(valid) - 1)
    no_global = zeros(Float64, count(valid), 0)

    @printf("threads=%d  blas=%d  mode=%s\n", Threads.nthreads(), BLAS.get_num_threads(), mode)
    @printf("train=%d  valid=%d  local cols=%d  tuning hours=%d  variables=%s\n",
        size(context.train_lonlat, 1), count(valid), size(Xlocal, 2),
        length(tuning), context.variables)
    println(repeat("-", 68))

    probe("_haversine_matrix (n x n)", 100, () ->
        DEMTerrainExperiment._haversine_matrix(train_lonlat, train_lonlat))
    probe("_local_hat (one group)", 100, () -> DEMTerrainExperiment._local_hat(
        Xlocal, Xlocal, distances, adjusted, kernel; adaptive=true, ridge=1e-8))
    probe("one hour-fit (mixed_gwr_predict)", 100, () ->
        DEMTerrainExperiment.mixed_gwr_predict(
            Xlocal, no_global, response, train_lonlat, Xlocal, no_global, train_lonlat,
            adjusted, kernel; adaptive=true, ridge=context.config.ridge,
            tolerance=context.config.tolerance,
            max_iterations=context.config.max_iterations, exclude_self=true,
            train_distances=distances, target_distances=distances))

    println(repeat("-", 68))
    loop = probe("dynamic_covariate_predict ($(length(tuning)) hours)", 5, () ->
        dynamic_covariate_predict(context, residuals, "mixed_gwr", [bandwidth], kernel;
            time_indices=tuning, leave_one_out=true))
    @printf("%-46s %9.1f %%\n", "  of which GC", 100 * loop.gctime / loop.time)
    println(repeat("-", 68))
    println("Compare the loop line across thread counts; lowest wins. See CLAUDE.md.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
