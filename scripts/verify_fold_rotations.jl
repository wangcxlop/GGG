#!/usr/bin/env julia

#=
Do the seed-free Hilbert rotations behave like the seeded partitions they replace?

`balanced_spatial` folds used to be initialised from an RNG: `_initial_spatial_centers` drew the
first center with `rand` and the rest by k-means++ weighting, so every reported number depended on
which seed happened to draw the partition. `:hilbert` replaces that draw with a declared rotation
of the station frame - repeat `i` uses rotation `i-1`, and rotation 0, the canonical partition, has
no free parameter at all.

Swapping the initialiser is only safe if two things hold, and neither is safe to assume:

  geometry     Fold *shape* is load-bearing. `tuning_geometry = :inner_spatial` is justified by the
               outer folds' nearest-training-gauge distances (p10/median/p90 = 8.6/23.9/48.9 km);
               if the new partition is less compact, selection and reporting stop being the same
               estimand and the justification in `selection_folds` no longer applies.

  distinctness Repeated cross-validation only measures split sensitivity if the partitions actually
               differ. The `:farthest` initialisation failed exactly here - 20 seeds collapsed into
               10 partitions differing by 1-5 stations - and capacity-constrained Lloyd can wash out
               the initialisation the same way.

This script measures both on the real network and prints the seeded incumbent beside them, so the
comparison is a measurement rather than a claim. It fits no models: it needs station coordinates
only, and runs in seconds.
=#

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "scripts", "run_interpolation_benchmark.jl"))

using CSV, DataFrames, Statistics

const OUTDIR = joinpath(ROOT, "output", "fold_rotation_verification")
const ROTATIONS = 0:19
# The geometry the `:inner_spatial` tuning argument was calibrated on.
const REFERENCE_MEDIAN_KM = 23.9

"""Nearest training-gauge distance for every station, pooled over the folds it is held out of.

Mirrors what `run_interpolation_benchmark` records per fold at its `haversine_distance_matrix`
call, so the numbers here are the same quantity the benchmark reports."""
function nearest_train_distances(ids, lonlat, folds, k)
    position = Dict(id => index for (index, id) in enumerate(ids))
    distances = Float64[]
    for fold in 1:k
        val = [position[id] for id in folds[fold]]
        train = [position[id] for id in reduce(vcat, (folds[j] for j in 1:k if j != fold))]
        matrix = haversine_distance_matrix(lonlat[train, :], lonlat[val, :])
        append!(distances, vec(minimum(matrix, dims=1)))
    end
    return distances
end

_signature(folds) = sort([sort(fold) for fold in folds])

function summarise(label, partition_id, ids, lonlat, folds, k)
    distances = nearest_train_distances(ids, lonlat, folds, k)
    return (;
        scheme=label, partition=partition_id,
        min_fold_size=minimum(length.(folds)), max_fold_size=maximum(length.(folds)),
        balanced=maximum(length.(folds)) - minimum(length.(folds)) <= 1,
        nearest_train_p10=round(quantile(distances, 0.1), digits=2),
        nearest_train_median=round(quantile(distances, 0.5), digits=2),
        nearest_train_p90=round(quantile(distances, 0.9), digits=2),
    )
end

function verify()
    mkpath(OUTDIR)
    cfg = benchmark_config(:full)
    k = cfg.k
    # The benchmark folds the stations common to every product, so the fold geometry has to be
    # measured on that same set rather than on everything in the station table.
    _, ids, _ = load_global_common_product_data(cfg.mger)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col,
        lat_col=cfg.mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    println("stations: $(length(ids))  k=$k")

    rows = NamedTuple[]
    rotations = [benchmark_folds(:balanced_spatial, ids, lonlat; k=k, seed=cfg.seed,
        center_init=:hilbert, rotation=r) for r in ROTATIONS]
    for (index, folds) in enumerate(rotations)
        push!(rows, summarise("hilbert", ROTATIONS[index], ids, lonlat, folds, k))
    end
    seeds = [20260627 + 1000 * i for i in 0:(length(ROTATIONS) - 1)]
    incumbent = [benchmark_folds(:balanced_spatial, ids, lonlat; k=k, seed=s,
        center_init=:kmeanspp) for s in seeds]
    for (index, folds) in enumerate(incumbent)
        push!(rows, summarise("kmeanspp", seeds[index], ids, lonlat, folds, k))
    end
    table = DataFrame(rows)
    CSV.write(joinpath(OUTDIR, "fold_rotation_geometry.csv"), table)

    distinct(partitions, n) = length(unique(_signature.(partitions[1:n])))
    canonical = summarise("hilbert", 0, ids, lonlat, rotations[1], k)
    seeded = summarise("kmeanspp", seeds[1], ids, lonlat, incumbent[1], k)

    println()
    println(table)
    println()
    println("canonical rotation 0   p10/median/p90 = $(canonical.nearest_train_p10) / " *
        "$(canonical.nearest_train_median) / $(canonical.nearest_train_p90) km")
    println("seeded incumbent       p10/median/p90 = $(seeded.nearest_train_p10) / " *
        "$(seeded.nearest_train_median) / $(seeded.nearest_train_p90) km")
    println("rotation 0 is the incumbent partition: " *
        string(_signature(rotations[1]) == _signature(incumbent[1])))
    println()
    println("distinct hilbert  partitions: $(distinct(rotations, 10))/10  " *
        "$(distinct(rotations, length(ROTATIONS)))/$(length(ROTATIONS))")
    println("distinct kmeanspp partitions: $(distinct(incumbent, 10))/10  " *
        "$(distinct(incumbent, length(ROTATIONS)))/$(length(ROTATIONS))")

    # The two gates. Geometry is the one that would invalidate `selection_folds`' argument;
    # distinctness only has to match what the seeded scheme already delivered.
    geometry_ok = abs(canonical.nearest_train_median - REFERENCE_MEDIAN_KM) <= 1.0
    distinct_ok = distinct(rotations, 10) >= distinct(incumbent, 10)
    balanced_ok = all(row -> row.balanced, rows)
    println()
    println("GATE geometry (rotation 0 median within 1 km of $REFERENCE_MEDIAN_KM): " *
        (geometry_ok ? "PASS" : "FAIL"))
    println("GATE distinctness (hilbert >= kmeanspp over 10 partitions): " *
        (distinct_ok ? "PASS" : "FAIL"))
    println("GATE balance (every partition within one station): " * (balanced_ok ? "PASS" : "FAIL"))
    println()
    println("wrote $(joinpath(OUTDIR, "fold_rotation_geometry.csv"))")
    return geometry_ok && distinct_ok && balanced_ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    verify() || exit(1)
end
