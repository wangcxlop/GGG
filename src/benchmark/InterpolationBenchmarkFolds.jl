"""Draw an index with probability proportional to `weights` (k-means++ sampling)."""
function _sample_weighted(weights::Vector{Float64}, rng::AbstractRNG)
    total = sum(weights)
    (isfinite(total) && total > 0) || return argmax(weights)
    threshold = rand(rng) * total
    cumulative = 0.0
    for index in eachindex(weights)
        cumulative += weights[index]
        cumulative >= threshold && return index
    end
    return lastindex(weights)
end

"""
Initial fold centers.

`:kmeanspp` draws each center with probability proportional to its squared distance to the
nearest existing center, so different seeds give genuinely different partitions. `:farthest`
is the original deterministic farthest-point traversal, where the seed only chose the first
center and therefore produced near-identical partitions across seeds.
"""
function _initial_spatial_centers(
    xy::Matrix{Float64}, k::Int, rng::AbstractRNG; center_init::Symbol=:kmeanspp,
)
    n = size(xy, 1)
    centers = Matrix{Float64}(undef, k, 2)
    first_index = rand(rng, 1:n)
    centers[1, :] = xy[first_index, :]
    nearest2 = fill(Inf, n)
    for center_index in 2:k
        previous = @view centers[center_index - 1, :]
        for station in 1:n
            d2 = sum(abs2, @view(xy[station, :]) .- previous)
            nearest2[station] = min(nearest2[station], d2)
        end
        next_index = center_init === :kmeanspp ? _sample_weighted(nearest2, rng) :
            argmax(nearest2)
        centers[center_index, :] = xy[next_index, :]
    end
    return centers
end

function _capacity_assignment(xy::Matrix{Float64}, centers::Matrix{Float64}, capacities::Vector{Int})
    n, k = size(xy, 1), size(centers, 1)
    distances = Matrix{Float64}(undef, n, k)
    for station in 1:n, cluster in 1:k
        distances[station, cluster] = sum(abs2, @view(xy[station, :]) .- @view(centers[cluster, :]))
    end
    certainty = Vector{Float64}(undef, n)
    for station in 1:n
        ordered = partialsort(@view(distances[station, :]), 1:min(2, k))
        certainty[station] = length(ordered) == 1 ? Inf : ordered[2] - ordered[1]
    end
    order = sortperm(1:n; by=station -> (-certainty[station], station))
    remaining = copy(capacities)
    assignments = zeros(Int, n)
    for station in order
        cluster_order = sortperm(1:k; by=cluster -> (distances[station, cluster], cluster))
        chosen = findfirst(cluster -> remaining[cluster] > 0, cluster_order)
        chosen === nothing && error("capacity assignment failed")
        cluster = cluster_order[chosen]
        assignments[station] = cluster
        remaining[cluster] -= 1
    end
    return assignments
end

"""Bits per axis in the Hilbert grid. 16 puts the grid resolution far below station spacing."""
const HILBERT_ORDER = 16
# Golden-ratio step over a full turn. Successive rotations stay well spread without the repeat
# count having to be known in advance, and rotation 0 leaves the frame unrotated, so the canonical
# partition is the one with no free parameter at all.
#
# The span is a full turn, not the 90 degrees the square grid's symmetry would suggest: the curve
# has fixed start and end corners, so a quarter turn does not map it to itself and the ordering
# genuinely differs. Measured on the 237-station network, a full turn yields 8 distinct partitions
# in the first 10 rotations and 14 in the first 20, against 8 and 12 for a quarter turn.
const HILBERT_ROTATION_STEP = 0.6180339887498949

"""
Hilbert-curve index of integer grid cell `(x, y)` on a `2^order` square grid.

The standard `xy2d` bit-interleave with a per-quadrant rotation. Cells adjacent along the curve are
adjacent in space, which is what makes a contiguous run of the sorted stations a compact block.
"""
function _hilbert_index(x::Int, y::Int, order::Int)
    n = 1 << order
    d = 0
    s = n >> 1
    while s > 0
        rx = (x & s) > 0 ? 1 : 0
        ry = (y & s) > 0 ? 1 : 0
        d += s * s * ((3 * rx) ⊻ ry)
        if ry == 0
            if rx == 1
                x = n - 1 - x
                y = n - 1 - y
            end
            x, y = y, x
        end
        s >>= 1
    end
    return d
end

"""Map `values` onto the `0:2^HILBERT_ORDER-1` integer grid; a degenerate axis collapses to 0."""
function _hilbert_axis(values::Vector{Float64})
    low, high = extrema(values)
    high > low || return zeros(Int, length(values))
    return round.(Int, (values .- low) .* (((1 << HILBERT_ORDER) - 1) / (high - low)))
end

"""
Deterministic initial fold centers, from a Hilbert-curve ordering of the stations.

This is what replaces the random seed. Rotating the cloud by `360° * frac(rotation * φ⁻¹)` before
indexing makes every `rotation` a declared, reproducible partition rather than an RNG draw, and
rotation 0 is the canonical one. Cutting the ordering into `k` equal runs already gives compact and
exactly balanced blocks; the caller's capacity-constrained Lloyd loop then refines them into the
Voronoi-like cells the reported fold geometry was measured on.

Centroids are returned in the unrotated frame, because that is the frame Lloyd iterates in.
"""
function _hilbert_centers(xy::Matrix{Float64}, ids::Vector{String}, k::Int, rotation::Int)
    n = size(xy, 1)
    angle = 2 * pi * mod(rotation * HILBERT_ROTATION_STEP, 1.0)
    cosine, sine = cos(angle), sin(angle)
    grid_x = _hilbert_axis(cosine .* xy[:, 1] .- sine .* xy[:, 2])
    grid_y = _hilbert_axis(sine .* xy[:, 1] .+ cosine .* xy[:, 2])
    index = [_hilbert_index(grid_x[i], grid_y[i], HILBERT_ORDER) for i in 1:n]
    # Station id breaks ties so coincident gauges cannot reorder between runs.
    order = sortperm(1:n; by=station -> (index[station], ids[station]))
    centers = Matrix{Float64}(undef, k, 2)
    base, extra = div(n, k), rem(n, k)
    start = 1
    for cluster in 1:k
        len = base + (cluster <= extra ? 1 : 0)
        centers[cluster, :] = vec(mean(xy[order[start:start + len - 1], :], dims=1))
        start += len
    end
    return centers
end

"""
Reproducible, compact, capacity-balanced two-dimensional spatial folds.

`center_init=:hilbert` is seed-free: `seed` is ignored and `rotation` selects the partition. The
`:kmeanspp` and `:farthest` initialisations keep reading `seed` so earlier runs still reproduce.
"""
function split_stations_balanced_spatial_kfold(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, seed::Int=20260627,
    center_init::Symbol=:kmeanspp, rotation::Int=0,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat, 1) == n || throw(DimensionMismatch("lonlat rows must match station_ids"))
    xy = local_km_coordinates(lonlat)
    capacities = fill(div(n, k), k)
    capacities[1:rem(n, k)] .+= 1
    centers = center_init === :hilbert ? _hilbert_centers(xy, station_ids, k, rotation) :
        _initial_spatial_centers(xy, k, MersenneTwister(seed); center_init)
    assignments = zeros(Int, n)
    for _ in 1:100
        updated = _capacity_assignment(xy, centers, capacities)
        updated == assignments && break
        assignments = updated
        for cluster in 1:k
            idx = findall(==(cluster), assignments)
            centers[cluster, :] = vec(mean(xy[idx, :], dims=1))
        end
    end
    folds = [String[] for _ in 1:k]
    for station in 1:n
        push!(folds[assignments[station]], station_ids[station])
    end
    return folds
end

function benchmark_folds(
    scheme::Symbol, ids::Vector{String}, lonlat::Matrix{Float64}; k::Int, seed::Int,
    center_init::Symbol=:kmeanspp, rotation::Int=0,
)
    if scheme == :balanced_spatial
        return split_stations_balanced_spatial_kfold(
            ids, lonlat; k=k, seed=seed, center_init, rotation)
    elseif scheme == :random
        return split_stations_kfold(ids; k=k, rng=MersenneTwister(seed))
    elseif scheme == :strip
        return split_stations_spatial_block_kfold(ids, lonlat; k=k)
    end
    throw(ArgumentError("unsupported CV scheme: $scheme"))
end

"""
Inner split of one fold's training stations, as positions within that training set.

Hyperparameters have to be chosen by the same kind of prediction the benchmark reports. The
historical criterion was leave-one-out at training stations, which leaves the held-out station
sitting inside its own neighbourhood: median 5.5 km to the nearest remaining gauge, against
23.9 km from an outer `balanced_spatial` fold station to its nearest training gauge. The two are
not the same estimand, and on the joint-covariate path the reported RMSE curve falls
monotonically across the whole bandwidth grid while the leave-one-out curve is U-shaped with a
minimum at 12-30 — over most of the grid they are anti-correlated, and the leave-one-out pick
costs +19% to +101% of reported RMSE.

Splitting the training stations with the same splitter, the same scheme and `k_inner = cfg.k`
reproduces the outer geometry closely: inner-val nearest-train km p10/median/p90 = 8.1/20.4/46.1
under the default `:hilbert` rotations, and 8.1/22.7/45.6 under the `:kmeanspp` seeds it replaced,
against an outer 8.6/23.9/48.9 that both produce identically. `k_inner=3` overshoots to a 32.2 km
median and `k_inner=8` undershoots to 17.5, so matching `cfg.k` is both the knob-free and the
accurate choice, and the gap between the two initialisations is small next to that spread.

Using `scheme` rather than always splitting spatially is deliberate: under `:random` the
reporting geometry (median 6.3 km) already matches leave-one-out, so a random inner split
correctly leaves that scheme's numbers essentially where they were.

Returns `nothing` when the fold is too small to split — holding out a whole group would leave
too few stations to fit the widest design (`mixed_gwr`/`mgwr` need three local columns plus a
global block). The caller falls back to leave-one-out and records that in `benchmark_scope.csv`,
so a degenerate run never silently reports one criterion while having used the other. At the real
station count this never triggers: 189 training stations split five ways leave 151.
"""
const MIN_SELECTION_TRAIN_STATIONS = 10

function selection_folds(
    cfg::InterpolationBenchmarkConfig, scheme::Symbol, train_ids::Vector{String},
    train_lonlat::Matrix{Float64}, fold::Int, repeat_seed::Int; repeat_index::Int=1,
)
    k_inner = cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k
    n_train = length(train_ids)
    k_inner = min(k_inner, n_train)
    k_inner >= 2 || return nothing
    n_train - cld(n_train, k_inner) >= MIN_SELECTION_TRAIN_STATIONS || return nothing
    # Vary by outer fold and repeat so the inner partition is not the same one every time. Not by
    # product: the station geometry is identical across products, and holding the split fixed
    # keeps the per-product comparison clean. Under `:hilbert` the rotation carries that variation
    # instead of the seed; `100 * repeat_index + fold` cannot collide with the outer rotations
    # (`0:n_repeats-1`) at any realistic repeat count, so an inner split is never the outer one.
    groups = benchmark_folds(scheme, train_ids, train_lonlat;
        k=k_inner, seed=repeat_seed + 7919 * fold, center_init=cfg.fold_center_init,
        rotation=100 * repeat_index + fold)
    position = Dict(id => index for (index, id) in enumerate(train_ids))
    return [[position[id] for id in group] for group in groups if !isempty(group)]
end

"""
Assemble out-of-fold predictions over the training stations.

`predict(inner_train_rows, inner_val_rows)` returns a `length(inner_val_rows) × n_time` matrix.
The result has the same shape leave-one-out produced, so `_candidate_metrics` and every scan row
downstream are unchanged.
"""
function _selection_oof(
    groups::Vector{Vector{Int}}, n_train::Int, n_time::Int, predict,
)
    out_of_fold = fill(NaN, n_train, n_time)
    for group in groups
        inner_train = setdiff(1:n_train, group)
        isempty(inner_train) && continue
        out_of_fold[group, :] = predict(inner_train, group)
    end
    return out_of_fold
end

function _write_split(path::String, ids::Vector{String}, folds::Vector{Vector{String}}, scheme::Symbol)
    fold_map = Dict(id => fold for (fold, fold_ids) in enumerate(folds) for id in fold_ids)
    CSV.write(path, DataFrame(
        station_id=ids, fold=[fold_map[id] for id in ids], scheme=fill(string(scheme), length(ids)),
    ))
end
