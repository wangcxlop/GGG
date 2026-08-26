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

"""Reproducible, compact, capacity-balanced two-dimensional spatial folds."""
function split_stations_balanced_spatial_kfold(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, seed::Int=20260627,
    center_init::Symbol=:kmeanspp,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat, 1) == n || throw(DimensionMismatch("lonlat rows must match station_ids"))
    xy = local_km_coordinates(lonlat)
    capacities = fill(div(n, k), k)
    capacities[1:rem(n, k)] .+= 1
    centers = _initial_spatial_centers(xy, k, MersenneTwister(seed); center_init)
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
    center_init::Symbol=:kmeanspp,
)
    if scheme == :balanced_spatial
        return split_stations_balanced_spatial_kfold(ids, lonlat; k=k, seed=seed, center_init)
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
reproduces the outer geometry closely (inner-val nearest-train km p10/median/p90 = 8.1/22.7/45.6
against the outer 8.6/23.9/48.9); `k_inner=3` overshoots to a 32.2 km median and `k_inner=8`
undershoots to 17.5, so matching `cfg.k` is both the knob-free and the accurate choice.

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
    train_lonlat::Matrix{Float64}, fold::Int, repeat_seed::Int,
)
    k_inner = cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k
    n_train = length(train_ids)
    k_inner = min(k_inner, n_train)
    k_inner >= 2 || return nothing
    n_train - cld(n_train, k_inner) >= MIN_SELECTION_TRAIN_STATIONS || return nothing
    # Vary by outer fold and repeat so the inner partition is not the same one every time. Not by
    # product: the station geometry is identical across products, and holding the split fixed
    # keeps the per-product comparison clean.
    groups = benchmark_folds(scheme, train_ids, train_lonlat;
        k=k_inner, seed=repeat_seed + 7919 * fold, center_init=cfg.fold_center_init)
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
