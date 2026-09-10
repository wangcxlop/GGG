"""Monte Carlo coefficient-surface variability test with a joint aspect statistic."""
function spatial_variability_test(
    terrain::DataFrame, lonlat::Matrix{Float64}, response::Vector{Float64}, groups::Vector{String};
    product::String="product", bandwidth_candidates::Vector{Int}=[30, 50, 80, 120, 160],
    permutations::Int=999, q_threshold::Float64=0.05, seed::Int=20260815,
    ridge::Float64=1e-8,
)
    isempty(groups) && return DataFrame(), DataFrame()
    valid = isfinite.(response)
    count(valid) >= 12 || throw(ArgumentError("at least twelve stations are required for spatial testing"))
    y = response[valid]
    coords = lonlat[valid, :]
    terrain_valid = terrain[valid, :]
    terrain_X, terrain_columns, terrain_column_groups = _terrain_matrix(terrain_valid, groups)
    spatial_X = hcat(coords[:, 1], coords[:, 2], terrain_X)
    spatial_Z, _, _ = _standardize_train(spatial_X)
    X = hcat(ones(size(spatial_Z, 1)), spatial_Z)
    design_groups = vcat(["intercept", "longitude", "latitude"], terrain_column_groups)
    distances = _haversine_matrix(coords, coords)
    best, scan = _select_gwr_bandwidth(X, y, distances, bandwidth_candidates; ridge)
    smoothers = _gwr_smoothers(
        X, distances, best.bandwidth; ridge, exclude_self=false,
    )
    observed_beta = hcat((smoother * y for smoother in smoothers)...)

    rng = MersenneTwister(seed)
    group_indices = [findall(==(group), design_groups) for group in groups]
    observed_statistics = [sum(var(@view observed_beta[:, index]) for index in indices)
        for indices in group_indices]
    # The permutations are drawn up front, from the same `rng` in the same order, so
    # `permutations_drawn[i]` is exactly what iteration `i` used to draw for itself. That leaves
    # the bodies independent, and each one is a full O(n^2 p) smoother rebuild - the dominant cost
    # of the whole test. `hits` is written one row per permutation and reduced afterwards;
    # summing booleans is exact and order-free, so `exceed` is unchanged.
    permutations_drawn = [randperm(rng, length(y)) for _ in 1:permutations]
    # Weights depend on the geometry and bandwidth only, never on the permuted design, so they are
    # built once here instead of once per permutation.
    local_weights = _gwr_local_weights(distances, best.bandwidth, size(X, 2); exclude_self=false)
    # `Matrix{Bool}`, not a `BitMatrix`: a BitArray packs 64 entries into one word, so two threads
    # writing different rows of the same column would read-modify-write the same word and lose
    # updates. One byte per entry makes the concurrent writes independent.
    hits = fill(false, permutations, length(groups))
    Threads.@threads :greedy for permutation_index in 1:permutations
        permutation = permutations_drawn[permutation_index]
        Xpermuted = X[permutation, :]
        ypermuted = y[permutation]
        permuted_smoothers = _gwr_smoothers(
            Xpermuted, distances, best.bandwidth; ridge, exclude_self=false, local_weights,
        )
        for group_index in eachindex(groups)
            statistic = sum(var(permuted_smoothers[index] * ypermuted)
                for index in group_indices[group_index])
            hits[permutation_index, group_index] = statistic >= observed_statistics[group_index]
        end
    end
    exceed = vec(sum(hits, dims=1))
    rows = NamedTuple[]
    pvalues = [(count + 1) / (permutations + 1) for count in exceed]
    for group in groups
        group_index = findfirst(==(group), groups)
        observed = observed_statistics[group_index]
        pvalue = pvalues[group_index]
        push!(rows, (;
            product, variable_group=group, bandwidth=best.bandwidth,
            variability_statistic=observed, pvalue, qvalue=NaN, role="uncertain",
        ))
    end
    adjusted = bh_adjust(pvalues)
    for i in eachindex(rows)
        rows[i] = merge(rows[i], (;
            qvalue=adjusted[i], role=adjusted[i] < q_threshold ? "local" : "global",
        ))
    end
    scan[!, :product] = fill(product, nrow(scan))
    return DataFrame(rows), scan
end
