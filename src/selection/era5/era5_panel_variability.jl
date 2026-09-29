function _haversine(lon1, lat1, lon2, lat2)
    dlat = deg2rad(lat2 - lat1); dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon / 2)^2
    return 6371.0088 * 2 * asin(sqrt(clamp(a, 0.0, 1.0)))
end

function _distance_matrix(lonlat::Matrix{Float64})
    n = size(lonlat, 1); result = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        d = _haversine(lonlat[i, 1], lonlat[i, 2], lonlat[j, 1], lonlat[j, 2])
        result[i, j] = result[j, i] = d
    end
    return result
end

function _adaptive_weights(distances::AbstractVector{<:Real}, neighbors::Int; exclude::Int=0)
    usable = [i for i in eachindex(distances) if i != exclude && isfinite(distances[i])]
    isempty(usable) && return zeros(length(distances))
    k = min(neighbors, length(usable))
    bandwidth = max(sort(distances[usable])[k], eps(Float64))
    weights = zeros(length(distances))
    for i in usable
        ratio = distances[i] / bandwidth
        ratio < 1 && (weights[i] = (1 - ratio^2)^2)
    end
    return weights
end

function _panel_sufficient_statistics(panel, selected::Vector{Symbol}, lonlat::Matrix{Float64}, cfg)
    eligible = findall(panel.eligible_station)
    isempty(eligible) && return nothing
    coords = lonlat[eligible, :]
    coord_mean = vec(mean(coords, dims=1)); coord_scale = vec(std(coords, dims=1))
    any(x -> !isfinite(x) || x <= 0, coord_scale) && return nothing
    p = 3 + length(selected)
    A = [zeros(p, p) for _ in eligible]
    b = [zeros(p) for _ in eligible]
    rows = [Vector{Vector{Float64}}() for _ in eligible]
    ys = [Float64[] for _ in eligible]
    for (position, i) in enumerate(eligible)
        valid = [t for t in axes(panel.y, 2) if isfinite(panel.y[i, t]) &&
            all(isfinite(panel.x[v][i, t]) for v in selected)]
        isempty(valid) && continue
        weight = 1 / length(valid)
        for t in valid
            x = [1.0, (coords[position, 1] - coord_mean[1]) / coord_scale[1],
                (coords[position, 2] - coord_mean[2]) / coord_scale[2],
                (panel.x[v][i, t] for v in selected)...]
            A[position] .+= weight .* (x * x')
            b[position] .+= weight .* x .* panel.y[i, t]
            push!(rows[position], x); push!(ys[position], panel.y[i, t])
        end
    end
    return (; A, b, rows, ys, coords, eligible, p)
end

"""One `_adaptive_weights` row per target, for reuse across a permutation test.

The weights are a function of `distances` and `neighbors` only — the block permutation reorders
which sufficient statistic each spatial slot contributes, never the weight on that slot. Building
them once removes a sort and two allocations per target from every one of the ~1000 permutations.
"""
_local_weight_rows(distances, neighbors, n::Int) =
    [_adaptive_weights(view(distances, target, :), neighbors) for target in 1:n]

function _local_coefficients(stats, distances, neighbors, cfg;
    block_order=collect(eachindex(stats.A)), weight_rows=nothing)
    n = length(stats.A); coefficients = fill(NaN, n, stats.p)
    rows = weight_rows === nothing ? _local_weight_rows(distances, neighbors, n) : weight_rows
    lhs = Matrix{Float64}(undef, stats.p, stats.p); rhs = Vector{Float64}(undef, stats.p)
    for target in 1:n
        weights = rows[target]
        # Reused across targets rather than reallocated; refilled to exactly what
        # `cfg.ridge .* Matrix{Float64}(I, p, p)` and `zeros(p)` produced.
        fill!(lhs, 0.0); fill!(rhs, 0.0)
        for i in 1:stats.p
            lhs[i, i] = cfg.ridge
        end
        for spatial_index in 1:n
            source = block_order[spatial_index]
            lhs .+= weights[spatial_index] .* stats.A[source]
            rhs .+= weights[spatial_index] .* stats.b[source]
        end
        try
            coefficients[target, :] .= lhs \ rhs
        catch
        end
    end
    return coefficients
end

function _loocv_panel_rmse(stats, distances, neighbors, cfg)
    n = length(stats.A); sse = 0.0; count_station = 0
    for target in 1:n
        weights = _adaptive_weights(view(distances, target, :), neighbors; exclude=target)
        lhs, rhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p), zeros(stats.p)
        for i in 1:n
            lhs .+= weights[i] .* stats.A[i]; rhs .+= weights[i] .* stats.b[i]
        end
        beta = try lhs \ rhs catch; continue end
        isempty(stats.rows[target]) && continue
        station_sse = mean((stats.ys[target][j] - dot(stats.rows[target][j], beta))^2
            for j in eachindex(stats.ys[target]))
        sse += station_sse; count_station += 1
    end
    return count_station == n ? sqrt(sse / n) : Inf
end

function panel_spatial_variability_test(
    panel, selected::Vector{Symbol}, lonlat::Matrix{Float64}, cfg::ERA5SelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed + 1),
)
    bandwidth_scan = DataFrame(neighbors=Int[], loocv_rmse=Float64[], available=Bool[])
    variability = DataFrame(
        variable=String[], statistic=Float64[], pvalue=Float64[], qvalue=Float64[],
        role=String[], neighbors=Union{Missing,Int}[], status=String[],
    )
    isempty(selected) && return (; bandwidth_scan, variability, bandwidth=missing)
    stats = _panel_sufficient_statistics(panel, selected, lonlat, cfg)
    if stats === nothing || length(stats.eligible) <= maximum([3, minimum(cfg.bandwidth_candidates)])
        for variable in selected
            push!(variability, (String(variable), NaN, NaN, NaN, "uncertain", missing, "insufficient_stations"))
        end
        return (; bandwidth_scan, variability, bandwidth=missing)
    end
    distances = _distance_matrix(stats.coords)
    candidates = filter(c -> max(4, stats.p + 1) <= c <= length(stats.eligible) - 1,
        unique(cfg.bandwidth_candidates))
    for candidate in candidates
        rmse = _loocv_panel_rmse(stats, distances, candidate, cfg)
        push!(bandwidth_scan, (candidate, rmse, isfinite(rmse)))
    end
    available = filter(:available => identity, bandwidth_scan)
    if nrow(available) == 0
        for variable in selected
            push!(variability, (String(variable), NaN, NaN, NaN, "uncertain", missing, "no_bandwidth"))
        end
        return (; bandwidth_scan, variability, bandwidth=missing)
    end
    bandwidth = available.neighbors[argmin(available.loocv_rmse)]
    weight_rows = _local_weight_rows(distances, bandwidth, length(stats.A))
    coefficients = _local_coefficients(stats, distances, bandwidth, cfg; weight_rows)
    observed = [var(coefficients[:, 3 + j]) for j in eachindex(selected)]
    # Drawn up front, from the same `rng` in the same order, so `orders[i]` is exactly what
    # iteration `i` drew for itself; the bodies are then independent and each is a full
    # O(n^2 p^2) coefficient-surface rebuild. `hits` is one row per permutation, reduced
    # afterwards - summing booleans is exact and order-free, so `exceed` is unchanged.
    orders = [station_block_permutation(length(stats.A), rng)
        for _ in 1:cfg.spatial_permutations]
    # `Matrix{Bool}`, not a `BitMatrix`: a BitArray packs 64 entries into one word, so two threads
    # writing different rows of the same column would read-modify-write the same word and lose
    # updates. One byte per entry makes the concurrent writes independent.
    hits = fill(false, cfg.spatial_permutations, length(selected))
    Threads.@threads :greedy for permutation_index in 1:cfg.spatial_permutations
        permuted = _local_coefficients(stats, distances, bandwidth, cfg;
            block_order=orders[permutation_index], weight_rows)
        for j in eachindex(selected)
            value = var(permuted[:, 3 + j])
            hits[permutation_index, j] = isfinite(value) && value >= observed[j]
        end
    end
    exceed = vec(sum(hits, dims=1))
    pvalues = cfg.spatial_permutations > 0 ?
        (exceed .+ 1) ./ (cfg.spatial_permutations + 1) : fill(NaN, length(selected))
    qvalues = all(isfinite, pvalues) ? _bh_adjust(pvalues) : fill(NaN, length(selected))
    for j in eachindex(selected)
        role = isfinite(qvalues[j]) ? (qvalues[j] < cfg.q_threshold ? "local" : "global") : "uncertain"
        push!(variability, (String(selected[j]), observed[j], pvalues[j], qvalues[j],
            role, bandwidth, role == "uncertain" ? "test_failed" : "ok"))
    end
    return (; bandwidth_scan, variability, bandwidth)
end
