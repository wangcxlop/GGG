# The per-row `deg2rad` and `cos(lat)` are hoisted out of the double loop: they depend on `i`
# (or `j`) alone, and the old body recomputed each row's pair once per column. Same values, same
# order of operations per entry - `a` is still `sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2`
# built from the identical `Float64`s - so the matrix is bit-for-bit what it was.
function _haversine_matrix(train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64})
    ntrain, ntarget = size(train_lonlat, 1), size(target_lonlat, 1)
    result = Matrix{Float64}(undef, ntrain, ntarget)
    radius = 6378.388
    train_lon = [deg2rad(train_lonlat[i, 1]) for i in 1:ntrain]
    train_lat = [deg2rad(train_lonlat[i, 2]) for i in 1:ntrain]
    train_coslat = cos.(train_lat)
    @inbounds for j in 1:ntarget
        lon2, lat2 = deg2rad(target_lonlat[j, 1]), deg2rad(target_lonlat[j, 2])
        coslat2 = cos(lat2)
        for i in 1:ntrain
            lat1 = train_lat[i]
            dlon, dlat = lon2 - train_lon[i], lat2 - lat1
            a = sin(dlat / 2)^2 + train_coslat[i] * coslat2 * sin(dlon / 2)^2
            result[i, j] = 2radius * asin(sqrt(clamp(a, 0.0, 1.0)))
        end
    end
    return result
end

"""
General local GWR weights for one target, written into `weights` and using `buffer` as scratch
when `adaptive`. `kernel(dist, bw)` is the weight formula (e.g. bisquare, gaussian, ...).

For an adaptive bandwidth, only the k-th smallest finite distance is needed, so `partialsort!`
on a reused buffer replaces `sort(filter(isfinite, distances))`; it returns the same value in
linear time and without allocating. For a fixed bandwidth, `bw` is applied directly with no
neighbor-count conversion. `distances` is not modified.
"""
# Reimplements the library's GWR weighting (`gw_weight!` in src/core/gw_weight.jl, kernels in
# src/core/kernel.jl) instead of calling it, and not identically: the adaptive bandwidth here is the
# distance to the k-th nearest neighbour via `partialsort!`, where `gw_weight!` uses the `dn = bw/n`
# ratio. Both are defensible, but they are different weights for the same nominal bandwidth, so
# the DEM path's `bw` is not directly comparable to the rest of the GWR family's. Unifying them
# changes fitted weights and therefore results - flagged, not fixed.
function _gw_local_weights!(
    weights::Vector{Float64}, distances::Vector{Float64}, bw::Float64,
    kernel::Function, buffer::Vector{Float64}; adaptive::Bool=true,
)
    n = length(distances)
    if !adaptive
        @inbounds for i in 1:n
            d = distances[i]
            weights[i] = isfinite(d) ? kernel(d, bw) : 0.0
        end
        return weights
    end
    finite_count = 0
    @inbounds for i in 1:n
        d = distances[i]
        isfinite(d) || continue
        finite_count += 1
        buffer[finite_count] = d
    end
    if finite_count == 0
        fill!(weights, 0.0)
        return weights
    end
    finite_distances = view(buffer, 1:finite_count)
    # Kept exactly as it was, quirk included: with a single finite distance `clamp` returns
    # `hi = 1` for the usual `neighbors >= 2`, but `lo = 2` for `neighbors <= 1`, which then
    # asks for an element that is not there. Both the old indexing and `partialsort!` below
    # raise `BoundsError` on that input, and no caller can reach it — the bandwidth grids
    # start at 8 and callers clamp with `min(bandwidth, n - 1)`.
    k = clamp(Int(round(bw)), 2, finite_count)
    # `partialsort!` permutes `finite_distances` without dropping elements, so the
    # `maximum` fallback below still sees the same multiset.
    bandwidth = partialsort!(finite_distances, k)
    bandwidth > 0 || (bandwidth = maximum(finite_distances))
    if !(bandwidth > 0)
        @inbounds for i in 1:n
            weights[i] = distances[i] == 0 ? 1.0 : 0.0
        end
        return weights
    end
    # `kernel` (not a hardcoded `d < bandwidth` cutoff) decides whether/how a distance beyond
    # `bandwidth` still contributes: that's a no-op for bisquare/tricube/boxcar, which are zero
    # there anyway, but gaussian/exponential have no hard cutoff and must keep tapering.
    @inbounds for i in 1:n
        d = distances[i]
        weights[i] = isfinite(d) ? kernel(d, bandwidth) : 0.0
    end
    return weights
end

_bisquare_kernel(dist::Float64, bw::Float64) = dist > bw ? 0.0 : (1 - (dist / bw)^2)^2

"""Adaptive bisquare weights; a thin, behaviour-preserving wrapper over [`_gw_local_weights!`](@ref)."""
function _adaptive_bisquare!(
    weights::Vector{Float64}, distances::Vector{Float64}, neighbors::Int,
    buffer::Vector{Float64},
)
    _gw_local_weights!(weights, distances, Float64(neighbors), _bisquare_kernel, buffer; adaptive=true)
end

function _adaptive_bisquare(distances::Vector{Float64}, neighbors::Int)
    n = length(distances)
    return _adaptive_bisquare!(
        Vector{Float64}(undef, n), distances, neighbors, Vector{Float64}(undef, n),
    )
end

"""
Per-target local weights for [`_gwr_smoothers`](@ref): the surviving neighbour indices and their
weights, for every target.

Split out because these depend on `distances`, `neighbors` and `exclude_self` alone — never on
the design matrix. `spatial_variability_test` rebuilds a smoother a thousand times over the same
geometry with only `X` permuted, and used to recompute this, including a sort per target, every
single time.

The `p + 1` sufficiency check lives here too, so it is raised from the caller's own task rather
than from inside a threaded permutation loop, where it would surface wrapped.
"""
function _gwr_local_weights(
    distances::Matrix{Float64}, neighbors::Int, p::Int; exclude_self::Bool=true,
)
    n = size(distances, 1)
    size(distances) == (n, n) || throw(DimensionMismatch("GWR smoother requires square distances"))
    entries = Vector{Tuple{Vector{Int},Vector{Float64}}}(undef, n)
    d = Vector{Float64}(undef, n)
    buffer = Vector{Float64}(undef, n)
    w = Vector{Float64}(undef, n)
    for target in 1:n
        copyto!(d, view(distances, :, target))
        exclude_self && (d[target] = Inf)
        _adaptive_bisquare!(w, d, neighbors, buffer)
        valid = w .> 0
        count(valid) >= p + 1 || throw(ArgumentError("insufficient local observations for bandwidth $neighbors"))
        indices = findall(valid)
        entries[target] = (indices, w[indices])
    end
    return entries
end

function _gwr_smoothers(
    X::Matrix{Float64}, distances::Matrix{Float64}, neighbors::Int;
    ridge::Float64=1e-8, exclude_self::Bool=true, local_weights=nothing,
)
    n, p = size(X)
    size(distances) == (n, n) || throw(DimensionMismatch("GWR smoother requires square distances"))
    entries = local_weights === nothing ?
        _gwr_local_weights(distances, neighbors, p; exclude_self) : local_weights
    smoothers = [zeros(Float64, n, n) for _ in 1:p]
    for target in 1:n
        indices, wv = entries[target]
        Xv = X[indices, :]
        A = Xv' * (wv .* Xv) + ridge * I
        B = A \ (Xv' .* wv')
        for coefficient in 1:p
            smoothers[coefficient][target, indices] = B[coefficient, :]
        end
    end
    return smoothers
end

function _select_gwr_bandwidth(
    X::Matrix{Float64}, y::Vector{Float64}, distances::Matrix{Float64}, candidates::Vector{Int};
    ridge::Float64=1e-8,
)
    rows = NamedTuple[]
    best = nothing
    for bandwidth in candidates
        bandwidth < size(X, 1) || continue
        try
            smoothers = _gwr_smoothers(X, distances, bandwidth; ridge, exclude_self=true)
            prediction = zeros(Float64, length(y))
            for j in axes(X, 2)
                prediction .+= X[:, j] .* (smoothers[j] * y)
            end
            rmse = sqrt(mean(abs2, prediction - y))
            push!(rows, (; bandwidth, RMSE=rmse, status="success", error=""))
            if best === nothing || rmse < best.RMSE
                best = (; bandwidth, RMSE=rmse, smoothers)
            end
        catch error
            push!(rows, (; bandwidth, RMSE=Inf, status="failed", error=sprint(showerror, error)))
        end
    end
    best === nothing && throw(ArgumentError("all GWR bandwidth candidates failed"))
    return best, DataFrame(rows)
end
