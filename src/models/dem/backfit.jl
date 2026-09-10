"""
Rows `rows` and columns `columns` of a supplied distance matrix, without copying when the
selection is the whole matrix.

The gather is a 286 KB allocation at fold scale and it is by far the most common case that
nothing is actually being dropped - every station is valid for most hours. `dynamic_covariate_predict`
fits one hour per task across every core, and at that width the large-object allocations, not the
arithmetic, are what stop the loop scaling: it runs fastest on four threads and gets *slower* on
twenty-four. A view would still allocate a wrapper and would push a non-`Matrix` type through
`_local_hat`; returning the matrix itself keeps both the type and the values exactly as they were.
"""
function _distance_subset(distances::Matrix{Float64}, rows, columns)
    length(rows) == size(distances, 1) && length(columns) == size(distances, 2) &&
        return distances
    return distances[rows, columns]
end

function _local_hat(
    Xtrain::Matrix{Float64}, Xtarget::Matrix{Float64}, distances::Matrix{Float64},
    bw::Float64, kernel::Function; adaptive::Bool=true, ridge::Float64=1e-8,
    exclude_self::Bool=false,
)
    ntrain, p = size(Xtrain)
    ntarget = size(Xtarget, 1)
    size(Xtarget, 2) == p || throw(DimensionMismatch("local design columns differ"))
    size(distances) == (ntrain, ntarget) || throw(DimensionMismatch("distance dimensions differ"))
    exclude_self && ntrain != ntarget &&
        throw(DimensionMismatch("exclude_self requires matching train and target rows"))
    H = zeros(Float64, ntarget, ntrain)
    # Allocated once per call rather than once per target: the old body allocated ten
    # arrays inside the loop and spent most of its time in the allocator, not in BLAS.
    # All buffers are call-local, so the `Threads.@threads` callers stay safe.
    d = Vector{Float64}(undef, ntrain)
    w = Vector{Float64}(undef, ntrain)
    buffer = Vector{Float64}(undef, ntrain)
    indices = Vector{Int}(undef, ntrain)
    wv = Vector{Float64}(undef, ntrain)
    hv = Vector{Float64}(undef, ntrain)
    Xv = Matrix{Float64}(undef, ntrain, p)
    WX = Matrix{Float64}(undef, ntrain, p)
    A = Matrix{Float64}(undef, p, p)
    xt = Vector{Float64}(undef, p)
    @inbounds for target in 1:ntarget
        copyto!(d, view(distances, :, target))
        exclude_self && (d[target] = Inf)
        _gw_local_weights!(w, d, bw, kernel, buffer; adaptive)
        n_valid = 0
        for i in 1:ntrain
            w[i] > 0 || continue
            n_valid += 1
            indices[n_valid] = i
            wv[n_valid] = w[i]
            for j in 1:p
                Xv[n_valid, j] = Xtrain[i, j]
                WX[n_valid, j] = w[i] * Xtrain[i, j]
            end
        end
        n_valid >= p + 1 || continue
        Xvalid = view(Xv, 1:n_valid, :)
        mul!(A, transpose(Xvalid), view(WX, 1:n_valid, :))
        for j in 1:p
            A[j, j] += ridge
        end
        for j in 1:p
            xt[j] = Xtarget[target, j]
        end
        # `A` is symmetric, so `xt' * (A \ (Xv' .* wv')) == (A \ xt)' * (Xv' .* wv')`.
        # Solving for the p-vector first and contracting afterwards replaces a solve with
        # `n_valid` right-hand sides by one with a single right-hand side, which is the bulk
        # of the saving. The plain `Matrix` backslash is kept deliberately (rather than a
        # Cholesky) so behaviour on ill-conditioned or non-finite `A` is unchanged.
        coefficients = A \ xt
        hvalid = view(hv, 1:n_valid)
        mul!(hvalid, Xvalid, coefficients)
        for t in 1:n_valid
            H[target, indices[t]] = hvalid[t] * wv[t]
        end
    end
    return H
end

function _global_projection(X::Matrix{Float64}; ridge::Float64=1e-8)
    isempty(X) && return zeros(Float64, size(X, 1), size(X, 1))
    return X * ((X' * X + ridge * I) \ X')
end

function _backfit_components(
    y::Vector{Float64}, local_hats::Vector{Matrix{Float64}}, global_hat::Matrix{Float64};
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_components = [zeros(Float64, length(y)) for _ in local_hats]
    global_component = global_hat * y
    previous = Inf
    for iteration in 1:max_iterations
        for j in eachindex(local_hats)
            partial = y - global_component
            for k in eachindex(local_components)
                k == j || (partial .-= local_components[k])
            end
            local_components[j] = local_hats[j] * partial
        end
        local_sum = isempty(local_components) ? zeros(Float64, length(y)) : reduce(+, local_components)
        global_component = global_hat * (y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, y - fitted)
        change = isfinite(previous) ? abs(previous - rss) / max(previous, eps()) : Inf
        change < tolerance && return (; local_components, global_component, fitted, converged=true, iterations=iteration)
        previous = rss
    end
    local_sum = isempty(local_components) ? zeros(Float64, length(y)) : reduce(+, local_components)
    fitted = local_sum + global_component
    return (; local_components, global_component, fitted, converged=false, iterations=max_iterations)
end

function _mixed_fit_complete(
    Xlocal::Matrix{Float64}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    distances::Matrix{Float64}, bw::Float64, kernel::Function; adaptive::Bool=true,
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_hat = _local_hat(Xlocal, Xlocal, distances, bw, kernel; adaptive, ridge)
    global_hat = _global_projection(Xglobal; ridge)
    result = _backfit_components(
        y, [local_hat], global_hat; tolerance, max_iterations,
    )
    return merge(result, (; local_hat, global_hat))
end

function _linear_loocv_rmse(hat::Matrix{Float64}, y::Vector{Float64})
    fitted = hat * y
    denominator = 1 .- diag(hat)
    valid = abs.(denominator) .> sqrt(eps())
    any(valid) || return Inf
    errors = (y[valid] - fitted[valid]) ./ denominator[valid]
    return sqrt(mean(abs2, errors))
end

function select_mixed_bandwidth(
    Xlocal::Matrix{Float64}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    lonlat::Matrix{Float64}, candidates::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    distances = _haversine_matrix(lonlat, lonlat)
    global_hat = _global_projection(Xglobal; ridge)
    rows = NamedTuple[]
    best_bandwidth = 0.0
    best_rmse = Inf
    for bandwidth in candidates
        # A fixed-km bandwidth has no relationship to the station count, so the sanity check
        # against `length(y)` only applies to the adaptive neighbor-count case.
        (!adaptive || bandwidth < length(y)) || continue
        try
            local_hat = _local_hat(Xlocal, Xlocal, distances, bandwidth, kernel; adaptive, ridge)
            identity_fit = _backfit_components(
                Matrix{Float64}(I, length(y), length(y)), [local_hat], global_hat;
                tolerance, max_iterations,
            )
            hat = identity_fit.fitted
            rmse = _linear_loocv_rmse(hat, y)
            push!(rows, (; bandwidth, RMSE=rmse, converged=identity_fit.converged,
                iterations=identity_fit.iterations, status="success", error=""))
            if identity_fit.converged && rmse < best_rmse
                best_bandwidth, best_rmse = bandwidth, rmse
            end
        catch error
            push!(rows, (; bandwidth, RMSE=Inf, converged=false, iterations=0,
                status="failed", error=sprint(showerror, error)))
        end
    end
    best_bandwidth > 0 || throw(ArgumentError("all Mixed GWR bandwidth candidates failed"))
    return best_bandwidth, DataFrame(rows)
end

function _backfit_components(
    Y::Matrix{Float64}, local_hats::Vector{Matrix{Float64}}, global_hat::Matrix{Float64};
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_components = [zeros(Float64, size(Y)) for _ in local_hats]
    global_component = global_hat * Y
    previous = Inf
    for iteration in 1:max_iterations
        for j in eachindex(local_hats)
            partial = Y - global_component
            for k in eachindex(local_components)
                k == j || (partial .-= local_components[k])
            end
            local_components[j] = local_hats[j] * partial
        end
        local_sum = isempty(local_components) ? zeros(Float64, size(Y)) : reduce(+, local_components)
        global_component = global_hat * (Y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, Y - fitted)
        change = isfinite(previous) ? abs(previous - rss) / max(previous, eps()) : Inf
        change < tolerance && return (; local_components, global_component, fitted, converged=true, iterations=iteration)
        previous = rss
    end
    local_sum = isempty(local_components) ? zeros(Float64, size(Y)) : reduce(+, local_components)
    fitted = local_sum + global_component
    return (; local_components, global_component, fitted, converged=false, iterations=max_iterations)
end

function _multiscale_fit_complete(
    local_groups::Vector{Matrix{Float64}}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    distances::Matrix{Float64}, bandwidths::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    length(local_groups) == length(bandwidths) || throw(DimensionMismatch("bandwidth count differs"))
    hats = [_local_hat(X, X, distances, bandwidth, kernel; adaptive, ridge)
        for (X, bandwidth) in zip(local_groups, bandwidths)]
    global_hat = _global_projection(Xglobal; ridge)
    result = _backfit_components(y, hats, global_hat; tolerance, max_iterations)
    return merge(result, (; hats, global_hat))
end

function select_multiscale_bandwidths(
    local_groups::Vector{Matrix{Float64}}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    lonlat::Matrix{Float64}, candidates::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    isempty(local_groups) && return Float64[], DataFrame(), true
    # A fixed-km bandwidth has no relationship to the station count, so the sanity check
    # against `length(y)` only applies to the adaptive neighbor-count case.
    usable = adaptive ? filter(<(length(y)), candidates) : candidates
    isempty(usable) && throw(ArgumentError("no usable multiscale bandwidth candidates"))
    distances = _haversine_matrix(lonlat, lonlat)
    global_hat = _global_projection(Xglobal; ridge)
    bandwidths = fill(last(usable), length(local_groups))
    components = [zeros(Float64, length(y)) for _ in local_groups]
    global_component = global_hat * y
    rows = NamedTuple[]
    previous_rss = Inf
    for iteration in 1:max_iterations
        changed = false
        for group_index in eachindex(local_groups)
            partial = y - global_component
            for other in eachindex(components)
                other == group_index || (partial .-= components[other])
            end
            best_bw, best_rmse, best_hat = 0.0, Inf, nothing
            for bandwidth in usable
                try
                    hat = _local_hat(
                        local_groups[group_index], local_groups[group_index], distances,
                        bandwidth, kernel; adaptive, ridge,
                    )
                    rmse = _linear_loocv_rmse(hat, partial)
                    push!(rows, (; iteration, group_index, bandwidth, RMSE=rmse,
                        status="success", error="", selected=false))
                    if rmse < best_rmse
                        best_bw, best_rmse, best_hat = bandwidth, rmse, hat
                    end
                catch error
                    push!(rows, (; iteration, group_index, bandwidth, RMSE=Inf,
                        status="failed", error=sprint(showerror, error), selected=false))
                end
            end
            best_hat === nothing && return bandwidths, DataFrame(rows), false
            selected_index = findlast(row -> row.iteration == iteration &&
                row.group_index == group_index && row.bandwidth == best_bw && row.status == "success", rows)
            rows[selected_index] = merge(rows[selected_index], (; selected=true))
            changed |= bandwidths[group_index] != best_bw
            bandwidths[group_index] = best_bw
            components[group_index] = best_hat * partial
        end
        local_sum = reduce(+, components)
        global_component = global_hat * (y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, y - fitted)
        relative_change = isfinite(previous_rss) ?
            abs(previous_rss - rss) / max(previous_rss, eps()) : Inf
        (!changed && relative_change < tolerance) && return bandwidths, DataFrame(rows), true
        previous_rss = rss
    end
    return bandwidths, DataFrame(rows), false
end
