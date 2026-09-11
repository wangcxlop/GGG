"""
`train_distances`/`target_distances` let a caller that fits many response columns over the same
station geometry hand in the haversine matrices instead of having them rebuilt here.

The contract is "the distances of the coordinate arguments exactly as passed":
`train_distances[i, j]` is the distance between `train_lonlat[i, :]` and `train_lonlat[j, :]`, and
`target_distances[i, j]` between `train_lonlat[i, :]` and `target_lonlat[j, :]`. The per-column
`indices`/`target_indices` subsetting below is then applied to the matrix the same way it is
applied to the coordinates, so the submatrix is elementwise what `_haversine_matrix` would have
returned for the subset. `nothing` (the default) computes them here, as before.

`unsupported` is forwarded to the *target* hat only, never to the back-fit hat - see
[`_local_hat`](@ref) for why that distinction matters. `:missing` is what lets a target the local
fit could not support surface as a `NaN` prediction instead of a silent zero correction.
"""
function mixed_gwr_predict(
    Xlocal_train::Matrix{Float64}, Xglobal_train::Matrix{Float64}, Ytrain::Matrix{Float64},
    train_lonlat::Matrix{Float64}, Xlocal_target::Matrix{Float64},
    Xglobal_target::Matrix{Float64}, target_lonlat::Matrix{Float64}, bandwidth::Float64,
    kernel::Function; adaptive::Bool=true,
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
    exclude_self::Bool=false, unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
)
    if exclude_self
        size(train_lonlat) == size(target_lonlat) ||
            throw(DimensionMismatch("exclude_self requires matching train and target coordinates"))
        train_lonlat == target_lonlat ||
            throw(ArgumentError("exclude_self requires train and target coordinates in the same order"))
    end
    prediction = fill(NaN, size(Xlocal_target, 1), size(Ytrain, 2))
    converged = trues(size(Ytrain, 2))
    # Keyed on the index vector itself. `Vector{Int}` hashes and compares by content, so the
    # cache behaves exactly as it did on the joined string - without building a ~1.5 KB
    # `String` per time column, which on the joint path (one column per call) was the whole
    # cost of a cache that could never hit.
    cache = Dict{Vector{Int},NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        if count(valid) <= size(Xlocal_train, 2) + size(Xglobal_train, 2) + 2
            converged[time] = false
            continue
        end
        indices = findall(valid)
        operators = get!(cache, indices) do
            local_train_valid = Xlocal_train[indices, :]
            global_train_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(Xlocal_target, 1))
            local_target_valid = Xlocal_target[target_indices, :]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            # Clamping to the valid station count only makes sense for an adaptive neighbor
            # count; a fixed-km bandwidth must be applied as given.
            adjusted_bandwidth = adaptive ? min(bandwidth, length(indices) - 1) : bandwidth
            train_distances_valid = train_distances === nothing ?
                _haversine_matrix(lonlat_valid, lonlat_valid) :
                _distance_subset(train_distances, indices, indices)
            target_distances_valid = target_distances === nothing ?
                _haversine_matrix(lonlat_valid, target_lonlat_valid) :
                _distance_subset(target_distances, indices, target_indices)
            local_hat = _local_hat(
                local_train_valid, local_train_valid,
                train_distances_valid, adjusted_bandwidth, kernel;
                adaptive, ridge,
            )
            target_hat = _local_hat(
                local_train_valid, local_target_valid,
                target_distances_valid, adjusted_bandwidth, kernel;
                adaptive, ridge, exclude_self, unsupported,
            )
            global_hat = _global_projection(global_train_valid; ridge)
            (; local_train_valid, global_train_valid, global_target_valid,
                target_indices, local_hat, target_hat, global_hat)
        end
        y = Vector{Float64}(Ytrain[indices, time])
        fitted = _backfit_components(
            y, [operators.local_hat], operators.global_hat; tolerance, max_iterations,
        )
        if !fitted.converged
            converged[time] = false
            continue
        end
        global_beta = isempty(operators.global_train_valid) ? Float64[] :
            (operators.global_train_valid' * operators.global_train_valid + ridge * I) \
            (operators.global_train_valid' * (y - fitted.local_components[1]))
        partial = y - (isempty(operators.global_train_valid) ? zeros(length(y)) :
            operators.global_train_valid * global_beta)
        target_prediction = operators.target_hat * partial
        isempty(operators.global_target_valid) ||
            (target_prediction .+= operators.global_target_valid * global_beta)
        prediction[operators.target_indices, time] = target_prediction
    end
    return prediction, converged
end

function multiscale_gwr_predict(
    local_train::Vector{Matrix{Float64}}, Xglobal_train::Matrix{Float64},
    Ytrain::Matrix{Float64}, train_lonlat::Matrix{Float64},
    local_target::Vector{Matrix{Float64}}, Xglobal_target::Matrix{Float64},
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200, exclude_self::Bool=false,
    unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
)
    if exclude_self
        isempty(Xglobal_train) ||
            throw(ArgumentError("multiscale exclude_self does not support global variables"))
        size(train_lonlat) == size(target_lonlat) ||
            throw(DimensionMismatch("exclude_self requires matching train and target coordinates"))
        train_lonlat == target_lonlat ||
            throw(ArgumentError("exclude_self requires train and target coordinates in the same order"))
    end
    prediction = fill(NaN, size(target_lonlat, 1), size(Ytrain, 2))
    converged = trues(size(Ytrain, 2))
    # Keyed on the index vector itself. `Vector{Int}` hashes and compares by content, so the
    # cache behaves exactly as it did on the joined string - without building a ~1.5 KB
    # `String` per time column, which on the joint path (one column per call) was the whole
    # cost of a cache that could never hit.
    cache = Dict{Vector{Int},NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        min_required = sum(size(X, 2) for X in local_train) + size(Xglobal_train, 2) + 2
        if count(valid) <= min_required
            converged[time] = false
            continue
        end
        indices = findall(valid)
        operators = get!(cache, indices) do
            local_valid = [X[indices, :] for X in local_train]
            global_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(target_lonlat, 1))
            local_target_valid = [X[target_indices, :] for X in local_target]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            # Clamping to the valid station count only makes sense for an adaptive neighbor
            # count; a fixed-km bandwidth must be applied as given.
            adjusted_bandwidths = adaptive ? min.(bandwidths, length(indices) - 1) : bandwidths
            # Same `train_distances`/`target_distances` contract as `mixed_gwr_predict`; see its
            # docstring. Local names differ from the keywords so the `do` block does not capture.
            train_distances_valid = train_distances === nothing ?
                _haversine_matrix(lonlat_valid, lonlat_valid) :
                _distance_subset(train_distances, indices, indices)
            train_hats = [_local_hat(X, X, train_distances_valid, bw, kernel; adaptive, ridge)
                for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            global_hat = _global_projection(global_valid; ridge)
            target_hats = if exclude_self
                [_local_hat(X, X, train_distances_valid, bw, kernel;
                    adaptive, ridge, exclude_self=true, unsupported)
                    for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            else
                target_distances_valid = target_distances === nothing ?
                    _haversine_matrix(lonlat_valid, target_lonlat_valid) :
                    _distance_subset(target_distances, indices, target_indices)
                [_local_hat(X, Xt, target_distances_valid, bw, kernel;
                    adaptive, ridge, unsupported)
                    for (X, Xt, bw) in zip(local_valid, local_target_valid, adjusted_bandwidths)]
            end
            (; local_valid, global_valid, global_target_valid, target_indices,
                train_hats, target_hats, global_hat)
        end
        y = Vector{Float64}(Ytrain[indices, time])
        fitted = _backfit_components(
            y, operators.train_hats, operators.global_hat; tolerance, max_iterations,
        )
        if !fitted.converged
            converged[time] = false
            continue
        end
        local_sum = reduce(+, fitted.local_components)
        global_beta = isempty(operators.global_valid) ? Float64[] :
            (operators.global_valid' * operators.global_valid + ridge * I) \
            (operators.global_valid' * (y - local_sum))
        global_train_component = isempty(operators.global_valid) ? zeros(length(y)) :
            operators.global_valid * global_beta
        target_prediction = isempty(operators.global_target_valid) ?
            zeros(length(operators.target_indices)) : operators.global_target_valid * global_beta
        for group_index in eachindex(operators.local_valid)
            partial = y - global_train_component
            for other in eachindex(fitted.local_components)
                other == group_index || (partial .-= fitted.local_components[other])
            end
            target_prediction .+= operators.target_hats[group_index] * partial
        end
        prediction[operators.target_indices, time] = target_prediction
    end
    return prediction, converged
end
