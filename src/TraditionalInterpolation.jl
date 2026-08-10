module TraditionalInterpolation

using LinearAlgebra
using Statistics

export haversine_distance_matrix, local_km_coordinates
export idw_predict, adw_predict, tps_predict, tps_loo_predict

const EARTH_RADIUS_KM = 6378.388

function _check_inputs(
    train_lonlat::AbstractMatrix, values::AbstractMatrix, target_lonlat::AbstractMatrix,
)
    size(train_lonlat, 2) == 2 || throw(DimensionMismatch("train_lonlat must have two columns"))
    size(target_lonlat, 2) == 2 || throw(DimensionMismatch("target_lonlat must have two columns"))
    size(values, 1) == size(train_lonlat, 1) ||
        throw(DimensionMismatch("values rows must match training stations"))
    return nothing
end

"""Pairwise great-circle distances in km, with rows=train and columns=target."""
function haversine_distance_matrix(
    train_lonlat::AbstractMatrix{<:Real}, target_lonlat::AbstractMatrix{<:Real},
)
    n_train = size(train_lonlat, 1)
    n_target = size(target_lonlat, 1)
    distances = Matrix{Float64}(undef, n_train, n_target)
    @inbounds for j in 1:n_target
        lon2 = deg2rad(Float64(target_lonlat[j, 1]))
        lat2 = deg2rad(Float64(target_lonlat[j, 2]))
        for i in 1:n_train
            lon1 = deg2rad(Float64(train_lonlat[i, 1]))
            lat1 = deg2rad(Float64(train_lonlat[i, 2]))
            dlon = lon2 - lon1
            dlat = lat2 - lat1
            a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
            distances[i, j] = 2 * EARTH_RADIUS_KM * asin(sqrt(clamp(a, 0.0, 1.0)))
        end
    end
    return distances
end

"""Equirectangular local coordinates in km, using a common lon/lat centre."""
function local_km_coordinates(
    lonlat::AbstractMatrix{<:Real};
    center::Tuple{<:Real,<:Real}=(mean(lonlat[:, 1]), mean(lonlat[:, 2])),
)
    lon0, lat0 = Float64.(center)
    x = (Float64.(lonlat[:, 1]) .- lon0) .* cosd(lat0) .* (pi / 180) .* EARTH_RADIUS_KM
    y = (Float64.(lonlat[:, 2]) .- lat0) .* (pi / 180) .* EARTH_RADIUS_KM
    return hcat(x, y)
end

function _candidate_indices(distances::AbstractVector{<:Real}, neighbors::Union{Nothing,Int})
    valid = findall(isfinite, distances)
    isempty(valid) && return valid
    if neighbors === nothing || neighbors >= length(valid)
        return valid
    end
    neighbors >= 1 || throw(ArgumentError("neighbors must be positive or nothing"))
    order = partialsortperm(@view(distances[valid]), 1:neighbors)
    return valid[order]
end

function _weighted_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real},
    target_lonlat::AbstractMatrix{<:Real}; power::Real=2.0,
    neighbors::Union{Nothing,Int}=nothing, angular::Bool=false,
    exclude_self::Bool=false,
)
    _check_inputs(train_lonlat, values, target_lonlat)
    power > 0 || throw(ArgumentError("power must be positive"))
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("exclude_self requires matching training and target rows"))
        @inbounds for i in 1:size(distances, 1)
            distances[i, i] = Inf
        end
    end

    n_target = size(target_lonlat, 1)
    n_time = size(values, 2)
    prediction = fill(NaN, n_target, n_time)
    @inbounds for target in 1:n_target
        geom_idx = _candidate_indices(@view(distances[:, target]), neighbors)
        isempty(geom_idx) && continue
        target_xy = local_km_coordinates(
            vcat(Float64.(train_lonlat[geom_idx, :]), Float64.(target_lonlat[target:target, :]));
            center=(target_lonlat[target, 1], target_lonlat[target, 2]),
        )
        vectors = @view target_xy[1:end-1, :]

        zero_local = findall(
            local_index -> distances[geom_idx[local_index], target] <= sqrt(eps(Float64)),
            eachindex(geom_idx),
        )
        nonzero_local = [local_index for local_index in eachindex(geom_idx) if !(local_index in zero_local)]
        base_all = Float64[
            distances[geom_idx[local_index], target]^(-power) for local_index in nonzero_local
        ]

        function angular_weights(local_positions::Vector{Int}, base::Vector{Float64})
            weights = copy(base)
            if angular && length(local_positions) > 1
                unit_vectors = Matrix{Float64}(undef, length(local_positions), 2)
                for (index, local_index) in enumerate(local_positions)
                    vx = vectors[local_index, 1]
                    vy = vectors[local_index, 2]
                    vector_norm = hypot(vx, vy)
                    unit_vectors[index, 1] = vx / vector_norm
                    unit_vectors[index, 2] = vy / vector_norm
                end
                cosine = clamp.(unit_vectors * transpose(unit_vectors), -1.0, 1.0)
                total = sum(base)
                numerator = total .- cosine * base
                denominator = total .- base
                correction = ifelse.(denominator .> 0, numerator ./ denominator, 0.0)
                weights .*= 1 .+ correction
            end
            return weights
        end

        full_weights = angular_weights(nonzero_local, base_all)
        if !isempty(nonzero_local)
            selected_values = Float64.(values[geom_idx[nonzero_local], :])
            valid = isfinite.(selected_values)
            clean_values = ifelse.(valid, selected_values, 0.0)
            numerator = vec(transpose(full_weights) * clean_values)
            denominator = vec(transpose(full_weights) * Float64.(valid))
            usable = isfinite.(denominator) .& (denominator .> 0)
            prediction[target, usable] = numerator[usable] ./ denominator[usable]
        end

        if !isempty(zero_local)
            zero_values = Float64.(values[geom_idx[zero_local], :])
            for time in 1:n_time
                available = filter(isfinite, @view(zero_values[:, time]))
                if !isempty(available)
                    prediction[target, time] = mean(available)
                end
            end
        end
    end
    return prediction
end

"""Inverse-distance interpolation with NaN-aware values."""
function idw_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real},
    target_lonlat::AbstractMatrix{<:Real}; power::Real=2.0,
    neighbors::Union{Nothing,Int}=nothing, exclude_self::Bool=false,
)
    return _weighted_predict(
        train_lonlat, values, target_lonlat;
        power=power, neighbors=neighbors, angular=false, exclude_self=exclude_self,
    )
end

"""
Shepard angular-distance weighting. The IDW weight is multiplied by a directional
correction based on the horizontal angles between the contributing stations.
"""
function adw_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real},
    target_lonlat::AbstractMatrix{<:Real}; power::Real=2.0,
    neighbors::Union{Nothing,Int}=nothing, exclude_self::Bool=false,
)
    return _weighted_predict(
        train_lonlat, values, target_lonlat;
        power=power, neighbors=neighbors, angular=true, exclude_self=exclude_self,
    )
end

@inline function _tps_kernel(r2::Float64)
    return r2 <= eps(Float64) ? 0.0 : 0.5 * r2 * log(r2)
end

function _tps_kernel_matrix(a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real})
    out = Matrix{Float64}(undef, size(a, 1), size(b, 1))
    @inbounds for j in 1:size(b, 1), i in 1:size(a, 1)
        dx = Float64(a[i, 1] - b[j, 1])
        dy = Float64(a[i, 2] - b[j, 2])
        out[i, j] = _tps_kernel(dx * dx + dy * dy)
    end
    return out
end

function _tps_system(xy::Matrix{Float64}, smooth::Float64)
    n = size(xy, 1)
    n >= 3 || throw(ArgumentError("TPS requires at least three valid stations"))
    P = hcat(ones(Float64, n), xy)
    rank(P) == 3 || throw(ArgumentError("TPS stations must not be collinear"))
    K = _tps_kernel_matrix(xy, xy)
    positive = abs.(K[.!iszero.(K)])
    scale = isempty(positive) ? 1.0 : median(positive)
    lambda = smooth * scale
    A = [K + lambda * I P; transpose(P) zeros(Float64, 3, 3)]
    return A, K, P
end

function _valid_groups(values::AbstractMatrix{<:Real})
    groups = Dict{Tuple,Vector{Int}}()
    for time in axes(values, 2)
        key = Tuple(findall(isfinite, @view(values[:, time])))
        push!(get!(groups, key, Int[]), time)
    end
    return groups
end

"""
Two-dimensional thin-plate smoothing spline. `smooth` is dimensionless and is
scaled by the median non-zero TPS kernel magnitude for the current station set.
"""
function tps_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real},
    target_lonlat::AbstractMatrix{<:Real}; smooth::Real=0.01,
)
    _check_inputs(train_lonlat, values, target_lonlat)
    smooth >= 0 || throw(ArgumentError("smooth must be non-negative"))
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    train_xy_all = local_km_coordinates(train_lonlat; center=center)
    target_xy = local_km_coordinates(target_lonlat; center=center)
    prediction = fill(NaN, size(target_lonlat, 1), size(values, 2))

    for (key, times) in _valid_groups(values)
        idx = collect(Int, key)
        length(idx) >= 3 || continue
        try
            train_xy = Matrix{Float64}(train_xy_all[idx, :])
            A, _, _ = _tps_system(train_xy, Float64(smooth))
            rhs = vcat(Float64.(values[idx, times]), zeros(Float64, 3, length(times)))
            coefficients = A \ rhs
            K_target = transpose(_tps_kernel_matrix(train_xy, target_xy))
            P_target = hcat(ones(Float64, size(target_xy, 1)), target_xy)
            prediction[:, times] = hcat(K_target, P_target) * coefficients
        catch
            # The caller records missing output as a failed/insufficient candidate.
        end
    end
    return prediction
end

"""Efficient leave-one-station-out predictions for TPS smoothing selection."""
function tps_loo_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real};
    smooth::Real=0.01,
)
    _check_inputs(train_lonlat, values, train_lonlat)
    smooth > 0 || throw(ArgumentError("TPS LOOCV requires positive smoothing"))
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    train_xy_all = local_km_coordinates(train_lonlat; center=center)
    prediction = fill(NaN, size(values))

    for (key, times) in _valid_groups(values)
        idx = collect(Int, key)
        length(idx) >= 4 || continue
        try
            xy = Matrix{Float64}(train_xy_all[idx, :])
            A, K, P = _tps_system(xy, Float64(smooth))
            selector = vcat(Matrix{Float64}(I, length(idx), length(idx)), zeros(3, length(idx)))
            mapping = A \ selector
            H = hcat(K, P) * mapping
            fitted = H * Float64.(values[idx, times])
            leverage_den = 1 .- diag(H)
            for (local_i, global_i) in enumerate(idx)
                abs(leverage_den[local_i]) <= 1e-10 && continue
                prediction[global_i, times] = Float64.(values[global_i, times]) .-
                    (Float64.(values[global_i, times]) .- fitted[local_i, :]) ./
                    leverage_den[local_i]
            end
        catch
            # Leave this missing pattern as NaN; coverage checks reject unusable candidates.
        end
    end
    return prediction
end

end
