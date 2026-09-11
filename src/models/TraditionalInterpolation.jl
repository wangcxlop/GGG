module TraditionalInterpolation

using LinearAlgebra
using Statistics

export haversine_distance_matrix, local_km_coordinates, nearest_training_distance
export idw_predict, adw_predict, tps_predict, tps_loo_predict

# Station distances are NOT computed against one radius project-wide. This value (the WGS72
# equatorial radius) is also used by MGERPipeline, DEMTerrainExperiment and
# PrecipitationCorrection, but JointCovariateModels and ERA5VariableSelection use the mean-Earth
# 6371.0088 instead. Each component is internally consistent, so nothing is wrong within a given
# model; the two families are simply not on the same metric, which matters when comparing a
# GWR-family bandwidth in km against a joint-covariate one. Left as-is deliberately: unifying
# them moves every published number. Fix it in a commit that re-runs the benchmark, not in a
# refactor.
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

"""
Distance from each validation station to its nearest training station, in km.

The quantity behind the benchmark's `nearest_train_km` stratum: how far a held-out station is
from any gauge the model actually saw. Written out inline in the benchmark's fold loop and again,
with two different surrounding interfaces, in `run_claim_reassessment.jl` and
`verify_fold_rotations.jl`.
"""
nearest_training_distance(train_lonlat::AbstractMatrix, val_lonlat::AbstractMatrix) =
    vec(minimum(haversine_distance_matrix(train_lonlat, val_lonlat), dims=1))


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

"""
Rows of `values` grouped by which stations are finite, as `mask => times`.

Keyed by a `BitVector` rather than a tuple of indices: the key is hashed once per hour, and a
tuple of up to a few hundred `Int`s is both type-unstable and expensive to hash. Groups, and so
every caller's output, are identical either way.
"""
function _valid_groups(values::AbstractMatrix{<:Real})
    groups = Dict{BitVector,Vector{Int}}()
    for time in axes(values, 2)
        mask = BitVector(isfinite(value) for value in @view(values[:, time]))
        push!(get!(groups, mask, Int[]), time)
    end
    return groups
end

"""
Distance-decay weights for one set of contributing stations, as
`(rows, weights, total, coincident_rows)` where `rows` indexes into the training stations.

`selection` holds local indices into `candidates`, and is expected to already be restricted to
the stations that report at the hours these weights will be used for.
"""
function _selection_weights(selection::Vector{Int}, geometry, angular::Bool, coincident::Float64)
    base, unit_x, unit_y = geometry.base, geometry.unit_x, geometry.unit_y
    weights = [base[i] for i in selection]
    if angular && length(selection) > 1
        # Shepard's correction needs only the pairwise cosines, and `cos_ij = u_i . u_j` for unit
        # vectors, so the k x k cosine matrix never has to be formed: `sum_j base_j cos_ij` is
        # `u_i . (U' base)`, two components. That keeps the cost linear in the selection size.
        decay = sum(weights)
        sum_x = sum(weights[j] * unit_x[selection[j]] for j in eachindex(selection))
        sum_y = sum(weights[j] * unit_y[selection[j]] for j in eachindex(selection))
        for j in eachindex(selection)
            denominator = decay - weights[j]
            denominator > 0 || continue
            i = selection[j]
            weights[j] *= 1 + (decay - (unit_x[i] * sum_x + unit_y[i] * sum_y)) / denominator
        end
    end
    return (
        [geometry.candidates[i] for i in selection],
        weights,
        sum(weights),
        [geometry.candidates[i] for i in selection if geometry.distances[i] <= coincident],
    )
end

"""
Shared IDW/ADW kernel: a distance-decay weighted mean, optionally with Shepard's directional
correction.

Both the neighbour selection and the angular correction are formed **per availability group** -
over the stations that actually report at those hours - rather than over the geometric
neighbourhood with the missing values filtered out afterwards. The latter is wrong in a way
renormalising the weights does not repair: a station contributing no value would still consume
a slot in the `neighbors` budget, and would still shadow (or still boost) its directional
neighbours in the correction, leaving the surviving weights shaped against a station geometry
that is not the one being used. `_valid_groups` is what keeps this affordable - the weights are
rebuilt once per distinct missing pattern, not once per hour.
"""
function _weighted_predict(
    train_lonlat::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real},
    target_lonlat::AbstractMatrix{<:Real}; power::Real=2.0,
    neighbors::Union{Nothing,Int}=nothing, angular::Bool=false,
    exclude_self::Bool=false,
)
    _check_inputs(train_lonlat, values, target_lonlat)
    power > 0 || throw(ArgumentError("power must be positive"))
    neighbors === nothing || neighbors >= 1 ||
        throw(ArgumentError("neighbors must be positive or nothing"))
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("exclude_self requires matching training and target rows"))
        @inbounds for i in 1:size(distances, 1)
            distances[i, i] = Inf
        end
    end

    n_target = size(target_lonlat, 1)
    prediction = fill(NaN, n_target, size(values, 2))
    # Missing entries are zeroed rather than skipped term by term: every station carrying a
    # weight reports over its own group's hours, so these zeros never reach a weighted sum.
    clean = [isfinite(value) ? Float64(value) : 0.0 for value in values]
    groups = _valid_groups(values)
    coincident = sqrt(eps(Float64))

    for target in 1:n_target
        target_distances = @view distances[:, target]
        candidates = findall(isfinite, target_distances)
        isempty(candidates) && continue
        candidate_distances = Float64[target_distances[i] for i in candidates]
        # A station on top of the target carries no distance-decay weight; it takes over through
        # the exact-interpolation branch below instead.
        base = [d <= coincident ? 0.0 : d^(-power) for d in candidate_distances]
        # Directions to each candidate, in km about the target - only ADW needs them.
        unit_x, unit_y = if angular
            offsets = local_km_coordinates(
                vcat(Float64.(train_lonlat[candidates, :]), Float64.(target_lonlat[target:target, :]));
                center=(target_lonlat[target, 1], target_lonlat[target, 2]),
            )
            norms = [hypot(offsets[i, 1], offsets[i, 2]) for i in eachindex(candidates)]
            ([norms[i] > 0 ? offsets[i, 1] / norms[i] : 0.0 for i in eachindex(candidates)],
             [norms[i] > 0 ? offsets[i, 2] / norms[i] : 0.0 for i in eachindex(candidates)])
        else
            (Float64[], Float64[])
        end
        geometry = (; candidates, distances=candidate_distances, base, unit_x, unit_y)
        # Nearest-first once per target, so a bounded neighbourhood is taken by walking this
        # until enough stations report rather than by re-ranking the candidates per group.
        nearest_first = neighbors === nothing ? Int[] : sortperm(candidate_distances)
        # A bounded neighbourhood collapses many missing patterns onto the same few stations, so
        # its weights are worth caching. An unbounded one selects every reporting station, which
        # is as distinct as the pattern itself, and caching would only retain what it rebuilds.
        cache = Dict{Vector{Int},Tuple{Vector{Int},Vector{Float64},Float64,Vector{Int}}}()
        selection = Int[]

        for (mask, times) in groups
            # Kept in candidate order either way, so the weighted sum walks rows ascending.
            empty!(selection)
            if neighbors === nothing
                for i in eachindex(candidates)
                    mask[candidates[i]] && push!(selection, i)
                end
            else
                for i in nearest_first
                    mask[candidates[i]] || continue
                    push!(selection, i)
                    length(selection) == neighbors && break
                end
                sort!(selection)
            end
            isempty(selection) && continue

            rows, weights, total, exact = if neighbors === nothing
                _selection_weights(selection, geometry, angular, coincident)
            else
                get!(() -> _selection_weights(selection, geometry, angular, coincident),
                    cache, copy(selection))
            end

            if !isempty(exact)
                # A station on top of the target wins outright, averaged when several coincide.
                for time in times
                    prediction[target, time] = sum(clean[row, time] for row in exact) / length(exact)
                end
            elseif total > 0 && isfinite(total)
                prediction[target, times] =
                    vec(transpose(weights) * @view(clean[rows, times])) ./ total
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

"""
The magnitude the dimensionless `smooth` is measured against: the median non-zero TPS kernel
value over a station network. Since `K = 0.5 r^2 ln r^2` is monotone in `r` past `r = e^-0.5`
km, this is the kernel evaluated at the median pairwise station distance - a proxy for the
network's extent and density.
"""
function _tps_scale(xy::AbstractMatrix{<:Real})
    K = _tps_kernel_matrix(xy, xy)
    positive = abs.(K[.!iszero.(K)])
    return isempty(positive) ? 1.0 : median(positive)
end

"""
Stations reporting at some hour, and the `_tps_scale` of the network they form.

Taken over the whole call rather than per missing-value group. Derived per group - as this was
until the hour-to-hour flicker was measured - the effective smoothing depended on which gauges
happened to report in a given hour, so the same `smooth` meant a different lambda from one hour
to the next. It also re-sorted a fresh n^2/2 vector for every group to do it.
"""
function _tps_reporting_scale(values::AbstractMatrix{<:Real}, xy::AbstractMatrix{<:Real})
    reporting = findall(station -> any(isfinite, @view(values[station, :])), axes(values, 1))
    return length(reporting) >= 2 ? _tps_scale(@view(xy[reporting, :])) : 1.0
end

function _tps_system(xy::Matrix{Float64}, smooth::Float64, scale::Float64)
    n = size(xy, 1)
    n >= 3 || throw(ArgumentError("TPS requires at least three valid stations"))
    P = hcat(ones(Float64, n), xy)
    rank(P) == 3 || throw(ArgumentError("TPS stations must not be collinear"))
    K = _tps_kernel_matrix(xy, xy)
    A = [K + smooth * scale * I P; transpose(P) zeros(Float64, 3, 3)]
    return A, K, P
end

"""
Two-dimensional thin-plate smoothing spline. `smooth` is dimensionless, scaled by the median
non-zero TPS kernel magnitude over the stations that report at least one hour.

That scale is fixed once per call, so the effective smoothing does not depend on which gauges
reported in any given hour. It does still differ between training sets of different extent -
the inner selection split spans a smaller network than the fold it is selected for, so the
applied lambda is not exactly the one validated. That is deliberate, and matches every other
method here: GWR's adaptive bandwidth is a neighbour count whose radius likewise grows when
stations are withheld.
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
    scale = _tps_reporting_scale(values, train_xy_all)
    prediction = fill(NaN, size(target_lonlat, 1), size(values, 2))

    for (mask, times) in _valid_groups(values)
        idx = findall(mask)
        length(idx) >= 3 || continue
        try
            train_xy = Matrix{Float64}(train_xy_all[idx, :])
            A, _, _ = _tps_system(train_xy, Float64(smooth), scale)
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
    scale = _tps_reporting_scale(values, train_xy_all)
    prediction = fill(NaN, size(values))

    for (mask, times) in _valid_groups(values)
        idx = findall(mask)
        length(idx) >= 4 || continue
        try
            xy = Matrix{Float64}(train_xy_all[idx, :])
            A, K, P = _tps_system(xy, Float64(smooth), scale)
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
