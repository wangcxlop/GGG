function _balanced_spatial_folds(
    ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, seed::Int=20260815,
)
    n = length(ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    center_lat = mean(lonlat[:, 2])
    xy = hcat(
        (lonlat[:, 1] .- mean(lonlat[:, 1])) .* (111.32cosd(center_lat)),
        (lonlat[:, 2] .- mean(lonlat[:, 2])) .* 110.57,
    )
    rng = MersenneTwister(seed)
    centers = Matrix{Float64}(undef, k, 2)
    centers[1, :] = xy[rand(rng, 1:n), :]
    nearest = fill(Inf, n)
    for cluster in 2:k
        for station in 1:n
            nearest[station] = min(nearest[station],
                sum(abs2, xy[station, :] - centers[cluster - 1, :]))
        end
        centers[cluster, :] = xy[argmax(nearest), :]
    end
    capacities = fill(div(n, k), k)
    capacities[1:rem(n, k)] .+= 1
    assignments = zeros(Int, n)
    for _ in 1:100
        distances = [sum(abs2, xy[i, :] - centers[j, :]) for i in 1:n, j in 1:k]
        certainty = [let ordered = sort(@view distances[i, :])
            length(ordered) > 1 ? ordered[2] - ordered[1] : Inf
        end for i in 1:n]
        order = sortperm(1:n; by=i -> (-certainty[i], i))
        remaining = copy(capacities)
        updated = zeros(Int, n)
        for station in order
            for cluster in sortperm(@view distances[station, :])
                if remaining[cluster] > 0
                    updated[station] = cluster
                    remaining[cluster] -= 1
                    break
                end
            end
        end
        updated == assignments && break
        assignments = updated
        for cluster in 1:k
            indices = findall(==(cluster), assignments)
            centers[cluster, :] = vec(mean(xy[indices, :], dims=1))
        end
    end
    return [findall(==(fold), assignments) for fold in 1:k]
end

function _mean_wet_residual(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}; threshold::Float64=0.1,
    min_hours::Int=100, time_indices=axes(Yobs, 2),
)
    size(Yobs) == size(Ysat) || throw(DimensionMismatch("observation and satellite matrices differ"))
    response = fill(NaN, size(Yobs, 1))
    counts = zeros(Int, size(Yobs, 1))
    for station in axes(Yobs, 1)
        valid = [time for time in time_indices if
            isfinite(Yobs[station, time]) && isfinite(Ysat[station, time]) &&
            Yobs[station, time] >= threshold]
        counts[station] = length(valid)
        length(valid) >= min_hours &&
            (response[station] = mean(Yobs[station, valid] - Ysat[station, valid]))
    end
    return response, counts
end

mean_wet_residual(args...; kwargs...) = _mean_wet_residual(args...; kwargs...)

function _monthly_rows(
    product::String, times::Vector{DateTime}, terrain::DataFrame,
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, cfg::DEMExperimentConfig,
)
    rows = NamedTuple[]
    for month_value in 6:9
        indices = findall(==(month_value), month.(times))
        isempty(indices) && continue
        response, counts = _mean_wet_residual(
            Yobs, Ysat; threshold=cfg.wet_threshold,
            min_hours=max(10, div(cfg.min_wet_hours, 4)), time_indices=indices,
        )
        valid = isfinite.(response)
        for group in ("elevation", "slope")
            x = Float64.(terrain[!, only(TERRAIN_COLUMNS[group])])
            value = count(valid) >= 3 ? _correlation(x[valid], response[valid]) : NaN
            push!(rows, (; product, month=month_value, variable_group=group,
                correlation=value, direction=isnan(value) ? "unavailable" : value >= 0 ? "positive" : "negative",
                valid_station_count=count(valid), min_wet_count=isempty(counts) ? 0 : minimum(counts)))
        end
    end
    return rows
end

monthly_correlation_rows(args...; kwargs...) = _monthly_rows(args...; kwargs...)

function _designs(
    terrain_train::DataFrame, terrain_target::DataFrame,
    lonlat_train::Matrix{Float64}, lonlat_target::Matrix{Float64},
    role_map::Dict{String,String},
)
    spatial_train = copy(lonlat_train)
    spatial_target = copy(lonlat_target)
    spatial_z, spatial_means, spatial_scales = _standardize_train(spatial_train)
    spatial_target_z = (spatial_target .- spatial_means') ./ spatial_scales'

    local_train = hcat(ones(size(spatial_z, 1)), spatial_z)
    local_target = hcat(ones(size(spatial_target_z, 1)), spatial_target_z)
    global_train = zeros(Float64, size(local_train, 1), 0)
    global_target = zeros(Float64, size(local_target, 1), 0)
    multiscale_train = Matrix{Float64}[
        ones(size(spatial_z, 1), 1), spatial_z[:, 1:1], spatial_z[:, 2:2],
    ]
    multiscale_target = Matrix{Float64}[
        ones(size(spatial_target_z, 1), 1),
        spatial_target_z[:, 1:1], spatial_target_z[:, 2:2],
    ]
    local_group_names = ["intercept", "longitude", "latitude"]

    selected = [group for group in TERRAIN_GROUPS if haskey(role_map, group)]
    if !isempty(selected)
        raw_train, _, column_groups = _terrain_matrix(terrain_train, selected)
        raw_target, _, _ = _terrain_matrix(terrain_target, selected)
        ztrain, means, scales = _standardize_train(raw_train)
        ztarget = (raw_target .- means') ./ scales'
        first_column = 1
        for group in selected
            width = length(TERRAIN_COLUMNS[group])
            columns = first_column:(first_column + width - 1)
            if role_map[group] == "local"
                local_train = hcat(local_train, ztrain[:, columns])
                local_target = hcat(local_target, ztarget[:, columns])
                push!(multiscale_train, ztrain[:, columns])
                push!(multiscale_target, ztarget[:, columns])
                push!(local_group_names, group)
            elseif role_map[group] == "global"
                global_train = hcat(global_train, ztrain[:, columns])
                global_target = hcat(global_target, ztarget[:, columns])
            else
                throw(ArgumentError("uncertain terrain roles cannot be fitted"))
            end
            first_column += width
        end
    end
    return (;
        mixed_local_train=local_train, mixed_local_target=local_target,
        global_train, global_target,
        multiscale_train, multiscale_target, local_group_names,
    )
end

"""Build leakage-safe terrain designs using training-fold centering and scaling."""
terrain_model_designs(args...; kwargs...) = _designs(args...; kwargs...)
