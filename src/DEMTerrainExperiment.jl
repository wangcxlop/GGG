module DEMTerrainExperiment

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

export DEMExperimentConfig, bh_adjust, terrain_screen, spatial_variability_test
export mixed_gwr_predict, multiscale_gwr_predict, select_mixed_bandwidth
export select_multiscale_bandwidths, run_dem_experiment
export mean_wet_residual, monthly_correlation_rows, terrain_model_designs
export terrain_groups, terrain_columns

const TERRAIN_GROUPS = ("elevation", "slope", "aspect")
const TERRAIN_COLUMNS = Dict(
    "elevation" => [:elevation_m],
    "slope" => [:slope_deg],
    "aspect" => [:aspect_sin, :aspect_cos],
)

terrain_groups() = collect(TERRAIN_GROUPS)
terrain_columns(group::String) = copy(TERRAIN_COLUMNS[group])

Base.@kwdef struct DEMExperimentConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    k::Int = 5
    seed::Int = 20260815
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    screen_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
    tolerance::Float64 = 1e-5
    max_iterations::Int = 200
end

function bh_adjust(pvalues::AbstractVector{<:Real})
    p = Float64.(pvalues)
    m = length(p)
    order = sortperm(p)
    adjusted = fill(NaN, m)
    running = 1.0
    for rank in m:-1:1
        index = order[rank]
        running = min(running, p[index] * m / rank)
        adjusted[index] = min(running, 1.0)
    end
    return adjusted
end

function _rank_ties(values::AbstractVector{<:Real})
    order = sortperm(values)
    ranks = zeros(Float64, length(values))
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && values[order[j + 1]] == values[order[i]]
            j += 1
        end
        rank = (i + j) / 2
        ranks[order[i:j]] .= rank
        i = j + 1
    end
    return ranks
end

function _correlation(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(DimensionMismatch("correlation vectors must match"))
    length(x) >= 3 || return NaN
    std(x) > 0 && std(y) > 0 || return NaN
    return cor(x, y)
end

_spearman(x, y) = _correlation(_rank_ties(x), _rank_ties(y))

function _permutation_pvalue(
    statistic::Function, x, y::Vector{Float64}, observed::Float64,
    permutations::Int, rng::AbstractRNG,
)
    permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
    permutations == 0 && return NaN
    isnan(observed) && return 1.0
    exceed = 0
    permuted = copy(y)
    for _ in 1:permutations
        shuffle!(rng, permuted)
        value = statistic(x, permuted)
        exceed += isfinite(value) && value >= observed
    end
    return (exceed + 1) / (permutations + 1)
end

function _joint_f_statistic(X::Matrix{Float64}, y::Vector{Float64}; ridge::Float64=1e-10)
    n, p = size(X)
    n > p + 1 || return NaN
    centered = y .- mean(y)
    rss0 = sum(abs2, centered)
    rss0 > eps() || return NaN
    design = hcat(ones(n), X)
    beta = (design' * design + ridge * I) \ (design' * y)
    rss1 = sum(abs2, y - design * beta)
    rss1 > 0 || return Inf
    return max((rss0 - rss1) / p, 0.0) / (rss1 / (n - p - 1))
end

function _vif_values(X::Matrix{Float64})
    p = size(X, 2)
    p == 0 && return Float64[]
    p == 1 && return [1.0]
    result = fill(Inf, p)
    for j in 1:p
        others = setdiff(1:p, j)
        y = X[:, j]
        design = hcat(ones(size(X, 1)), X[:, others])
        beta = design \ y
        rss = sum(abs2, y - design * beta)
        tss = sum(abs2, y .- mean(y))
        r2 = tss > 0 ? clamp(1 - rss / tss, 0.0, 1.0) : 1.0
        result[j] = r2 < 1 ? 1 / (1 - r2) : Inf
    end
    return result
end

function _standardize_train(X::Matrix{Float64})
    means = vec(mean(X, dims=1))
    scales = vec(std(X, dims=1))
    all(scales .> 0) || throw(ArgumentError("constant predictor in terrain design"))
    return (X .- means') ./ scales', means, scales
end

function _terrain_matrix(terrain::DataFrame, groups::Vector{String})
    columns = Symbol[]
    column_groups = String[]
    for group in groups
        for column in TERRAIN_COLUMNS[group]
            push!(columns, column)
            push!(column_groups, group)
        end
    end
    return Matrix{Float64}(terrain[:, columns]), columns, column_groups
end

"""Traditional correlation, joint-aspect, BH, and grouped VIF screening."""
function terrain_screen(
    terrain::DataFrame, response::Vector{Float64}; product::String="product",
    permutations::Int=999, q_threshold::Float64=0.05, vif_threshold::Float64=5.0,
    seed::Int=20260815,
)
    length(response) == nrow(terrain) || throw(DimensionMismatch("terrain and response rows differ"))
    valid = isfinite.(response)
    count(valid) >= 8 || throw(ArgumentError("at least eight stations are required for screening"))
    y = response[valid]
    rows = NamedTuple[]
    raw_p = Float64[]
    rng = MersenneTwister(seed)

    for group in TERRAIN_GROUPS
        columns = TERRAIN_COLUMNS[group]
        X = Matrix{Float64}(terrain[valid, columns])
        if group == "aspect"
            observed = _joint_f_statistic(X, y)
            pvalue = _permutation_pvalue(
                (a, b) -> _joint_f_statistic(a, b), X, y, observed, permutations, rng,
            )
            push!(raw_p, pvalue)
            push!(rows, (;
                product, variable_group=group, test="joint_F", statistic=observed,
                pearson=NaN, spearman=NaN, pvalue, qvalue=NaN,
                direction_stable=true, selected_pre_vif=false, selected=false,
                exclusion_reason="",
            ))
        else
            x = vec(X)
            pearson = _correlation(x, y)
            spearman = _spearman(x, y)
            observed = abs(pearson)
            pvalue = _permutation_pvalue(
                (a, b) -> abs(_correlation(vec(a), b)), X, y, observed, permutations, rng,
            )
            direction_stable = isfinite(pearson) && isfinite(spearman) && pearson * spearman >= 0
            push!(raw_p, pvalue)
            push!(rows, (;
                product, variable_group=group, test="pearson", statistic=pearson,
                pearson, spearman, pvalue, qvalue=NaN,
                direction_stable, selected_pre_vif=false, selected=false,
                exclusion_reason="",
            ))
        end
    end

    adjusted = bh_adjust(raw_p)
    selected_groups = String[]
    for i in eachindex(rows)
        selected = adjusted[i] < q_threshold && rows[i].direction_stable
        reason = selected ? "" :
            (!rows[i].direction_stable ? "pearson_spearman_direction_conflict" : "q_not_significant")
        rows[i] = merge(rows[i], (;
            qvalue=adjusted[i], selected_pre_vif=selected, selected, exclusion_reason=reason,
        ))
        selected && push!(selected_groups, rows[i].variable_group)
    end

    vif_rows = NamedTuple[]
    active = copy(selected_groups)
    while !isempty(active)
        X, columns, column_groups = _terrain_matrix(terrain[valid, :], active)
        Xz, _, _ = _standardize_train(X)
        vifs = _vif_values(Xz)
        for (column, group, value) in zip(columns, column_groups, vifs)
            push!(vif_rows, (; product, variable_group=group, variable=String(column), VIF=value,
                iteration=length(selected_groups) - length(active) + 1))
        end
        maximum(vifs) < vif_threshold && break
        remove_group = column_groups[argmax(vifs)]
        filter!(!=(remove_group), active)
        row_index = findfirst(row -> row.variable_group == remove_group, rows)
        rows[row_index] = merge(rows[row_index], (;
            selected=false, exclusion_reason="VIF_at_or_above_$(vif_threshold)",
        ))
    end
    return DataFrame(rows), DataFrame(vif_rows), active
end

function _haversine_matrix(train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64})
    result = Matrix{Float64}(undef, size(train_lonlat, 1), size(target_lonlat, 1))
    radius = 6378.388
    for j in axes(target_lonlat, 1), i in axes(train_lonlat, 1)
        lon1, lat1 = deg2rad(train_lonlat[i, 1]), deg2rad(train_lonlat[i, 2])
        lon2, lat2 = deg2rad(target_lonlat[j, 1]), deg2rad(target_lonlat[j, 2])
        dlon, dlat = lon2 - lon1, lat2 - lat1
        a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
        result[i, j] = 2radius * asin(sqrt(clamp(a, 0.0, 1.0)))
    end
    return result
end

function _adaptive_bisquare(distances::Vector{Float64}, neighbors::Int)
    finite_distances = sort(filter(isfinite, distances))
    isempty(finite_distances) && return zeros(Float64, length(distances))
    k = clamp(neighbors, 2, length(finite_distances))
    bandwidth = finite_distances[k]
    bandwidth > 0 || (bandwidth = maximum(finite_distances))
    bandwidth > 0 || return Float64.(distances .== 0)
    return [isfinite(d) && d < bandwidth ? (1 - (d / bandwidth)^2)^2 : 0.0 for d in distances]
end

function _gwr_smoothers(
    X::Matrix{Float64}, distances::Matrix{Float64}, neighbors::Int;
    ridge::Float64=1e-8, exclude_self::Bool=true,
)
    n, p = size(X)
    size(distances) == (n, n) || throw(DimensionMismatch("GWR smoother requires square distances"))
    smoothers = [zeros(Float64, n, n) for _ in 1:p]
    for target in 1:n
        d = copy(@view distances[:, target])
        exclude_self && (d[target] = Inf)
        w = _adaptive_bisquare(d, neighbors)
        valid = w .> 0
        count(valid) >= p + 1 || throw(ArgumentError("insufficient local observations for bandwidth $neighbors"))
        Xv = X[valid, :]
        wv = w[valid]
        A = Xv' * (wv .* Xv) + ridge * I
        B = A \ (Xv' .* wv')
        indices = findall(valid)
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
    exceed = zeros(Int, length(groups))
    for _ in 1:permutations
        permutation = randperm(rng, length(y))
        Xpermuted = X[permutation, :]
        ypermuted = y[permutation]
        permuted_smoothers = _gwr_smoothers(
            Xpermuted, distances, best.bandwidth; ridge, exclude_self=false,
        )
        for group_index in eachindex(groups)
            statistic = sum(var(permuted_smoothers[index] * ypermuted)
                for index in group_indices[group_index])
            exceed[group_index] += statistic >= observed_statistics[group_index]
        end
    end
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

function _local_hat(
    Xtrain::Matrix{Float64}, Xtarget::Matrix{Float64}, distances::Matrix{Float64},
    neighbors::Int; ridge::Float64=1e-8, exclude_self::Bool=false,
)
    ntrain, p = size(Xtrain)
    ntarget = size(Xtarget, 1)
    size(Xtarget, 2) == p || throw(DimensionMismatch("local design columns differ"))
    size(distances) == (ntrain, ntarget) || throw(DimensionMismatch("distance dimensions differ"))
    exclude_self && ntrain != ntarget &&
        throw(DimensionMismatch("exclude_self requires matching train and target rows"))
    H = zeros(Float64, ntarget, ntrain)
    for target in 1:ntarget
        d = copy(@view distances[:, target])
        exclude_self && (d[target] = Inf)
        w = _adaptive_bisquare(d, neighbors)
        valid = w .> 0
        count(valid) >= p + 1 || continue
        Xv = Xtrain[valid, :]
        wv = w[valid]
        A = Xv' * (wv .* Xv) + ridge * I
        B = A \ (Xv' .* wv')
        H[target, findall(valid)] = vec(Xtarget[target, :]' * B)
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
    distances::Matrix{Float64}, bandwidth::Int;
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_hat = _local_hat(Xlocal, Xlocal, distances, bandwidth; ridge)
    global_hat = _global_projection(Xglobal; ridge)
    result = _backfit_components(
        y, [local_hat], global_hat; tolerance, max_iterations,
    )
    return merge(result, (; local_hat, global_hat))
end

function _mixed_predict_complete(
    Xlocal_train::Matrix{Float64}, Xglobal_train::Matrix{Float64}, y::Vector{Float64},
    train_lonlat::Matrix{Float64}, Xlocal_target::Matrix{Float64},
    Xglobal_target::Matrix{Float64}, target_lonlat::Matrix{Float64}, bandwidth::Int;
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
)
    train_distances = _haversine_matrix(train_lonlat, train_lonlat)
    fitted = _mixed_fit_complete(
        Xlocal_train, Xglobal_train, y, train_distances, bandwidth;
        ridge, tolerance, max_iterations,
    )
    global_beta = isempty(Xglobal_train) ? Float64[] :
        (Xglobal_train' * Xglobal_train + ridge * I) \ (Xglobal_train' * (y - fitted.local_components[1]))
    partial = y - (isempty(Xglobal_train) ? zeros(length(y)) : Xglobal_train * global_beta)
    target_hat = _local_hat(
        Xlocal_train, Xlocal_target,
        _haversine_matrix(train_lonlat, target_lonlat), bandwidth; ridge,
    )
    prediction = target_hat * partial
    isempty(Xglobal_target) || (prediction .+= Xglobal_target * global_beta)
    return prediction, fitted.converged, fitted.iterations
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
    lonlat::Matrix{Float64}, candidates::Vector{Int}; ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    distances = _haversine_matrix(lonlat, lonlat)
    global_hat = _global_projection(Xglobal; ridge)
    rows = NamedTuple[]
    best_bandwidth = 0
    best_rmse = Inf
    for bandwidth in candidates
        bandwidth < length(y) || continue
        try
            local_hat = _local_hat(Xlocal, Xlocal, distances, bandwidth; ridge)
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
    distances::Matrix{Float64}, bandwidths::Vector{Int}; ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    length(local_groups) == length(bandwidths) || throw(DimensionMismatch("bandwidth count differs"))
    hats = [_local_hat(X, X, distances, bandwidth; ridge)
        for (X, bandwidth) in zip(local_groups, bandwidths)]
    global_hat = _global_projection(Xglobal; ridge)
    result = _backfit_components(y, hats, global_hat; tolerance, max_iterations)
    return merge(result, (; hats, global_hat))
end

function select_multiscale_bandwidths(
    local_groups::Vector{Matrix{Float64}}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    lonlat::Matrix{Float64}, candidates::Vector{Int}; ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    isempty(local_groups) && return Int[], DataFrame(), true
    usable = filter(<(length(y)), candidates)
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
            best_bw, best_rmse, best_hat = 0, Inf, nothing
            for bandwidth in usable
                try
                    hat = _local_hat(
                        local_groups[group_index], local_groups[group_index], distances,
                        bandwidth; ridge,
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

function _multiscale_predict_complete(
    local_train::Vector{Matrix{Float64}}, Xglobal_train::Matrix{Float64}, y::Vector{Float64},
    train_lonlat::Matrix{Float64}, local_target::Vector{Matrix{Float64}},
    Xglobal_target::Matrix{Float64}, target_lonlat::Matrix{Float64}, bandwidths::Vector{Int};
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
)
    train_distances = _haversine_matrix(train_lonlat, train_lonlat)
    fitted = _multiscale_fit_complete(
        local_train, Xglobal_train, y, train_distances, bandwidths;
        ridge, tolerance, max_iterations,
    )
    fitted.converged || return fill(NaN, size(target_lonlat, 1)), false, fitted.iterations
    local_sum = reduce(+, fitted.local_components)
    global_beta = isempty(Xglobal_train) ? Float64[] :
        (Xglobal_train' * Xglobal_train + ridge * I) \ (Xglobal_train' * (y - local_sum))
    global_train = isempty(Xglobal_train) ? zeros(length(y)) : Xglobal_train * global_beta
    prediction = isempty(Xglobal_target) ? zeros(size(target_lonlat, 1)) : Xglobal_target * global_beta
    target_distances = _haversine_matrix(train_lonlat, target_lonlat)
    for group_index in eachindex(local_train)
        partial = y - global_train
        for other in eachindex(fitted.local_components)
            other == group_index || (partial .-= fitted.local_components[other])
        end
        target_hat = _local_hat(
            local_train[group_index], local_target[group_index], target_distances,
            bandwidths[group_index]; ridge,
        )
        prediction .+= target_hat * partial
    end
    return prediction, true, fitted.iterations
end

function mixed_gwr_predict(
    Xlocal_train::Matrix{Float64}, Xglobal_train::Matrix{Float64}, Ytrain::Matrix{Float64},
    train_lonlat::Matrix{Float64}, Xlocal_target::Matrix{Float64},
    Xglobal_target::Matrix{Float64}, target_lonlat::Matrix{Float64}, bandwidth::Int;
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
    exclude_self::Bool=false,
)
    if exclude_self
        size(train_lonlat) == size(target_lonlat) ||
            throw(DimensionMismatch("exclude_self requires matching train and target coordinates"))
        train_lonlat == target_lonlat ||
            throw(ArgumentError("exclude_self requires train and target coordinates in the same order"))
    end
    prediction = fill(NaN, size(Xlocal_target, 1), size(Ytrain, 2))
    converged = trues(size(Ytrain, 2))
    cache = Dict{String,NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        if count(valid) <= size(Xlocal_train, 2) + size(Xglobal_train, 2) + 2
            converged[time] = false
            continue
        end
        indices = findall(valid)
        key = join(indices, ',')
        operators = get!(cache, key) do
            local_train_valid = Xlocal_train[indices, :]
            global_train_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(Xlocal_target, 1))
            local_target_valid = Xlocal_target[target_indices, :]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            adjusted_bandwidth = min(bandwidth, length(indices) - 1)
            local_hat = _local_hat(
                local_train_valid, local_train_valid,
                _haversine_matrix(lonlat_valid, lonlat_valid), adjusted_bandwidth; ridge,
            )
            target_hat = _local_hat(
                local_train_valid, local_target_valid,
                _haversine_matrix(lonlat_valid, target_lonlat_valid), adjusted_bandwidth;
                ridge, exclude_self,
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
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Int}; ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200, exclude_self::Bool=false,
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
    cache = Dict{String,NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        min_required = sum(size(X, 2) for X in local_train) + size(Xglobal_train, 2) + 2
        if count(valid) <= min_required
            converged[time] = false
            continue
        end
        indices = findall(valid)
        key = join(indices, ',')
        operators = get!(cache, key) do
            local_valid = [X[indices, :] for X in local_train]
            global_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(target_lonlat, 1))
            local_target_valid = [X[target_indices, :] for X in local_target]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            adjusted_bandwidths = min.(bandwidths, length(indices) - 1)
            train_distances = _haversine_matrix(lonlat_valid, lonlat_valid)
            train_hats = [_local_hat(X, X, train_distances, bw; ridge)
                for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            global_hat = _global_projection(global_valid; ridge)
            target_hats = if exclude_self
                [_local_hat(X, X, train_distances, bw; ridge, exclude_self=true)
                    for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            else
                target_distances = _haversine_matrix(lonlat_valid, target_lonlat_valid)
                [_local_hat(X, Xt, target_distances, bw; ridge)
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

function _metric_row(
    product::String, fold::Int, method::String, quantity::String, stratum::String,
    truth::Matrix{Float64}, prediction::Matrix{Float64}, mask::BitMatrix,
    eligible_n::Int,
)
    n = count(mask)
    n == 0 && return (; product, fold, method, quantity, stratum, n=0,
        coverage=0.0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN)
    actual = truth[mask]
    estimated = prediction[mask]
    error = estimated - actual
    return (;
        product, fold, method, quantity, stratum, n,
        coverage=eligible_n > 0 ? n / eligible_n : 0.0, RMSE=sqrt(mean(abs2, error)),
        MAE=mean(abs, error), Bias=mean(error),
        r=n > 1 && std(actual) > 0 && std(estimated) > 0 ? cor(actual, estimated) : NaN,
    )
end

function _append_metrics!(
    rows::Vector{NamedTuple}, product::String, fold::Int, method::String,
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, residual_prediction::Matrix{Float64},
    threshold::Float64,
)
    residual_truth = Yobs - Ysat
    corrected = max.(Ysat + residual_prediction, 0.0)
    eligible = isfinite.(Yobs) .& isfinite.(Ysat)
    strata = (
        ("all", eligible),
        ("wet", eligible .& (Yobs .>= threshold)),
        ("no_rain", eligible .& (Yobs .< 0.1)),
        ("light", eligible .& (Yobs .>= 0.1) .& (Yobs .< 2.5)),
        ("moderate", eligible .& (Yobs .>= 2.5) .& (Yobs .< 8.0)),
        ("heavy", eligible .& (Yobs .>= 8.0)),
    )
    for (name, eligible_mask) in strata
        bitmask = BitMatrix(eligible_mask .& isfinite.(residual_prediction))
        eligible_n = count(eligible_mask)
        push!(rows, _metric_row(
            product, fold, method, "residual", name,
            residual_truth, residual_prediction, bitmask, eligible_n,
        ))
        push!(rows, _metric_row(
            product, fold, method, "corrected_precipitation", name,
            Yobs, corrected, bitmask, eligible_n,
        ))
    end
    return rows
end

function _write_wide(path::String, times::Vector{DateTime}, ids::Vector{String}, values::Matrix{Float64})
    table = DataFrame(time=Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS"))
    for (index, id) in enumerate(ids)
        table[!, Symbol(id)] = values[index, :]
    end
    CSV.write(path, table)
end

function _consensus_table(role_tables::Vector{DataFrame})
    rows = NamedTuple[]
    products = [nrow(table) == 0 ? "" : String(table.product[1]) for table in role_tables]
    for group in TERRAIN_GROUPS
        selected_products = String[]
        roles = String[]
        for (product, table) in zip(products, role_tables)
            matches = filter(:variable_group => ==(group), table)
            nrow(matches) == 1 || continue
            push!(selected_products, product)
            push!(roles, String(matches.role[1]))
        end
        consensus_selected = length(selected_products) >= 2
        local_count = count(==("local"), roles)
        global_count = count(==("global"), roles)
        consensus_role = local_count >= 2 ? "local" : global_count >= 2 ? "global" :
            consensus_selected ? "product_specific" : "not_selected"
        push!(rows, (;
            variable_group=group, selected_product_count=length(selected_products),
            selected_products=join(selected_products, ","), local_product_count=local_count,
            global_product_count=global_count, consensus_selected, consensus_role,
        ))
    end
    return DataFrame(rows)
end

"""Run the standalone two-step DEM screening and spatial-validation experiment."""
function run_dem_experiment(
    cfg::DEMExperimentConfig, products::Vector{String}, ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64},
    satellite::Dict{String,Matrix{Float64}}, terrain::DataFrame,
    lonlat::Matrix{Float64},
)
    mkpath(cfg.outdir)
    length(ids) == nrow(terrain) == size(lonlat, 1) == size(Yobs, 1) ||
        throw(DimensionMismatch("station dimensions are not aligned"))
    allunique(ids) || throw(ArgumentError("duplicate station IDs"))
    terrain_ids = String.(terrain.station_id)
    terrain_ids == ids || throw(ArgumentError("terrain rows must follow common station order"))
    all(isfinite, Matrix{Float64}(terrain[:, [:elevation_m, :slope_deg, :aspect_sin, :aspect_cos]])) ||
        throw(ArgumentError("non-finite terrain values"))
    2 <= cfg.k <= length(ids) || throw(ArgumentError("invalid spatial fold count"))

    qc = DataFrame(
        key=["station_count", "time_count", "products", "wet_threshold", "min_wet_hours",
            "screen_permutations", "spatial_permutations", "selection_rule"],
        value=[string(length(ids)), string(length(times)), join(products, ","),
            string(cfg.wet_threshold), string(cfg.min_wet_hours),
            string(cfg.screen_permutations), string(cfg.spatial_permutations),
            "correlation_and_VIF_then_Monte_Carlo_spatial_variability"],
    )
    CSV.write(joinpath(cfg.outdir, "data_qc.csv"), qc)

    all_screen = DataFrame[]
    all_vif = DataFrame[]
    all_spatial = DataFrame[]
    all_spatial_scans = DataFrame[]
    all_monthly = NamedTuple[]
    role_tables = DataFrame[]
    responses = Dict{String,Vector{Float64}}()
    counts_by_product = Dict{String,Vector{Int}}()

    for product in products
        Ysat = satellite[product]
        response, counts = _mean_wet_residual(
            Yobs, Ysat; threshold=cfg.wet_threshold, min_hours=cfg.min_wet_hours,
        )
        responses[product] = response
        counts_by_product[product] = counts
        screen, vif, selected = terrain_screen(
            terrain, response; product, permutations=cfg.screen_permutations,
            q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
            seed=cfg.seed + Int(sum(codeunits(product))),
        )
        push!(all_screen, screen)
        push!(all_vif, vif)
        append!(all_monthly, _monthly_rows(product, times, terrain, Yobs, Ysat, cfg))
        if isempty(selected)
            push!(role_tables, DataFrame(
                product=String[], variable_group=String[], bandwidth=Int[],
                variability_statistic=Float64[], pvalue=Float64[], qvalue=Float64[], role=String[],
            ))
            continue
        end
        spatial, scan = spatial_variability_test(
            terrain, lonlat, response, selected; product,
            bandwidth_candidates=cfg.bandwidth_candidates,
            permutations=cfg.spatial_permutations, q_threshold=cfg.q_threshold,
            seed=cfg.seed + 10_000 + Int(sum(codeunits(product))), ridge=cfg.ridge,
        )
        push!(all_spatial, spatial)
        push!(all_spatial_scans, scan)
        push!(role_tables, spatial)
    end

    screen_table = isempty(all_screen) ? DataFrame() : vcat(all_screen...; cols=:union)
    vif_table = isempty(all_vif) ? DataFrame() : vcat(all_vif...; cols=:union)
    spatial_table = isempty(all_spatial) ? DataFrame() : vcat(all_spatial...; cols=:union)
    spatial_scan = isempty(all_spatial_scans) ? DataFrame() : vcat(all_spatial_scans...; cols=:union)
    CSV.write(joinpath(cfg.outdir, "correlation_screen.csv"), screen_table)
    CSV.write(joinpath(cfg.outdir, "vif.csv"), vif_table)
    CSV.write(joinpath(cfg.outdir, "monthly_correlation_direction.csv"), DataFrame(all_monthly))
    CSV.write(joinpath(cfg.outdir, "spatial_variability.csv"), spatial_table)
    CSV.write(joinpath(cfg.outdir, "spatial_bandwidth_scan.csv"), spatial_scan)
    CSV.write(joinpath(cfg.outdir, "cross_product_consensus.csv"), _consensus_table(role_tables))

    count_table = DataFrame(station_id=ids)
    for product in products
        count_table[!, Symbol(lowercase(product), "_wet_hours")] = counts_by_product[product]
        count_table[!, Symbol(lowercase(product), "_mean_wet_residual")] = responses[product]
    end
    CSV.write(joinpath(cfg.outdir, "station_residual_summary.csv"), count_table)

    folds = _balanced_spatial_folds(ids, lonlat; k=cfg.k, seed=cfg.seed)
    fold_map = zeros(Int, length(ids))
    for fold in 1:cfg.k
        fold_map[folds[fold]] .= fold
    end
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), DataFrame(station_id=ids, fold=fold_map))

    metric_rows = NamedTuple[]
    tuning_rows = NamedTuple[]
    status_rows = NamedTuple[]
    method_names = ["residual_gwr", "mixed_gwr", "multiscale_gwr"]
    oof = Dict((product, method) => fill(NaN, size(Yobs))
        for product in products for method in method_names)

    fold_role_rows = NamedTuple[]
    for (product_index, product) in enumerate(products)
        Ysat = satellite[product]
        for fold in 1:cfg.k
            val_idx = folds[fold]
            train_idx = setdiff(1:length(ids), val_idx)

            # Terrain-role selection (which variables are used, and whether each is
            # spatially-varying) is recomputed from training-fold stations only, so a
            # fold's held-out validation stations never influence its own model structure.
            _, _, fold_selected_groups = terrain_screen(
                terrain[train_idx, :], responses[product][train_idx]; product,
                permutations=cfg.screen_permutations, q_threshold=cfg.q_threshold,
                vif_threshold=cfg.vif_threshold,
                seed=cfg.seed + Int(sum(codeunits(product))) + 100 * fold,
            )
            role_map = Dict{String,String}()
            if !isempty(fold_selected_groups)
                fold_spatial, _ = spatial_variability_test(
                    terrain[train_idx, :], lonlat[train_idx, :], responses[product][train_idx],
                    fold_selected_groups; product, bandwidth_candidates=cfg.bandwidth_candidates,
                    permutations=cfg.spatial_permutations, q_threshold=cfg.q_threshold,
                    seed=cfg.seed + 10_000 + Int(sum(codeunits(product))) + 100 * fold,
                    ridge=cfg.ridge,
                )
                role_map = Dict(String(row.variable_group) => String(row.role) for row in eachrow(fold_spatial))
                for row in eachrow(fold_spatial)
                    push!(fold_role_rows, (; product, fold, variable_group=row.variable_group,
                        bandwidth=row.bandwidth, variability_statistic=row.variability_statistic,
                        pvalue=row.pvalue, qvalue=row.qvalue, role=row.role))
                end
            end
            can_fit_dem = !isempty(role_map) && all(!=("uncertain"), values(role_map))

            designs = _designs(
                terrain[train_idx, :], terrain[val_idx, :], lonlat[train_idx, :],
                lonlat[val_idx, :], can_fit_dem ? role_map : Dict{String,String}(),
            )
            aggregate = responses[product][train_idx]
            valid_aggregate = isfinite.(aggregate)
            train_count = count(valid_aggregate)
            usable_candidates = filter(<(train_count), cfg.bandwidth_candidates)
            if train_count < 12 || isempty(usable_candidates)
                push!(status_rows, (; product, fold, method="all", status="failed",
                    error="insufficient aggregate training stations"))
                continue
            end
            aggregate_lonlat = lonlat[train_idx[valid_aggregate], :]
            baseline_X = designs.mixed_local_train[valid_aggregate, 1:3]
            baseline_bw, baseline_scan = select_mixed_bandwidth(
                baseline_X, zeros(Float64, train_count, 0), aggregate[valid_aggregate],
                aggregate_lonlat, usable_candidates; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(baseline_scan)
                push!(tuning_rows, (; product, fold, method="residual_gwr", group="all_local",
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.bandwidth == baseline_bw,
                    status=row.status, error=row.error))
            end
            baseline_prediction, baseline_converged = mixed_gwr_predict(
                designs.mixed_local_train[:, 1:3], zeros(Float64, length(train_idx), 0),
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.mixed_local_target[:, 1:3], zeros(Float64, length(val_idx), 0),
                lonlat[val_idx, :], baseline_bw; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            oof[(product, "residual_gwr")][val_idx, :] = baseline_prediction
            _append_metrics!(metric_rows, product, fold, "residual_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], baseline_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="residual_gwr",
                status=all(baseline_converged) ? "success" : "partial",
                error=all(baseline_converged) ? "" : "one or more hourly fits failed"))

            can_fit_dem || continue
            mixed_bw, mixed_scan = select_mixed_bandwidth(
                designs.mixed_local_train[valid_aggregate, :],
                designs.global_train[valid_aggregate, :], aggregate[valid_aggregate],
                aggregate_lonlat, usable_candidates; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(mixed_scan)
                push!(tuning_rows, (; product, fold, method="mixed_gwr", group="all_local",
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.bandwidth == mixed_bw,
                    status=row.status, error=row.error))
            end
            mixed_prediction, mixed_converged = mixed_gwr_predict(
                designs.mixed_local_train, designs.global_train,
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.mixed_local_target, designs.global_target, lonlat[val_idx, :], mixed_bw;
                ridge=cfg.ridge, tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            oof[(product, "mixed_gwr")][val_idx, :] = mixed_prediction
            _append_metrics!(metric_rows, product, fold, "mixed_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], mixed_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="mixed_gwr",
                status=all(mixed_converged) ? "success" : "partial",
                error=all(mixed_converged) ? "" : "one or more hourly fits failed"))

            local_valid = [X[valid_aggregate, :] for X in designs.multiscale_train]
            multiscale_bw, multiscale_scan, bandwidth_converged = select_multiscale_bandwidths(
                local_valid, designs.global_train[valid_aggregate, :], aggregate[valid_aggregate],
                aggregate_lonlat, usable_candidates; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(multiscale_scan)
                group_name = designs.local_group_names[row.group_index]
                push!(tuning_rows, (; product, fold, method="multiscale_gwr", group=group_name,
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.selected,
                    status=row.status, error=row.error))
            end
            if !bandwidth_converged
                push!(status_rows, (; product, fold, method="multiscale_gwr", status="failed",
                    error="bandwidth backfitting did not converge"))
                continue
            end
            multiscale_prediction, multiscale_converged = multiscale_gwr_predict(
                designs.multiscale_train, designs.global_train,
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.multiscale_target, designs.global_target, lonlat[val_idx, :], multiscale_bw;
                ridge=cfg.ridge, tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            oof[(product, "multiscale_gwr")][val_idx, :] = multiscale_prediction
            _append_metrics!(metric_rows, product, fold, "multiscale_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], multiscale_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="multiscale_gwr",
                status=all(multiscale_converged) ? "success" : "partial",
                error=all(multiscale_converged) ? "" : "one or more hourly fits failed"))
        end
    end

    metrics = DataFrame(metric_rows)
    tuning = DataFrame(tuning_rows)
    status = DataFrame(status_rows)
    fold_roles = DataFrame(fold_role_rows)
    CSV.write(joinpath(cfg.outdir, "validation_metrics.csv"), metrics)
    CSV.write(joinpath(cfg.outdir, "bandwidth_scan.csv"), tuning)
    CSV.write(joinpath(cfg.outdir, "run_status.csv"), status)
    CSV.write(joinpath(cfg.outdir, "spatial_variability_by_fold.csv"), fold_roles)
    for product in products
        product_dir = joinpath(cfg.outdir, lowercase(product))
        mkpath(product_dir)
        for method in method_names
            values = oof[(product, method)]
            any(isfinite, values) || continue
            _write_wide(joinpath(product_dir, "oof_residual_$(method).csv"), times, ids, values)
            corrected = max.(satellite[product] + values, 0.0)
            _write_wide(joinpath(product_dir, "oof_corrected_$(method).csv"), times, ids, corrected)
        end
    end
    return (; screen=screen_table, vif=vif_table, spatial=spatial_table,
        consensus=_consensus_table(role_tables), metrics, tuning, status, fold_roles)
end

end
