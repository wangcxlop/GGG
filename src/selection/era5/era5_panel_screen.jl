# One of four different spatial-fold builders in the project, and the simplest: a sort-and-slice
# along the longer axis. The others are `split_stations_balanced_spatial_kfold`
# (InterpolationBenchmarkFolds.jl - Hilbert-seeded capacity-constrained Lloyd),
# `DEMTerrainExperiment._balanced_spatial_folds` (farthest-point-seeded Lloyd) and
# `MGERPipeline.split_stations_spatial_block_kfold` (contiguous lon/lat blocks). They are genuinely
# different partitioning algorithms rather than copies, so they are not merged; but which one a
# given selection path uses is historical, not reasoned, and that is worth knowing before reading
# across their results.
function balanced_spatial_folds(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("lonlat must have n × 2 rows"))
    lat_mid = mean(lonlat[:, 2])
    lon_span = (maximum(lonlat[:, 1]) - minimum(lonlat[:, 1])) * cosd(lat_mid) * 111.32
    lat_span = (maximum(lonlat[:, 2]) - minimum(lonlat[:, 2])) * 110.57
    column = lon_span >= lat_span ? 1 : 2
    order = sortperm(1:n; by=i -> (lonlat[i, column], lonlat[i, 3 - column], station_ids[i]))
    assignment = zeros(Int, n)
    base, extra = divrem(n, k)
    left = 1
    for fold in 1:k
        count_in_fold = base + (fold <= extra ? 1 : 0)
        indices = order[left:(left + count_in_fold - 1)]
        assignment[indices] .= fold
        left += count_in_fold
    end
    return assignment
end

function _weighted_mean(x::Vector{Float64}, w::Vector{Float64})
    total = sum(w)
    total > 0 || return NaN
    return dot(x, w) / total
end

function _weighted_correlation(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64})
    length(x) >= 3 || return NaN
    mx, my = _weighted_mean(x, w), _weighted_mean(y, w)
    dx, dy = x .- mx, y .- my
    vx, vy = dot(w, dx .* dx), dot(w, dy .* dy)
    vx > 0 && vy > 0 || return NaN
    return dot(w, dx .* dy) / sqrt(vx * vy)
end

function _tie_ranks(values::Vector{Float64})
    order = sortperm(values)
    ranks = zeros(Float64, length(values))
    left = 1
    while left <= length(order)
        right = left
        while right < length(order) && values[order[right + 1]] == values[order[left]]
            right += 1
        end
        ranks[order[left:right]] .= (left + right) / 2
        left = right + 1
    end
    return ranks
end

function _flatten_balanced(x::Matrix{Float64}, y::Matrix{Float64})
    xs, ys, ws = Float64[], Float64[], Float64[]
    for i in axes(y, 1)
        valid = findall(j -> isfinite(x[i, j]) && isfinite(y[i, j]), axes(y, 2))
        isempty(valid) && continue
        weight = 1 / length(valid)
        for j in valid
            push!(xs, x[i, j]); push!(ys, y[i, j]); push!(ws, weight)
        end
    end
    return xs, ys, ws
end

"""Return a row permutation; indexing a panel with it preserves every station block."""
station_block_permutation(nstations::Int, rng::AbstractRNG) = randperm(rng, nstations)

function _station_pair_moments(x::Matrix{Float64}, y::Matrix{Float64})
    n = size(y, 1)
    size(x) == size(y) || throw(DimensionMismatch("panel matrices must match"))
    valid_x = [findall(isfinite, view(x, i, :)) for i in 1:n]
    valid_y = [findall(isfinite, view(y, i, :)) for i in 1:n]
    count_pair = zeros(Int, n, n)
    sum_x = zeros(n, n); sum_y = zeros(n, n)
    sum_xx = zeros(n, n); sum_yy = zeros(n, n); sum_xy = zeros(n, n)
    for source in 1:n, target in 1:n
        xs, ys = valid_x[source], valid_y[target]
        ix = 1; iy = 1
        while ix <= length(xs) && iy <= length(ys)
            tx, ty = xs[ix], ys[iy]
            if tx == ty
                xv, yv = x[source, tx], y[target, ty]
                count_pair[source, target] += 1
                sum_x[source, target] += xv; sum_y[source, target] += yv
                sum_xx[source, target] += xv^2; sum_yy[source, target] += yv^2
                sum_xy[source, target] += xv * yv
                ix += 1; iy += 1
            elseif tx < ty
                ix += 1
            else
                iy += 1
            end
        end
    end
    return (; count=count_pair, sum_x, sum_y, sum_xx, sum_yy, sum_xy)
end

function _correlation_from_pair_moments(moments, order::Vector{Int})
    total_x = 0.0; total_y = 0.0; total_xx = 0.0
    total_yy = 0.0; total_xy = 0.0; station_count = 0
    for target in eachindex(order)
        source = order[target]
        count_pair = moments.count[source, target]
        count_pair > 0 || continue
        scale = 1 / count_pair
        total_x += moments.sum_x[source, target] * scale
        total_y += moments.sum_y[source, target] * scale
        total_xx += moments.sum_xx[source, target] * scale
        total_yy += moments.sum_yy[source, target] * scale
        total_xy += moments.sum_xy[source, target] * scale
        station_count += 1
    end
    station_count >= 3 || return NaN
    mean_x, mean_y = total_x / station_count, total_y / station_count
    variance_x = total_xx / station_count - mean_x^2
    variance_y = total_yy / station_count - mean_y^2
    variance_x > 0 && variance_y > 0 || return NaN
    covariance = total_xy / station_count - mean_x * mean_y
    return covariance / sqrt(variance_x * variance_y)
end

function prepare_dynamic_panel(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, era5::AbstractDict,
    train_indices::Vector{Int}, times::Vector{DateTime}, cfg::ERA5SelectionConfig,
)
    size(Yobs) == size(Ysat) || throw(DimensionMismatch("observation and satellite matrices differ"))
    size(Yobs, 2) == length(times) || throw(DimensionMismatch("time dimension differs"))
    variables = collect(ERA5_VARIABLES)
    n, nt = length(train_indices), length(times)
    residual = Yobs[train_indices, :] .- Ysat[train_indices, :]
    raw = Dict(v => Float64.(era5[v][train_indices, :]) for v in variables)
    mask = falses(n, nt)
    for i in 1:n, t in 1:nt
        mask[i, t] = isfinite(Yobs[train_indices[i], t]) &&
            isfinite(Ysat[train_indices[i], t]) && Yobs[train_indices[i], t] >= cfg.wet_threshold &&
            all(isfinite(raw[v][i, t]) for v in variables)
    end
    previous_count = -1
    while count(mask) != previous_count
        previous_count = count(mask)
        eligible_station_now = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
        mask[.!eligible_station_now, :] .= false
        eligible_time_now = vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time
        mask[:, .!eligible_time_now] .= false
    end
    eligible_station = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours

    y = fill(NaN, n, nt)
    centered = Dict(v => fill(NaN, n, nt) for v in variables)
    for t in 1:nt
        indices = findall(mask[:, t])
        length(indices) >= cfg.min_stations_per_time || continue
        ymean = mean(residual[indices, t])
        y[indices, t] .= residual[indices, t] .- ymean
        for v in variables
            xmean = mean(raw[v][indices, t])
            centered[v][indices, t] .= raw[v][indices, t] .- xmean
        end
    end
    scales = DataFrame(variable=String[], mean=Float64[], scale=Float64[])
    for v in variables
        x, _, w = _flatten_balanced(centered[v], y)
        μ = _weighted_mean(x, w)
        σ = sqrt(_weighted_mean((x .- μ) .^ 2, w))
        isfinite(σ) && σ > 0 || (σ = NaN)
        if isfinite(σ)
            for index in eachindex(centered[v])
                isfinite(centered[v][index]) && (centered[v][index] = (centered[v][index] - μ) / σ)
            end
        end
        push!(scales, (String(v), μ, σ))
    end
    qc = DataFrame(
        training_station_count=[n], eligible_station_count=[count(eligible_station)],
        eligible_time_count=[count(vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time)],
        wet_station_hours=[count(mask)], min_wet_hours=[cfg.min_wet_hours],
        min_stations_per_time=[cfg.min_stations_per_time],
        status=[count(eligible_station) >= cfg.min_stations_per_time ? "ok" : "insufficient_stations"],
    )
    return (; y, x=centered, mask, eligible_station, scales, qc)
end

function _association(panel, variable::Symbol, permutations::Int, rng::AbstractRNG)
    x, y, w = _flatten_balanced(panel.x[variable], panel.y)
    moments = _station_pair_moments(panel.x[variable], panel.y)
    pearson = _correlation_from_pair_moments(moments, collect(1:size(panel.y, 1)))
    spearman = _weighted_correlation(_tie_ranks(x), _tie_ranks(y), w)
    if permutations <= 0 || !isfinite(pearson)
        return pearson, spearman, isfinite(pearson) ? NaN : 1.0
    end
    exceed = 0
    n = size(panel.y, 1)
    for _ in 1:permutations
        order = station_block_permutation(n, rng)
        value = _correlation_from_pair_moments(moments, order)
        exceed += isfinite(value) && abs(value) >= abs(pearson)
    end
    return pearson, spearman, (exceed + 1) / (permutations + 1)
end

function _weighted_vifs(panel, active::Vector{Symbol})
    length(active) <= 1 && return fill(1.0, length(active))
    rows = Vector{Vector{Float64}}()
    weights = Float64[]
    for i in axes(panel.y, 1)
        valid = [t for t in axes(panel.y, 2) if isfinite(panel.y[i, t]) &&
            all(isfinite(panel.x[v][i, t]) for v in active)]
        isempty(valid) && continue
        weight = 1 / length(valid)
        for t in valid
            push!(rows, [panel.x[v][i, t] for v in active]); push!(weights, weight)
        end
    end
    isempty(rows) && return fill(Inf, length(active))
    X = reduce(vcat, permutedims.(rows))
    sw = sqrt.(weights ./ sum(weights))
    result = fill(Inf, length(active))
    for j in eachindex(active)
        other = setdiff(eachindex(active), [j])
        design = hcat(ones(size(X, 1)), X[:, other])
        beta = (design .* sw) \ (X[:, j] .* sw)
        prediction = design * beta
        μ = _weighted_mean(X[:, j], weights)
        rss = sum(weights .* (X[:, j] .- prediction) .^ 2)
        tss = sum(weights .* (X[:, j] .- μ) .^ 2)
        r2 = tss > 0 ? clamp(1 - rss / tss, 0.0, 1.0) : 1.0
        result[j] = r2 < 1 ? 1 / (1 - r2) : Inf
    end
    return result
end

function dynamic_panel_screen(
    panel, times::Vector{DateTime}, cfg::ERA5SelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed),
)
    association = DataFrame(
        variable=String[], pearson=Float64[], spearman=Float64[],
        direction_stable=Bool[], pvalue=Float64[], qvalue=Float64[], selected_bh=Bool[],
    )
    raw = NamedTuple[]
    for variable in ERA5_VARIABLES
        pearson, spearman, pvalue = _association(
            panel, variable, cfg.association_permutations, rng,
        )
        stable = isfinite(pearson) && isfinite(spearman) && sign(pearson) == sign(spearman)
        push!(raw, (; variable, pearson, spearman, pvalue, stable))
    end
    qvalues = _bh_adjust([row.pvalue for row in raw])
    for (row, qvalue) in zip(raw, qvalues)
        push!(association, (
            String(row.variable), row.pearson, row.spearman, row.stable,
            row.pvalue, qvalue, row.stable && qvalue < cfg.q_threshold,
        ))
    end
    active = Symbol.(association.variable[association.selected_bh])
    vif = DataFrame(iteration=Int[], variable=String[], vif=Float64[], removed=Bool[], reason=String[])
    iteration = 1
    while !isempty(active)
        values = _weighted_vifs(panel, active)
        maxvif = maximum(values)
        if maxvif < cfg.vif_threshold
            for (variable, value) in zip(active, values)
                push!(vif, (iteration, String(variable), value, false, "retained"))
            end
            break
        end
        tied = findall(v -> isapprox(v, maxvif; rtol=1e-10, atol=1e-10) || (isinf(v) && isinf(maxvif)), values)
        qmap = Dict(Symbol(row.variable) => row.qvalue for row in eachrow(association))
        remove_index = sort(tied; by=j -> (qmap[active[j]], findfirst(==(active[j]), ERA5_VARIABLES)), rev=true)[1]
        for (j, (variable, value)) in enumerate(zip(active, values))
            push!(vif, (iteration, String(variable), value, j == remove_index,
                j == remove_index ? "vif_ge_threshold" : "pending"))
        end
        deleteat!(active, remove_index)
        iteration += 1
    end
    monthly = DataFrame(variable=String[], month=Int[], pearson=Float64[], spearman=Float64[],
        pearson_direction=String[], spearman_direction=String[])
    for variable in ERA5_VARIABLES, month_value in 6:9
        columns = findall(t -> month(t) == month_value, times)
        if isempty(columns)
            p, s = NaN, NaN
        else
            x, y, w = _flatten_balanced(panel.x[variable][:, columns], panel.y[:, columns])
            p = _weighted_correlation(x, y, w)
            s = _weighted_correlation(_tie_ranks(x), _tie_ranks(y), w)
        end
        direction(value) = !isfinite(value) ? "uncertain" : value > 0 ? "positive" :
            value < 0 ? "negative" : "zero"
        push!(monthly, (String(variable), month_value, p, s, direction(p), direction(s)))
    end
    return (; association, vif, monthly, selected=active)
end
