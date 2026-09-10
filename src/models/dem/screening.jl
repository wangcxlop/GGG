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
