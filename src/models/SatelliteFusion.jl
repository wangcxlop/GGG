"""
Combine several satellite precipitation products into one anchor.

The benchmark corrects one product at a time, and F4 reads the GWR family's loss as a property of
the correction: a false alarm is handed to the model as a fixed offset it cannot argue with. That
is true, and it makes the anchor itself the thing worth changing. The three products are not
redundant estimates of the same error - measured over the 3,080,258 cells where the gauge and all
three are finite, anchor RMSE is 1.317 (FY4B), 0.965 (GPM) and 1.223 (GSMaP), while a
least-squares combination of them is 0.888, better than every GWR-family method's own result on
FY4B and better than `tps` and direct `gwr`.

Two variants, and the difference between them is not a detail:

- `:ols` fits `y_obs ~ 1 + p1 + p2 + ...` and is the better *pooled* anchor. Its coefficients sum
  to well under one, so it shrinks - which is what a squared-error fit should do, and which costs
  it the heavy stratum (−4.3% against the traditional baselines as a bare anchor).
- `:mean` weights the products equally, fits nothing, and does the reverse: worse pooled (0.988),
  better on heavy. Carrying it, `mgwr` clears the pre-registered ≥5% heavy gate on all three
  products, which no shipped configuration does.

Both are offered because that trade-off cannot be resolved by tuning. A shrinkage ladder selected
on inner-split RMSE would collapse to the pooled optimum every time and could never reach the
heavy gate, so the choice is exposed as two products rather than hidden in a selection step.

Leakage is the caller's responsibility and the whole point of `station_rows`: the benchmark fits
these coefficients on a fold's *training* stations only and then applies them everywhere, so the
held-out stations' gauge values never reach the anchor their own predictions are built on.
"""
module SatelliteFusion

using Dates: Hour
using LinearAlgebra

export FUSION_VARIANTS, fit_satellite_fusion, apply_satellite_fusion
export equal_weight_coefficients, fusion_coefficients
export GROUPED_FUSION_VARIANTS, fusion_features, fit_grouped_fusion, apply_grouped_fusion

"""The fusion variants a caller may ask for. See the module docstring for why both exist."""
const FUSION_VARIANTS = (:ols, :mean)

"""
Fusion variants whose design is wider than the products at one hour, fitted per group of cells.

`:ols_lagnbr` regresses the gauge on each product at `t`, `t-1` and `t+1` and on each product's
mean over the nearest other stations, with one coefficient set per agreement-envelope band. The
lags absorb a product's timing error and the neighbour means its position error - the two things a
single-hour, single-point fit cannot see. Screened on the benchmark's own balanced_spatial folds it
cut the bare anchor's RMSE from 0.888 to 0.876, with the inner-split and in-sample scores within
0.003 of that, so the extra coefficients are not what buys it. It does not rescue the heavy
stratum: a squared-error anchor keeps a heavy bias near -12 mm/h whatever it is given.

Kept apart from `FUSION_VARIANTS` so `--fused-anchor` still means exactly the two products it
always meant.
"""
const GROUPED_FUSION_VARIANTS = (:ols_lagnbr,)

"""
Equal weights over `n` products, with a zero intercept.

Same layout as a fitted coefficient vector - intercept first, then one slope per product - so the
two are interchangeable everywhere downstream.
"""
equal_weight_coefficients(n::Int) = vcat(0.0, fill(1 / n, n))

"""
Least-squares fusion coefficients from one set of station rows: `y_obs ~ 1 + p1 + p2 + ...`.

Accumulated as normal equations rather than as a tall design matrix. The training half of a full
run is ~2.5M cells, and a `(p+1)×(p+1)` Gram matrix says the same thing without allocating 80 MB
per fold - which matters because this is refitted once per fold per derived product.

Only cells where the gauge *and* every product are finite contribute, because a partial row would
silently change which design the coefficients belong to.

Returns `(nothing, used)` when no cell is usable or the system is singular, which the caller is
expected to report as a fall back to `equal_weight_coefficients` rather than treat as a fit.
"""
function fit_satellite_fusion(
    y_obs::AbstractMatrix{Float64}, sources::Vector{<:AbstractMatrix{Float64}},
    station_rows::AbstractVector{Int},
)
    isempty(sources) && throw(ArgumentError("satellite fusion needs at least one source product"))
    for source in sources
        size(source) == size(y_obs) || throw(DimensionMismatch(
            "every source product must have the shape of the observation matrix",
        ))
    end
    p = length(sources) + 1
    gram = zeros(Float64, p, p)
    rhs = zeros(Float64, p)
    row = ones(Float64, p)
    used = 0
    @inbounds for station in station_rows, time in axes(y_obs, 2)
        observed = y_obs[station, time]
        isnan(observed) && continue
        usable = true
        for (index, source) in enumerate(sources)
            value = source[station, time]
            if isnan(value)
                usable = false
                break
            end
            row[index + 1] = value
        end
        usable || continue
        used += 1
        for i in 1:p
            rhs[i] += row[i] * observed
            for j in i:p
                gram[i, j] += row[i] * row[j]
            end
        end
    end
    used == 0 && return nothing, 0
    @inbounds for i in 1:p, j in 1:(i - 1)
        gram[i, j] = gram[j, i]
    end
    coefficients = try
        cholesky(Symmetric(gram)) \ rhs
    catch
        return nothing, used
    end
    all(isfinite, coefficients) || return nothing, used
    return coefficients, used
end

"""
Coefficients for one variant, and whether the fit fell back.

`:mean` fits nothing. `:ols` falls back to the equal-weight vector when its normal equations are
singular, so a fold can never end up without an anchor - but the fallback is reported rather than
silent, because a fold whose fusion did not fit is not the same experiment as one whose did.
"""
function fusion_coefficients(
    variant::Symbol, y_obs::AbstractMatrix{Float64},
    sources::Vector{<:AbstractMatrix{Float64}}, station_rows::AbstractVector{Int},
)
    variant in FUSION_VARIANTS ||
        throw(ArgumentError("unknown fusion variant $variant; expected one of $FUSION_VARIANTS"))
    equal = equal_weight_coefficients(length(sources))
    variant === :mean && return equal, 0, false
    coefficients, used = fit_satellite_fusion(y_obs, sources, station_rows)
    coefficients === nothing && return equal, used, true
    return coefficients, used, false
end

"""
Apply fusion coefficients over the whole grid, clipped at zero.

A cell where any source is NaN comes out NaN. That is the honest answer rather than a convenience:
the fusion is a function of all of the products and there is no partial version of it. Those cells
drop out of the derived product's evaluation mask, which is a real difference in cell population
between it and the products it is built from, and one the run records.

The clip matches what every method's prediction already gets (`predict_selected`'s
`max.(..., 0.0)`), so the anchor and the predictions built on it live on the same half-line.
"""
function apply_satellite_fusion(
    sources::Vector{<:AbstractMatrix{Float64}}, coefficients::AbstractVector{Float64},
)
    isempty(sources) && throw(ArgumentError("satellite fusion needs at least one source product"))
    length(coefficients) == length(sources) + 1 || throw(DimensionMismatch(
        "expected one intercept plus one coefficient per source product",
    ))
    out = fill(NaN, size(first(sources)))
    @inbounds for index in eachindex(out)
        total = coefficients[1]
        usable = true
        for (position, source) in enumerate(sources)
            value = source[index]
            if isnan(value)
                usable = false
                break
            end
            total += coefficients[position + 1] * value
        end
        usable && (out[index] = max(total, 0.0))
    end
    return out
end

"""
The design of `:ols_lagnbr`, one matrix per term, each of the sources' shape.

Order: every source at `t`, then every source at `t-1`, then at `t+1`, then every source's mean
over `neighbours[s][1:k]`. `neighbours` holds each station's other stations nearest first; it is
passed in so this module needs no geometry.

A lag whose adjacent hour is not `step` away on `times` (a gap in the common grid) or is NaN
falls back to the hour's own value, and a neighbour mean with no finite neighbour falls back to the
station's own value. Both only ever read satellite values, so a held-out gauge cannot reach them.
"""
function fusion_features(
    sources::Vector{<:AbstractMatrix{Float64}}, times::AbstractVector,
    neighbours::AbstractVector{<:AbstractVector{Int}}; k::Int=8, step=Hour(1),
)
    isempty(sources) && throw(ArgumentError("satellite fusion needs at least one source product"))
    n_station, n_time = size(first(sources))
    length(times) == n_time || throw(DimensionMismatch("times must match the sources' columns"))
    length(neighbours) == n_station ||
        throw(DimensionMismatch("need one neighbour list per station"))
    function lagged(source, shift)
        out = Matrix{Float64}(source)
        @inbounds for time in 1:n_time
            other = time + shift
            (1 <= other <= n_time && times[other] - times[time] == shift * step) || continue
            for station in 1:n_station
                value = source[station, other]
                isnan(value) || (out[station, time] = value)
            end
        end
        return out
    end
    function neighbour_mean(source)
        out = Matrix{Float64}(source)
        @inbounds for station in 1:n_station
            nearest = neighbours[station][1:min(k, length(neighbours[station]))]
            for time in 1:n_time
                total = 0.0
                count = 0
                for other in nearest
                    value = source[other, time]
                    isnan(value) || (total += value; count += 1)
                end
                count > 0 && (out[station, time] = total / count)
            end
        end
        return out
    end
    return vcat(
        [Matrix{Float64}(source) for source in sources],
        [lagged(source, -1) for source in sources],
        [lagged(source, 1) for source in sources],
        [neighbour_mean(source) for source in sources],
    )
end

"""
Per-group least squares `y_obs ~ 1 + features...` over `station_rows`, one column per group.

`groups` assigns every cell a group in `1:n_group`. The accumulation is `fit_satellite_fusion`'s,
once per group. The design is the cells where the gauge and every one of `sources` is finite -
the same cells MERGED_OLS fits on - rather than every cell the features happen to cover.

A group with fewer than `min_cells` cells, or a singular system, takes the pooled fit over all
groups and is flagged in `fell_back`; a pooled fit that is itself singular returns `nothing`.
"""
function fit_grouped_fusion(
    y_obs::AbstractMatrix{Float64}, sources::Vector{<:AbstractMatrix{Float64}},
    features::Vector{<:AbstractMatrix{Float64}}, groups::AbstractMatrix{Int}, n_group::Int,
    station_rows::AbstractVector{Int}; min_cells::Int=2000,
)
    for matrix in vcat(sources, features, [groups])
        size(matrix) == size(y_obs) || throw(DimensionMismatch(
            "every source, feature and the group matrix must have the observation's shape",
        ))
    end
    p = length(features) + 1
    gram = zeros(Float64, p, p, n_group)
    rhs = zeros(Float64, p, n_group)
    used = zeros(Int, n_group)
    row = ones(Float64, p)
    @inbounds for station in station_rows, time in axes(y_obs, 2)
        observed = y_obs[station, time]
        (isnan(observed) || any(source -> isnan(source[station, time]), sources)) && continue
        group = groups[station, time]
        1 <= group <= n_group || throw(ArgumentError("group $group is outside 1:$n_group"))
        for (index, feature) in enumerate(features)
            row[index + 1] = feature[station, time]
        end
        used[group] += 1
        for i in 1:p
            rhs[i, group] += row[i] * observed
            for j in i:p
                gram[i, j, group] += row[i] * row[j]
            end
        end
    end
    solve(g, r) = try
        solution = cholesky(Symmetric(g, :U)) \ r
        all(isfinite, solution) ? solution : nothing
    catch
        nothing
    end
    pooled = solve(dropdims(sum(gram; dims=3); dims=3), vec(sum(rhs; dims=2)))
    pooled === nothing && return nothing, used, trues(n_group)
    coefficients = zeros(Float64, p, n_group)
    fell_back = falses(n_group)
    for group in 1:n_group
        fitted = used[group] >= min_cells ? solve(gram[:, :, group], rhs[:, group]) : nothing
        fell_back[group] = fitted === nothing
        coefficients[:, group] = something(fitted, pooled)
    end
    return coefficients, used, fell_back
end

"""
Apply per-group coefficients over the whole grid, clipped at zero.

NaN wherever any of `sources` is NaN - the same cells `apply_satellite_fusion` leaves NaN, so the
derived product covers exactly what MERGED_OLS covers.
"""
function apply_grouped_fusion(
    sources::Vector{<:AbstractMatrix{Float64}}, features::Vector{<:AbstractMatrix{Float64}},
    coefficients::AbstractMatrix{Float64}, groups::AbstractMatrix{Int},
)
    size(coefficients, 1) == length(features) + 1 || throw(DimensionMismatch(
        "expected one intercept plus one coefficient per feature",
    ))
    out = fill(NaN, size(first(sources)))
    @inbounds for index in eachindex(out)
        any(source -> isnan(source[index]), sources) && continue
        group = groups[index]
        total = coefficients[1, group]
        for (position, feature) in enumerate(features)
            total += coefficients[position + 1, group] * feature[index]
        end
        out[index] = max(total, 0.0)
    end
    return out
end

end # module
