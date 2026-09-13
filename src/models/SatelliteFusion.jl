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

using LinearAlgebra

export FUSION_VARIANTS, fit_satellite_fusion, apply_satellite_fusion
export equal_weight_coefficients, fusion_coefficients

"""The fusion variants a caller may ask for. See the module docstring for why both exist."""
const FUSION_VARIANTS = (:ols, :mean)

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

end # module
