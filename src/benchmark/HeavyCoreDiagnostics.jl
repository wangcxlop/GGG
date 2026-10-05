"""
Why a satellite product underestimates the core of a heavy-rain event.

`HeavyRainEvents` scores an event across the whole gauge network. This module looks only at the
*core* gauges - those with the largest event totals - and splits each one's deficit
(gauge - satellite) into three parts that can each be measured:

- *representativeness*: a 0.1 deg pixel is an areal mean, so even a perfect product sits below a
  point maximum. Measured by the mean of the gauges within a radius (`upscaled_gauge_truth`).
- *displacement*: the product has the core, but in a neighbouring pixel. Measured by the largest
  product value within a radius (`neighbourhood_max`).
- *magnitude*: what is left - the product has nowhere near the areal-mean amount anywhere nearby.

`deficit_attribution` combines the three; `core_hourly_diagnostics` says whether the missing rain
is concentrated in the gauge's peak hours or spread over the event; `neighbour_consistency` checks
the core gauge against its nearest neighbours, so that a faulty gauge is not mistaken for a
product error.

Distances are great-circle km between `(lon, lat)` rows, as in `HeavyRainEvents`.
"""
module HeavyCoreDiagnostics

using Statistics
using Main.TraditionalInterpolation: haversine_distance_matrix

export upscaled_gauge_truth, neighbourhood_max, core_hourly_diagnostics, deficit_attribution
export neighbour_consistency

"""Distances in km from `center` `(lon, lat)` to every row of `lonlat`."""
_distances_km(lonlat::AbstractMatrix{<:Real}, center) =
    vec(haversine_distance_matrix(lonlat, [center[1] center[2]]))

"""
Areal gauge truth around gauge `i`: the mean of the finite `totals` of the gauges within each
radius (gauge `i` included), one `NamedTuple` per radius with `radius_km`, `n` and `mean`.

A product that reproduced the areal mean exactly would still miss `totals[i] - mean`; that gap is
the representativeness part of the deficit, not a product error.
"""
function upscaled_gauge_truth(
    lonlat::AbstractMatrix{<:Real}, totals::AbstractVector{<:Real}, i::Integer;
    radii_km=(5.0, 10.0, 15.0),
)
    size(lonlat, 1) == length(totals) || throw(DimensionMismatch("lonlat and totals must describe the same gauges"))
    distances = _distances_km(lonlat, (lonlat[i, 1], lonlat[i, 2]))
    return map(radii_km) do radius
        inside = [k for k in eachindex(totals) if distances[k] <= radius && isfinite(totals[k])]
        (; radius_km=Float64(radius), n=length(inside), mean=isempty(inside) ? NaN : mean(totals[inside]))
    end
end

"""
The largest finite `values` within each radius of `center` `(lon, lat)`, from points `lonlat`
(gauge pixels, or every cell of a flattened grid). One `NamedTuple` per radius: `radius_km`, `n`,
`max`, and where that maximum sits - `max_lon`, `max_lat`, `offset_km` from `center`.
"""
function neighbourhood_max(
    lonlat::AbstractMatrix{<:Real}, values::AbstractVector{<:Real}, center;
    radii_km=(10.0, 25.0, 50.0),
)
    size(lonlat, 1) == length(values) || throw(DimensionMismatch("lonlat and values must describe the same points"))
    distances = _distances_km(lonlat, center)
    return map(radii_km) do radius
        inside = [k for k in eachindex(values) if distances[k] <= radius && isfinite(values[k])]
        isempty(inside) && return (; radius_km=Float64(radius), n=0, max=NaN, max_lon=NaN, max_lat=NaN, offset_km=NaN)
        best = inside[argmax(values[inside])]
        (; radius_km=Float64(radius), n=length(inside), max=Float64(values[best]),
            max_lon=Float64(lonlat[best, 1]), max_lat=Float64(lonlat[best, 2]), offset_km=distances[best])
    end
end

"""
Hourly behaviour of one core gauge against a product over an event.

`obs` and `sat` are hourly series over the event window padded by `pad` hours on each side, so
`obs[pad+1:end-pad]` is the event itself. Returns a `NamedTuple`:

- `obs_total`, `sat_total` over the event, and `sat_total_padded` over the padded window: rain the
  product placed just outside the event window, which a timing error would move there
- `obs_peak`, `sat_peak` (largest hourly values), `peak_ratio = sat_peak / obs_peak` and
  `peak_lag_h` (product peak hour minus gauge peak hour; positive = product late)
- `obs_wet_hours`, `sat_wet_hours`: hours with at least `wet_mm`
- `top_share_obs`: share of the gauge total that falls in its `top_k` wettest hours
- `top_deficit_share`: share of the event deficit that falls in those same hours. Near 1 means the
  product misses the bursts and gets the rest right; near `top_share_obs` means it is uniformly low.
- `best_lag_h`, `best_lag_r`: the lag in `-max_lag:max_lag` that maximises the correlation of
  `obs[t]` with `sat[t + lag]`, and `r_lag0` for comparison
"""
function core_hourly_diagnostics(
    obs::AbstractVector{<:Real}, sat::AbstractVector{<:Real};
    pad::Integer=3, max_lag::Integer=pad, top_k::Integer=3, wet_mm::Real=0.1,
)
    length(obs) == length(sat) || throw(DimensionMismatch("obs and sat must cover the same hours"))
    0 <= max_lag <= pad || throw(ArgumentError("max_lag must lie in 0:pad"))
    event = (pad + 1):(length(obs) - pad)
    o, s = Float64.(obs[event]), Float64.(sat[event])
    all(isfinite, o) && all(isfinite, s) || throw(ArgumentError("the event window must be complete"))

    obs_total, sat_total = sum(o), sum(s)
    padded = filter(isfinite, Float64.(sat))
    top = partialsortperm(o, 1:min(top_k, length(o)); rev=true)
    deficit = obs_total - sat_total
    lagged_r(lag) = begin
        shifted = Float64.(sat[event .+ lag])
        all(isfinite, shifted) && std(o) > 0 && std(shifted) > 0 ? cor(o, shifted) : NaN
    end
    lags = collect(-max_lag:max_lag)
    rs = lagged_r.(lags)
    finite = findall(isfinite, rs)
    best = isempty(finite) ? 0 : finite[argmax(rs[finite])]
    return (;
        obs_total, sat_total, sat_total_padded=sum(padded),
        obs_peak=maximum(o), sat_peak=maximum(s), peak_ratio=maximum(s) / maximum(o),
        peak_lag_h=argmax(s) - argmax(o),
        obs_wet_hours=count(>=(wet_mm), o), sat_wet_hours=count(>=(wet_mm), s),
        top_share_obs=obs_total > 0 ? sum(o[top]) / obs_total : NaN,
        top_deficit_share=deficit > 0 ? sum(o[top] .- s[top]) / deficit : NaN,
        best_lag_h=best == 0 ? 0 : lags[best], best_lag_r=best == 0 ? NaN : rs[best],
        r_lag0=rs[max_lag + 1],
    )
end

"""
Split the deficit `obs - sat` at a core gauge into representativeness, displacement and magnitude.

- `representativeness = obs - areal`, where `areal` is the upscaled gauge truth: what a perfect
  areal product would still miss.
- `displacement = clamp(nearby - sat, 0, areal - sat)`, where `nearby` is the product's maximum
  near the gauge: the part of the areal amount the product has, but elsewhere. Capped at
  `areal - sat` so that a wet pixel further out cannot explain more than the pixel-scale deficit.
- `magnitude = areal - sat - displacement`: the areal amount the product has nowhere nearby.

The three sum to the deficit. When the product exceeds the areal truth (`sat > areal`) the
displacement and magnitude parts are zero and negative, and `representativeness` exceeds the
deficit. Shares are each part over the deficit, `NaN` when there is no deficit.
"""
function deficit_attribution(obs::Real, sat::Real, areal::Real, nearby::Real)
    deficit = obs - sat
    representativeness = obs - areal
    displacement = clamp(nearby - sat, 0.0, max(areal - sat, 0.0))
    magnitude = areal - sat - displacement
    share(x) = deficit > 0 ? x / deficit : NaN
    return (;
        deficit, representativeness, displacement, magnitude,
        representativeness_share=share(representativeness), displacement_share=share(displacement),
        magnitude_share=share(magnitude),
    )
end

"""
Gauge `i` against its `k` nearest gauges with a finite total: their mean distance, the mean and
maximum of their totals, and `ratio = totals[i] / neighbour mean`. A large ratio over short
distances, with no neighbour close to gauge `i`, is the signature of a faulty gauge rather than a
real core.
"""
function neighbour_consistency(lonlat::AbstractMatrix{<:Real}, totals::AbstractVector{<:Real}, i::Integer; k::Integer=3)
    distances = _distances_km(lonlat, (lonlat[i, 1], lonlat[i, 2]))
    candidates = [j for j in eachindex(totals) if j != i && isfinite(totals[j])]
    nearest = candidates[partialsortperm(distances[candidates], 1:min(k, length(candidates)))]
    neighbour_mean = mean(totals[nearest])
    return (;
        k=length(nearest), mean_distance_km=mean(distances[nearest]),
        neighbour_mean, neighbour_max=maximum(totals[nearest]), ratio=totals[i] / neighbour_mean,
    )
end

end # module
