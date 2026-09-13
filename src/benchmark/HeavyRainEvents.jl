"""
Heavy-rainfall event evaluation: does a satellite product reproduce the spatial pattern of
gauge-observed rain on a heavy-rain day?

The interpolation benchmark scores every station-hour, where ~93% of cells are dry. This module
looks only at representative heavy-rain days and compares, across the gauge network, each
product's daily accumulation with the gauges' - spatial correlation, error, and whether the rain
centre and the heavy-rain footprint land in the right place.

Conventions every function here relies on:

- Time labels are hour-*ending* Beijing time, the convention every study-area table shares:
  gauges natively, GPM/GSMaP via `MGERDataPrep.prepare_satellite_wide` (GEE UTC hour-start + 9 h),
  FY4B via `FY4BPreprocessing.aligned_hour_end`. The label `D 09:00` is the first hour of the
  08-08 BJT day `D`; see `met_day`.
- Matrices are `[station, time]`, as `read_hourly_wide` returns them, with `NaN` for missing.
- Every product is scored over the *same* hours and the *same* stations. FY4B systematically
  lacks the 23:00 label and often 00:00-01:00, so a full-day FY4B total would be biased low next
  to a full-day GPM total - a gap, not a product error. `event_hours` therefore builds one common
  hour set per day and reports how much gauge rain it leaves out.
- Metrics come only from the paired gauge-pixel values. `idw_surface` exists for visualization
  and is never scored.
"""
module HeavyRainEvents

using DataFrames, Dates, Statistics
using MixedGWR: metric_continuous, metric_event
using Main.TraditionalInterpolation: haversine_distance_matrix, idw_predict

export met_day, align_to_reference, daily_totals, screen_heavy_rain_days, tag_rain_processes!
export event_hours, select_representative_events, event_station_totals, spatial_metrics
export idw_surface

"""
The meteorological day an hour-ending label belongs to.

With `day_start_hour=8` (the 08-08 BJT day) the labels `D 09:00` through `D+1 08:00` all map to
`D`. `day_start_hour=20` gives the 20-20 BJT day instead.
"""
met_day(time::DateTime; day_start_hour::Int=8) = Date(time - Hour(day_start_hour + 1))

"""
Re-index a `[station, time]` matrix onto a reference station order and time axis.

A reference hour absent from `times` becomes a `NaN` column, which is how FY4B's missing hours
enter the event window. A reference *station* absent from `ids` is an error rather than a `NaN`
row: a product silently missing a gauge would shrink every metric's sample without a trace.
"""
function align_to_reference(
    ref_times::AbstractVector{DateTime}, ref_ids::AbstractVector{<:AbstractString},
    times::AbstractVector{DateTime}, ids::AbstractVector{<:AbstractString},
    Y::AbstractMatrix{<:Real},
)
    size(Y) == (length(ids), length(times)) ||
        throw(DimensionMismatch("Y must be [station, time] matching ids and times"))
    allunique(ids) || throw(ArgumentError("ids contain duplicates"))
    allunique(times) || throw(ArgumentError("times contain duplicates"))
    row_by_id = Dict(id => row for (row, id) in enumerate(ids))
    missing_ids = filter(id -> !haskey(row_by_id, id), ref_ids)
    isempty(missing_ids) ||
        throw(ArgumentError("$(length(missing_ids)) reference stations are absent, e.g. $(first(missing_ids))"))
    column_by_time = Dict(time => column for (column, time) in enumerate(times))
    rows = [row_by_id[id] for id in ref_ids]
    aligned = fill(NaN, length(ref_ids), length(ref_times))
    for (j, time) in enumerate(ref_times)
        column = get(column_by_time, time, 0)
        column == 0 && continue
        for (i, row) in enumerate(rows)
            aligned[i, j] = Y[row, column]
        end
    end
    return aligned
end

"""
Daily accumulations `(days, totals[station, day])`.

A station-day total is `NaN` unless all 24 of its hours are finite, so a heavy-rain count never
rests on a partial day. Days with fewer than 24 time rows (the ends of the record) are skipped.
"""
function daily_totals(times::AbstractVector{DateTime}, Y::AbstractMatrix{<:Real}; day_start_hour::Int=8)
    size(Y, 2) == length(times) || throw(DimensionMismatch("Y columns must match times"))
    allunique(times) || throw(ArgumentError("times contain duplicates"))
    days = met_day.(times; day_start_hour)
    columns = Dict{Date,Vector{Int}}()
    for (j, day) in enumerate(days)
        push!(get!(columns, day, Int[]), j)
    end
    kept = sort([day for (day, idx) in columns if length(idx) == 24])
    totals = fill(NaN, size(Y, 1), length(kept))
    for (k, day) in enumerate(kept), i in axes(Y, 1)
        total = 0.0
        complete = true
        for j in columns[day]
            value = Y[i, j]
            isfinite(value) || (complete = false; break)
            total += value
        end
        complete && (totals[i, k] = total)
    end
    return kept, totals
end

"""
Days on which at least `min_stations` gauges exceed `threshold` mm - by default the heavy-rain
definition "more than three stations above 50 mm/day", with the comparison strict.

`n_wide`/`frac_wide` count gauges above `wide_threshold`, the areal extent measure
`select_representative_events` separates widespread from localized days with.
"""
function screen_heavy_rain_days(
    days::AbstractVector{Date}, obs_daily::AbstractMatrix{<:Real},
    ids::AbstractVector{<:AbstractString};
    threshold::Real=50.0, min_stations::Int=4, wide_threshold::Real=25.0,
)
    size(obs_daily) == (length(ids), length(days)) ||
        throw(DimensionMismatch("obs_daily must be [station, day] matching ids and days"))
    table = DataFrame(
        day=Date[], n_valid=Int[], n_heavy=Int[], n_wide=Int[], frac_wide=Float64[],
        max_mm=Float64[], mean_mm=Float64[], sd_mm=Float64[], cv=Float64[],
        peak_station_id=String[],
    )
    for (k, day) in enumerate(days)
        valid = findall(isfinite, @view obs_daily[:, k])
        isempty(valid) && continue
        values = obs_daily[valid, k]
        n_heavy = count(>(threshold), values)
        n_heavy >= min_stations || continue
        n_wide = count(>(wide_threshold), values)
        push!(table, (
            day, length(valid), n_heavy, n_wide, n_wide / length(valid),
            maximum(values), mean(values), std(values), std(values) / mean(values),
            String(ids[valid[argmax(values)]]),
        ))
    end
    table.threshold_mm = fill(Float64(threshold), nrow(table))
    table.wide_threshold_mm = fill(Float64(wide_threshold), nrow(table))
    return table
end

"""
Number runs of consecutive qualifying days as one rain process (`process_id`), so a multi-day
system is not selected twice under two dates.
"""
function tag_rain_processes!(table::DataFrame)
    issorted(table.day) || throw(ArgumentError("table must be sorted by day"))
    process = zeros(Int, nrow(table))
    for i in eachindex(process)
        process[i] = i == 1 ? 1 :
            (table.day[i] - table.day[i - 1] == Day(1) ? process[i - 1] : process[i - 1] + 1)
    end
    table.process_id = process
    return table
end

"""
The hours of `day` every source can be scored over.

An hour is available to a source when at least `min_station_fraction` of its stations are finite
there. `hour_idx` is the set available to the gauges *and* every product at once - the common
window all sources are accumulated over. Two shares say whether that window still describes the
event:

- `excluded_rain_share`: the fraction of the day's gauge rain, over the whole network, that falls
  outside the window.
- `peak_excluded_rain_share`: the same fraction at the storm maximum alone - the gauge with the
  largest complete 24-h total. A network-wide share can look harmless while the storm core falls
  in the missing hours: on 2024-07-30 FY4B's two gaps removed 17% of the network's rain but 48% of
  the peak gauge's 144 mm. Pooling over every gauge above 50 mm dilutes it the same way (19%),
  and a percentile over those gauges is too strict for widespread days, where some gauge on the
  edge of a rain band always loses a large share. `NaN` when no gauge has rain.
"""
function event_hours(
    times::AbstractVector{DateTime}, Y_obs::AbstractMatrix{<:Real},
    products::AbstractDict{<:AbstractString,<:AbstractMatrix}, day::Date;
    day_start_hour::Int=8, min_station_fraction::Real=0.9,
)
    size(Y_obs, 2) == length(times) || throw(DimensionMismatch("Y_obs columns must match times"))
    all(Y -> size(Y) == size(Y_obs), values(products)) ||
        throw(DimensionMismatch("every product must be aligned to Y_obs"))
    day_hour_idx = findall(time -> met_day(time; day_start_hour) == day, times)
    length(day_hour_idx) == 24 ||
        throw(ArgumentError("$day has $(length(day_hour_idx)) hours on the time axis, expected 24"))

    n_station = size(Y_obs, 1)
    available(Y, j) = count(isfinite, @view Y[:, j]) >= min_station_fraction * n_station
    product_hours = Dict(String(p) => count(j -> available(Y, j), day_hour_idx) for (p, Y) in products)
    hour_idx = filter(j -> available(Y_obs, j) && all(Y -> available(Y, j), values(products)), day_hour_idx)

    excluded_hour_idx = setdiff(day_hour_idx, hour_idx)
    hour_rain(j) = sum(value for value in @view(Y_obs[:, j]) if isfinite(value); init=0.0)
    day_rain = sum(hour_rain, day_hour_idx)
    excluded_rain = sum(hour_rain, excluded_hour_idx; init=0.0)
    excluded_rain_share = day_rain > 0 ? excluded_rain / day_rain : NaN

    # The same gauge `screen_heavy_rain_days` reports as `peak_station_id`: only complete days count.
    totals = [all(isfinite, @view Y_obs[i, day_hour_idx]) ? sum(@view Y_obs[i, day_hour_idx]) : -Inf
        for i in 1:n_station]
    peak = argmax(totals)
    peak_excluded_rain_share = totals[peak] > 0 ?
        sum((Y_obs[peak, j] for j in excluded_hour_idx); init=0.0) / totals[peak] : NaN
    return (;
        day_hour_idx, hour_idx, common_hours=length(hour_idx), product_hours,
        excluded_rain_share, peak_excluded_rain_share,
    )
end

"""
Pick the representative events: `n_widespread` widespread days, then `n_localized` localized ones.

- *Eligible*: the common window keeps at least `min_common_hours` hours and leaves out at most
  `max_excluded_rain_share` of the day's gauge rain - both over the whole network and at the peak
  gauge (`peak_excluded_rain_share`, see `event_hours`). A product with no data that day
  has no common hours, so this also requires every product to be present.
- *Widespread* (`frac_wide >= localized_frac`): ranked by gauges above the heavy-rain threshold.
- *Localized* (`frac_wide < localized_frac`): ranked by the single-gauge maximum - intense rain
  over a small footprint, the hardest case for a satellite's resolution.
- One day per `process_id` across both groups; ties go to the earlier date.

Returns a copy of `table` with `eligible`, `selected` and `event_type` columns.
"""
function select_representative_events(
    table::DataFrame; n_widespread::Int=4, n_localized::Int=2, localized_frac::Real=0.25,
    min_common_hours::Int=20, max_excluded_rain_share::Real=0.20,
)
    required = [
        :day, :n_heavy, :frac_wide, :max_mm, :process_id, :common_hours,
        :excluded_rain_share, :peak_excluded_rain_share,
    ]
    absent = setdiff(required, propertynames(table))
    isempty(absent) || throw(ArgumentError("table is missing columns: $(join(absent, ", "))"))

    selection = copy(table)
    selection.eligible = (selection.common_hours .>= min_common_hours) .&
        (selection.excluded_rain_share .<= max_excluded_rain_share) .&
        (selection.peak_excluded_rain_share .<= max_excluded_rain_share)
    selection.selected = falses(nrow(selection))
    selection.event_type = fill("", nrow(selection))
    used_processes = Set{Int}()

    function pick!(candidates, n, event_type)
        taken = 0
        for i in candidates
            taken == n && break
            selection.process_id[i] in used_processes && continue
            selection.selected[i] = true
            selection.event_type[i] = event_type
            push!(used_processes, selection.process_id[i])
            taken += 1
        end
        taken < n && @warn "only $taken of $n $event_type events are eligible"
        return nothing
    end

    eligible = findall(selection.eligible)
    widespread = filter(i -> selection.frac_wide[i] >= localized_frac, eligible)
    pick!(sort(widespread; by=i -> (-selection.n_heavy[i], selection.day[i])), n_widespread, "widespread")
    localized = filter(i -> selection.frac_wide[i] < localized_frac, eligible)
    pick!(sort(localized; by=i -> (-selection.max_mm[i], selection.day[i])), n_localized, "localized")
    return selection
end

"""
Per-station accumulations over `hour_idx` for the gauges and every product.

A station is kept only when the gauge and *every* product are finite over all of `hour_idx`, so
all products are scored on one identical station set; dropped stations are `NaN` in every vector.
"""
function event_station_totals(
    Y_obs::AbstractMatrix{<:Real}, products::AbstractDict{<:AbstractString,<:AbstractMatrix},
    hour_idx::AbstractVector{<:Integer},
)
    all(Y -> size(Y) == size(Y_obs), values(products)) ||
        throw(DimensionMismatch("every product must be aligned to Y_obs"))
    complete(Y) = vec(all(isfinite, Y[:, hour_idx]; dims=2))
    keep = complete(Y_obs)
    for Y in values(products)
        keep .&= complete(Y)
    end
    total(Y) = [keep[i] ? sum(@view Y[i, hour_idx]) : NaN for i in axes(Y, 1)]
    return (; keep, obs=total(Y_obs), sat=Dict(String(p) => total(Y) for (p, Y) in products))
end

"""Ranks with ties given their average rank, for Spearman correlation."""
function _tied_ranks(x::AbstractVector{<:Real})
    order = sortperm(x)
    ranks = Vector{Float64}(undef, length(x))
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && x[order[j + 1]] == x[order[i]]
            j += 1
        end
        ranks[order[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return ranks
end

"""Great-circle distance in km between two `(lon, lat)` points."""
_distance_km(a, b) = haversine_distance_matrix([a[1] a[2]], [b[1] b[2]])[1, 1]

"""Precipitation-weighted centre of the network, `(lon, lat)`; `NaN`s when nothing fell."""
function _rain_centroid(lonlat::AbstractMatrix{<:Real}, values::AbstractVector{<:Real})
    weights = max.(values, 0.0)
    total = sum(weights)
    total > 0 || return (NaN, NaN)
    return (sum(weights .* lonlat[:, 1]) / total, sum(weights .* lonlat[:, 2]) / total)
end

"""Metric-name suffix for a threshold: `25.0` -> `"25"`, `2.5` -> `"2p5"`."""
_threshold_label(threshold::Real) =
    isinteger(threshold) ? string(Int(threshold)) : replace(string(threshold), "." => "p")

"""
Spatial agreement between gauge accumulations `obs` and a product's `est` at the same stations.

Pairs with either side non-finite are dropped. Returns a `NamedTuple`:

- `n`, `obs_mean`, `est_mean`
- `RMSE`, `MAE`, `Bias`, `r` from `metric_continuous`; `RB_pct = 100 Bias / obs_mean`
- `rho`: Spearman correlation - robust to the skew of heavy-rain totals
- `CRMSE`: centred RMSE, the pattern error left once the mean bias is removed
- `sd_ratio` (sigma_est / sigma_obs, population std, so `CRMSE`, `sd_ratio` and `r` satisfy the
  Taylor-diagram identity) and `cv_ratio`: whether the product reproduces the *amplitude* of the
  spatial heterogeneity, not only its arrangement
- `KGE = 1 - sqrt((r-1)^2 + (sd_ratio-1)^2 + (est_mean/obs_mean-1)^2)`
- `POD_<t>`, `FAR_<t>`, `CSI_<t>` per threshold, from `metric_event`, which counts `>= t` - unlike
  the strict `>` of `screen_heavy_rain_days`
- `centroid_shift_km`: distance between the gauges' and the product's rain-weighted centres
- `peak_shift_km`: distance between the gauge with the largest observed total and the gauge where
  the product puts its maximum
"""
function spatial_metrics(
    obs::AbstractVector{<:Real}, est::AbstractVector{<:Real}, lonlat::AbstractMatrix{<:Real};
    thresholds=(25.0, 50.0),
)
    length(obs) == length(est) == size(lonlat, 1) ||
        throw(DimensionMismatch("obs, est and lonlat must describe the same stations"))
    pairs = findall(i -> isfinite(obs[i]) && isfinite(est[i]), eachindex(obs))
    length(pairs) >= 2 || throw(ArgumentError("need at least two finite station pairs"))
    o = Float64.(obs[pairs])
    e = Float64.(est[pairs])
    xy = Float64.(lonlat[pairs, :])

    continuous = metric_continuous(o, e)
    obs_mean, est_mean = mean(o), mean(e)
    sd_obs, sd_est = std(o; corrected=false), std(e; corrected=false)
    sd_ratio = sd_est / sd_obs
    cv_ratio = (sd_est / est_mean) / (sd_obs / obs_mean)
    crmse = sqrt(mean(((e .- est_mean) .- (o .- obs_mean)) .^ 2))
    kge = 1 - sqrt((continuous.r - 1)^2 + (sd_ratio - 1)^2 + (est_mean / obs_mean - 1)^2)

    event_names = Symbol[]
    event_values = Float64[]
    for threshold in thresholds
        scores = metric_event(o, e; thr=Float64(threshold))
        label = _threshold_label(threshold)
        append!(event_names, Symbol.(["POD_", "FAR_", "CSI_"], label))
        append!(event_values, [scores.POD, scores.FAR, scores.CSI])
    end

    obs_centre = _rain_centroid(xy, o)
    est_centre = _rain_centroid(xy, e)
    centroid_shift_km = all(isfinite, (obs_centre..., est_centre...)) ?
        _distance_km(obs_centre, est_centre) : NaN
    peak_shift_km = _distance_km(xy[argmax(o), :], xy[argmax(e), :])

    head = (;
        n=length(pairs), obs_mean, est_mean,
        continuous.RMSE, continuous.MAE, continuous.Bias, RB_pct=100 * continuous.Bias / obs_mean,
        continuous.r, rho=cor(_tied_ranks(o), _tied_ranks(e)), CRMSE=crmse,
        sd_ratio, cv_ratio, KGE=kge,
    )
    return merge(head, NamedTuple{Tuple(event_names)}(Tuple(event_values)), (; centroid_shift_km, peak_shift_km))
end

"""
Inverse-distance surface of station `values` on a regular lon/lat grid, as a long
`(lon, lat, value)` table.

For visualization only: it draws a product's values *sampled at the gauges*, not the product's
native field, and no metric is computed from it. Every gauge with a finite value contributes, so
calling it with the same settings for the gauges and each product gives surfaces that differ only
in their input values. `bounds` is `(west, east, south, north)`.
"""
function idw_surface(
    lonlat::AbstractMatrix{<:Real}, values::AbstractVector{<:Real};
    bounds::NTuple{4,Real}, step_deg::Real=0.025, power::Real=2.0,
)
    size(lonlat, 1) == length(values) ||
        throw(DimensionMismatch("lonlat rows must match values"))
    west, east, south, north = Float64.(bounds)
    step_deg > 0 || throw(ArgumentError("step_deg must be positive"))
    lons = range(west, east; length=round(Int, (east - west) / step_deg) + 1)
    lats = range(south, north; length=round(Int, (north - south) / step_deg) + 1)
    grid = [repeat(collect(lons), inner=length(lats)) repeat(collect(lats), outer=length(lons))]
    valid = findall(isfinite, values)
    isempty(valid) && throw(ArgumentError("no finite station values"))
    surface = idw_predict(lonlat[valid, :], reshape(Float64.(values[valid]), :, 1), grid; power)
    return DataFrame(lon=grid[:, 1], lat=grid[:, 2], value=surface[:, 1])
end

end # module
