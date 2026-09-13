"""
Temporal evaluation of satellite precipitation against gauges: does a product capture how rain evolves
hour by hour, and does that change with rain intensity, season and landform region?

`HeavyRainEvents` scores the spatial pattern of daily totals on selected heavy-rain days. This module
scores the time dimension over the whole record, in six tables:

- `intensity_class_table`: gauge-wet hours stratified by the gauge's hourly intensity class, plus the
  gauge-class x satellite-class confusion counts.
- `false_alarm_table`: gauge-dry hours on their own - how often, and how hard, the satellite rains
  when the gauge does not.
- `event_tables` / `event_summary`: station rain events and the timing, detection and volume errors
  of each product over each event.
- `diurnal_cycle_table` / `diurnal_summary`: composite diurnal cycles and their phase errors.
- `station_correlation_table` / `station_correlation_summary`: hourly correlation at each gauge.
- `regional_series_table`: correlation and error of regional-mean series at 1, 3 and 24 h.

Conventions every function relies on:

- Matrices are `[station, time]` on one contiguous hour-ending Beijing-time axis, `NaN` for missing,
  every product already aligned to the gauges (`HeavyRainEvents.align_to_reference`).
- A `mask` names the cells a *sample* scores: those where the gauge and every product in the sample are
  finite, so every product in a sample is scored on identical cells.
- Hourly classes (`INTENSITY_BOUNDS`) are lower-bound inclusive: no rain < 0.1, light [0.1, 2),
  moderate [2, 4), heavy [4, 8), rainstorm [8, 20), severe rainstorm >= 20 mm/h.
- The gauges are 0.5 mm tipping buckets: their "light" hours are 0.5, 1.0 or 1.5 mm and they cannot
  see 0.1-0.5 mm. Where a satellite threshold matters, a `_res` variant repeats the calculation at the
  gauge resolution `GAUGE_RESOLUTION_MM`.
- Seasons are meteorological (MAM, JJA, SON, DJF) by the month of the 08-08 BJT day (`met_day`);
  day-block bootstraps resample those days. Every stratum is also reported as `all`.
"""
module SatelliteTemporalEvaluation

using DataFrames, Dates, Random, Statistics
using MixedGWR: metric_continuous, metric_event
using Main.HeavyRainEvents: met_day, _tied_ranks

export WET_MM, GAUGE_RESOLUTION_MM, INTENSITY_BOUNDS, INTENSITY_CLASSES, SEASONS, FALSE_ALARM_BINS
export intensity_class, season_of, evaluation_context, sample_mask, day_block_bootstrap
export intensity_class_table, false_alarm_table, sample_coverage_table
export station_rain_events, score_event, best_lag, event_tables, event_summary
export harmonic_phase, circular_hour_difference, diurnal_cycle_table, diurnal_summary
export station_correlation_table, station_correlation_summary
export regional_mean_series, aggregate_series, series_metrics, regional_series_table

const WET_MM = 0.1
const GAUGE_RESOLUTION_MM = 0.5
const INTENSITY_BOUNDS = [0.1, 2.0, 4.0, 8.0, 20.0]
const INTENSITY_CLASSES = ["no_rain", "light", "moderate", "heavy", "rainstorm", "severe_rainstorm"]
const SEASONS = ["MAM", "JJA", "SON", "DJF"]
const SEASON_LABELS = ["all"; SEASONS]
const FALSE_ALARM_BINS = ["no_rain", "below_gauge_resolution", "light_resolved", "moderate", "heavy",
    "rainstorm", "severe_rainstorm"]

"""Hourly intensity class index: 0 no rain, 1 light, ..., 5 severe rainstorm; -1 for a non-finite value."""
intensity_class(value::Real) = isfinite(value) ? searchsortedlast(INTENSITY_BOUNDS, value) : -1

"""Meteorological season of a date."""
function season_of(day::Date)
    m = month(day)
    return m in 3:5 ? "MAM" : m in 6:8 ? "JJA" : m in 9:11 ? "SON" : "DJF"
end

"""
The strata every table shares, computed once: each column's met-day, season, hour of day and hour
within the met-day, and each station's region index.

`regions` fixes the region order (index 1, 2, ...); index 0 means "all" everywhere below, as does
season index 0. The axis must be contiguous hourly, which events, lags and 3-h blocks depend on.
"""
function evaluation_context(
    times::AbstractVector{DateTime}, station_region::AbstractVector{<:AbstractString};
    regions::AbstractVector{<:AbstractString}=sort(unique(station_region)), day_start_hour::Int=8,
)
    length(times) >= 2 && all(==(Hour(1)), diff(times)) ||
        throw(ArgumentError("times must be a contiguous hourly axis"))
    unknown = setdiff(unique(station_region), regions)
    isempty(unknown) || throw(ArgumentError("stations carry regions not listed in `regions`: $(join(unknown, ", "))"))
    days = met_day.(times; day_start_hour)
    unique_days = unique(days)
    day_number = Dict(day => k for (k, day) in enumerate(unique_days))
    column_day = [day_number[day] for day in days]
    day_season = [findfirst(==(season_of(day)), SEASONS) for day in unique_days]
    column_hour = hour.(times)
    return (;
        times=collect(times), days=unique_days, column_day, day_season,
        column_season=day_season[column_day], column_hour,
        column_offset=mod.(column_hour .- (day_start_hour + 1), 24),
        regions=String.(regions), station_region_index=[findfirst(==(r), regions) for r in station_region],
    )
end

_region_labels(ctx) = ["all"; ctx.regions]

"""Cells where every matrix is finite."""
function sample_mask(matrices::AbstractMatrix...)
    mask = trues(size(first(matrices)))
    for Y in matrices
        size(Y) == size(mask) || throw(DimensionMismatch("every matrix must have the same size"))
        mask .&= isfinite.(Y)
    end
    return mask
end

function _check_inputs(ctx, matrices...)
    for Y in matrices
        size(Y, 2) == length(ctx.times) || throw(DimensionMismatch("matrix columns must match the time axis"))
        size(Y, 1) == length(ctx.station_region_index) || throw(DimensionMismatch("matrix rows must match the stations"))
    end
end

"""Pearson r over the pairs where both are finite; `NaN` below `min_n` pairs or with a constant side."""
function _finite_cor(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; min_n::Int=2)
    keep = [i for i in eachindex(x, y) if isfinite(x[i]) && isfinite(y[i])]
    length(keep) >= max(min_n, 2) || return NaN
    a, b = Float64.(x[keep]), Float64.(y[keep])
    (std(a) > 0 && std(b) > 0) || return NaN
    return cor(a, b)
end

_spearman(x, y; min_n::Int=2) = length(x) >= max(min_n, 2) ? _finite_cor(_tied_ranks(x), _tied_ranks(y); min_n) : NaN
_median(v) = isempty(v) ? NaN : median(v)
_quantile(v, p) = isempty(v) ? NaN : quantile(v, p)
_share(v) = isempty(v) ? NaN : mean(v)
_row(pairs) = NamedTuple{Tuple(first.(pairs))}(Tuple(last.(pairs)))

"""
Percentile confidence interval of `statistic` under a day-block bootstrap.

`sums[cell, k]` are additive per-cell quantities; they are summed per `day` (positive integers), days
are resampled with replacement `reps` times, and `statistic(totals)` maps the summed vector to a tuple
of statistics. Resampling whole days keeps the hours and stations of one weather system together.
Returns `(; lo, hi)`, or `nothing` when there are no cells.
"""
function day_block_bootstrap(
    day::AbstractVector{<:Integer}, sums::AbstractMatrix{<:Real}, statistic::Function;
    reps::Int, rng::AbstractRNG, level::Real=0.95,
)
    length(day) == size(sums, 1) || throw(DimensionMismatch("one day per row of sums"))
    isempty(day) && return nothing
    per_day = zeros(maximum(day), size(sums, 2))
    for i in eachindex(day)
        @views per_day[day[i], :] .+= sums[i, :]
    end
    per_day = per_day[sort(unique(day)), :]
    n_days = size(per_day, 1)
    totals = zeros(size(sums, 2))
    draws = map(1:reps) do _
        fill!(totals, 0.0)
        for d in rand(rng, 1:n_days, n_days)
            @views totals .+= per_day[d, :]
        end
        collect(Float64, statistic(totals))
    end
    alpha = (1 - level) / 2
    columns = [filter(isfinite, getindex.(draws, j)) for j in eachindex(first(draws))]
    return (; lo=[_quantile(c, alpha) for c in columns], hi=[_quantile(c, 1 - alpha) for c in columns])
end

"""Cells scored per sample, season and region: the denominators behind every other table."""
function sample_coverage_table(ctx, Y_obs::AbstractMatrix, mask::AbstractMatrix; sample::AbstractString)
    _check_inputs(ctx, Y_obs, mask)
    G = length(ctx.regions) + 1
    n_cells, n_wet = zeros(Int, 5, G), zeros(Int, 5, G)
    hours, days = [Set{Int}() for _ in 1:5, _ in 1:G], [Set{Int}() for _ in 1:5, _ in 1:G]
    for j in axes(Y_obs, 2), i in axes(Y_obs, 1)
        mask[i, j] || continue
        s, g = ctx.column_season[j] + 1, ctx.station_region_index[i] + 1
        for (a, b) in ((1, 1), (1, g), (s, 1), (s, g))
            n_cells[a, b] += 1
            n_wet[a, b] += Y_obs[i, j] >= WET_MM
            push!(hours[a, b], j)
            push!(days[a, b], ctx.column_day[j])
        end
    end
    return DataFrame([(; sample, season=SEASON_LABELS[a], region=_region_labels(ctx)[b], n_hours=length(hours[a, b]),
        n_days=length(days[a, b]), n_cells=n_cells[a, b], n_gauge_wet=n_wet[a, b],
        gauge_wet_share=n_cells[a, b] > 0 ? n_wet[a, b] / n_cells[a, b] : NaN) for a in 1:5 for b in 1:G])
end

"""Metrics of one stratum of gauge-wet cells; see `intensity_class_table`."""
function _wet_metrics(obs, sat, obs_class, sat_class, day, station; with_r::Bool, reps::Int, rng::AbstractRNG)
    n = length(obs)
    err = sat .- obs
    base = (;
        n, n_stations=length(unique(station)), n_days=length(unique(day)),
        obs_mean=mean(obs), sat_mean=mean(sat), Bias=mean(err),
        RB_pct=100 * (sum(sat) - sum(obs)) / sum(obs), MAE=mean(abs.(err)), RMSE=sqrt(mean(abs2.(err))),
        POD_rain=_share(sat .>= WET_MM), POD_rain_res=_share(sat .>= GAUGE_RESOLUTION_MM),
        class_hit=_share(sat_class .== obs_class), under_class=_share(sat_class .< obs_class),
        over_class=_share(sat_class .> obs_class), median_ratio=_median(sat ./ obs),
        r=with_r ? _finite_cor(obs, sat) : NaN,
        low_sample=n < 100 || length(unique(day)) < 10,
    )
    sums = [ones(n) obs sat abs2.(err) (sat .>= WET_MM) (sat_class .== obs_class)]
    ci = day_block_bootstrap(day, sums, t -> (
        (t[3] - t[2]) / t[1], 100 * (t[3] - t[2]) / t[2], sqrt(t[4] / t[1]), t[5] / t[1], t[6] / t[1],
    ); reps, rng)
    names = (:Bias, :RB_pct, :RMSE, :POD_rain, :class_hit)
    bounds = ci === nothing ? fill(NaN, 2 * length(names)) : vcat(([ci.lo[k], ci.hi[k]] for k in eachindex(names))...)
    labels = vcat(([Symbol(name, "_lo"), Symbol(name, "_hi")] for name in names)...)
    return merge(base, NamedTuple{Tuple(labels)}(Tuple(bounds)))
end

"""
Hourly metrics by gauge intensity class, over gauge-wet cells only (`gauge >= 0.1`), and the
gauge-class x satellite-class confusion counts over the same cells.

Per sample x product x season x region x gauge class (`all_wet` pools the five classes): n, stations,
days, means, Bias, RB_pct, MAE, RMSE, `POD_rain` (satellite >= 0.1; `_res` at 0.5), `class_hit` /
`under_class` / `over_class` (satellite class equal / lower / higher, "lower" including no rain), the
median satellite/gauge ratio, and r on `all_wet` only - within one class the gauge range is truncated,
so r would mostly measure that truncation. `_lo`/`_hi` are 95% day-block bootstrap bounds;
`low_sample` flags n < 100 or fewer than 10 days.
"""
function intensity_class_table(
    ctx, Y_obs::AbstractMatrix, Y_sat::AbstractMatrix, mask::AbstractMatrix;
    sample::AbstractString, product::AbstractString, reps::Int=1000, seed::Int=20260913,
)
    _check_inputs(ctx, Y_obs, Y_sat, mask)
    n_station = size(Y_obs, 1)
    cells = [k for k in eachindex(Y_obs) if mask[k] && Y_obs[k] >= WET_MM]
    station = [(k - 1) % n_station + 1 for k in cells]
    column = [(k - 1) ÷ n_station + 1 for k in cells]
    obs, sat = Float64[Y_obs[k] for k in cells], Float64[Y_sat[k] for k in cells]
    obs_class, sat_class = intensity_class.(obs), intensity_class.(sat)
    day, season = ctx.column_day[column], ctx.column_season[column]
    region = ctx.station_region_index[station]

    metrics, confusion = NamedTuple[], NamedTuple[]
    stratum = 0
    for (s, season_label) in enumerate(SEASON_LABELS), (g, region_label) in enumerate(_region_labels(ctx))
        base = findall(i -> (s == 1 || season[i] == s - 1) && (g == 1 || region[i] == g - 1), eachindex(obs))
        for c in 0:5
            selected = c == 0 ? base : filter(i -> obs_class[i] == c, base)
            stratum += 1
            scores = _wet_metrics(obs[selected], sat[selected], obs_class[selected], sat_class[selected],
                day[selected], station[selected]; with_r=c == 0, reps, rng=MersenneTwister(seed + stratum))
            label = c == 0 ? "all_wet" : INTENSITY_CLASSES[c + 1]
            push!(metrics, merge((; sample, product, season=season_label, region=region_label, gauge_class=label), scores))
        end
        counts = zeros(Int, 5, 6)
        for i in base
            counts[obs_class[i], sat_class[i] + 1] += 1
        end
        for c in 1:5, sc in 0:5
            row_total = sum(counts[c, :])
            push!(confusion, (; sample, product, season=season_label, region=region_label,
                gauge_class=INTENSITY_CLASSES[c + 1], sat_class=INTENSITY_CLASSES[sc + 1], n=counts[c, sc + 1],
                row_fraction=row_total > 0 ? counts[c, sc + 1] / row_total : NaN))
        end
    end
    return DataFrame(metrics), DataFrame(confusion)
end

"""Bin of a satellite value on a gauge-dry hour, splitting light rain at the gauge resolution."""
_false_alarm_bin(value::Real) = value < WET_MM ? 1 : value < GAUGE_RESOLUTION_MM ? 2 : intensity_class(value) + 2

"""
The side table for gauge-dry hours (`gauge < 0.1`), which the intensity classes exclude.

Per sample x product x season x region:
- `false_alarm_rate`: P(satellite >= 0.1 | gauge dry); `false_alarm_rate_res` at 0.5
- `dry_share_<bin>`: the satellite's value class on gauge-dry hours, with light rain split at the 0.5 mm
  gauge resolution (`below_gauge_resolution` hours cannot be called false alarms with confidence)
- `FAR` / `FAR_res`: false-alarm hours over all satellite-wet hours
- `dry_volume_share`: fraction of the satellite's rain volume that falls on gauge-dry hours
- `gauge_dry_given_<class>`: for each satellite rain class, the fraction of its hours the gauge is dry
"""
function false_alarm_table(
    ctx, Y_obs::AbstractMatrix, Y_sat::AbstractMatrix, mask::AbstractMatrix;
    sample::AbstractString, product::AbstractString,
)
    _check_inputs(ctx, Y_obs, Y_sat, mask)
    G = length(ctx.regions) + 1
    n_all, dry_n = zeros(Int, 5, G), zeros(Int, 5, G)
    dry_bin = zeros(Int, 5, G, length(FALSE_ALARM_BINS))
    dry_volume, volume = zeros(5, G), zeros(5, G)
    sat_wet, sat_wet_res = zeros(Int, 5, G), zeros(Int, 5, G)
    class_total, class_dry = zeros(Int, 5, G, 6), zeros(Int, 5, G, 6)
    for j in axes(Y_obs, 2)
        s = ctx.column_season[j] + 1
        for i in axes(Y_obs, 1)
            mask[i, j] || continue
            obs, sat = Y_obs[i, j], Y_sat[i, j]
            g = ctx.station_region_index[i] + 1
            dry = obs < WET_MM
            c = intensity_class(sat) + 1
            b = _false_alarm_bin(sat)
            for (a, r) in ((1, 1), (1, g), (s, 1), (s, g))
                n_all[a, r] += 1
                volume[a, r] += sat
                sat_wet[a, r] += sat >= WET_MM
                sat_wet_res[a, r] += sat >= GAUGE_RESOLUTION_MM
                class_total[a, r, c] += 1
                if dry
                    dry_n[a, r] += 1
                    dry_bin[a, r, b] += 1
                    dry_volume[a, r] += sat
                    class_dry[a, r, c] += 1
                end
            end
        end
    end
    ratio(a, b) = b > 0 ? a / b : NaN
    rows = NamedTuple[]
    for a in 1:5, r in 1:G
        false_alarms = dry_n[a, r] - dry_bin[a, r, 1]
        false_alarms_res = sum(dry_bin[a, r, 3:end])
        pairs = Pair{Symbol,Any}[
            :sample => sample, :product => product, :season => SEASON_LABELS[a], :region => _region_labels(ctx)[r],
            :n_cells => n_all[a, r], :n_gauge_dry => dry_n[a, r], :gauge_dry_share => ratio(dry_n[a, r], n_all[a, r]),
            :false_alarm_rate => ratio(false_alarms, dry_n[a, r]),
            :false_alarm_rate_res => ratio(false_alarms_res, dry_n[a, r]),
        ]
        for (k, bin) in enumerate(FALSE_ALARM_BINS)
            k == 1 && continue
            push!(pairs, Symbol("dry_share_", bin) => ratio(dry_bin[a, r, k], dry_n[a, r]))
        end
        append!(pairs, [
            :FAR => ratio(false_alarms, sat_wet[a, r]), :FAR_res => ratio(false_alarms_res, sat_wet_res[a, r]),
            :dry_volume_share => ratio(dry_volume[a, r], volume[a, r]),
        ])
        for c in 2:6
            push!(pairs, Symbol("gauge_dry_given_", INTENSITY_CLASSES[c]) => ratio(class_dry[a, r, c], class_total[a, r, c]))
        end
        push!(rows, _row(pairs))
    end
    return DataFrame(rows)
end

"""
Rain events in one gauge's hourly series.

A wet hour is `>= wet`. Consecutive wet hours belong to one event unless at least `min_dry_gap`
non-wet hours separate them. Each event is `(; start, stop, peak_index, peak_mm, total_mm, complete)`,
indices into `y`; the peak is the first hour holding the maximum. `complete` is false when a missing
hour falls inside the event or within `min_dry_gap` hours of either end, or the axis ends there - the
event's extent cannot then be confirmed, so it is not scored.
"""
function station_rain_events(y::AbstractVector{<:Real}; wet::Real=WET_MM, min_dry_gap::Int=3)
    min_dry_gap >= 1 || throw(ArgumentError("min_dry_gap must be at least 1"))
    T = length(y)
    events = NamedTuple{(:start, :stop, :peak_index, :peak_mm, :total_mm, :complete),
        Tuple{Int,Int,Int,Float64,Float64,Bool}}[]
    wet_hours = findall(v -> isfinite(v) && v >= wet, y)
    isempty(wet_hours) && return events
    function close_event(start, stop)
        values = @view y[start:stop]
        peak_index = start - 1 + argmax([isfinite(v) ? v : -Inf for v in values])
        lo, hi = start - min_dry_gap, stop + min_dry_gap
        complete = lo >= 1 && hi <= T && all(isfinite, @view y[lo:hi])
        push!(events, (start, stop, peak_index, Float64(y[peak_index]), sum(v for v in values if isfinite(v)), complete))
    end
    start = previous = wet_hours[1]
    for t in @view wet_hours[2:end]
        if t - previous - 1 >= min_dry_gap
            close_event(start, previous)
            start = t
        end
        previous = t
    end
    close_event(start, previous)
    return events
end

"""
The lag in hours (`-max_lag..max_lag`) at which `y` correlates best with `x`; positive means `y` is late.

At lag `L`, `x[t]` is paired with `y[t + L]` over the overlap, so a satellite series that repeats the
gauge two hours later peaks at `L = +2`. Ties go to the smaller `|L|`, then to the earlier lag; `NaN`
when no lag has `min_overlap` finite, non-constant pairs.
"""
function best_lag(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; max_lag::Int, min_overlap::Int=4)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    best, best_r = NaN, -Inf
    for lag in sort(collect(-max_lag:max_lag); by=l -> (abs(l), l))
        t = max(1, 1 - lag):min(length(x), length(x) - lag)
        isempty(t) && continue
        r = _finite_cor(view(x, t), view(y, t .+ lag); min_n=min_overlap)
        if isfinite(r) && r > best_r
            best, best_r = Float64(lag), r
        end
    end
    return best
end

"""
Scores of one satellite series against one gauge event, over the event's scoring window.

`obs` and `sat` are the window (finite), and `start`/`stop` the event's first and last wet hour inside it.
Timing errors are satellite minus gauge in hours (positive = satellite late):
- `peak_error_h`: first satellite maximum minus first gauge maximum (tipping-bucket ties resolve to the
  earliest hour)
- `centroid_error_h`: difference of the rain-weighted mean hours - robust where peaks are flat
- `onset_error_h` / `end_error_h` / `duration_sat_h`: first and last satellite hour `>= wet`, censored by
  the window (an onset cannot be earlier than the window start); `_res` repeats them at `resolution`
Timing is `NaN` when the satellite has no hour `>= wet` (`detected` false). `r_event` needs 4 hours and
`best_lag_h` 4 overlapping hours.
"""
function score_event(
    obs::AbstractVector{<:Real}, sat::AbstractVector{<:Real}, start::Int, stop::Int;
    wet::Real=WET_MM, resolution::Real=GAUGE_RESOLUTION_MM, max_lag::Int=3,
)
    length(obs) == length(sat) || throw(DimensionMismatch("obs and sat must have the same length"))
    1 <= start <= stop <= length(obs) || throw(ArgumentError("event must lie inside the window"))
    all(isfinite, obs) && all(isfinite, sat) || throw(ArgumentError("the scoring window must be finite"))
    hours = 1:length(obs)
    obs_total, sat_total = sum(obs), sum(sat)
    function timing(threshold)
        wet_hours = findall(>=(threshold), sat)
        isempty(wet_hours) && return (false, NaN, NaN, NaN)
        return (true, Float64(first(wet_hours) - start), Float64(last(wet_hours) - stop),
            Float64(last(wet_hours) - first(wet_hours) + 1))
    end
    detected, onset, finish, duration = timing(wet)
    detected_res, onset_res, finish_res, duration_res = timing(resolution)
    return (;
        window_hours=length(obs), duration_obs_h=stop - start + 1,
        obs_total_mm=obs_total, sat_total_mm=sat_total, volume_rel_bias=(sat_total - obs_total) / obs_total,
        obs_peak_mm=maximum(obs), sat_peak_mm=maximum(sat), peak_ratio=maximum(sat) / maximum(obs),
        detected, peak_error_h=detected ? Float64(argmax(sat) - argmax(obs)) : NaN,
        centroid_error_h=detected ? sum(hours .* sat) / sat_total - sum(hours .* obs) / obs_total : NaN,
        onset_error_h=onset, end_error_h=finish, duration_sat_h=duration,
        detected_res, onset_error_h_res=onset_res, end_error_h_res=finish_res, duration_sat_h_res=duration_res,
        r_event=_finite_cor(obs, sat; min_n=4), best_lag_h=best_lag(obs, sat; max_lag, min_overlap=4),
    )
end

"""
Every gauge event, and each product's scores over it.

Events come from `station_rain_events` with `min_dry_gap`. The scoring window is `pad` hours either
side of the event, cut back so it never reaches another event's wet hours. An event is scored for a
product when the gauge event is `complete`, the gauge is finite over the window, and the product is
too. `samples` maps a sample name to its products; `in_<sample>` marks the events every product of
that sample can score, which is the set that sample summarises - so within a sample all products are
scored on the same events.

Returns `(; events, scores)`: one row per gauge event, and one per scored event x product carrying the
event's strata.
"""
function event_tables(
    ctx, Y_obs::AbstractMatrix, products::AbstractDict{<:AbstractString,<:AbstractMatrix},
    samples::AbstractVector{<:Pair}, station_ids::AbstractVector{<:AbstractString};
    min_dry_gap::Int=3, pad::Int=3, max_lag::Int=3,
)
    _check_inputs(ctx, Y_obs, values(products)...)
    T = size(Y_obs, 2)
    product_names = sort(collect(keys(products)))
    event_rows, score_rows = NamedTuple[], NamedTuple[]
    event_id = 0
    for i in axes(Y_obs, 1)
        y = @view Y_obs[i, :]
        events = station_rain_events(y; min_dry_gap)
        region = _region_labels(ctx)[ctx.station_region_index[i] + 1]
        for (k, event) in enumerate(events)
            event_id += 1
            lo = max(1, event.start - pad, k > 1 ? events[k - 1].stop + 1 : 1)
            hi = min(T, event.stop + pad, k < length(events) ? events[k + 1].start - 1 : T)
            gauge_ok = event.complete && all(isfinite, @view y[lo:hi])
            available = Dict(p => gauge_ok && all(isfinite, @view products[p][i, lo:hi]) for p in product_names)
            in_sample = [Symbol("in_", name) => all(p -> available[p], members) for (name, members) in samples]
            strata = (; event_id, station_id=String(station_ids[i]), region,
                season=SEASONS[ctx.column_season[event.start]], day=ctx.column_day[event.start],
                event_class=INTENSITY_CLASSES[intensity_class(event.peak_mm) + 1])
            push!(event_rows, merge(strata, (;
                start_time=ctx.times[event.start], stop_time=ctx.times[event.stop],
                duration_h=event.stop - event.start + 1, total_mm=event.total_mm, peak_mm=event.peak_mm,
                gauge_complete=gauge_ok, window_start=ctx.times[lo], window_stop=ctx.times[hi],
            ), _row(in_sample)))
            gauge_ok || continue
            for p in product_names
                available[p] || continue
                scores = score_event(Y_obs[i, lo:hi], products[p][i, lo:hi], event.start - lo + 1, event.stop - lo + 1; max_lag)
                push!(score_rows, merge(strata, (; product=p), _row(in_sample), scores))
            end
        end
    end
    return (; events=DataFrame(event_rows), scores=DataFrame(score_rows))
end

"""
Event scores summarised per sample x product x season x region x event class (`all` pools classes).

`n_gauge_events` counts the stratum's complete gauge events and `retained_share` the fraction the
sample could score - where FY4B's missing hours bite. Among scored events: `POD_event` (`_res` at the
gauge resolution); medians of peak, centroid, onset and end errors over detected events, the peak
error's IQR and `peak_within_1h` share; median duration ratio, volume bias, peak ratio and event r;
the pooled volume bias (sum over events); the shares of best lags that are early, zero or late.
`_lo`/`_hi` are 95% day-block bootstrap bounds for `POD_event`, `peak_within_1h` and
`pooled_volume_rel_bias` (additive statistics; medians get none). `low_sample` flags < 30 events.
"""
function event_summary(
    ctx, events::DataFrame, scores::DataFrame, samples::AbstractVector{<:Pair};
    reps::Int=1000, seed::Int=20260913, bootstrap::Bool=true,
)
    region_labels = _region_labels(ctx)
    classes = ["all"; INTENSITY_CLASSES[2:end]]
    in_stratum(season, region, class, s, g, c) = (s == "all" || season == s) && (g == "all" || region == g) &&
        (c == "all" || class == c)
    # Typed copies of the stratum columns: DataFrame column access inside these loops is not inferred.
    event_complete, event_season = Vector{Bool}(events.gauge_complete), String.(events.season)
    event_region, event_class = String.(events.region), String.(events.event_class)
    score_season, score_region, score_class = String.(scores.season), String.(scores.region), String.(scores.event_class)
    rows = NamedTuple[]
    stratum = 0
    for (sample, members) in samples, product in members
        product_rows = findall((scores.product .== product) .& scores[!, Symbol("in_", sample)])
        for s in SEASON_LABELS, g in region_labels, c in classes
            stratum += 1
            n_gauge = count(i -> event_complete[i] && in_stratum(event_season[i], event_region[i], event_class[i], s, g, c),
                eachindex(event_complete))
            sel = filter(i -> in_stratum(score_season[i], score_region[i], score_class[i], s, g, c), product_rows)
            n = length(sel)
            detected = scores.detected[sel]
            hit = sel[detected]
            peak = scores.peak_error_h[hit]
            lags = filter(isfinite, scores.best_lag_h[sel])
            base = (;
                sample, product, season=s, region=g, event_class=c,
                n_gauge_events=n_gauge, n_events=n, retained_share=n_gauge > 0 ? n / n_gauge : NaN,
                n_stations=length(unique(scores.station_id[sel])),
                POD_event=_share(detected), POD_event_res=_share(scores.detected_res[sel]),
                median_peak_error_h=_median(peak), peak_error_q25_h=_quantile(peak, 0.25),
                peak_error_q75_h=_quantile(peak, 0.75), peak_within_1h=_share(abs.(peak) .<= 1),
                median_centroid_error_h=_median(scores.centroid_error_h[hit]),
                median_onset_error_h=_median(scores.onset_error_h[hit]),
                median_end_error_h=_median(scores.end_error_h[hit]),
                median_onset_error_h_res=_median(filter(isfinite, scores.onset_error_h_res[sel])),
                median_end_error_h_res=_median(filter(isfinite, scores.end_error_h_res[sel])),
                median_duration_ratio=_median(scores.duration_sat_h[hit] ./ scores.duration_obs_h[hit]),
                median_volume_rel_bias=_median(scores.volume_rel_bias[sel]),
                pooled_volume_rel_bias=n > 0 ? (sum(scores.sat_total_mm[sel]) - sum(scores.obs_total_mm[sel])) / sum(scores.obs_total_mm[sel]) : NaN,
                median_peak_ratio=_median(scores.peak_ratio[sel]),
                median_r_event=_median(filter(isfinite, scores.r_event[sel])),
                lag_early_share=_share(lags .< 0), lag_zero_share=_share(lags .== 0), lag_late_share=_share(lags .> 0),
                low_sample=n < 30,
            )
            ci = nothing
            if bootstrap && n > 0
                within = detected .& (abs.(scores.peak_error_h[sel]) .<= 1)
                sums = [ones(n) detected within scores.obs_total_mm[sel] scores.sat_total_mm[sel]]
                ci = day_block_bootstrap(scores.day[sel], sums, t -> (t[2] / t[1], t[3] / t[2], (t[5] - t[4]) / t[4]);
                    reps, rng=MersenneTwister(seed + stratum))
            end
            bounds = ci === nothing ? fill(NaN, 6) : [ci.lo[1], ci.hi[1], ci.lo[2], ci.hi[2], ci.lo[3], ci.hi[3]]
            push!(rows, merge(base, NamedTuple{(:POD_event_lo, :POD_event_hi, :peak_within_1h_lo, :peak_within_1h_hi,
                :pooled_volume_rel_bias_lo, :pooled_volume_rel_bias_hi)}(Tuple(bounds))))
        end
    end
    return DataFrame(rows)
end

"""
First-harmonic fit `mean + amplitude * cos(2π(h - phase_hour)/24)` of a diurnal curve by least squares,
so hours with no data can simply be left out.
"""
function harmonic_phase(hours::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    length(hours) == length(values) || throw(DimensionMismatch("hours and values must have the same length"))
    length(hours) >= 3 || throw(ArgumentError("need at least three hours"))
    omega = 2π .* Float64.(hours) ./ 24
    beta = [ones(length(omega)) cos.(omega) sin.(omega)] \ Float64.(values)
    return (; mean=beta[1], amplitude=hypot(beta[2], beta[3]), phase_hour=mod(atan(beta[3], beta[2]) * 24 / 2π, 24))
end

"""`a - b` on the 24-hour clock, in `[-12, 12)`."""
circular_hour_difference(a::Real, b::Real) = mod(a - b + 12, 24) - 12

"""
Composite diurnal cycles per sample x source x season x region x hour of day (hour-ending BJT).

Over the sample's cells: `mean_amount` (every hour, dry ones included), `wet_freq` (share of hours
`>= 0.1`), `wet_intensity` (mean over wet hours) and `freq_<class>` for each rain class. `sources` is
a vector of `name => matrix`, the gauge first by convention.
"""
function diurnal_cycle_table(ctx, sources::AbstractVector{<:Pair}, mask::AbstractMatrix; sample::AbstractString)
    G = length(ctx.regions) + 1
    rows = NamedTuple[]
    for (source, Y) in sources
        _check_inputs(ctx, Y, mask)
        n, n_wet = zeros(Int, 5, G, 24), zeros(Int, 5, G, 24)
        amount, wet_amount = zeros(5, G, 24), zeros(5, G, 24)
        class_n = zeros(Int, 5, G, 24, 5)
        for j in axes(Y, 2)
            s, h = ctx.column_season[j] + 1, ctx.column_hour[j] + 1
            for i in axes(Y, 1)
                mask[i, j] || continue
                v = Y[i, j]
                g = ctx.station_region_index[i] + 1
                c = intensity_class(v)
                for (a, b) in ((1, 1), (1, g), (s, 1), (s, g))
                    n[a, b, h] += 1
                    amount[a, b, h] += v
                    if c >= 1
                        n_wet[a, b, h] += 1
                        wet_amount[a, b, h] += v
                        class_n[a, b, h, c] += 1
                    end
                end
            end
        end
        for a in 1:5, b in 1:G, h in 1:24
            count_ = n[a, b, h]
            pairs = Pair{Symbol,Any}[
                :sample => sample, :source => String(source), :season => SEASON_LABELS[a],
                :region => _region_labels(ctx)[b], :hour => h - 1, :n => count_,
                :mean_amount => count_ > 0 ? amount[a, b, h] / count_ : NaN,
                :wet_freq => count_ > 0 ? n_wet[a, b, h] / count_ : NaN,
                :wet_intensity => n_wet[a, b, h] > 0 ? wet_amount[a, b, h] / n_wet[a, b, h] : NaN,
            ]
            for c in 1:5
                push!(pairs, Symbol("freq_", INTENSITY_CLASSES[c + 1]) => count_ > 0 ? class_n[a, b, h, c] / count_ : NaN)
            end
            push!(rows, _row(pairs))
        end
    end
    return DataFrame(rows)
end

"""
Each product's diurnal cycle against the gauge's, per sample x season x region.

Over the hours with data: r of the amount and frequency curves, peak hours and their circular
difference, first-harmonic phases and their difference, the relative amplitude (amplitude / mean) of
each and its ratio, and the ratio of mean amounts. `NaN` below 12 hours or 100 cells an hour.
"""
function diurnal_summary(diurnal::DataFrame, products::AbstractVector{<:AbstractString}; gauge::AbstractString="Gauge")
    rows = NamedTuple[]
    for group in groupby(diurnal, [:sample, :season, :region]; sort=true)
        reference = sort(group[group.source .== gauge, :], :hour)
        for product in products
            estimate = sort(group[group.source .== product, :], :hour)
            nrow(estimate) == nrow(reference) || continue
            keep = findall((reference.n .>= 100) .& (estimate.n .>= 100))
            hours = reference.hour[keep]
            enough = length(keep) >= 12 && sum(reference.mean_amount[keep]) > 0 && sum(estimate.mean_amount[keep]) > 0
            fit(values) = enough ? harmonic_phase(hours, values) : (; mean=NaN, amplitude=NaN, phase_hour=NaN)
            obs_amount, sat_amount = fit(reference.mean_amount[keep]), fit(estimate.mean_amount[keep])
            obs_freq, sat_freq = fit(reference.wet_freq[keep]), fit(estimate.wet_freq[keep])
            peak(values) = enough ? Float64(hours[argmax(values)]) : NaN
            obs_peak, sat_peak = peak(reference.mean_amount[keep]), peak(estimate.mean_amount[keep])
            push!(rows, (;
                sample=group.sample[1], product, season=group.season[1], region=group.region[1], n_hours=length(keep),
                r_amount=enough ? _finite_cor(reference.mean_amount[keep], estimate.mean_amount[keep]) : NaN,
                r_freq=enough ? _finite_cor(reference.wet_freq[keep], estimate.wet_freq[keep]) : NaN,
                peak_hour_gauge=obs_peak, peak_hour_sat=sat_peak,
                peak_hour_diff_h=enough ? circular_hour_difference(sat_peak, obs_peak) : NaN,
                phase_hour_gauge=obs_amount.phase_hour, phase_hour_sat=sat_amount.phase_hour,
                phase_diff_h=enough ? circular_hour_difference(sat_amount.phase_hour, obs_amount.phase_hour) : NaN,
                freq_phase_diff_h=enough ? circular_hour_difference(sat_freq.phase_hour, obs_freq.phase_hour) : NaN,
                rel_amplitude_gauge=obs_amount.amplitude / obs_amount.mean,
                rel_amplitude_sat=sat_amount.amplitude / sat_amount.mean,
                rel_amplitude_ratio=(sat_amount.amplitude / sat_amount.mean) / (obs_amount.amplitude / obs_amount.mean),
                mean_amount_ratio=obs_amount.mean > 0 ? sat_amount.mean / obs_amount.mean : NaN,
            ))
        end
    end
    return DataFrame(rows)
end

"""
Hourly agreement at each gauge, per sample x product x season x station.

`r_union_wet` / `rho_union_wet` use only hours where the gauge or the satellite is `>= 0.1`, so the
mass of shared dry hours does not inflate them; `r_all_hours` is the conventional all-hours value, for
reference. `NaN` below `min_n` hours.
"""
function station_correlation_table(
    ctx, Y_obs::AbstractMatrix, Y_sat::AbstractMatrix, mask::AbstractMatrix;
    sample::AbstractString, product::AbstractString, station_ids::AbstractVector{<:AbstractString}, min_n::Int=50,
)
    _check_inputs(ctx, Y_obs, Y_sat, mask)
    rows = NamedTuple[]
    for i in axes(Y_obs, 1), s in 0:4
        columns = [j for j in axes(Y_obs, 2) if mask[i, j] && (s == 0 || ctx.column_season[j] == s)]
        obs, sat = Y_obs[i, columns], Y_sat[i, columns]
        wet = (obs .>= WET_MM) .| (sat .>= WET_MM)
        n_wet = count(wet)
        push!(rows, (;
            sample, product, season=SEASON_LABELS[s + 1], station_id=String(station_ids[i]),
            region=_region_labels(ctx)[ctx.station_region_index[i] + 1],
            n_hours=length(columns), n_union_wet=n_wet,
            r_union_wet=n_wet >= min_n ? _finite_cor(obs[wet], sat[wet]; min_n) : NaN,
            rho_union_wet=n_wet >= min_n ? _spearman(obs[wet], sat[wet]; min_n) : NaN,
            r_all_hours=length(columns) >= min_n ? _finite_cor(obs, sat; min_n) : NaN,
        ))
    end
    return DataFrame(rows)
end

"""Median and IQR of the station correlations per sample x product x season x region (`all` included)."""
function station_correlation_summary(stations::DataFrame)
    rows = NamedTuple[]
    table = vcat(stations, transform(stations, :region => (r -> fill("all", length(r))) => :region))
    for group in groupby(table, [:sample, :product, :season, :region]; sort=true)
        r = filter(isfinite, group.r_union_wet)
        rho = filter(isfinite, group.rho_union_wet)
        r_all = filter(isfinite, group.r_all_hours)
        push!(rows, (;
            sample=group.sample[1], product=group.product[1], season=group.season[1], region=group.region[1],
            n_stations=length(r), median_r_union_wet=_median(r), r_union_wet_q25=_quantile(r, 0.25),
            r_union_wet_q75=_quantile(r, 0.75), median_rho_union_wet=_median(rho), median_r_all_hours=_median(r_all),
            low_sample=length(r) < 10,
        ))
    end
    return DataFrame(rows)
end

"""
Hourly mean over a region's stations (`region = 0` for all), from the cells in `mask`.

An hour is kept only when at least `min_station_fraction` of the region's stations are in the mask, so
the average is not taken over a handful of gauges; `NaN` otherwise. Because a sample's mask is shared by
its sources, every source is averaged over the same stations each hour.
"""
function regional_mean_series(ctx, Y::AbstractMatrix, mask::AbstractMatrix, region::Int; min_station_fraction::Real=0.9)
    _check_inputs(ctx, Y, mask)
    stations = region == 0 ? collect(axes(Y, 1)) : findall(==(region), ctx.station_region_index)
    series = fill(NaN, size(Y, 2))
    isempty(stations) && return series
    for j in axes(Y, 2)
        total, n = 0.0, 0
        for i in stations
            mask[i, j] || continue
            total += Y[i, j]
            n += 1
        end
        n > 0 && n >= min_station_fraction * length(stations) && (series[j] = total / n)
    end
    return series
end

"""
Accumulate an hourly series to `timescale_hours` = 1, 3 or 24.

3-h blocks and 24-h totals follow the met-day (the first block starts at the day's first hour); a block
missing any hour is `NaN`. Returns `(; values, day)` with each block's met-day index.
"""
function aggregate_series(ctx, x::AbstractVector{<:Real}, timescale_hours::Int)
    timescale_hours in (1, 3, 24) || throw(ArgumentError("timescale must be 1, 3 or 24 hours"))
    length(x) == length(ctx.times) || throw(DimensionMismatch("series must match the time axis"))
    timescale_hours == 1 && return (; values=Float64.(x), day=copy(ctx.column_day))
    key(j) = (ctx.column_day[j], ctx.column_offset[j] ÷ timescale_hours)
    values, day = Float64[], Int[]
    j = 1
    while j <= length(x)
        stop = j
        while stop < length(x) && key(stop + 1) == key(j)
            stop += 1
        end
        block = @view x[j:stop]
        push!(values, stop - j + 1 == timescale_hours && all(isfinite, block) ? sum(block) : NaN)
        push!(day, ctx.column_day[j])
        j = stop + 1
    end
    return (; values, day)
end

"""
Agreement of two series over their finite pairs: n, means, RMSE, MAE, Bias, RB_pct, r, rho, sd ratio,
KGE, POD/FAR/CSI at `wet`, r over pairs where either side is wet, and - with `max_lag > 0` - the best
lag (positive = `sat` late). All `NaN` below 10 pairs.
"""
function series_metrics(obs::AbstractVector{<:Real}, sat::AbstractVector{<:Real}; wet::Real=WET_MM, max_lag::Int=0)
    length(obs) == length(sat) || throw(DimensionMismatch("obs and sat must have the same length"))
    pairs = findall(i -> isfinite(obs[i]) && isfinite(sat[i]), eachindex(obs))
    n = length(pairs)
    names = (:n, :obs_mean, :sat_mean, :RMSE, :MAE, :Bias, :RB_pct, :r, :rho, :sd_ratio, :KGE, :POD, :FAR, :CSI,
        :n_wet, :r_wet, :best_lag_h)
    n < 10 && return NamedTuple{names}((n, fill(NaN, length(names) - 1)...))
    o, s = Float64.(obs[pairs]), Float64.(sat[pairs])
    continuous = metric_continuous(o, s)
    events = metric_event(o, s; thr=Float64(wet))
    sd_ratio = std(s; corrected=false) / std(o; corrected=false)
    kge = 1 - sqrt((continuous.r - 1)^2 + (sd_ratio - 1)^2 + (mean(s) / mean(o) - 1)^2)
    wet_pairs = (o .>= wet) .| (s .>= wet)
    return NamedTuple{names}((
        n, mean(o), mean(s), continuous.RMSE, continuous.MAE, continuous.Bias, 100 * (sum(s) - sum(o)) / sum(o),
        continuous.r, _spearman(o, s), sd_ratio, kge, events.POD, events.FAR, events.CSI,
        count(wet_pairs), _finite_cor(o[wet_pairs], s[wet_pairs]; min_n=10),
        max_lag > 0 ? best_lag(obs, sat; max_lag, min_overlap=10) : NaN,
    ))
end

"""
Regional-mean series metrics per sample x product x region x season x timescale (1, 3, 24 h).

Series come from `regional_mean_series` over the sample's mask. A season keeps the periods whose
met-day falls in it (for the 1-h lag search the gauge side is masked and the satellite side left
whole, so a lag can reach across the edge of a season by at most `max_lag` hours). Rows with no scored
period are omitted - FY4B never completes a met-day, so `all_products` has no 24-h rows.
"""
function regional_series_table(
    ctx, Y_obs::AbstractMatrix, products::AbstractDict{<:AbstractString,<:AbstractMatrix}, mask::AbstractMatrix;
    sample::AbstractString, max_lag::Int=6, min_station_fraction::Real=0.9,
)
    rows = NamedTuple[]
    for (g, region) in enumerate(_region_labels(ctx))
        obs_hourly = regional_mean_series(ctx, Y_obs, mask, g - 1; min_station_fraction)
        for product in sort(collect(keys(products)))
            sat_hourly = regional_mean_series(ctx, products[product], mask, g - 1; min_station_fraction)
            for timescale in (1, 3, 24)
                obs = aggregate_series(ctx, obs_hourly, timescale)
                sat = aggregate_series(ctx, sat_hourly, timescale)
                for (s, season) in enumerate(SEASON_LABELS)
                    in_season = [s == 1 || ctx.day_season[d] == s - 1 for d in obs.day]
                    masked_obs = [in_season[k] ? obs.values[k] : NaN for k in eachindex(obs.values)]
                    scores = series_metrics(masked_obs, sat.values; max_lag=timescale == 1 ? max_lag : 0)
                    scores.n == 0 && continue
                    push!(rows, merge((; sample, product, region, season, timescale_h=timescale),
                        scores, (; low_sample=scores.n < 100)))
                end
            end
        end
    end
    return DataFrame(rows)
end

end # module
