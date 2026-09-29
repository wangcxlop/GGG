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
- `peak_error_h`: first satellite maximum minus first gauge maximum (ties, common in the 0.5 mm-quantised
  gauge record, resolve to the earliest hour)
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
