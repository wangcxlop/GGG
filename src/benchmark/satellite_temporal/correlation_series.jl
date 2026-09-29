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
