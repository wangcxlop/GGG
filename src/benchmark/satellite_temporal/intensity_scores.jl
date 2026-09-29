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
