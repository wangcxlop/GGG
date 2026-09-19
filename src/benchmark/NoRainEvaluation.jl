"""
No-rain evaluation: what the satellite products, and the interpolation and fusion methods built on
them, report on the station-hours where the gauge records no rain.

About 93% of gauge station-hours are dry, so dry-hour behaviour decides most pooled scores, yet
`SatelliteTemporalEvaluation.false_alarm_table` is the only place it has been looked at on its own.
This module widens that view:

- `occurrence_table`: the full occurrence contingency and the dry-hour amount scores of one estimate,
  per level of any stratifier, with day-block bootstrap intervals. Every stratification below (season,
  hour, station, distance to gauge rain, neighbouring gauges, the network's state, meteorology, product
  agreement) goes through this one function, so they are all scored the same way.
- Stratifiers of the gauge-dry hours: `hours_to_rain` / `proximity_labels` / `side_labels` (where a dry
  hour sits relative to the station's own rain), `neighbour_wet_share` / `neighbour_labels` (whether
  nearby gauges are wet), `network_wet_share` / `network_labels` (whether the whole network is dry),
  and `bin_labels` for continuous covariates.
- `neighbour_gauge_baseline`: how often a *neighbouring gauge* reports rain when a station is dry, by
  separation - the false-alarm rate a perfect areal estimate would still show against a point gauge.
- `run_lengths` / `dry_spell_table` / `daily_sums`: dry spells, hourly and daily.
- `silent_gauge_spells`: the other side of a false alarm - long runs of gauge zeros during which the
  neighbouring gauges were clearly wet, i.e. a gauge reporting 0.0 while it was not measuring.
- `paired_error_table` / `zero_threshold_sensitivity`: dry-hour error of one prediction against a
  reference, and what zeroing small values would change.
- `long_to_wide_hourly` / `lagged_correlation`: bring a long-format product onto the gauge grid and
  check its clock.

Conventions, shared with `SatelliteTemporalEvaluation`:

- Matrices are `[station, time]`, `NaN` for missing, every estimate aligned to the gauges.
- `mask` names the cells a sample scores; a cell outside it is ignored everywhere.
- A gauge-dry cell is `gauge < WET_MM` (0.1 mm/h). The gauges are 0.5 mm tipping buckets, so that is a
  0.0 reading, and a satellite value between 0.1 and 0.5 mm/h may be rain the bucket had not yet tipped.
- A stratifier is an integer label per cell - a matrix, or a function `(i, j) -> label` for labels that
  only depend on the station or the hour. Labels run 1..L, and 0 leaves a cell out.
- `day[j]` is column `j`'s met-day number (positive integers); intervals resample whole days.
"""
module NoRainEvaluation

using DataFrames, Dates, Random, Statistics
using Main.SatelliteTemporalEvaluation: WET_MM, GAUGE_RESOLUTION_MM, day_block_bootstrap
using Main.TraditionalInterpolation: haversine_distance_matrix

export DRY_THRESHOLDS, EXCEEDANCE_THRESHOLDS, HOURS_PER_YEAR
export PROXIMITY_EDGES, PROXIMITY_LEVELS, SIDE_LEVELS, NEIGHBOUR_LEVELS, NETWORK_EDGES, NETWORK_LEVELS
export SPELL_BIN_EDGES
export hours_to_rain, proximity_labels, side_labels, neighbour_wet_share, neighbour_labels
export network_wet_share, network_labels, bin_labels, phantom_extent
export occurrence_table, paired_error_table, neighbour_gauge_baseline
export run_lengths, dry_spell_table, daily_sums, silent_gauge_spells
export long_to_wide_hourly, lagged_correlation, zero_threshold_sensitivity

"""Estimate thresholds for occurrence: the wet threshold and the gauge's 0.5 mm resolution."""
const DRY_THRESHOLDS = [WET_MM, GAUGE_RESOLUTION_MM]
"""Thresholds of the exceedance curve P(estimate >= t | gauge dry)."""
const EXCEEDANCE_THRESHOLDS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]
const HOURS_PER_YEAR = 8766.0
"""Upper bounds (h) of the distance-to-rain classes; the last class is open."""
const PROXIMITY_EDGES = [1, 3, 6, 24]
const PROXIMITY_LEVELS = ["1 h", "2-3 h", "4-6 h", "7-24 h", "> 24 h"]
const SIDE_LEVELS = ["inside an event gap", "before rain", "after rain", "away from rain"]
const NEIGHBOUR_LEVELS = ["all neighbours dry", "under half wet", "half or more wet"]
"""Upper bounds of the network wet share for the classes after `network dry`."""
const NETWORK_EDGES = [0.05, 0.25]
const NETWORK_LEVELS = ["network dry", "isolated (<= 5% wet)", "scattered (5-25% wet)", "widespread (> 25% wet)"]
"""Lower bounds (h, or days on a daily axis) of the spell-length histogram bins."""
const SPELL_BIN_EDGES = [1, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 72, 96, 144, 192, 288, 384, 576, 768]

_ratio(a, b) = b > 0 ? a / b : NaN
_label(labels::AbstractMatrix, i, j) = labels[i, j]
_label(labels::Function, i, j) = labels(i, j)

# ---------------------------------------------------------------------------------------------
# Where a dry hour sits

"""
Hours from each hour of one gauge series to its nearest wet hour (`>= wet`), looking back (`before`)
and ahead (`after`).

A distance is positive when a wet hour is reached, `Inf` when the series ends first, and negative
(minus the distance) when a missing hour is reached first - the rain could be hiding there, so the
true distance is only bounded. Wet hours are at distance 0; missing hours carry `NaN`. The matrix
method applies this to every row.
"""
function hours_to_rain(y::AbstractVector{<:Real}; wet::Real=WET_MM)
    n = length(y)
    before, after = fill(NaN, n), fill(NaN, n)
    for (range, out) in ((1:n, before), (n:-1:1, after))
        last, last_missing = 0, false
        for t in range
            v = y[t]
            if !isfinite(v)
                last, last_missing = t, true
            elseif v >= wet
                out[t] = 0.0
                last, last_missing = t, false
            else
                out[t] = last == 0 ? Inf : last_missing ? -abs(t - last) : abs(t - last)
            end
        end
    end
    return before, after
end

function hours_to_rain(Y::AbstractMatrix{<:Real}; wet::Real=WET_MM)
    before, after = fill(NaN, size(Y)), fill(NaN, size(Y))
    for i in axes(Y, 1)
        b, a = hours_to_rain(view(Y, i, :); wet)
        before[i, :] = b
        after[i, :] = a
    end
    return before, after
end

"""
Class of the distance from a dry hour to the station's nearest rain (`PROXIMITY_LEVELS`), 0 for wet or
missing hours and for hours whose class a nearby missing stretch leaves undecided.

With the nearest known wet hour at `dw` and the nearest missing hour at `dm`, the true distance lies in
`[min(dw, dm), dw]`; the class is kept only when both ends fall in the same class.
"""
function proximity_labels(before::AbstractMatrix, after::AbstractMatrix; edges=PROXIMITY_EDGES)
    size(before) == size(after) || throw(DimensionMismatch("before and after must have the same size"))
    labels = zeros(Int8, size(before))
    for k in eachindex(before)
        b, a = before[k], after[k]
        (isnan(b) || isnan(a) || b == 0 || a == 0) && continue
        dw = min(b > 0 ? b : Inf, a > 0 ? a : Inf)
        dm = min(b < 0 ? -b : Inf, a < 0 ? -a : Inf)
        known = searchsortedfirst(edges, dw)
        (dm >= dw || searchsortedfirst(edges, dm) == known) && (labels[k] = known)
    end
    return labels
end

"""
Which side of the station's rain a dry hour is on (`SIDE_LEVELS`), looking `window` hours each way: rain
on both sides (a gap inside an event), rain only ahead (before the onset), rain only behind (after the
end), or neither. 0 for wet or missing hours, and when a missing hour within the window leaves a side
undecided.
"""
function side_labels(before::AbstractMatrix, after::AbstractMatrix; window::Int=3)
    size(before) == size(after) || throw(DimensionMismatch("before and after must have the same size"))
    labels = zeros(Int8, size(before))
    side(v) = 0 < v <= window ? 1 : -window <= v < 0 ? -1 : 0   # rain, undecided, no rain in the window
    for k in eachindex(before)
        b, a = before[k], after[k]
        (isnan(b) || isnan(a) || b == 0 || a == 0) && continue
        sb, sa = side(b), side(a)
        (sb == -1 || sa == -1) && continue
        labels[k] = sb == 1 && sa == 1 ? 1 : sa == 1 ? 2 : sb == 1 ? 3 : 4
    end
    return labels
end

"""
Share of each station's reporting neighbours within `radius_km` that are wet at each hour; `NaN` where
a station has no neighbour in range or none of them reports.
"""
function neighbour_wet_share(Y::AbstractMatrix{<:Real}, lonlat::AbstractMatrix{<:Real};
        radius_km::Real=15.0, wet::Real=WET_MM)
    size(lonlat, 1) == size(Y, 1) || throw(DimensionMismatch("one lon/lat row per station"))
    distance = haversine_distance_matrix(lonlat, lonlat)
    neighbours = [[k for k in axes(Y, 1) if k != i && distance[k, i] <= radius_km] for i in axes(Y, 1)]
    share = fill(NaN, size(Y))
    for j in axes(Y, 2), i in axes(Y, 1)
        reporting, raining = 0, 0
        for k in neighbours[i]
            v = Y[k, j]
            isfinite(v) || continue
            reporting += 1
            raining += v >= wet
        end
        reporting > 0 && (share[i, j] = raining / reporting)
    end
    return share, neighbours
end

"""`NEIGHBOUR_LEVELS` of the gauge-dry cells from `neighbour_wet_share`; 0 elsewhere."""
function neighbour_labels(share::AbstractMatrix, Y_obs::AbstractMatrix; wet::Real=WET_MM)
    labels = zeros(Int8, size(share))
    for k in eachindex(share)
        s, o = share[k], Y_obs[k]
        (isfinite(o) && o < wet && isfinite(s)) || continue
        labels[k] = s == 0 ? 1 : s < 0.5 ? 2 : 3
    end
    return labels
end

"""
Share of reporting gauges that are wet, per hour; `NaN` for hours where fewer than `min_reporting` of
the stations report, which cannot be called network-dry.
"""
function network_wet_share(Y::AbstractMatrix{<:Real}; min_reporting::Real=0.9, wet::Real=WET_MM)
    share = fill(NaN, size(Y, 2))
    for j in axes(Y, 2)
        reporting = count(isfinite, view(Y, :, j))
        reporting >= min_reporting * size(Y, 1) || continue
        share[j] = count(v -> isfinite(v) && v >= wet, view(Y, :, j)) / reporting
    end
    return share
end

"""`NETWORK_LEVELS` of the gauge-dry cells from a per-hour `network_wet_share`; 0 elsewhere."""
function network_labels(share::AbstractVector, Y_obs::AbstractMatrix; wet::Real=WET_MM, edges=NETWORK_EDGES)
    length(share) == size(Y_obs, 2) || throw(DimensionMismatch("one share per hour"))
    labels = zeros(Int8, size(Y_obs))
    for j in axes(Y_obs, 2)
        s = share[j]
        isfinite(s) || continue
        level = s == 0 ? 1 : searchsortedfirst(edges, s) + 1
        for i in axes(Y_obs, 1)
            o = Y_obs[i, j]
            isfinite(o) && o < wet && (labels[i, j] = level)
        end
    end
    return labels
end

"""
How widely an estimate rains on the hours every reporting gauge is dry (`network_wet_share == 0`).

For each such hour, the share of the estimate's scored cells (inside `mask`) that are `>= wet`; returns
`n_hours` and, for each cutoff `c`, the share of those hours on which that share exceeds `c` (`c = 0`:
at least one gauge pixel rains). With no gauge wet anywhere, rain at a gauge pixel is either rain that
fell between the gauges or rain that never reached the ground.
"""
function phantom_extent(network_share::AbstractVector, Y_est::AbstractMatrix, mask::AbstractMatrix;
        wet::Real=WET_MM, cutoffs=[0.0, 0.05, 0.25, 0.5], key::NamedTuple=(;))
    length(network_share) == size(Y_est, 2) || throw(DimensionMismatch("one share per hour"))
    size(Y_est) == size(mask) || throw(DimensionMismatch("estimate and mask must match"))
    shares = Float64[]
    for j in axes(Y_est, 2)
        network_share[j] == 0 || continue
        scored = [Y_est[i, j] for i in axes(Y_est, 1) if mask[i, j] && isfinite(Y_est[i, j])]
        isempty(scored) || push!(shares, count(>=(wet), scored) / length(scored))
    end
    exceed = (; (Symbol("share_hours_gt_", round(Int, 100c)) => _ratio(count(>(c), shares), length(shares))
        for c in cutoffs)...)
    return merge(key, (; n_hours=length(shares), mean_wet_share=isempty(shares) ? NaN : mean(shares)), exceed)
end

"""Label `k` where `edges[k] <= x < edges[k + 1]`; 0 for `NaN` or out of range."""
function bin_labels(X::AbstractArray{<:Real}, edges::AbstractVector{<:Real})
    issorted(edges) || throw(ArgumentError("edges must be sorted"))
    labels = zeros(Int8, size(X))
    for k in eachindex(X)
        x = X[k]
        isnan(x) && continue
        b = searchsortedlast(edges, x)
        1 <= b < length(edges) && (labels[k] = b)
    end
    return labels
end

# ---------------------------------------------------------------------------------------------
# Occurrence and dry-hour amounts

# Per-day accumulator columns; threshold t adds hits at BASE + 2t - 1 and false alarms at BASE + 2t.
const N_DRY, N_WET, VOL_DRY, VOL_ALL, ZERO_DRY, SSE_DRY, SSE_ALL, VOL_OBS, VOL_OBS_DRY = 1:9
const BASE = 9

"""Scores from one accumulator total; wet-side scores are `NaN` where the level has no gauge-wet cell."""
function _occurrence_scores(t::AbstractVector, n_thresholds::Int)
    n_dry, n_wet = t[N_DRY], t[N_WET]
    out = [_ratio(t[VOL_DRY], n_dry), sqrt(_ratio(t[SSE_DRY], n_dry)), _ratio(t[SSE_DRY], t[SSE_ALL])]
    for k in 1:n_thresholds
        a, b = t[BASE + 2k - 1], t[BASE + 2k]
        c, d = n_wet - a, n_dry - b
        pofd = _ratio(b, n_dry)
        if n_wet > 0
            hss_den = (a + c) * (c + d) + (a + b) * (b + d)
            append!(out, [pofd, _ratio(b, a + b), _ratio(a + b, n_wet), _ratio(d, c + d), _ratio(a, n_wet) - pofd,
                hss_den > 0 ? 2 * (a * d - b * c) / hss_den : NaN])
        else
            append!(out, [pofd, NaN, NaN, NaN, NaN, NaN])
        end
    end
    return out
end

const _THRESHOLD_SCORES = (:POFD, :FAR, :freq_bias, :NPV, :PSS, :HSS)

"""Rows of an occurrence table from per-level, per-day accumulators `acc[level, day, column]`."""
function _occurrence_rows(acc::Array{Float64,3}, levels, thresholds; reps::Int, seed::Int, key::NamedTuple)
    T = length(thresholds)
    rows = NamedTuple[]
    for (l, level) in enumerate(levels)
        per_day = acc[l, :, :]
        used = [d for d in axes(per_day, 1) if per_day[d, N_DRY] + per_day[d, N_WET] > 0]
        total = vec(sum(per_day; dims=1))
        n_dry, n_wet = total[N_DRY], total[N_WET]
        n_cells = n_dry + n_wet
        point = _occurrence_scores(total, T)
        ci = reps > 0 && !isempty(used) ?
            day_block_bootstrap(collect(eachindex(used)), per_day[used, :], t -> _occurrence_scores(t, T);
                reps, rng=MersenneTwister(seed + l)) : nothing
        lo(k) = ci === nothing ? NaN : ci.lo[k]
        hi(k) = ci === nothing ? NaN : ci.hi[k]
        base = (;
            level, n_cells=Int(n_cells), n_dry=Int(n_dry), n_wet=Int(n_wet), dry_share=_ratio(n_dry, n_cells),
            n_days=length(used),
            mean_dry=point[1], mean_dry_lo=lo(1), mean_dry_hi=hi(1),
            RMSE_dry=point[2], RMSE_dry_lo=lo(2), RMSE_dry_hi=hi(2),
            mean_obs_dry=_ratio(total[VOL_OBS_DRY], n_dry),
            zero_share_dry=_ratio(total[ZERO_DRY], n_dry), dry_volume_share=_ratio(total[VOL_DRY], total[VOL_ALL]),
            dry_sse_share=point[3],
            spurious_mm_per_year=_ratio(total[VOL_DRY], n_cells) * HOURS_PER_YEAR,
            estimate_mm_per_year=_ratio(total[VOL_ALL], n_cells) * HOURS_PER_YEAR,
            gauge_mm_per_year=_ratio(total[VOL_OBS], n_cells) * HOURS_PER_YEAR,
        )
        for (k, threshold) in enumerate(thresholds)
            offset = 3 + 6 * (k - 1)
            scores = Pair{Symbol,Any}[:threshold => Float64(threshold), :POD => _ratio(total[BASE + 2k - 1], n_wet),
                :specificity => 1 - point[offset + 1]]
            for (s, name) in enumerate(_THRESHOLD_SCORES)
                append!(scores, [name => point[offset + s], Symbol(name, "_lo") => lo(offset + s),
                    Symbol(name, "_hi") => hi(offset + s)])
            end
            push!(rows, merge(key, base, NamedTuple{Tuple(first.(scores))}(Tuple(last.(scores)))))
        end
    end
    return DataFrame(rows)
end

"""
Occurrence and dry-hour skill of an estimate against the gauges, per level of a stratifier.

`labels` (a matrix, or a function `(i, j) -> label`) puts cell `(i, j)` in level 1..length(levels); 0,
or a cell outside `mask` or with a non-finite value, leaves it out. A cell is gauge-wet at
`obs >= wet`, otherwise gauge-dry; the estimate is wet at `est >= threshold` for each of `thresholds`.
One row per level and threshold, prefixed by the fields of `key`:

- counts: `n_cells`, `n_dry`, `n_wet`, `dry_share`, `n_days`
- occurrence at the threshold: `POFD` = P(est wet | gauge dry), the false-alarm *rate*;
  `specificity` = 1 - POFD; `FAR` = P(gauge dry | est wet), the false-alarm *ratio*;
  `NPV` = P(gauge dry | est dry); `POD`; `freq_bias` = est-wet over gauge-wet cells; `PSS` = POD - POFD;
  `HSS`. Wet-side scores are `NaN` for a level without gauge-wet cells.
- dry-hour amounts, the same on every threshold row: `mean_dry` (what the estimate reports on gauge-dry
  cells, mm/h), `mean_obs_dry` (the gauge's own mean there: a few readings are non-zero but below 0.1),
  `RMSE_dry`, `zero_share_dry` (estimate exactly 0), `dry_volume_share` (estimate volume on gauge-dry
  cells over its total), `dry_sse_share` (squared error on gauge-dry cells over the level's total), and
  per year of record `spurious_mm_per_year` (estimate volume on gauge-dry cells), `estimate_mm_per_year`
  and `gauge_mm_per_year`.

`_lo` / `_hi` are 95% day-block bootstrap bounds; `reps = 0` leaves them `NaN`.
"""
function occurrence_table(
    labels, levels::AbstractVector, Y_obs::AbstractMatrix, Y_est::AbstractMatrix, mask::AbstractMatrix,
    day::AbstractVector{<:Integer}; thresholds=DRY_THRESHOLDS, wet::Real=WET_MM, reps::Int=1000,
    seed::Int=20260919, key::NamedTuple=(;),
)
    size(Y_obs) == size(Y_est) == size(mask) || throw(DimensionMismatch("obs, estimate and mask must match"))
    length(day) == size(Y_obs, 2) || throw(DimensionMismatch("one day per column"))
    labels isa AbstractMatrix && size(labels) != size(Y_obs) && throw(DimensionMismatch("labels must match obs"))
    T = length(thresholds)
    acc = zeros(length(levels), maximum(day), BASE + 2T)
    for j in axes(Y_obs, 2)
        d = day[j]
        for i in axes(Y_obs, 1)
            mask[i, j] || continue
            l = _label(labels, i, j)
            l == 0 && continue
            o, e = Y_obs[i, j], Y_est[i, j]
            (isfinite(o) && isfinite(e)) || continue
            err2 = (e - o)^2
            acc[l, d, VOL_ALL] += e
            acc[l, d, SSE_ALL] += err2
            acc[l, d, VOL_OBS] += o
            if o >= wet
                acc[l, d, N_WET] += 1
                for k in 1:T
                    e >= thresholds[k] && (acc[l, d, BASE + 2k - 1] += 1)
                end
            else
                acc[l, d, N_DRY] += 1
                acc[l, d, VOL_DRY] += e
                acc[l, d, VOL_OBS_DRY] += o
                acc[l, d, ZERO_DRY] += e == 0
                acc[l, d, SSE_DRY] += err2
                for k in 1:T
                    e >= thresholds[k] && (acc[l, d, BASE + 2k] += 1)
                end
            end
        end
    end
    return _occurrence_rows(acc, levels, thresholds; reps, seed, key)
end

"""
Squared error of `Y_method` against `Y_ref` on the same cells, per level of a stratifier.

Per level: `n`, `RMSE_method`, `RMSE_ref`, `improvement` = 1 - RMSE_method / RMSE_ref (positive: the
method is better), the mean of each (`mean_method`, `mean_ref`, `mean_obs`), and `sse_gap` = SSE_method
- SSE_ref with `sse_gap_share`, its share of the summed gap over all levels - meaningful when the levels
partition the cells, as the gauge x anchor quadrants do. `_lo` / `_hi`: paired day-block bootstrap
bounds of `improvement` and of `mean_diff` = mean_method - mean_ref.
"""
function paired_error_table(
    labels, levels::AbstractVector, Y_obs::AbstractMatrix, Y_method::AbstractMatrix, Y_ref::AbstractMatrix,
    mask::AbstractMatrix, day::AbstractVector{<:Integer}; reps::Int=1000, seed::Int=20260919,
    key::NamedTuple=(;),
)
    size(Y_obs) == size(Y_method) == size(Y_ref) == size(mask) || throw(DimensionMismatch("matrices must match"))
    acc = zeros(length(levels), maximum(day), 6)     # n, SSE method, SSE ref, sum method, sum ref, sum obs
    for j in axes(Y_obs, 2)
        d = day[j]
        for i in axes(Y_obs, 1)
            mask[i, j] || continue
            l = _label(labels, i, j)
            l == 0 && continue
            o, m, r = Y_obs[i, j], Y_method[i, j], Y_ref[i, j]
            (isfinite(o) && isfinite(m) && isfinite(r)) || continue
            acc[l, d, 1] += 1
            acc[l, d, 2] += (m - o)^2
            acc[l, d, 3] += (r - o)^2
            acc[l, d, 4] += m
            acc[l, d, 5] += r
            acc[l, d, 6] += o
        end
    end
    scores(t) = (1 - sqrt(_ratio(t[2], t[1])) / sqrt(_ratio(t[3], t[1])), _ratio(t[4] - t[5], t[1]))
    gaps = [sum(acc[l, :, 2]) - sum(acc[l, :, 3]) for l in eachindex(levels)]
    rows = NamedTuple[]
    for (l, level) in enumerate(levels)
        per_day = acc[l, :, :]
        used = [d for d in axes(per_day, 1) if per_day[d, 1] > 0]
        t = vec(sum(per_day; dims=1))
        point = scores(t)
        ci = reps > 0 && !isempty(used) ?
            day_block_bootstrap(collect(eachindex(used)), per_day[used, :], scores; reps, rng=MersenneTwister(seed + l)) :
            nothing
        push!(rows, merge(key, (;
            level, n=Int(t[1]), RMSE_method=sqrt(_ratio(t[2], t[1])), RMSE_ref=sqrt(_ratio(t[3], t[1])),
            improvement=point[1], improvement_lo=ci === nothing ? NaN : ci.lo[1],
            improvement_hi=ci === nothing ? NaN : ci.hi[1],
            mean_method=_ratio(t[4], t[1]), mean_ref=_ratio(t[5], t[1]), mean_obs=_ratio(t[6], t[1]),
            mean_diff=point[2], mean_diff_lo=ci === nothing ? NaN : ci.lo[2], mean_diff_hi=ci === nothing ? NaN : ci.hi[2],
            sse_gap=gaps[l], sse_gap_share=sum(gaps) != 0 ? gaps[l] / sum(gaps) : NaN,
        )))
    end
    return DataFrame(rows)
end

"""
The false-alarm rate of a *gauge* as an estimate for its neighbour, by separation.

Every ordered pair of stations `(i, k)` whose separation falls in `[edges_km[b], edges_km[b + 1])` scores
gauge `k` as the estimate at station `i`, on the hours both report inside `mask`; a last level pairs every
station with its nearest neighbour whatever the distance. Scored with `occurrence_table`'s columns at
`threshold = wet`, plus `n_pairs` and `median_km`. Two gauges a few km apart sample the same rain field
at two points, so this is the floor a point-to-area comparison cannot get under: part of what a
satellite pixel is charged as a false alarm is rain the gauge simply was not under.
"""
function neighbour_gauge_baseline(
    Y_obs::AbstractMatrix, lonlat::AbstractMatrix, mask::AbstractMatrix, day::AbstractVector{<:Integer};
    edges_km=[0.0, 5.0, 10.0, 15.0, 20.0, 30.0], wet::Real=WET_MM, reps::Int=1000, seed::Int=20260919,
    key::NamedTuple=(;),
)
    n = size(Y_obs, 1)
    distance = haversine_distance_matrix(lonlat, lonlat)
    for i in 1:n
        distance[i, i] = Inf
    end
    B = length(edges_km) - 1
    levels = vcat(["$(Int(edges_km[b]))-$(Int(edges_km[b + 1])) km" for b in 1:B], ["nearest gauge"])
    pairs = [Tuple{Int,Int}[] for _ in 1:(B + 1)]
    for i in 1:n, k in 1:n
        b = searchsortedlast(edges_km, distance[k, i])
        1 <= b <= B && push!(pairs[b], (i, k))
    end
    for i in 1:n
        push!(pairs[B + 1], (i, argmin(view(distance, :, i))))
    end
    acc = zeros(length(levels), maximum(day), BASE + 2)
    for (l, level_pairs) in enumerate(pairs), (i, k) in level_pairs, j in axes(Y_obs, 2)
        (mask[i, j] && mask[k, j]) || continue
        o, e = Y_obs[i, j], Y_obs[k, j]
        (isfinite(o) && isfinite(e)) || continue
        d = day[j]
        acc[l, d, VOL_ALL] += e
        acc[l, d, SSE_ALL] += (e - o)^2
        acc[l, d, VOL_OBS] += o
        if o >= wet
            acc[l, d, N_WET] += 1
            e >= wet && (acc[l, d, BASE + 1] += 1)
        else
            acc[l, d, N_DRY] += 1
            acc[l, d, VOL_DRY] += e
            acc[l, d, VOL_OBS_DRY] += o
            acc[l, d, ZERO_DRY] += e == 0
            acc[l, d, SSE_DRY] += (e - o)^2
            e >= wet && (acc[l, d, BASE + 2] += 1)
        end
    end
    table = _occurrence_rows(acc, levels, [wet]; reps, seed, key)
    table.n_pairs = length.(pairs)
    table.median_km = [isempty(p) ? NaN : median(distance[k, i] for (i, k) in p) for p in pairs]
    return table
end

# ---------------------------------------------------------------------------------------------
# Dry spells and dry days

"""
Runs of consecutive `true` in `flags` over the positions where `valid` holds.

`breaks[t]` marks a gap in time before position `t` (a non-contiguous axis), which ends a run as a
missing position does. Returns `(; lengths, censored)`: a run is censored when it touches the ends of
the axis, a missing position or a break, since its true length is then unknown.
"""
function run_lengths(flags::AbstractVector{Bool}, valid::AbstractVector{Bool};
        breaks::AbstractVector{Bool}=falses(length(flags)))
    length(flags) == length(valid) == length(breaks) || throw(DimensionMismatch("flags, valid and breaks must match"))
    lengths, censored = Int[], Bool[]
    n = length(flags)
    t = 1
    while t <= n
        if valid[t] && flags[t]
            start = t
            while t < n && valid[t + 1] && flags[t + 1] && !breaks[t + 1]
                t += 1
            end
            left_open = start == 1 || breaks[start] || !valid[start - 1]
            right_open = t == n || breaks[t + 1] || !valid[t + 1]
            push!(lengths, t - start + 1)
            push!(censored, left_open || right_open)
        end
        t += 1
    end
    return (; lengths, censored)
end

"""
Dry spells of one source: runs of cells below `threshold` inside `valid`.

Returns `(; summary, histogram)`. `summary` has a pooled `all` row and one row per station: dry share
of the valid cells, the number of complete and censored spells, and the mean, median, 90th percentile
and maximum length of the complete ones, plus the share of all valid cells lying in complete spells of
at least 24 and 72 steps. `histogram` bins the pooled complete spells by `bin_edges` (lower bounds):
spell count and the cells they hold. On a daily axis the same numbers are in days.
"""
function dry_spell_table(
    Y::AbstractMatrix, valid::AbstractMatrix; threshold::Real, station_ids::AbstractVector,
    breaks::AbstractVector{Bool}=falses(size(Y, 2)), bin_edges=SPELL_BIN_EDGES, key::NamedTuple=(;),
)
    size(Y) == size(valid) || throw(DimensionMismatch("Y and valid must match"))
    summarize(level, lengths, censored, n_valid, n_dry) = begin
        complete = lengths[.!censored]
        (; level, n_valid, dry_share=_ratio(n_dry, n_valid), n_spells=length(complete), n_censored=count(censored),
            mean_length=isempty(complete) ? NaN : mean(complete), median_length=isempty(complete) ? NaN : median(complete),
            p90_length=isempty(complete) ? NaN : quantile(complete, 0.9), max_length=isempty(complete) ? 0 : maximum(complete),
            share_in_ge24=_ratio(sum(filter(>=(24), complete); init=0), n_valid),
            share_in_ge72=_ratio(sum(filter(>=(72), complete); init=0), n_valid))
    end
    rows = NamedTuple[]
    pooled_lengths, pooled_censored = Int[], Bool[]
    pooled_valid, pooled_dry = 0, 0
    for i in axes(Y, 1)
        v = BitVector(view(valid, i, :) .& isfinite.(view(Y, i, :)))
        flags = BitVector([v[j] && Y[i, j] < threshold for j in axes(Y, 2)])
        runs = run_lengths(flags, v; breaks)
        n_valid, n_dry = count(v), count(flags)
        push!(rows, merge(key, summarize(String(station_ids[i]), runs.lengths, runs.censored, n_valid, n_dry)))
        append!(pooled_lengths, runs.lengths)
        append!(pooled_censored, runs.censored)
        pooled_valid += n_valid
        pooled_dry += n_dry
    end
    pushfirst!(rows, merge(key, summarize("all", pooled_lengths, pooled_censored, pooled_valid, pooled_dry)))
    complete = pooled_lengths[.!pooled_censored]
    hist = NamedTuple[]
    for (b, lower) in enumerate(bin_edges)
        upper = b < length(bin_edges) ? bin_edges[b + 1] : typemax(Int)
        in_bin = filter(x -> lower <= x < upper, complete)
        push!(hist, merge(key, (; bin_lower=lower, bin_upper=upper == typemax(Int) ? -1 : upper,
            n_spells=length(in_bin), cells_in_spells=sum(in_bin; init=0),
            spell_share=_ratio(length(in_bin), length(complete)),
            cell_share=_ratio(sum(in_bin; init=0), sum(complete; init=0)))))
    end
    return (; summary=DataFrame(rows), histogram=DataFrame(hist))
end

"""
Daily sums over the `valid` cells of each met-day: `totals[station, day]` for `day = 1:maximum(day)`,
`NaN` where fewer than `min_hours` of that station-day's hours are valid. Summing gauge and estimate over
one shared `valid` mask keeps a partial day fair: both sides then cover the same hours.
"""
function daily_sums(day::AbstractVector{<:Integer}, Y::AbstractMatrix, valid::AbstractMatrix; min_hours::Int=24)
    size(Y) == size(valid) || throw(DimensionMismatch("Y and valid must match"))
    length(day) == size(Y, 2) || throw(DimensionMismatch("one day per column"))
    D = maximum(day)
    totals, hours = zeros(size(Y, 1), D), zeros(Int, size(Y, 1), D)
    for j in axes(Y, 2), i in axes(Y, 1)
        (valid[i, j] && isfinite(Y[i, j])) || continue
        totals[i, day[j]] += Y[i, j]
        hours[i, day[j]] += 1
    end
    totals[hours .< min_hours] .= NaN
    return totals
end

"""
Runs of dry days at a gauge during which its neighbours were clearly wet: the signature of a gauge that
reads 0.0 while it is not measuring (a blocked funnel, a stopped logger), not of a dry spell.

`D[station, day]` holds daily totals (`NaN` when incomplete). Every run of at least `min_days` consecutive
complete days below `dry_mm` is compared, day by day, with the median of the `k` nearest stations reporting
that day; the run is flagged when that median reaches `wet_mm` on at least `min_wet_days` of its days.
`p_dry_given_wet` - how often any gauge stays below `dry_mm` on a day its neighbours' median reaches
`wet_mm` - is returned alongside, as the single-day rate the rule rests on: several such days inside one
dry run are very unlikely to be weather.

Returns `(; spells, flagged, p_dry_given_wet)`: one row per run (`station_id`, `first_day`, `last_day`
as indices into the day axis, `n_days`, `gauge_total`, `neighbour_total`, `neighbour_wet_days`,
`flagged`), and `flagged[station, day]`, true on every day of a flagged run.
"""
function silent_gauge_spells(D::AbstractMatrix, lonlat::AbstractMatrix, station_ids::AbstractVector;
        dry_mm::Real=1.0, wet_mm::Real=5.0, min_days::Int=7, min_wet_days::Int=3, k::Int=3)
    n, n_days = size(D)
    size(lonlat, 1) == n == length(station_ids) || throw(DimensionMismatch("one lon/lat row and id per station"))
    distance = haversine_distance_matrix(lonlat, lonlat)
    order = [filter(!=(i), sortperm(view(distance, :, i))) for i in 1:n]
    neighbour = fill(NaN, n, n_days)
    for d in 1:n_days, i in 1:n
        values = Float64[]
        for m in order[i]
            isfinite(D[m, d]) && push!(values, D[m, d])
            length(values) == k && break
        end
        isempty(values) || (neighbour[i, d] = median(values))
    end
    paired = [(D[c] < dry_mm) for c in eachindex(D) if isfinite(D[c]) && isfinite(neighbour[c]) && neighbour[c] >= wet_mm]
    flagged = falses(n, n_days)
    rows = NamedTuple[]
    for i in 1:n
        d = 1
        while d <= n_days
            if isfinite(D[i, d]) && D[i, d] < dry_mm
                first_day = d
                while d < n_days && isfinite(D[i, d + 1]) && D[i, d + 1] < dry_mm
                    d += 1
                end
                if d - first_day + 1 >= min_days
                    around = filter(isfinite, neighbour[i, first_day:d])
                    wet_days = count(>=(wet_mm), around)
                    flag = wet_days >= min_wet_days
                    flag && (flagged[i, first_day:d] .= true)
                    push!(rows, (; station_id=String(station_ids[i]), first_day, last_day=d, n_days=d - first_day + 1,
                        gauge_total=sum(D[i, first_day:d]), neighbour_total=sum(around; init=0.0),
                        neighbour_wet_days=wet_days, flagged=flag))
                end
            end
            d += 1
        end
    end
    return (; spells=DataFrame(rows), flagged, p_dry_given_wet=isempty(paired) ? NaN : mean(paired))
end

# ---------------------------------------------------------------------------------------------
# Inputs and sensitivity

"""
A long table (`time`, `station`, `value`) as a `[station, time]` matrix on `ref_ids` x `ref_times`.

`time + Hour(offset_hours)` is matched to `ref_times`: 9 turns a UTC hour-start label into the
hour-ending Beijing-time label the gauges use. Rows for other stations or hours are skipped; a
duplicate cell is an error. Returns `(; Y, matched)`.
"""
function long_to_wide_hourly(
    time::AbstractVector{DateTime}, station::AbstractVector{<:AbstractString}, value::AbstractVector{<:Real},
    ref_times::AbstractVector{DateTime}, ref_ids::AbstractVector{<:AbstractString}; offset_hours::Int=9,
)
    length(time) == length(station) == length(value) || throw(DimensionMismatch("time, station and value must match"))
    row_of = Dict(String(id) => i for (i, id) in enumerate(ref_ids))
    column_of = Dict(t => j for (j, t) in enumerate(ref_times))
    Y = fill(NaN, length(ref_ids), length(ref_times))
    seen = falses(size(Y))
    matched = 0
    for k in eachindex(time)
        i = get(row_of, String(station[k]), 0)
        j = get(column_of, time[k] + Hour(offset_hours), 0)
        (i == 0 || j == 0) && continue
        seen[i, j] && throw(ArgumentError("duplicate value for $(station[k]) at $(time[k])"))
        seen[i, j] = true
        Y[i, j] = value[k]
        matched += 1
    end
    return (; Y, matched)
end

"""
Pooled Pearson r between `A[:, j]` and `B[:, j + lag]` for each lag, over cells both finite and inside
`mask` at their own column. A clock offset between two series of the same field shows as a peak away
from lag 0.
"""
function lagged_correlation(A::AbstractMatrix, B::AbstractMatrix, mask::AbstractMatrix; lags=-3:3)
    size(A) == size(B) == size(mask) || throw(DimensionMismatch("A, B and mask must match"))
    rows = NamedTuple[]
    for lag in lags
        sx = sy = sxx = syy = sxy = 0.0
        n = 0
        for j in axes(A, 2)
            jj = j + lag
            1 <= jj <= size(B, 2) || continue
            for i in axes(A, 1)
                (mask[i, j] && mask[i, jj]) || continue
                x, y = A[i, j], B[i, jj]
                (isfinite(x) && isfinite(y)) || continue
                n += 1
                sx += x; sy += y; sxx += x^2; syy += y^2; sxy += x * y
            end
        end
        cov = sxy / n - (sx / n) * (sy / n)
        push!(rows, (; lag, n, r=cov / sqrt((sxx / n - (sx / n)^2) * (syy / n - (sy / n)^2))))
    end
    return DataFrame(rows)
end

"""
What setting every value below `tau` to zero would change, for each `tau` in `taus` (0 is the estimate
as it is): overall and dry-hour RMSE, the mean on gauge-dry cells, the volume ratio to the gauges, and
POD / POFD / FAR / CSI at `wet`. Descriptive only - a `tau` read off this table is chosen on held-out
cells, so it is not an out-of-sample result.
"""
function zero_threshold_sensitivity(Y_obs::AbstractMatrix, Y_est::AbstractMatrix, mask::AbstractMatrix;
        taus=[0.0, 0.05, 0.1, 0.2, 0.3, 0.5], wet::Real=WET_MM, key::NamedTuple=(;))
    size(Y_obs) == size(Y_est) == size(mask) || throw(DimensionMismatch("obs, estimate and mask must match"))
    rows = NamedTuple[]
    for tau in taus
        n = n_dry = hits = misses = false_alarms = 0
        sse = sse_dry = vol = vol_dry = vol_obs = 0.0
        for k in eachindex(Y_obs)
            mask[k] || continue
            o, e = Y_obs[k], Y_est[k]
            (isfinite(o) && isfinite(e)) || continue
            e = e < tau ? 0.0 : e
            n += 1
            sse += (e - o)^2
            vol += e
            vol_obs += o
            if o >= wet
                e >= wet ? (hits += 1) : (misses += 1)
            else
                n_dry += 1
                sse_dry += (e - o)^2
                vol_dry += e
                e >= wet && (false_alarms += 1)
            end
        end
        push!(rows, merge(key, (; tau=Float64(tau), n, RMSE=sqrt(_ratio(sse, n)), RMSE_dry=sqrt(_ratio(sse_dry, n_dry)),
            mean_dry=_ratio(vol_dry, n_dry), volume_ratio=_ratio(vol, vol_obs), POD=_ratio(hits, hits + misses),
            POFD=_ratio(false_alarms, n_dry), FAR=_ratio(false_alarms, hits + false_alarms),
            CSI=_ratio(hits, hits + misses + false_alarms))))
    end
    return DataFrame(rows)
end

end
