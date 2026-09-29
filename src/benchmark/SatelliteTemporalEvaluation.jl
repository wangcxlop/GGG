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
- The gauge record is quantised to 0.5 mm: its "light" hours are 0.5, 1.0 or 1.5 mm and it cannot
  show 0.1-0.5 mm. Where a satellite threshold matters, a `_res` variant repeats the calculation at the
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

include("satellite_temporal/intensity_scores.jl")
include("satellite_temporal/rain_events.jl")
include("satellite_temporal/diurnal_cycle.jl")
include("satellite_temporal/correlation_series.jl")

end # module
