function _empty_metric_row(; scheme, product, method, fold=missing, repeat=1, seed=0,
    group="overall", level="all",
    threshold=NaN, n=0, coverage=0.0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN,
    POD=NaN, FAR=NaN, CSI=NaN)
    return (;
        scheme=String(scheme), product=String(product), method=String(method), fold,
        repeat=Int(repeat), seed=Int(seed),
        group=String(group), level=String(level), threshold=Float64(threshold), n=Int(n),
        coverage=Float64(coverage), RMSE=Float64(RMSE), MAE=Float64(MAE),
        Bias=Float64(Bias), r=Float64(r), POD=Float64(POD), FAR=Float64(FAR), CSI=Float64(CSI),
    )
end

function _continuous_row(y_obs, prediction, mask; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs...)
    metric = metric_continuous(y_obs, prediction; mask=mask)
    return _empty_metric_row(;
        kwargs..., n=metric.n, coverage=metric.n / length(mask), RMSE=metric.RMSE,
        MAE=metric.MAE, Bias=metric.Bias, r=metric.r,
    )
end

function _event_row(y_obs, prediction, mask, threshold; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs..., threshold)
    metric = metric_event(y_obs, prediction; mask=mask, thr=threshold)
    return _empty_metric_row(;
        kwargs..., threshold, n, coverage=n / length(mask),
        POD=metric.POD, FAR=metric.FAR, CSI=metric.CSI,
    )
end

function append_stratified_metrics!(
    rows::Vector{NamedTuple}, scheme::String, product::String, method::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, prediction::Matrix{Float64},
    common_mask::BitMatrix, nearest_distance::Vector{Float64}, thresholds::Vector{Float64};
    fold=missing, repeat::Int=1, seed::Int=0,
)
    base = (; scheme, product, method, fold, repeat, seed)
    # `common_mask` is pinned to MASK_METHODS, so a method outside that set can still be NaN
    # inside it and would otherwise score NaN rather than "worse". Drop its own gaps and let
    # the `coverage` column report the shortfall instead. This is a no-op for the MASK_METHODS
    # themselves, whose NaNs already defined the mask.
    scored_mask = common_mask .& .!isnan.(prediction)
    push!(rows, _continuous_row(y_obs, prediction, scored_mask; base..., group="overall", level="all"))

    rain_groups = (
        ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
        ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
    )
    for (name, lower, upper) in rain_groups
        stratum = scored_mask .& (y_obs .>= lower) .& (y_obs .< upper)
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="rain_intensity", level=name))
    end

    for year_value in sort(unique(year.(times)))
        time_mask = reshape(year.(times) .== year_value, 1, :)
        stratum = scored_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="year", level=string(year_value)))
    end
    for month_value in sort(unique(month.(times)))
        time_mask = reshape(month.(times) .== month_value, 1, :)
        stratum = scored_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="month", level=lpad(month_value, 2, '0')))
    end

    distance_groups = (
        ("0_20", -Inf, 20.0), ("20_50", 20.0, 50.0),
        ("50_100", 50.0, 100.0), ("100_plus", 100.0, Inf),
    )
    for (name, lower, upper) in distance_groups
        station_mask = reshape((nearest_distance .>= lower) .& (nearest_distance .< upper), :, 1)
        stratum = scored_mask .& station_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="nearest_train_km", level=name))
    end

    for threshold in thresholds
        push!(rows, _event_row(y_obs, prediction, scored_mask, threshold;
            base..., group="event_threshold", level=string(threshold)))
    end
    return rows
end

function _common_method_mask(y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}})
    mask = .!isnan.(y_obs)
    for method in MASK_METHODS
        mask .&= .!isnan.(predictions[method])
    end
    return BitMatrix(mask)
end

const REPEAT_SUMMARY_METRICS = (:RMSE, :MAE, :Bias, :r, :POD, :FAR, :CSI)

_finite(values) = filter(isfinite, collect(skipmissing(values)))

"""
Spread of each pooled metric across repeated cross-validation partitions.

One row per (scheme, product, method, group, level) carrying mean/std/min/max over repeats, so a
reported value is read as a distribution rather than a single split's point estimate. `mean_n` and
`mean_coverage` are kept because `_common_method_mask` requires every method to be finite, so the
evaluated sample size itself shifts between partitions.
"""
function summarize_repeats(metrics::DataFrame)
    nrow(metrics) == 0 && return DataFrame()
    pooled = filter(row -> ismissing(row.fold), metrics)
    nrow(pooled) == 0 && return DataFrame()
    rows = NamedTuple[]
    for subset in groupby(pooled, [:scheme, :product, :method, :group, :level])
        entry = (;
            scheme=first(subset.scheme), product=first(subset.product),
            method=first(subset.method), group=first(subset.group), level=first(subset.level),
            n_repeats=nrow(subset), mean_n=mean(subset.n), mean_coverage=mean(subset.coverage),
        )
        for name in REPEAT_SUMMARY_METRICS
            values = _finite(subset[!, name])
            entry = merge(entry, (;
                Symbol(name, :_mean) => isempty(values) ? NaN : mean(values),
                Symbol(name, :_std) => length(values) < 2 ? NaN : std(values),
                Symbol(name, :_min) => isempty(values) ? NaN : minimum(values),
                Symbol(name, :_max) => isempty(values) ? NaN : maximum(values),
            ))
        end
        push!(rows, entry)
    end
    return sort!(DataFrame(rows), [:scheme, :product, :group, :level, :method])
end

"""
Across-fold spread of each method's per-fold metrics, within a single partition.

`metrics_pooled.csv` stitches every fold's out-of-fold predictions into one matrix and scores that
once. It is the right point estimate, but it carries no dispersion at all: a method that is steady
across the region and one that is excellent in three folds and poor in two report the same number.
This summarises `metrics_folds.csv` so the two can be read side by side.

Read `_std` as a **dispersion diagnostic, not a standard error**. The folds share training stations
and the same satellite record, so they are not independent replicates and `_std / sqrt(n_folds)` is
not a valid uncertainty for the pooled estimate — `paired_comparisons.csv` is where uncertainty is
quantified. What the spread does say is whether a pooled gap between two methods is large or small
next to how much either one moves between folds, which on this data it very often is not.

Grouped by `repeat` as well, so repeated partitions summarise separately rather than having their
fold-to-fold and partition-to-partition variation silently pooled into one number;
`summarize_repeats` is the across-partition view.
"""
function summarize_fold_spread(metrics::DataFrame)
    nrow(metrics) == 0 && return DataFrame()
    per_fold = filter(row -> !ismissing(row.fold), metrics)
    nrow(per_fold) == 0 && return DataFrame()
    rows = NamedTuple[]
    for subset in groupby(per_fold, [:scheme, :product, :method, :group, :level, :repeat])
        entry = (;
            scheme=first(subset.scheme), product=first(subset.product),
            method=first(subset.method), group=first(subset.group), level=first(subset.level),
            repeat=first(subset.repeat), n_folds=nrow(subset), total_n=sum(subset.n),
        )
        for name in REPEAT_SUMMARY_METRICS
            values = _finite(subset[!, name])
            entry = merge(entry, (;
                Symbol(name, :_mean) => isempty(values) ? NaN : mean(values),
                Symbol(name, :_std) => length(values) < 2 ? NaN : std(values),
                Symbol(name, :_min) => isempty(values) ? NaN : minimum(values),
                Symbol(name, :_max) => isempty(values) ? NaN : maximum(values),
            ))
        end
        push!(rows, entry)
    end
    return sort!(DataFrame(rows), [:scheme, :product, :repeat, :group, :level, :method])
end

"""
How often each method wins, across repeated cross-validation partitions.

Methods are ranked by pooled RMSE within every repeat; a `win_count` well below `n_repeats` means
the ranking is an artifact of the particular split rather than a stable result.
"""
function method_rank_stability(metrics::DataFrame)
    nrow(metrics) == 0 && return DataFrame()
    pooled = filter(row -> ismissing(row.fold) && isfinite(row.RMSE), metrics)
    nrow(pooled) == 0 && return DataFrame()
    rows = NamedTuple[]
    for cell in groupby(pooled, [:scheme, :product, :group, :level])
        ranks = Dict{String,Vector{Int}}()
        for repeat_value in sort(unique(cell.repeat))
            ordered = sort(filter(:repeat => ==(repeat_value), DataFrame(cell)), :RMSE)
            for (position, row) in enumerate(eachrow(ordered))
                push!(get!(ranks, row.method, Int[]), position)
            end
        end
        for method in sort(collect(keys(ranks)))
            positions = ranks[method]
            push!(rows, (;
                scheme=first(cell.scheme), product=first(cell.product),
                group=first(cell.group), level=first(cell.level), method,
                n_repeats=length(positions), mean_rank=mean(positions),
                best_rank=minimum(positions), worst_rank=maximum(positions),
                win_count=count(==(1), positions),
            ))
        end
    end
    return sort!(DataFrame(rows), [:scheme, :product, :group, :level, :mean_rank])
end

"""
Pooled balanced-spatial metric for one repeat.

`repeat` used to be pinned to `1` here and in `assess_gwr_claim`'s bootstrap filter, which meant a
repeated run could produce every partition's metrics and still report the claim from the first one
— the extra partitions were computed and then discarded at the last step. The repeat is now the
caller's choice, and `assess_gwr_claim` walks all of them.
"""
function _one_metric(
    metrics::DataFrame, product::String, method::String, group::String, level::String;
    repeat::Int=1,
)
    rows = filter(row ->
        row.scheme == "balanced_spatial" && row.product == product && row.method == method &&
        ismissing(row.fold) && row.repeat == repeat && row.group == group && row.level == level,
        metrics,
    )
    nrow(rows) == 1 || throw(ArgumentError(
        "expected one pooled metric for $product/$method/$group/$level at repeat $repeat, " *
        "got $(nrow(rows))",
    ))
    return rows[1, :]
end

"""
Assess whether the residual-GWR correction is supported for each product.

Every comparison against the traditional baselines (`idw`/`adw`/`tps`) is required to hold
against **all three**, not just whichever looks friendliest on this pooled sample. The previous
version picked a single baseline per stratum via `argmin`/`argmax` on the same data it reported
— an unprotected selection-on-test step — then tested or compared only against that one.
Requiring all three removes the selection entirely rather than correcting for it: the three
per-baseline `overall` comparisons are already Holm-corrected against each other by
`paired_bootstrap_rows`, and no further correction is needed to require all three to pass.
`best_traditional`/`RMSE_best_traditional`/`overall_relative_improvement` are kept as
descriptive-only fields (how much better than the strongest competitor) and do not gate
`product_supported`; the `*_win_count` fields show how many of the 3 baselines were actually
beaten in each stratum.

`method` selects which method is assessed (default `residual_gwr`, the historical behaviour).
`own_coverage` maps product to the method's *own* prediction coverage; when supplied, the
coverage gate uses it instead of `coverage` on the shared evaluation mask. The shared mask is
pinned to `MASK_METHODS`; historically `mgwr` back-fit non-convergence dropped whole hours and
dragged it down to ~0.8 (not `gwr`, which has always had ~1.0 own coverage), so gating on the
shared mask asked every method to answer for mgwr's failures, which no method could pass. That
non-convergence is now fixed and the shared mask sits at ~0.96-0.98 (see
scripts/run_mgwr_diagnostics.jl), but `own_coverage` remains available since it's still the more
correct number to gate on. Both numbers are reported side by side. Naming either argument
switches the output to a self-describing schema (`method` + `RMSE_method`) so rows for different
methods cannot be silently mixed.
"""
function assess_gwr_claim(
    metrics::DataFrame, bootstrap::DataFrame, products::Vector{String};
    method::String=DEFAULT_CLAIM_METHOD, own_coverage::Union{Nothing,AbstractDict}=nothing,
)
    all_rows = NamedTuple[]
    tagged = method != DEFAULT_CLAIM_METHOD || own_coverage !== nothing
    bootstrap_has_method = :method in propertynames(bootstrap)
    for repeat_value in claim_repeats(metrics)
    rows = NamedTuple[]
    for product in products
        # `baseline` rather than `method` as the comprehension variable: the kwarg of the same
        # name is read two lines down, and Julia scoping made the shadow harmless only by luck.
        overall_traditional = [
            _one_metric(metrics, product, baseline, "overall", "all"; repeat=repeat_value)
            for baseline in TRADITIONAL_METHODS]
        best_index = argmin([row.RMSE for row in overall_traditional])
        best_method = TRADITIONAL_METHODS[best_index]
        best_overall = overall_traditional[best_index]
        gwr_overall = _one_metric(metrics, product, method, "overall", "all"; repeat=repeat_value)

        # Overall: GWR must be significantly better than every baseline (each row already
        # Holm-corrected against the other two by paired_bootstrap_rows), not just whichever
        # one was easiest to beat.
        overall_bootstrap = filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.stratum == "overall" && row.repeat == repeat_value &&
            (!bootstrap_has_method || row.method == method),
            bootstrap,
        )
        significant_per_baseline = Dict(String(row.baseline) =>
            row.ci_low > 0 && row.pvalue_holm < 0.05 for row in eachrow(overall_bootstrap))
        overall_win_count = count(values(significant_per_baseline))
        significant = length(significant_per_baseline) == length(TRADITIONAL_METHODS) &&
            overall_win_count == length(TRADITIONAL_METHODS)

        # Heavy: GWR must improve on every baseline by at least 5%.
        gwr_heavy = _one_metric(metrics, product, method, "rain_intensity", "heavy"; repeat=repeat_value).RMSE
        heavy_improvements = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "rain_intensity", "heavy"; repeat=repeat_value).RMSE
                isfinite(base) && base > 0 ? (base - gwr_heavy) / base : NaN
            end for baseline in TRADITIONAL_METHODS)
        heavy_win_count = count(>=(0.05), filter(isfinite, collect(values(heavy_improvements))))
        heavy_ok = heavy_win_count == length(TRADITIONAL_METHODS)
        # Worst case across baselines, matching what the gate above requires.
        heavy_relative_improvement = minimum(values(heavy_improvements))

        # Moderate: GWR is allowed to be marginally worse here (non-inferiority), but not by
        # more than 2% against any baseline.
        gwr_moderate = _one_metric(metrics, product, method, "rain_intensity", "moderate"; repeat=repeat_value).RMSE
        moderate_degradations = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "rain_intensity", "moderate"; repeat=repeat_value).RMSE
                isfinite(base) && base > 0 ? (gwr_moderate - base) / base : NaN
            end for baseline in TRADITIONAL_METHODS)
        moderate_win_count = count(<=(0.02), filter(isfinite, collect(values(moderate_degradations))))
        moderate_ok = moderate_win_count == length(TRADITIONAL_METHODS)
        moderate_relative_degradation = maximum(values(moderate_degradations))

        # Year: a "win" requires beating every baseline that year, not just one.
        # `String.` because a `metrics` table read back from CSV carries InlineString levels,
        # which `_one_metric`'s `::String` signature would reject.
        year_levels = String.(unique(filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.method == method && ismissing(row.fold) && row.repeat == repeat_value &&
            row.group == "year",
            metrics,
        ).level))
        year_wins = 0
        for year_level in year_levels
            gwr_year = _one_metric(metrics, product, method, "year", year_level; repeat=repeat_value).RMSE
            year_wins += all(baseline ->
                gwr_year < _one_metric(metrics, product, baseline, "year", year_level; repeat=repeat_value).RMSE,
                TRADITIONAL_METHODS)
        end
        majority_years = isempty(year_levels) ? false : year_wins > length(year_levels) / 2

        # Event threshold: must not degrade CSI/FAR against any baseline.
        gwr_event = _one_metric(metrics, product, method, "event_threshold", "0.1"; repeat=repeat_value)
        event_checks = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "event_threshold", "0.1"; repeat=repeat_value)
                gwr_event.CSI >= base.CSI - 0.02 && gwr_event.FAR <= base.FAR + 0.02
            end for baseline in TRADITIONAL_METHODS)
        event_win_count = count(values(event_checks))
        event_not_degraded = event_win_count == length(TRADITIONAL_METHODS)

        gated_coverage = own_coverage === nothing ?
            gwr_overall.coverage : Float64(own_coverage[product])
        coverage_acceptable = gated_coverage >= 0.95
        direction_improved = gwr_overall.RMSE < best_overall.RMSE
        product_supported = significant && heavy_ok && moderate_ok && majority_years &&
            event_not_degraded && coverage_acceptable
        head = (;
            repeat=repeat_value, product, best_traditional=best_method,
            RMSE_best_traditional=best_overall.RMSE,
        )
        tail = (;
            overall_relative_improvement=(best_overall.RMSE - gwr_overall.RMSE) / best_overall.RMSE,
            paired_significant=significant, overall_win_count=overall_win_count,
            heavy_relative_improvement=heavy_relative_improvement, heavy_win_count=heavy_win_count,
            moderate_relative_degradation=moderate_relative_degradation,
            moderate_win_count=moderate_win_count,
            year_win_count=year_wins, year_count=length(year_levels), majority_years,
            event_not_degraded, event_win_count, common_coverage=gwr_overall.coverage,
            coverage_acceptable, direction_improved, product_supported,
        )
        push!(rows, if tagged
            merge((; method), head,
                (; RMSE_method=gwr_overall.RMSE, own_coverage=gated_coverage), tail)
        else
            merge(head, (; RMSE_residual_gwr=gwr_overall.RMSE), tail)
        end)
    end
    # Counted within the partition, not across partitions. Two products supported in different
    # repeats is not the same evidence as two products supported in the same one.
    supported_products = count(row -> row.product_supported, rows)
    improved_products = count(row -> row.direction_improved, rows)
    overall_supported = supported_products >= 2 && improved_products >= 2
    append!(all_rows, [merge(row, (;
        supported_product_count=supported_products,
        improved_product_count=improved_products,
        overall_claim_supported=overall_supported,
    )) for row in rows])
    end
    return DataFrame(all_rows)
end

"""
Partitions the claim can be assessed on: every repeat present in the pooled balanced-spatial rows.

A single-partition run yields `[1]`, which is what the claim path always used to assume.
"""
claim_repeats(metrics::DataFrame) = nrow(metrics) == 0 ? Int[] :
    sort(unique(Int.(filter(row ->
        row.scheme == "balanced_spatial" && ismissing(row.fold), metrics,
    ).repeat)))

"""
How consistently each product's claim survives across repeated partitions.

`claim_assessment.csv` now carries one row per (repeat, product); this collapses them so the
question "does this hold, or did it hold once" is answerable without reading the whole table. With
a single partition every `*_repeats` count is 1 and `unanimous` is trivially true — the table is
only informative once a run has more than one seed, which is the point of keeping it cheap.
"""
function claim_agreement(claim::DataFrame)
    (nrow(claim) == 0 || !(:repeat in propertynames(claim))) && return DataFrame()
    rows = NamedTuple[]
    for subset in groupby(claim, [:product])
        supported = count(subset.product_supported)
        push!(rows, (;
            product=first(subset.product), n_repeats=nrow(subset),
            supported_repeats=supported, improved_repeats=count(subset.direction_improved),
            unanimous=supported == nrow(subset) || supported == 0,
            mean_overall_relative_improvement=mean(subset.overall_relative_improvement),
            min_overall_relative_improvement=minimum(subset.overall_relative_improvement),
            max_overall_relative_improvement=maximum(subset.overall_relative_improvement),
        ))
    end
    return sort!(DataFrame(rows), :product)
end
