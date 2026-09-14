function _daily_bootstrap_delta(
    rng::AbstractRNG, times::Vector{DateTime}, y_obs::Matrix{Float64},
    baseline::Matrix{Float64}, gwr::Matrix{Float64}, mask::BitMatrix, reps::Int,
)
    days = Date.(times)
    unique_days = sort(unique(days))
    daily_n = zeros(Int, length(unique_days))
    daily_sse_baseline = zeros(Float64, length(unique_days))
    daily_sse_gwr = zeros(Float64, length(unique_days))
    for (day_index, day) in enumerate(unique_days)
        time_idx = findall(==(day), days)
        day_mask = mask[:, time_idx]
        daily_n[day_index] = count(day_mask)
        daily_sse_baseline[day_index] = sum(abs2, (baseline[:, time_idx] .- y_obs[:, time_idx])[day_mask])
        daily_sse_gwr[day_index] = sum(abs2, (gwr[:, time_idx] .- y_obs[:, time_idx])[day_mask])
    end
    keep = findall(>(0), daily_n)
    isempty(keep) && return Float64[]
    deltas = Vector{Float64}(undef, reps)
    for rep in 1:reps
        sampled = rand(rng, keep, length(keep))
        n = sum(daily_n[sampled])
        rmse_baseline = sqrt(sum(daily_sse_baseline[sampled]) / n)
        rmse_gwr = sqrt(sum(daily_sse_gwr[sampled]) / n)
        deltas[rep] = rmse_baseline - rmse_gwr
    end
    return deltas
end

"""
Day-block bootstrap of the difference in detection skill, `treatment - baseline`, at one threshold.

The categorical counterpart of `_daily_bootstrap_delta`. POD, FAR and CSI are ratios of pooled
contingency counts, so a day is resampled as its (hits, misses, false alarms) for both predictions
together: the pairing is kept, the hours and stations of one weather system stay in one block, and
every replicate's scores are recomputed from the summed counts exactly as `metric_event` computes
them from the cells (`>=` on both sides).

Only days with at least one masked cell are resampled, matching `_daily_bootstrap_delta`.

Returns `baseline` and `treatment` (POD, FAR, CSI on the full sample) and `deltas`, a `reps x 3`
matrix of `treatment - baseline` for POD, FAR and CSI. A replicate whose denominator is empty scores
NaN, so callers should summarise the finite values.
"""
function daily_bootstrap_event_delta(
    rng::AbstractRNG, times::Vector{DateTime}, y_obs::Matrix{Float64},
    baseline::Matrix{Float64}, treatment::Matrix{Float64}, mask::BitMatrix,
    threshold::Float64, reps::Int,
)
    size(baseline) == size(y_obs) == size(treatment) == size(mask) || throw(DimensionMismatch(
        "observations, predictions and mask must have the same shape",
    ))
    days = Date.(times)
    unique_days = sort(unique(days))
    # Columns: hits, misses, false alarms for the baseline, then the same for the treatment.
    counts = zeros(Int, length(unique_days), 6)
    cells = zeros(Int, length(unique_days))
    for (day_index, day) in enumerate(unique_days)
        for time in findall(==(day), days), station in axes(y_obs, 1)
            mask[station, time] || continue
            cells[day_index] += 1
            observed = y_obs[station, time] >= threshold
            for (offset, prediction) in ((0, baseline), (3, treatment))
                estimated = prediction[station, time] >= threshold
                if observed && estimated
                    counts[day_index, offset + 1] += 1
                elseif observed
                    counts[day_index, offset + 2] += 1
                elseif estimated
                    counts[day_index, offset + 3] += 1
                end
            end
        end
    end

    scores(hit, miss, false_alarm) = (
        hit + miss > 0 ? hit / (hit + miss) : NaN,
        hit + false_alarm > 0 ? false_alarm / (hit + false_alarm) : NaN,
        hit + miss + false_alarm > 0 ? hit / (hit + miss + false_alarm) : NaN,
    )
    function paired_scores(rows)
        total = vec(sum(@view(counts[rows, :]); dims=1))
        return scores(total[1], total[2], total[3]), scores(total[4], total[5], total[6])
    end

    keep = findall(>(0), cells)
    deltas = fill(NaN, reps, 3)
    isempty(keep) && return (; baseline=(NaN, NaN, NaN), treatment=(NaN, NaN, NaN), deltas)
    base_all, treat_all = paired_scores(keep)
    for rep in 1:reps
        base, treat = paired_scores(rand(rng, keep, length(keep)))
        for metric in 1:3
            deltas[rep, metric] = treat[metric] - base[metric]
        end
    end
    return (; baseline=base_all, treatment=treat_all, deltas)
end

function _holm_adjust(pvalues::AbstractVector)
    pvalues = Float64.(pvalues)
    m = length(pvalues)
    order = sortperm(pvalues)
    adjusted = fill(NaN, m)
    running = 0.0
    for (rank_index, original_index) in enumerate(order)
        value = min(1.0, (m - rank_index + 1) * pvalues[original_index])
        running = max(running, value)
        adjusted[original_index] = running
    end
    return adjusted
end

"""
Method under test when the caller does not name one. Kept as the default everywhere so the
benchmark's own `paired_comparisons.csv` / `claim_assessment.csv` are unchanged by the
method-aware parameters added for `scripts/run_claim_reassessment.jl`.
"""
const DEFAULT_CLAIM_METHOD = "residual_gwr"

function paired_bootstrap_rows(
    cfg::InterpolationBenchmarkConfig, scheme::String, product::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    common_mask::BitMatrix; repeat::Int=1, seed::Int=cfg.seed,
    method::String=DEFAULT_CLAIM_METHOD, pairwise_mask::Bool=false,
)
    rows = NamedTuple[]
    cfg.bootstrap_reps == 0 && return rows
    treatment = predictions[method]
    # A `method` column only appears once the caller asks for something other than the historical
    # single-method behaviour: the row already carries `baseline` but nothing naming the treatment,
    # so rows for two methods appended to one table would be indistinguishable.
    tagged = method != DEFAULT_CLAIM_METHOD || pairwise_mask
    strata = (
        ("overall", common_mask),
        ("moderate", common_mask .& (y_obs .>= 2.5) .& (y_obs .< 8.0)),
        ("heavy", common_mask .& (y_obs .>= 8.0)),
    )
    for (stratum, mask) in strata
        local_rows = NamedTuple[]
        for (baseline_index, baseline_method) in enumerate(TRADITIONAL_METHODS)
            baseline = predictions[baseline_method]
            # `common_mask` is pinned to MASK_METHODS, so a method outside that set can still be
            # NaN inside it and would score NaN rather than "worse". Restrict each test to the
            # cells both sides actually predicted.
            test_mask = pairwise_mask ?
                BitMatrix(mask .& .!isnan.(baseline) .& .!isnan.(treatment)) : BitMatrix(mask)
            deltas = _daily_bootstrap_delta(
                MersenneTwister(seed + 1000 * baseline_index + sum(codeunits(product))),
                times, y_obs, baseline, treatment, test_mask, cfg.bootstrap_reps,
            )
            isempty(deltas) && continue
            observed_baseline = metric_continuous(y_obs, baseline; mask=test_mask).RMSE
            observed_gwr = metric_continuous(y_obs, treatment; mask=test_mask).RMSE
            delta = observed_baseline - observed_gwr
            pvalue = min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0)))
            row = (;
                scheme, product, stratum, baseline=baseline_method, repeat, seed,
                n=count(test_mask), n_day=length(unique(Date.(times))), reps=cfg.bootstrap_reps,
                RMSE_baseline=observed_baseline, RMSE_gwr=observed_gwr,
                delta_RMSE=delta, relative_improvement=delta / observed_baseline,
                ci_low=quantile(deltas, 0.025), ci_high=quantile(deltas, 0.975),
                pvalue, pvalue_holm=NaN,
            )
            push!(local_rows, tagged ? merge((; method), row) : row)
        end
        adjusted = _holm_adjust([row.pvalue for row in local_rows])
        for (index, row) in enumerate(local_rows)
            push!(rows, merge(row, (; pvalue_holm=adjusted[index])))
        end
    end
    return rows
end
