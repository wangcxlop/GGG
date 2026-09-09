# `hurdle_gwr`: kept, but not part of the reported comparison.
#
# The model is absent from `BENCHMARK_METHODS`, `BENCHMARK_RUNS`, `MASK_METHODS` and
# `AUTO_CANDIDATE_RUNS` (see `InterpolationBenchmarkConfig.jl`), so nothing in a normal run
# reaches any of this. It is retained so re-enabling it stays the one-line change that comment
# promises, and collected here so a path that never executes stops padding out `predict_selected`
# and `select_interpolation_parameter!`, which are read constantly.
#
# Nothing below changed when it moved; only its address did.
#
# Be aware that none of this has automated coverage - the test suite cannot reach it, because
# `select_interpolation_parameter!` rejects any (mode, method) pair outside `BENCHMARK_RUNS`. To
# exercise it, push `("direct", "hurdle_gwr")` onto `BENCHMARK_RUNS` and call
# `build_hurdle_context` -> `select_interpolation_parameter!` -> `predict_selected`. That is how
# this move was checked: the scan rows, the selected candidate, the prediction and the diagnostic
# row all came out identical to the pre-move code.

"""
Per-fold context for `hurdle_gwr`.

`fit_hurdle_gwr`/`predict_hurdle_gwr` need the fold's `times` vector, which the rest of the
benchmark's method dispatch never has to carry; this mirrors the existing `dem_context` /
`joint_context` keyword pattern rather than widening every signature. `diagnostics` is a shared
accumulator (same pattern as `scan_rows`) for the degeneracy report described at
`_hurdle_predict`.
"""
function build_hurdle_context(
    times::Vector{DateTime}, cfg::InterpolationBenchmarkConfig,
    diagnostics::Vector{NamedTuple}; scheme::String, product::String, fold::Int,
    repeat::Int=1,
)
    wet_threshold = cfg.mger.rain_threshold
    # `HurdleGWRConfig` requires `first(intensity_breaks) == wet_threshold`, sorted and unique.
    # Deriving the breaks from the benchmark's own event thresholds keeps the occurrence table
    # binned the same way the categorical scores are.
    breaks = sort(unique(vcat(wet_threshold, filter(>(wet_threshold), cfg.event_thresholds))))
    return (; times, wet_threshold, intensity_breaks=breaks, diagnostics,
        scheme, product, fold, repeat)
end

"""Adapt the benchmark's `(adaptive, bw)` convention to HurdleGWR's `SpatialBandwidth`."""
_hurdle_bandwidth(adaptive::Bool, bw::Float64) =
    adaptive ? AdaptiveBandwidth(Int(round(bw))) : FixedBandwidth(bw)

_hurdle_config(context, kernel::Int, adaptive::Bool, bw::Float64) = HurdleGWRConfig(
    wet_threshold=context.wet_threshold,
    intensity_breaks=context.intensity_breaks,
    kernel=kernel,
    bandwidth=_hurdle_bandwidth(adaptive, bw),
)

"""
Fit the hurdle model on the training fold and predict `P(wet) * E[amount | wet]`.

Unlike every other method here this is a time-pooled calibration, not a per-hour spatial fit:
one local regression per target station, with time entering through annual/diurnal harmonics and
a calendar-month occurrence bin. The satellite enters the design as a *fitted* `log1p_satellite`
coefficient rather than as an offset with its coefficient forced to 1.

`predict_hurdle_gwr` takes no gauge observations at the target at all, so held-out leakage is
impossible by construction. `exclude_self` drives its `exclude_train_indices`, which zeroes a
station's own spatial weight - the same leave-one-out convention `_gwr_predict` uses.

Returns the full prediction object: the caller needs `used_global_amount` and
`occurrence_fallback` as well as `corrected`, because when a local fit misses
`min_amount_samples` / `min_effective_stations` the model *silently* falls back to global
coefficients, and an undiagnosed fallback would report a global model as a local one.
"""
function _hurdle_predict(
    train_lonlat::Matrix{Float64}, y_obs_train::Matrix{Float64},
    y_sat_train::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_sat_target::Matrix{Float64}, times::Vector{DateTime};
    config::HurdleGWRConfig, exclude_self::Bool=false,
)
    model = fit_hurdle_gwr(y_obs_train, y_sat_train, train_lonlat, times; config)
    exclusions = exclude_self ? collect(1:size(train_lonlat, 1)) : nothing
    return predict_hurdle_gwr(
        model, y_sat_target, target_lonlat, times; exclude_train_indices=exclusions,
    )
end

"""Summarise how much of a hurdle prediction actually came from a *local* fit."""
function _hurdle_diagnostic_row(prediction, context; kernel::Int, adaptive::Bool, bw::Float64)
    fallback = prediction.occurrence_fallback
    total = length(fallback)
    share(code) = total == 0 ? NaN : count(==(code), fallback) / total
    return (;
        scheme=context.scheme, product=context.product, fold=context.fold,
        repeat=context.repeat, kernel, adaptive, bw,
        n_target=length(prediction.used_global_amount),
        used_global_amount_fraction=mean(prediction.used_global_amount),
        mean_effective_stations=mean(prediction.effective_stations),
        occurrence_monthly=share(0x00), occurrence_all_month=share(0x01),
        occurrence_global=share(0x02), occurrence_skipped=share(0x03),
    )
end

"""Hurdle GWR's candidate scan: the shared kernel x bandwidth grid, minus the global
candidate, which `FixedBandwidth` cannot express."""
function _scan_hurdle_gwr!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig, row::NamedTuple, scored,
    train_lonlat::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64},
    time_weights, bandwidth_families, hurdle_context, tuning_times,
)
    for kernel in cfg.mger.kernels
        for (adaptive, candidates) in bandwidth_families
            for bw in candidates
                # `FixedBandwidth` rejects a non-finite bandwidth by construction, so the
                # family's global candidate cannot be expressed as a `HurdleGWRConfig`.
                # Skip it rather than logging a guaranteed-failure row per kernel; giving
                # `hurdle_gwr` a global option means teaching `SpatialBandwidth` about one.
                isfinite(bw) || continue
                _scan_candidate!(scan_rows; row..., kernel, adaptive, bw) do
                    hurdle_config = _hurdle_config(hurdle_context, kernel, adaptive, bw)
                    corrected = scored(
                        (tr, va) -> _hurdle_predict(
                            train_lonlat[tr, :], y_obs[tr, :], y_sat[tr, :],
                            train_lonlat[va, :], y_sat[va, :], something(tuning_times);
                            config=hurdle_config,
                        ).corrected,
                        () -> _hurdle_predict(
                            train_lonlat, y_obs, y_sat, train_lonlat, y_sat,
                            something(tuning_times);
                            config=hurdle_config, exclude_self=true,
                        ).corrected,
                    )
                    _candidate_metrics(y_obs, y_sat, corrected; time_weights)
                end
            end
        end
    end
    return scan_rows
end

"""
Predict with the hyperparameters `select_interpolation_parameter!` chose, and record how local the
fit actually was. Lifted out of `predict_selected`'s dispatch chain along with the rest of the
hurdle adapters.
"""
function _predict_hurdle_selected(
    selected, train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64}, y_sat_target::Matrix{Float64},
    hurdle_context,
)
    hurdle_context === nothing &&
        throw(ArgumentError("hurdle_gwr requires a hurdle_context carrying the fold times"))
    prediction = _hurdle_predict(
        train_lonlat, y_obs_train, y_sat_train, target_lonlat, y_sat_target,
        hurdle_context.times;
        config=_hurdle_config(
            hurdle_context, selected.kernel, selected.adaptive, selected.bw,
        ),
    )
    push!(hurdle_context.diagnostics, _hurdle_diagnostic_row(
        prediction, hurdle_context;
        kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw,
    ))
    # `_hurdle_predict` returns the full model output; the benchmark wants the corrected field.
    return prediction.corrected
end
