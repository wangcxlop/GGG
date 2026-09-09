# ------------------------------------------- D5: which cells does a method fail to predict?

"""
Rebuild `common_evaluation_mask.csv` from the out-of-fold matrices: a cell counts only when
the gauge value and every mask-defining method's prediction are all finite.

Mirrors `InterpolationBenchmark._common_method_mask`. Pass a shorter `methods` list to ask
what the mask would have been had a given method not been allowed to define it.
"""
function rebuild_common_mask(
    y_obs::Matrix{Float64}, predictions::AbstractDict{String,Matrix{Float64}};
    methods::Vector{String}=MASK_METHODS,
)
    mask = .!isnan.(y_obs)
    for method in methods
        haskey(predictions, method) ||
            throw(ArgumentError("mask-defining method is absent from the run: $method"))
        mask .&= .!isnan.(predictions[method])
    end
    return BitMatrix(mask)
end

"""
Per fold, how many cells a method failed to predict, and whether the failures come in whole
hours or scattered station-by-station.

The denominator is the *evaluable* cell count — gauge and satellite both finite — not the
shared mask, because the shared mask has already had these very failures removed from it.
A whole-hour dropout (every evaluable validation station NaN in the same hour) is the
signature of a fit-level failure such as back-fitting non-convergence, which discards the
hour for all targets at once; scattered dropouts instead mean individual local fits were
underdetermined. The split is what distinguishes the two, and they need different fixes.
"""
function dropout_table(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64},
    predictions::AbstractDict{String,Matrix{Float64}}, ids::Vector{String},
    fold_map::Dict{String,Int}; scheme::String, product::String,
    methods::Vector{String}=["mgwr", "residual_gwr", "mixed_gwr"],
)
    fold_of = [fold_map[id] for id in ids]
    n_time = size(y_obs, 2)
    rows = NamedTuple[]
    for method in methods
        haskey(predictions, method) || continue
        prediction = predictions[method]
        for fold in sort(unique(fold_of))
            val = findall(==(fold), fold_of)
            isempty(val) && continue
            evaluable = 0
            dropped = 0
            hours_evaluable = 0
            hours_all_dropped = 0
            cells_in_dropped_hours = 0
            @inbounds for time in 1:n_time
                hour_evaluable = 0
                hour_dropped = 0
                for station in val
                    (isnan(y_obs[station, time]) || isnan(y_sat[station, time])) && continue
                    hour_evaluable += 1
                    isnan(prediction[station, time]) && (hour_dropped += 1)
                end
                hour_evaluable == 0 && continue
                hours_evaluable += 1
                evaluable += hour_evaluable
                dropped += hour_dropped
                if hour_dropped == hour_evaluable
                    hours_all_dropped += 1
                    cells_in_dropped_hours += hour_dropped
                end
            end
            push!(rows, (;
                scheme, product, fold, method,
                n_evaluable=evaluable, n_dropped=dropped,
                own_coverage=evaluable > 0 ? 1 - dropped / evaluable : NaN,
                hours_evaluable, hours_all_dropped,
                hour_dropout_fraction=hours_evaluable > 0 ?
                    hours_all_dropped / hours_evaluable : NaN,
                cells_in_all_dropped_hours=cells_in_dropped_hours,
                cells_scattered=dropped - cells_in_dropped_hours,
            ))
        end
    end
    return DataFrame(rows)
end

"""
What one mask-defining method's failures cost every other method.

The shared mask keeps only cells where *every* method in `MASK_METHODS` is finite, so an
unstable method silently shrinks the denominator the whole field is scored on. This rebuilds
the mask with `excluded` removed from the defining set and rescores everyone on both, so the
cost is a number rather than an argument. Cells recovered are by construction the harder ones
(the excluded method could not fit them), so a method's RMSE is expected to rise on the wider
mask — the comparison of interest is between methods, not against zero.
"""
function mask_cost_table(
    y_obs::Matrix{Float64}, predictions::AbstractDict{String,Matrix{Float64}};
    scheme::String, product::String, excluded::String="mgwr",
    methods::Vector{String}=MASK_METHODS,
)
    excluded in methods ||
        throw(ArgumentError("$excluded does not define the mask, so removing it changes nothing"))
    full = rebuild_common_mask(y_obs, predictions; methods)
    reduced = rebuild_common_mask(y_obs, predictions; methods=filter(!=(excluded), methods))
    rows = NamedTuple[]
    for method in sort(collect(keys(predictions)))
        on_full = _metrics(y_obs, predictions[method], full)
        on_reduced = _metrics(y_obs, predictions[method], reduced)
        push!(rows, (;
            scheme, product, excluded_method=excluded, method,
            mask_cells_full=count(full), mask_cells_reduced=count(reduced),
            cells_recovered=count(reduced) - count(full),
            n_full=on_full.n, RMSE_full=on_full.RMSE,
            n_reduced=on_reduced.n, RMSE_reduced=on_reduced.RMSE,
            RMSE_delta=on_reduced.RMSE - on_full.RMSE,
            relative_change=(isnan(on_full.RMSE) || on_full.RMSE == 0) ? NaN :
                (on_reduced.RMSE - on_full.RMSE) / on_full.RMSE,
        ))
    end
    return DataFrame(rows)
end

