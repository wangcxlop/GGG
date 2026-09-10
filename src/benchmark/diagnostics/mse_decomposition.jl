# ------------------------------------------------------- D3: where the MSE gap lives

"""
Per-method MSE split into bias² + variance within each rain class, plus each class's share
of the total MSE gap against `reference` (the best traditional interpolator).
"""
function mse_decomposition_table(
    y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    mask::AbstractMatrix; scheme::String, product::String, reference::String="adw",
)
    overall_n = count(mask)
    reference_overall = _metrics(y_obs, predictions[reference], mask)
    rows = NamedTuple[]
    for method in sort(collect(keys(predictions)))
        method == reference && continue
        overall = _metrics(y_obs, predictions[method], mask)
        # Total MSE gap, weighted the way the pooled RMSE weights it.
        total_gap = overall.MSE - reference_overall.MSE
        for (name, low, high) in RAIN_CLASSES
            stratum = _rain_mask(y_obs, mask, low, high)
            score = _metrics(y_obs, predictions[method], stratum)
            reference_score = _metrics(y_obs, predictions[reference], stratum)
            score.n == 0 && continue
            weight = score.n / overall_n
            contribution = weight * (score.MSE - reference_score.MSE)
            push!(rows, (;
                scheme, product, method, reference, rain_class=name, n=score.n,
                sample_share=weight, RMSE=score.RMSE, MSE=score.MSE,
                bias_squared=score.Bias^2, variance=score.variance,
                reference_MSE=reference_score.MSE,
                mse_gap_contribution=contribution,
                gap_share=total_gap != 0 ? contribution / total_gap : NaN,
            ))
        end
    end
    return DataFrame(rows)
end

