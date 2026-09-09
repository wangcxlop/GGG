# ---------------------------------------------------------------------------- scoring

function _metrics(y_true::AbstractMatrix, y_pred::AbstractMatrix, mask::AbstractMatrix)
    a = Float64[]
    b = Float64[]
    @inbounds for index in eachindex(mask)
        mask[index] || continue
        observed = y_true[index]
        predicted = y_pred[index]
        (isnan(observed) || isnan(predicted)) && continue
        push!(a, observed)
        push!(b, predicted)
    end
    n = length(a)
    n == 0 && return (; n=0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN, MSE=NaN, variance=NaN)
    e = b .- a
    bias = mean(e)
    mse = mean(e .^ 2)
    return (;
        n, RMSE=sqrt(mse), MAE=mean(abs.(e)), Bias=bias,
        r=n > 1 ? cor(a, b) : NaN, MSE=mse, variance=mse - bias^2,
    )
end

"""`station × time` mask selecting the cells of `mask` whose gauge value falls in `[low, high)`."""
function _rain_mask(y_obs::AbstractMatrix, mask::AbstractMatrix, low::Float64, high::Float64)
    out = falses(size(mask))
    @inbounds for index in eachindex(mask)
        mask[index] || continue
        observed = y_obs[index]
        isnan(observed) && continue
        out[index] = low <= observed < high
    end
    return out
end

