"""
Accuracy metrics shared by every evaluation path in the project.

These lived in `MGERPipeline.jl`, which is a top-level script rather than a module, so nothing
that *is* a module could reach them - `BenchmarkDiagnostics` says as much where it reimplements
`metric_continuous` (see its module docstring). Moved into `MixedGWR` so there is one definition
to call. The bodies are unchanged, including their edge cases: `metric_continuous` asserts on an
empty sample rather than returning a NaN row.
"""

function metric_continuous(y_true::AbstractArray{<:Real}, y_pred::AbstractArray{<:Real}; mask=nothing)
	a = Float64.(vec(y_true))
	b = Float64.(vec(y_pred))
	use_mask = mask === nothing ? (.!isnan.(a) .& .!isnan.(b)) : vec(mask)
	aa = a[use_mask]
	bb = b[use_mask]
	n = length(aa)
	@assert n > 0 "没有可用样本用于统计"
	e = bb .- aa
	rmse = sqrt(mean(e .^ 2))
	mae = mean(abs.(e))
	bias = mean(e)
	r = n > 1 ? cor(aa, bb) : NaN
	return (; n, RMSE=rmse, MAE=mae, Bias=bias, r)
end


function metric_event(y_true::AbstractArray{<:Real}, y_pred::AbstractArray{<:Real}; thr::Float64=0.1, mask=nothing)
	a = Float64.(vec(y_true))
	b = Float64.(vec(y_pred))
	use_mask = mask === nothing ? (.!isnan.(a) .& .!isnan.(b)) : vec(mask)
	aa = a[use_mask]
	bb = b[use_mask]

	obs = aa .>= thr
	est = bb .>= thr
	hit = sum(obs .& est)
	miss = sum(obs .& .!est)
	fal = sum((.!obs) .& est)

	pod_den = hit + miss
	far_den = hit + fal
	csi_den = hit + miss + fal

	pod = pod_den > 0 ? hit / pod_den : NaN
	far = far_den > 0 ? fal / far_den : NaN
	csi = csi_den > 0 ? hit / csi_den : NaN
	return (; POD=pod, FAR=far, CSI=csi)
end


function common_valid_mask(arrays::AbstractArray...)
	mask = trues(size(arrays[1]))
	for a in arrays
		size(a) == size(mask) ||
			throw(DimensionMismatch("all arrays must have the same size for a common metric mask"))
		mask .&= .!isnan.(a)
	end
	return mask
end


function complete_time_mask(arrays::AbstractMatrix...)
	mask = trues(size(arrays[1], 2))
	for a in arrays
		size(a, 2) == length(mask) ||
			throw(DimensionMismatch("all matrices must have the same time dimension"))
		mask .&= vec(all(.!isnan.(a), dims=1))
	end
	return mask
end
