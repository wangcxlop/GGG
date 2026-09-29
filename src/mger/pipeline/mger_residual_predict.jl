function pairwise_haversine_km(X_lonlat::Matrix{Float64})
	points = map(x -> (x[1], x[2]), eachrow(X_lonlat))
	fun_dist = Haversine(6378.388)
	return pairwise(fun_dist, points)
end


function make_loocv_dist(dMat::Matrix{Float64})
	d = copy(dMat)
	@inbounds for i in 1:min(size(d)...)
		d[i, i] = Inf
	end
	return d
end


# `metric_continuous`, `metric_event`, `common_valid_mask` and `complete_time_mask` now live in
# `src/core/metrics.jl` inside the `MixedGWR` module, which this file already imports at the top. They
# are reused by the benchmark, the diagnostics and the tests, so they are library code rather than
# part of this pipeline.




function st_gwr_predict_nanaware(
	X::Matrix{Float64}, Y::Matrix{Float64}, wMat::Matrix{Float64};
	Xpred::Matrix{Float64}, min_obs::Int=size(X, 2), lambda::Float64=1e-8,
)
	n_control, p = size(X)
	size(Y, 1) == n_control ||
		throw(DimensionMismatch("Y must have the same number of rows as X"))
	size(wMat, 1) == n_control ||
		throw(DimensionMismatch("weight matrix rows must match X rows"))
	size(Xpred, 2) == p ||
		throw(DimensionMismatch("Xpred must have the same number of columns as X"))
	size(Xpred, 1) == size(wMat, 2) ||
		throw(DimensionMismatch("Xpred rows must match weight matrix columns"))
	p == 3 ||
		throw(ArgumentError("nan-aware ST-GWR currently expects X = [1, lon, lat] with 3 columns"))

	n_target = size(wMat, 2)
	ntime = size(Y, 2)
	Ypred = fill(NaN, n_target, ntime)

	# `:greedy` rather than the default chunked schedule: `n_target` is a fold's held-out station
	# count (tens), so equal chunks leave threads idle on the ragged last round. Each iteration
	# writes only `Ypred[i, :]`, so the schedule cannot change a value.
	@inbounds Threads.@threads :greedy for i in 1:n_target
		xp1 = Xpred[i, 1]
		xp2 = Xpred[i, 2]
		xp3 = Xpred[i, 3]
		for t in 1:ntime
			a11 = lambda
			a12 = 0.0
			a13 = 0.0
			a22 = lambda
			a23 = 0.0
			a33 = lambda
			b1 = 0.0
			b2 = 0.0
			b3 = 0.0
			n_eff = 0

			for r in 1:n_control
				y = Y[r, t]
				w = wMat[r, i]
				if !isnan(y) && isfinite(w) && w > 0
					x1 = X[r, 1]
					x2 = X[r, 2]
					x3 = X[r, 3]
					wy = w * y
					a11 += w * x1 * x1
					a12 += w * x1 * x2
					a13 += w * x1 * x3
					a22 += w * x2 * x2
					a23 += w * x2 * x3
					a33 += w * x3 * x3
					b1 += wy * x1
					b2 += wy * x2
					b3 += wy * x3
					n_eff += 1
				end
			end

			if n_eff >= min_obs
				det = a11 * (a22 * a33 - a23 * a23) -
					a12 * (a12 * a33 - a13 * a23) +
					a13 * (a12 * a23 - a13 * a22)
				if abs(det) > eps(Float64)
					beta1 = (b1 * (a22 * a33 - a23 * a23) -
						a12 * (b2 * a33 - a23 * b3) +
						a13 * (b2 * a23 - a22 * b3)) / det
					beta2 = (a11 * (b2 * a33 - a23 * b3) -
						b1 * (a12 * a33 - a13 * a23) +
						a13 * (a12 * b3 - b2 * a13)) / det
					beta3 = (a11 * (a22 * b3 - b2 * a23) -
						a12 * (a12 * b3 - b2 * a13) +
						b1 * (a12 * a23 - a13 * a22)) / det
					Ypred[i, t] = xp1 * beta1 + xp2 * beta2 + xp3 * beta3
				end
			end
		end
	end
	return Ypred
end


"""Return target-centred east/north offsets in kilometres for training stations."""
function target_centered_offsets_km(
	train_lonlat::Matrix{Float64}, target_lon::Float64, target_lat::Float64;
	earth_radius_km::Float64=6378.388,
)
	size(train_lonlat, 2) == 2 ||
		throw(DimensionMismatch("train_lonlat must have columns [lon, lat]"))
	isfinite(target_lon) && isfinite(target_lat) ||
		throw(ArgumentError("target longitude and latitude must be finite"))

	n = size(train_lonlat, 1)
	east = Vector{Float64}(undef, n)
	north = Vector{Float64}(undef, n)
	cos_lat = cosd(target_lat)
	@inbounds for i in 1:n
		lon = train_lonlat[i, 1]
		lat = train_lonlat[i, 2]
		if isfinite(lon) && isfinite(lat)
			dlon = mod(lon - target_lon + 180.0, 360.0) - 180.0
			east[i] = earth_radius_km * cos_lat * deg2rad(dlon)
			north[i] = earth_radius_km * deg2rad(lat - target_lat)
		else
			east[i] = NaN
			north[i] = NaN
		end
	end
	return east, north
end


"""
Predict a time-varying residual at each target with target-centred local linear GWR.

For each target and time, weights are normalised over valid training residuals and
east/north offsets are scaled by their weighted RMS distance. Ridge regularisation
is applied only to the two spatial slopes, so the target prediction is the local
intercept. Invalid or underdetermined fits remain `NaN`.
"""
function local_linear_residual_predict(
	train_lonlat::Matrix{Float64}, residuals::Matrix{Float64},
	target_lonlat::Matrix{Float64}, wMat::Matrix{Float64};
	slope_ridge::Float64, min_obs::Int=3,
)
	n_control = size(train_lonlat, 1)
	size(train_lonlat, 2) == 2 ||
		throw(DimensionMismatch("train_lonlat must have columns [lon, lat]"))
	size(target_lonlat, 2) == 2 ||
		throw(DimensionMismatch("target_lonlat must have columns [lon, lat]"))
	size(residuals, 1) == n_control ||
		throw(DimensionMismatch("residual rows must match training stations"))
	size(wMat, 1) == n_control ||
		throw(DimensionMismatch("weight matrix rows must match training stations"))
	size(wMat, 2) == size(target_lonlat, 1) ||
		throw(DimensionMismatch("weight matrix columns must match target stations"))
	min_obs >= 3 || throw(ArgumentError("min_obs must be at least 3"))
	isfinite(slope_ridge) && slope_ridge >= 0 ||
		throw(ArgumentError("slope_ridge must be finite and non-negative"))

	n_target = size(target_lonlat, 1)
	ntime = size(residuals, 2)
	prediction = fill(NaN, n_target, ntime)

	@inbounds Threads.@threads for j in 1:n_target
		target_lon = target_lonlat[j, 1]
		target_lat = target_lonlat[j, 2]
		if !isfinite(target_lon) || !isfinite(target_lat)
			continue
		end
		east, north = target_centered_offsets_km(
			train_lonlat, target_lon, target_lat,
		)

		for t in 1:ntime
			n_eff = 0
			sum_w = 0.0
			for i in 1:n_control
				y = residuals[i, t]
				w = wMat[i, j]
				if isfinite(y) && isfinite(w) && w > 0 &&
					isfinite(east[i]) && isfinite(north[i])
					n_eff += 1
					sum_w += w
				end
			end
			n_eff >= min_obs && isfinite(sum_w) && sum_w > 0 || continue

			scale2 = 0.0
			for i in 1:n_control
				y = residuals[i, t]
				w = wMat[i, j]
				if isfinite(y) && isfinite(w) && w > 0 &&
					isfinite(east[i]) && isfinite(north[i])
					wn = w / sum_w
					scale2 += wn * (east[i]^2 + north[i]^2)
				end
			end
			isfinite(scale2) && scale2 > eps(Float64) || continue
			scale = sqrt(scale2)

			a11 = 0.0
			a12 = 0.0
			a13 = 0.0
			a22 = slope_ridge
			a23 = 0.0
			a33 = slope_ridge
			b1 = 0.0
			b2 = 0.0
			b3 = 0.0
			for i in 1:n_control
				y = residuals[i, t]
				w = wMat[i, j]
				if isfinite(y) && isfinite(w) && w > 0 &&
					isfinite(east[i]) && isfinite(north[i])
					wn = w / sum_w
					x2 = east[i] / scale
					x3 = north[i] / scale
					a11 += wn
					a12 += wn * x2
					a13 += wn * x3
					a22 += wn * x2 * x2
					a23 += wn * x2 * x3
					a33 += wn * x3 * x3
					wy = wn * y
					b1 += wy
					b2 += wy * x2
					b3 += wy * x3
				end
			end

			det = a11 * (a22 * a33 - a23 * a23) -
				a12 * (a12 * a33 - a13 * a23) +
				a13 * (a12 * a23 - a13 * a22)
			if isfinite(det) && abs(det) > eps(Float64)
				beta0 = (b1 * (a22 * a33 - a23 * a23) -
					a12 * (b2 * a33 - a23 * b3) +
					a13 * (b2 * a23 - a22 * b3)) / det
				isfinite(beta0) && (prediction[j, t] = beta0)
			end
		end
	end
	return prediction
end


function negative_output_stats(values::AbstractArray{<:Real})
	finite_mask = isfinite.(values)
	n_finite = count(finite_mask)
	if n_finite == 0
		return (; negative_n=0, negative_fraction=NaN, min_corrected=NaN)
	end
	negative_n = count(finite_mask .& (values .< 0))
	return (;
		negative_n,
		negative_fraction=negative_n / n_finite,
		min_corrected=minimum(values[finite_mask]),
	)
end


# Core correction step

function bias_correct_stgwr(
	lonlat::Matrix{Float64}, Y_obs::Matrix{Float64}, Y_sat::Matrix{Float64},
	dMat::Matrix{Float64}; kernel::Int=BISQUARE, adaptive::Bool=true, bw::Float64=50.0,
	use_loocv::Bool=true, slope_ridge::Float64=1e-6,
)
	dist = use_loocv ? make_loocv_dist(dMat) : dMat
	wMat = gw_weight(dist, bw; kernel, adaptive)
	R = Y_obs .- Y_sat
	Rhat = local_linear_residual_predict(
		lonlat, R, lonlat, wMat; slope_ridge,
	)
	Y_corr = Y_sat .+ Rhat
	return Y_corr, Rhat, wMat
end


function scan_params(
	lonlat::Matrix{Float64}, Y_obs::Matrix{Float64}, Y_sat::Matrix{Float64}, dMat::Matrix{Float64};
	kernels::Vector{Int}, bw_adaptive::Vector{Float64}, bw_fixed_km::Vector{Float64},
	slope_ridge_candidates::Vector{Float64},
	min_scan_coverage::Float64=0.95,
	rain_threshold::Float64=0.1, use_loocv::Bool=true, fail_path::Union{Nothing,String}=nothing,
)
	isempty(slope_ridge_candidates) &&
		throw(ArgumentError("slope_ridge_candidates must not be empty"))
	all(isfinite(ridge) && ridge >= 0 for ridge in slope_ridge_candidates) ||
		throw(ArgumentError("slope_ridge_candidates must be finite and non-negative"))
	isfinite(min_scan_coverage) && 0.0 <= min_scan_coverage <= 1.0 ||
		throw(ArgumentError("min_scan_coverage must be between 0 and 1"))
	rows = NamedTuple[]
	fail_df = DataFrame(
		kernel=Int[], adaptive=Bool[], bw=Float64[], slope_ridge=Float64[], error=String[],
	)
	for k in kernels
		for adaptive in (true, false)
			bandwidths = adaptive ? bw_adaptive : bw_fixed_km
			for bw in bandwidths, slope_ridge in slope_ridge_candidates
				try
					Yc, _, _ = bias_correct_stgwr(
						lonlat, Y_obs, Y_sat, dMat;
						kernel=k, adaptive, bw, use_loocv, slope_ridge,
					)
					eval_mask = common_valid_mask(Y_obs, Y_sat, Yc)
					mc = metric_continuous(Y_obs, Yc; mask=eval_mask)
					me = metric_event(Y_obs, Yc; thr=rain_threshold, mask=eval_mask)
					neg = negative_output_stats(Yc)
					push!(rows, (;
						kernel=k, adaptive, bw, slope_ridge,
						n=mc.n, coverage=mc.n / length(eval_mask),
						RMSE=mc.RMSE, AICc=NaN, criterion="LOOCV_RMSE",
						MAE=mc.MAE, Bias=mc.Bias, r=mc.r,
						POD=me.POD, FAR=me.FAR, CSI=me.CSI,
						neg...,
					))
				catch e
					push!(fail_df, (;
						kernel=k, adaptive, bw, slope_ridge,
						error=sprint(showerror, e),
					))
				end
			end
		end
	end
	if fail_path !== nothing
		CSV.write(fail_path, fail_df)
	end
	@assert !isempty(rows) "every parameter combination failed (usually a singular local regression matrix). Try a larger bandwidth or adaptive bandwidths only."
	df = DataFrame(rows)
	if nrow(fail_df) > 0
		println("[scan_params] failed combinations skipped: ", nrow(fail_df))
	end
	eligible = df[df.coverage .>= min_scan_coverage, :]
	if isempty(eligible)
		max_coverage = maximum(df.coverage)
		throw(ArgumentError(
			"no parameter candidate reached min_scan_coverage=$min_scan_coverage; " *
			"maximum coverage was $max_coverage",
		))
	end
	sort!(df, [:RMSE, :MAE, :slope_ridge, :bw, :adaptive])
	sort!(eligible, [:RMSE, :MAE, :slope_ridge, :bw, :adaptive])
	best = eligible[1, :]
	return df, best
end
