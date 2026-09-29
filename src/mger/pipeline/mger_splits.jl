"""
    split_stations_train_val(station_ids::Vector{String}; train_frac::Float64=0.8, rng::AbstractRNG=Random.GLOBAL_RNG)

Split stations into training and validation sets by a random shuffle.
Returns: (train_ids, val_ids)

The held-out stations are new/unseen *locations*, so this tests generalisation to stations that
were not fitted. It is NOT a spatial split: the shuffle never reads `lonlat`, so held-out
stations sit interleaved with training ones and a validation station typically has a training
neighbour a few km away. For held-out *area*, use `split_stations_spatial_block_kfold`.
"""

# Station-wise train/validation splits

function split_stations_train_val(station_ids::Vector{String}; train_frac::Float64=0.8, rng::AbstractRNG=Random.GLOBAL_RNG)
    n = length(station_ids)
    n_train = max(1, floor(Int, n * train_frac))

    shuffled = shuffle(rng, station_ids)
    train_ids = shuffled[1:n_train]
    val_ids = shuffled[n_train+1:end]

    return train_ids, val_ids
end


function split_stations_kfold(station_ids::Vector{String}; k::Int=5, rng::AbstractRNG=Random.GLOBAL_RNG)
	n = length(station_ids)
	2 <= k <= n || throw(ArgumentError("k must be between 2 and the number of stations"))

	shuffled = shuffle(rng, station_ids)
	folds = [String[] for _ in 1:k]
	for (i, station_id) in enumerate(shuffled)
		push!(folds[mod1(i, k)], station_id)
	end
	return folds
end


function split_stations_spatial_block_kfold(
	station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, axis::Symbol=:auto,
)
	n = length(station_ids)
	2 <= k <= n || throw(ArgumentError("k must be between 2 and the number of stations"))
	size(lonlat, 1) == n || throw(DimensionMismatch("lonlat rows must match station_ids"))

	chosen_axis = if axis == :auto
		lat_mid = mean(lonlat[:, 2])
		lon_span_km = (maximum(lonlat[:, 1]) - minimum(lonlat[:, 1])) * cosd(lat_mid) * 111.32
		lat_span_km = (maximum(lonlat[:, 2]) - minimum(lonlat[:, 2])) * 110.57
		lon_span_km >= lat_span_km ? :lon : :lat
	elseif axis in (:lon, :lat)
		axis
	else
		throw(ArgumentError("axis must be :auto, :lon, or :lat"))
	end

	col = chosen_axis == :lon ? 1 : 2
	order = sortperm(1:n; by=i -> (lonlat[i, col], lonlat[i, 3 - col], station_ids[i]))
	base = div(n, k)
	extra = rem(n, k)
	folds = Vector{String}[]
	start = 1
	for fold_idx in 1:k
		len = base + (fold_idx <= extra ? 1 : 0)
		idx = order[start:start+len-1]
		push!(folds, station_ids[idx])
		start += len
	end
	return folds
end


function evaluate_spatial_holdout(
	product::AbstractString, ids::Vector{String}, times::Vector{DateTime},
	Y_obs::Matrix{Float64}, Y_sat::Matrix{Float64}, lonlat::Matrix{Float64},
	train_ids::Vector{String}, val_ids::Vector{String}, cfg::MGERConfig;
	scan_path::AbstractString, fail_path::AbstractString,
)
	id_map = Dict(id => i for (i, id) in enumerate(ids))
	train_idx = [id_map[id] for id in train_ids]
	val_idx = [id_map[id] for id in val_ids]

	lonlat_train = lonlat[train_idx, :]
	lonlat_val = lonlat[val_idx, :]
	Y_obs_train = Y_obs[train_idx, :]
	Y_sat_train = Y_sat[train_idx, :]
	Y_obs_val = Y_obs[val_idx, :]
	Y_sat_val = Y_sat[val_idx, :]

	dMat_train = pairwise_haversine_km(lonlat_train)
	dMat_train_val = let
		points_train = map(x -> (x[1], x[2]), eachrow(lonlat_train))
		points_val = map(x -> (x[1], x[2]), eachrow(lonlat_val))
		pairwise(Haversine(6378.388), points_train, points_val)
	end

	scan_df, best = scan_params(
		lonlat_train, Y_obs_train, Y_sat_train, dMat_train;
		kernels=cfg.kernels,
		bw_adaptive=cfg.bw_adaptive,
		bw_fixed_km=cfg.bw_fixed_km,
		slope_ridge_candidates=cfg.slope_ridge_candidates,
		min_scan_coverage=cfg.min_scan_coverage,
		rain_threshold=cfg.rain_threshold,
		use_loocv=cfg.use_loocv_eval,
		fail_path=fail_path,
	)
	CSV.write(scan_path, scan_df)

	wMat_train_val = gw_weight(dMat_train_val, Float64(best.bw); kernel=Int(best.kernel), adaptive=Bool(best.adaptive))
	R_train = Y_obs_train .- Y_sat_train
	Rhat_val = local_linear_residual_predict(
		lonlat_train, R_train, lonlat_val, wMat_train_val;
		slope_ridge=Float64(best.slope_ridge),
	)
	Y_corr_val = Y_sat_val .+ Rhat_val
	negative = negative_output_stats(Y_corr_val)

	pre_mask = common_valid_mask(Y_obs_val, Y_sat_val)
	post_mask = common_valid_mask(Y_obs_val, Y_corr_val)
	eval_mask = common_valid_mask(Y_obs_val, Y_sat_val, Y_corr_val)
	@assert any(eval_mask) "[$product] validation set has no common valid samples; cannot compare before/after correction"

	pre_c = metric_continuous(Y_obs_val, Y_sat_val; mask=eval_mask)
	pre_e = metric_event(Y_obs_val, Y_sat_val; thr=cfg.rain_threshold, mask=eval_mask)
	post_c = metric_continuous(Y_obs_val, Y_corr_val; mask=eval_mask)
	post_e = metric_event(Y_obs_val, Y_corr_val; thr=cfg.rain_threshold, mask=eval_mask)

	return (;
		val_idx, Y_corr_val,
		n_train=length(train_ids),
		n_val=length(val_ids),
		n=post_c.n,
		pre_n=count(pre_mask),
		post_n=count(post_mask),
		common_n=count(eval_mask),
		total_n=length(eval_mask),
		coverage=count(eval_mask) / length(eval_mask),
		scan_n=Int(best.n),
		scan_coverage=Float64(best.coverage),
		kernel=Int(best.kernel),
		adaptive=Bool(best.adaptive),
		bw=Float64(best.bw),
		slope_ridge=Float64(best.slope_ridge),
		RMSE_pre=pre_c.RMSE, RMSE_post=post_c.RMSE,
		MAE_pre=pre_c.MAE, MAE_post=post_c.MAE,
		Bias_pre=pre_c.Bias, Bias_post=post_c.Bias,
		r_pre=pre_c.r, r_post=post_c.r,
		POD_pre=pre_e.POD, POD_post=post_e.POD,
		FAR_pre=pre_e.FAR, FAR_post=post_e.FAR,
		CSI_pre=pre_e.CSI, CSI_post=post_e.CSI,
		negative...,
	)
end


function summarize_fold_metrics(fold_df::DataFrame)
	metric_cols = [
		:RMSE_pre, :RMSE_post, :MAE_pre, :MAE_post, :Bias_pre, :Bias_post,
		:r_pre, :r_post, :POD_pre, :POD_post, :FAR_pre, :FAR_post, :CSI_pre, :CSI_post,
	]
	out = DataFrame()
	for sdf in groupby(fold_df, :product)
		row = Dict{Symbol, Any}(
			:product => first(sdf.product),
			:k => nrow(sdf),
			:n_mean => mean(Float64.(sdf.n)),
			:n_sum => sum(Int.(sdf.n)),
		)
		for col in metric_cols
			vals = Float64.(sdf[!, col])
			row[Symbol(String(col), "_mean")] = mean(vals)
			row[Symbol(String(col), "_std")] = length(vals) > 1 ? std(vals) : NaN
		end
		push!(out, row; cols=:union)
	end
	return out
end
