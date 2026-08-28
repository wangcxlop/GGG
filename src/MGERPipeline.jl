using MixedGWR
using CSV, DataFrames, Dates, Statistics, Random
using Distances


"""
配置参数
- 残差模型使用目标中心的局部东西/南北公里坐标
- 同时扫描 adaptive=true/false 两类带宽
- 同时扫描局部空间斜率的岭正则强度
"""
Base.@kwdef struct MGERConfig
	station_meta_path::String
	obs_hourly_wide_path::String
	sat_paths::Dict{String, String} # "FY4B"=>path, "GPM"=>path, "GSMaP"=>path
	outdir::String = "output/mger_final"

	time_col::Symbol = :time
	station_id_col::Symbol = :station_id
	lon_col::Symbol = :lon
	lat_col::Symbol = :lat

	kernels::Vector{Int} = [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR]
	bw_adaptive::Vector{Float64} = [30.0, 50.0, 80.0, 120.0]
	bw_fixed_km::Vector{Float64} = [10.0, 20.0, 30.0, 50.0]
	slope_ridge_candidates::Vector{Float64} = [1e-8, 1e-6, 1e-4, 1e-2]
	min_scan_coverage::Float64 = 0.95
	rain_threshold::Float64 = 0.1

	use_loocv_eval::Bool = true
	analysis_start::Union{Nothing, DateTime} = nothing
	analysis_end::Union{Nothing, DateTime} = nothing
	expected_common_time_count::Union{Nothing, Int} = nothing
end


function config_for_kernel(cfg::MGERConfig, kernel::Int, outdir::AbstractString)
	return MGERConfig(
		station_meta_path=cfg.station_meta_path,
		obs_hourly_wide_path=cfg.obs_hourly_wide_path,
		sat_paths=copy(cfg.sat_paths),
		outdir=String(outdir),
		time_col=cfg.time_col,
		station_id_col=cfg.station_id_col,
		lon_col=cfg.lon_col,
		lat_col=cfg.lat_col,
		kernels=[kernel],
		bw_adaptive=copy(cfg.bw_adaptive),
		bw_fixed_km=copy(cfg.bw_fixed_km),
		slope_ridge_candidates=copy(cfg.slope_ridge_candidates),
		min_scan_coverage=cfg.min_scan_coverage,
		rain_threshold=cfg.rain_threshold,
		use_loocv_eval=cfg.use_loocv_eval,
		analysis_start=cfg.analysis_start,
		analysis_end=cfg.analysis_end,
		expected_common_time_count=cfg.expected_common_time_count,
	)
end


function parse_time_utc(s::AbstractString)
	t = replace(strip(s), "Z" => "")
	fmts = (
		dateformat"yyyy-mm-ddTHH:MM:SS",
		dateformat"yyyy-mm-dd HH:MM:SS",
		dateformat"yyyy/mm/dd HH:MM:SS",
		dateformat"yyyy-mm-ddTHH:MM",
		dateformat"yyyy-mm-dd HH:MM",
	)
	for f in fmts
		try
			return DateTime(t, f)
		catch
		end
	end
	error("无法解析时间: $s")
end


"""读取宽表: 行=时间，列=站点，返回 Y[站点, 时间]"""
function read_hourly_wide(path::AbstractString; time_col::Symbol=:time)
	df = CSV.read(path, DataFrame)
	names_lower = Dict(Symbol(lowercase(String(c))) => c for c in names(df))
	c_time = get(names_lower, Symbol(lowercase(String(time_col))), nothing)
	@assert c_time !== nothing "缺少 time 列: $time_col"

	times = parse_time_utc.(string.(df[!, c_time]))
	st_cols = [c for c in names(df) if c != c_time]
	station_ids = string.(st_cols)

	raw = Matrix(df[:, st_cols]) # [ntime, nstation]
	vals = Array{Float64}(undef, size(raw))
	@inbounds for i in eachindex(raw)
		v = raw[i]
		vals[i] = ismissing(v) ? NaN : Float64(v)
	end
	Y = permutedims(vals) # [nstation, ntime]
	return times, station_ids, Y
end

function subset_time_window(
	times::Vector{DateTime},
	Y::AbstractMatrix,
	analysis_start::Union{Nothing, DateTime},
	analysis_end::Union{Nothing, DateTime},
)
	@assert size(Y, 2) == length(times) "时间轴长度与数据列数不一致"
	if analysis_start !== nothing && analysis_end !== nothing
		@assert analysis_start <= analysis_end "分析开始时间不得晚于结束时间"
	end
	indices = findall(time ->
		(analysis_start === nothing || time >= analysis_start) &&
		(analysis_end === nothing || time <= analysis_end),
		times,
	)
	return times[indices], Y[:, indices]
end

function write_global_time_qc(
	path::AbstractString,
	obs_times::Vector{DateTime},
	sat_inputs::Dict{String, Any};
	reference_product::Union{Nothing, String}=nothing,
)
	reference_product = reference_product === nothing ?
		(haskey(sat_inputs, "FY4B") ? "FY4B" : first(sort(collect(keys(sat_inputs))))) :
		reference_product
	@assert haskey(sat_inputs, reference_product) "缺少全局时间质控参考产品 $reference_product"
	reference_times = sort(unique(sat_inputs[reference_product].times))
	obs_set = Set(obs_times)
	qc = DataFrame(time=Dates.format.(reference_times, dateformat"yyyy-mm-ddTHH:MM:SS") .* "Z")
	qc[!, :observation_available] = [time in obs_set for time in reference_times]
	availability_columns = Symbol[:observation_available]
	for product in sort(collect(keys(sat_inputs)))
		column = Symbol(lowercase(product), "_available")
		product_times = Set(sat_inputs[product].times)
		qc[!, column] = [time in product_times for time in reference_times]
		push!(availability_columns, column)
	end
	qc[!, :all_available] = [all(qc[row, column] for column in availability_columns) for row in 1:nrow(qc)]
	CSV.write(path, qc)
	return qc, reference_times[qc.all_available]
end

function load_global_common_product_data(cfg::MGERConfig)
	obs_times, obs_ids, Y_obs0 = read_hourly_wide(cfg.obs_hourly_wide_path; time_col=cfg.time_col)
	obs_times, Y_obs0 = subset_time_window(obs_times, Y_obs0, cfg.analysis_start, cfg.analysis_end)
	@assert length(unique(obs_times)) == length(obs_times) "观测存在重复时间戳"

	products = sort(collect(keys(cfg.sat_paths)))
	sat_inputs = Dict{String, Any}()
	common_product_ids = copy(obs_ids)
	for product in products
		sat_times, sat_ids, Y_sat0 = read_hourly_wide(cfg.sat_paths[product]; time_col=cfg.time_col)
		sat_times, Y_sat0 = subset_time_window(sat_times, Y_sat0, cfg.analysis_start, cfg.analysis_end)
		@assert length(unique(sat_times)) == length(sat_times) "[$product] 存在重复时间戳"
		id_set = Set(sat_ids)
		common_product_ids = [id for id in common_product_ids if id in id_set]
		sat_inputs[product] = (; times=sat_times, ids=sat_ids, Y=Y_sat0)
	end
	@assert !isempty(common_product_ids) "三产品没有共同站点，无法进行公平横向评估"

	time_qc_path = joinpath(cfg.outdir, "global_common_time_qc.csv")
	_, common_product_times = write_global_time_qc(time_qc_path, obs_times, sat_inputs)
	if cfg.expected_common_time_count !== nothing && length(common_product_times) != cfg.expected_common_time_count
		error("全局共同时间数 $(length(common_product_times)) 与预期 $(cfg.expected_common_time_count) 不一致。差异清单已写入 $time_qc_path")
	end
	@assert !isempty(common_product_times) "三产品与观测没有全局共同时间"

	obs_id_map = Dict(id => index for (index, id) in enumerate(obs_ids))
	obs_time_map = Dict(time => index for (index, time) in enumerate(obs_times))
	obs_id_idx = [obs_id_map[id] for id in common_product_ids]
	obs_time_idx = [obs_time_map[time] for time in common_product_times]
	Y_obs_common = Y_obs0[obs_id_idx, obs_time_idx]

	product_data = Dict{String, Any}()
	for product in products
		input = sat_inputs[product]
		id_map = Dict(id => index for (index, id) in enumerate(input.ids))
		time_map = Dict(time => index for (index, time) in enumerate(input.times))
		id_idx = [id_map[id] for id in common_product_ids]
		time_idx = [time_map[time] for time in common_product_times]
		product_data[product] = (;
			times=common_product_times,
			ids=common_product_ids,
			Y_obs=Y_obs_common,
			Y_sat=input.Y[id_idx, time_idx],
		)
	end
	return products, common_product_ids, product_data
end


function load_station_meta(path::AbstractString;
	station_id_col::Symbol=:station_id, lon_col::Symbol=:lon, lat_col::Symbol=:lat)
	st = CSV.read(path, DataFrame)
	names_lower = Dict(Symbol(lowercase(String(c))) => c for c in names(st))
	c_sid = get(names_lower, Symbol(lowercase(String(station_id_col))), nothing)
	c_lon = get(names_lower, Symbol(lowercase(String(lon_col))), nothing)
	c_lat = get(names_lower, Symbol(lowercase(String(lat_col))), nothing)

	@assert c_sid !== nothing "缺少站号列: $station_id_col"
	@assert c_lon !== nothing "缺少经度列: $lon_col"
	@assert c_lat !== nothing "缺少纬度列: $lat_col"

	rename!(st, c_sid => :station_id, c_lon => :lon, c_lat => :lat)
	st[!, :station_id] = string.(st[!, :station_id])
	return st
end


# 目的是对两个已加载的数据集进行同步/重新索引，使它们在站点和时间维度上匹配

function align_obs_sat(
	obs_times::Vector{DateTime}, obs_ids::Vector{String}, Y_obs::Matrix{Float64},
	sat_times::Vector{DateTime}, sat_ids::Vector{String}, Y_sat::Matrix{Float64}
)
	common_ids = intersect(obs_ids, sat_ids)
	common_times = intersect(obs_times, sat_times)
	sort!(common_times)

	@assert !isempty(common_ids) "站号没有交集"
	@assert !isempty(common_times) "时间没有交集"

	obs_id_map = Dict(v => i for (i, v) in pairs(obs_ids))
	sat_id_map = Dict(v => i for (i, v) in pairs(sat_ids))
	obs_t_map = Dict(v => i for (i, v) in pairs(obs_times))
	sat_t_map = Dict(v => i for (i, v) in pairs(sat_times))

	i_obs = [obs_id_map[x] for x in common_ids]
	i_sat = [sat_id_map[x] for x in common_ids]
	t_obs = [obs_t_map[x] for x in common_times]
	t_sat = [sat_t_map[x] for x in common_times]

	Yo = Y_obs[i_obs, t_obs]
	Ys = Y_sat[i_sat, t_sat]
	return common_times, common_ids, Yo, Ys
end

"""按 common_ids 顺序构建经纬度坐标矩阵，用于距离计算。"""
function build_X_lonlat(st::DataFrame, common_ids::Vector{String};
	station_id_col::Symbol=:station_id, lon_col::Symbol=:lon, lat_col::Symbol=:lat)
	st_map = Dict(string(st[i, station_id_col]) => i for i in 1:nrow(st))
	idx = [st_map[sid] for sid in common_ids if haskey(st_map, sid)]
	@assert length(idx) == length(common_ids) "站点元数据不完整，部分站号缺失"
	lon = Float64.(st[idx, lon_col])
	lat = Float64.(st[idx, lat_col])
	X = hcat(lon, lat)
	return X
end


"""构建 ST-GWR 回归设计矩阵 X=[1, lon_centered, lat_centered]。"""
function build_X_intercept_centered(lonlat::Matrix{Float64};
	center::Tuple{Float64,Float64}=(mean(lonlat[:, 1]), mean(lonlat[:, 2])))
	lon = lonlat[:, 1]
	lat = lonlat[:, 2]
	lon_center, lat_center = center
	return hcat(ones(Float64, size(lonlat, 1)), lon .- lon_center, lat .- lat_center)
end


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
# `src/metrics.jl` inside the `MixedGWR` module, which this file already imports at the top. They
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


# 核心矫正步骤

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


"""
计算 AICc，带 trace 有效性检查。
如果 trace 过大导致分母非正，则使用 GCV（广义交叉验证）作为替代。

GCV 公式：GCV = n * RSS / (n - trace)^2
"""

# 计算迹（有效复杂度）、残差平方和，并根据条件选择 AICc 或 GCV 作为模型选择标准

function calc_aicc_or_gcv(
	X::Matrix{Float64}, Ytrue::Matrix{Float64}, Ypred::Matrix{Float64}, wMat::Matrix{Float64}
)
	size(Ytrue) == size(Ypred) ||
		throw(DimensionMismatch("Ytrue and Ypred must have the same size"))

	mask = .!isnan.(Ytrue) .& .!isnan.(Ypred)
	n = count(mask)
	if n <= 3
		return Inf, :INSUFFICIENT
	end

	size(wMat, 2) == size(X, 1) ||
		throw(DimensionMismatch("AICc/GCV requires a calibration weight matrix whose columns match X rows"))

	# The hat-matrix trace depends only on X and weights, not on the response values.
	y_trace = zeros(Float64, size(X, 1))
	trace_one_time = GWR_calib(X, y_trace, wMat).trace
	trace = trace_one_time * (n / size(X, 1))
	resid = Ytrue[mask] .- Ypred[mask]
	RSS = max(sum(abs2, resid), eps(Float64))

	# Check if AICc is computable
	den1 = n - trace
	den2 = n - 2 - trace

	if den1 <= 0 || den2 <= 0
		# Fall back to GCV when AICc is unstable
		# Use approximate effective degrees of freedom
		effective_df = min(trace, n - 3)
		gcv = n * RSS / (n - effective_df)^2
		return isfinite(gcv) ? gcv : Inf, :GCV
	end

	aicc = (log(RSS / den1) + log(2pi) + (n + trace) / den2) * n
	return isfinite(aicc) ? aicc : Inf, :AICc
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
	@assert !isempty(rows) "所有参数组合都失败（常见原因：局部回归矩阵奇异）。可尝试增大带宽或仅用 adaptive 带宽。"
	df = DataFrame(rows)
	if nrow(fail_df) > 0
		println("[scan_params] 跳过失败组合数: ", nrow(fail_df))
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

# 将矫正后的结果写回宽表格式，列=站点，行=时间   回宽CSV格式

function write_wide(path::AbstractString, times::Vector{DateTime}, station_ids::Vector{String}, Y::Matrix{Float64})
        n_station, n_time = size(Y)
        @assert n_station == length(station_ids)
        @assert n_time == length(times)

	df = DataFrame(time = Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS"))
	Yt = permutedims(Y) # [ntime, nstation]
        for (j, sid) in enumerate(station_ids)
                df[!, Symbol(sid)] = Yt[:, j]
        end
        tmp_path = string(path, ".tmp-", getpid(), ".csv")
        try
                CSV.write(tmp_path, df)
                mv(tmp_path, path; force=true)
        catch e
                isfile(tmp_path) && rm(tmp_path; force=true)
                rethrow(e)
        end
end


function write_fullfit_product(
        outdir::AbstractString, product::AbstractString, times::Vector{DateTime}, ids::Vector{String},
		lonlat::Matrix{Float64}, Y_obs::Matrix{Float64}, Y_sat::Matrix{Float64}, dMat::Matrix{Float64}, best,
)
	wMat = gw_weight(dMat, Float64(best.bw); kernel=Int(best.kernel), adaptive=Bool(best.adaptive))
	R = Y_obs .- Y_sat
	Rhat = local_linear_residual_predict(
		lonlat, R, lonlat, wMat; slope_ridge=Float64(best.slope_ridge),
	)
	Y_corr = Y_sat .+ Rhat
        path = joinpath(outdir, "corr_$(product)_fullfit_insample.csv")
        try
                write_wide(path, times, ids, Y_corr)
                return path
        catch e
                stamp = Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
                fallback_path = joinpath(outdir, "corr_$(product)_fullfit_insample_$(stamp).csv")
                @warn "Could not replace full-fit product file; writing timestamped fallback instead." product path fallback_path exception=(e, catch_backtrace())
                write_wide(fallback_path, times, ids, Y_corr)
                return fallback_path
        end
end


function write_validation_scope(
	outdir::AbstractString; cv_scheme::AbstractString, use_loocv_eval::Bool,
	fold_scheme::Union{Nothing,Symbol}=nothing,
	slope_ridge_candidates::Vector{Float64},
	min_scan_coverage::Float64,
	kernel_candidates::Vector{Int}=Int[],
)
	# The held-out geometry decides what the run can claim. `:random` interleaves validation
	# stations with training ones, so it supports "stations not used in fitting" but NOT
	# spatial generalisation to gauge-free ground; only `:spatial_block` holds out contiguous
	# area. Asserting the spatial claim under `:random` is a scientific overclaim, not a
	# naming slip, so the string is derived rather than fixed.
	supported_claim = if fold_scheme === :spatial_block
		"gauge-network-assisted spatial bias correction/interpolation for stations in held-out spatial blocks"
	elseif fold_scheme === :random
		"gauge-network-assisted bias correction/interpolation for stations not used in fitting, " *
			"held out at random and therefore interleaved with training stations; generalisation " *
			"to gauge-free regions is NOT supported by this split"
	else
		"gauge-network-assisted bias correction/interpolation for stations not used in fitting"
	end
	rows = [
		(key="validation_target", value="held-out station locations at matched observation/satellite timestamps"),
		(key="training_signal", value="same-timestamp observed-minus-satellite residuals from training stations"),
		(key="residual_definition", value="R = P_obs - P_sat; P_corr = P_sat + R_hat"),
		(key="local_model", value="target-centred local linear weighted ridge regression"),
		(key="local_coordinates", value="east/north offsets in kilometres, scaled by weighted RMS distance"),
		(key="slope_ridge_candidates", value=join(slope_ridge_candidates, ",")),
		(key="min_scan_coverage", value=string(min_scan_coverage)),
		(key="kernel_selection", value=length(kernel_candidates) > 1 ?
			"selected jointly with bandwidth/ridge via training-fold LOOCV, candidates=$(join(kernel_candidates, ","))" :
			"fixed by caller, not selected (kernel_candidates=$(join(kernel_candidates, ",")))"),
		(key="parameter_selection", value="coverage threshold, then LOOCV RMSE, MAE, slope_ridge, bandwidth, fixed before adaptive"),
		(key="negative_precipitation_policy", value="retain raw negative corrected values without clipping"),
		(key="validation_station_observations_used_for_fit", value="false"),
		(key="same_timestamp_training_station_observations_required", value="true"),
		(key="temporal_holdout", value="false"),
		(key="parameter_scan_use_loocv_eval", value=string(use_loocv_eval)),
		(key="cv_scheme", value=cv_scheme),
		(key="fold_scheme", value=fold_scheme === nothing ? "none" : string(fold_scheme)),
		(key="supported_claim", value=supported_claim),
		(key="unsupported_claim", value="standalone satellite correction or temporal forecast without concurrent training-station observations"),
	]
	CSV.write(joinpath(outdir, "validation_scope.csv"), DataFrame(rows))
end


"""
    split_stations_train_val(station_ids::Vector{String}; train_frac::Float64=0.8, rng::AbstractRNG=Random.GLOBAL_RNG)

Split stations into training and validation sets by a random shuffle.
Returns: (train_ids, val_ids)

The held-out stations are new/unseen *locations*, so this tests generalisation to stations that
were not fitted. It is NOT a spatial split: the shuffle never reads `lonlat`, so held-out
stations sit interleaved with training ones and a validation station typically has a training
neighbour a few km away. For held-out *area*, use `split_stations_spatial_block_kfold`.
"""

# 按站点进行空间训练/验证集划分 

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
	@assert any(eval_mask) "[$product] 验证集没有共同有效样本，无法比较校正前后"

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


function run_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, rng::AbstractRNG=Random.GLOBAL_RNG, fold_scheme::Symbol=:random,
)
	mkpath(cfg.outdir)

	st = load_station_meta(cfg.station_meta_path;
		station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
	products, common_product_ids, product_data = load_global_common_product_data(cfg)

	common_lonlat = build_X_lonlat(st, common_product_ids; station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
	folds = if fold_scheme == :random
		split_stations_kfold(common_product_ids; k=k, rng=rng)
	elseif fold_scheme == :spatial_block
		split_stations_spatial_block_kfold(common_product_ids, common_lonlat; k=k)
	else
		throw(ArgumentError("fold_scheme must be :random or :spatial_block"))
	end
	cv_scheme = fold_scheme == :spatial_block ? "station_spatial_block_$(k)fold" : "station_$(k)fold"
	write_validation_scope(
		cfg.outdir;
		cv_scheme="$(cv_scheme)_pooled",
		use_loocv_eval=cfg.use_loocv_eval,
		fold_scheme,
		slope_ridge_candidates=cfg.slope_ridge_candidates,
		min_scan_coverage=cfg.min_scan_coverage,
		kernel_candidates=cfg.kernels,
	)

	common_split = DataFrame(station_id=String[], fold=Int[], fold_scheme=String[])
	for (fold_idx, fold_ids) in enumerate(folds)
		for station_id in fold_ids
			push!(common_split, (station_id=station_id, fold=fold_idx, fold_scheme=string(fold_scheme)))
		end
	end
	CSV.write(joinpath(cfg.outdir, "split_common_spatial5fold.csv"), common_split)

	fold_summary_rows = NamedTuple[]
	pooled_summary_rows = NamedTuple[]
	for product in products
		data = product_data[product]
		times = data.times
		id_map_all = Dict(id => i for (i, id) in enumerate(data.ids))
		use_idx = [id_map_all[id] for id in common_product_ids]
		ids = common_product_ids
		Y_obs = data.Y_obs[use_idx, :]
		Y_sat = data.Y_sat[use_idx, :]
		lonlat = build_X_lonlat(st, ids; station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)
		Y_corr_pooled = fill(NaN, size(Y_obs))
		selected_slope_ridges = Float64[]

		product_split = DataFrame(station_id=String[], fold=Int[], fold_scheme=String[])
		for (fold_idx, fold_ids) in enumerate(folds)
			for station_id in fold_ids
				push!(product_split, (station_id=station_id, fold=fold_idx, fold_scheme=string(fold_scheme)))
			end
		end
		CSV.write(joinpath(cfg.outdir, "split_$(product)_spatial5fold.csv"), product_split)

		for fold_idx in 1:k
			val_ids = folds[fold_idx]
			train_ids = String[]
			for other_fold in 1:k
				other_fold == fold_idx && continue
				append!(train_ids, folds[other_fold])
			end

			fold_dir = joinpath(cfg.outdir, "fold_$(fold_idx)")
			mkpath(fold_dir)
			split_df = DataFrame(
				station_id = vcat(train_ids, val_ids),
				split = vcat(fill("train", length(train_ids)), fill("validation", length(val_ids))),
				fold = fill(fold_idx, length(train_ids) + length(val_ids)),
				fold_scheme = fill(string(fold_scheme), length(train_ids) + length(val_ids)),
			)
			CSV.write(joinpath(fold_dir, "split_$(product)_spatialcv.csv"), split_df)

			println("[$product] $cv_scheme fold $fold_idx/$k: $(length(train_ids)) train sites, $(length(val_ids)) val sites")
			result = evaluate_spatial_holdout(
				product, ids, times, Y_obs, Y_sat, lonlat, train_ids, val_ids, cfg;
				scan_path=joinpath(fold_dir, "scan_$(product)_spatialcv.csv"),
				fail_path=joinpath(fold_dir, "scan_$(product)_spatialcv_failures.csv"),
			)
			Y_corr_pooled[result.val_idx, :] = result.Y_corr_val
			push!(selected_slope_ridges, result.slope_ridge)
			write_wide(joinpath(fold_dir, "corr_$(product)_spatialcv_val.csv"), times, val_ids, result.Y_corr_val)

			push!(fold_summary_rows, (
				product=product,
				fold=fold_idx,
				spatial_cv=true,
				cv_scheme=cv_scheme,
				n_train=result.n_train,
				n_val=result.n_val,
				n=result.n,
				pre_n=result.pre_n,
				post_n=result.post_n,
				common_n=result.common_n,
				total_n=result.total_n,
				coverage=result.coverage,
				scan_n=result.scan_n,
				scan_coverage=result.scan_coverage,
				kernel=result.kernel,
				adaptive=result.adaptive,
				bw=result.bw,
				slope_ridge=result.slope_ridge,
				RMSE_pre=result.RMSE_pre, RMSE_post=result.RMSE_post,
				MAE_pre=result.MAE_pre, MAE_post=result.MAE_post,
				Bias_pre=result.Bias_pre, Bias_post=result.Bias_post,
				r_pre=result.r_pre, r_post=result.r_post,
				POD_pre=result.POD_pre, POD_post=result.POD_post,
				FAR_pre=result.FAR_pre, FAR_post=result.FAR_post,
				CSI_pre=result.CSI_pre, CSI_post=result.CSI_post,
				negative_n=result.negative_n,
				negative_fraction=result.negative_fraction,
				min_corrected=result.min_corrected,
			))
		end

		write_wide(joinpath(cfg.outdir, "corr_$(product)_spatial5fold_val.csv"), times, ids, Y_corr_pooled)
		pre_mask = common_valid_mask(Y_obs, Y_sat)
		post_mask = common_valid_mask(Y_obs, Y_corr_pooled)
		eval_mask = common_valid_mask(Y_obs, Y_sat, Y_corr_pooled)
		@assert any(eval_mask) "[$product] 五折验证没有共同有效样本，无法比较校正前后"

		pre_c = metric_continuous(Y_obs, Y_sat; mask=eval_mask)
		pre_e = metric_event(Y_obs, Y_sat; thr=cfg.rain_threshold, mask=eval_mask)
		post_c = metric_continuous(Y_obs, Y_corr_pooled; mask=eval_mask)
		post_e = metric_event(Y_obs, Y_corr_pooled; thr=cfg.rain_threshold, mask=eval_mask)
		negative = negative_output_stats(Y_corr_pooled)
		push!(pooled_summary_rows, (
			product=product,
			spatial_cv=true,
			cv_scheme="$(cv_scheme)_pooled",
			k=k,
			n_station=length(ids),
			n=post_c.n,
			pre_n=count(pre_mask),
			post_n=count(post_mask),
			common_n=count(eval_mask),
			total_n=length(eval_mask),
			coverage=count(eval_mask) / length(eval_mask),
			slope_ridge_by_fold=join(selected_slope_ridges, ","),
			RMSE_pre=pre_c.RMSE, RMSE_post=post_c.RMSE,
			MAE_pre=pre_c.MAE, MAE_post=post_c.MAE,
			Bias_pre=pre_c.Bias, Bias_post=post_c.Bias,
			r_pre=pre_c.r, r_post=post_c.r,
			POD_pre=pre_e.POD, POD_post=post_e.POD,
			FAR_pre=pre_e.FAR, FAR_post=post_e.FAR,
			CSI_pre=pre_e.CSI, CSI_post=post_e.CSI,
			negative...,
		))
	end

	fold_summary_df = DataFrame(fold_summary_rows)
	pooled_summary_df = DataFrame(pooled_summary_rows)
	fold_stats_df = summarize_fold_metrics(fold_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_folds.csv"), fold_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_fold_stats.csv"), fold_stats_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products_pooled.csv"), pooled_summary_df)
	CSV.write(joinpath(cfg.outdir, "summary_three_products.csv"), pooled_summary_df)
	return pooled_summary_df
end


function add_kernel_metadata!(df::DataFrame, kernel::Int, name::AbstractString)
	if :kernel in propertynames(df)
		df[!, :kernel] .= kernel
	else
		insertcols!(df, 2, :kernel => fill(kernel, nrow(df)))
	end
	if :kernel_name in propertynames(df)
		df[!, :kernel_name] .= String(name)
	else
		insertcols!(df, 3, :kernel_name => fill(String(name), nrow(df)))
	end
	return df
end


"""
    run_multikernel_spatial_kfold_pipeline(
        cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
    )

Run one independent spatial K-fold experiment per configured kernel. Each kernel
uses the same station split and selects its bandwidth from the configured grids.
"""
function run_multikernel_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
)
	isempty(cfg.kernels) && throw(ArgumentError("cfg.kernels must not be empty"))
	length(unique(cfg.kernels)) == length(cfg.kernels) ||
		throw(ArgumentError("cfg.kernels contains duplicate kernel values"))
	isempty(cfg.bw_adaptive) && isempty(cfg.bw_fixed_km) &&
		throw(ArgumentError("at least one adaptive or fixed bandwidth is required"))
	isempty(cfg.slope_ridge_candidates) &&
		throw(ArgumentError("cfg.slope_ridge_candidates must not be empty"))
	all(isfinite(ridge) && ridge >= 0 for ridge in cfg.slope_ridge_candidates) ||
		throw(ArgumentError("cfg.slope_ridge_candidates must be finite and non-negative"))
	isfinite(cfg.min_scan_coverage) && 0.0 <= cfg.min_scan_coverage <= 1.0 ||
		throw(ArgumentError("cfg.min_scan_coverage must be between 0 and 1"))
	for kernel in cfg.kernels
		kernel_name(kernel) # Validate before creating any output.
	end

	mkpath(cfg.outdir)
	status_df = DataFrame(
		kernel=Int[], kernel_name=String[], status=String[],
		error=String[], outdir=String[],
	)
	pooled_tables = DataFrame[]
	fold_tables = DataFrame[]
	fold_stats_tables = DataFrame[]

	for kernel in cfg.kernels
		name = kernel_name(kernel)
		kernel_outdir = joinpath(cfg.outdir, name)
		kernel_cfg = config_for_kernel(cfg, kernel, kernel_outdir)
		try
			pooled = run_spatial_kfold_pipeline(
				kernel_cfg;
				k=k,
				rng=MersenneTwister(seed),
				fold_scheme=fold_scheme,
			)
			add_kernel_metadata!(pooled, kernel, name)
			folds = CSV.read(joinpath(kernel_outdir, "summary_three_products_folds.csv"), DataFrame)
			fold_stats = CSV.read(joinpath(kernel_outdir, "summary_three_products_fold_stats.csv"), DataFrame)
			add_kernel_metadata!(folds, kernel, name)
			add_kernel_metadata!(fold_stats, kernel, name)
			push!(pooled_tables, pooled)
			push!(fold_tables, folds)
			push!(fold_stats_tables, fold_stats)
			push!(status_df, (kernel, name, "success", "", kernel_outdir))
		catch e
			message = sprint(showerror, e)
			push!(status_df, (kernel, name, "failed", message, kernel_outdir))
			@error "Kernel experiment failed" kernel name exception=(e, catch_backtrace())
		end
	end

	CSV.write(joinpath(cfg.outdir, "kernel_run_status.csv"), status_df)
	if !isempty(pooled_tables)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_pooled.csv"),
			vcat(pooled_tables...; cols=:union),
		)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_folds.csv"),
			vcat(fold_tables...; cols=:union),
		)
		CSV.write(
			joinpath(cfg.outdir, "summary_five_kernels_fold_stats.csv"),
			vcat(fold_stats_tables...; cols=:union),
		)
	end

	failed = status_df[status_df.status .== "failed", :]
	if nrow(failed) > 0
		error("$(nrow(failed)) kernel experiment(s) failed; see $(joinpath(cfg.outdir, "kernel_run_status.csv"))")
	end
	return vcat(pooled_tables...; cols=:union)
end

"""
Per-product count of how many of the k folds picked each kernel, from a
`summary_three_products_folds.csv`-shaped table (must carry `product`, `fold`, `kernel` columns).
A kernel winning every fold means the choice is stable; a split vote means it is fold-dependent.
"""
function _kernel_selection_stability(folds::DataFrame)
	nrow(folds) == 0 && return DataFrame()
	rows = NamedTuple[]
	for group in groupby(folds, [:product, :kernel])
		push!(rows, (;
			product=String(group.product[1]), kernel=Int(group.kernel[1]),
			kernel_name=kernel_name(Int(group.kernel[1])), win_count=nrow(group),
		))
	end
	stability = DataFrame(rows)
	fold_counts = combine(groupby(stability, :product), :win_count => sum => :fold_count)
	stability = leftjoin(stability, fold_counts; on=:product)
	stability.selection_frequency = stability.win_count ./ stability.fold_count
	return sort!(stability, [:product, order(:win_count; rev=true)])
end

"""
    run_nested_kernel_spatial_kfold_pipeline(
        cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
    )

Select the kernel via nested cross-validation: every training fold picks its own best
`(kernel, bw, adaptive, slope_ridge)` by that fold's own LOOCV, using only training-fold data,
and the held-out fold stations validate the winning combination — reusing `scan_params`'s
existing joint scan over `cfg.kernels` unrestricted, rather than forcing one kernel per run.

Contrast with `run_multikernel_spatial_kfold_pipeline`, which forces one kernel per run (via
`config_for_kernel`) to give each kernel its own paired head-to-head comparison on the same
fold split; its pooled output answers "how does each kernel perform," not "which kernel should
be used," and reading it as a selection recommendation is exactly the leak this function closes.
"""
function run_nested_kernel_spatial_kfold_pipeline(
	cfg::MGERConfig; k::Int=5, seed::Integer, fold_scheme::Symbol=:random,
)
	length(cfg.kernels) > 1 || throw(ArgumentError(
		"run_nested_kernel_spatial_kfold_pipeline needs more than one kernel in cfg.kernels; " *
		"for a single fixed kernel use run_spatial_kfold_pipeline directly",
	))
	pooled = run_spatial_kfold_pipeline(cfg; k, rng=MersenneTwister(seed), fold_scheme)
	folds = CSV.read(joinpath(cfg.outdir, "summary_three_products_folds.csv"), DataFrame)
	stability = _kernel_selection_stability(folds)
	CSV.write(joinpath(cfg.outdir, "kernel_selection_stability.csv"), stability)
	return pooled
end


"""
    run_pipeline(cfg::MGERConfig; spatial_cv::Bool=true, train_frac::Float64=0.8, rng::AbstractRNG=Random.GLOBAL_RNG)

Run MGER_FINAL pipeline with an optional held-out station split.

# Arguments
- cfg: MGERConfig with pipeline parameters
- spatial_cv: If true, hold out a random `1 - train_frac` share of stations (see
  `split_stations_train_val` - a random station holdout, not a spatial one)
- train_frac: Fraction of stations to use for training (default 0.8)
- rng: Random number generator for reproducible splits

# Returns
- summary_df: DataFrame with evaluation metrics
- If spatial_cv=true, also returns train/val station lists for each product
"""
function run_pipeline(cfg::MGERConfig; spatial_cv::Bool=true, train_frac::Float64=0.8, rng::AbstractRNG=Random.GLOBAL_RNG)
	mkpath(cfg.outdir)
	if !spatial_cv && !cfg.use_loocv_eval
		error("Refusing spatial_cv=false with use_loocv_eval=false: this would report full-fit in-sample metrics as if they were validation results. Run spatial_cv=true for independent evaluation; that path also writes corr_*_fullfit_insample.csv product files.")
	end
	# `spatial_cv=true` is a random 80/20 station holdout - `split_stations_train_val` shuffles
	# the id list and never reads coordinates - so the label must not say "spatial".
	cv_scheme = spatial_cv ? "random_station_holdout" : "loocv_eval_all_stations"
	write_validation_scope(
		cfg.outdir;
		cv_scheme=cv_scheme,
		use_loocv_eval=cfg.use_loocv_eval,
		fold_scheme=nothing,
		slope_ridge_candidates=cfg.slope_ridge_candidates,
		min_scan_coverage=cfg.min_scan_coverage,
		kernel_candidates=cfg.kernels,
	)

	st = load_station_meta(cfg.station_meta_path;
		station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)

	products, common_product_ids, product_data = load_global_common_product_data(cfg)

	common_train_ids = String[]
	common_val_ids = String[]
	if spatial_cv
		common_train_ids, common_val_ids = split_stations_train_val(common_product_ids; train_frac=train_frac, rng=rng)
		split_df = DataFrame(
			station_id = vcat(common_train_ids, common_val_ids),
			split = vcat(fill("train", length(common_train_ids)), fill("validation", length(common_val_ids))),
		)
		CSV.write(joinpath(cfg.outdir, "split_common_spatialcv.csv"), split_df)
	end

	summary_rows = NamedTuple[]
	for product in products
		data = product_data[product]
		times = data.times
		id_map_all = Dict(id => i for (i, id) in enumerate(data.ids))
		use_idx = [id_map_all[id] for id in common_product_ids]
		ids = common_product_ids
		Y_obs = data.Y_obs[use_idx, :]
		Y_sat = data.Y_sat[use_idx, :]
		lonlat = build_X_lonlat(st, ids; station_id_col=cfg.station_id_col, lon_col=cfg.lon_col, lat_col=cfg.lat_col)

		if spatial_cv
			train_ids = common_train_ids
			val_ids = common_val_ids
			split_df = DataFrame(
				station_id = vcat(train_ids, val_ids),
				split = vcat(fill("train", length(train_ids)), fill("validation", length(val_ids))),
			)
			CSV.write(joinpath(cfg.outdir, "split_$(product)_spatialcv.csv"), split_df)

			# Get indices for train and val
			id_map = Dict(id => i for (i, id) in enumerate(ids))
			train_idx = [id_map[id] for id in train_ids]
			val_idx = [id_map[id] for id in val_ids]

			# Build train/val matrices
			lonlat_train = lonlat[train_idx, :]
			lonlat_val = lonlat[val_idx, :]
			Y_obs_train = Y_obs[train_idx, :]
			Y_sat_train = Y_sat[train_idx, :]
			Y_obs_val = Y_obs[val_idx, :]
			Y_sat_val = Y_sat[val_idx, :]

			# Distance matrix among training sites (for fitting)
			dMat_train = pairwise_haversine_km(lonlat_train)
			# Distance matrix from train to val sites (for prediction)
			dMat_train_val = let
				points_train = map(x -> (x[1], x[2]), eachrow(lonlat_train))
				points_val = map(x -> (x[1], x[2]), eachrow(lonlat_val))
				pairwise(Haversine(6378.388), points_train, points_val)
			end

			println("[$product] $cv_scheme: $(length(train_ids)) train sites, $(length(val_ids)) val sites")

			# Parameter scan on training data only
			scan_df, best = scan_params(
				lonlat_train, Y_obs_train, Y_sat_train, dMat_train;
				kernels=cfg.kernels,
				bw_adaptive=cfg.bw_adaptive,
				bw_fixed_km=cfg.bw_fixed_km,
				slope_ridge_candidates=cfg.slope_ridge_candidates,
				min_scan_coverage=cfg.min_scan_coverage,
				rain_threshold=cfg.rain_threshold,
				use_loocv=cfg.use_loocv_eval,
				fail_path=joinpath(cfg.outdir, "scan_$(product)_spatialcv_failures.csv"),
			)
			CSV.write(joinpath(cfg.outdir, "scan_$(product)_spatialcv.csv"), scan_df)

			# Train on train sites and predict validation residuals with the
			# target-centred local linear model.
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
			@assert any(eval_mask) "[$product] 验证集没有共同有效样本，无法比较校正前后"

			pre_c = metric_continuous(Y_obs_val, Y_sat_val; mask=eval_mask)
			pre_e = metric_event(Y_obs_val, Y_sat_val; thr=cfg.rain_threshold, mask=eval_mask)
			post_c = metric_continuous(Y_obs_val, Y_corr_val; mask=eval_mask)
			post_e = metric_event(Y_obs_val, Y_corr_val; thr=cfg.rain_threshold, mask=eval_mask)

			# Save corrected data for validation sites
			write_wide(joinpath(cfg.outdir, "corr_$(product)_spatialcv_val.csv"), times, val_ids, Y_corr_val)

			# Product output: fit on all common stations and write a clearly labeled
			# in-sample file. This file is not used for validation metrics.
			dMat_full = pairwise_haversine_km(lonlat)
			fullfit_path = write_fullfit_product(
				cfg.outdir, product, times, ids, lonlat, Y_obs, Y_sat, dMat_full, best,
			)
			println("[$product] Full-fit product saved for product use only: $fullfit_path")

			push!(summary_rows, (
				product=product,
				spatial_cv=true,
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
			))
		else
			# Standard pipeline: all sites
			dMat = pairwise_haversine_km(lonlat)

			scan_df, best = scan_params(
				lonlat, Y_obs, Y_sat, dMat;
				kernels=cfg.kernels,
				bw_adaptive=cfg.bw_adaptive,
				bw_fixed_km=cfg.bw_fixed_km,
				slope_ridge_candidates=cfg.slope_ridge_candidates,
				min_scan_coverage=cfg.min_scan_coverage,
				rain_threshold=cfg.rain_threshold,
				use_loocv=cfg.use_loocv_eval,
				fail_path=joinpath(cfg.outdir, "scan_$(product)_failures.csv"),
			)
			CSV.write(joinpath(cfg.outdir, "scan_$(product).csv"), scan_df)

			dist = cfg.use_loocv_eval ? make_loocv_dist(dMat) : dMat
			wMat = gw_weight(dist, Float64(best.bw); kernel=Int(best.kernel), adaptive=Bool(best.adaptive))
			R = Y_obs .- Y_sat
			Rhat = local_linear_residual_predict(
				lonlat, R, lonlat, wMat; slope_ridge=Float64(best.slope_ridge),
			)
			Y_corr = Y_sat .+ Rhat
			negative = negative_output_stats(Y_corr)

			pre_mask = common_valid_mask(Y_obs, Y_sat)
			post_mask = common_valid_mask(Y_obs, Y_corr)
			eval_mask = common_valid_mask(Y_obs, Y_sat, Y_corr)
			@assert any(eval_mask) "[$product] 没有共同有效样本，无法比较校正前后"

			pre_c = metric_continuous(Y_obs, Y_sat; mask=eval_mask)
			pre_e = metric_event(Y_obs, Y_sat; thr=cfg.rain_threshold, mask=eval_mask)
			post_c = metric_continuous(Y_obs, Y_corr; mask=eval_mask)
			post_e = metric_event(Y_obs, Y_corr; thr=cfg.rain_threshold, mask=eval_mask)

			write_wide(joinpath(cfg.outdir, "corr_$(product)_loocv_eval.csv"), times, ids, Y_corr)

			push!(summary_rows, (
				product=product,
				spatial_cv=false,
				n_train=length(ids),
				n_val=missing,
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
			))
		end
	end

	summary_df = DataFrame(summary_rows)
	CSV.write(joinpath(cfg.outdir, "summary_three_products.csv"), summary_df)
	return summary_df
end


"""
Main audited entry points:

    julia --project=. scripts/run_mger_three_products_timealigned.jl
    julia --project=. scripts/run_mger_three_products_timealigned_loocv.jl
"""
