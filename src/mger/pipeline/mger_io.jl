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
	error("cannot parse time: $s")
end


"""Read a wide table (rows = times, columns = stations) and return Y[station, time]."""
function read_hourly_wide(path::AbstractString; time_col::Symbol=:time)
	df = CSV.read(path, DataFrame)
	names_lower = Dict(Symbol(lowercase(String(c))) => c for c in names(df))
	c_time = get(names_lower, Symbol(lowercase(String(time_col))), nothing)
	@assert c_time !== nothing "missing time column: $time_col"

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
	@assert size(Y, 2) == length(times) "time axis length does not match the number of data columns"
	if analysis_start !== nothing && analysis_end !== nothing
		@assert analysis_start <= analysis_end "analysis start must not be later than analysis end"
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
	@assert haskey(sat_inputs, reference_product) "missing reference product for global time QC: $reference_product"
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
	@assert length(unique(obs_times)) == length(obs_times) "observations contain duplicate timestamps"

	products = sort(collect(keys(cfg.sat_paths)))
	sat_inputs = Dict{String, Any}()
	common_product_ids = copy(obs_ids)
	for product in products
		sat_times, sat_ids, Y_sat0 = read_hourly_wide(cfg.sat_paths[product]; time_col=cfg.time_col)
		sat_times, Y_sat0 = subset_time_window(sat_times, Y_sat0, cfg.analysis_start, cfg.analysis_end)
		@assert length(unique(sat_times)) == length(sat_times) "[$product] contains duplicate timestamps"
		id_set = Set(sat_ids)
		common_product_ids = [id for id in common_product_ids if id in id_set]
		sat_inputs[product] = (; times=sat_times, ids=sat_ids, Y=Y_sat0)
	end
	@assert !isempty(common_product_ids) "the three products share no stations; a fair cross-product evaluation is impossible"

	time_qc_path = joinpath(cfg.outdir, "global_common_time_qc.csv")
	_, common_product_times = write_global_time_qc(time_qc_path, obs_times, sat_inputs)
	if cfg.expected_common_time_count !== nothing && length(common_product_times) != cfg.expected_common_time_count
		error("global common time count $(length(common_product_times)) does not match the expected $(cfg.expected_common_time_count). The difference list was written to $time_qc_path")
	end
	@assert !isempty(common_product_times) "the three products and the observations share no common times"

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

	@assert c_sid !== nothing "missing station id column: $station_id_col"
	@assert c_lon !== nothing "missing longitude column: $lon_col"
	@assert c_lat !== nothing "missing latitude column: $lat_col"

	rename!(st, c_sid => :station_id, c_lon => :lon, c_lat => :lat)
	st[!, :station_id] = string.(st[!, :station_id])
	return st
end


"""Build the lon/lat coordinate matrix in `common_ids` order, for distance calculations."""
function build_X_lonlat(st::DataFrame, common_ids::Vector{String};
	station_id_col::Symbol=:station_id, lon_col::Symbol=:lon, lat_col::Symbol=:lat)
	st_map = Dict(string(st[i, station_id_col]) => i for i in 1:nrow(st))
	idx = [st_map[sid] for sid in common_ids if haskey(st_map, sid)]
	@assert length(idx) == length(common_ids) "station metadata is incomplete; some station ids are missing"
	lon = Float64.(st[idx, lon_col])
	lat = Float64.(st[idx, lat_col])
	X = hcat(lon, lat)
	return X
end


"""Build the ST-GWR design matrix X=[1, lon_centered, lat_centered]."""
function build_X_intercept_centered(lonlat::Matrix{Float64};
	center::Tuple{Float64,Float64}=(mean(lonlat[:, 1]), mean(lonlat[:, 2])))
	lon = lonlat[:, 1]
	lat = lonlat[:, 2]
	lon_center, lat_center = center
	return hcat(ones(Float64, size(lonlat, 1)), lon .- lon_center, lat .- lat_center)
end
