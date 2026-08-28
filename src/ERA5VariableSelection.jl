module ERA5VariableSelection

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

export ERA5SelectionConfig, ERA5_VARIABLES, align_feature_time, load_era5_panel
export balanced_spatial_folds, prepare_dynamic_panel, station_block_permutation
export dynamic_panel_screen, panel_spatial_variability_test, run_era5_variable_selection

const ERA5_VARIABLES = (:t2m_c, :d2m_c, :u10, :v10, :sp_hpa)

Base.@kwdef struct ERA5SelectionConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    min_stations_per_time::Int = 12
    k::Int = 5
    seed::Int = 20260816
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    association_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
    time_offset_hours::Int = 9
end

align_feature_time(time_utc::DateTime; offset_hours::Int=9) = time_utc + Hour(offset_hours)

function _parse_datetime(value)
    value isa DateTime && return value
    text = String(value)
    try
        return DateTime(text)
    catch
        for format in (dateformat"yyyy-mm-dd HH:MM:SS", dateformat"yyyy-mm-ddTHH:MM:SS")
            try
                return DateTime(text, format)
            catch
            end
        end
    end
    throw(ArgumentError("invalid ERA5 timestamp: $text"))
end

function _finite_float(value)
    ismissing(value) && return NaN
    value isa Real && return Float64(value)
    parsed = tryparse(Float64, String(value))
    return parsed === nothing ? NaN : parsed
end

function _csv_field(text::AbstractString)
    value = strip(text)
    if ncodeunits(value) >= 2 && startswith(value, '"') && endswith(value, '"')
        return replace(value[2:end-1], "\"\"" => "\"")
    end
    return value
end

function _bh_adjust(pvalues::AbstractVector{<:Real})
    p = Float64.(pvalues)
    m = length(p)
    m == 0 && return Float64[]
    order = sortperm(p)
    adjusted = fill(NaN, m)
    running = 1.0
    for rank in m:-1:1
        index = order[rank]
        running = min(running, p[index] * m / rank)
        adjusted[index] = min(running, 1.0)
    end
    return adjusted
end

"""Load requested station-hour ERA5 values using `time_utc + offset`, never `time_bjt`.

The function streams annual files and returns station × feature-time matrices plus a
quality-control table. Missing, duplicate, non-finite, or incorrectly aligned records
are reported and cause an error by default; no temporal filling is performed.
"""
function load_era5_panel(
    annual_paths::AbstractDict{<:Integer,<:AbstractString}, station_ids::Vector{String},
    feature_times::Vector{DateTime}; variables=collect(ERA5_VARIABLES),
    offset_hours::Int=9, strict::Bool=true, validate_annual_complete::Bool=true,
)
    length(unique(station_ids)) == length(station_ids) ||
        throw(ArgumentError("station IDs must be unique"))
    length(unique(feature_times)) == length(feature_times) ||
        throw(ArgumentError("feature times must be unique"))
    vars = Symbol.(variables)
    all(v -> v in ERA5_VARIABLES, vars) || throw(ArgumentError("unsupported ERA5 variable"))
    id_index = Dict(id => i for (i, id) in enumerate(station_ids))
    time_index = Dict(time => j for (j, time) in enumerate(feature_times))
    matrices = Dict(v => fill(NaN, length(station_ids), length(feature_times)) for v in vars)
    seen = falses(length(station_ids), length(feature_times))
    duplicates = 0
    offset_mismatch = 0
    matched = 0
    annual_expected = 0
    annual_matched = 0
    annual_duplicates = 0
    annual_nonfinite = 0
    annual_invalid_times = 0
    annual_unknown_stations = 0
    needed = vcat([:station_id, :time_utc, :time_bjt], vars)
    target_years = Set(year(t - Hour(offset_hours)) for t in feature_times)
    missing_years = sort(collect(setdiff(target_years, Set(Int.(keys(annual_paths))))))
    isempty(missing_years) || throw(ArgumentError("missing ERA5 annual paths: $(join(missing_years, ", "))"))

    for yr in sort(collect(target_years))
        path = String(annual_paths[yr])
        isfile(path) || throw(ArgumentError("ERA5 annual file does not exist: $path"))
        year_start = DateTime(yr, 1, 1)
        year_end = DateTime(yr + 1, 1, 1)
        hours_in_year = Int(Dates.value(year_end - year_start) ÷ 3_600_000)
        year_seen = falses(length(station_ids), hours_in_year)
        annual_expected += length(station_ids) * hours_in_year
        open(path, "r") do io
            eof(io) && throw(ArgumentError("empty ERA5 annual file: $path"))
            header = Symbol.(_csv_field.(split(chomp(readline(io)), ','; keepempty=true)))
            positions = Dict(name => findfirst(==(name), header) for name in needed)
            missing_columns = [name for (name, position) in positions if position === nothing]
            isempty(missing_columns) || throw(ArgumentError(
                "ERA5 file is missing columns: $(join(missing_columns, ", "))",
            ))
            for line in eachline(io)
                fields = split(chomp(line), ','; keepempty=true)
                sid = _csv_field(fields[positions[:station_id]])
                if !haskey(id_index, sid)
                    annual_unknown_stations += 1
                    continue
                end
                utc = _parse_datetime(_csv_field(fields[positions[:time_utc]]))
                slot = Int(Dates.value(utc - year_start) ÷ 3_600_000) + 1
                exact_hour = year_start <= utc < year_end &&
                    utc == year_start + Hour(slot - 1) && 1 <= slot <= hours_in_year
                if !exact_hour
                    annual_invalid_times += 1
                    continue
                end
                i = id_index[sid]
                if year_seen[i, slot]
                    annual_duplicates += 1
                else
                    year_seen[i, slot] = true
                    annual_matched += 1
                end
                bjt = _parse_datetime(_csv_field(fields[positions[:time_bjt]]))
                bjt == utc + Hour(8) || (offset_mismatch += 1)
                feature_time = align_feature_time(utc; offset_hours=offset_hours)
                is_target = haskey(time_index, feature_time)
                j = is_target ? time_index[feature_time] : 0
                target_duplicate = is_target && seen[i, j]
                if target_duplicate
                    duplicates += 1
                elseif is_target
                    seen[i, j] = true
                    matched += 1
                end
                for variable in vars
                    value = _finite_float(_csv_field(fields[positions[variable]]))
                    !isfinite(value) && (annual_nonfinite += 1)
                    is_target && !target_duplicate && (matrices[variable][i, j] = value)
                end
            end
        end
    end
    missing_cells = count(!, seen)
    nonfinite_cells = sum(count(x -> !isfinite(x), matrix) for matrix in values(matrices))
    qc = DataFrame(
        station_count=[length(station_ids)], target_time_count=[length(feature_times)],
        expected_station_hours=[length(station_ids) * length(feature_times)],
        matched_station_hours=[matched], duplicate_keys=[duplicates],
        missing_station_hours=[missing_cells], nonfinite_values=[nonfinite_cells],
        annual_expected_station_hours=[annual_expected], annual_matched_station_hours=[annual_matched],
        annual_duplicate_keys=[annual_duplicates],
        annual_missing_station_hours=[annual_expected - annual_matched],
        annual_nonfinite_values=[annual_nonfinite], annual_invalid_times=[annual_invalid_times],
        annual_unknown_station_rows=[annual_unknown_stations],
        time_bjt_offset_mismatches=[offset_mismatch], feature_time_offset_hours=[offset_hours],
        annual_complete=[annual_duplicates == 0 && annual_expected == annual_matched &&
            annual_nonfinite == 0 && annual_invalid_times == 0 && annual_unknown_stations == 0],
        complete=[duplicates == 0 && missing_cells == 0 && nonfinite_cells == 0 &&
            (!validate_annual_complete || (annual_duplicates == 0 && annual_expected == annual_matched &&
                annual_nonfinite == 0 && annual_invalid_times == 0 && annual_unknown_stations == 0)) &&
            offset_mismatch == 0],
    )
    if strict && !qc.complete[1]
        error("ERA5 station-hour quality control failed: $(NamedTuple(qc[1, :]))")
    end
    return (; values=matrices, qc)
end

# One of four different spatial-fold builders in the project, and the simplest: a sort-and-slice
# along the longer axis. The others are `split_stations_balanced_spatial_kfold`
# (InterpolationBenchmarkFolds.jl - Hilbert-seeded capacity-constrained Lloyd),
# `DEMTerrainExperiment._balanced_spatial_folds` (farthest-point-seeded Lloyd) and
# `MGERPipeline.split_stations_spatial_block_kfold` (contiguous lon/lat blocks). They are genuinely
# different partitioning algorithms rather than copies, so they are not merged; but which one a
# given selection path uses is historical, not reasoned, and that is worth knowing before reading
# across their results.
function balanced_spatial_folds(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("lonlat must have n × 2 rows"))
    lat_mid = mean(lonlat[:, 2])
    lon_span = (maximum(lonlat[:, 1]) - minimum(lonlat[:, 1])) * cosd(lat_mid) * 111.32
    lat_span = (maximum(lonlat[:, 2]) - minimum(lonlat[:, 2])) * 110.57
    column = lon_span >= lat_span ? 1 : 2
    order = sortperm(1:n; by=i -> (lonlat[i, column], lonlat[i, 3 - column], station_ids[i]))
    assignment = zeros(Int, n)
    base, extra = divrem(n, k)
    left = 1
    for fold in 1:k
        count_in_fold = base + (fold <= extra ? 1 : 0)
        indices = order[left:(left + count_in_fold - 1)]
        assignment[indices] .= fold
        left += count_in_fold
    end
    return assignment
end

function _weighted_mean(x::Vector{Float64}, w::Vector{Float64})
    total = sum(w)
    total > 0 || return NaN
    return dot(x, w) / total
end

function _weighted_correlation(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64})
    length(x) >= 3 || return NaN
    mx, my = _weighted_mean(x, w), _weighted_mean(y, w)
    dx, dy = x .- mx, y .- my
    vx, vy = dot(w, dx .* dx), dot(w, dy .* dy)
    vx > 0 && vy > 0 || return NaN
    return dot(w, dx .* dy) / sqrt(vx * vy)
end

function _tie_ranks(values::Vector{Float64})
    order = sortperm(values)
    ranks = zeros(Float64, length(values))
    left = 1
    while left <= length(order)
        right = left
        while right < length(order) && values[order[right + 1]] == values[order[left]]
            right += 1
        end
        ranks[order[left:right]] .= (left + right) / 2
        left = right + 1
    end
    return ranks
end

function _flatten_balanced(x::Matrix{Float64}, y::Matrix{Float64})
    xs, ys, ws = Float64[], Float64[], Float64[]
    for i in axes(y, 1)
        valid = findall(j -> isfinite(x[i, j]) && isfinite(y[i, j]), axes(y, 2))
        isempty(valid) && continue
        weight = 1 / length(valid)
        for j in valid
            push!(xs, x[i, j]); push!(ys, y[i, j]); push!(ws, weight)
        end
    end
    return xs, ys, ws
end

"""Return a row permutation; indexing a panel with it preserves every station block."""
station_block_permutation(nstations::Int, rng::AbstractRNG) = randperm(rng, nstations)

function _station_pair_moments(x::Matrix{Float64}, y::Matrix{Float64})
    n = size(y, 1)
    size(x) == size(y) || throw(DimensionMismatch("panel matrices must match"))
    valid_x = [findall(isfinite, view(x, i, :)) for i in 1:n]
    valid_y = [findall(isfinite, view(y, i, :)) for i in 1:n]
    count_pair = zeros(Int, n, n)
    sum_x = zeros(n, n); sum_y = zeros(n, n)
    sum_xx = zeros(n, n); sum_yy = zeros(n, n); sum_xy = zeros(n, n)
    for source in 1:n, target in 1:n
        xs, ys = valid_x[source], valid_y[target]
        ix = 1; iy = 1
        while ix <= length(xs) && iy <= length(ys)
            tx, ty = xs[ix], ys[iy]
            if tx == ty
                xv, yv = x[source, tx], y[target, ty]
                count_pair[source, target] += 1
                sum_x[source, target] += xv; sum_y[source, target] += yv
                sum_xx[source, target] += xv^2; sum_yy[source, target] += yv^2
                sum_xy[source, target] += xv * yv
                ix += 1; iy += 1
            elseif tx < ty
                ix += 1
            else
                iy += 1
            end
        end
    end
    return (; count=count_pair, sum_x, sum_y, sum_xx, sum_yy, sum_xy)
end

function _correlation_from_pair_moments(moments, order::Vector{Int})
    total_x = 0.0; total_y = 0.0; total_xx = 0.0
    total_yy = 0.0; total_xy = 0.0; station_count = 0
    for target in eachindex(order)
        source = order[target]
        count_pair = moments.count[source, target]
        count_pair > 0 || continue
        scale = 1 / count_pair
        total_x += moments.sum_x[source, target] * scale
        total_y += moments.sum_y[source, target] * scale
        total_xx += moments.sum_xx[source, target] * scale
        total_yy += moments.sum_yy[source, target] * scale
        total_xy += moments.sum_xy[source, target] * scale
        station_count += 1
    end
    station_count >= 3 || return NaN
    mean_x, mean_y = total_x / station_count, total_y / station_count
    variance_x = total_xx / station_count - mean_x^2
    variance_y = total_yy / station_count - mean_y^2
    variance_x > 0 && variance_y > 0 || return NaN
    covariance = total_xy / station_count - mean_x * mean_y
    return covariance / sqrt(variance_x * variance_y)
end

function prepare_dynamic_panel(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, era5::AbstractDict,
    train_indices::Vector{Int}, times::Vector{DateTime}, cfg::ERA5SelectionConfig,
)
    size(Yobs) == size(Ysat) || throw(DimensionMismatch("observation and satellite matrices differ"))
    size(Yobs, 2) == length(times) || throw(DimensionMismatch("time dimension differs"))
    variables = collect(ERA5_VARIABLES)
    n, nt = length(train_indices), length(times)
    residual = Yobs[train_indices, :] .- Ysat[train_indices, :]
    raw = Dict(v => Float64.(era5[v][train_indices, :]) for v in variables)
    mask = falses(n, nt)
    for i in 1:n, t in 1:nt
        mask[i, t] = isfinite(Yobs[train_indices[i], t]) &&
            isfinite(Ysat[train_indices[i], t]) && Yobs[train_indices[i], t] >= cfg.wet_threshold &&
            all(isfinite(raw[v][i, t]) for v in variables)
    end
    previous_count = -1
    while count(mask) != previous_count
        previous_count = count(mask)
        eligible_station_now = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
        mask[.!eligible_station_now, :] .= false
        eligible_time_now = vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time
        mask[:, .!eligible_time_now] .= false
    end
    eligible_station = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours

    y = fill(NaN, n, nt)
    centered = Dict(v => fill(NaN, n, nt) for v in variables)
    for t in 1:nt
        indices = findall(mask[:, t])
        length(indices) >= cfg.min_stations_per_time || continue
        ymean = mean(residual[indices, t])
        y[indices, t] .= residual[indices, t] .- ymean
        for v in variables
            xmean = mean(raw[v][indices, t])
            centered[v][indices, t] .= raw[v][indices, t] .- xmean
        end
    end
    scales = DataFrame(variable=String[], mean=Float64[], scale=Float64[])
    for v in variables
        x, _, w = _flatten_balanced(centered[v], y)
        μ = _weighted_mean(x, w)
        σ = sqrt(_weighted_mean((x .- μ) .^ 2, w))
        isfinite(σ) && σ > 0 || (σ = NaN)
        if isfinite(σ)
            for index in eachindex(centered[v])
                isfinite(centered[v][index]) && (centered[v][index] = (centered[v][index] - μ) / σ)
            end
        end
        push!(scales, (String(v), μ, σ))
    end
    qc = DataFrame(
        training_station_count=[n], eligible_station_count=[count(eligible_station)],
        eligible_time_count=[count(vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time)],
        wet_station_hours=[count(mask)], min_wet_hours=[cfg.min_wet_hours],
        min_stations_per_time=[cfg.min_stations_per_time],
        status=[count(eligible_station) >= cfg.min_stations_per_time ? "ok" : "insufficient_stations"],
    )
    return (; y, x=centered, mask, eligible_station, scales, qc)
end

function _association(panel, variable::Symbol, permutations::Int, rng::AbstractRNG)
    x, y, w = _flatten_balanced(panel.x[variable], panel.y)
    moments = _station_pair_moments(panel.x[variable], panel.y)
    pearson = _correlation_from_pair_moments(moments, collect(1:size(panel.y, 1)))
    spearman = _weighted_correlation(_tie_ranks(x), _tie_ranks(y), w)
    if permutations <= 0 || !isfinite(pearson)
        return pearson, spearman, isfinite(pearson) ? NaN : 1.0
    end
    exceed = 0
    n = size(panel.y, 1)
    for _ in 1:permutations
        order = station_block_permutation(n, rng)
        value = _correlation_from_pair_moments(moments, order)
        exceed += isfinite(value) && abs(value) >= abs(pearson)
    end
    return pearson, spearman, (exceed + 1) / (permutations + 1)
end

function _weighted_vifs(panel, active::Vector{Symbol})
    length(active) <= 1 && return fill(1.0, length(active))
    rows = Vector{Vector{Float64}}()
    weights = Float64[]
    for i in axes(panel.y, 1)
        valid = [t for t in axes(panel.y, 2) if isfinite(panel.y[i, t]) &&
            all(isfinite(panel.x[v][i, t]) for v in active)]
        isempty(valid) && continue
        weight = 1 / length(valid)
        for t in valid
            push!(rows, [panel.x[v][i, t] for v in active]); push!(weights, weight)
        end
    end
    isempty(rows) && return fill(Inf, length(active))
    X = reduce(vcat, permutedims.(rows))
    sw = sqrt.(weights ./ sum(weights))
    result = fill(Inf, length(active))
    for j in eachindex(active)
        other = setdiff(eachindex(active), [j])
        design = hcat(ones(size(X, 1)), X[:, other])
        beta = (design .* sw) \ (X[:, j] .* sw)
        prediction = design * beta
        μ = _weighted_mean(X[:, j], weights)
        rss = sum(weights .* (X[:, j] .- prediction) .^ 2)
        tss = sum(weights .* (X[:, j] .- μ) .^ 2)
        r2 = tss > 0 ? clamp(1 - rss / tss, 0.0, 1.0) : 1.0
        result[j] = r2 < 1 ? 1 / (1 - r2) : Inf
    end
    return result
end

function dynamic_panel_screen(
    panel, times::Vector{DateTime}, cfg::ERA5SelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed),
)
    association = DataFrame(
        variable=String[], pearson=Float64[], spearman=Float64[],
        direction_stable=Bool[], pvalue=Float64[], qvalue=Float64[], selected_bh=Bool[],
    )
    raw = NamedTuple[]
    for variable in ERA5_VARIABLES
        pearson, spearman, pvalue = _association(
            panel, variable, cfg.association_permutations, rng,
        )
        stable = isfinite(pearson) && isfinite(spearman) && sign(pearson) == sign(spearman)
        push!(raw, (; variable, pearson, spearman, pvalue, stable))
    end
    qvalues = _bh_adjust([row.pvalue for row in raw])
    for (row, qvalue) in zip(raw, qvalues)
        push!(association, (
            String(row.variable), row.pearson, row.spearman, row.stable,
            row.pvalue, qvalue, row.stable && qvalue < cfg.q_threshold,
        ))
    end
    active = Symbol.(association.variable[association.selected_bh])
    vif = DataFrame(iteration=Int[], variable=String[], vif=Float64[], removed=Bool[], reason=String[])
    iteration = 1
    while !isempty(active)
        values = _weighted_vifs(panel, active)
        maxvif = maximum(values)
        if maxvif < cfg.vif_threshold
            for (variable, value) in zip(active, values)
                push!(vif, (iteration, String(variable), value, false, "retained"))
            end
            break
        end
        tied = findall(v -> isapprox(v, maxvif; rtol=1e-10, atol=1e-10) || (isinf(v) && isinf(maxvif)), values)
        qmap = Dict(Symbol(row.variable) => row.qvalue for row in eachrow(association))
        remove_index = sort(tied; by=j -> (qmap[active[j]], findfirst(==(active[j]), ERA5_VARIABLES)), rev=true)[1]
        for (j, (variable, value)) in enumerate(zip(active, values))
            push!(vif, (iteration, String(variable), value, j == remove_index,
                j == remove_index ? "vif_ge_threshold" : "pending"))
        end
        deleteat!(active, remove_index)
        iteration += 1
    end
    monthly = DataFrame(variable=String[], month=Int[], pearson=Float64[], spearman=Float64[],
        pearson_direction=String[], spearman_direction=String[])
    for variable in ERA5_VARIABLES, month_value in 6:9
        columns = findall(t -> month(t) == month_value, times)
        if isempty(columns)
            p, s = NaN, NaN
        else
            x, y, w = _flatten_balanced(panel.x[variable][:, columns], panel.y[:, columns])
            p = _weighted_correlation(x, y, w)
            s = _weighted_correlation(_tie_ranks(x), _tie_ranks(y), w)
        end
        direction(value) = !isfinite(value) ? "uncertain" : value > 0 ? "positive" :
            value < 0 ? "negative" : "zero"
        push!(monthly, (String(variable), month_value, p, s, direction(p), direction(s)))
    end
    return (; association, vif, monthly, selected=active)
end

function _haversine(lon1, lat1, lon2, lat2)
    dlat = deg2rad(lat2 - lat1); dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon / 2)^2
    return 6371.0088 * 2 * asin(sqrt(clamp(a, 0.0, 1.0)))
end

function _distance_matrix(lonlat::Matrix{Float64})
    n = size(lonlat, 1); result = zeros(n, n)
    for i in 1:n, j in (i + 1):n
        d = _haversine(lonlat[i, 1], lonlat[i, 2], lonlat[j, 1], lonlat[j, 2])
        result[i, j] = result[j, i] = d
    end
    return result
end

function _adaptive_weights(distances::AbstractVector{<:Real}, neighbors::Int; exclude::Int=0)
    usable = [i for i in eachindex(distances) if i != exclude && isfinite(distances[i])]
    isempty(usable) && return zeros(length(distances))
    k = min(neighbors, length(usable))
    bandwidth = max(sort(distances[usable])[k], eps(Float64))
    weights = zeros(length(distances))
    for i in usable
        ratio = distances[i] / bandwidth
        ratio < 1 && (weights[i] = (1 - ratio^2)^2)
    end
    return weights
end

function _panel_sufficient_statistics(panel, selected::Vector{Symbol}, lonlat::Matrix{Float64}, cfg)
    eligible = findall(panel.eligible_station)
    isempty(eligible) && return nothing
    coords = lonlat[eligible, :]
    coord_mean = vec(mean(coords, dims=1)); coord_scale = vec(std(coords, dims=1))
    any(x -> !isfinite(x) || x <= 0, coord_scale) && return nothing
    p = 3 + length(selected)
    A = [zeros(p, p) for _ in eligible]
    b = [zeros(p) for _ in eligible]
    rows = [Vector{Vector{Float64}}() for _ in eligible]
    ys = [Float64[] for _ in eligible]
    for (position, i) in enumerate(eligible)
        valid = [t for t in axes(panel.y, 2) if isfinite(panel.y[i, t]) &&
            all(isfinite(panel.x[v][i, t]) for v in selected)]
        isempty(valid) && continue
        weight = 1 / length(valid)
        for t in valid
            x = [1.0, (coords[position, 1] - coord_mean[1]) / coord_scale[1],
                (coords[position, 2] - coord_mean[2]) / coord_scale[2],
                (panel.x[v][i, t] for v in selected)...]
            A[position] .+= weight .* (x * x')
            b[position] .+= weight .* x .* panel.y[i, t]
            push!(rows[position], x); push!(ys[position], panel.y[i, t])
        end
    end
    return (; A, b, rows, ys, coords, eligible, p)
end

function _local_coefficients(stats, distances, neighbors, cfg; block_order=collect(eachindex(stats.A)))
    n = length(stats.A); coefficients = fill(NaN, n, stats.p)
    for target in 1:n
        weights = _adaptive_weights(view(distances, target, :), neighbors)
        lhs, rhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p), zeros(stats.p)
        for spatial_index in 1:n
            source = block_order[spatial_index]
            lhs .+= weights[spatial_index] .* stats.A[source]
            rhs .+= weights[spatial_index] .* stats.b[source]
        end
        try
            coefficients[target, :] .= lhs \ rhs
        catch
        end
    end
    return coefficients
end

function _loocv_panel_rmse(stats, distances, neighbors, cfg)
    n = length(stats.A); sse = 0.0; count_station = 0
    for target in 1:n
        weights = _adaptive_weights(view(distances, target, :), neighbors; exclude=target)
        lhs, rhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p), zeros(stats.p)
        for i in 1:n
            lhs .+= weights[i] .* stats.A[i]; rhs .+= weights[i] .* stats.b[i]
        end
        beta = try lhs \ rhs catch; continue end
        isempty(stats.rows[target]) && continue
        station_sse = mean((stats.ys[target][j] - dot(stats.rows[target][j], beta))^2
            for j in eachindex(stats.ys[target]))
        sse += station_sse; count_station += 1
    end
    return count_station == n ? sqrt(sse / n) : Inf
end

function panel_spatial_variability_test(
    panel, selected::Vector{Symbol}, lonlat::Matrix{Float64}, cfg::ERA5SelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed + 1),
)
    bandwidth_scan = DataFrame(neighbors=Int[], loocv_rmse=Float64[], available=Bool[])
    variability = DataFrame(
        variable=String[], statistic=Float64[], pvalue=Float64[], qvalue=Float64[],
        role=String[], neighbors=Union{Missing,Int}[], status=String[],
    )
    isempty(selected) && return (; bandwidth_scan, variability, bandwidth=missing)
    stats = _panel_sufficient_statistics(panel, selected, lonlat, cfg)
    if stats === nothing || length(stats.eligible) <= maximum([3, minimum(cfg.bandwidth_candidates)])
        for variable in selected
            push!(variability, (String(variable), NaN, NaN, NaN, "uncertain", missing, "insufficient_stations"))
        end
        return (; bandwidth_scan, variability, bandwidth=missing)
    end
    distances = _distance_matrix(stats.coords)
    candidates = filter(c -> max(4, stats.p + 1) <= c <= length(stats.eligible) - 1,
        unique(cfg.bandwidth_candidates))
    for candidate in candidates
        rmse = _loocv_panel_rmse(stats, distances, candidate, cfg)
        push!(bandwidth_scan, (candidate, rmse, isfinite(rmse)))
    end
    available = filter(:available => identity, bandwidth_scan)
    if nrow(available) == 0
        for variable in selected
            push!(variability, (String(variable), NaN, NaN, NaN, "uncertain", missing, "no_bandwidth"))
        end
        return (; bandwidth_scan, variability, bandwidth=missing)
    end
    bandwidth = available.neighbors[argmin(available.loocv_rmse)]
    coefficients = _local_coefficients(stats, distances, bandwidth, cfg)
    observed = [var(coefficients[:, 3 + j]) for j in eachindex(selected)]
    exceed = zeros(Int, length(selected))
    for _ in 1:cfg.spatial_permutations
        order = station_block_permutation(length(stats.A), rng)
        permuted = _local_coefficients(stats, distances, bandwidth, cfg; block_order=order)
        for j in eachindex(selected)
            value = var(permuted[:, 3 + j])
            exceed[j] += isfinite(value) && value >= observed[j]
        end
    end
    pvalues = cfg.spatial_permutations > 0 ?
        (exceed .+ 1) ./ (cfg.spatial_permutations + 1) : fill(NaN, length(selected))
    qvalues = all(isfinite, pvalues) ? _bh_adjust(pvalues) : fill(NaN, length(selected))
    for j in eachindex(selected)
        role = isfinite(qvalues[j]) ? (qvalues[j] < cfg.q_threshold ? "local" : "global") : "uncertain"
        push!(variability, (String(selected[j]), observed[j], pvalues[j], qvalues[j],
            role, bandwidth, role == "uncertain" ? "test_failed" : "ok"))
    end
    return (; bandwidth_scan, variability, bandwidth)
end

function _annotate!(table::DataFrame, product::String, scheme::String, fold::Int)
    insertcols!(table, 1, :product => fill(product, nrow(table)),
        :scheme => fill(scheme, nrow(table)), :fold => fill(fold, nrow(table)))
    return table
end

function _append_table!(target::DataFrame, source::DataFrame)
    if ncol(target) == 0
        for column in propertynames(source)
            target[!, column] = similar(source[!, column], 0)
        end
    end
    nrow(source) > 0 && append!(target, source; cols=:union)
    return target
end

function run_era5_variable_selection(
    cfg::ERA5SelectionConfig, products::Vector{String}, station_ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64}, satellite::AbstractDict,
    era5::AbstractDict, lonlat::Matrix{Float64}; data_qc::Union{Nothing,DataFrame}=nothing,
)
    mkpath(cfg.outdir)
    n = length(station_ids)
    size(Yobs) == (n, length(times)) || throw(DimensionMismatch("observation dimensions differ"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("coordinate dimensions differ"))
    folds = balanced_spatial_folds(station_ids, lonlat; k=cfg.k)
    fold_table = DataFrame(station_id=station_ids, fold=folds, lon=lonlat[:, 1], lat=lonlat[:, 2])
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), fold_table)
    data_qc !== nothing && CSV.write(joinpath(cfg.outdir, "era5_data_qc.csv"), data_qc)

    qc_all, association_all, monthly_all, vif_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    bandwidth_all, variability_all, role_all, status_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    specifications = DataFrame(product=String[], variable=String[], selected=Bool[], role=String[],
        association_qvalue=Float64[], spatial_qvalue=Union{Missing,Float64}[], status=String[])
    schemes = [("full_data", 0, collect(1:n))]
    append!(schemes, [("spatial_cv", fold, findall(!=(fold), folds)) for fold in 1:cfg.k])
    for product in products
        Ysat = Float64.(satellite[product])
        for (scheme, fold, train_indices) in schemes
            run_seed = cfg.seed + 1000 * findfirst(==(product), products) + 10 * fold
            panel = prepare_dynamic_panel(Yobs, Ysat, era5, train_indices, times, cfg)
            qc = copy(panel.qc); _annotate!(qc, product, scheme, fold); _append_table!(qc_all, qc)
            if panel.qc.status[1] != "ok"
                for variable in ERA5_VARIABLES
                    push!(role_all, (product=product, scheme=scheme, fold=fold,
                        variable=String(variable), selected=false, role="uncertain", status="insufficient_stations"))
                    if scheme == "full_data"
                        push!(specifications, (product, String(variable), false, "uncertain",
                            NaN, missing, "insufficient_stations"))
                    end
                end
                push!(status_all, (product=product, scheme=scheme, fold=fold,
                    status="failed", reason="insufficient_stations"))
                continue
            end
            screen = dynamic_panel_screen(panel, times, cfg; rng=MersenneTwister(run_seed))
            for table in (screen.association, screen.monthly, screen.vif)
                _annotate!(table, product, scheme, fold)
            end
            _append_table!(association_all, screen.association)
            _append_table!(monthly_all, screen.monthly)
            _append_table!(vif_all, screen.vif)
            spatial = panel_spatial_variability_test(panel, screen.selected, lonlat[train_indices, :], cfg;
                rng=MersenneTwister(run_seed + 1))
            _annotate!(spatial.bandwidth_scan, product, scheme, fold)
            _annotate!(spatial.variability, product, scheme, fold)
            _append_table!(bandwidth_all, spatial.bandwidth_scan)
            _append_table!(variability_all, spatial.variability)
            role_map = Dict(Symbol(row.variable) => (role=row.role, q=row.qvalue, status=row.status)
                for row in eachrow(spatial.variability))
            association_q = Dict(Symbol(row.variable) => row.qvalue for row in eachrow(screen.association))
            for variable in ERA5_VARIABLES
                selected = variable in screen.selected
                info = selected ? get(role_map, variable, (role="uncertain", q=NaN, status="test_failed")) :
                    (role="not_selected", q=NaN, status="not_selected")
                push!(role_all, (product=product, scheme=scheme, fold=fold,
                    variable=String(variable), selected=selected, role=info.role, status=info.status))
                if scheme == "full_data"
                    push!(specifications, (product, String(variable), selected, info.role,
                        association_q[variable], isfinite(info.q) ? info.q : missing, info.status))
                end
            end
            push!(status_all, (product=product, scheme=scheme, fold=fold,
                status="ok", reason=isempty(screen.selected) ? "no_variable_selected" : ""))
        end
    end
    stability = DataFrame(product=String[], variable=String[], selection_rate=Float64[],
        local_rate=Float64[], global_rate=Float64[], uncertain_rate=Float64[])
    for product in products, variable in ERA5_VARIABLES
        rows = filter([:product, :scheme, :variable] =>
            (p, s, v) -> p == product && s == "spatial_cv" && v == String(variable), role_all)
        denominator = max(nrow(rows), 1)
        push!(stability, (product, String(variable), count(rows.selected) / denominator,
            count(==("local"), rows.role) / denominator, count(==("global"), rows.role) / denominator,
            count(==("uncertain"), rows.role) / denominator))
    end
    consensus = DataFrame(variable=String[], selected_product_count=Int[], consensus_selected=Bool[],
        local_product_count=Int[], global_product_count=Int[], consensus_role=String[])
    for variable in ERA5_VARIABLES
        rows = filter(:variable => ==(String(variable)), specifications)
        selected_count = count(rows.selected)
        local_count = count(==("local"), rows.role); global_count = count(==("global"), rows.role)
        role = selected_count == 0 ? "not_selected" : local_count >= 2 ? "local" :
            global_count >= 2 ? "global" : selected_count >= 2 ? "role_unstable" : "product_specific"
        push!(consensus, (String(variable), selected_count, selected_count >= 2,
            local_count, global_count, role))
    end
    outputs = Dict(
        "era5_fold_quality_control.csv" => qc_all,
        "era5_fold_association.csv" => association_all,
        "era5_fold_monthly_direction.csv" => monthly_all,
        "era5_fold_vif.csv" => vif_all,
        "era5_fold_bandwidth_scan.csv" => bandwidth_all,
        "era5_fold_spatial_variability.csv" => variability_all,
        "era5_fold_roles.csv" => role_all,
        "era5_role_stability.csv" => stability,
        "era5_final_full_data_spec.csv" => specifications,
        "era5_cross_product_consensus.csv" => consensus,
        "run_status.csv" => status_all,
    )
    for (filename, table) in outputs
        CSV.write(joinpath(cfg.outdir, filename), table)
    end
    return (; specifications, stability, consensus, status=status_all)
end

end
