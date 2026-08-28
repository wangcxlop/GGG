module JointVariableSelection

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

# All three are loaded through src/load_modules.jl before this file, so each is shared
# rather than recompiled into this module. See src/load_modules.jl.
using Main.DEMTerrainExperiment
using Main.ERA5VariableSelection
using Main.NDVIVariableSelection

export JointSelectionConfig, JOINT_GROUP_ORDER, prepare_joint_panel, joint_vif
export joint_spatial_variability_test, run_joint_variable_selection, select_joint_covariates

const JOINT_GROUP_ORDER = [
    "elevation", "slope", "aspect", "t2m_c", "d2m_c",
    "u10", "v10", "sp_hpa", "ndvi",
]
const GROUP_COLUMNS = Dict(
    "elevation" => [:elevation_m],
    "slope" => [:slope_deg],
    "aspect" => [:aspect_sin, :aspect_cos],
    "t2m_c" => [:t2m_c],
    "d2m_c" => [:d2m_c],
    "u10" => [:u10],
    "v10" => [:v10],
    "sp_hpa" => [:sp_hpa],
    "ndvi" => [:ndvi],
)
const GROUP_FAMILY = Dict(
    group => (group in ("elevation", "slope", "aspect") ? "dem" :
        group == "ndvi" ? "ndvi" : "era5") for group in JOINT_GROUP_ORDER
)

_predictor_label(group::String) = group == "ndvi" ? "ndvi_qc" :
    join(String.(GROUP_COLUMNS[group]), "+")

Base.@kwdef struct JointSelectionConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    min_stations_per_time::Int = 12
    min_ndvi_wet_hours_per_period::Int = 5
    min_ndvi_periods_per_station::Int = 8
    max_ndvi_age_days::Int = 32
    k::Int = 5
    seed::Int = 20260818
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    independent_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
end

function _family_configs(cfg::JointSelectionConfig)
    era = ERA5SelectionConfig(
        outdir=cfg.outdir, wet_threshold=cfg.wet_threshold,
        min_wet_hours=cfg.min_wet_hours,
        min_stations_per_time=cfg.min_stations_per_time,
        k=cfg.k, seed=20260816,
        bandwidth_candidates=cfg.bandwidth_candidates,
        association_permutations=cfg.independent_permutations,
        spatial_permutations=cfg.spatial_permutations,
        q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
        ridge=cfg.ridge,
    )
    ndvi = NDVISelectionConfig(
        outdir=cfg.outdir, wet_threshold=cfg.wet_threshold,
        min_wet_hours_per_period=cfg.min_ndvi_wet_hours_per_period,
        min_periods_per_station=cfg.min_ndvi_periods_per_station,
        min_stations_per_period=cfg.min_stations_per_time,
        max_age_days=cfg.max_ndvi_age_days,
        k=cfg.k, seed=20260817,
        bandwidth_candidates=cfg.bandwidth_candidates,
        association_permutations=cfg.independent_permutations,
        spatial_permutations=cfg.spatial_permutations,
        q_threshold=cfg.q_threshold, ridge=cfg.ridge,
    )
    return era, ndvi
end

function _candidate_rows(
    product::String, Yobs::Matrix{Float64}, Ysat::Matrix{Float64},
    terrain::DataFrame, era5::AbstractDict, ndvi_aligned,
    train_indices::Vector{Int}, station_ids::Vector{String},
    times::Vector{DateTime}, cfg::JointSelectionConfig,
    run_seed::Int, dem_seed::Int,
)
    rows = NamedTuple[]
    active = String[]
    # DEM screening uses station mean wet residuals and grouped terrain VIF.
    response, _ = mean_wet_residual(
        Yobs[train_indices, :], Ysat[train_indices, :];
        threshold=cfg.wet_threshold, min_hours=cfg.min_wet_hours,
    )
    dem_screen, _, dem_selected = terrain_screen(
        terrain[train_indices, :], response; product,
        permutations=cfg.independent_permutations,
        q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
        seed=dem_seed,
    )
    for group in ("elevation", "slope", "aspect")
        row = only(eachrow(filter(:variable_group => ==(group), dem_screen)))
        selected = group in dem_selected
        push!(rows, (;
            family="dem", variable_group=group,
            predictor_columns=_predictor_label(group),
            independent_qvalue=Float64(row.qvalue),
            direction_stable=Bool(row.direction_stable),
            independent_selected=selected,
            independent_status=selected ? "selected" : String(row.exclusion_reason),
        ))
        selected && push!(active, group)
    end

    era_cfg, ndvi_cfg = _family_configs(cfg)
    era_panel = prepare_dynamic_panel(Yobs, Ysat, era5, train_indices, times, era_cfg)
    if era_panel.qc.status[1] == "ok"
        era_screen = dynamic_panel_screen(
            era_panel, times, era_cfg; rng=MersenneTwister(20260816 + run_seed),
        )
        for variable in ERA5_VARIABLES
            row = only(eachrow(filter(:variable => ==(String(variable)), era_screen.association)))
            selected = variable in era_screen.selected
            push!(rows, (;
                family="era5", variable_group=String(variable),
                predictor_columns=String(variable), independent_qvalue=Float64(row.qvalue),
                direction_stable=Bool(row.direction_stable), independent_selected=selected,
                independent_status=selected ? "selected" :
                    (!row.direction_stable ? "direction_conflict" : "q_not_significant_or_family_vif"),
            ))
            selected && push!(active, String(variable))
        end
    else
        for variable in ERA5_VARIABLES
            push!(rows, (;
                family="era5", variable_group=String(variable), predictor_columns=String(variable),
                independent_qvalue=NaN, direction_stable=false, independent_selected=false,
                independent_status="insufficient_stations",
            ))
        end
    end

    ndvi_panel = prepare_ndvi_panel(
        Yobs, Ysat, ndvi_aligned, train_indices,
        station_ids, ndvi_cfg,
    )
    if ndvi_panel.qc.status[1] == "ok"
        ndvi_screen = screen_ndvi(
            ndvi_panel, ndvi_aligned.periods, ndvi_cfg;
            rng=MersenneTwister(20260817 + run_seed),
        )
        row = ndvi_screen.association[1, :]
        selected = !isempty(ndvi_screen.selected)
        push!(rows, (;
            family="ndvi", variable_group="ndvi", predictor_columns="ndvi_qc",
            independent_qvalue=Float64(row.qvalue),
            direction_stable=Bool(row.direction_stable), independent_selected=selected,
            independent_status=selected ? "selected" :
                (!row.direction_stable ? "direction_conflict" : "q_not_significant"),
        ))
        selected && push!(active, "ndvi")
    else
        push!(rows, (;
            family="ndvi", variable_group="ndvi", predictor_columns="ndvi_qc",
            independent_qvalue=NaN, direction_stable=false, independent_selected=false,
            independent_status="insufficient_stations_or_variation",
        ))
    end
    active = [group for group in JOINT_GROUP_ORDER if group in Set(active)]
    return DataFrame(rows), active
end

function _raw_predictors(
    groups::Vector{String}, terrain::DataFrame, era5::AbstractDict,
    ndvi_aligned, train_indices::Vector{Int}, time_count::Int,
)
    predictors = Dict{Symbol,Matrix{Float64}}()
    column_groups = String[]
    columns = Symbol[]
    for group in groups
        for column in GROUP_COLUMNS[group]
            matrix = if GROUP_FAMILY[group] == "dem"
                repeat(reshape(Float64.(terrain[train_indices, column]), :, 1), 1, time_count)
            elseif GROUP_FAMILY[group] == "era5"
                Float64.(era5[column][train_indices, :])
            else
                Float64.(ndvi_aligned.values[train_indices, :])
            end
            predictors[column] = matrix
            push!(columns, column); push!(column_groups, group)
        end
    end
    return predictors, columns, column_groups
end

function prepare_joint_panel(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned, train_indices::Vector{Int},
    groups::Vector{String}, cfg::JointSelectionConfig,
)
    n, nt = length(train_indices), size(Yobs, 2)
    predictors, columns, column_groups = _raw_predictors(
        groups, terrain, era5, ndvi_aligned, train_indices, nt,
    )
    residual = Yobs[train_indices, :] .- Ysat[train_indices, :]
    mask = falses(n, nt)
    for i in 1:n, time in 1:nt
        global_station = train_indices[i]
        mask[i, time] = isfinite(Yobs[global_station, time]) &&
            isfinite(Ysat[global_station, time]) &&
            Yobs[global_station, time] >= cfg.wet_threshold &&
            all(isfinite(predictors[column][i, time]) for column in columns)
    end
    previous = -1
    while count(mask) != previous
        previous = count(mask)
        if "ndvi" in groups
            for i in 1:n
                global_station = train_indices[i]
                period_values = ndvi_aligned.source_period[global_station, :]
                for period in unique(filter(>(0), period_values[mask[i, :]]))
                    period_times = findall(mask[i, :] .& (period_values .== period))
                    length(period_times) >= cfg.min_ndvi_wet_hours_per_period ||
                        (mask[i, period_times] .= false)
                end
                retained_periods = unique(filter(>(0), period_values[mask[i, :]]))
                length(retained_periods) >= cfg.min_ndvi_periods_per_station ||
                    (mask[i, :] .= false)
            end
        end
        station_ok = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
        mask[.!station_ok, :] .= false
        time_ok = vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time
        mask[:, .!time_ok] .= false
    end
    eligible_station = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
    y = fill(NaN, n, nt)
    x = Dict(column => fill(NaN, n, nt) for column in columns)
    for time in 1:nt
        stations = findall(mask[:, time])
        length(stations) >= cfg.min_stations_per_time || continue
        y[stations, time] .= residual[stations, time] .- mean(residual[stations, time])
        for column in columns
            values = predictors[column][stations, time]
            x[column][stations, time] .= values .- mean(values)
        end
    end
    weights = zeros(n, nt)
    uses_ndvi = "ndvi" in groups
    for i in 1:n
        valid_times = findall(mask[i, :])
        isempty(valid_times) && continue
        if uses_ndvi
            global_station = train_indices[i]
            period_values = ndvi_aligned.source_period[global_station, valid_times]
            valid_periods = sort(unique(filter(>(0), period_values)))
            isempty(valid_periods) && continue
            for period in valid_periods
                period_times = valid_times[period_values .== period]
                weights[i, period_times] .= 1 / (length(valid_periods) * length(period_times))
            end
        else
            weights[i, valid_times] .= 1 / length(valid_times)
        end
    end
    scale_rows = NamedTuple[]
    for column in columns
        values = Float64[]; value_weights = Float64[]
        for index in eachindex(y)
            isfinite(y[index]) && isfinite(x[column][index]) && weights[index] > 0 || continue
            push!(values, x[column][index]); push!(value_weights, weights[index])
        end
        total = sum(value_weights)
        mean_value = total > 0 ? dot(values, value_weights) / total : NaN
        scale = total > 0 ? sqrt(dot((values .- mean_value) .^ 2, value_weights) / total) : NaN
        if isfinite(scale) && scale > 0
            for index in eachindex(x[column])
                isfinite(x[column][index]) &&
                    (x[column][index] = (x[column][index] - mean_value) / scale)
            end
        end
        push!(scale_rows, (; column=String(column), mean=mean_value, scale))
    end
    station_weight_sums = vec(sum(weights, dims=2))
    positive = station_weight_sums[station_weight_sums .> 0]
    qc = DataFrame(
        training_station_count=[n], eligible_station_count=[count(eligible_station)],
        eligible_time_count=[count(vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time)],
        used_station_hours=[count(mask)], group_count=[length(groups)],
        predictor_column_count=[length(columns)], uses_ndvi_hierarchical_weight=[uses_ndvi],
        minimum_station_weight=[isempty(positive) ? NaN : minimum(positive)],
        maximum_station_weight=[isempty(positive) ? NaN : maximum(positive)],
        status=[count(eligible_station) >= cfg.min_stations_per_time &&
            all(isfinite(row.scale) && row.scale > 0 for row in scale_rows) ? "ok" :
            "insufficient_stations_or_variation"],
    )
    weight_rows = NamedTuple[]
    for i in 1:n
        global_station = train_indices[i]
        valid_times = findall(weights[i, :] .> 0)
        period_weights = Float64[]
        period_count = 0
        if uses_ndvi && !isempty(valid_times)
            periods = ndvi_aligned.source_period[global_station, valid_times]
            valid_periods = sort(unique(filter(>(0), periods)))
            period_count = length(valid_periods)
            period_weights = [sum(weights[i, valid_times[periods .== period]])
                for period in valid_periods]
        end
        push!(weight_rows, (;
            local_station_index=i, global_station_index=global_station,
            eligible=eligible_station[i], used_hour_count=length(valid_times),
            ndvi_period_count=period_count, total_weight=station_weight_sums[i],
            minimum_period_weight=isempty(period_weights) ? NaN : minimum(period_weights),
            maximum_period_weight=isempty(period_weights) ? NaN : maximum(period_weights),
        ))
    end
    return (; y, x, weights, mask, eligible_station, columns, column_groups,
        groups, scales=DataFrame(scale_rows), weight_audit=DataFrame(weight_rows), qc)
end

function _weighted_rows(panel, active::Vector{String})
    columns = Symbol[]; column_groups = String[]
    for group in active, column in GROUP_COLUMNS[group]
        push!(columns, column); push!(column_groups, group)
    end
    rows = Vector{Vector{Float64}}(); weights = Float64[]
    for i in axes(panel.y, 1), time in axes(panel.y, 2)
        panel.weights[i, time] > 0 && isfinite(panel.y[i, time]) &&
            all(isfinite(panel.x[column][i, time]) for column in columns) || continue
        push!(rows, [panel.x[column][i, time] for column in columns])
        push!(weights, panel.weights[i, time])
    end
    X = isempty(rows) ? zeros(0, length(columns)) : reduce(vcat, permutedims.(rows))
    return X, weights, columns, column_groups
end

function _column_vifs(X::Matrix{Float64}, weights::Vector{Float64})
    p = size(X, 2)
    p == 0 && return Float64[]
    p == 1 && return [1.0]
    sw = sqrt.(weights ./ sum(weights))
    result = fill(Inf, p)
    for j in 1:p
        others = setdiff(1:p, [j])
        design = hcat(ones(size(X, 1)), X[:, others])
        beta = (design .* sw) \ (X[:, j] .* sw)
        prediction = design * beta
        mean_value = dot(X[:, j], weights) / sum(weights)
        rss = sum(weights .* (X[:, j] .- prediction) .^ 2)
        tss = sum(weights .* (X[:, j] .- mean_value) .^ 2)
        r2 = tss > 0 ? clamp(1 - rss / tss, 0.0, 1.0) : 1.0
        result[j] = r2 < 1 ? 1 / (1 - r2) : Inf
    end
    return result
end

function joint_vif(
    panel, candidates::Vector{String}, independent_q::Dict{String,Float64},
    cfg::JointSelectionConfig,
)
    active = copy(candidates)
    audit = DataFrame(
        iteration=Int[], variable_group=String[], predictor_columns=String[],
        group_vif=Float64[], removed=Bool[], reason=String[],
    )
    iteration = 1
    while !isempty(active)
        X, weights, columns, column_groups = _weighted_rows(panel, active)
        column_vifs = _column_vifs(X, weights)
        group_values = Dict(group => maximum(column_vifs[column_groups .== group]) for group in active)
        maximum_value = maximum(Base.values(group_values))
        if maximum_value < cfg.vif_threshold
            for group in active
                push!(audit, (iteration, group, join(String.(GROUP_COLUMNS[group]), "+"),
                    group_values[group], false, "retained"))
            end
            break
        end
        tied = [group for group in active if isapprox(group_values[group], maximum_value;
            rtol=1e-10, atol=1e-10) || (isinf(group_values[group]) && isinf(maximum_value))]
        remove_group = sort(tied; by=group ->
            (-get(independent_q, group, Inf), findfirst(==(group), JOINT_GROUP_ORDER)))[1]
        for group in active
            push!(audit, (iteration, group, join(String.(GROUP_COLUMNS[group]), "+"),
                group_values[group], group == remove_group,
                group == remove_group ? "joint_vif_ge_threshold" : "pending"))
        end
        filter!(!=(remove_group), active)
        iteration += 1
    end
    return audit, active
end

function _panel_stats(panel, groups::Vector{String}, lonlat::Matrix{Float64}, cfg)
    eligible = findall(panel.eligible_station)
    isempty(eligible) && return nothing
    coords = lonlat[eligible, :]
    means = vec(mean(coords, dims=1)); scales = vec(std(coords, dims=1))
    any(x -> !isfinite(x) || x <= 0, scales) && return nothing
    columns = reduce(vcat, [GROUP_COLUMNS[group] for group in groups]; init=Symbol[])
    p = 3 + length(columns)
    A = [zeros(p, p) for _ in eligible]
    b = [zeros(p) for _ in eligible]
    rows = [Vector{Vector{Float64}}() for _ in eligible]
    ys = [Float64[] for _ in eligible]
    row_weights = [Float64[] for _ in eligible]
    for (position, i) in enumerate(eligible), time in axes(panel.y, 2)
        weight = panel.weights[i, time]
        weight > 0 && isfinite(panel.y[i, time]) &&
            all(isfinite(panel.x[column][i, time]) for column in columns) || continue
        design = [1.0, (coords[position, 1] - means[1]) / scales[1],
            (coords[position, 2] - means[2]) / scales[2],
            (panel.x[column][i, time] for column in columns)...]
        A[position] .+= weight .* (design * design')
        b[position] .+= weight .* design .* panel.y[i, time]
        push!(rows[position], design); push!(ys[position], panel.y[i, time])
        push!(row_weights[position], weight)
    end
    group_indices = Dict{String,Vector{Int}}()
    offset = 3
    for group in groups
        width = length(GROUP_COLUMNS[group])
        group_indices[group] = collect((offset + 1):(offset + width))
        offset += width
    end
    return (; A, b, rows, ys, row_weights, coords, p, group_indices, eligible)
end

function _coefficients(stats, distances, neighbors, cfg; order=collect(eachindex(stats.A)))
    n = length(stats.A); result = fill(NaN, n, stats.p)
    for target in 1:n
        spatial_weights = ERA5VariableSelection._adaptive_weights(view(distances, target, :), neighbors)
        lhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p); rhs = zeros(stats.p)
        for spatial_index in 1:n
            source = order[spatial_index]
            lhs .+= spatial_weights[spatial_index] .* stats.A[source]
            rhs .+= spatial_weights[spatial_index] .* stats.b[source]
        end
        try result[target, :] .= lhs \ rhs catch end
    end
    return result
end

function _loocv(stats, distances, neighbors, cfg)
    total = 0.0; station_count = 0; n = length(stats.A)
    for target in 1:n
        spatial_weights = ERA5VariableSelection._adaptive_weights(
            view(distances, target, :), neighbors; exclude=target,
        )
        lhs = cfg.ridge .* Matrix{Float64}(I, stats.p, stats.p); rhs = zeros(stats.p)
        for i in 1:n
            lhs .+= spatial_weights[i] .* stats.A[i]; rhs .+= spatial_weights[i] .* stats.b[i]
        end
        beta = try lhs \ rhs catch; continue end
        isempty(stats.rows[target]) && continue
        station_sse = sum(stats.row_weights[target][j] *
            (stats.ys[target][j] - dot(stats.rows[target][j], beta))^2
            for j in eachindex(stats.ys[target])) / sum(stats.row_weights[target])
        total += station_sse; station_count += 1
    end
    return station_count == n ? sqrt(total / n) : Inf
end

function joint_spatial_variability_test(
    panel, groups::Vector{String}, lonlat::Matrix{Float64}, cfg::JointSelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed),
)
    scan = DataFrame(neighbors=Int[], loocv_rmse=Float64[], available=Bool[])
    result = DataFrame(
        variable_group=String[], statistic=Float64[], pvalue=Float64[], qvalue=Float64[],
        role=String[], neighbors=Union{Missing,Int}[], status=String[],
    )
    isempty(groups) && return (; scan, result, bandwidth=missing)
    stats = _panel_stats(panel, groups, lonlat, cfg)
    if stats === nothing
        for group in groups
            push!(result, (group, NaN, NaN, NaN, "uncertain", missing, "insufficient_stations"))
        end
        return (; scan, result, bandwidth=missing)
    end
    distances = ERA5VariableSelection._distance_matrix(stats.coords)
    candidates = filter(c -> max(4, stats.p + 1) <= c <= length(stats.eligible) - 1,
        unique(cfg.bandwidth_candidates))
    for candidate in candidates
        value = _loocv(stats, distances, candidate, cfg)
        push!(scan, (candidate, value, isfinite(value)))
    end
    available = filter(:available => identity, scan)
    if nrow(available) == 0
        for group in groups
            push!(result, (group, NaN, NaN, NaN, "uncertain", missing, "no_bandwidth"))
        end
        return (; scan, result, bandwidth=missing)
    end
    bandwidth = available.neighbors[argmin(available.loocv_rmse)]
    coefficients = _coefficients(stats, distances, bandwidth, cfg)
    observed = [sum(var(coefficients[:, index]) for index in stats.group_indices[group]) for group in groups]
    exceed = zeros(Int, length(groups))
    for _ in 1:cfg.spatial_permutations
        order = station_block_permutation(length(stats.A), rng)
        permuted = _coefficients(stats, distances, bandwidth, cfg; order)
        for (group_index, group) in enumerate(groups)
            value = sum(var(permuted[:, index]) for index in stats.group_indices[group])
            exceed[group_index] += isfinite(value) && value >= observed[group_index]
        end
    end
    pvalues = (exceed .+ 1) ./ (cfg.spatial_permutations + 1)
    qvalues = DEMTerrainExperiment.bh_adjust(pvalues)
    for (index, group) in enumerate(groups)
        role = qvalues[index] < cfg.q_threshold ? "local" : "global"
        push!(result, (group, observed[index], pvalues[index], qvalues[index],
            role, bandwidth, "ok"))
    end
    return (; scan, result, bandwidth)
end

"""
Select joint covariates for one (product, `train_indices`) cell: independent per-family
screening → cross-family joint VIF pruning → local/global spatial-role test. Used both by
`run_joint_variable_selection`'s own full-data/stability reporting and, per training fold, by
the interpolation benchmark's nested selection mode — the only difference between the two
callers is which `train_indices` they pass in.
"""
function select_joint_covariates(
    product::String, Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned, train_indices::Vector{Int}, station_ids::Vector{String},
    times::Vector{DateTime}, lonlat::Matrix{Float64}, cfg::JointSelectionConfig,
    run_seed::Int, dem_seed::Int,
)
    candidates, active = _candidate_rows(
        product, Yobs, Ysat, terrain, era5, ndvi_aligned,
        train_indices, station_ids, times, cfg, run_seed, dem_seed,
    )
    if isempty(active)
        qc = DataFrame(
            training_station_count=[length(train_indices)], eligible_station_count=[0],
            eligible_time_count=[0], used_station_hours=[0], group_count=[0],
            predictor_column_count=[0], uses_ndvi_hierarchical_weight=[false],
            minimum_station_weight=[NaN], maximum_station_weight=[NaN],
            status=["no_independent_candidates"],
        )
        return (;
            candidates, active, qc, scales=DataFrame(), weight_audit=DataFrame(),
            vif=DataFrame(), retained=String[], spatial_scan=DataFrame(),
            spatial_result=DataFrame(), role_map=Dict{String,String}(),
            status="ok", status_reason="no_independent_candidates",
        )
    end
    independent_q = Dict(String(row.variable_group) => Float64(row.independent_qvalue)
        for row in eachrow(candidates))
    panel = prepare_joint_panel(Yobs, Ysat, terrain, era5, ndvi_aligned, train_indices, active, cfg)
    if panel.qc.status[1] != "ok"
        retained = String[]; vif = DataFrame()
    else
        vif, retained = joint_vif(panel, active, independent_q, cfg)
    end
    spatial = panel.qc.status[1] == "ok" ? joint_spatial_variability_test(
        panel, retained, lonlat[train_indices, :], cfg; rng=MersenneTwister(cfg.seed + run_seed),
    ) : joint_spatial_variability_test(panel, String[], lonlat[train_indices, :], cfg)
    role_map = Dict(String(row.variable_group) => String(row.role)
        for row in eachrow(spatial.result) if row.role in ("local", "global"))
    status = panel.qc.status[1] == "ok" ? "ok" : "failed"
    status_reason = status == "ok" ? (isempty(retained) ? "all_removed_by_joint_vif" : "") :
        String(panel.qc.status[1])
    return (;
        candidates, active, qc=panel.qc, scales=panel.scales, weight_audit=panel.weight_audit,
        vif, retained, spatial_scan=spatial.scan, spatial_result=spatial.result,
        role_map, status, status_reason,
    )
end

function _annotate!(table::DataFrame, product::String, scheme::String, fold::Int)
    insertcols!(table, 1, :product => fill(product, nrow(table)),
        :scheme => fill(scheme, nrow(table)), :fold => fill(fold, nrow(table)))
    return table
end

function _append!(target::DataFrame, source::DataFrame)
    if ncol(target) == 0
        for column in propertynames(source)
            target[!, column] = similar(source[!, column], 0)
        end
    end
    nrow(source) > 0 && append!(target, source; cols=:union)
    return target
end

function run_joint_variable_selection(
    cfg::JointSelectionConfig, products::Vector{String}, station_ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64}, satellite::AbstractDict,
    terrain::DataFrame, era5::AbstractDict, ndvi_table::DataFrame,
    lonlat::Matrix{Float64}; prerequisite_audit::Union{Nothing,DataFrame}=nothing,
)
    mkpath(cfg.outdir)
    n = length(station_ids)
    length(times) == size(Yobs, 2) || throw(DimensionMismatch("observation time dimension differs"))
    size(Yobs, 1) == n || throw(DimensionMismatch("observation station dimension differs"))
    nrow(terrain) == n || throw(DimensionMismatch("terrain station dimension differs"))
    String.(terrain.station_id) == station_ids ||
        throw(ArgumentError("terrain rows must follow common station order"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("coordinate dimensions differ"))
    all(size(satellite[product]) == size(Yobs) for product in products) ||
        throw(DimensionMismatch("satellite dimensions differ"))
    all(size(era5[variable]) == size(Yobs) for variable in ERA5_VARIABLES) ||
        throw(DimensionMismatch("ERA5 dimensions differ"))
    folds = balanced_spatial_folds(station_ids, lonlat; k=cfg.k)
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), DataFrame(
        station_id=station_ids, fold=folds, lon=lonlat[:, 1], lat=lonlat[:, 2],
    ))
    prerequisite_audit !== nothing && CSV.write(
        joinpath(cfg.outdir, "prerequisite_full_audit.csv"), prerequisite_audit,
    )
    ndvi_aligned = align_ndvi_asof(
        ndvi_table, station_ids, times; max_age_days=cfg.max_ndvi_age_days,
    )
    candidate_all, qc_all, scale_all, weight_all, vif_all =
        DataFrame(), DataFrame(), DataFrame(), DataFrame(), DataFrame()
    bandwidth_all, spatial_all, role_all, status_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    final_spec = DataFrame(
        product=String[], family=String[], variable_group=String[], predictor_columns=String[],
        independent_selected=Bool[], independent_qvalue=Float64[], joint_vif_retained=Bool[],
        role=String[], final_included=Bool[], exclusion_reason=String[],
    )
    schemes = [("full_data", 0, collect(1:n))]
    append!(schemes, [("spatial_cv", fold, findall(!=(fold), folds)) for fold in 1:cfg.k])
    for (product_index, product) in enumerate(products)
        Ysat = Float64.(satellite[product])
        for (scheme, fold, train_indices) in schemes
            run_seed = 1000 * product_index + 10 * fold
            dem_seed = 20260815 + Int(sum(codeunits(product))) +
                (fold == 0 ? 0 : 10_000 + 10 * fold)
            selection = select_joint_covariates(
                product, Yobs, Ysat, terrain, era5, ndvi_aligned,
                train_indices, station_ids, times, lonlat, cfg, run_seed, dem_seed,
            )
            candidates = selection.candidates
            _annotate!(candidates, product, scheme, fold); _append!(candidate_all, candidates)
            if isempty(selection.active)
                qc = selection.qc
                _annotate!(qc, product, scheme, fold); _append!(qc_all, qc)
                push!(status_all, (product=product, scheme=scheme, fold=fold,
                    status="ok", reason="no_independent_candidates"))
                for row in eachrow(candidates)
                    push!(role_all, (product=product, scheme=scheme, fold=fold,
                        family=row.family, variable_group=row.variable_group,
                        independent_selected=false, joint_vif_retained=false,
                        role="not_selected", final_included=false,
                        exclusion_reason=row.independent_status))
                    if scheme == "full_data"
                        push!(final_spec, (product, row.family, row.variable_group,
                            row.predictor_columns, false, row.independent_qvalue, false,
                            "not_selected", false, row.independent_status))
                    end
                end
                continue
            end
            qc = copy(selection.qc); _annotate!(qc, product, scheme, fold); _append!(qc_all, qc)
            scales = copy(selection.scales)
            _annotate!(scales, product, scheme, fold); _append!(scale_all, scales)
            weight_audit = copy(selection.weight_audit)
            weight_audit.station_id = station_ids[weight_audit.global_station_index]
            _annotate!(weight_audit, product, scheme, fold); _append!(weight_all, weight_audit)
            vif = selection.vif; retained = selection.retained
            _annotate!(vif, product, scheme, fold); _append!(vif_all, vif)
            _annotate!(selection.spatial_scan, product, scheme, fold)
            _append!(bandwidth_all, selection.spatial_scan)
            _annotate!(selection.spatial_result, product, scheme, fold)
            _append!(spatial_all, selection.spatial_result)
            role_status_map = Dict(String(row.variable_group) => (role=String(row.role), status=String(row.status))
                for row in eachrow(selection.spatial_result))
            for row in eachrow(candidates)
                group = String(row.variable_group)
                selected = Bool(row.independent_selected)
                vif_retained = group in retained
                if !selected
                    role, included, reason = "not_selected", false, String(row.independent_status)
                elseif !vif_retained
                    role, included, reason = "not_selected", false,
                        selection.status == "ok" ? "removed_by_joint_vif" : selection.status_reason
                else
                    info = get(role_status_map, group, (role="uncertain", status="test_failed"))
                    role = info.role; included = role in ("local", "global")
                    reason = included ? "" : info.status
                end
                push!(role_all, (product=product, scheme=scheme, fold=fold,
                    family=row.family, variable_group=group, independent_selected=selected,
                    joint_vif_retained=vif_retained, role, final_included=included,
                    exclusion_reason=reason))
                if scheme == "full_data"
                    push!(final_spec, (product, String(row.family), group,
                        String(row.predictor_columns), selected, Float64(row.independent_qvalue),
                        vif_retained, role, included, reason))
                end
            end
            push!(status_all, (product=product, scheme=scheme, fold=fold,
                status=selection.status, reason=selection.status_reason))
        end
    end
    stability = DataFrame(
        product=String[], family=String[], variable_group=String[], fold_count=Int[],
        independent_selected_count=Int[], vif_retained_count=Int[],
        local_count=Int[], global_count=Int[], uncertain_count=Int[], final_included_count=Int[],
    )
    for group in groupby(filter(:scheme => ==("spatial_cv"), role_all),
        [:product, :family, :variable_group])
        push!(stability, (
            first(group.product), first(group.family), first(group.variable_group), nrow(group),
            count(group.independent_selected), count(group.joint_vif_retained),
            count(==("local"), group.role), count(==("global"), group.role),
            count(==("uncertain"), group.role), count(group.final_included),
        ))
    end
    consensus = DataFrame(
        family=String[], variable_group=String[], included_product_count=Int[],
        local_product_count=Int[], global_product_count=Int[],
        consensus_included=Bool[], consensus_role=String[],
    )
    for variable_group in JOINT_GROUP_ORDER
        rows = filter(:variable_group => ==(variable_group), final_spec)
        included_count = count(rows.final_included)
        local_count = count(==("local"), rows.role); global_count = count(==("global"), rows.role)
        role = included_count == 0 ? "not_selected" : local_count >= 2 ? "local" :
            global_count >= 2 ? "global" : included_count >= 2 ? "role_unstable" : "product_specific"
        push!(consensus, (GROUP_FAMILY[variable_group], variable_group, included_count,
            local_count, global_count, included_count >= 2, role))
    end
    outputs = Dict(
        "joint_independent_candidates.csv" => candidate_all,
        "joint_panel_quality_control.csv" => qc_all,
        "joint_scaling.csv" => scale_all,
        "joint_weight_audit.csv" => weight_all,
        "joint_vif.csv" => vif_all,
        "joint_bandwidth_scan.csv" => bandwidth_all,
        "joint_spatial_variability.csv" => spatial_all,
        "joint_fold_roles.csv" => role_all,
        "joint_role_stability.csv" => stability,
        "joint_final_full_data_spec.csv" => final_spec,
        "joint_cross_product_consensus.csv" => consensus,
        "run_status.csv" => status_all,
    )
    for (filename, table) in outputs
        CSV.write(joinpath(cfg.outdir, filename), table)
    end
    return (; final_spec, stability, consensus, status=status_all)
end

end
