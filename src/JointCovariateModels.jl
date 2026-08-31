module JointCovariateModels

using CSV
using DataFrames
using Dates
using LinearAlgebra
using SHA
using Statistics

# Loaded through src/load_modules.jl by whoever loads this file. Previously included here,
# which compiled a second DEMTerrainExperiment whose types were incompatible with the one
# the benchmark and JointVariableSelection used.
# `using Main.X: a, b` imports only `a` and `b`; it does not bind `X` itself, and this module
# also reaches DEMTerrainExperiment by qualification (`_haversine_matrix`, `_local_hat`,
# `_global_projection`). Bind the module name explicitly so both forms resolve, without
# pulling in the rest of its exported surface.
using Main: DEMTerrainExperiment
using Main.DEMTerrainExperiment: mixed_gwr_predict, multiscale_gwr_predict

export JointCovariateBenchmarkConfig, JointFoldContext
export load_joint_covariate_spec, joint_spec_sha256, build_joint_fold_context
export dynamic_covariate_predict, joint_effective_roles, joint_group_names
export joint_models_coincide

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
const GROUP_FAMILY = Dict(group =>
    group in ("elevation", "slope", "aspect") ? "dem" :
    group == "ndvi" ? "ndvi" : "era5" for group in JOINT_GROUP_ORDER)

Base.@kwdef struct JointCovariateBenchmarkConfig
    # `nothing` when the benchmark is configured for nested per-fold selection instead of a
    # fixed full-data specification (see `InterpolationBenchmarkConfig.joint_selection`).
    spec_path::Union{Nothing,String} = nothing
    terrain_path::String
    era5_paths::Dict{Int,String}
    ndvi_path::Union{Nothing,String} = nothing
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    feature_time_offset_hours::Int = 9
    max_ndvi_age_days::Int = 32
    wet_threshold::Float64 = 0.1
    ridge::Float64 = 1e-8
    tolerance::Float64 = 1e-5
    # The back-fit routinely needs several hundred sweeps once smooth covariate groups are in
    # play; at the old cap of 200 the median sweep count sat at 150-194 and MGWR discarded whole
    # hours for want of a few more iterations. Measured on the real `_local_hat`/`_global_projection`
    # operators (150 stations, 93%-dry residuals, two smooth ERA5-like covariate groups): the
    # `:split` layout converges on 60% of hours at 200 sweeps and 100% at 500.
    max_iterations::Int = 1000
    # Successive over-relaxation factor for the back-fitting sweeps. It converges to the same
    # fixed point for any admissible value - `relaxation` cancels at the fixed point - so this
    # trades convergence *rate* only.
    #
    # DO NOT "fix" this to 1.0. The name `_multiscale_predict_damped` suggests damping and a
    # value above 1 looks like a bug, but plain Gauss-Seidel is far worse here: on the fixture
    # above, the share of hours converging within 200 sweeps runs
    #
    #   layout            w=1.0  w=1.1  w=1.2  w=1.3  w=1.5  w=1.7
    #   :split               4%    12%    40%    98%    60%     1%
    #   :shared              3%     8%    22%   100%   100%    52%
    #   :intercept_only    100%   100%   100%   100%   100%   100%
    #
    # i.e. over-relaxation is what makes the coupled layouts usable at all. 1.5 is kept as the
    # historical value; 1.3 looked better on that one synthetic fixture, which is not enough
    # evidence to move a default, and the layout that wins on accuracy is insensitive to it.
    relaxation::Float64 = 1.5
    # How MGWR splits the spatial trend into back-fitting groups. Covariate groups always get
    # their own bandwidth, which is what makes the model multiscale; this controls only the
    # `[1, z_lon, z_lat]` block.
    #
    #   :split           - one single-column group each for the intercept, longitude and latitude
    #                      (the historical layout).
    #   :shared          - one three-column group, solved by a single weighted least squares the
    #                      way `mixed_gwr` solves its spatial block. That is a statement about the
    #                      block, not about the model: the covariates stay in their own groups, so
    #                      `:shared` is not `mixed_gwr` with extra bandwidths either.
    #   :intercept_only  - drop the coordinate columns; a locally varying intercept plus the
    #                      covariate groups.
    #
    # The coordinate groups are near-constant inside any local window, so they are near-collinear
    # with the intercept group *and* with any spatially smooth covariate. That coupling, not the
    # relaxation factor, is what governs how many sweeps the back-fit needs.
    #
    # `:intercept_only` is the default on measured evidence. Against a common-cell paired
    # comparison of the smoke run (`scripts/compare_benchmark_runs.jl --smoke`, where the
    # traditional methods move +0.0% and so confirm the pairing), MGWR's RMSE change is
    #
    #   layout            FY4B     GPM   GSMaP   coverage at 1.000   smoke wall clock
    #   :split           +0.1%   -6.8%   -5.1%             8 / 15             30-58 min
    #   :shared          +0.9%   -5.7%   -3.5%            13 / 15                16 min
    #   :intercept_only  -4.3%   -6.5%   -6.4%            15 / 15                14 min
    #
    # i.e. it wins on two products, is 0.3pp behind on the third against a layout scored on
    # fewer cells, and is the only one that never discards an hour. Dropping the coordinate plane
    # also frees the intercept to go genuinely local (8-80 neighbours) instead of pinning near the
    # grid maximum, which is the multiscale behaviour the method is supposed to have.
    #
    # What none of the layouts is, is `mixed_gwr` with the single-bandwidth constraint relaxed.
    # `:intercept_only` cannot be: it has two fewer design columns at the same group count, so no
    # bandwidth vector recovers `mixed_gwr`. Even `:shared`, which carries the identical columns,
    # only approximates it — driving every group to the same bandwidth and comparing against
    # `mixed_gwr` gives machine precision (7.8e-16) when there are no local covariates and a
    # single group on each side, but with covariate groups the difference plateaus under a
    # tightening tolerance rather than vanishing: 1.6e-7 relative, worst case 1.3e-5 at the
    # narrowest bandwidth, falling to 4e-10 by bw=24. `_local_hat` is a local-linear smoother and
    # is not idempotent, so back-fitting converges to an additive fixed point rather than to the
    # joint weighted least squares `mixed_gwr` solves; the two agree only in the one-group case
    # and as the bandwidth widens toward a projection. Four orders below reporting precision, so
    # `:shared` is a usable nested comparison in practice - but it is not an identity, and
    # `:intercept_only` is not nested at all.
    mgwr_spatial_grouping::Symbol = :intercept_only
    # How a local group's fit at a target with too few supporting stations is recorded.
    #
    # `:missing` marks that target NaN, so the coverage machinery (`min_tuning_coverage`,
    # `_prediction_coverage`, `_common_method_mask`) can see it and reject the candidate.
    # `:zero` is the historical all-zero hat row, which reports a fabricated "no local
    # correction" as a successful prediction and is invisible to every one of those gates:
    # measured on the full nested run, the gate rejected 195/675 of direct `gwr`'s candidates
    # and 0/675, 0/675, 0/3939 of `residual_gwr`/`mixed_gwr`/`mgwr`'s on the identical grid,
    # and `mgwr`'s reported GPM own-coverage of 0.9985 was 0.793 once the unfitted targets
    # were counted honestly. Only a fixed-km bandwidth with a compact kernel
    # (bisquare/tricube/boxcar) can reach the branch at all.
    #
    # `:zero` is kept so a run can reproduce the pre-fix numbers exactly.
    unsupported_local_target::Symbol = :missing
end

struct JointFoldContext
    product::String
    variables::Vector{String}
    roles::Dict{String,String}
    spatial_train::Matrix{Float64}
    spatial_target::Matrix{Float64}
    predictor_train::Dict{String,Array{Float64,3}}
    predictor_target::Dict{String,Array{Float64,3}}
    train_lonlat::Matrix{Float64}
    target_lonlat::Matrix{Float64}
    bandwidth_candidates::Vector{Int}
    scaling::DataFrame
    quality_control::DataFrame
    config::JointCovariateBenchmarkConfig
end

function joint_spec_sha256(path::String)
    isfile(path) || throw(ArgumentError("joint variable specification does not exist: $path"))
    return bytes2hex(open(sha256, path))
end

function load_joint_covariate_spec(path::String, products::Vector{String})
    isfile(path) || throw(ArgumentError("joint variable specification does not exist: $path"))
    table = CSV.read(path, DataFrame)
    required = [
        :product, :family, :variable_group, :predictor_columns,
        :independent_selected, :joint_vif_retained, :role, :final_included,
    ]
    missing_columns = setdiff(required, propertynames(table))
    isempty(missing_columns) || throw(ArgumentError(
        "joint variable specification is missing columns: $(join(string.(missing_columns), ", "))",
    ))
    allunique(zip(String.(table.product), String.(table.variable_group))) ||
        throw(ArgumentError("joint variable specification contains duplicate product/group rows"))
    Set(String.(table.product)) == Set(products) ||
        throw(ArgumentError("joint variable specification products differ from benchmark products"))
    all(in(Set(JOINT_GROUP_ORDER)), String.(table.variable_group)) ||
        throw(ArgumentError("joint variable specification contains an unsupported variable group"))
    included = filter(:final_included => identity, table)
    for row in eachrow(included)
        Bool(row.independent_selected) && Bool(row.joint_vif_retained) ||
            throw(ArgumentError("included variable bypassed independent screening or joint VIF"))
        String(row.role) in ("local", "global") ||
            throw(ArgumentError("included variable has unresolved role"))
        expected_family = GROUP_FAMILY[String(row.variable_group)]
        String(row.family) == expected_family ||
            throw(ArgumentError("variable family differs for $(row.variable_group)"))
    end
    role_maps = Dict{String,Dict{String,String}}()
    for product in products
        rows = filter([:product, :final_included] =>
            (p, included) -> String(p) == product && included, table)
        role_maps[product] = Dict(String(row.variable_group) => String(row.role)
            for row in eachrow(rows))
    end
    return (; table, included, role_maps)
end

function _raw_group(
    group::String, indices::Vector{Int}, nt::Int, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned,
)
    columns = GROUP_COLUMNS[group]
    values = Array{Float64}(undef, length(indices), length(columns), nt)
    for (column_index, column) in enumerate(columns)
        if GROUP_FAMILY[group] == "dem"
            values[:, column_index, :] .= reshape(
                Float64.(terrain[indices, column]), length(indices), 1,
            )
        elseif GROUP_FAMILY[group] == "era5"
            values[:, column_index, :] .= Float64.(era5[column][indices, :])
        else
            ndvi_aligned === nothing && throw(ArgumentError("NDVI is selected but not loaded"))
            values[:, column_index, :] .= Float64.(ndvi_aligned.values[indices, :])
        end
    end
    return values
end

function _spatial_design(train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64})
    center = vec(mean(train_lonlat, dims=1))
    scale = vec(std(train_lonlat, dims=1))
    all(isfinite, scale) && all(>(0), scale) ||
        throw(ArgumentError("training coordinates have insufficient variation"))
    train = hcat(ones(size(train_lonlat, 1)),
        (train_lonlat[:, 1] .- center[1]) ./ scale[1],
        (train_lonlat[:, 2] .- center[2]) ./ scale[2])
    target = hcat(ones(size(target_lonlat, 1)),
        (target_lonlat[:, 1] .- center[1]) ./ scale[1],
        (target_lonlat[:, 2] .- center[2]) ./ scale[2])
    return train, target
end

function _hierarchical_weights(mask::BitMatrix, train_indices::Vector{Int}, ndvi_aligned)
    n, nt = size(mask)
    weights = zeros(n, nt)
    uses_ndvi = ndvi_aligned !== nothing
    for i in 1:n
        valid_times = findall(mask[i, :])
        isempty(valid_times) && continue
        if uses_ndvi
            periods = ndvi_aligned.source_period[train_indices[i], valid_times]
            valid_periods = sort(unique(filter(>(0), periods)))
            isempty(valid_periods) && continue
            for period in valid_periods
                period_times = valid_times[periods .== period]
                weights[i, period_times] .= 1 / (length(valid_periods) * length(period_times))
            end
        else
            weights[i, valid_times] .= 1 / length(valid_times)
        end
    end
    return weights
end

function build_joint_fold_context(
    product::String, roles::Dict{String,String}, train_indices::Vector{Int},
    target_indices::Vector{Int}, lonlat::Matrix{Float64}, Yobs::Matrix{Float64},
    Ysat::Matrix{Float64}, terrain::DataFrame, era5::AbstractDict, ndvi_aligned,
    cfg::JointCovariateBenchmarkConfig,
)
    variables = [group for group in JOINT_GROUP_ORDER if haskey(roles, group)]
    nt = size(Yobs, 2)
    train_lonlat = Matrix{Float64}(lonlat[train_indices, :])
    target_lonlat = Matrix{Float64}(lonlat[target_indices, :])
    spatial_train, spatial_target = _spatial_design(train_lonlat, target_lonlat)
    raw_train = Dict(group => _raw_group(
        group, train_indices, nt, terrain, era5,
        group == "ndvi" ? ndvi_aligned : nothing,
    ) for group in variables)
    raw_target = Dict(group => _raw_group(
        group, target_indices, nt, terrain, era5,
        group == "ndvi" ? ndvi_aligned : nothing,
    ) for group in variables)

    mask = isfinite.(Yobs[train_indices, :]) .& isfinite.(Ysat[train_indices, :]) .&
        (Yobs[train_indices, :] .>= cfg.wet_threshold)
    for group in variables, column in axes(raw_train[group], 2)
        mask .&= isfinite.(raw_train[group][:, column, :])
    end
    weights = _hierarchical_weights(
        mask, train_indices, "ndvi" in variables ? ndvi_aligned : nothing,
    )
    predictor_train = Dict{String,Array{Float64,3}}()
    predictor_target = Dict{String,Array{Float64,3}}()
    scaling_rows = NamedTuple[]
    for group in variables
        train = copy(raw_train[group]); target = copy(raw_target[group])
        for column in axes(train, 2)
            for time in 1:nt
                finite_train = isfinite.(@view train[:, column, time])
                if any(finite_train)
                    hour_mean = mean(@view train[finite_train, column, time])
                    train[finite_train, column, time] .-= hour_mean
                    finite_target = isfinite.(@view target[:, column, time])
                    target[finite_target, column, time] .-= hour_mean
                end
            end
            weighted_values = Float64[]; value_weights = Float64[]
            for i in axes(train, 1), time in 1:nt
                weight = weights[i, time]
                weight > 0 && isfinite(train[i, column, time]) || continue
                push!(weighted_values, train[i, column, time])
                push!(value_weights, weight)
            end
            total = sum(value_weights)
            center = total > 0 ? dot(weighted_values, value_weights) / total : NaN
            scale = total > 0 ? sqrt(dot(
                (weighted_values .- center) .^ 2, value_weights,
            ) / total) : NaN
            isfinite(scale) && scale > 0 || throw(ArgumentError(
                "training-fold scale is unavailable for $product/$group column $column",
            ))
            for index in eachindex(train)
                isfinite(train[index]) && (train[index] = (train[index] - center) / scale)
            end
            for index in eachindex(target)
                isfinite(target[index]) && (target[index] = (target[index] - center) / scale)
            end
            push!(scaling_rows, (;
                product, variable_group=group,
                predictor_column=String(GROUP_COLUMNS[group][column]),
                center, scale, wet_station_hours=count(weights .> 0),
            ))
        end
        predictor_train[group] = train
        predictor_target[group] = target
    end
    valid_station_weights = vec(sum(weights, dims=2))
    positive = valid_station_weights[valid_station_weights .> 0]
    qc = DataFrame([(
        product=product, training_station_count=length(train_indices),
        target_station_count=length(target_indices), time_count=nt,
        variable_count=length(variables), variables=join(variables, ","),
        roles=join(["$group=$(roles[group])" for group in variables], ";"),
        wet_station_hours=count(mask),
        minimum_station_weight=isempty(positive) ? NaN : minimum(positive),
        maximum_station_weight=isempty(positive) ? NaN : maximum(positive),
        status="ok",
    )])
    return JointFoldContext(
        product, variables, copy(roles), spatial_train, spatial_target,
        predictor_train, predictor_target, train_lonlat, target_lonlat,
        filter(<(length(train_indices)), unique(cfg.bandwidth_candidates)),
        DataFrame(scaling_rows), qc, cfg,
    )
end

function joint_effective_roles(context::JointFoldContext, method::String)
    method == "residual_gwr" && return Dict(group => "local" for group in context.variables)
    method in ("mixed_gwr", "mgwr") || throw(ArgumentError("unsupported joint model $method"))
    return copy(context.roles)
end

"""
Whether two joint models reduce to the same fit on this fold's covariate roles.

`residual_gwr` forces every covariate local; `mixed_gwr` honours the fold's role map. When that
map contains no `"global"` the two build the identical design — `_design_at` then differs only in
the group *name* (`"shared_all_local"` vs `"shared_local"`) — and they search the identical grid
through the identical scorer, so they are one model reported under two names. On the full nested
run that held at 19 of 30 fold-cells, with bit-identical selected kernel, bandwidth, shrinkage and
inner RMSE; `joint_role_stability.csv` shows why, a `"global"` role appearing in only 6 of 54
group-cells.

`mgwr` is never a duplicate of either: it gives every group its own bandwidth whatever the roles
say, and under `:intercept_only` it also drops the coordinate columns.

This is what `run_status.csv`'s `duplicate_of` column and `select_auto_method`'s contender
collapsing are keyed on, so a reader does not have to infer "same model" from matching RMSEs.
"""
function joint_models_coincide(
    context::JointFoldContext, first_method::AbstractString, second_method::AbstractString,
)
    first_method == second_method && return true
    ("mgwr" in (first_method, second_method)) && return false
    return joint_effective_roles(context, String(first_method)) ==
        joint_effective_roles(context, String(second_method))
end

function joint_group_names(context::JointFoldContext, method::String)
    design = _design_at(context, method, 1; target=false)
    return design.group_names
end

function _design_at(context::JointFoldContext, method::String, time::Int; target::Bool=false)
    spatial = target ? context.spatial_target : context.spatial_train
    predictors = target ? context.predictor_target : context.predictor_train
    roles = joint_effective_roles(context, method)
    local_groups = Matrix{Float64}[]
    group_names = String[]
    if method == "mgwr"
        # `mgwr_spatial_grouping` controls only the spatial block; covariates always get their
        # own group, and therefore their own bandwidth, in every variant.
        grouping = context.config.mgwr_spatial_grouping
        if grouping === :split
            for index in axes(spatial, 2)
                push!(local_groups, spatial[:, index:index])
                push!(group_names, ("intercept", "longitude", "latitude")[index])
            end
        elseif grouping === :shared
            push!(local_groups, copy(spatial))
            # Deliberately not "shared_local" (mixed_gwr's name) or "all" (`_scan_row`'s default
            # for ungrouped methods): the scan-row lookups match on the group name.
            push!(group_names, "shared_spatial")
        elseif grouping === :intercept_only
            push!(local_groups, spatial[:, 1:1])
            push!(group_names, "intercept")
        else
            throw(ArgumentError("unsupported mgwr_spatial_grouping: $grouping"))
        end
        for group in context.variables
            roles[group] == "local" || continue
            push!(local_groups, Matrix{Float64}(predictors[group][:, :, time]))
            push!(group_names, group)
        end
    else
        local_design = copy(spatial)
        for group in context.variables
            roles[group] == "local" || continue
            local_design = hcat(local_design, predictors[group][:, :, time])
        end
        push!(local_groups, local_design)
        push!(group_names, method == "residual_gwr" ? "shared_all_local" : "shared_local")
    end
    global_design = zeros(Float64, size(spatial, 1), 0)
    for group in context.variables
        roles[group] == "global" || continue
        global_design = hcat(global_design, predictors[group][:, :, time])
    end
    return (; local_groups, global_design, group_names)
end

function _haversine_matrix(train::Matrix{Float64}, target::Matrix{Float64})
    radius = 6371.0088
    result = Matrix{Float64}(undef, size(train, 1), size(target, 1))
    for i in axes(train, 1), j in axes(target, 1)
        lat1, lat2 = deg2rad(train[i, 2]), deg2rad(target[j, 2])
        dlat = lat2 - lat1
        dlon = deg2rad(target[j, 1] - train[i, 1])
        value = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
        result[i, j] = 2radius * asin(sqrt(clamp(value, 0.0, 1.0)))
    end
    return result
end

function _adaptive_weights(distances::AbstractVector{Float64}, neighbors::Int)
    finite = sort(filter(isfinite, distances))
    isempty(finite) && return zeros(length(distances))
    k = clamp(neighbors, 1, length(finite))
    bandwidth = finite[k]
    bandwidth <= 0 && (bandwidth = maximum(finite[finite .> 0]; init=1.0))
    weights = zeros(length(distances))
    for i in eachindex(distances)
        ratio = distances[i] / bandwidth
        isfinite(ratio) && ratio < 1 && (weights[i] = (1 - ratio^2)^2)
    end
    return weights
end

function _weight_matrix(
    distances::Matrix{Float64}, neighbors::Int; exclude_self::Bool=false,
)
    weights = zeros(size(distances, 2), size(distances, 1))
    for target in axes(distances, 2)
        d = copy(@view distances[:, target])
        exclude_self && target <= length(d) && (d[target] = Inf)
        weights[target, :] .= _adaptive_weights(d, neighbors)
    end
    return weights
end

function _local_predict(
    Xtrain::Matrix{Float64}, y::Vector{Float64}, Xtarget::Matrix{Float64},
    weights::Matrix{Float64}, ridge::Float64,
)
    prediction = fill(NaN, size(Xtarget, 1))
    p = size(Xtrain, 2)
    for target in axes(Xtarget, 1)
        w = @view weights[target, :]
        valid = (w .> 0) .& isfinite.(y) .& vec(all(isfinite, Xtrain; dims=2))
        count(valid) >= p + 1 || continue
        X = Xtrain[valid, :]; values = y[valid]; selected_weights = w[valid]
        beta = try
            (X' * (selected_weights .* X) + ridge * I) \
                (X' * (selected_weights .* values))
        catch
            continue
        end
        all(isfinite, @view Xtarget[target, :]) &&
            (prediction[target] = dot(@view(Xtarget[target, :]), beta))
    end
    return prediction
end

function _global_predict(
    Xtrain::Matrix{Float64}, y::Vector{Float64}, Xtarget::Matrix{Float64},
    ridge::Float64; leave_one_out::Bool=false,
)
    isempty(Xtrain) && return zeros(size(Xtarget, 1))
    if !leave_one_out
        beta = (Xtrain' * Xtrain + ridge * I) \ (Xtrain' * y)
        return Xtarget * beta
    end
    size(Xtrain, 1) == size(Xtarget, 1) ||
        throw(DimensionMismatch("global leave-one-out requires matching rows"))
    prediction = fill(NaN, size(Xtrain, 1))
    for target in axes(Xtrain, 1)
        keep = collect(axes(Xtrain, 1)) .!= target
        X = Xtrain[keep, :]; values = y[keep]
        beta = try (X' * X + ridge * I) \ (X' * values) catch; continue end
        prediction[target] = dot(@view(Xtrain[target, :]), beta)
    end
    return prediction
end

function _predict_time(
    local_train::Vector{Matrix{Float64}}, global_train::Matrix{Float64},
    y::Vector{Float64}, train_lonlat::Matrix{Float64},
    local_target::Vector{Matrix{Float64}}, global_target::Matrix{Float64},
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Int}, cfg;
    leave_one_out::Bool=false,
)
    length(local_train) == length(local_target) == length(bandwidths) ||
        throw(DimensionMismatch("local group dimensions differ"))
    train_distances = _haversine_matrix(train_lonlat, train_lonlat)
    target_distances = leave_one_out ? train_distances :
        _haversine_matrix(train_lonlat, target_lonlat)
    train_weights = [_weight_matrix(train_distances, bandwidth;
        exclude_self=leave_one_out) for bandwidth in bandwidths]
    target_weights = leave_one_out ? train_weights :
        [_weight_matrix(target_distances, bandwidth) for bandwidth in bandwidths]
    local_components = [zeros(length(y)) for _ in local_train]
    global_component = _global_predict(
        global_train, y, global_train, cfg.ridge; leave_one_out,
    )
    previous_fitted = fill(NaN, length(y))
    converged = false
    for _ in 1:cfg.max_iterations
        for group in eachindex(local_train)
            partial = y - global_component
            for other in eachindex(local_components)
                other == group || (partial .-= local_components[other])
            end
            local_components[group] = _local_predict(
                local_train[group], partial, local_train[group],
                train_weights[group], cfg.ridge,
            )
        end
        local_sum = isempty(local_components) ? zeros(length(y)) : reduce(+, local_components)
        global_component = _global_predict(
            global_train, y - local_sum, global_train, cfg.ridge; leave_one_out,
        )
        fitted = local_sum + global_component
        all(isfinite, fitted) || return fill(NaN, size(target_lonlat, 1)), false
        change = all(isfinite, previous_fitted) ?
            norm(fitted - previous_fitted) / max(norm(previous_fitted), eps()) : Inf
        if change < cfg.tolerance
            converged = true
            break
        end
        previous_fitted .= fitted
    end
    converged || return fill(NaN, size(target_lonlat, 1)), false
    leave_one_out && return reduce(+, local_components; init=global_component), true
    local_sum = reduce(+, local_components; init=zeros(length(y)))
    global_target_component = _global_predict(
        global_train, y - local_sum, global_target, cfg.ridge,
    )
    prediction = global_target_component
    for group in eachindex(local_train)
        partial = y - global_component
        for other in eachindex(local_components)
            other == group || (partial .-= local_components[other])
        end
        prediction .+= _local_predict(
            local_train[group], partial, local_target[group],
            target_weights[group], cfg.ridge,
        )
    end
    return prediction, true
end

function _multiscale_predict_damped(
    local_train::Vector{Matrix{Float64}}, global_train::Matrix{Float64},
    response::Matrix{Float64}, train_lonlat::Matrix{Float64},
    local_target::Vector{Matrix{Float64}}, global_target::Matrix{Float64},
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Float64}, kernel::Function, cfg;
    adaptive::Bool=true, exclude_self::Bool=false, unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
)
    length(local_train) == length(local_target) == length(bandwidths) ||
        throw(DimensionMismatch(
            "local group count ($(length(local_train)) train, $(length(local_target)) target) " *
            "must match the bandwidth count ($(length(bandwidths)))",
        ))
    y = vec(response)
    # Supplied by `dynamic_covariate_predict`, which holds the station geometry fixed across every
    # hour it fits and so builds these once instead of once per hour. Same contract as
    # `mixed_gwr_predict`: the distances of `train_lonlat`/`target_lonlat` exactly as passed.
    train_d::Matrix{Float64} = train_distances === nothing ?
        DEMTerrainExperiment._haversine_matrix(train_lonlat, train_lonlat) : train_distances
    train_hats = [DEMTerrainExperiment._local_hat(
        X, X, train_d, bandwidth, kernel; adaptive, ridge=cfg.ridge, exclude_self,
    ) for (X, bandwidth) in zip(local_train, bandwidths)]
    global_hat = DEMTerrainExperiment._global_projection(global_train; ridge=cfg.ridge)
    local_components = [zeros(length(y)) for _ in local_train]
    global_component = global_hat * y
    previous = fill(NaN, length(y))
    converged = false
    relaxation = cfg.relaxation
    # Scratch for the sweep loop, allocated once instead of once per group per sweep. The sweep
    # count runs to the hundreds and this loop is the benchmark's innermost hot spot, so the
    # allocations were thousands of short-lived vectors per hour-fit, inside a threaded region
    # where the resulting GC pauses stop every thread at once.
    #
    # The arithmetic is deliberately left in exactly its original order and grouping, including
    # `component_sum` starting from zero (which is what `reduce(+, ...; init=zeros(...))` did) and
    # the `norm` of a materialised difference, so every value is bit-for-bit unchanged.
    n_obs = length(y)
    partial = Vector{Float64}(undef, n_obs)
    update = Vector{Float64}(undef, n_obs)
    component_sum = Vector{Float64}(undef, n_obs)
    global_update = Vector{Float64}(undef, n_obs)
    fitted = Vector{Float64}(undef, n_obs)
    difference = Vector{Float64}(undef, n_obs)
    for _ in 1:cfg.max_iterations
        for group in eachindex(local_train)
            partial .= y .- global_component
            for other in eachindex(local_components)
                other == group || (partial .-= local_components[other])
            end
            mul!(update, train_hats[group], partial)
            local_components[group] .= relaxation .* update .+
                (1 - relaxation) .* local_components[group]
        end
        fill!(component_sum, 0.0)
        for component in local_components
            component_sum .+= component
        end
        partial .= y .- component_sum
        mul!(global_update, global_hat, partial)
        global_component .= relaxation .* global_update .+
            (1 - relaxation) .* global_component
        fitted .= component_sum .+ global_component
        # A diverging sweep makes `change` NaN, which fails the tolerance test silently and then
        # burns the whole iteration budget. Bail as soon as the iterate leaves the reals, the way
        # the undamped sibling does.
        all(isfinite, fitted) || return fill(NaN, size(target_lonlat, 1), 1), falses(1)
        change = if all(isfinite, previous)
            difference .= fitted .- previous
            norm(difference) / max(norm(previous), eps())
        else
            Inf
        end
        previous .= fitted
        if change < cfg.tolerance
            converged = true
            break
        end
    end
    converged || return fill(NaN, size(target_lonlat, 1), 1), falses(1)
    if exclude_self
        # `unsupported` cannot be honoured here. Unlike `mixed_gwr_predict` and
        # `multiscale_gwr_predict`, this branch has no separate target hat to mark: the
        # leave-one-out prediction *is* the fitted vector of the self-excluding `train_hats`
        # above, and a NaN row in those would propagate into `fitted`, trip the
        # `all(isfinite, fitted)` guard and discard the whole hour rather than one station.
        # Only `tuning_geometry = :loocv` and the no-inner-split joint scoring path reach it;
        # the default `:inner_spatial` never does.
        return reshape(previous, :, 1), trues(1)
    end
    local_sum = reduce(+, local_components; init=zeros(length(y)))
    global_beta = isempty(global_train) ? Float64[] :
        (global_train' * global_train + cfg.ridge * I) \
            (global_train' * (y - local_sum))
    prediction = isempty(global_target) ? zeros(size(target_lonlat, 1)) :
        global_target * global_beta
    target_d::Matrix{Float64} = target_distances === nothing ?
        DEMTerrainExperiment._haversine_matrix(train_lonlat, target_lonlat) : target_distances
    for group in eachindex(local_train)
        partial = y - (isempty(global_train) ? zeros(length(y)) : global_train * global_beta)
        for other in eachindex(local_components)
            other == group || (partial .-= local_components[other])
        end
        target_hat = DEMTerrainExperiment._local_hat(
            local_train[group], local_target[group], target_d,
            bandwidths[group], kernel; adaptive, ridge=cfg.ridge, unsupported,
        )
        prediction .+= target_hat * partial
    end
    return reshape(prediction, :, 1), trues(1)
end

function dynamic_covariate_predict(
    context::JointFoldContext, residuals::Matrix{Float64}, method::String,
    bandwidths::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, time_indices::Vector{Int}=collect(axes(residuals, 2)),
    leave_one_out::Bool=false,
)
    method in ("residual_gwr", "mixed_gwr", "mgwr") ||
        throw(ArgumentError("unsupported joint model $method"))
    # Forwarded to every *target* hat below. It is a documented no-op on
    # `_multiscale_predict_damped`'s `exclude_self` path, which has no separate target hat;
    # see the comment at that early return.
    unsupported = context.config.unsupported_local_target
    target_count = leave_one_out ? size(residuals, 1) : size(context.target_lonlat, 1)
    prediction = fill(NaN, target_count, length(time_indices))
    # `Vector{Bool}`, not a `BitVector`. The loop below writes one entry per iteration from
    # many threads, and a BitVector packs 64 entries per word, so those writes would be
    # read-modify-writes of a shared word and would drop each other. That race predates the
    # scheduling change below but was mostly hidden by contiguous per-thread chunks; interleaved
    # scheduling would make it routine. One byte per entry makes each write independent.
    converged = fill(false, length(time_indices))
    # The station geometry is a property of the fold context, not of the hour: only which rows are
    # `valid` changes below. Building the two haversine matrices here rather than inside the loop
    # replaces ~2*N^2 transcendental evaluations per hour with one row/column gather. The
    # submatrix of a distance matrix is elementwise the distance matrix of the submatrix, so the
    # values handed to `_local_hat` are unchanged.
    full_train_distances =
        DEMTerrainExperiment._haversine_matrix(context.train_lonlat, context.train_lonlat)
    full_target_distances = leave_one_out ? full_train_distances :
        DEMTerrainExperiment._haversine_matrix(context.train_lonlat, context.target_lonlat)
    # `:greedy`, not the default chunked schedule. Hour cost here spans orders of magnitude - a
    # dry or underdetermined hour falls straight through the `count(valid) > min_required` guard
    # below, while a wet one runs hundreds of back-fitting sweeps - and precipitation is strongly
    # autocorrelated in time, so splitting the range into contiguous per-thread chunks drops whole
    # storms onto single threads. Iterations write only their own column of `prediction` and
    # `converged` and read nothing another iteration writes, so the schedule cannot change a value.
    Threads.@threads :greedy for output_time in eachindex(time_indices)
        time = time_indices[output_time]
        train_design = _design_at(context, method, time; target=false)
        target_design = leave_one_out ? train_design :
            _design_at(context, method, time; target=true)
        y = Vector{Float64}(residuals[:, time])
        valid = isfinite.(y)
        for X in train_design.local_groups
            valid .&= vec(all(isfinite, X; dims=2))
        end
        isempty(train_design.global_design) ||
            (valid .&= vec(all(isfinite, train_design.global_design; dims=2)))
        min_required = sum(size(X, 2) for X in train_design.local_groups) +
            size(train_design.global_design, 2) + 2
        count(valid) > min_required || continue
        local_train = [Matrix{Float64}(X[valid, :]) for X in train_design.local_groups]
        global_train = Matrix{Float64}(train_design.global_design[valid, :])
        train_lonlat = context.train_lonlat[valid, :]
        # Subset the fold's distance matrices the same way the coordinates are subset, so the
        # models below receive exactly the matrices they would have rebuilt for themselves - and
        # skip the copy entirely on the common hour where nothing is dropped. See
        # `DEMTerrainExperiment._distance_subset` for why these copies are worth avoiding.
        all_valid = all(valid)
        train_distances = all_valid ? full_train_distances : full_train_distances[valid, valid]
        if leave_one_out
            local_target = local_train
            global_target = global_train
            target_lonlat = train_lonlat
            target_distances = train_distances
        else
            local_target = target_design.local_groups
            global_target = target_design.global_design
            target_lonlat = context.target_lonlat
            target_distances = all_valid ? full_target_distances : full_target_distances[valid, :]
        end
        # Clamping to the valid station count only makes sense for an adaptive neighbor count; a
        # fixed-km bandwidth must be applied as given.
        adjusted = adaptive ? min.(bandwidths, count(valid) - 1) : bandwidths
        response = reshape(y[valid], :, 1)
        values_matrix, ok_vector = if leave_one_out && !isempty(global_train)
            # The held station is excluded analytically from the global coefficient;
            # local components are then tuned with a self-excluding spatial fit.
            global_loo = _global_predict(
                global_train, vec(response), global_train,
                context.config.ridge; leave_one_out=true,
            )
            partial = reshape(vec(response) - global_loo, :, 1)
            if method == "mgwr"
                local_prediction, ok = _multiscale_predict_damped(
                    local_train, zeros(Float64, count(valid), 0), partial,
                    train_lonlat, local_train, zeros(Float64, count(valid), 0),
                    train_lonlat, adjusted, kernel, context.config;
                    adaptive, exclude_self=true, unsupported,
                    train_distances, target_distances=train_distances,
                )
                local_prediction .+= global_loo
                local_prediction, ok
            else
                local_prediction, ok = mixed_gwr_predict(
                    only(local_train), zeros(Float64, count(valid), 0), partial,
                    train_lonlat, only(local_train), zeros(Float64, count(valid), 0),
                    train_lonlat, only(adjusted), kernel; adaptive, ridge=context.config.ridge,
                    tolerance=context.config.tolerance,
                    max_iterations=context.config.max_iterations,
                    exclude_self=true, unsupported,
                    train_distances, target_distances=train_distances,
                )
                local_prediction .+= global_loo
                local_prediction, ok
            end
        elseif method == "mgwr"
            _multiscale_predict_damped(
                local_train, global_train, response, train_lonlat,
                local_target, global_target, target_lonlat, adjusted, kernel, context.config;
                adaptive, exclude_self=leave_one_out, unsupported,
                train_distances, target_distances,
            )
        else
            mixed_gwr_predict(
                only(local_train), global_train, response, train_lonlat,
                only(local_target), global_target, target_lonlat, only(adjusted), kernel;
                adaptive, ridge=context.config.ridge, tolerance=context.config.tolerance,
                max_iterations=context.config.max_iterations,
                exclude_self=leave_one_out, unsupported,
                train_distances, target_distances,
            )
        end
        values = vec(values_matrix)
        ok = Bool(ok_vector[1])
        if leave_one_out
            prediction[findall(valid), output_time] = values
        else
            prediction[:, output_time] = values
        end
        converged[output_time] = ok
    end
    return prediction, converged
end

end
