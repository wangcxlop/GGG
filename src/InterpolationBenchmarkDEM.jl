function load_aligned_terrain(path::String, ids::Vector{String})
    terrain = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    required = [:station_id, :elevation_m, :slope_deg, :aspect_sin, :aspect_cos]
    missing_columns = setdiff(required, propertynames(terrain))
    isempty(missing_columns) || throw(ArgumentError(
        "terrain table is missing columns: $(join(string.(missing_columns), ", "))",
    ))
    allunique(terrain.station_id) || throw(ArgumentError("duplicate station IDs in terrain table"))
    row_by_id = Dict(String(id) => row for (row, id) in enumerate(terrain.station_id))
    missing_ids = filter(id -> !haskey(row_by_id, id), ids)
    isempty(missing_ids) || throw(ArgumentError(
        "terrain table is missing $(length(missing_ids)) benchmark stations",
    ))
    aligned = terrain[[row_by_id[id] for id in ids], :]
    values = Matrix{Float64}(aligned[:, required[2:end]])
    all(isfinite, values) || throw(ArgumentError("terrain predictors contain non-finite values"))
    return aligned
end

function _tag_dem_table!(
    table::DataFrame, scheme::String, product::String, fold::Int, phase::String,
)
    table[!, :scheme] = fill(scheme, nrow(table))
    :product in propertynames(table) || (table[!, :product] = fill(product, nrow(table)))
    table[!, :fold] = fill(fold, nrow(table))
    table[!, :phase] = fill(phase, nrow(table))
    select!(table, :scheme, :product, :fold, :phase, Not([:scheme, :product, :fold, :phase]))
    return table
end

function _failed_screen_table(product::String, error::String)
    return DataFrame([(
        product=product, variable_group=group, test=group == "aspect" ? "joint_F" : "pearson",
        statistic=NaN, pearson=NaN, spearman=NaN, pvalue=NaN, qvalue=NaN,
        direction_stable=false, selected_pre_vif=false, selected=false,
        exclusion_reason=error,
    ) for group in terrain_groups()])
end

function _dem_role_audit(
    screen::DataFrame, spatial::DataFrame; scheme::String, product::String,
    fold::Int, phase::String, selection_status::String, error::String="",
)
    rows = NamedTuple[]
    for group in terrain_groups()
        screen_row = only(eachrow(filter(:variable_group => ==(group), screen)))
        spatial_match = :variable_group in propertynames(spatial) ?
            filter(:variable_group => ==(group), spatial) : DataFrame()
        selected = Bool(screen_row.selected)
        role = selected ? (nrow(spatial_match) == 1 ? String(spatial_match.role[1]) : "uncertain") :
            "not_selected"
        model_included = selected && role in ("local", "global")
        push!(rows, (;
            scheme, product, fold, phase, variable_group=group,
            predictor_columns=join(string.(terrain_columns(group)), "+"),
            selected_pre_vif=Bool(screen_row.selected_pre_vif), selected,
            model_included,
            selection_qvalue=Float64(screen_row.qvalue), role,
            variability_pvalue=nrow(spatial_match) == 1 ? Float64(spatial_match.pvalue[1]) : NaN,
            variability_qvalue=nrow(spatial_match) == 1 ? Float64(spatial_match.qvalue[1]) : NaN,
            selection_status, error,
        ))
    end
    return DataFrame(rows)
end

function screen_dem_subset(
    terrain::DataFrame, lonlat::Matrix{Float64}, times::Vector{DateTime},
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, dem::DEMExperimentConfig;
    scheme::String, product::String, fold::Int, phase::String, seed::Int,
)
    response, counts = mean_wet_residual(
        y_obs, y_sat; threshold=dem.wet_threshold, min_hours=dem.min_wet_hours,
    )
    valid_count = count(isfinite, response)
    finite_counts = isempty(counts) ? [0] : counts
    qc = DataFrame([(
        scheme=scheme, product=product, fold=fold, phase=phase,
        training_station_count=nrow(terrain), valid_residual_station_count=valid_count,
        wet_threshold=dem.wet_threshold, min_wet_hours=dem.min_wet_hours,
        min_station_wet_hours=minimum(finite_counts),
        median_station_wet_hours=median(finite_counts),
        max_station_wet_hours=maximum(finite_counts),
    )])

    screen = DataFrame()
    vif = DataFrame()
    selected = String[]
    selection_status = "selected"
    selection_error = ""
    try
        screen, vif, selected = terrain_screen(
            terrain, response; product, permutations=dem.screen_permutations,
            q_threshold=dem.q_threshold, vif_threshold=dem.vif_threshold, seed,
        )
        isempty(selected) && (selection_status = "no_dem_selected")
    catch error
        selection_status = "screen_failed"
        selection_error = sprint(showerror, error)
        screen = _failed_screen_table(product, selection_error)
    end

    monthly = monthly_correlation_rows(product, times, terrain, y_obs, y_sat, dem)
    monthly_table = DataFrame(monthly)
    spatial = DataFrame()
    spatial_scan = DataFrame()
    if !isempty(selected)
        try
            spatial, spatial_scan = spatial_variability_test(
                terrain, lonlat, response, selected; product,
                bandwidth_candidates=dem.bandwidth_candidates,
                permutations=dem.spatial_permutations, q_threshold=dem.q_threshold,
                seed=seed + 10_000, ridge=dem.ridge,
            )
        catch error
            selection_status = "spatial_test_failed"
            selection_error = sprint(showerror, error)
        end
    end
    role_table = _dem_role_audit(
        screen, spatial; scheme, product, fold, phase, selection_status,
        error=selection_error,
    )
    role_map = Dict{String,String}()
    for row in eachrow(role_table)
        row.selected && row.role in ("local", "global") || continue
        role_map[String(row.variable_group)] = String(row.role)
    end
    if !isempty(selected) && isempty(role_map) && selection_status == "selected"
        selection_status = "no_role_resolved"
        role_table.selection_status .= selection_status
    end
    _tag_dem_table!(screen, scheme, product, fold, phase)
    _tag_dem_table!(vif, scheme, product, fold, phase)
    _tag_dem_table!(monthly_table, scheme, product, fold, phase)
    _tag_dem_table!(spatial, scheme, product, fold, phase)
    _tag_dem_table!(spatial_scan, scheme, product, fold, phase)
    return (;
        response, counts, qc, screen, vif, monthly=monthly_table, spatial,
        spatial_scan, roles=role_table, role_map, selection_status,
        error=selection_error,
    )
end

function build_dem_fold_context(
    selection, terrain_train::DataFrame, terrain_target::DataFrame,
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    dem::DEMExperimentConfig,
)
    role_map = selection.role_map
    all_local_map = Dict{String,String}(group => "local" for group in keys(role_map))
    mixed = terrain_model_designs(
        terrain_train, terrain_target, train_lonlat, target_lonlat, role_map,
    )
    all_local = terrain_model_designs(
        terrain_train, terrain_target, train_lonlat, target_lonlat, all_local_map,
    )
    valid_aggregate = isfinite.(selection.response)
    usable_candidates = filter(<(count(valid_aggregate)), dem.bandwidth_candidates)
    isempty(usable_candidates) && throw(ArgumentError(
        "no DEM bandwidth candidate is smaller than the valid training-station count",
    ))
    variables = sort(collect(keys(role_map)))
    role_text = join(["$group=$(role_map[group])" for group in variables], ";")
    return (;
        mixed, all_local, aggregate=selection.response, valid_aggregate,
        bandwidth_candidates=usable_candidates, variables,
        dem_variable_count=length(variables), dem_roles=role_text,
        selection_status=selection.selection_status, train_lonlat, dem,
    )
end

function select_dem_parameter!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    method::String, mode::String, scheme::Symbol, product::String, fold::Int,
    dem_context,
)
    mode == "residual" || throw(ArgumentError("DEM models only support residual mode"))
    dem = something(cfg.dem)
    valid = dem_context.valid_aggregate
    response = dem_context.aggregate[valid]
    n = length(response)

    if method in ("gwr", "mixed_gwr")
        designs = method == "gwr" ? dem_context.all_local : dem_context.mixed
        local_design = designs.mixed_local_train[valid, :]
        global_design = method == "gwr" ? zeros(Float64, n, 0) :
            designs.global_train[valid, :]
        train_lonlat = dem_context.train_lonlat[valid, :]
        bandwidth, scan = select_mixed_bandwidth(
            local_design, global_design, response, train_lonlat,
            dem_context.bandwidth_candidates; ridge=dem.ridge,
            tolerance=dem.tolerance, max_iterations=dem.max_iterations,
        )
        for row in eachrow(scan)
            converged = Bool(row.converged)
            push!(scan_rows, _scan_row(;
                scheme, product, fold, mode, method,
                group=method == "gwr" ? "all_dem_local" : "shared_local",
                kernel=BISQUARE, adaptive=true, bw=row.bandwidth, n,
                coverage=1.0, RMSE=row.RMSE, MAE=NaN,
                status=row.status == "success" && converged ? "success" : "failed",
                error=converged ? String(row.error) : "backfitting did not converge",
                selected=row.bandwidth == bandwidth,
            ))
        end
        return (; bw=Float64(bandwidth), dem_model=true)
    elseif method == "mgwr"
        designs = dem_context.mixed
        local_designs = [X[valid, :] for X in designs.multiscale_train]
        global_design = designs.global_train[valid, :]
        train_lonlat = dem_context.train_lonlat[valid, :]
        bandwidths, scan, converged = select_multiscale_bandwidths(
            local_designs, global_design, response, train_lonlat,
            dem_context.bandwidth_candidates; ridge=dem.ridge,
            tolerance=dem.tolerance, max_iterations=dem.max_iterations,
        )
        final_iteration = nrow(scan) == 0 ? 0 : maximum(scan.iteration)
        for row in eachrow(scan)
            group_name = designs.local_group_names[row.group_index]
            selected = converged && row.iteration == final_iteration &&
                row.bandwidth == bandwidths[row.group_index] && row.status == "success"
            push!(scan_rows, _scan_row(;
                scheme, product, fold, mode, method, group=group_name,
                iteration=row.iteration, kernel=BISQUARE, adaptive=true,
                bw=row.bandwidth, n, coverage=1.0, RMSE=row.RMSE, MAE=NaN,
                status=String(row.status), error=String(row.error), selected,
            ))
        end
        converged || throw(ArgumentError("MGWR bandwidth backfitting did not converge"))
        return (; bandwidths, dem_model=true)
    end
    throw(ArgumentError("unsupported DEM residual method: $method"))
end

function _empty_dem_store()
    return Dict(key => DataFrame[] for key in
        (:qc, :screen, :vif, :monthly, :spatial, :spatial_scan, :roles))
end

function _store_dem_selection!(store::Dict, selection)
    for key in keys(store)
        push!(store[key], getproperty(selection, key))
    end
    return store
end

function _combine_dem_tables(store::Dict, key::Symbol)
    tables = store[key]
    isempty(tables) && return DataFrame()
    return vcat(tables...; cols=:union)
end

function _dem_role_stability(roles::DataFrame)
    nrow(roles) == 0 && return DataFrame()
    cv_roles = filter(:phase => ==("cv"), roles)
    rows = NamedTuple[]
    for group in groupby(cv_roles, [:scheme, :product, :variable_group])
        resolved = group.role .∈ Ref(["local", "global"])
        selected_count = count(resolved)
        local_count = count(==("local"), group.role)
        global_count = count(==("global"), group.role)
        uncertain_count = count(==("uncertain"), group.role)
        total = nrow(group)
        push!(rows, (;
            scheme=String(group.scheme[1]), product=String(group.product[1]),
            variable_group=String(group.variable_group[1]), fold_count=total,
            selected_count, local_count, global_count, uncertain_count,
            not_selected_count=count(==("not_selected"), group.role),
            selection_frequency=total > 0 ? selected_count / total : NaN,
            local_frequency_selected=selected_count > 0 ? local_count / selected_count : NaN,
            global_frequency_selected=selected_count > 0 ? global_count / selected_count : NaN,
        ))
    end
    return DataFrame(rows)
end

function _dem_bandwidth_table(scans::DataFrame)
    nrow(scans) == 0 && return DataFrame()
    rows = NamedTuple[]
    for row in eachrow(scans)
        row.mode == "residual" || continue
        row.method in ("gwr", "mixed_gwr", "mgwr") || continue
        row.selected || continue
        push!(rows, (;
            scheme=String(row.scheme), product=String(row.product), fold=Int(row.fold),
            method=row.method == "gwr" ? "residual_gwr" : String(row.method),
            variable_group=String(row.group), bandwidth_neighbors=Int(round(row.bw)),
            iteration=Int(row.iteration), selection_RMSE=Float64(row.RMSE),
        ))
    end
    return DataFrame(rows)
end

function _write_dem_outputs(
    outdir::String, store::Dict, scans::DataFrame, status::DataFrame,
)
    qc = _combine_dem_tables(store, :qc)
    screen = _combine_dem_tables(store, :screen)
    vif = _combine_dem_tables(store, :vif)
    monthly = _combine_dem_tables(store, :monthly)
    spatial = _combine_dem_tables(store, :spatial)
    spatial_scan = _combine_dem_tables(store, :spatial_scan)
    roles = _combine_dem_tables(store, :roles)
    correlation = filter(:variable_group => !=("aspect"), screen)
    aspect = filter(:variable_group => ==("aspect"), screen)
    final_spec = filter(:phase => ==("full_data"), roles)
    CSV.write(joinpath(outdir, "dem_fold_quality_control.csv"), qc)
    CSV.write(joinpath(outdir, "dem_fold_correlation.csv"), correlation)
    CSV.write(joinpath(outdir, "dem_fold_aspect_joint_test.csv"), aspect)
    CSV.write(joinpath(outdir, "dem_fold_vif.csv"), vif)
    CSV.write(joinpath(outdir, "dem_fold_monthly_correlation.csv"), monthly)
    CSV.write(joinpath(outdir, "dem_fold_spatial_variability.csv"), spatial)
    CSV.write(joinpath(outdir, "dem_spatial_bandwidth_scan.csv"), spatial_scan)
    CSV.write(joinpath(outdir, "dem_fold_roles.csv"), roles)
    CSV.write(joinpath(outdir, "dem_role_stability.csv"), _dem_role_stability(roles))
    CSV.write(joinpath(outdir, "dem_final_full_data_spec.csv"), final_spec)
    CSV.write(joinpath(outdir, "dem_bandwidths.csv"), _dem_bandwidth_table(scans))
    CSV.write(joinpath(outdir, "model_status.csv"), status)
    return nothing
end
