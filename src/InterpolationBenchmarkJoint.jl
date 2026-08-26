function load_joint_benchmark_inputs(
    joint::JointCovariateBenchmarkConfig, products::Vector{String}, ids::Vector{String},
    times::Vector{DateTime},
)
    specification = joint.spec_path === nothing ? nothing :
        load_joint_covariate_spec(joint.spec_path, products)
    terrain = load_aligned_terrain(joint.terrain_path, ids)
    required_years = sort(unique(year(time - Hour(joint.feature_time_offset_hours)) for time in times))
    missing_years = filter(year_value -> !haskey(joint.era5_paths, year_value), required_years)
    isempty(missing_years) || throw(ArgumentError(
        "joint ERA5 paths are missing years: $(join(missing_years, ", "))",
    ))
    era5 = ERA5VariableSelection.load_era5_panel(
        joint.era5_paths, ids, times; offset_hours=joint.feature_time_offset_hours,
    )
    Bool(era5.qc.complete[1]) || throw(ArgumentError("joint ERA5 input failed quality control"))
    # Nested mode can't know ahead of time whether any training fold will select NDVI, so load
    # it whenever it's configured; fixed mode keeps the narrower "only if the spec uses it" check.
    uses_ndvi = specification === nothing ? joint.ndvi_path !== nothing :
        any("ndvi" in keys(role_map) for role_map in values(specification.role_maps))
    ndvi = if uses_ndvi
        joint.ndvi_path === nothing && throw(ArgumentError("NDVI is selected but ndvi_path is absent"))
        loaded = NDVIVariableSelection.load_ndvi_covariates(something(joint.ndvi_path), ids)
        aligned = NDVIVariableSelection.align_ndvi_asof(
            loaded.table, ids, times; max_age_days=joint.max_ndvi_age_days,
        )
        (; loaded, aligned)
    else
        nothing
    end
    return (; specification, terrain, era5, ndvi,
        spec_sha256=joint.spec_path === nothing ? nothing : joint_spec_sha256(joint.spec_path))
end

"""
Score one joint-covariate candidate over the training stations.

`selection_contexts`, when given, is one `(context, target_positions)` pair per inner selection
group: a `JointFoldContext` built with that group's stations as its *target* and the rest of the
fold's training stations as its *train*. Scoring then walks those with `leave_one_out=false`,
which predicts at `context.target_lonlat`, and stitches the results back into one out-of-fold
matrix over all training stations — the same shape and meaning the leave-one-out path produced.

`build_joint_fold_context` fits every scaled quantity on its `train_indices` alone, so each inner
context is leak-free by construction. The `roles` map is the one selected on the whole outer
training fold and so did see the inner-target stations; that cannot bias the *bandwidth*
comparison, since the covariate set is held fixed across every candidate, and re-running role
selection per inner group would multiply the permutation cost for no effect on the ranking.
"""
function _joint_candidate_metrics(
    context, residuals::Matrix{Float64}, y_obs::Matrix{Float64},
    y_sat::Matrix{Float64}, method::String, bandwidths::Vector{Int},
    time_indices::Vector{Int}, time_weights::Union{Nothing,Vector{Float64}}=nothing,
    selection_contexts=nothing, shrink_candidates::Vector{Float64}=[1.0],
)
    isempty(shrink_candidates) &&
        throw(ArgumentError("shrink_candidates must not be empty"))
    residual_prediction, converged = if selection_contexts === nothing
        dynamic_covariate_predict(
            context, residuals, method, bandwidths; time_indices, leave_one_out=true,
        )
    else
        out_of_fold = fill(NaN, size(residuals, 1), length(time_indices))
        any_converged = falses(length(time_indices))
        for entry in selection_contexts
            # `train_residuals` is the inner-training row slice, materialised once per group by
            # `select_joint_parameter!` rather than re-sliced for every candidate — at full size
            # that copy is ~14 MB and the scan walks dozens of candidates.
            group_prediction, group_converged = dynamic_covariate_predict(
                entry.context, entry.train_residuals, method, bandwidths;
                time_indices, leave_one_out=false,
            )
            out_of_fold[entry.target_positions, :] = group_prediction
            any_converged .|= group_converged
        end
        (out_of_fold, any_converged)
    end
    # The fit is the expensive part and does not depend on the shrinkage scale, so the ladder is
    # rescored against the same `residual_prediction`. Scoring the clipped correction (rather than
    # solving `E[g*r]/E[g^2]` in closed form) keeps the criterion identical to the reported metric,
    # which the `max(., 0)` clip would otherwise break.
    best = nothing
    for shrink in shrink_candidates
        corrected = max.(y_sat[:, time_indices] .+ shrink .* residual_prediction, 0.0)
        metrics = _candidate_metrics(
            y_obs[:, time_indices], y_sat[:, time_indices], corrected; time_weights,
        )
        if best === nothing || (metrics.RMSE, metrics.MAE) < (best.RMSE, best.MAE)
            best = merge(metrics, (; shrink))
        end
    end
    return merge(best, (;
        converged_times=count(converged), total_times=length(converged),
    ))
end

function select_joint_parameter!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    method::String, mode::String, scheme::Symbol, product::String, fold::Int,
    joint_context, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    joint_selection_contexts=nothing,
)
    mode == "residual" || throw(ArgumentError("joint covariates only support residual mode"))
    joint_method = method == "gwr" ? "residual_gwr" : method
    residuals = y_obs .- y_sat
    # Materialise each inner group's training rows once, not once per candidate.
    selection_entries = joint_selection_contexts === nothing ? nothing : [(
        context=entry.context, target_positions=entry.target_positions,
        train_residuals=residuals[
            setdiff(1:size(residuals, 1), entry.target_positions), :],
    ) for entry in joint_selection_contexts]
    times, time_weights = _tuning_time_sample(
        y_obs, cfg.tuning_max_times, cfg.tuning_time_weighting,
    )
    candidates = joint_context.bandwidth_candidates
    isempty(candidates) && throw(ArgumentError("no joint bandwidth candidate is usable"))
    group_names = joint_group_names(joint_context, joint_method)
    # Shrinkage rides along with every candidate: the model fit does not depend on it, so the
    # ladder is rescored inside `_joint_candidate_metrics` for free and the winner is reported
    # in the scan row's `shrink` column alongside its bandwidth.
    shrink_candidates = cfg.residual_shrinkage_candidates
    if joint_method in ("residual_gwr", "mixed_gwr")
        first_row = length(scan_rows) + 1
        for bandwidth in candidates
            try
                metrics = _joint_candidate_metrics(
                    joint_context, residuals, y_obs, y_sat, joint_method,
                    [bandwidth], times, time_weights, selection_entries, shrink_candidates,
                )
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method,
                    group=only(group_names), kernel=BISQUARE, adaptive=true,
                    bw=bandwidth, shrink=metrics.shrink, n=metrics.n, coverage=metrics.coverage,
                    RMSE=metrics.RMSE, MAE=metrics.MAE,
                    status=metrics.coverage >= cfg.min_tuning_coverage ? "success" : "failed",
                    error=metrics.coverage >= cfg.min_tuning_coverage ? "" :
                        "LOOCV coverage below minimum",
                ))
            catch error
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method,
                    group=only(group_names), kernel=BISQUARE, adaptive=true,
                    bw=bandwidth, status="failed", error=sprint(showerror, error),
                ))
            end
        end
        selected = _select_candidate!(scan_rows, first_row, cfg.min_tuning_coverage)
        return (; bandwidths=[Int(round(selected.bw))], shrink=selected.shrink,
            joint_model=true, joint_method)
    end

    bandwidths = fill(last(candidates), length(group_names))
    converged = false
    final_iteration = 0
    for iteration in 1:cfg.mgwr_max_tuning_iterations
        previous = copy(bandwidths)
        final_iteration = iteration
        for group_index in eachindex(group_names)
            candidate_rows = Int[]
            for bandwidth in candidates
                trial = copy(bandwidths); trial[group_index] = bandwidth
                try
                    metrics = _joint_candidate_metrics(
                        joint_context, residuals, y_obs, y_sat, joint_method,
                        trial, times, time_weights, selection_entries, shrink_candidates,
                    )
                    push!(scan_rows, _scan_row(;
                        scheme, product, fold, mode, method, iteration,
                        group=group_names[group_index], kernel=BISQUARE, adaptive=true,
                        bw=bandwidth, shrink=metrics.shrink,
                        n=metrics.n, coverage=metrics.coverage,
                        RMSE=metrics.RMSE, MAE=metrics.MAE,
                        status=metrics.coverage >= cfg.min_tuning_coverage ? "success" : "failed",
                        error=metrics.coverage >= cfg.min_tuning_coverage ? "" :
                            "LOOCV coverage below minimum",
                    ))
                catch error
                    push!(scan_rows, _scan_row(;
                        scheme, product, fold, mode, method, iteration,
                        group=group_names[group_index], kernel=BISQUARE, adaptive=true,
                        bw=bandwidth, status="failed", error=sprint(showerror, error),
                    ))
                end
                push!(candidate_rows, length(scan_rows))
            end
            valid = filter(candidate_rows) do index
                row = scan_rows[index]
                row.status == "success" && row.coverage >= cfg.min_tuning_coverage &&
                    isfinite(row.RMSE)
            end
            isempty(valid) && throw(ArgumentError(
                "all MGWR candidates failed for $(group_names[group_index])",
            ))
            selected_index = sort(valid; by=index ->
                (scan_rows[index].RMSE, scan_rows[index].MAE))[1]
            bandwidths[group_index] = Int(round(scan_rows[selected_index].bw))
        end
        if bandwidths == previous
            converged = true
            break
        end
    end
    converged || throw(ArgumentError(
        "joint MGWR bandwidth search did not converge in $(cfg.mgwr_max_tuning_iterations) iterations",
    ))
    selected_shrink = NaN
    for group_index in eachindex(group_names)
        matching = filter(eachindex(scan_rows)) do index
            row = scan_rows[index]
            row.scheme == string(scheme) && row.product == product && row.fold == fold &&
                row.method == method && row.iteration == final_iteration &&
                row.group == group_names[group_index] &&
                Int(round(row.bw)) == bandwidths[group_index] && row.status == "success"
        end
        isempty(matching) && error("selected joint MGWR bandwidth row is absent")
        selected_index = last(matching)
        scan_rows[selected_index] = merge(scan_rows[selected_index], (; selected=true))
        # The search only stops once a whole sweep leaves `bandwidths` unchanged, so every group's
        # winning row in the final iteration was scored against the same bandwidth vector and
        # therefore the same residual prediction. Any of their shrink values is the right one.
        selected_shrink = scan_rows[selected_index].shrink
    end
    return (; bandwidths, shrink=selected_shrink, joint_model=true, joint_method)
end

function _tag_joint_table!(
    table::DataFrame, scheme::String, product::String, fold::Int, repeat::Int, seed::Int,
)
    table[!, :scheme] = fill(scheme, nrow(table))
    :product in propertynames(table) || (table[!, :product] = fill(product, nrow(table)))
    table[!, :fold] = fill(fold, nrow(table))
    table[!, :repeat] = fill(repeat, nrow(table))
    table[!, :seed] = fill(seed, nrow(table))
    select!(table, :scheme, :product, :fold, :repeat, :seed,
        Not([:scheme, :product, :fold, :repeat, :seed]))
    return table
end

"""
Nested per-fold joint-covariate selection: screens, VIF-prunes, and role-tests covariates using
only `terrain[train_idx,:]`-equivalent training data (`select_joint_covariates` receives
`train_idx` and slices internally), then tags the diagnostic tables for storage — the same
relationship `screen_dem_subset` has to `terrain_screen`/`spatial_variability_test`.
"""
function _screen_joint_subset(
    product::String, y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned, train_idx::Vector{Int}, ids::Vector{String},
    times::Vector{DateTime}, lonlat::Matrix{Float64}, selection_cfg::JointSelectionConfig;
    scheme::String, fold::Int, repeat::Int, seed::Int, run_seed::Int, dem_seed::Int,
)
    # `_candidate_rows` unconditionally screens NDVI (unlike DEM/ERA5, it has no "absent" branch),
    # so when no NDVI data is configured we hand it an empty-period alignment instead of `nothing`
    # — `prepare_ndvi_panel` degrades that to `status="insufficient_stations_or_variation"` and
    # NDVI is simply never selected, rather than erroring.
    effective_ndvi = ndvi_aligned !== nothing ? ndvi_aligned : (;
        values=fill(NaN, size(y_obs)), source_period=zeros(Int, size(y_obs)), periods=Date[],
    )
    selection = select_joint_covariates(
        product, y_obs, y_sat, terrain, era5, effective_ndvi, train_idx, ids, times, lonlat,
        selection_cfg, run_seed, dem_seed,
    )
    role_rows = DataFrame([(
        variable_group=group,
        independent_selected=group in selection.active,
        joint_vif_retained=group in selection.retained,
        role=get(selection.role_map, group, "not_selected"),
        final_included=haskey(selection.role_map, group),
    ) for group in JointCovariateModels.JOINT_GROUP_ORDER])
    candidates = copy(selection.candidates)
    vif = copy(selection.vif)
    spatial_result = copy(selection.spatial_result)
    _tag_joint_table!(candidates, scheme, product, fold, repeat, seed)
    _tag_joint_table!(vif, scheme, product, fold, repeat, seed)
    _tag_joint_table!(spatial_result, scheme, product, fold, repeat, seed)
    _tag_joint_table!(role_rows, scheme, product, fold, repeat, seed)
    return (; role_map=selection.role_map, candidates, vif, spatial_result, roles=role_rows)
end

function _empty_joint_store()
    return Dict(key => DataFrame[] for key in (:candidates, :vif, :spatial, :roles))
end

function _store_joint_selection!(store::Dict, selection)
    push!(store[:candidates], selection.candidates)
    push!(store[:vif], selection.vif)
    push!(store[:spatial], selection.spatial_result)
    push!(store[:roles], selection.roles)
    return store
end

function _joint_role_stability(roles::DataFrame)
    nrow(roles) == 0 && return DataFrame()
    rows = NamedTuple[]
    for group in groupby(roles, [:scheme, :product, :variable_group])
        total = nrow(group)
        local_count = count(==("local"), group.role)
        global_count = count(==("global"), group.role)
        included_count = count(group.final_included)
        push!(rows, (;
            scheme=String(group.scheme[1]), product=String(group.product[1]),
            variable_group=String(group.variable_group[1]), fold_count=total,
            included_count, local_count, global_count,
            selection_frequency=total > 0 ? included_count / total : NaN,
            local_frequency_selected=included_count > 0 ? local_count / included_count : NaN,
            global_frequency_selected=included_count > 0 ? global_count / included_count : NaN,
        ))
    end
    return DataFrame(rows)
end

function _write_joint_selection_outputs(outdir::String, store::Dict)
    combine(key) = isempty(store[key]) ? DataFrame() : vcat(store[key]...; cols=:union)
    roles = combine(:roles)
    CSV.write(joinpath(outdir, "joint_fold_candidates.csv"), combine(:candidates))
    CSV.write(joinpath(outdir, "joint_fold_vif.csv"), combine(:vif))
    CSV.write(joinpath(outdir, "joint_fold_spatial_variability.csv"), combine(:spatial))
    CSV.write(joinpath(outdir, "joint_fold_roles.csv"), roles)
    CSV.write(joinpath(outdir, "joint_role_stability.csv"), _joint_role_stability(roles))
    return nothing
end
