function _metric_row(
    product::String, fold::Int, method::String, quantity::String, stratum::String,
    truth::Matrix{Float64}, prediction::Matrix{Float64}, mask::BitMatrix,
    eligible_n::Int,
)
    n = count(mask)
    n == 0 && return (; product, fold, method, quantity, stratum, n=0,
        coverage=0.0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN)
    actual = truth[mask]
    estimated = prediction[mask]
    error = estimated - actual
    return (;
        product, fold, method, quantity, stratum, n,
        coverage=eligible_n > 0 ? n / eligible_n : 0.0, RMSE=sqrt(mean(abs2, error)),
        MAE=mean(abs, error), Bias=mean(error),
        r=n > 1 && std(actual) > 0 && std(estimated) > 0 ? cor(actual, estimated) : NaN,
    )
end

function _append_metrics!(
    rows::Vector{NamedTuple}, product::String, fold::Int, method::String,
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, residual_prediction::Matrix{Float64},
    threshold::Float64,
)
    residual_truth = Yobs - Ysat
    corrected = max.(Ysat + residual_prediction, 0.0)
    eligible = isfinite.(Yobs) .& isfinite.(Ysat)
    strata = (
        ("all", eligible),
        ("wet", eligible .& (Yobs .>= threshold)),
        ("no_rain", eligible .& (Yobs .< 0.1)),
        ("light", eligible .& (Yobs .>= 0.1) .& (Yobs .< 2.5)),
        ("moderate", eligible .& (Yobs .>= 2.5) .& (Yobs .< 8.0)),
        ("heavy", eligible .& (Yobs .>= 8.0)),
    )
    for (name, eligible_mask) in strata
        bitmask = BitMatrix(eligible_mask .& isfinite.(residual_prediction))
        eligible_n = count(eligible_mask)
        push!(rows, _metric_row(
            product, fold, method, "residual", name,
            residual_truth, residual_prediction, bitmask, eligible_n,
        ))
        push!(rows, _metric_row(
            product, fold, method, "corrected_precipitation", name,
            Yobs, corrected, bitmask, eligible_n,
        ))
    end
    return rows
end

function _write_wide(path::String, times::Vector{DateTime}, ids::Vector{String}, values::Matrix{Float64})
    table = DataFrame(time=Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS"))
    for (index, id) in enumerate(ids)
        table[!, Symbol(id)] = values[index, :]
    end
    CSV.write(path, table)
end

function _consensus_table(role_tables::Vector{DataFrame})
    rows = NamedTuple[]
    products = [nrow(table) == 0 ? "" : String(table.product[1]) for table in role_tables]
    for group in TERRAIN_GROUPS
        selected_products = String[]
        roles = String[]
        for (product, table) in zip(products, role_tables)
            matches = filter(:variable_group => ==(group), table)
            nrow(matches) == 1 || continue
            push!(selected_products, product)
            push!(roles, String(matches.role[1]))
        end
        consensus_selected = length(selected_products) >= 2
        local_count = count(==("local"), roles)
        global_count = count(==("global"), roles)
        consensus_role = local_count >= 2 ? "local" : global_count >= 2 ? "global" :
            consensus_selected ? "product_specific" : "not_selected"
        push!(rows, (;
            variable_group=group, selected_product_count=length(selected_products),
            selected_products=join(selected_products, ","), local_product_count=local_count,
            global_product_count=global_count, consensus_selected, consensus_role,
        ))
    end
    return DataFrame(rows)
end

"""Run the standalone two-step DEM screening and spatial-validation experiment."""
function run_dem_experiment(
    cfg::DEMExperimentConfig, products::Vector{String}, ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64},
    satellite::Dict{String,Matrix{Float64}}, terrain::DataFrame,
    lonlat::Matrix{Float64},
)
    mkpath(cfg.outdir)
    length(ids) == nrow(terrain) == size(lonlat, 1) == size(Yobs, 1) ||
        throw(DimensionMismatch("station dimensions are not aligned"))
    allunique(ids) || throw(ArgumentError("duplicate station IDs"))
    terrain_ids = String.(terrain.station_id)
    terrain_ids == ids || throw(ArgumentError("terrain rows must follow common station order"))
    all(isfinite, Matrix{Float64}(terrain[:, [:elevation_m, :slope_deg, :aspect_sin, :aspect_cos]])) ||
        throw(ArgumentError("non-finite terrain values"))
    2 <= cfg.k <= length(ids) || throw(ArgumentError("invalid spatial fold count"))

    qc = DataFrame(
        key=["station_count", "time_count", "products", "wet_threshold", "min_wet_hours",
            "screen_permutations", "spatial_permutations", "selection_rule"],
        value=[string(length(ids)), string(length(times)), join(products, ","),
            string(cfg.wet_threshold), string(cfg.min_wet_hours),
            string(cfg.screen_permutations), string(cfg.spatial_permutations),
            "correlation_and_VIF_then_Monte_Carlo_spatial_variability"],
    )
    CSV.write(joinpath(cfg.outdir, "data_qc.csv"), qc)

    all_screen = DataFrame[]
    all_vif = DataFrame[]
    all_spatial = DataFrame[]
    all_spatial_scans = DataFrame[]
    all_monthly = NamedTuple[]
    role_tables = DataFrame[]
    responses = Dict{String,Vector{Float64}}()
    counts_by_product = Dict{String,Vector{Int}}()

    for product in products
        Ysat = satellite[product]
        response, counts = _mean_wet_residual(
            Yobs, Ysat; threshold=cfg.wet_threshold, min_hours=cfg.min_wet_hours,
        )
        responses[product] = response
        counts_by_product[product] = counts
        screen, vif, selected = terrain_screen(
            terrain, response; product, permutations=cfg.screen_permutations,
            q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
            seed=cfg.seed + Int(sum(codeunits(product))),
        )
        push!(all_screen, screen)
        push!(all_vif, vif)
        append!(all_monthly, _monthly_rows(product, times, terrain, Yobs, Ysat, cfg))
        if isempty(selected)
            push!(role_tables, DataFrame(
                product=String[], variable_group=String[], bandwidth=Int[],
                variability_statistic=Float64[], pvalue=Float64[], qvalue=Float64[], role=String[],
            ))
            continue
        end
        spatial, scan = spatial_variability_test(
            terrain, lonlat, response, selected; product,
            bandwidth_candidates=cfg.bandwidth_candidates,
            permutations=cfg.spatial_permutations, q_threshold=cfg.q_threshold,
            seed=cfg.seed + 10_000 + Int(sum(codeunits(product))), ridge=cfg.ridge,
        )
        push!(all_spatial, spatial)
        push!(all_spatial_scans, scan)
        push!(role_tables, spatial)
    end

    screen_table = isempty(all_screen) ? DataFrame() : vcat(all_screen...; cols=:union)
    vif_table = isempty(all_vif) ? DataFrame() : vcat(all_vif...; cols=:union)
    spatial_table = isempty(all_spatial) ? DataFrame() : vcat(all_spatial...; cols=:union)
    spatial_scan = isempty(all_spatial_scans) ? DataFrame() : vcat(all_spatial_scans...; cols=:union)
    CSV.write(joinpath(cfg.outdir, "correlation_screen.csv"), screen_table)
    CSV.write(joinpath(cfg.outdir, "vif.csv"), vif_table)
    CSV.write(joinpath(cfg.outdir, "monthly_correlation_direction.csv"), DataFrame(all_monthly))
    CSV.write(joinpath(cfg.outdir, "spatial_variability.csv"), spatial_table)
    CSV.write(joinpath(cfg.outdir, "spatial_bandwidth_scan.csv"), spatial_scan)
    CSV.write(joinpath(cfg.outdir, "cross_product_consensus.csv"), _consensus_table(role_tables))

    count_table = DataFrame(station_id=ids)
    for product in products
        count_table[!, Symbol(lowercase(product), "_wet_hours")] = counts_by_product[product]
        count_table[!, Symbol(lowercase(product), "_mean_wet_residual")] = responses[product]
    end
    CSV.write(joinpath(cfg.outdir, "station_residual_summary.csv"), count_table)

    folds = _balanced_spatial_folds(ids, lonlat; k=cfg.k, seed=cfg.seed)
    fold_map = zeros(Int, length(ids))
    for fold in 1:cfg.k
        fold_map[folds[fold]] .= fold
    end
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), DataFrame(station_id=ids, fold=fold_map))

    metric_rows = NamedTuple[]
    tuning_rows = NamedTuple[]
    status_rows = NamedTuple[]
    method_names = ["residual_gwr", "mixed_gwr", "multiscale_gwr"]
    oof = Dict((product, method) => fill(NaN, size(Yobs))
        for product in products for method in method_names)

    fold_role_rows = NamedTuple[]
    for (product_index, product) in enumerate(products)
        Ysat = satellite[product]
        for fold in 1:cfg.k
            val_idx = folds[fold]
            train_idx = setdiff(1:length(ids), val_idx)

            # Terrain-role selection (which variables are used, and whether each is
            # spatially-varying) is recomputed from training-fold stations only, so a
            # fold's held-out validation stations never influence its own model structure.
            _, _, fold_selected_groups = terrain_screen(
                terrain[train_idx, :], responses[product][train_idx]; product,
                permutations=cfg.screen_permutations, q_threshold=cfg.q_threshold,
                vif_threshold=cfg.vif_threshold,
                seed=cfg.seed + Int(sum(codeunits(product))) + 100 * fold,
            )
            role_map = Dict{String,String}()
            if !isempty(fold_selected_groups)
                fold_spatial, _ = spatial_variability_test(
                    terrain[train_idx, :], lonlat[train_idx, :], responses[product][train_idx],
                    fold_selected_groups; product, bandwidth_candidates=cfg.bandwidth_candidates,
                    permutations=cfg.spatial_permutations, q_threshold=cfg.q_threshold,
                    seed=cfg.seed + 10_000 + Int(sum(codeunits(product))) + 100 * fold,
                    ridge=cfg.ridge,
                )
                role_map = Dict(String(row.variable_group) => String(row.role) for row in eachrow(fold_spatial))
                for row in eachrow(fold_spatial)
                    push!(fold_role_rows, (; product, fold, variable_group=row.variable_group,
                        bandwidth=row.bandwidth, variability_statistic=row.variability_statistic,
                        pvalue=row.pvalue, qvalue=row.qvalue, role=row.role))
                end
            end
            can_fit_dem = !isempty(role_map) && all(!=("uncertain"), values(role_map))

            designs = _designs(
                terrain[train_idx, :], terrain[val_idx, :], lonlat[train_idx, :],
                lonlat[val_idx, :], can_fit_dem ? role_map : Dict{String,String}(),
            )
            aggregate = responses[product][train_idx]
            valid_aggregate = isfinite.(aggregate)
            train_count = count(valid_aggregate)
            usable_candidates = filter(<(train_count), cfg.bandwidth_candidates)
            if train_count < 12 || isempty(usable_candidates)
                push!(status_rows, (; product, fold, method="all", status="failed",
                    error="insufficient aggregate training stations"))
                continue
            end
            aggregate_lonlat = lonlat[train_idx[valid_aggregate], :]
            # This experiment predates kernel selection; it keeps the historical adaptive-bisquare
            # behaviour explicitly rather than sweeping kernels, since it is not the benchmark's
            # `mixed_gwr`/`mgwr` tuning path.
            candidates_f = Float64.(usable_candidates)
            baseline_X = designs.mixed_local_train[valid_aggregate, 1:3]
            baseline_bw, baseline_scan = select_mixed_bandwidth(
                baseline_X, zeros(Float64, train_count, 0), aggregate[valid_aggregate],
                aggregate_lonlat, candidates_f, _bisquare_kernel; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(baseline_scan)
                push!(tuning_rows, (; product, fold, method="residual_gwr", group="all_local",
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.bandwidth == baseline_bw,
                    status=row.status, error=row.error))
            end
            baseline_prediction, baseline_converged = mixed_gwr_predict(
                designs.mixed_local_train[:, 1:3], zeros(Float64, length(train_idx), 0),
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.mixed_local_target[:, 1:3], zeros(Float64, length(val_idx), 0),
                lonlat[val_idx, :], baseline_bw, _bisquare_kernel; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            oof[(product, "residual_gwr")][val_idx, :] = baseline_prediction
            _append_metrics!(metric_rows, product, fold, "residual_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], baseline_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="residual_gwr",
                status=all(baseline_converged) ? "success" : "partial",
                error=all(baseline_converged) ? "" : "one or more hourly fits failed"))

            can_fit_dem || continue
            mixed_bw, mixed_scan = select_mixed_bandwidth(
                designs.mixed_local_train[valid_aggregate, :],
                designs.global_train[valid_aggregate, :], aggregate[valid_aggregate],
                aggregate_lonlat, candidates_f, _bisquare_kernel; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(mixed_scan)
                push!(tuning_rows, (; product, fold, method="mixed_gwr", group="all_local",
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.bandwidth == mixed_bw,
                    status=row.status, error=row.error))
            end
            mixed_prediction, mixed_converged = mixed_gwr_predict(
                designs.mixed_local_train, designs.global_train,
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.mixed_local_target, designs.global_target, lonlat[val_idx, :], mixed_bw,
                _bisquare_kernel; ridge=cfg.ridge, tolerance=cfg.tolerance,
                max_iterations=cfg.max_iterations,
            )
            oof[(product, "mixed_gwr")][val_idx, :] = mixed_prediction
            _append_metrics!(metric_rows, product, fold, "mixed_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], mixed_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="mixed_gwr",
                status=all(mixed_converged) ? "success" : "partial",
                error=all(mixed_converged) ? "" : "one or more hourly fits failed"))

            local_valid = [X[valid_aggregate, :] for X in designs.multiscale_train]
            multiscale_bw, multiscale_scan, bandwidth_converged = select_multiscale_bandwidths(
                local_valid, designs.global_train[valid_aggregate, :], aggregate[valid_aggregate],
                aggregate_lonlat, candidates_f, _bisquare_kernel; ridge=cfg.ridge,
                tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            for row in eachrow(multiscale_scan)
                group_name = designs.local_group_names[row.group_index]
                push!(tuning_rows, (; product, fold, method="multiscale_gwr", group=group_name,
                    bandwidth=row.bandwidth, RMSE=row.RMSE, selected=row.selected,
                    status=row.status, error=row.error))
            end
            if !bandwidth_converged
                push!(status_rows, (; product, fold, method="multiscale_gwr", status="failed",
                    error="bandwidth backfitting did not converge"))
                continue
            end
            multiscale_prediction, multiscale_converged = multiscale_gwr_predict(
                designs.multiscale_train, designs.global_train,
                Yobs[train_idx, :] - Ysat[train_idx, :], lonlat[train_idx, :],
                designs.multiscale_target, designs.global_target, lonlat[val_idx, :],
                multiscale_bw, _bisquare_kernel;
                ridge=cfg.ridge, tolerance=cfg.tolerance, max_iterations=cfg.max_iterations,
            )
            oof[(product, "multiscale_gwr")][val_idx, :] = multiscale_prediction
            _append_metrics!(metric_rows, product, fold, "multiscale_gwr",
                Yobs[val_idx, :], Ysat[val_idx, :], multiscale_prediction, cfg.wet_threshold)
            push!(status_rows, (; product, fold, method="multiscale_gwr",
                status=all(multiscale_converged) ? "success" : "partial",
                error=all(multiscale_converged) ? "" : "one or more hourly fits failed"))
        end
    end

    metrics = DataFrame(metric_rows)
    tuning = DataFrame(tuning_rows)
    status = DataFrame(status_rows)
    fold_roles = DataFrame(fold_role_rows)
    CSV.write(joinpath(cfg.outdir, "validation_metrics.csv"), metrics)
    CSV.write(joinpath(cfg.outdir, "bandwidth_scan.csv"), tuning)
    CSV.write(joinpath(cfg.outdir, "run_status.csv"), status)
    CSV.write(joinpath(cfg.outdir, "spatial_variability_by_fold.csv"), fold_roles)
    for product in products
        product_dir = joinpath(cfg.outdir, lowercase(product))
        mkpath(product_dir)
        for method in method_names
            values = oof[(product, method)]
            any(isfinite, values) || continue
            _write_wide(joinpath(product_dir, "oof_residual_$(method).csv"), times, ids, values)
            corrected = max.(satellite[product] + values, 0.0)
            _write_wide(joinpath(product_dir, "oof_corrected_$(method).csv"), times, ids, corrected)
        end
    end
    return (; screen=screen_table, vif=vif_table, spatial=spatial_table,
        consensus=_consensus_table(role_tables), metrics, tuning, status, fold_roles)
end
