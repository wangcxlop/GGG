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
    schemes = selection_schemes(n, folds, cfg.k)
    for product in products
        Ysat = Float64.(satellite[product])
        for (scheme, fold, train_indices) in schemes
            run_seed = cfg.seed + 1000 * findfirst(==(product), products) + 10 * fold
            panel = prepare_dynamic_panel(Yobs, Ysat, era5, train_indices, times, cfg)
            qc = copy(panel.qc); annotate_selection!(qc, product, scheme, fold); append_selection!(qc_all, qc)
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
                annotate_selection!(table, product, scheme, fold)
            end
            append_selection!(association_all, screen.association)
            append_selection!(monthly_all, screen.monthly)
            append_selection!(vif_all, screen.vif)
            spatial = panel_spatial_variability_test(panel, screen.selected, lonlat[train_indices, :], cfg;
                rng=MersenneTwister(run_seed + 1))
            annotate_selection!(spatial.bandwidth_scan, product, scheme, fold)
            annotate_selection!(spatial.variability, product, scheme, fold)
            append_selection!(bandwidth_all, spatial.bandwidth_scan)
            append_selection!(variability_all, spatial.variability)
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
