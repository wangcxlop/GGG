"""
Products derived by fusing the ones that have files, keyed by product name.

Each entry carries the variant and the source matrices, which are the real products' own `Y_sat`
in `products` order - so the fusion's design is the same one everywhere, and
`fused_anchor_selection.csv` can name its coefficients.

The derived product shares the source products' `times`, `ids` and `Y_obs`:
`load_global_common_product_data` has already put every product on one common station set and one
common time grid, so there is nothing to align. Its `Y_sat` starts as NaN and is filled per fold,
because the anchor depends on which stations were held out - which is exactly what keeps a held-out
gauge out of the anchor its own prediction is built on.
"""
function _build_fused_products!(
    cfg::InterpolationBenchmarkConfig, products::Vector{String}, product_data::Dict{String,Any},
    lonlat::AbstractMatrix,
)
    isempty(cfg.fused_anchor_variants) && return String[], Dict{String,Any}()
    source_names = copy(products)
    sources = [Matrix{Float64}(product_data[name].Y_sat) for name in source_names]
    template = product_data[first(source_names)]
    derived = String[]
    specifications = Dict{String,Any}()
    for variant in sort(cfg.fused_anchor_variants)
        name = _fused_product_name(variant)
        haskey(product_data, name) &&
            throw(ArgumentError("derived product $name collides with a configured product"))
        product_data[name] = (;
            template.times, template.ids, template.Y_obs,
            Y_sat=fill(NaN, size(template.Y_obs)),
        )
        # A grouped variant's design and bands depend only on the source products, so they are
        # built once here; only the coefficients are refitted per fold.
        design = if variant in GROUPED_FUSION_VARIANTS
            distance = haversine_distance_matrix(lonlat, lonlat)
            neighbours = [filter(!=(station), sortperm(distance[station, :]))
                          for station in axes(distance, 1)]
            groups = blend_band_matrix(
                :agreement_envelope, first(sources), sources, cfg.mger.rain_threshold,
            ) .+ 1
            (; features=fusion_features(sources, template.times, neighbours), groups,
                n_group=blend_band_count(:agreement_envelope) + 1)
        else
            nothing
        end
        specifications[name] = (; variant, source_names, sources, design)
        push!(derived, name)
    end
    return derived, specifications
end

function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    source_products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))
    # Derived products are appended after validation, and the joint inputs below are loaded for the
    # source products only: `load_joint_benchmark_inputs` reads a fixed specification keyed by
    # product, which a product with no file cannot appear in, and which the config validation has
    # already ruled out alongside a fused anchor.
    _, fused_specifications = _build_fused_products!(cfg, source_products, product_data, lonlat)
    products = benchmark_products(cfg, source_products)
    terrain = _dem_enabled(cfg) ? load_aligned_terrain(something(cfg.terrain_path), ids) : nothing
    dem_store = _empty_dem_store()
    common_times = product_data[first(products)].times
    joint_inputs = _joint_enabled(cfg) ? load_joint_benchmark_inputs(
        something(cfg.joint_covariates), source_products, ids, common_times,
    ) : nothing
    nested_joint = _joint_enabled(cfg) && cfg.joint_selection !== nothing
    # Band matrices for the non-constant blend axes. They are read off the source products alone,
    # so they depend on neither the fold nor the scheme nor which product is being corrected, and
    # are built once here rather than per fold. The constant axis keeps its own code path in
    # `_run_fold_blend!` so that turning another axis on cannot perturb what it reports.
    blend_bands = Dict{Symbol,Matrix{Int}}()
    if cfg.satellite_wet_blend
        band_sources = [Matrix{Float64}(product_data[name].Y_sat) for name in source_products]
        for axis in cfg.blend_axes
            axis === :constant && continue
            blend_bands[axis] = blend_band_matrix(
                axis, first(band_sources), band_sources, cfg.mger.rain_threshold,
            )
        end
    end
    # The stack's agreement bands come from the same source products, so they too are built once.
    stack_agreement = if cfg.stacked_blend
        stack_sources = [Matrix{Float64}(product_data[name].Y_sat) for name in source_products]
        blend_band_matrix(
            :agreement_envelope, first(stack_sources), stack_sources, cfg.mger.rain_threshold,
        )
    else
        nothing
    end
    joint_store = _empty_joint_store()
    joint_scaling_tables = DataFrame[]
    joint_qc_tables = DataFrame[]
    joint_inputs === nothing || _write_joint_provenance(cfg, joint_inputs, nested_joint)

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]
    auto_selection_rows = NamedTuple[]
    blend_selection_rows = NamedTuple[]
    stack_selection_rows = NamedTuple[]
    fused_anchor_rows = NamedTuple[]
    fused_grouped_rows = NamedTuple[]
    hurdle_rows = NamedTuple[]
    # Cells where the fold was too small for an inner selection split and fell back to
    # leave-one-out. Reported in `benchmark_scope.csv` so the fallback is never silent.
    selection_fallback_cells = String[]

    _dem_enabled(cfg) && _screen_dem_full_data!(dem_store, cfg, terrain, lonlat, products, product_data)

    seeds = benchmark_seeds(cfg)
    # Repeated cross-validation: one independent fold partition per repeat.
    #
    # `repeat_seed` and the partition index are deliberately separate. Under the default
    # `:hilbert` initialisation the partition comes from `rotation = repeat_index - 1` and no RNG
    # is involved at all; `repeat_seed` still seeds the genuinely stochastic sub-processes below
    # (paired bootstrap, DEM permutation tests, joint variable selection) and the `:random`
    # scheme, which is random by definition.
    for (repeat_index, repeat_seed) in enumerate(seeds)
        repeat_root = _repeat_dir(cfg.mger.outdir, repeat_index, length(seeds))
        for scheme_symbol in cfg.cv_schemes
            scheme = string(scheme_symbol)
            scheme_dir = joinpath(repeat_root, scheme)
            mkpath(scheme_dir)
            folds = benchmark_folds(
                scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
                center_init=cfg.fold_center_init, rotation=repeat_index - 1,
            )
            _write_split(joinpath(scheme_dir, "split_common.csv"), ids, folds, scheme_symbol)
            id_map = Dict(id => index for (index, id) in enumerate(ids))

            for product in products
                data = product_data[product]
                y_obs = data.Y_obs
                y_sat = data.Y_sat
                fusion = get(fused_specifications, product, nothing)
                predictions = Dict(
                    method => fill(NaN, size(y_obs)) for method in benchmark_methods(cfg)
                )
                # A derived product has no anchor until a fold says which stations may be used to
                # fit one, so `raw` is filled per fold below instead of in one assignment here.
                fusion === nothing && (predictions["raw"] .= y_sat)
                nearest_train_distance = fill(NaN, length(ids))

                for fold in 1:cfg.k
                    fold_sat = y_sat
                    if fusion !== nothing && fusion.design !== nothing
                        _, _, fusion_train_idx, fusion_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        design = fusion.design
                        coefficients, used, fell_back = fit_grouped_fusion(
                            y_obs, fusion.sources, design.features, design.groups,
                            design.n_group, fusion_train_idx,
                        )
                        coefficients === nothing && error(
                            "$product fold $fold: the pooled fusion fit is singular",
                        )
                        fold_sat = apply_grouped_fusion(
                            fusion.sources, design.features, coefficients, design.groups,
                        )
                        predictions["raw"][fusion_val_idx, :] = fold_sat[fusion_val_idx, :]
                        lower = lowercase.(fusion.source_names)
                        terms = vcat("intercept", lower, lower .* "_lag1",
                            lower .* "_lead1", lower .* "_nbr8")
                        for group in 1:design.n_group, (row, term) in enumerate(terms)
                            push!(fused_grouped_rows, (;
                                scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                                variant=String(fusion.variant), group=group - 1,
                                n_cell=used[group], fell_back=fell_back[group], term,
                                coefficient=coefficients[row, group],
                            ))
                        end
                    elseif fusion !== nothing
                        _, _, fusion_train_idx, fusion_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        coefficients, used, fell_back = fusion_coefficients(
                            fusion.variant, y_obs, fusion.sources, fusion_train_idx,
                        )
                        # Fitted on the training stations, applied everywhere: the fold's own
                        # training rows need the same anchor its held-out rows get, or the
                        # residual target would mean two different things inside one fit.
                        fold_sat = apply_satellite_fusion(fusion.sources, coefficients)
                        predictions["raw"][fusion_val_idx, :] = fold_sat[fusion_val_idx, :]
                        push!(fused_anchor_rows, merge(
                            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                                variant=String(fusion.variant),
                                n_train_station=length(fusion_train_idx), n_train_cell=used,
                                fell_back, intercept=coefficients[1]),
                            NamedTuple{Tuple(Symbol("beta_", lowercase(name))
                                             for name in fusion.source_names)}(
                                Tuple(coefficients[2:end])),
                        ))
                    end
                    _run_benchmark_fold!(
                        cfg, fold, folds, id_map, ids, products, product,
                        data, lonlat, y_obs, fold_sat, terrain, joint_inputs, nested_joint,
                        scheme, scheme_symbol, repeat_index, repeat_seed,
                        predictions, nearest_train_distance, blend_bands,
                        dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
                        all_metric_rows, all_scan_rows, run_status_rows, auto_selection_rows,
                        blend_selection_rows, hurdle_rows, selection_fallback_cells,
                        stack_agreement, stack_selection_rows,
                    )
                end

                common_mask = _common_method_mask(y_obs, predictions)
                if !any(common_mask)
                    failures = ["fold=$(row.fold) method=$(row.method): $(row.error)" for
                        row in run_status_rows if row.scheme == scheme && row.product == product &&
                        row.status != "success"]
                    error("[$scheme/$product] no common valid OOF samples across all methods; " *
                        join(failures, " | "))
                end
                product_dir = joinpath(scheme_dir, lowercase(product))
                mkpath(product_dir)
                for method in benchmark_methods(cfg)
                    # The per-station OOF tables are large; only the first repeat writes them.
                    repeat_index == 1 && write_wide(
                        joinpath(product_dir, "oof_$(method).csv"), data.times, ids, predictions[method],
                    )
                    append_stratified_metrics!(
                        all_metric_rows, scheme, product, method, data.times, y_obs,
                        predictions[method], common_mask, nearest_train_distance, cfg.event_thresholds;
                        repeat=repeat_index, seed=repeat_seed,
                    )
                end
                if repeat_index == 1
                    mask_df = DataFrame(time=Dates.format.(data.times, dateformat"yyyy-mm-ddTHH:MM:SS"))
                    for (station_index, station_id) in enumerate(ids)
                        mask_df[!, Symbol(station_id)] = common_mask[station_index, :]
                    end
                    CSV.write(joinpath(product_dir, "common_evaluation_mask.csv"), mask_df)
                end
                if scheme_symbol == :balanced_spatial && cfg.bootstrap_reps > 0
                    append!(all_bootstrap_rows, paired_bootstrap_rows(
                        cfg, scheme, product, data.times, y_obs, predictions, common_mask;
                        repeat=repeat_index, seed=repeat_seed,
                    ))
                end
            end
        end
    end # repeat loop

    # One row per fold x band x term; `group` is the agreement-envelope band (0 = no product
    # wet). `fell_back` marks a band too sparse or singular to fit alone, which took the pooled fit.
    isempty(fused_grouped_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "fused_anchor_grouped_coefficients.csv"),
            DataFrame(fused_grouped_rows))
    # One row per fold x stack band: the band's inner-split cell count, whether it took the pooled
    # fit, and its weight on every `STACK_FEATURES` member.
    isempty(stack_selection_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "stack_selection.csv"), DataFrame(stack_selection_rows))

    return _write_benchmark_outputs(
        cfg, products, seeds, nested_joint, joint_inputs,
        dem_store, joint_store, joint_scaling_tables, joint_qc_tables,
        all_metric_rows, all_scan_rows, all_bootstrap_rows, run_status_rows,
        auto_selection_rows, blend_selection_rows, fused_anchor_rows, hurdle_rows,
        selection_fallback_cells,
    )
end
