module NDVIVariableSelection

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

include(joinpath(@__DIR__, "ERA5VariableSelection.jl"))
using .ERA5VariableSelection

export NDVISelectionConfig, load_ndvi_covariates, align_ndvi_asof
export prepare_ndvi_panel, screen_ndvi, ndvi_spatial_variability_test
export run_ndvi_variable_selection

Base.@kwdef struct NDVISelectionConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours_per_period::Int = 5
    min_periods_per_station::Int = 8
    min_stations_per_period::Int = 12
    max_age_days::Int = 32
    k::Int = 5
    seed::Int = 20260817
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    association_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    ridge::Float64 = 1e-8
end

function _era_config(cfg::NDVISelectionConfig)
    return ERA5SelectionConfig(
        outdir=cfg.outdir,
        wet_threshold=cfg.wet_threshold,
        min_wet_hours=cfg.min_periods_per_station,
        min_stations_per_time=cfg.min_stations_per_period,
        k=cfg.k,
        seed=cfg.seed,
        bandwidth_candidates=cfg.bandwidth_candidates,
        association_permutations=cfg.association_permutations,
        spatial_permutations=cfg.spatial_permutations,
        q_threshold=cfg.q_threshold,
        ridge=cfg.ridge,
    )
end

function _direction(value::Real)
    !isfinite(value) && return "uncertain"
    value > 0 && return "positive"
    value < 0 && return "negative"
    return "zero"
end

"""Read and validate the 237-station MOD13A2 table used by the NDVI experiment."""
function load_ndvi_covariates(path::AbstractString, station_ids::Vector{String})
    isfile(path) || throw(ArgumentError("NDVI table does not exist: $path"))
    table = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    required = [
        :station_id, :lon, :lat, :composite_start, :observation_date,
        :ndvi_qc, :pixel_reliability, :quality_class, :land_water_class,
    ]
    missing_columns = setdiff(required, propertynames(table))
    isempty(missing_columns) || throw(ArgumentError(
        "NDVI table is missing columns: $(join(missing_columns, ", "))",
    ))
    allunique(table[:, [:station_id, :composite_start]]) ||
        error("Duplicate station/composite-period rows in NDVI table")
    index = Set(station_ids)
    selected = filter(:station_id => in(index), table)
    found = Set(selected.station_id)
    missing_stations = filter(id -> !(id in found), station_ids)
    isempty(missing_stations) || error("NDVI table is missing $(length(missing_stations)) stations")
    extra_stations = count(id -> !(id in index), unique(table.station_id))
    valid_quality = [
        !ismissing(value) && quality in ("good", "marginal")
        for (value, quality) in zip(selected.ndvi_qc, selected.quality_class)
    ]
    all((ismissing(value) || -0.2 <= value <= 1.0) for value in selected.ndvi_qc) ||
        error("ndvi_qc contains values outside [-0.2, 1.0]")
    all(valid_quality .== .!ismissing.(selected.ndvi_qc)) ||
        error("ndvi_qc disagrees with pixel reliability quality class")
    all(!ismissing(selected.observation_date[i]) for i in eachindex(valid_quality) if valid_quality[i]) ||
        error("A quality-valid NDVI row is missing observation_date")
    sort!(selected, [:station_id, :composite_start])
    periods = sort(unique(selected.composite_start))
    station_qc = combine(groupby(selected, :station_id),
        nrow => :period_count,
        :ndvi_qc => (x -> count(!ismissing, x)) => :quality_valid_count,
        :land_water_class => (x -> count(!=("land"), x)) => :nonland_period_count,
    )
    station_qc.quality_valid_fraction = station_qc.quality_valid_count ./ station_qc.period_count
    station_qc.nonland_pixel_flag = station_qc.nonland_period_count .> 0
    data_qc = DataFrame(
        station_count=[length(found)], period_count=[length(periods)], row_count=[nrow(selected)],
        expected_rows=[length(found) * length(periods)], duplicate_keys=[0],
        quality_valid_count=[count(valid_quality)], quality_missing_count=[count(!, valid_quality)],
        nonland_row_count=[count(!=("land"), selected.land_water_class)],
        nonland_station_count=[count(station_qc.nonland_pixel_flag)],
        extra_station_count=[extra_stations], all_237_stations=[length(found) == length(station_ids)],
        complete=[nrow(selected) == length(found) * length(periods)],
    )
    data_qc.complete[1] || error("NDVI station-period coverage is incomplete")
    return (; table=selected, periods, data_qc, station_qc)
end

"""As-of align quality-valid NDVI to precipitation hours without future values.

An observation becomes available at 00:00 Beijing time on the day after
`observation_date`. It remains valid until the next quality-valid observation or
for `max_age_days`, whichever comes first. Missing composites are never interpolated.
"""
function align_ndvi_asof(
    table::DataFrame, station_ids::Vector{String}, times::Vector{DateTime};
    max_age_days::Int=32,
)
    max_age_days > 0 || throw(ArgumentError("max_age_days must be positive"))
    periods = sort(unique(table.composite_start))
    period_index = Dict(period => i for (i, period) in enumerate(periods))
    groups = Dict(first(group.station_id) => group for group in groupby(table, :station_id))
    source_period = zeros(Int, length(station_ids), length(times))
    values = fill(NaN, length(station_ids), length(times))
    effective_times = fill(DateTime(1), length(station_ids), length(times))
    matched_hours = zeros(Int, length(station_ids))
    intervals = DataFrame(
        station_id=String[], composite_start=Date[], observation_date=Date[],
        effective_time=DateTime[], valid_until_exclusive=DateTime[],
        ndvi_qc=Float64[], land_water_class=String[], matched_hour_count=Int[],
    )
    for (station, station_id) in enumerate(station_ids)
        haskey(groups, station_id) || continue
        rows = groups[station_id]
        valid = filter(:ndvi_qc => !ismissing, rows)
        nrow(valid) == 0 && continue
        effective = DateTime.(valid.observation_date) .+ Day(1)
        order = sortperm(effective)
        effective = effective[order]
        valid = valid[order, :]
        pointer = 0
        time_order = sortperm(times)
        for time_index in time_order
            time = times[time_index]
            while pointer < length(effective) && effective[pointer + 1] <= time
                pointer += 1
            end
            pointer == 0 && continue
            start_time = effective[pointer]
            time < start_time + Day(max_age_days) || continue
            value = valid.ndvi_qc[pointer]
            ismissing(value) && continue
            source_period[station, time_index] = period_index[valid.composite_start[pointer]]
            values[station, time_index] = Float64(value)
            effective_times[station, time_index] = start_time
            matched_hours[station] += 1
        end
        for row_index in eachindex(effective)
            valid_until = min(
                effective[row_index] + Day(max_age_days),
                row_index < length(effective) ? effective[row_index + 1] :
                    effective[row_index] + Day(max_age_days),
            )
            source = period_index[valid.composite_start[row_index]]
            push!(intervals, (
                station_id, valid.composite_start[row_index], valid.observation_date[row_index],
                effective[row_index], valid_until, Float64(valid.ndvi_qc[row_index]),
                String(valid.land_water_class[row_index]),
                count(==(source), view(source_period, station, :)),
            ))
        end
    end
    qc = DataFrame(
        station_count=[length(station_ids)], time_count=[length(times)],
        station_hours=[length(station_ids) * length(times)],
        matched_station_hours=[count(>(0), source_period)],
        unmatched_station_hours=[count(==(0), source_period)],
        stations_with_match=[count(>(0), matched_hours)], max_age_days=[max_age_days],
        effective_time_rule=["observation_date_plus_1_day_at_00_BJT"],
        interpolation=["none"],
    )
    return (; source_period, values, effective_times, periods, intervals, qc)
end

function _standardize_panel!(x::Matrix{Float64}, y::Matrix{Float64})
    values, _, weights = ERA5VariableSelection._flatten_balanced(x, y)
    isempty(values) && return NaN, NaN
    mean_value = ERA5VariableSelection._weighted_mean(values, weights)
    scale = sqrt(ERA5VariableSelection._weighted_mean((values .- mean_value) .^ 2, weights))
    if !isfinite(scale) || scale <= 0
        return mean_value, NaN
    end
    for index in eachindex(x)
        isfinite(x[index]) && (x[index] = (x[index] - mean_value) / scale)
    end
    return mean_value, scale
end

function prepare_ndvi_panel(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, aligned,
    train_indices::Vector{Int}, station_ids::Vector{String}, cfg::NDVISelectionConfig,
)
    size(Yobs) == size(Ysat) || throw(DimensionMismatch("observation and satellite matrices differ"))
    size(Yobs) == size(aligned.source_period) || throw(DimensionMismatch("NDVI alignment dimensions differ"))
    n, period_count = length(train_indices), length(aligned.periods)
    residual_sum = zeros(n, period_count)
    ndvi_sum = zeros(n, period_count)
    wet_count = zeros(Int, n, period_count)
    for (local_station, global_station) in enumerate(train_indices), time in axes(Yobs, 2)
        period = aligned.source_period[global_station, time]
        period > 0 || continue
        observed, satellite = Yobs[global_station, time], Ysat[global_station, time]
        ndvi = aligned.values[global_station, time]
        isfinite(observed) && isfinite(satellite) && isfinite(ndvi) || continue
        observed >= cfg.wet_threshold || continue
        residual_sum[local_station, period] += observed - satellite
        ndvi_sum[local_station, period] += ndvi
        wet_count[local_station, period] += 1
    end
    yraw = fill(NaN, n, period_count)
    xraw = fill(NaN, n, period_count)
    mask = wet_count .>= cfg.min_wet_hours_per_period
    for i in 1:n, period in 1:period_count
        mask[i, period] || continue
        yraw[i, period] = residual_sum[i, period] / wet_count[i, period]
        xraw[i, period] = ndvi_sum[i, period] / wet_count[i, period]
    end
    previous = -1
    while count(mask) != previous
        previous = count(mask)
        station_ok = vec(sum(mask, dims=2)) .>= cfg.min_periods_per_station
        mask[.!station_ok, :] .= false
        period_ok = vec(sum(mask, dims=1)) .>= cfg.min_stations_per_period
        mask[:, .!period_ok] .= false
    end
    eligible_station = vec(sum(mask, dims=2)) .>= cfg.min_periods_per_station
    y = fill(NaN, n, period_count)
    x = fill(NaN, n, period_count)
    for period in 1:period_count
        stations = findall(mask[:, period])
        length(stations) >= cfg.min_stations_per_period || continue
        y[stations, period] .= yraw[stations, period] .- mean(yraw[stations, period])
        x[stations, period] .= xraw[stations, period] .- mean(xraw[stations, period])
    end
    mean_value, scale = _standardize_panel!(x, y)
    cells = DataFrame(
        station_id=String[], composite_start=Date[], wet_hour_count=Int[],
        ndvi_raw=Float64[], mean_residual=Float64[], used=Bool[],
    )
    for (local_station, global_station) in enumerate(train_indices), period in 1:period_count
        wet_count[local_station, period] > 0 || continue
        count_value = wet_count[local_station, period]
        push!(cells, (
            station_ids[global_station], aligned.periods[period], count_value,
            ndvi_sum[local_station, period] / count_value,
            residual_sum[local_station, period] / count_value,
            mask[local_station, period],
        ))
    end
    qc = DataFrame(
        training_station_count=[n], eligible_station_count=[count(eligible_station)],
        eligible_period_count=[count(vec(sum(mask, dims=1)) .>= cfg.min_stations_per_period)],
        used_station_periods=[count(mask)],
        min_wet_hours_per_period=[cfg.min_wet_hours_per_period],
        min_periods_per_station=[cfg.min_periods_per_station],
        min_stations_per_period=[cfg.min_stations_per_period],
        ndvi_train_mean=[mean_value], ndvi_train_scale=[scale],
        status=[count(eligible_station) >= cfg.min_stations_per_period && isfinite(scale) ?
            "ok" : "insufficient_stations_or_variation"],
    )
    return (; y, x=Dict(:ndvi => x), mask, eligible_station, cells, qc)
end

function screen_ndvi(
    panel, periods::Vector{Date}, cfg::NDVISelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed),
)
    pearson, spearman, pvalue = ERA5VariableSelection._association(
        panel, :ndvi, cfg.association_permutations, rng,
    )
    stable = isfinite(pearson) && isfinite(spearman) && sign(pearson) == sign(spearman)
    qvalue = pvalue
    selected = stable && isfinite(qvalue) && qvalue < cfg.q_threshold
    association = DataFrame(
        variable=["ndvi"], pearson=[pearson], spearman=[spearman],
        direction_stable=[stable], pvalue=[pvalue], qvalue=[qvalue], selected=[selected],
    )
    vif = DataFrame(
        variable=["ndvi"], vif=Union{Missing,Float64}[missing],
        status=["not_applicable_single_predictor"],
    )
    monthly = DataFrame(
        variable=String[], month=Int[], pearson=Float64[], spearman=Float64[],
        pearson_direction=String[], spearman_direction=String[],
    )
    for month_value in 6:9
        columns = findall(period -> month(period) == month_value, periods)
        if isempty(columns)
            p, s = NaN, NaN
        else
            xv, yv, weights = ERA5VariableSelection._flatten_balanced(
                panel.x[:ndvi][:, columns], panel.y[:, columns],
            )
            p = ERA5VariableSelection._weighted_correlation(xv, yv, weights)
            s = ERA5VariableSelection._weighted_correlation(
                ERA5VariableSelection._tie_ranks(xv), ERA5VariableSelection._tie_ranks(yv), weights,
            )
        end
        push!(monthly, ("ndvi", month_value, p, s, _direction(p), _direction(s)))
    end
    return (; association, vif, monthly, selected=selected ? [:ndvi] : Symbol[])
end

function ndvi_spatial_variability_test(
    panel, selected::Vector{Symbol}, lonlat::Matrix{Float64}, cfg::NDVISelectionConfig;
    rng::AbstractRNG=MersenneTwister(cfg.seed + 1),
)
    return ERA5VariableSelection.panel_spatial_variability_test(
        panel, selected, lonlat, _era_config(cfg); rng,
    )
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

function run_ndvi_variable_selection(
    cfg::NDVISelectionConfig, products::Vector{String}, station_ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64}, satellite::AbstractDict,
    ndvi_table::DataFrame, lonlat::Matrix{Float64};
    data_qc::Union{Nothing,DataFrame}=nothing,
    station_qc::Union{Nothing,DataFrame}=nothing,
)
    mkpath(cfg.outdir)
    n = length(station_ids)
    size(Yobs) == (n, length(times)) || throw(DimensionMismatch("observation dimensions differ"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("coordinate dimensions differ"))
    aligned = align_ndvi_asof(ndvi_table, station_ids, times; max_age_days=cfg.max_age_days)
    folds = ERA5VariableSelection.balanced_spatial_folds(station_ids, lonlat; k=cfg.k)
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), DataFrame(
        station_id=station_ids, fold=folds, lon=lonlat[:, 1], lat=lonlat[:, 2],
    ))
    data_qc !== nothing && CSV.write(joinpath(cfg.outdir, "ndvi_data_qc.csv"), data_qc)
    station_qc !== nothing && CSV.write(joinpath(cfg.outdir, "ndvi_station_qc.csv"), station_qc)
    CSV.write(joinpath(cfg.outdir, "ndvi_alignment_qc.csv"), aligned.qc)
    CSV.write(joinpath(cfg.outdir, "ndvi_alignment_intervals.csv"), aligned.intervals)

    qc_all, cells_all, association_all, monthly_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    vif_all, bandwidth_all, variability_all, roles_all, status_all =
        DataFrame(), DataFrame(), DataFrame(), DataFrame(), DataFrame()
    specifications = DataFrame(
        product=String[], variable=String[], selected=Bool[], role=String[],
        association_qvalue=Float64[], spatial_qvalue=Union{Missing,Float64}[], status=String[],
    )
    schemes = [("full_data", 0, collect(1:n))]
    append!(schemes, [("spatial_cv", fold, findall(!=(fold), folds)) for fold in 1:cfg.k])
    for product in products
        Ysat = Float64.(satellite[product])
        for (scheme, fold, train_indices) in schemes
            run_seed = cfg.seed + 1000 * findfirst(==(product), products) + 10 * fold
            panel = prepare_ndvi_panel(Yobs, Ysat, aligned, train_indices, station_ids, cfg)
            qc = copy(panel.qc); _annotate!(qc, product, scheme, fold); _append_table!(qc_all, qc)
            cells = copy(panel.cells); _annotate!(cells, product, scheme, fold); _append_table!(cells_all, cells)
            if panel.qc.status[1] != "ok"
                push!(roles_all, (product=product, scheme=scheme, fold=fold,
                    variable="ndvi", selected=false, role="uncertain",
                    status="insufficient_stations_or_variation"))
                if scheme == "full_data"
                    push!(specifications, (product, "ndvi", false, "uncertain", NaN,
                        missing, "insufficient_stations_or_variation"))
                end
                push!(status_all, (product=product, scheme=scheme, fold=fold,
                    status="failed", reason="insufficient_stations_or_variation"))
                continue
            end
            screen = screen_ndvi(panel, aligned.periods, cfg; rng=MersenneTwister(run_seed))
            for table in (screen.association, screen.monthly, screen.vif)
                _annotate!(table, product, scheme, fold)
            end
            _append_table!(association_all, screen.association)
            _append_table!(monthly_all, screen.monthly)
            _append_table!(vif_all, screen.vif)
            spatial = ndvi_spatial_variability_test(
                panel, screen.selected, lonlat[train_indices, :], cfg;
                rng=MersenneTwister(run_seed + 1),
            )
            _annotate!(spatial.bandwidth_scan, product, scheme, fold)
            _annotate!(spatial.variability, product, scheme, fold)
            _append_table!(bandwidth_all, spatial.bandwidth_scan)
            _append_table!(variability_all, spatial.variability)
            if isempty(screen.selected)
                role, spatial_q, role_status = "not_selected", missing, "not_selected"
            elseif nrow(spatial.variability) == 0
                role, spatial_q, role_status = "uncertain", missing, "test_failed"
            else
                row = spatial.variability[1, :]
                role = row.role
                spatial_q = isfinite(row.qvalue) ? row.qvalue : missing
                role_status = row.status
            end
            push!(roles_all, (product=product, scheme=scheme, fold=fold,
                variable="ndvi", selected=!isempty(screen.selected), role, status=role_status))
            if scheme == "full_data"
                push!(specifications, (product, "ndvi", !isempty(screen.selected), role,
                    screen.association.qvalue[1], spatial_q, role_status))
            end
            push!(status_all, (product=product, scheme=scheme, fold=fold,
                status="ok", reason=isempty(screen.selected) ? "no_variable_selected" : ""))
        end
    end
    stability = DataFrame(
        product=String[], variable=String[], selection_rate=Float64[],
        local_rate=Float64[], global_rate=Float64[], uncertain_rate=Float64[],
    )
    for product in products
        rows = filter([:product, :scheme] => (p, s) -> p == product && s == "spatial_cv", roles_all)
        denominator = max(nrow(rows), 1)
        push!(stability, (product, "ndvi", count(rows.selected) / denominator,
            count(==("local"), rows.role) / denominator,
            count(==("global"), rows.role) / denominator,
            count(==("uncertain"), rows.role) / denominator))
    end
    selected_count = count(specifications.selected)
    local_count = count(==("local"), specifications.role)
    global_count = count(==("global"), specifications.role)
    consensus_role = selected_count == 0 ? "not_selected" : local_count >= 2 ? "local" :
        global_count >= 2 ? "global" : selected_count >= 2 ? "role_unstable" : "product_specific"
    consensus = DataFrame(
        variable=["ndvi"], selected_product_count=[selected_count],
        consensus_selected=[selected_count >= 2], local_product_count=[local_count],
        global_product_count=[global_count], consensus_role=[consensus_role],
    )
    outputs = Dict(
        "ndvi_fold_quality_control.csv" => qc_all,
        "ndvi_station_period_residual.csv" => cells_all,
        "ndvi_fold_association.csv" => association_all,
        "ndvi_fold_monthly_direction.csv" => monthly_all,
        "ndvi_fold_vif.csv" => vif_all,
        "ndvi_fold_bandwidth_scan.csv" => bandwidth_all,
        "ndvi_fold_spatial_variability.csv" => variability_all,
        "ndvi_fold_roles.csv" => roles_all,
        "ndvi_role_stability.csv" => stability,
        "ndvi_final_full_data_spec.csv" => specifications,
        "ndvi_cross_product_consensus.csv" => consensus,
        "run_status.csv" => status_all,
    )
    for (filename, table) in outputs
        CSV.write(joinpath(cfg.outdir, filename), table)
    end
    return (; specifications, stability, consensus, status=status_all)
end

end
