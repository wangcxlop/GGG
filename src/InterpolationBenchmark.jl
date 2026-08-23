using CSV, DataFrames, Dates, LinearAlgebra, Random, Statistics
using MixedGWR

if !isdefined(@__MODULE__, :MGERConfig)
    include(joinpath(@__DIR__, "MGERPipeline.jl"))
end
if !isdefined(@__MODULE__, :TraditionalInterpolation)
    include(joinpath(@__DIR__, "TraditionalInterpolation.jl"))
end
using .TraditionalInterpolation
if !isdefined(@__MODULE__, :DEMTerrainExperiment)
    include(joinpath(@__DIR__, "DEMTerrainExperiment.jl"))
end
using .DEMTerrainExperiment: DEMExperimentConfig, terrain_screen, spatial_variability_test
using .DEMTerrainExperiment: mean_wet_residual, monthly_correlation_rows
using .DEMTerrainExperiment: terrain_model_designs, terrain_groups, terrain_columns
using .DEMTerrainExperiment: mixed_gwr_predict, multiscale_gwr_predict
using .DEMTerrainExperiment: select_mixed_bandwidth, select_multiscale_bandwidths
if !isdefined(@__MODULE__, :JointCovariateModels)
    include(joinpath(@__DIR__, "JointCovariateModels.jl"))
end
using .JointCovariateModels
if !isdefined(@__MODULE__, :ERA5VariableSelection)
    include(joinpath(@__DIR__, "ERA5VariableSelection.jl"))
end
if !isdefined(@__MODULE__, :NDVIVariableSelection)
    include(joinpath(@__DIR__, "NDVIVariableSelection.jl"))
end
if !isdefined(@__MODULE__, :JointVariableSelection)
    include(joinpath(@__DIR__, "JointVariableSelection.jl"))
end
using .JointVariableSelection: JointSelectionConfig, select_joint_covariates

const BENCHMARK_METHODS = [
    "raw", "zero", "train_clim", "hour_field_mean",
    "idw", "adw", "tps", "gwr", "gwr_const",
    "residual_gwr", "residual_gwr_const", "mixed_gwr", "mgwr", "hurdle_gwr",
]
const TRADITIONAL_METHODS = ["idw", "adw", "tps"]
# Hyperparameter-free reference predictors, estimated from training stations only. They are
# reported so a method's RMSE can be read as skill rather than as a bare number: on hourly
# precipitation ~93% of station-hours are dry, so a large part of any RMSE is simply how hard
# a method shrinks toward zero. Deliberately not in TRADITIONAL_METHODS - they are context for
# the reader, not baselines the GWR claim is assessed against.
const NULL_METHODS = ["zero", "train_clim", "hour_field_mean"]
# The common evaluation mask stays pinned to the original comparison set. Every method must be
# finite for a cell to count, so letting a newly added method into the mask would silently
# re-score all the others and break comparability with earlier runs - the same coverage
# confound that already makes the `balanced_spatial` numbers hard to read (direct `gwr` fails
# often enough to drag the shared mask down to ~0.76). Diagnostic and reference methods are
# therefore scored *on* the mask without being allowed to *define* it.
const MASK_METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]
const BENCHMARK_RUNS = [
    ("direct", "idw"), ("direct", "adw"), ("direct", "tps"),
    ("direct", "gwr"), ("direct", "gwr_const"), ("direct", "hurdle_gwr"),
    ("residual", "gwr"), ("residual", "gwr_const"), ("residual", "mixed_gwr"),
    ("residual", "mgwr"),
]

Base.@kwdef struct InterpolationBenchmarkConfig
    mger::MGERConfig
    terrain_path::Union{Nothing,String} = nothing
    dem::Union{Nothing,DEMExperimentConfig} = nothing
    joint_covariates::Union{Nothing,JointCovariateBenchmarkConfig} = nothing
    # Set to run joint-covariate variable selection fresh inside every training fold instead of
    # once on the full station set (`joint_covariates.spec_path`). Exactly one of the two must
    # be set when `joint_covariates` is enabled.
    joint_selection::Union{Nothing,JointSelectionConfig} = nothing
    k::Int = 5
    seed::Int = 20260627
    # Repeated cross-validation: one independent fold partition per seed. Empty means a
    # single repeat using `seed`, which reproduces the historical single-split behaviour.
    seeds::Vector{Int} = Int[]
    fold_center_init::Symbol = :kmeanspp
    cv_schemes::Vector{Symbol} = [:balanced_spatial, :random]
    idw_powers::Vector{Float64} = [1.0, 1.5, 2.0, 2.5, 3.0]
    neighbor_candidates::Vector{Union{Nothing,Int}} = Union{Nothing,Int}[8, 16, 32, 64, nothing]
    tps_smooth_candidates::Vector{Float64} = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    min_tuning_coverage::Float64 = 0.95
    tuning_max_times::Int = 336
    # `:stratified` weights the tuning subsample so its RMSE estimates the pooled RMSE the
    # benchmark reports, instead of the wet-hour RMSE the unweighted subsample scores.
    #
    # It is NOT the default, despite being the more correct estimator, because measurement says
    # the level error it removes was very nearly a constant multiplier and therefore cancelled
    # in the ranking: over 30 product/fold/method cells it costs 0.586% mean regret against the
    # exact full-record curve where `:uniform` costs 0.096%, and on the joint-covariate path the
    # two are a wash (63.4% vs 63.7%). Switching the default would move published numbers for no
    # measured gain. See `scripts/verify_tuning_time_weighting.jl`, and the header of
    # `_tuning_time_sample` for what actually dominates selection error.
    tuning_time_weighting::Symbol = :uniform
    # How candidates are scored during selection. `:inner_spatial` predicts out-of-fold onto an
    # inner split of the training stations, built by the same splitter and scheme as the outer
    # partition, so the selection criterion is the same estimand the benchmark reports.
    # `:loocv` is the historical leave-one-out-at-training-stations criterion, which scores
    # interpolation skill (median 5.5 km to the nearest remaining gauge) while the results table
    # reports extrapolation skill (median 23.9 km under `balanced_spatial`). See
    # `selection_folds` for the measurements.
    tuning_geometry::Symbol = :inner_spatial
    # Groups in the inner selection split; 0 means "use `cfg.k`", which is what reproduces the
    # outer fold geometry. Ignored when `tuning_geometry === :loocv`.
    tuning_inner_k::Int = 0
    mgwr_max_tuning_iterations::Int = 5
    # Scale applied to the joint dynamic models' residual correction before it is added back to
    # the satellite field. GWR is unbiased and never shrinks, so the correction it produces is
    # the right shape at the wrong magnitude: `satellite_offset.csv` measures the optimal rescale
    # at 0.84 / 0.47 / 0.53 for FY4B / GPM / GSMaP, i.e. the models over-correct by up to 2x.
    #
    # The scale is selected jointly with the bandwidth on the inner spatial split, and costs
    # nothing to add: `_joint_candidate_metrics` rescores the residual prediction it has already
    # computed, so the ladder needs no extra model fits. Scoring the clipped objective rather than
    # the closed-form `E[g*r]/E[g^2]` keeps selection consistent with the reported metric, which
    # the `max(., 0)` clip would otherwise break.
    #
    # The ladder starts at 0.1 rather than 0.3 because a smoke run with a 0.3 floor put 31 of 75
    # selected joint rows *on* that floor - the grid, not the data, would have been picking the
    # scale, which is the same failure the bandwidth grids hit before `--local-grid`.
    #
    # `[1.0]` disables shrinkage and reproduces the historical unshrunk behaviour exactly.
    residual_shrinkage_candidates::Vector{Float64} =
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    event_thresholds::Vector{Float64} = [0.1, 2.5, 8.0, 16.0]
    bootstrap_reps::Int = 2000
end

"""Fold seeds for the run: one repeat per seed, defaulting to a single repeat on `cfg.seed`."""
benchmark_seeds(cfg::InterpolationBenchmarkConfig) =
    isempty(cfg.seeds) ? [cfg.seed] : cfg.seeds

function _validate_benchmark_config(cfg::InterpolationBenchmarkConfig, n_station::Int)
    2 <= cfg.k <= n_station || throw(ArgumentError("k must be between 2 and station count"))
    length(unique(benchmark_seeds(cfg))) == length(benchmark_seeds(cfg)) ||
        throw(ArgumentError("seeds must be unique"))
    cfg.fold_center_init in (:kmeanspp, :farthest) ||
        throw(ArgumentError("fold_center_init must be :kmeanspp or :farthest"))
    all(s -> s in (:balanced_spatial, :random, :strip), cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must contain :balanced_spatial, :random, or :strip"))
    length(unique(cfg.cv_schemes)) == length(cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must be unique"))
    all(>(0), cfg.idw_powers) || throw(ArgumentError("IDW/ADW powers must be positive"))
    all(x -> x === nothing || x > 0, cfg.neighbor_candidates) ||
        throw(ArgumentError("neighbor candidates must be positive or nothing"))
    all(>(0), cfg.tps_smooth_candidates) ||
        throw(ArgumentError("TPS smoothing candidates must be positive"))
    0 < cfg.min_tuning_coverage <= 1 ||
        throw(ArgumentError("min_tuning_coverage must be in (0, 1]"))
    cfg.bootstrap_reps >= 0 || throw(ArgumentError("bootstrap_reps must be non-negative"))
    cfg.tuning_max_times >= 0 || throw(ArgumentError("tuning_max_times must be non-negative"))
    cfg.tuning_time_weighting in (:stratified, :uniform) ||
        throw(ArgumentError("tuning_time_weighting must be :stratified or :uniform"))
    cfg.tuning_geometry in (:inner_spatial, :loocv) ||
        throw(ArgumentError("tuning_geometry must be :inner_spatial or :loocv"))
    cfg.tuning_inner_k == 0 || cfg.tuning_inner_k >= 2 ||
        throw(ArgumentError("tuning_inner_k must be 0 (use cfg.k) or at least 2"))
    cfg.mgwr_max_tuning_iterations > 0 ||
        throw(ArgumentError("mgwr_max_tuning_iterations must be positive"))
    isempty(cfg.residual_shrinkage_candidates) &&
        throw(ArgumentError("residual_shrinkage_candidates must not be empty"))
    all(s -> 0 < s <= 1, cfg.residual_shrinkage_candidates) ||
        throw(ArgumentError("residual shrinkage candidates must be in (0, 1]"))
    xor(cfg.terrain_path === nothing, cfg.dem === nothing) && throw(ArgumentError(
        "terrain_path and dem must either both be configured or both be omitted",
    ))
    cfg.joint_covariates !== nothing && cfg.terrain_path !== nothing && throw(ArgumentError(
        "joint covariates and legacy DEM screening cannot be enabled together",
    ))
    if cfg.dem !== nothing
        cfg.dem.min_wet_hours > 0 || throw(ArgumentError("DEM min_wet_hours must be positive"))
        cfg.dem.screen_permutations > 0 ||
            throw(ArgumentError("DEM screen_permutations must be positive"))
        cfg.dem.spatial_permutations > 0 ||
            throw(ArgumentError("DEM spatial_permutations must be positive"))
        all(>(1), cfg.dem.bandwidth_candidates) ||
            throw(ArgumentError("DEM bandwidth candidates must exceed one neighbor"))
    end
    if cfg.joint_covariates !== nothing
        joint = cfg.joint_covariates
        xor(joint.spec_path === nothing, cfg.joint_selection === nothing) || throw(ArgumentError(
            "joint covariates require exactly one of spec_path (fixed full-data selection) or " *
            "joint_selection (nested per-fold selection)",
        ))
        joint.spec_path === nothing || isfile(joint.spec_path) ||
            throw(ArgumentError("joint specification does not exist"))
        isfile(joint.terrain_path) || throw(ArgumentError("joint terrain table does not exist"))
        all(>(1), joint.bandwidth_candidates) ||
            throw(ArgumentError("joint bandwidth candidates must exceed one neighbor"))
        joint.max_iterations > 0 || throw(ArgumentError("joint max_iterations must be positive"))
        joint.tolerance > 0 || throw(ArgumentError("joint tolerance must be positive"))
        # Successive over-relaxation diverges outside (0, 2); see the field comment in
        # `JointCovariateBenchmarkConfig` for why the default sits above 1.
        0 < joint.relaxation < 2 ||
            throw(ArgumentError("joint relaxation must be in (0, 2)"))
        joint.mgwr_spatial_grouping in (:split, :shared, :intercept_only) || throw(ArgumentError(
            "mgwr_spatial_grouping must be :split, :shared, or :intercept_only",
        ))
    end
    return nothing
end

"""
Per-repeat output root. A single-repeat run keeps the historical `outdir/<scheme>/...` layout so
existing verify scripts and downstream readers are unaffected; repeated runs nest under
`outdir/repeat_<i>/`.
"""
_repeat_dir(outdir::String, repeat_index::Int, n_repeats::Int) =
    n_repeats == 1 ? outdir : joinpath(outdir, "repeat_$(lpad(repeat_index, 2, '0'))")

_dem_enabled(cfg::InterpolationBenchmarkConfig) = cfg.terrain_path !== nothing
_joint_enabled(cfg::InterpolationBenchmarkConfig) = cfg.joint_covariates !== nothing

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

"""Draw an index with probability proportional to `weights` (k-means++ sampling)."""
function _sample_weighted(weights::Vector{Float64}, rng::AbstractRNG)
    total = sum(weights)
    (isfinite(total) && total > 0) || return argmax(weights)
    threshold = rand(rng) * total
    cumulative = 0.0
    for index in eachindex(weights)
        cumulative += weights[index]
        cumulative >= threshold && return index
    end
    return lastindex(weights)
end

"""
Initial fold centers.

`:kmeanspp` draws each center with probability proportional to its squared distance to the
nearest existing center, so different seeds give genuinely different partitions. `:farthest`
is the original deterministic farthest-point traversal, where the seed only chose the first
center and therefore produced near-identical partitions across seeds.
"""
function _initial_spatial_centers(
    xy::Matrix{Float64}, k::Int, rng::AbstractRNG; center_init::Symbol=:kmeanspp,
)
    n = size(xy, 1)
    centers = Matrix{Float64}(undef, k, 2)
    first_index = rand(rng, 1:n)
    centers[1, :] = xy[first_index, :]
    nearest2 = fill(Inf, n)
    for center_index in 2:k
        previous = @view centers[center_index - 1, :]
        for station in 1:n
            d2 = sum(abs2, @view(xy[station, :]) .- previous)
            nearest2[station] = min(nearest2[station], d2)
        end
        next_index = center_init === :kmeanspp ? _sample_weighted(nearest2, rng) :
            argmax(nearest2)
        centers[center_index, :] = xy[next_index, :]
    end
    return centers
end

function _capacity_assignment(xy::Matrix{Float64}, centers::Matrix{Float64}, capacities::Vector{Int})
    n, k = size(xy, 1), size(centers, 1)
    distances = Matrix{Float64}(undef, n, k)
    for station in 1:n, cluster in 1:k
        distances[station, cluster] = sum(abs2, @view(xy[station, :]) .- @view(centers[cluster, :]))
    end
    certainty = Vector{Float64}(undef, n)
    for station in 1:n
        ordered = partialsort(@view(distances[station, :]), 1:min(2, k))
        certainty[station] = length(ordered) == 1 ? Inf : ordered[2] - ordered[1]
    end
    order = sortperm(1:n; by=station -> (-certainty[station], station))
    remaining = copy(capacities)
    assignments = zeros(Int, n)
    for station in order
        cluster_order = sortperm(1:k; by=cluster -> (distances[station, cluster], cluster))
        chosen = findfirst(cluster -> remaining[cluster] > 0, cluster_order)
        chosen === nothing && error("capacity assignment failed")
        cluster = cluster_order[chosen]
        assignments[station] = cluster
        remaining[cluster] -= 1
    end
    return assignments
end

"""Reproducible, compact, capacity-balanced two-dimensional spatial folds."""
function split_stations_balanced_spatial_kfold(
    station_ids::Vector{String}, lonlat::Matrix{Float64}; k::Int=5, seed::Int=20260627,
    center_init::Symbol=:kmeanspp,
)
    n = length(station_ids)
    2 <= k <= n || throw(ArgumentError("k must be between 2 and station count"))
    size(lonlat, 1) == n || throw(DimensionMismatch("lonlat rows must match station_ids"))
    xy = local_km_coordinates(lonlat)
    capacities = fill(div(n, k), k)
    capacities[1:rem(n, k)] .+= 1
    centers = _initial_spatial_centers(xy, k, MersenneTwister(seed); center_init)
    assignments = zeros(Int, n)
    for _ in 1:100
        updated = _capacity_assignment(xy, centers, capacities)
        updated == assignments && break
        assignments = updated
        for cluster in 1:k
            idx = findall(==(cluster), assignments)
            centers[cluster, :] = vec(mean(xy[idx, :], dims=1))
        end
    end
    folds = [String[] for _ in 1:k]
    for station in 1:n
        push!(folds[assignments[station]], station_ids[station])
    end
    return folds
end

function benchmark_folds(
    scheme::Symbol, ids::Vector{String}, lonlat::Matrix{Float64}; k::Int, seed::Int,
    center_init::Symbol=:kmeanspp,
)
    if scheme == :balanced_spatial
        return split_stations_balanced_spatial_kfold(ids, lonlat; k=k, seed=seed, center_init)
    elseif scheme == :random
        return split_stations_kfold(ids; k=k, rng=MersenneTwister(seed))
    elseif scheme == :strip
        return split_stations_spatial_block_kfold(ids, lonlat; k=k)
    end
    throw(ArgumentError("unsupported CV scheme: $scheme"))
end

"""
Inner split of one fold's training stations, as positions within that training set.

Hyperparameters have to be chosen by the same kind of prediction the benchmark reports. The
historical criterion was leave-one-out at training stations, which leaves the held-out station
sitting inside its own neighbourhood: median 5.5 km to the nearest remaining gauge, against
23.9 km from an outer `balanced_spatial` fold station to its nearest training gauge. The two are
not the same estimand, and on the joint-covariate path the reported RMSE curve falls
monotonically across the whole bandwidth grid while the leave-one-out curve is U-shaped with a
minimum at 12-30 — over most of the grid they are anti-correlated, and the leave-one-out pick
costs +19% to +101% of reported RMSE.

Splitting the training stations with the same splitter, the same scheme and `k_inner = cfg.k`
reproduces the outer geometry closely (inner-val nearest-train km p10/median/p90 = 8.1/22.7/45.6
against the outer 8.6/23.9/48.9); `k_inner=3` overshoots to a 32.2 km median and `k_inner=8`
undershoots to 17.5, so matching `cfg.k` is both the knob-free and the accurate choice.

Using `scheme` rather than always splitting spatially is deliberate: under `:random` the
reporting geometry (median 6.3 km) already matches leave-one-out, so a random inner split
correctly leaves that scheme's numbers essentially where they were.

Returns `nothing` when the fold is too small to split — holding out a whole group would leave
too few stations to fit the widest design (`mixed_gwr`/`mgwr` need three local columns plus a
global block). The caller falls back to leave-one-out and records that in `benchmark_scope.csv`,
so a degenerate run never silently reports one criterion while having used the other. At the real
station count this never triggers: 189 training stations split five ways leave 151.
"""
const MIN_SELECTION_TRAIN_STATIONS = 10

function selection_folds(
    cfg::InterpolationBenchmarkConfig, scheme::Symbol, train_ids::Vector{String},
    train_lonlat::Matrix{Float64}, fold::Int, repeat_seed::Int,
)
    k_inner = cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k
    n_train = length(train_ids)
    k_inner = min(k_inner, n_train)
    k_inner >= 2 || return nothing
    n_train - cld(n_train, k_inner) >= MIN_SELECTION_TRAIN_STATIONS || return nothing
    # Vary by outer fold and repeat so the inner partition is not the same one every time. Not by
    # product: the station geometry is identical across products, and holding the split fixed
    # keeps the per-product comparison clean.
    groups = benchmark_folds(scheme, train_ids, train_lonlat;
        k=k_inner, seed=repeat_seed + 7919 * fold, center_init=cfg.fold_center_init)
    position = Dict(id => index for (index, id) in enumerate(train_ids))
    return [[position[id] for id in group] for group in groups if !isempty(group)]
end

"""
Assemble out-of-fold predictions over the training stations.

`predict(inner_train_rows, inner_val_rows)` returns a `length(inner_val_rows) × n_time` matrix.
The result has the same shape leave-one-out produced, so `_candidate_metrics` and every scan row
downstream are unchanged.
"""
function _selection_oof(
    groups::Vector{Vector{Int}}, n_train::Int, n_time::Int, predict,
)
    out_of_fold = fill(NaN, n_train, n_time)
    for group in groups
        inner_train = setdiff(1:n_train, group)
        isempty(inner_train) && continue
        out_of_fold[group, :] = predict(inner_train, group)
    end
    return out_of_fold
end

function _gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64}, target_lonlat::Matrix{Float64};
    kernel::Int, adaptive::Bool, bw::Float64, exclude_self::Bool=false,
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center=center)
    X_target = build_X_intercept_centered(target_lonlat; center=center)
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("GWR exclude_self requires matching rows"))
        for i in 1:size(distances, 1)
            distances[i, i] = Inf
        end
    end
    weights = gw_weight(distances, bw; kernel=kernel, adaptive=adaptive)
    return st_gwr_predict_nanaware(X_train, values, weights; Xpred=X_target, min_obs=3)
end

"""
Local-constant GWR: the same kernel and bandwidth machinery as `_gwr_predict`, but fitting a
local mean instead of a local plane.

This exists to separate two explanations that are otherwise confounded. `_gwr_predict` fits
`[1, lon, lat]`, so it needs at least three effectively-weighted neighbours and cannot be as
local as IDW/ADW, which need two. Running both over the same bandwidth grid says whether GWR
loses because of the bandwidth it was given or because of the order of the local model.

Unlike `_gwr_predict` this cannot go through `st_gwr_predict_nanaware`, which hard-requires a
three-column design. The NaN handling mirrors `TraditionalInterpolation._weighted_predict`:
renormalize by the weights of the stations that actually reported each hour.
"""
function _gwr_const_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64}, target_lonlat::Matrix{Float64};
    kernel::Int, adaptive::Bool, bw::Float64, exclude_self::Bool=false, min_obs::Int=1,
)
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("GWR exclude_self requires matching rows"))
        for i in 1:size(distances, 1)
            distances[i, i] = Inf
        end
    end
    weights = gw_weight(distances, bw; kernel=kernel, adaptive=adaptive)
    n_target = size(target_lonlat, 1)
    n_time = size(values, 2)
    prediction = fill(NaN, n_target, n_time)
    Threads.@threads for target in 1:n_target
        @inbounds for time in 1:n_time
            total = 0.0
            weight_sum = 0.0
            effective = 0
            for station in axes(values, 1)
                weight = weights[station, target]
                value = values[station, time]
                (weight > 0 && isfinite(weight) && !isnan(value)) || continue
                total += weight * value
                weight_sum += weight
                effective += 1
            end
            if effective >= min_obs && weight_sum > 0
                prediction[target, time] = total / weight_sum
            end
        end
    end
    return prediction
end

"""
Per-fold context for `hurdle_gwr`.

`fit_hurdle_gwr`/`predict_hurdle_gwr` need the fold's `times` vector, which the rest of the
benchmark's method dispatch never has to carry; this mirrors the existing `dem_context` /
`joint_context` keyword pattern rather than widening every signature. `diagnostics` is a shared
accumulator (same pattern as `scan_rows`) for the degeneracy report described at
`_hurdle_predict`.
"""
function build_hurdle_context(
    times::Vector{DateTime}, cfg::InterpolationBenchmarkConfig,
    diagnostics::Vector{NamedTuple}; scheme::String, product::String, fold::Int,
)
    wet_threshold = cfg.mger.rain_threshold
    # `HurdleGWRConfig` requires `first(intensity_breaks) == wet_threshold`, sorted and unique.
    # Deriving the breaks from the benchmark's own event thresholds keeps the occurrence table
    # binned the same way the categorical scores are.
    breaks = sort(unique(vcat(wet_threshold, filter(>(wet_threshold), cfg.event_thresholds))))
    return (; times, wet_threshold, intensity_breaks=breaks, diagnostics,
        scheme, product, fold)
end

"""Adapt the benchmark's `(adaptive, bw)` convention to HurdleGWR's `SpatialBandwidth`."""
_hurdle_bandwidth(adaptive::Bool, bw::Float64) =
    adaptive ? AdaptiveBandwidth(Int(round(bw))) : FixedBandwidth(bw)

_hurdle_config(context, kernel::Int, adaptive::Bool, bw::Float64) = HurdleGWRConfig(
    wet_threshold=context.wet_threshold,
    intensity_breaks=context.intensity_breaks,
    kernel=kernel,
    bandwidth=_hurdle_bandwidth(adaptive, bw),
)

"""
Fit the hurdle model on the training fold and predict `P(wet) * E[amount | wet]`.

Unlike every other method here this is a time-pooled calibration, not a per-hour spatial fit:
one local regression per target station, with time entering through annual/diurnal harmonics and
a calendar-month occurrence bin. The satellite enters the design as a *fitted* `log1p_satellite`
coefficient rather than as an offset with its coefficient forced to 1.

`predict_hurdle_gwr` takes no gauge observations at the target at all, so held-out leakage is
impossible by construction. `exclude_self` drives its `exclude_train_indices`, which zeroes a
station's own spatial weight - the same leave-one-out convention `_gwr_predict` uses.

Returns the full prediction object: the caller needs `used_global_amount` and
`occurrence_fallback` as well as `corrected`, because when a local fit misses
`min_amount_samples` / `min_effective_stations` the model *silently* falls back to global
coefficients, and an undiagnosed fallback would report a global model as a local one.
"""
function _hurdle_predict(
    train_lonlat::Matrix{Float64}, y_obs_train::Matrix{Float64},
    y_sat_train::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_sat_target::Matrix{Float64}, times::Vector{DateTime};
    config::HurdleGWRConfig, exclude_self::Bool=false,
)
    model = fit_hurdle_gwr(y_obs_train, y_sat_train, train_lonlat, times; config)
    exclusions = exclude_self ? collect(1:size(train_lonlat, 1)) : nothing
    return predict_hurdle_gwr(
        model, y_sat_target, target_lonlat, times; exclude_train_indices=exclusions,
    )
end

"""Summarise how much of a hurdle prediction actually came from a *local* fit."""
function _hurdle_diagnostic_row(prediction, context; kernel::Int, adaptive::Bool, bw::Float64)
    fallback = prediction.occurrence_fallback
    total = length(fallback)
    share(code) = total == 0 ? NaN : count(==(code), fallback) / total
    return (;
        scheme=context.scheme, product=context.product, fold=context.fold,
        kernel, adaptive, bw,
        n_target=length(prediction.used_global_amount),
        used_global_amount_fraction=mean(prediction.used_global_amount),
        mean_effective_stations=mean(prediction.effective_stations),
        occurrence_monthly=share(0x00), occurrence_all_month=share(0x01),
        occurrence_global=share(0x02), occurrence_skipped=share(0x03),
    )
end

"""
Hyperparameter-free reference predictors for one fold, estimated from training stations only.

- `zero`: constant 0. Hourly precipitation is ~93% dry, so this is a surprisingly strong
  RMSE competitor and the floor any method must clear.
- `train_clim`: the training stations' overall mean. A per-station climatology is not
  estimable for a held-out station, so the pooled mean is the station-free analogue.
- `hour_field_mean`: each hour's spatial mean over training stations - "is it raining
  anywhere in the domain right now", carrying no spatial structure at all. An interpolator
  that does not beat this is not doing spatial work.

Reported alongside the real methods so `metrics_*.csv` can be read as skill rather than as a
bare RMSE. See `NULL_METHODS`.
"""
function _null_fold_predictions(y_obs_train::Matrix{Float64}, n_val::Int)
    n_time = size(y_obs_train, 2)
    finite = filter(isfinite, vec(y_obs_train))
    climatology = isempty(finite) ? NaN : mean(finite)
    field_mean = Matrix{Float64}(undef, n_val, n_time)
    @inbounds for time in 1:n_time
        total = 0.0
        count = 0
        for station in axes(y_obs_train, 1)
            value = y_obs_train[station, time]
            if !isnan(value)
                total += value
                count += 1
            end
        end
        field_mean[:, time] .= count > 0 ? total / count : NaN
    end
    return Dict{String,Matrix{Float64}}(
        "zero" => zeros(Float64, n_val, n_time),
        "train_clim" => fill(climatology, n_val, n_time),
        "hour_field_mean" => field_mean,
    )
end

"""Build the provisional mixed-GWR design; replace this when variable roles are finalized."""
function build_mixed_gwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    local_train = build_X_intercept_centered(train_lonlat; center)
    local_target = build_X_intercept_centered(target_lonlat; center)
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    return (; local_train, local_target, global_train, global_target)
end

function _mixed_gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bw::Float64, exclude_self::Bool=false,
)
    designs = build_mixed_gwr_designs(train_lonlat, target_lonlat)
    prediction, _ = mixed_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, Int(round(bw));
        exclude_self,
    )
    return prediction
end

"""Build the provisional MGWR groups; replace these when covariates are finalized."""
function build_mgwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center)
    X_target = build_X_intercept_centered(target_lonlat; center)
    local_train = [X_train[:, index:index] for index in axes(X_train, 2)]
    local_target = [X_target[:, index:index] for index in axes(X_target, 2)]
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    group_names = ["intercept", "longitude", "latitude"]
    return (; local_train, local_target, global_train, global_target, group_names)
end

function _mgwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bandwidths::Vector{Int}, exclude_self::Bool=false,
)
    designs = build_mgwr_designs(train_lonlat, target_lonlat)
    prediction, _ = multiscale_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, bandwidths;
        exclude_self,
    )
    return prediction
end

"""
Score one tuning candidate.

`time_weights`, when given, holds one weight per column of `y_obs` and makes RMSE/MAE weighted
means, so a stratified tuning subsample (see [`_tuning_time_sample`](@ref)) estimates the metric
over the full hour set rather than over the subsample's own wet-heavy mixture. `n` and
`coverage` stay raw cell counts either way: they only gate whether the candidate produced
predictions at all (`min_tuning_coverage`), which reweighting would silently redefine.
"""
function _candidate_metrics(
    y_obs::Matrix{Float64}, y_sat::Matrix{Float64}, prediction::Matrix{Float64};
    require_satellite::Bool=true, time_weights::Union{Nothing,Vector{Float64}}=nothing,
)
    base_mask = require_satellite ? (.!isnan.(y_obs) .& .!isnan.(y_sat)) : .!isnan.(y_obs)
    mask = base_mask .& .!isnan.(prediction)
    base_n = count(base_mask)
    n = count(mask)
    n == 0 && return (; n=0, coverage=0.0, RMSE=Inf, MAE=Inf)
    if time_weights === nothing
        metric = metric_continuous(y_obs, prediction; mask=mask)
        return (; n, coverage=n / base_n, RMSE=metric.RMSE, MAE=metric.MAE)
    end
    length(time_weights) == size(y_obs, 2) || throw(ArgumentError(
        "time_weights must have one entry per tuning hour",
    ))
    weighted_square = 0.0
    weighted_absolute = 0.0
    weight_total = 0.0
    for time in axes(y_obs, 2)
        weight = time_weights[time]
        for station in axes(y_obs, 1)
            mask[station, time] || continue
            error = prediction[station, time] - y_obs[station, time]
            weighted_square += weight * error^2
            weighted_absolute += weight * abs(error)
            weight_total += weight
        end
    end
    weight_total > 0 || return (; n, coverage=n / base_n, RMSE=Inf, MAE=Inf)
    return (; n, coverage=n / base_n,
        RMSE=sqrt(weighted_square / weight_total), MAE=weighted_absolute / weight_total)
end

function _scan_row(; scheme, product, fold, mode, method, group="all", iteration=0,
    repeat=1, seed=0,
    power=NaN, neighbors=missing,
    smooth=NaN, kernel=missing, adaptive=missing, bw=NaN, shrink=NaN, n=0, coverage=0.0,
    RMSE=Inf, MAE=Inf, status="success", error="", selected=false)
    stored_neighbors = neighbors === nothing ? 0 : neighbors
    return (;
        scheme=string(scheme), product=String(product), fold=Int(fold), mode=String(mode),
        method=String(method), group=String(group), iteration=Int(iteration),
        repeat=Int(repeat), seed=Int(seed),
        power=Float64(power), neighbors=stored_neighbors, smooth=Float64(smooth),
        kernel, adaptive, bw=Float64(bw), shrink=Float64(shrink),
        n=Int(n), coverage=Float64(coverage),
        RMSE=Float64(RMSE), MAE=Float64(MAE), status=String(status), error=String(error),
        selected=Bool(selected),
    )
end

function _select_candidate!(rows::Vector{NamedTuple}, first_row::Int, min_coverage::Float64)
    valid = [i for i in first_row:length(rows) if
        rows[i].status == "success" && rows[i].coverage >= min_coverage && isfinite(rows[i].RMSE)]
    isempty(valid) && throw(ArgumentError("all parameter candidates failed or had insufficient coverage"))
    best_index = sort(valid; by=i -> (rows[i].RMSE, rows[i].MAE, -rows[i].coverage))[1]
    rows[best_index] = merge(rows[best_index], (; selected=true))
    return rows[best_index]
end

function select_mgwr_bandwidths!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    scheme::Symbol, product::String, fold::Int, train_lonlat::Matrix{Float64},
    residuals::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    time_weights::Union{Nothing,Vector{Float64}}=nothing, selection_groups=nothing,
)
    n_train, n_time = size(residuals)
    candidates = sort(unique(Int(round(bw)) for bw in cfg.mger.bw_adaptive))
    isempty(candidates) && throw(ArgumentError("MGWR requires adaptive bandwidth candidates"))
    designs = build_mgwr_designs(train_lonlat, train_lonlat)
    bandwidths = fill(last(candidates), length(designs.group_names))
    first_row = length(scan_rows) + 1
    final_iteration = 0
    converged = false

    for iteration in 1:cfg.mgwr_max_tuning_iterations
        previous = copy(bandwidths)
        for (group_index, group) in enumerate(designs.group_names)
            candidate_rows = Int[]
            for bw in candidates
                trial = copy(bandwidths)
                trial[group_index] = bw
                try
                    interpolated = selection_groups === nothing ?
                        _mgwr_predict(
                            train_lonlat, residuals, train_lonlat;
                            bandwidths=trial, exclude_self=true,
                        ) :
                        _selection_oof(selection_groups, n_train, n_time, (tr, va) ->
                            _mgwr_predict(
                                train_lonlat[tr, :], residuals[tr, :], train_lonlat[va, :];
                                bandwidths=trial,
                            ))
                    prediction = max.(y_sat .+ interpolated, 0.0)
                    metrics = _candidate_metrics(y_obs, y_sat, prediction; time_weights)
                    push!(scan_rows, _scan_row(;
                        scheme, product, fold, mode="residual", method="mgwr",
                        group, iteration, kernel=BISQUARE, adaptive=true, bw, metrics...,
                    ))
                catch e
                    push!(scan_rows, _scan_row(;
                        scheme, product, fold, mode="residual", method="mgwr",
                        group, iteration, kernel=BISQUARE, adaptive=true, bw,
                        status="failed", error=sprint(showerror, e),
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
                "all MGWR candidates failed for group $group or had insufficient coverage",
            ))
            best_index = sort(valid; by=index -> (
                scan_rows[index].RMSE, scan_rows[index].MAE, -scan_rows[index].coverage,
            ))[1]
            bandwidths[group_index] = Int(round(scan_rows[best_index].bw))
        end
        final_iteration = iteration
        if bandwidths == previous
            converged = true
            break
        end
    end
    converged || throw(ArgumentError(
        "MGWR bandwidth search did not converge in $(cfg.mgwr_max_tuning_iterations) iterations",
    ))

    for (group_index, group) in enumerate(designs.group_names)
        matching_rows = filter(first_row:length(scan_rows)) do index
            row = scan_rows[index]
            row.method == "mgwr" && row.group == group &&
                row.iteration == final_iteration && row.status == "success" &&
                Int(round(row.bw)) == bandwidths[group_index]
        end
        isempty(matching_rows) && error("MGWR selected bandwidth row is missing for $group")
        selected_index = last(matching_rows)
        scan_rows[selected_index] = merge(scan_rows[selected_index], (; selected=true))
    end
    return (; bandwidths)
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

"""
Two-stratum tuning-hour sample and its Horvitz-Thompson weights.

Hyperparameters are chosen on this subsample but the benchmark reports pooled RMSE over every
hour, and the two populations are nothing alike: on the 11426 common hours the wettest half of
the subsample drives the tuning set to 40% wet cells against 7% in the reported population, and
13x its mean squared observation. An unweighted RMSE over the subsample therefore scores
wet-hour skill while the results table scores mostly-dry-hour skill.

So the wettest `maximum_times ÷ 8` hours are taken with certainty and the rest of the hours are
systematically subsampled and up-weighted by the inverse of their inclusion probability. The
weighted MSE is then an unbiased estimator of the pooled MSE over all `size(y_obs, 2)` hours.

The certainty stratum is an eighth of the budget rather than the historical half because under
weighting it only carries its true population share of the weight (~1.5%); its job is to cover
the extreme-error tail, and every hour beyond that is budget taken away from the weighted
remainder, which is where almost all of the weight — and therefore the estimator's variance —
actually sits. Measured over 12 product/fold/method cells against the exact full-record curve,
an eighth reproduces the legacy objective's candidate choice in 12/12 cells at 0.0% regret while
cutting the estimator's bias from +197% to -2.5%; a half costs 0.7-0.9% mean regret, and
dropping the certainty stratum entirely costs 0.4%. See
`scripts/verify_tuning_time_weighting.jl`.

`weighting=:uniform` returns the historical wet-heavy index set with no weights, and is the
default. What this sampler corrects is real but second-order:

!!! note "The temporal mismatch is not what dominates selection error"
    Scored against the exact full-record curve, the unweighted subsample's RMSE is +201% too
    high — yet that error is nearly a constant multiplier, so it cancels in the *ranking* and
    costs only ~0.1% regret. Correcting the level trades that harmless bias for real variance
    at a fixed hour budget and costs ~0.6%.

    What dominates is a different mismatch entirely, and it is spatial. Candidates are scored by
    leave-one-out at *training* stations, where the held-out station sits inside its own
    neighbourhood, but the benchmark reports RMSE at *fold* stations 20+ km from any training
    gauge. On the joint-covariate path the held-out RMSE curve falls monotonically from bw=8 to
    bw=160 while both tuning curves are U-shaped with a minimum at 12-30: over most of the grid
    the two criteria are anti-correlated, and the tuner's pick costs +19% to +101%. This is why
    widening the grid (`--local-grid`, floor 30 -> 8) made `residual_gwr`/`mixed_gwr` worse — the
    old floor was accidentally shielding the tuner from its own preference. Fixing that needs an
    inner *spatial* split of the training fold, not a reweighted hour sample.
"""
function _tuning_time_sample(
    y_obs::Matrix{Float64}, maximum_times::Int, weighting::Symbol=:stratified,
)
    total = size(y_obs, 2)
    (maximum_times <= 0 || total <= maximum_times) && return collect(1:total), nothing
    wetness = [let values = filter(isfinite, @view(y_obs[:, time]))
        isempty(values) ? -Inf : mean(values)
    end for time in 1:total]
    if weighting === :uniform
        wet_count = div(maximum_times, 2)
        wet = partialsortperm(wetness, 1:wet_count; rev=true)
        spaced = unique(round.(Int, range(1, total; length=maximum_times - wet_count)))
        return sort(unique(vcat(wet, spaced))), nothing
    end
    wet_count = max(1, div(maximum_times, 8))
    wet = partialsortperm(wetness, 1:wet_count; rev=true)
    # `rest` excludes the certainty stratum, so the two strata are disjoint and the sample size
    # is exactly `maximum_times` (the historical union could collapse to fewer hours).
    rest = setdiff(1:total, wet)
    dry_count = min(maximum_times - wet_count, length(rest))
    # Sample the remainder systematically along the wetness ordering rather than the calendar
    # ordering. Inclusion probability is uniform either way, so the weights stay valid, but
    # spanning the wetness distribution stops the estimate from swinging on how many wet hours
    # that fell outside the certainty stratum happen to be picked up at ~67x weight. Measured
    # over 40 synthetic replicates this cuts the mean error against the full-sample RMSE from
    # 6.3% to 2.2%, and the worst case from 19.8% to 2.2%.
    frame = rest[sortperm(wetness[rest])]
    dry = frame[unique(round.(Int, range(1, length(frame); length=dry_count)))]
    indices = vcat(wet, dry)
    weights = vcat(ones(Float64, length(wet)), fill(length(rest) / length(dry), length(dry)))
    order = sortperm(indices)
    return indices[order], weights[order]
end

_tuning_time_indices(y_obs::Matrix{Float64}, maximum_times::Int) =
    first(_tuning_time_sample(y_obs, maximum_times, :uniform))

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

function select_interpolation_parameter!(
    scan_rows::Vector{NamedTuple}, cfg::InterpolationBenchmarkConfig,
    method::String, mode::String, scheme::Symbol, product::String, fold::Int,
    train_lonlat::Matrix{Float64}, y_obs::Matrix{Float64}, y_sat::Matrix{Float64};
    dem_context=nothing, joint_context=nothing, hurdle_context=nothing,
    selection_groups=nothing, joint_selection_contexts=nothing, repeat_seed::Int=cfg.seed,
)
    (mode, method) in BENCHMARK_RUNS ||
        throw(ArgumentError("unsupported benchmark method/mode pair: $method/$mode"))
    if dem_context !== nothing && mode == "residual" && method in ("gwr", "mixed_gwr", "mgwr")
        return select_dem_parameter!(
            scan_rows, cfg, method, mode, scheme, product, fold, dem_context,
        )
    end
    # One inner split per fold, shared by every method, so candidates from different methods are
    # scored against the same held-out stations. The caller normally supplies it (it also needs
    # it to build the joint inner contexts); computing it here keeps direct callers working.
    if selection_groups === nothing && cfg.tuning_geometry === :inner_spatial
        selection_groups = selection_folds(
            cfg, scheme, string.(1:size(train_lonlat, 1)), train_lonlat, fold, repeat_seed,
        )
    end
    if joint_context !== nothing && mode == "residual" && method in ("gwr", "mixed_gwr", "mgwr")
        return select_joint_parameter!(
            scan_rows, cfg, method, mode, scheme, product, fold,
            joint_context, y_obs, y_sat; joint_selection_contexts,
        )
    end
    # `hurdle_gwr` is the one method that needs the timestamps themselves, so the tuning subset
    # has to be applied to them as well as to the data matrices.
    tuning_times = hurdle_context === nothing ? nothing : hurdle_context.times
    time_weights = nothing
    if cfg.tuning_max_times > 0 && size(y_obs, 2) > cfg.tuning_max_times
        time_idx, time_weights = _tuning_time_sample(
            y_obs, cfg.tuning_max_times, cfg.tuning_time_weighting,
        )
        y_obs = y_obs[:, time_idx]
        y_sat = y_sat[:, time_idx]
        tuning_times = tuning_times === nothing ? nothing : tuning_times[time_idx]
    end
    target_values = mode == "direct" ? y_obs : y_obs .- y_sat
    n_train, n_time = size(target_values)
    # `:loocv` keeps the historical criterion — predict at every training station with only that
    # station excluded. `:inner_spatial` predicts out-of-fold onto whole spatial groups instead,
    # which is the geometry the benchmark reports. Both produce an `n_train × n_time` matrix, so
    # everything downstream is identical.
    scored(predict_pair, predict_loocv) = selection_groups === nothing ?
        predict_loocv() : _selection_oof(selection_groups, n_train, n_time, predict_pair)
    first_row = length(scan_rows) + 1
    if method in ("idw", "adw")
        predictor = method == "idw" ? idw_predict : adw_predict
        for power in cfg.idw_powers, neighbors in cfg.neighbor_candidates
            try
                interpolated = scored(
                    (tr, va) -> predictor(
                        train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                        power=power, neighbors=neighbors,
                    ),
                    () -> predictor(
                        train_lonlat, target_values, train_lonlat;
                        power=power, neighbors=neighbors, exclude_self=true,
                    ),
                )
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(
                    y_obs, y_sat, prediction;
                    require_satellite=mode == "residual", time_weights,
                )
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, power, neighbors, metrics...,
                ))
            catch e
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, power, neighbors,
                    status="failed", error=sprint(showerror, e),
                ))
            end
        end
    elseif method == "tps"
        for smooth in cfg.tps_smooth_candidates
            try
                interpolated = scored(
                    (tr, va) -> tps_predict(
                        train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                        smooth=smooth,
                    ),
                    () -> tps_loo_predict(train_lonlat, target_values; smooth=smooth),
                )
                prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(
                    y_obs, y_sat, prediction;
                    require_satellite=mode == "residual", time_weights,
                )
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, smooth, metrics...,
                ))
            catch e
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, smooth,
                    status="failed", error=sprint(showerror, e),
                ))
            end
        end
    elseif method in ("gwr", "gwr_const")
        predictor = method == "gwr" ? _gwr_predict : _gwr_const_predict
        for kernel in cfg.mger.kernels
            for (adaptive, candidates) in ((true, cfg.mger.bw_adaptive), (false, cfg.mger.bw_fixed_km))
                for bw in candidates
                    try
                        interpolated = scored(
                            (tr, va) -> predictor(
                                train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :];
                                kernel, adaptive, bw,
                            ),
                            () -> predictor(
                                train_lonlat, target_values, train_lonlat;
                                kernel, adaptive, bw, exclude_self=true,
                            ),
                        )
                        prediction = mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat .+ interpolated, 0.0)
                        metrics = _candidate_metrics(
                            y_obs, y_sat, prediction;
                            require_satellite=mode == "residual", time_weights,
                        )
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw, metrics...,
                        ))
                    catch e
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw,
                            status="failed", error=sprint(showerror, e),
                        ))
                    end
                end
            end
        end
    elseif method == "hurdle_gwr"
        mode == "direct" || throw(ArgumentError("hurdle_gwr only supports direct mode"))
        hurdle_context === nothing &&
            throw(ArgumentError("hurdle_gwr requires a hurdle_context carrying the fold times"))
        for kernel in cfg.mger.kernels
            for (adaptive, candidates) in ((true, cfg.mger.bw_adaptive), (false, cfg.mger.bw_fixed_km))
                for bw in candidates
                    try
                        hurdle_config = _hurdle_config(hurdle_context, kernel, adaptive, bw)
                        corrected = scored(
                            (tr, va) -> _hurdle_predict(
                                train_lonlat[tr, :], y_obs[tr, :], y_sat[tr, :],
                                train_lonlat[va, :], y_sat[va, :], something(tuning_times);
                                config=hurdle_config,
                            ).corrected,
                            () -> _hurdle_predict(
                                train_lonlat, y_obs, y_sat, train_lonlat, y_sat,
                                something(tuning_times);
                                config=hurdle_config, exclude_self=true,
                            ).corrected,
                        )
                        metrics = _candidate_metrics(y_obs, y_sat, corrected; time_weights)
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw, metrics...,
                        ))
                    catch e
                        push!(scan_rows, _scan_row(;
                            scheme, product, fold, mode, method, kernel, adaptive, bw,
                            status="failed", error=sprint(showerror, e),
                        ))
                    end
                end
            end
        end
    elseif method == "mixed_gwr"
        mode == "residual" || throw(ArgumentError("mixed_gwr only supports residual mode"))
        for bw in cfg.mger.bw_adaptive
            try
                interpolated = scored(
                    (tr, va) -> _mixed_gwr_predict(
                        train_lonlat[tr, :], target_values[tr, :], train_lonlat[va, :]; bw,
                    ),
                    () -> _mixed_gwr_predict(
                        train_lonlat, target_values, train_lonlat; bw, exclude_self=true,
                    ),
                )
                prediction = max.(y_sat .+ interpolated, 0.0)
                metrics = _candidate_metrics(y_obs, y_sat, prediction; time_weights)
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, kernel=BISQUARE,
                    adaptive=true, bw, metrics...,
                ))
            catch e
                push!(scan_rows, _scan_row(;
                    scheme, product, fold, mode, method, kernel=BISQUARE,
                    adaptive=true, bw, status="failed", error=sprint(showerror, e),
                ))
            end
        end
    elseif method == "mgwr"
        mode == "residual" || throw(ArgumentError("mgwr only supports residual mode"))
        return select_mgwr_bandwidths!(
            scan_rows, cfg, scheme, product, fold, train_lonlat,
            target_values, y_obs, y_sat; time_weights, selection_groups,
        )
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return _select_candidate!(scan_rows, first_row, cfg.min_tuning_coverage)
end

function predict_selected(
    selected, method::String, mode::String,
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64}, y_sat_target::Matrix{Float64};
    dem_context=nothing, joint_context=nothing, hurdle_context=nothing,
)
    (mode, method) in BENCHMARK_RUNS ||
        throw(ArgumentError("unsupported benchmark method/mode pair: $method/$mode"))
    values = mode == "direct" ? y_obs_train : y_obs_train .- y_sat_train
    is_dem_model = dem_context !== nothing && hasproperty(selected, :dem_model) && selected.dem_model
    is_joint_model = joint_context !== nothing &&
        hasproperty(selected, :joint_model) && selected.joint_model
    interpolated = if is_joint_model && mode == "residual"
        prediction, converged = dynamic_covariate_predict(
            joint_context, values, selected.joint_method, selected.bandwidths,
        )
        any(.!converged) && @warn(
            "joint dynamic model did not converge for some hours",
            product=joint_context.product, method=selected.joint_method,
            failed_hours=count(.!converged),
        )
        # Folded in here rather than at the shared `y_sat_target .+ interpolated` tail below,
        # which serves every other method and must keep its historical behaviour.
        shrink = hasproperty(selected, :shrink) ? selected.shrink : 1.0
        isfinite(shrink) ? prediction .* shrink : prediction
    elseif is_dem_model && method == "gwr" && mode == "residual"
        designs = dem_context.all_local
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, zeros(Float64, size(train_lonlat, 1), 0), values,
            train_lonlat, designs.mixed_local_target,
            zeros(Float64, size(target_lonlat, 1), 0), target_lonlat,
            Int(round(selected.bw)); ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mixed_gwr"
        designs = dem_context.mixed
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, designs.global_train, values, train_lonlat,
            designs.mixed_local_target, designs.global_target, target_lonlat,
            Int(round(selected.bw)); ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mgwr"
        designs = dem_context.mixed
        prediction, _ = multiscale_gwr_predict(
            designs.multiscale_train, designs.global_train, values, train_lonlat,
            designs.multiscale_target, designs.global_target, target_lonlat,
            selected.bandwidths; ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif method == "idw"
        selected_neighbors = ismissing(selected.neighbors) ? nothing :
            (selected.neighbors == 0 ? nothing : Int(selected.neighbors))
        idw_predict(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "adw"
        selected_neighbors = ismissing(selected.neighbors) ? nothing :
            (selected.neighbors == 0 ? nothing : Int(selected.neighbors))
        adw_predict(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "tps"
        tps_predict(train_lonlat, values, target_lonlat; smooth=selected.smooth)
    elseif method == "gwr"
        _gwr_predict(train_lonlat, values, target_lonlat;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw)
    elseif method == "gwr_const"
        _gwr_const_predict(train_lonlat, values, target_lonlat;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw)
    elseif method == "hurdle_gwr"
        hurdle_context === nothing &&
            throw(ArgumentError("hurdle_gwr requires a hurdle_context carrying the fold times"))
        prediction = _hurdle_predict(
            train_lonlat, y_obs_train, y_sat_train, target_lonlat, y_sat_target,
            hurdle_context.times;
            config=_hurdle_config(
                hurdle_context, selected.kernel, selected.adaptive, selected.bw,
            ),
        )
        push!(hurdle_context.diagnostics, _hurdle_diagnostic_row(
            prediction, hurdle_context;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw,
        ))
        prediction.corrected
    elseif method == "mixed_gwr"
        _mixed_gwr_predict(train_lonlat, values, target_lonlat; bw=selected.bw)
    elseif method == "mgwr"
        _mgwr_predict(
            train_lonlat, values, target_lonlat; bandwidths=selected.bandwidths,
        )
    else
        throw(ArgumentError("unknown method: $method"))
    end
    return mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat_target .+ interpolated, 0.0)
end

function _empty_metric_row(; scheme, product, method, fold=missing, repeat=1, seed=0,
    group="overall", level="all",
    threshold=NaN, n=0, coverage=0.0, RMSE=NaN, MAE=NaN, Bias=NaN, r=NaN,
    POD=NaN, FAR=NaN, CSI=NaN)
    return (;
        scheme=String(scheme), product=String(product), method=String(method), fold,
        repeat=Int(repeat), seed=Int(seed),
        group=String(group), level=String(level), threshold=Float64(threshold), n=Int(n),
        coverage=Float64(coverage), RMSE=Float64(RMSE), MAE=Float64(MAE),
        Bias=Float64(Bias), r=Float64(r), POD=Float64(POD), FAR=Float64(FAR), CSI=Float64(CSI),
    )
end

function _continuous_row(y_obs, prediction, mask; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs...)
    metric = metric_continuous(y_obs, prediction; mask=mask)
    return _empty_metric_row(;
        kwargs..., n=metric.n, coverage=metric.n / length(mask), RMSE=metric.RMSE,
        MAE=metric.MAE, Bias=metric.Bias, r=metric.r,
    )
end

function _event_row(y_obs, prediction, mask, threshold; kwargs...)
    n = count(mask)
    n == 0 && return _empty_metric_row(; kwargs..., threshold)
    metric = metric_event(y_obs, prediction; mask=mask, thr=threshold)
    return _empty_metric_row(;
        kwargs..., threshold, n, coverage=n / length(mask),
        POD=metric.POD, FAR=metric.FAR, CSI=metric.CSI,
    )
end

function append_stratified_metrics!(
    rows::Vector{NamedTuple}, scheme::String, product::String, method::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, prediction::Matrix{Float64},
    common_mask::BitMatrix, nearest_distance::Vector{Float64}, thresholds::Vector{Float64};
    fold=missing, repeat::Int=1, seed::Int=0,
)
    base = (; scheme, product, method, fold, repeat, seed)
    # `common_mask` is pinned to MASK_METHODS, so a method outside that set can still be NaN
    # inside it and would otherwise score NaN rather than "worse". Drop its own gaps and let
    # the `coverage` column report the shortfall instead. This is a no-op for the MASK_METHODS
    # themselves, whose NaNs already defined the mask.
    scored_mask = common_mask .& .!isnan.(prediction)
    push!(rows, _continuous_row(y_obs, prediction, scored_mask; base..., group="overall", level="all"))

    rain_groups = (
        ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
        ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
    )
    for (name, lower, upper) in rain_groups
        stratum = scored_mask .& (y_obs .>= lower) .& (y_obs .< upper)
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="rain_intensity", level=name))
    end

    for year_value in sort(unique(year.(times)))
        time_mask = reshape(year.(times) .== year_value, 1, :)
        stratum = scored_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="year", level=string(year_value)))
    end
    for month_value in sort(unique(month.(times)))
        time_mask = reshape(month.(times) .== month_value, 1, :)
        stratum = scored_mask .& time_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="month", level=lpad(month_value, 2, '0')))
    end

    distance_groups = (
        ("0_20", -Inf, 20.0), ("20_50", 20.0, 50.0),
        ("50_100", 50.0, 100.0), ("100_plus", 100.0, Inf),
    )
    for (name, lower, upper) in distance_groups
        station_mask = reshape((nearest_distance .>= lower) .& (nearest_distance .< upper), :, 1)
        stratum = scored_mask .& station_mask
        push!(rows, _continuous_row(y_obs, prediction, stratum;
            base..., group="nearest_train_km", level=name))
    end

    for threshold in thresholds
        push!(rows, _event_row(y_obs, prediction, scored_mask, threshold;
            base..., group="event_threshold", level=string(threshold)))
    end
    return rows
end

function _common_method_mask(y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}})
    mask = .!isnan.(y_obs)
    for method in MASK_METHODS
        mask .&= .!isnan.(predictions[method])
    end
    return BitMatrix(mask)
end

function _daily_bootstrap_delta(
    rng::AbstractRNG, times::Vector{DateTime}, y_obs::Matrix{Float64},
    baseline::Matrix{Float64}, gwr::Matrix{Float64}, mask::BitMatrix, reps::Int,
)
    days = Date.(times)
    unique_days = sort(unique(days))
    daily_n = zeros(Int, length(unique_days))
    daily_sse_baseline = zeros(Float64, length(unique_days))
    daily_sse_gwr = zeros(Float64, length(unique_days))
    for (day_index, day) in enumerate(unique_days)
        time_idx = findall(==(day), days)
        day_mask = mask[:, time_idx]
        daily_n[day_index] = count(day_mask)
        daily_sse_baseline[day_index] = sum(abs2, (baseline[:, time_idx] .- y_obs[:, time_idx])[day_mask])
        daily_sse_gwr[day_index] = sum(abs2, (gwr[:, time_idx] .- y_obs[:, time_idx])[day_mask])
    end
    keep = findall(>(0), daily_n)
    isempty(keep) && return Float64[]
    deltas = Vector{Float64}(undef, reps)
    for rep in 1:reps
        sampled = rand(rng, keep, length(keep))
        n = sum(daily_n[sampled])
        rmse_baseline = sqrt(sum(daily_sse_baseline[sampled]) / n)
        rmse_gwr = sqrt(sum(daily_sse_gwr[sampled]) / n)
        deltas[rep] = rmse_baseline - rmse_gwr
    end
    return deltas
end

function _holm_adjust(pvalues::AbstractVector)
    pvalues = Float64.(pvalues)
    m = length(pvalues)
    order = sortperm(pvalues)
    adjusted = fill(NaN, m)
    running = 0.0
    for (rank_index, original_index) in enumerate(order)
        value = min(1.0, (m - rank_index + 1) * pvalues[original_index])
        running = max(running, value)
        adjusted[original_index] = running
    end
    return adjusted
end

"""
Method under test when the caller does not name one. Kept as the default everywhere so the
benchmark's own `paired_comparisons.csv` / `claim_assessment.csv` are unchanged by the
method-aware parameters added for `scripts/run_claim_reassessment.jl`.
"""
const DEFAULT_CLAIM_METHOD = "residual_gwr"

function paired_bootstrap_rows(
    cfg::InterpolationBenchmarkConfig, scheme::String, product::String,
    times::Vector{DateTime}, y_obs::Matrix{Float64}, predictions::Dict{String,Matrix{Float64}},
    common_mask::BitMatrix; repeat::Int=1, seed::Int=cfg.seed,
    method::String=DEFAULT_CLAIM_METHOD, pairwise_mask::Bool=false,
)
    rows = NamedTuple[]
    cfg.bootstrap_reps == 0 && return rows
    treatment = predictions[method]
    # A `method` column only appears once the caller asks for something other than the historical
    # single-method behaviour: the row already carries `baseline` but nothing naming the treatment,
    # so rows for two methods appended to one table would be indistinguishable.
    tagged = method != DEFAULT_CLAIM_METHOD || pairwise_mask
    strata = (
        ("overall", common_mask),
        ("moderate", common_mask .& (y_obs .>= 2.5) .& (y_obs .< 8.0)),
        ("heavy", common_mask .& (y_obs .>= 8.0)),
    )
    for (stratum, mask) in strata
        local_rows = NamedTuple[]
        for (baseline_index, baseline_method) in enumerate(TRADITIONAL_METHODS)
            baseline = predictions[baseline_method]
            # `common_mask` is pinned to MASK_METHODS, so a method outside that set can still be
            # NaN inside it and would score NaN rather than "worse". Restrict each test to the
            # cells both sides actually predicted.
            test_mask = pairwise_mask ?
                BitMatrix(mask .& .!isnan.(baseline) .& .!isnan.(treatment)) : BitMatrix(mask)
            deltas = _daily_bootstrap_delta(
                MersenneTwister(seed + 1000 * baseline_index + sum(codeunits(product))),
                times, y_obs, baseline, treatment, test_mask, cfg.bootstrap_reps,
            )
            isempty(deltas) && continue
            observed_baseline = metric_continuous(y_obs, baseline; mask=test_mask).RMSE
            observed_gwr = metric_continuous(y_obs, treatment; mask=test_mask).RMSE
            delta = observed_baseline - observed_gwr
            pvalue = min(1.0, 2 * min(mean(deltas .<= 0), mean(deltas .>= 0)))
            row = (;
                scheme, product, stratum, baseline=baseline_method, repeat, seed,
                n=count(test_mask), n_day=length(unique(Date.(times))), reps=cfg.bootstrap_reps,
                RMSE_baseline=observed_baseline, RMSE_gwr=observed_gwr,
                delta_RMSE=delta, relative_improvement=delta / observed_baseline,
                ci_low=quantile(deltas, 0.025), ci_high=quantile(deltas, 0.975),
                pvalue, pvalue_holm=NaN,
            )
            push!(local_rows, tagged ? merge((; method), row) : row)
        end
        adjusted = _holm_adjust([row.pvalue for row in local_rows])
        for (index, row) in enumerate(local_rows)
            push!(rows, merge(row, (; pvalue_holm=adjusted[index])))
        end
    end
    return rows
end

const REPEAT_SUMMARY_METRICS = (:RMSE, :MAE, :Bias, :r, :POD, :FAR, :CSI)

_finite(values) = filter(isfinite, collect(skipmissing(values)))

"""
Spread of each pooled metric across repeated cross-validation partitions.

One row per (scheme, product, method, group, level) carrying mean/std/min/max over repeats, so a
reported value is read as a distribution rather than a single split's point estimate. `mean_n` and
`mean_coverage` are kept because `_common_method_mask` requires every method to be finite, so the
evaluated sample size itself shifts between partitions.
"""
function summarize_repeats(metrics::DataFrame)
    nrow(metrics) == 0 && return DataFrame()
    pooled = filter(row -> ismissing(row.fold), metrics)
    nrow(pooled) == 0 && return DataFrame()
    rows = NamedTuple[]
    for subset in groupby(pooled, [:scheme, :product, :method, :group, :level])
        entry = (;
            scheme=first(subset.scheme), product=first(subset.product),
            method=first(subset.method), group=first(subset.group), level=first(subset.level),
            n_repeats=nrow(subset), mean_n=mean(subset.n), mean_coverage=mean(subset.coverage),
        )
        for name in REPEAT_SUMMARY_METRICS
            values = _finite(subset[!, name])
            entry = merge(entry, (;
                Symbol(name, :_mean) => isempty(values) ? NaN : mean(values),
                Symbol(name, :_std) => length(values) < 2 ? NaN : std(values),
                Symbol(name, :_min) => isempty(values) ? NaN : minimum(values),
                Symbol(name, :_max) => isempty(values) ? NaN : maximum(values),
            ))
        end
        push!(rows, entry)
    end
    return sort!(DataFrame(rows), [:scheme, :product, :group, :level, :method])
end

"""
How often each method wins, across repeated cross-validation partitions.

Methods are ranked by pooled RMSE within every repeat; a `win_count` well below `n_repeats` means
the ranking is an artifact of the particular split rather than a stable result.
"""
function method_rank_stability(metrics::DataFrame)
    nrow(metrics) == 0 && return DataFrame()
    pooled = filter(row -> ismissing(row.fold) && isfinite(row.RMSE), metrics)
    nrow(pooled) == 0 && return DataFrame()
    rows = NamedTuple[]
    for cell in groupby(pooled, [:scheme, :product, :group, :level])
        ranks = Dict{String,Vector{Int}}()
        for repeat_value in sort(unique(cell.repeat))
            ordered = sort(filter(:repeat => ==(repeat_value), DataFrame(cell)), :RMSE)
            for (position, row) in enumerate(eachrow(ordered))
                push!(get!(ranks, row.method, Int[]), position)
            end
        end
        for method in sort(collect(keys(ranks)))
            positions = ranks[method]
            push!(rows, (;
                scheme=first(cell.scheme), product=first(cell.product),
                group=first(cell.group), level=first(cell.level), method,
                n_repeats=length(positions), mean_rank=mean(positions),
                best_rank=minimum(positions), worst_rank=maximum(positions),
                win_count=count(==(1), positions),
            ))
        end
    end
    return sort!(DataFrame(rows), [:scheme, :product, :group, :level, :mean_rank])
end

function _write_split(path::String, ids::Vector{String}, folds::Vector{Vector{String}}, scheme::Symbol)
    fold_map = Dict(id => fold for (fold, fold_ids) in enumerate(folds) for id in fold_ids)
    CSV.write(path, DataFrame(
        station_id=ids, fold=[fold_map[id] for id in ids], scheme=fill(string(scheme), length(ids)),
    ))
end

"""
Pooled balanced-spatial metric for the first repeat. The claim assessment is defined on a single
partition; across-partition evidence lives in `metrics_repeat_summary.csv` instead.
"""
function _one_metric(metrics::DataFrame, product::String, method::String, group::String, level::String)
    rows = filter(row ->
        row.scheme == "balanced_spatial" && row.product == product && row.method == method &&
        ismissing(row.fold) && row.repeat == 1 && row.group == group && row.level == level,
        metrics,
    )
    nrow(rows) == 1 || throw(ArgumentError(
        "expected one pooled metric for $product/$method/$group/$level, got $(nrow(rows))",
    ))
    return rows[1, :]
end

"""
Assess whether the residual-GWR correction is supported for each product.

Every comparison against the traditional baselines (`idw`/`adw`/`tps`) is required to hold
against **all three**, not just whichever looks friendliest on this pooled sample. The previous
version picked a single baseline per stratum via `argmin`/`argmax` on the same data it reported
— an unprotected selection-on-test step — then tested or compared only against that one.
Requiring all three removes the selection entirely rather than correcting for it: the three
per-baseline `overall` comparisons are already Holm-corrected against each other by
`paired_bootstrap_rows`, and no further correction is needed to require all three to pass.
`best_traditional`/`RMSE_best_traditional`/`overall_relative_improvement` are kept as
descriptive-only fields (how much better than the strongest competitor) and do not gate
`product_supported`; the `*_win_count` fields show how many of the 3 baselines were actually
beaten in each stratum.

`method` selects which method is assessed (default `residual_gwr`, the historical behaviour).
`own_coverage` maps product to the method's *own* prediction coverage; when supplied, the
coverage gate uses it instead of `coverage` on the shared evaluation mask. The shared mask is
pinned to `MASK_METHODS`, so it sits at ~0.8 because direct `gwr` fails ~20% of cells — gating on
it asks every method to answer for a different method's failures, which no method can pass. Both
numbers are reported side by side. Naming either argument switches the output to a self-describing
schema (`method` + `RMSE_method`) so rows for different methods cannot be silently mixed.
"""
function assess_gwr_claim(
    metrics::DataFrame, bootstrap::DataFrame, products::Vector{String};
    method::String=DEFAULT_CLAIM_METHOD, own_coverage::Union{Nothing,AbstractDict}=nothing,
)
    rows = NamedTuple[]
    tagged = method != DEFAULT_CLAIM_METHOD || own_coverage !== nothing
    bootstrap_has_method = :method in propertynames(bootstrap)
    for product in products
        overall_traditional = [_one_metric(metrics, product, method, "overall", "all")
            for method in TRADITIONAL_METHODS]
        best_index = argmin([row.RMSE for row in overall_traditional])
        best_method = TRADITIONAL_METHODS[best_index]
        best_overall = overall_traditional[best_index]
        gwr_overall = _one_metric(metrics, product, method, "overall", "all")

        # Overall: GWR must be significantly better than every baseline (each row already
        # Holm-corrected against the other two by paired_bootstrap_rows), not just whichever
        # one was easiest to beat.
        overall_bootstrap = filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.stratum == "overall" && row.repeat == 1 &&
            (!bootstrap_has_method || row.method == method),
            bootstrap,
        )
        significant_per_baseline = Dict(String(row.baseline) =>
            row.ci_low > 0 && row.pvalue_holm < 0.05 for row in eachrow(overall_bootstrap))
        overall_win_count = count(values(significant_per_baseline))
        significant = length(significant_per_baseline) == length(TRADITIONAL_METHODS) &&
            overall_win_count == length(TRADITIONAL_METHODS)

        # Heavy: GWR must improve on every baseline by at least 5%.
        gwr_heavy = _one_metric(metrics, product, method, "rain_intensity", "heavy").RMSE
        heavy_improvements = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "rain_intensity", "heavy").RMSE
                isfinite(base) && base > 0 ? (base - gwr_heavy) / base : NaN
            end for baseline in TRADITIONAL_METHODS)
        heavy_win_count = count(>=(0.05), filter(isfinite, collect(values(heavy_improvements))))
        heavy_ok = heavy_win_count == length(TRADITIONAL_METHODS)
        # Worst case across baselines, matching what the gate above requires.
        heavy_relative_improvement = minimum(values(heavy_improvements))

        # Moderate: GWR is allowed to be marginally worse here (non-inferiority), but not by
        # more than 2% against any baseline.
        gwr_moderate = _one_metric(metrics, product, method, "rain_intensity", "moderate").RMSE
        moderate_degradations = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "rain_intensity", "moderate").RMSE
                isfinite(base) && base > 0 ? (gwr_moderate - base) / base : NaN
            end for baseline in TRADITIONAL_METHODS)
        moderate_win_count = count(<=(0.02), filter(isfinite, collect(values(moderate_degradations))))
        moderate_ok = moderate_win_count == length(TRADITIONAL_METHODS)
        moderate_relative_degradation = maximum(values(moderate_degradations))

        # Year: a "win" requires beating every baseline that year, not just one.
        # `String.` because a `metrics` table read back from CSV carries InlineString levels,
        # which `_one_metric`'s `::String` signature would reject.
        year_levels = String.(unique(filter(row ->
            row.scheme == "balanced_spatial" && row.product == product &&
            row.method == method && ismissing(row.fold) && row.group == "year",
            metrics,
        ).level))
        year_wins = 0
        for year_level in year_levels
            gwr_year = _one_metric(metrics, product, method, "year", year_level).RMSE
            year_wins += all(baseline ->
                gwr_year < _one_metric(metrics, product, baseline, "year", year_level).RMSE,
                TRADITIONAL_METHODS)
        end
        majority_years = isempty(year_levels) ? false : year_wins > length(year_levels) / 2

        # Event threshold: must not degrade CSI/FAR against any baseline.
        gwr_event = _one_metric(metrics, product, method, "event_threshold", "0.1")
        event_checks = Dict(baseline => let
                base = _one_metric(metrics, product, baseline, "event_threshold", "0.1")
                gwr_event.CSI >= base.CSI - 0.02 && gwr_event.FAR <= base.FAR + 0.02
            end for baseline in TRADITIONAL_METHODS)
        event_win_count = count(values(event_checks))
        event_not_degraded = event_win_count == length(TRADITIONAL_METHODS)

        gated_coverage = own_coverage === nothing ?
            gwr_overall.coverage : Float64(own_coverage[product])
        coverage_acceptable = gated_coverage >= 0.95
        direction_improved = gwr_overall.RMSE < best_overall.RMSE
        product_supported = significant && heavy_ok && moderate_ok && majority_years &&
            event_not_degraded && coverage_acceptable
        head = (;
            product, best_traditional=best_method,
            RMSE_best_traditional=best_overall.RMSE,
        )
        tail = (;
            overall_relative_improvement=(best_overall.RMSE - gwr_overall.RMSE) / best_overall.RMSE,
            paired_significant=significant, overall_win_count=overall_win_count,
            heavy_relative_improvement=heavy_relative_improvement, heavy_win_count=heavy_win_count,
            moderate_relative_degradation=moderate_relative_degradation,
            moderate_win_count=moderate_win_count,
            year_win_count=year_wins, year_count=length(year_levels), majority_years,
            event_not_degraded, event_win_count, common_coverage=gwr_overall.coverage,
            coverage_acceptable, direction_improved, product_supported,
        )
        push!(rows, if tagged
            merge((; method), head,
                (; RMSE_method=gwr_overall.RMSE, own_coverage=gated_coverage), tail)
        else
            merge(head, (; RMSE_residual_gwr=gwr_overall.RMSE), tail)
        end)
    end
    supported_products = count(row -> row.product_supported, rows)
    improved_products = count(row -> row.direction_improved, rows)
    overall_supported = supported_products >= 2 && improved_products >= 2
    return DataFrame([merge(row, (;
        supported_product_count=supported_products,
        improved_product_count=improved_products,
        overall_claim_supported=overall_supported,
    )) for row in rows])
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

function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))
    terrain = _dem_enabled(cfg) ? load_aligned_terrain(something(cfg.terrain_path), ids) : nothing
    dem_store = _empty_dem_store()
    common_times = product_data[first(products)].times
    joint_inputs = _joint_enabled(cfg) ? load_joint_benchmark_inputs(
        something(cfg.joint_covariates), products, ids, common_times,
    ) : nothing
    nested_joint = _joint_enabled(cfg) && cfg.joint_selection !== nothing
    joint_store = _empty_joint_store()
    joint_scaling_tables = DataFrame[]
    joint_qc_tables = DataFrame[]
    if joint_inputs !== nothing
        if nested_joint
            CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
                source_path="", sha256="", selection_mode="nested_per_fold", confirmatory=true,
            )]))
        else
            CSV.write(
                joinpath(cfg.mger.outdir, "joint_variable_spec_used.csv"),
                joint_inputs.specification.table,
            )
            CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
                source_path=abspath(something(cfg.joint_covariates).spec_path),
                sha256=joint_inputs.spec_sha256,
                selection_mode="fixed_full_data",
                confirmatory=false,
            )]))
        end
        CSV.write(joinpath(cfg.mger.outdir, "joint_era5_input_qc.csv"), joint_inputs.era5.qc)
        if joint_inputs.ndvi !== nothing
            CSV.write(joinpath(cfg.mger.outdir, "joint_ndvi_alignment_qc.csv"),
                joint_inputs.ndvi.aligned.qc)
        end
    end

    all_metric_rows = NamedTuple[]
    all_scan_rows = NamedTuple[]
    all_bootstrap_rows = NamedTuple[]
    run_status_rows = NamedTuple[]
    hurdle_rows = NamedTuple[]
    # Cells where the fold was too small for an inner selection split and fell back to
    # leave-one-out. Reported in `benchmark_scope.csv` so the fallback is never silent.
    selection_fallback_cells = String[]

    if _dem_enabled(cfg)
        dem = something(cfg.dem)
        for product in products
            data = product_data[product]
            selection = screen_dem_subset(
                terrain, lonlat, data.times, data.Y_obs, data.Y_sat, dem;
                scheme="full_data", product, fold=0, phase="full_data",
                seed=cfg.seed + 100_000 + Int(sum(codeunits(product))),
            )
            _store_dem_selection!(dem_store, selection)
        end
    end

    seeds = benchmark_seeds(cfg)
    # Repeated cross-validation: one independent fold partition per seed. The body is left at its
    # original indentation to keep the diff reviewable.
    for (repeat_index, repeat_seed) in enumerate(seeds)
    repeat_root = _repeat_dir(cfg.mger.outdir, repeat_index, length(seeds))
    for scheme_symbol in cfg.cv_schemes
        scheme = string(scheme_symbol)
        scheme_dir = joinpath(repeat_root, scheme)
        mkpath(scheme_dir)
        folds = benchmark_folds(
            scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
            center_init=cfg.fold_center_init,
        )
        _write_split(joinpath(scheme_dir, "split_common.csv"), ids, folds, scheme_symbol)
        id_map = Dict(id => index for (index, id) in enumerate(ids))

        for product in products
            data = product_data[product]
            y_obs = data.Y_obs
            y_sat = data.Y_sat
            predictions = Dict(method => fill(NaN, size(y_obs)) for method in BENCHMARK_METHODS)
            predictions["raw"] .= y_sat
            nearest_train_distance = fill(NaN, length(ids))

            for fold in 1:cfg.k
                val_ids = folds[fold]
                train_ids = reduce(vcat, (folds[index] for index in 1:cfg.k if index != fold))
                train_idx = [id_map[id] for id in train_ids]
                val_idx = [id_map[id] for id in val_ids]
                train_lonlat = Matrix{Float64}(lonlat[train_idx, :])
                val_lonlat = Matrix{Float64}(lonlat[val_idx, :])
                y_obs_train = Matrix{Float64}(y_obs[train_idx, :])
                y_sat_train = Matrix{Float64}(y_sat[train_idx, :])
                y_sat_val = Matrix{Float64}(y_sat[val_idx, :])
                distance_train_val = haversine_distance_matrix(train_lonlat, val_lonlat)
                nearest_train_distance[val_idx] = vec(minimum(distance_train_val, dims=1))
                null_predictions = _null_fold_predictions(y_obs_train, length(val_idx))
                for (null_method, null_prediction) in null_predictions
                    predictions[null_method][val_idx, :] = null_prediction
                end
                # One inner split for the whole fold, so every method's candidates are scored
                # against the same held-out stations, and so the joint contexts below agree with
                # what the spatial-only methods use.
                selection_groups = cfg.tuning_geometry === :inner_spatial ?
                    selection_folds(cfg, scheme_symbol, train_ids, train_lonlat, fold, repeat_seed) :
                    nothing
                if cfg.tuning_geometry === :inner_spatial && selection_groups === nothing
                    push!(selection_fallback_cells, "$scheme/$product/fold$fold")
                end

                dem_context = nothing
                if _dem_enabled(cfg)
                    dem = something(cfg.dem)
                    selection = screen_dem_subset(
                        terrain[train_idx, :], train_lonlat, data.times,
                        y_obs_train, y_sat_train, dem; scheme, product, fold, phase="cv",
                        seed=repeat_seed + 10_000 * fold + Int(sum(codeunits(scheme * product))),
                    )
                    _store_dem_selection!(dem_store, selection)
                    dem_context = build_dem_fold_context(
                        selection, terrain[train_idx, :], terrain[val_idx, :],
                        train_lonlat, val_lonlat, dem,
                    )
                end
                joint_context = nothing
                joint_selection_contexts = nothing
                if joint_inputs !== nothing
                    joint = something(cfg.joint_covariates)
                    role_map = if nested_joint
                        product_index = findfirst(==(product), products)
                        run_seed = repeat_seed + 1000 * product_index + 10 * fold
                        dem_seed = repeat_seed + 20260815 + Int(sum(codeunits(product))) +
                            10_000 + 10 * fold
                        joint_selection = _screen_joint_subset(
                            # `select_joint_covariates` expects population-wide matrices and
                            # slices them via `train_idx` itself (unlike the DEM path's
                            # pre-sliced `terrain[train_idx,:]` convention above).
                            product, y_obs, y_sat, joint_inputs.terrain,
                            joint_inputs.era5.values,
                            joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
                            train_idx, ids, data.times, lonlat, something(cfg.joint_selection);
                            scheme, fold, repeat=repeat_index, seed=repeat_seed, run_seed, dem_seed,
                        )
                        _store_joint_selection!(joint_store, joint_selection)
                        joint_selection.role_map
                    else
                        joint_inputs.specification.role_maps[product]
                    end
                    joint_context = build_joint_fold_context(
                        product, role_map,
                        train_idx, val_idx, lonlat, y_obs, y_sat,
                        joint_inputs.terrain, joint_inputs.era5.values,
                        joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
                        joint,
                    )
                    scaling = copy(joint_context.scaling)
                    insertcols!(scaling, 1,
                        :scheme => fill(scheme, nrow(scaling)),
                        :fold => fill(fold, nrow(scaling)))
                    push!(joint_scaling_tables, scaling)
                    quality = copy(joint_context.quality_control)
                    insertcols!(quality, 1,
                        :scheme => fill(scheme, nrow(quality)),
                        :fold => fill(fold, nrow(quality)))
                    push!(joint_qc_tables, quality)
                    # One context per inner selection group, so joint candidates can be scored by
                    # predicting onto held-out stations instead of leave-one-out. Same builder,
                    # same role map, inner indices — all scaling refits on the inner training set.
                    if selection_groups !== nothing
                        joint_selection_contexts = [(
                            target_positions=group,
                            context=build_joint_fold_context(
                                product, role_map,
                                train_idx[setdiff(1:length(train_idx), group)],
                                train_idx[group], lonlat, y_obs, y_sat,
                                joint_inputs.terrain, joint_inputs.era5.values,
                                joint_inputs.ndvi === nothing ? nothing : joint_inputs.ndvi.aligned,
                                joint,
                            ),
                        ) for group in selection_groups]
                    end
                end

                hurdle_context = build_hurdle_context(
                    data.times, cfg, hurdle_rows; scheme, product, fold,
                )

                fold_predictions = merge(
                    Dict{String,Matrix{Float64}}("raw" => y_sat_val), null_predictions,
                )
                scan_start = length(all_scan_rows) + 1
                for (mode, method) in BENCHMARK_RUNS
                    output_method = method in ("mixed_gwr", "mgwr") ? method :
                        (mode == "direct" ? method : "residual_$(method)")
                    try
                        selected = select_interpolation_parameter!(
                            all_scan_rows, cfg, method, mode, scheme_symbol, product, fold,
                            train_lonlat, y_obs_train, y_sat_train;
                            dem_context, joint_context, hurdle_context,
                            selection_groups, joint_selection_contexts, repeat_seed,
                        )
                        fold_predictions[output_method] = predict_selected(
                            selected, method, mode, train_lonlat, val_lonlat,
                            y_obs_train, y_sat_train, y_sat_val;
                            dem_context, joint_context, hurdle_context,
                        )
                        predictions[output_method][val_idx, :] = fold_predictions[output_method]
                        uses_dem = dem_context !== nothing && mode == "residual" &&
                            method in ("gwr", "mixed_gwr", "mgwr")
                        uses_joint = joint_context !== nothing && mode == "residual" &&
                            method in ("gwr", "mixed_gwr", "mgwr")
                        selected_roles = uses_joint ? join([
                            "$group=$(joint_context.roles[group])" for group in joint_context.variables
                        ], ";") : ""
                        effective_role_map = uses_joint ? joint_effective_roles(
                            joint_context, output_method,
                        ) : Dict{String,String}()
                        effective_roles = uses_joint ? join([
                            "$group=$(effective_role_map[group])" for group in joint_context.variables
                        ], ";") : ""
                        eligible = .!isnan.(y_obs[val_idx, :]) .& .!isnan.(y_sat_val)
                        coverage = count(eligible) == 0 ? 0.0 :
                            count(eligible .& .!isnan.(fold_predictions[output_method])) / count(eligible)
                        push!(run_status_rows, (;
                            scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            method=output_method,
                            status=coverage >= cfg.min_tuning_coverage ? "success" : "partial",
                            error=coverage >= cfg.min_tuning_coverage ? "" :
                                "prediction coverage below minimum",
                            prediction_coverage=coverage,
                            dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
                            dem_selection_status=uses_dem ? dem_context.selection_status :
                                "not_applicable",
                            dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
                            dem_roles=uses_dem ? dem_context.dem_roles : "",
                            covariate_selection_mode=uses_joint ?
                                (nested_joint ? "nested_per_fold" : "fixed_full_data") :
                                "not_applicable",
                            covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
                            covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
                            covariate_selected_roles=selected_roles,
                            covariate_effective_roles=effective_roles,
                            covariate_spec_sha256=uses_joint ?
                                something(joint_inputs.spec_sha256, "") : "",
                        ))
                    catch e
                        fold_predictions[output_method] = fill(NaN, length(val_idx), size(y_obs, 2))
                        uses_dem = dem_context !== nothing && mode == "residual" &&
                            method in ("gwr", "mixed_gwr", "mgwr")
                        uses_joint = joint_context !== nothing && mode == "residual" &&
                            method in ("gwr", "mixed_gwr", "mgwr")
                        selected_roles = uses_joint ? join([
                            "$group=$(joint_context.roles[group])" for group in joint_context.variables
                        ], ";") : ""
                        effective_role_map = uses_joint ? joint_effective_roles(
                            joint_context, output_method,
                        ) : Dict{String,String}()
                        effective_roles = uses_joint ? join([
                            "$group=$(effective_role_map[group])" for group in joint_context.variables
                        ], ";") : ""
                        push!(run_status_rows, (;
                            scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            method=output_method, status="failed",
                            error=sprint(showerror, e), prediction_coverage=0.0,
                            dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
                            dem_selection_status=uses_dem ? dem_context.selection_status :
                                "not_applicable",
                            dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
                            dem_roles=uses_dem ? dem_context.dem_roles : "",
                            covariate_selection_mode=uses_joint ?
                                (nested_joint ? "nested_per_fold" : "fixed_full_data") :
                                "not_applicable",
                            covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
                            covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
                            covariate_selected_roles=selected_roles,
                            covariate_effective_roles=effective_roles,
                            covariate_spec_sha256=uses_joint ?
                                something(joint_inputs.spec_sha256, "") : "",
                        ))
                    end
                end

                for index in scan_start:length(all_scan_rows)
                    all_scan_rows[index] = merge(
                        all_scan_rows[index], (; repeat=repeat_index, seed=repeat_seed),
                    )
                end

                fold_mask = _common_method_mask(Matrix{Float64}(y_obs[val_idx, :]), fold_predictions)
                if any(fold_mask)
                    for method in BENCHMARK_METHODS
                        append_stratified_metrics!(
                            all_metric_rows, scheme, product, method, data.times,
                            Matrix{Float64}(y_obs[val_idx, :]), fold_predictions[method], fold_mask,
                            nearest_train_distance[val_idx], cfg.event_thresholds; fold,
                            repeat=repeat_index, seed=repeat_seed,
                        )
                    end
                end
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
            for method in BENCHMARK_METHODS
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

    metrics = DataFrame(all_metric_rows)
    scans = DataFrame(all_scan_rows)
    bootstrap = DataFrame(all_bootstrap_rows)
    status = DataFrame(run_status_rows)
    claim = if :balanced_spatial in cfg.cv_schemes && cfg.bootstrap_reps > 0
        assess_gwr_claim(metrics, bootstrap, products)
    else
        DataFrame()
    end
    repeat_summary = summarize_repeats(metrics)
    rank_stability = method_rank_stability(metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_stratified.csv"), metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_folds.csv"), filter(:fold => (x -> !ismissing(x)), metrics))
    CSV.write(joinpath(cfg.mger.outdir, "metrics_pooled.csv"), filter(:fold => ismissing, metrics))
    nrow(repeat_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_repeat_summary.csv"), repeat_summary)
    nrow(rank_stability) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "method_rank_stability.csv"), rank_stability)
    CSV.write(joinpath(cfg.mger.outdir, "parameter_scan.csv"), scans)
    # How much of each hurdle prediction actually came from a local fit. Without this a
    # globally-degenerate model is indistinguishable from a local one in the metrics.
    isempty(hurdle_rows) ||
        CSV.write(joinpath(cfg.mger.outdir, "hurdle_diagnostics.csv"), DataFrame(hurdle_rows))
    nrow(bootstrap) > 0 && CSV.write(joinpath(cfg.mger.outdir, "paired_comparisons.csv"), bootstrap)
    CSV.write(joinpath(cfg.mger.outdir, "run_status.csv"), status)
    _dem_enabled(cfg) && _write_dem_outputs(cfg.mger.outdir, dem_store, scans, status)
    if joint_inputs !== nothing
        scaling = isempty(joint_scaling_tables) ? DataFrame() :
            vcat(joint_scaling_tables...; cols=:union)
        quality = isempty(joint_qc_tables) ? DataFrame() :
            vcat(joint_qc_tables...; cols=:union)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_scaling.csv"), scaling)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_quality_control.csv"), quality)
        bandwidths = filter([:mode, :method, :selected] =>
            (mode, method, selected) -> mode == "residual" &&
                method in ("gwr", "mixed_gwr", "mgwr") && selected, scans)
        CSV.write(joinpath(cfg.mger.outdir, "joint_bandwidths.csv"), bandwidths)
        CSV.write(joinpath(cfg.mger.outdir, "covariate_model_status.csv"), filter(
            :covariate_selection_mode => in(("fixed_full_data", "nested_per_fold")), status,
        ))
        nested_joint && _write_joint_selection_outputs(cfg.mger.outdir, joint_store)
    end
    if ncol(claim) > 0
        CSV.write(joinpath(cfg.mger.outdir, "claim_assessment.csv"), claim)
    end
    scope = DataFrame(
        key=[
            "repeated_cv_partitions", "repeated_cv_seeds", "fold_center_init",
            "validation_target", "training_signal", "temporal_holdout", "primary_cv",
            "secondary_cv", "common_evaluation_mask", "tuning_time_limit",
            "tuning_time_weighting", "tuning_geometry",
            "dem_variable_selection", "dem_role_assignment", "dem_leakage_control",
            "dem_empty_selection", "supported_claim", "unsupported_claim",
        ],
        value=[
            string(length(seeds)), join(seeds, ","), string(cfg.fold_center_init),
            "held-out stations at matched observation/satellite timestamps",
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus fold-selected DEM variables for residual_gwr, mixed_gwr, and mgwr",
            "false", "balanced two-dimensional spatial 5-fold", "random station 5-fold",
            "true across all eight methods", string(cfg.tuning_max_times),
            cfg.tuning_time_weighting === :stratified ?
                "stratified: wettest eighth taken with certainty, remaining hours systematically subsampled and inverse-probability weighted so the tuning RMSE estimates the reported pooled RMSE" :
                "uniform (default): unweighted RMSE over a wet-oversampled subsample, so the tuning RMSE runs ~3x the reported pooled RMSE it stands in for; the level error is close to a constant multiplier and cancels in the candidate ranking",
            cfg.tuning_geometry === :inner_spatial ?
                "inner_spatial (default): candidates scored by out-of-fold prediction onto an inner $(cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k)-group split of the training stations, built with the same splitter and scheme as the outer partition, so selection and reporting are the same estimand" *
                (isempty(selection_fallback_cells) ? "" :
                    "; FELL BACK to leave-one-out in $(length(selection_fallback_cells)) cell(s) too small to split: $(join(selection_fallback_cells, ", "))") :
                "loocv (legacy): candidates scored by leave-one-out at training stations, which measures interpolation next to a retained gauge while the reported metric measures extrapolation to a station 20+ km from any gauge",
            _dem_enabled(cfg) ? "training-fold Pearson/Spearman direction check, joint aspect F test, BH q<0.05, and grouped VIF<5" : "disabled",
            _dem_enabled(cfg) ? "training-fold GWR Monte Carlo spatial nonstationarity test with BH q<0.05" : "disabled",
            _dem_enabled(cfg) ? "all DEM screening, scaling, role tests, and bandwidth selection use training stations only" : "not applicable",
            _dem_enabled(cfg) ? "fall back to the original spatial-only residual design and record dem_variable_count=0" : "not applicable",
            "gauge-network-assisted interpolation/correction at stations not used for fitting",
            "temporal forecast or correction without concurrent gauge observations",
        ],
    )
    if joint_inputs !== nothing
        joint = something(cfg.joint_covariates)
        append!(scope, DataFrame(
            key=["mgwr_spatial_grouping", "backfit_relaxation", "backfit_max_iterations",
                "residual_shrinkage"],
            value=[
                joint.mgwr_spatial_grouping === :split ?
                    "split (default): the intercept, longitude and latitude each form their own single-column back-fitting group" :
                    joint.mgwr_spatial_grouping === :shared ?
                        "shared: intercept, longitude and latitude solved as one weighted least squares, so only the covariates carry separate bandwidths" :
                        "intercept_only: the coordinate columns are dropped; a locally varying intercept plus one group per local covariate",
                "$(joint.relaxation) (successive over-relaxation on the back-fitting sweeps; converges to the same fixed point for any admissible value, so this trades rate only — plain Gauss-Seidel at 1.0 fails to converge on most hours)",
                string(joint.max_iterations),
                length(cfg.residual_shrinkage_candidates) == 1 &&
                    only(cfg.residual_shrinkage_candidates) == 1.0 ?
                    "disabled: the residual correction is added back unshrunk, as GWR produces it" :
                    "enabled: a scale in $(minimum(cfg.residual_shrinkage_candidates))-$(maximum(cfg.residual_shrinkage_candidates)) is selected with the bandwidth on the inner spatial split and applied to the joint dynamic models' residual correction",
            ],
        ))
        training_row = findfirst(==("training_signal"), scope.key)
        scope.value[training_row] = nested_joint ?
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a per-fold, training-station-only product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr" :
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a fixed full-data product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr"
        append!(scope, DataFrame(
            key=[
                "covariate_selection_mode", "covariate_selection_leakage",
                "covariate_parameter_leakage", "covariate_temporal_fit",
                "inference_policy",
            ],
            value=nested_joint ? [
                "nested_per_fold",
                "variable and role selection re-run inside every training fold; validation stations never seen by selection",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "eligible for the same paired significance test and claim assessment as the legacy DEM path",
            ] : [
                "fixed_full_data",
                "non-nested variable selection used all 237 stations before spatial cross-validation",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "exploratory performance metrics only; no paired significance or claim assessment",
            ],
        ))
    end
    CSV.write(joinpath(cfg.mger.outdir, "benchmark_scope.csv"), scope)
    return (; metrics, scans, bootstrap, status, claim, repeat_summary, rank_stability)
end
