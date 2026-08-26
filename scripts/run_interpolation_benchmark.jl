#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using DataFrames, Dates

include(joinpath(ROOT, "src", "InterpolationBenchmark.jl"))

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")

"""Deterministic seed list for repeated cross-validation, so the sweep itself is reproducible."""
repeat_seeds(repeats::Int) = repeats <= 1 ? Int[] : [20260627 + 1000 * i for i in 0:(repeats - 1)]

function benchmark_config(
    mode::Symbol; with_random::Bool=false, legacy_dem::Bool=false, repeats::Int=1,
    nested_covariates::Bool=false, local_grid::Bool=false,
    stratified_tuning_weights::Bool=false, legacy_tuning_geometry::Bool=false,
    mgwr_grouping::Symbol=:intercept_only, residual_shrinkage::Bool=true,
)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be :smoke or :full"))
    nested_covariates && legacy_dem && throw(ArgumentError(
        "--nested-covariates only applies to the joint-covariate path, not --legacy-dem",
    ))
    smoke = mode == :smoke
    # Repeated runs, and nested-selection runs, get their own directory so they never overwrite
    # the canonical single-partition/fixed-selection outputs that the verify scripts target.
    outdir = joinpath(
        ROOT, "output", (legacy_dem ?
            (smoke ? "interpolation_benchmark_smoke_dem" : "interpolation_benchmark_full_dem") :
            (smoke ? "interpolation_benchmark_smoke_joint_covariates" :
                "interpolation_benchmark_full_joint_covariates")) *
            (nested_covariates ? "_nested" : "") *
            (local_grid ? "_localgrid" : "") *
            (stratified_tuning_weights ? "_strattuning" : "") *
            (legacy_tuning_geometry ? "_loocvtuning" : "") *
            # Keyed on the historical layout rather than on the current default, so the pre-fix
            # run directories keep their names and a default run lands somewhere new instead of
            # overwriting the before-picture.
            (mgwr_grouping === :split ? "" : "_mgwr$(mgwr_grouping)") *
            (residual_shrinkage ? "" : "_noshrink") *
            (repeats > 1 ? "_repeats$(repeats)" : ""),
    )
    mkpath(outdir)
    # One adaptive grid for the whole GWR family - direct `gwr`, `hurdle_gwr`, residual `gwr`,
    # `mixed_gwr` and `mgwr` all search it, so a comparison between them is at equal search
    # budget. It used to be two: `mger.bw_adaptive` topped out at 120 while the joint models'
    # grid reached 160, which direct `gwr` could never select.
    #
    # The floor matters most. The canonical grid starts at 30 nearest neighbours (~37 km with
    # 237 Hubei stations) while IDW/ADW may use 8 (~18 km), and 72% of all selected bandwidths
    # landed on a grid endpoint with the CV curve still running into it - the grid, not the
    # data, was picking the bandwidth (`output/benchmark_diagnostics/*/bandwidth_saturation.csv`).
    # `--local-grid` drops the GWR floor to IDW/ADW's locality so the two compare at equal
    # bandwidth.
    #
    # 160 is gone from the ceiling. Its only real job was "effectively global": it exceeded the
    # inner selection split's training size, so `gw_weight` fell through its `dn > 1` branch to
    # a near-uniform fit during scoring and then refitted as the 160th-nearest-of-~190 at
    # prediction time - the same `bw` naming two different models.
    # `InterpolationBenchmarkConfig.bw_include_global` now offers global as an explicit
    # candidate instead, and the 120 ceiling stays below the inner training size so the
    # fallthrough branch is never reached.
    bw_adaptive = local_grid ? [8.0, 12.0, 16.0, 20.0, 30.0, 50.0, 80.0, 120.0] :
        (smoke ? [30.0, 80.0] : [30.0, 50.0, 80.0, 120.0])
    mger = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
                    "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv" :
                    "hubei_gpm_hourly_2022_2024_full_aligned.csv",
            ),
            "GSMaP" => joinpath(
                STUDY_DATA,
                smoke ? "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv" :
                    "hubei_gsmap_hourly_2022_2024_full_aligned.csv",
            ),
        ),
        outdir=outdir,
        kernels=smoke ? [GAUSSIAN, BISQUARE] :
            [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR],
        bw_adaptive=bw_adaptive,
        bw_fixed_km=local_grid ? [5.0, 10.0, 20.0, 30.0, 50.0] :
            (smoke ? [30.0, 50.0] : [10.0, 20.0, 30.0, 50.0]),
        rain_threshold=0.1,
        use_loocv_eval=true,
        # Full mode spans Jan 2022 - Dec 2024. FY4B's strict-completeness QC yields zero
        # usable hours in every November/December on record (2022-2024), so this window
        # is a no-op vs. a narrower Nov-ending window - stated explicitly rather than
        # relying on that gap.
        analysis_start=smoke ? DateTime(2022, 6, 1, 9) : DateTime(2022, 1, 1, 9),
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2025, 1, 1, 8),
        # FY4B's off-season NC files rarely form a strict-complete hour (near-zero in
        # Nov/Dec and Jan-Mar), so the common-time count grows from 8067 (Jun-Sep only)
        # to 11426, not to a full 3-year hourly count.
        expected_common_time_count=smoke ? nothing : 11426,
    )
    schemes = smoke && !with_random ? [:balanced_spatial] : [:balanced_spatial, :random]
    dem = legacy_dem ? DEMTerrainExperiment.DEMExperimentConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours=smoke ? 20 : 100,
        k=5,
        seed=20260815,
        bandwidth_candidates=[30, 50, 80, 120, 160],
        screen_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
        vif_threshold=5.0,
        tolerance=1e-5,
        max_iterations=200,
    ) : nothing
    joint = legacy_dem ? nothing : JointCovariateBenchmarkConfig(
        spec_path=nested_covariates ? nothing : joinpath(
            ROOT, "output", "joint_variable_selection", "full",
            "joint_final_full_data_spec.csv",
        ),
        terrain_path=joinpath(
            ROOT, "data", "processed", "covariates", "station_terrain.csv",
        ),
        era5_paths=Dict(year_value => joinpath(
            ROOT, "data", "processed", "covariates", "era5_land",
            "era5_land_station_hourly_utc_$(year_value).csv",
        ) for year_value in 2022:2024),
        ndvi_path=joinpath(
            ROOT, "data", "processed", "covariates",
            "station_ndvi_16day_2022_2024.csv",
        ),
        # Same grid as `mger.bw_adaptive` above: the joint models are part of the GWR family and
        # are compared against the rest of it, so they must not search a wider one.
        bandwidth_candidates=Int.(bw_adaptive),
        feature_time_offset_hours=9,
        max_ndvi_age_days=32,
        wet_threshold=0.1,
        tolerance=1e-5,
        # 200 was binding: the back-fit's median sweep count sat at 150-194 once smooth covariate
        # groups entered, so MGWR discarded whole hours for want of a few more iterations. See the
        # field comment in `JointCovariateBenchmarkConfig`.
        max_iterations=1000,
        mgwr_spatial_grouping=mgwr_grouping,
    )
    # Same permutation rigor as the per-fold nested selection is asked to match the legacy DEM
    # path's precedent (999 full / 99 smoke), rather than reducing it for speed.
    joint_selection = nested_covariates ? JointSelectionConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours=smoke ? 20 : 100,
        k=5,
        # Deliberately NOT extended by --local-grid. This grid feeds the ERA5/NDVI/DEM
        # local-vs-global role tests (JointVariableSelection._family_configs), i.e. which
        # covariates get selected. Holding it fixed keeps the same covariate set as the
        # baseline nested run, so a change in RMSE is attributable to the fitting bandwidth
        # alone rather than to a different set of variables.
        bandwidth_candidates=[30, 50, 80, 120, 160],
        independent_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
        vif_threshold=5.0,
    ) : nothing
    # The joint path's significance test was disabled only because fixed full-data selection
    # made it unreliable; nested selection removes that leak, so it gets the same bootstrap
    # budget as the legacy DEM path.
    joint_bootstrap_reps = nested_covariates ? (smoke ? 200 : 2000) : 0
    return InterpolationBenchmarkConfig(
        mger=mger,
        terrain_path=legacy_dem ? joinpath(
            ROOT, "data", "processed", "covariates", "station_terrain.csv",
        ) : nothing,
        dem=dem,
        joint_covariates=joint,
        joint_selection=joint_selection,
        k=5,
        seed=20260627,
        seeds=repeat_seeds(repeats),
        cv_schemes=schemes,
        idw_powers=smoke ? [1.5, 2.0, 2.5] : [1.0, 1.5, 2.0, 2.5, 3.0],
        neighbor_candidates=smoke ? Union{Nothing,Int}[16, 32] :
            Union{Nothing,Int}[8, 16, 32, 64, nothing],
        tps_smooth_candidates=smoke ? [1e-3, 1e-2, 1e-1] :
            [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
        min_tuning_coverage=0.95,
        tuning_max_times=smoke ? 72 : 336,
        tuning_time_weighting=stratified_tuning_weights ? :stratified : :uniform,
        tuning_geometry=legacy_tuning_geometry ? :loocv : :inner_spatial,
        # `[1.0]` is the historical unshrunk correction; the ladder costs no extra model fits.
        residual_shrinkage_candidates=residual_shrinkage ?
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] : [1.0],
        event_thresholds=[0.1, 2.5, 8.0, 16.0],
        bootstrap_reps=legacy_dem ? (smoke ? 200 : 2000) : joint_bootstrap_reps,
    )
end

"""
Parse `--mgwr-grouping <split|shared|intercept_only>`; absent means `:intercept_only`.
`--mgwr-grouping split` restores the historical three-single-column-group layout.
"""
function parse_mgwr_grouping(args)
    index = findfirst(startswith("--mgwr-grouping"), args)
    index === nothing && return :intercept_only
    text = args[index] == "--mgwr-grouping" ? get(args, index + 1, "intercept_only") :
        split(args[index], "=", limit=2)[2]
    grouping = Symbol(text)
    grouping in (:split, :shared, :intercept_only) || throw(ArgumentError(
        "--mgwr-grouping expects split, shared, or intercept_only, got $text",
    ))
    return grouping
end

"""Parse `--repeats N`; absent means a single partition (the historical behaviour)."""
function parse_repeats(args)
    index = findfirst(startswith("--repeats"), args)
    index === nothing && return 1
    text = args[index] == "--repeats" ? get(args, index + 1, "1") :
        split(args[index], "=", limit=2)[2]
    repeats = tryparse(Int, text)
    (repeats === nothing || repeats < 1) &&
        throw(ArgumentError("--repeats expects a positive integer, got $text"))
    return repeats
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    with_random = "--with-random" in args
    legacy_dem = "--legacy-dem" in args
    # :full defaults to the leak-free nested per-fold covariate selection (closes the
    # full-data-reused-across-every-fold leak); --no-nested-covariates reproduces the old
    # fixed full-data comparison path, and --legacy-dem never defaults into nested selection.
    nested_covariates = if "--nested-covariates" in args
        true
    elseif "--no-nested-covariates" in args
        false
    else
        mode == :full && !legacy_dem
    end
    repeats = parse_repeats(args)
    local_grid = "--local-grid" in args
    # Hyperparameters are tuned on a wet-oversampled subsample but scored on every hour, so the
    # tuning RMSE runs ~3x the metric it estimates. --stratified-tuning-weights reweights the
    # subsample to remove that. It is opt-in, not the default: the level error turned out to be
    # nearly a constant multiplier that cancels in the ranking, so correcting it changes
    # published numbers without improving selection. See scripts/verify_tuning_time_weighting.jl.
    stratified_tuning_weights = "--stratified-tuning-weights" in args
    # Candidates are now scored by predicting onto an inner spatial split of the training
    # stations, matching the geometry the results table reports. The old leave-one-out criterion
    # measured interpolation next to a retained gauge and picked bandwidths 19-101% worse on the
    # reported metric; --legacy-tuning-geometry restores it for reproducing pre-fix numbers.
    legacy_tuning_geometry = "--legacy-tuning-geometry" in args
    # How MGWR splits the spatial trend into back-fitting groups. The coordinate columns are
    # near-constant inside a local window, so they couple with the intercept and with any smooth
    # covariate; that coupling drives how many sweeps the back-fit needs and whether it converges
    # at all. See the field comment in `JointCovariateBenchmarkConfig`.
    mgwr_grouping = parse_mgwr_grouping(args)
    # GWR is unbiased and never shrinks, so its residual correction comes out the right shape at
    # the wrong magnitude - the measured optimal rescale is 0.47-0.84 depending on the product.
    # The scale is now selected with the bandwidth on the inner spatial split, at no extra fitting
    # cost; --no-residual-shrinkage pins it to 1 and reproduces the unshrunk numbers.
    residual_shrinkage = !("--no-residual-shrinkage" in args)
    cfg = benchmark_config(mode; with_random, legacy_dem, repeats, nested_covariates,
        local_grid, stratified_tuning_weights, legacy_tuning_geometry, mgwr_grouping,
        residual_shrinkage)
    !legacy_dem && Threads.nthreads() == 1 && @warn(
        "Joint dynamic models are compute intensive; use julia -t auto for parallel hourly fits",
    )
    nested_covariates && @warn(
        "Nested per-fold covariate selection reruns variable selection for every " *
        "(scheme, product, fold, repeat) cell — expect substantially longer runtime than " *
        "the fixed full-data selection path.",
    )
    println("Running interpolation benchmark: mode=$mode, repeats=$repeats, " *
        "nested_covariates=$nested_covariates, local_grid=$local_grid, " *
        "tuning_time_weighting=$(cfg.tuning_time_weighting), " *
        "tuning_geometry=$(cfg.tuning_geometry), " *
        "mgwr_grouping=$mgwr_grouping, " *
        "residual_shrinkage=$residual_shrinkage, " *
        "output=$(cfg.mger.outdir)")
    result = run_interpolation_benchmark(cfg)
    println("Finished: $(nrow(result.metrics)) metric rows, $(nrow(result.scans)) scan rows")
    if nrow(result.claim) > 0
        println(result.claim)
    end
    if repeats > 1
        println("Across-partition spread: metrics_repeat_summary.csv, method_rank_stability.csv")
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
