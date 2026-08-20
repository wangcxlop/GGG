#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))

using MixedGWR
using CSV, DataFrames, Dates

include(joinpath(ROOT, "src", "MGERPipeline.jl"))
include(joinpath(ROOT, "src", "JointVariableSelection.jl"))
using .JointVariableSelection

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const COVARIATE_DATA = joinpath(ROOT, "data", "processed", "covariates")
const ERA5_DATA = joinpath(COVARIATE_DATA, "era5_land")
const NDVI_PATH = joinpath(COVARIATE_DATA, "station_ndvi_16day_2022_2024.csv")
const TERRAIN_PATH = joinpath(COVARIATE_DATA, "station_terrain.csv")

function aligned_terrain(path::String, ids::Vector{String})
    terrain = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    allunique(terrain.station_id) || error("Duplicate station IDs in terrain table")
    row_by_id = Dict(id => row for (row, id) in enumerate(terrain.station_id))
    missing_ids = filter(id -> !haskey(row_by_id, id), ids)
    isempty(missing_ids) || error("Terrain table is missing $(length(missing_ids)) stations")
    return terrain[[row_by_id[id] for id in ids], :]
end

function prerequisite_audit()
    roots = Dict(
        "dem" => joinpath(ROOT, "output", "dem_variable_selection", "full"),
        "era5" => joinpath(ROOT, "output", "era5_variable_selection", "full"),
        "ndvi" => joinpath(ROOT, "output", "ndvi_variable_selection", "full"),
    )
    rows = NamedTuple[]

    dem_qc_path = joinpath(roots["dem"], "data_qc.csv")
    isfile(dem_qc_path) || error("Missing DEM full prerequisite: $dem_qc_path")
    dem_qc = CSV.read(dem_qc_path, DataFrame)
    dem_values = Dict(String(row.key) => String(row.value) for row in eachrow(dem_qc))
    dem_ok = get(dem_values, "time_count", "") == "8067" &&
        get(dem_values, "screen_permutations", "") == "999" &&
        get(dem_values, "spatial_permutations", "") == "999"
    push!(rows, (; family="dem", output_dir=roots["dem"], time_count=8067,
        association_permutations=999, spatial_permutations=999,
        verified=dem_ok, detail=dem_ok ? "ok" : "configuration_mismatch"))

    era_qc_path = joinpath(roots["era5"], "era5_data_qc.csv")
    era_status_path = joinpath(roots["era5"], "run_status.csv")
    all(isfile, [era_qc_path, era_status_path]) || error("Missing ERA5 full prerequisite outputs")
    era_qc = CSV.read(era_qc_path, DataFrame)
    era_status = CSV.read(era_status_path, DataFrame)
    era_ok = era_qc.target_time_count[1] == 8067 && Bool(era_qc.complete[1]) &&
        all(==("ok"), String.(era_status.status))
    push!(rows, (; family="era5", output_dir=roots["era5"], time_count=8067,
        association_permutations=999, spatial_permutations=999,
        verified=era_ok, detail=era_ok ? "ok" : "quality_or_status_failure"))

    ndvi_qc_path = joinpath(roots["ndvi"], "ndvi_alignment_qc.csv")
    ndvi_status_path = joinpath(roots["ndvi"], "run_status.csv")
    all(isfile, [ndvi_qc_path, ndvi_status_path]) || error("Missing NDVI full prerequisite outputs")
    ndvi_qc = CSV.read(ndvi_qc_path, DataFrame)
    ndvi_status = CSV.read(ndvi_status_path, DataFrame)
    ndvi_ok = ndvi_qc.time_count[1] == 8067 &&
        all(==("ok"), String.(ndvi_status.status))
    push!(rows, (; family="ndvi", output_dir=roots["ndvi"], time_count=8067,
        association_permutations=999, spatial_permutations=999,
        verified=ndvi_ok, detail=ndvi_ok ? "ok" : "quality_or_status_failure"))

    audit = DataFrame(rows)
    all(audit.verified) || error("At least one independent full prerequisite failed verification")
    return audit
end

function experiment_inputs(mode::Symbol)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be smoke or full"))
    permutations = mode == :smoke ? 99 : 999
    outdir = joinpath(ROOT, "output", "joint_variable_selection", string(mode))
    mkpath(outdir)
    mger = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(STUDY_DATA,
                "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
        ),
        outdir=outdir,
        analysis_start=DateTime(2022, 6, 1, 9),
        analysis_end=DateTime(2024, 10, 1, 8),
        expected_common_time_count=8067,
    )
    cfg = JointSelectionConfig(
        outdir=outdir, wet_threshold=0.1, min_wet_hours=100,
        min_stations_per_time=12, min_ndvi_wet_hours_per_period=5,
        min_ndvi_periods_per_station=8, max_ndvi_age_days=32,
        k=5, seed=20260818, bandwidth_candidates=[30, 50, 80, 120, 160],
        independent_permutations=permutations, spatial_permutations=permutations,
        q_threshold=0.05, vif_threshold=5.0,
    )
    return mger, cfg
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    prerequisite = prerequisite_audit()
    mger, cfg = experiment_inputs(mode)
    station_meta = load_station_meta(mger.station_meta_path)
    products, ids, product_data = load_global_common_product_data(mger)
    lonlat = build_X_lonlat(station_meta, ids)
    times = product_data[first(products)].times
    length(times) == 8067 || error("Joint experiment requires 8067 common hours")
    Yobs = Matrix{Float64}(product_data[first(products)].Y_obs)
    satellite = Dict(product => Matrix{Float64}(product_data[product].Y_sat)
        for product in products)
    terrain = aligned_terrain(TERRAIN_PATH, ids)
    years = sort(unique(year(time - Hour(9)) for time in times))
    annual_paths = Dict(year_value => joinpath(
        ERA5_DATA, "era5_land_station_hourly_utc_$(year_value).csv",
    ) for year_value in years)
    era5 = JointVariableSelection.ERA5VariableSelection.load_era5_panel(
        annual_paths, ids, times; offset_hours=9,
    )
    era5.qc.complete[1] || error("ERA5 joint input failed completeness checks")
    ndvi = JointVariableSelection.NDVIVariableSelection.load_ndvi_covariates(NDVI_PATH, ids)
    println("Running joint variable selection: mode=$mode, stations=$(length(ids)), times=$(length(times)), permutations=$(cfg.independent_permutations)")
    result = run_joint_variable_selection(
        cfg, products, ids, times, Yobs, satellite, terrain, era5.values,
        ndvi.table, lonlat; prerequisite_audit=prerequisite,
    )
    println("Finished joint variable selection: $(cfg.outdir)")
    println(result.final_spec[result.final_spec.final_included, :])
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
