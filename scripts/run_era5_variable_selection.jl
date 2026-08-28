#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates

include(joinpath(ROOT, "src", "MGERPipeline.jl"))
include(joinpath(ROOT, "src", "ERA5VariableSelection.jl"))
using .ERA5VariableSelection

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const ERA5_DATA = joinpath(ROOT, "data", "processed", "covariates", "era5_land")

function experiment_inputs(mode::Symbol)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be smoke or full"))
    smoke = mode == :smoke
    outdir = joinpath(ROOT, "output", "era5_variable_selection", string(mode))
    mkpath(outdir)
    mger = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(STUDY_DATA, smoke ?
                "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
                "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
        ),
        outdir=outdir,
        analysis_start=DateTime(2022, 6, 1, 9),
        analysis_end=smoke ? DateTime(2022, 7, 1, 8) : DateTime(2024, 10, 1, 8),
        expected_common_time_count=smoke ? nothing : 8067,
    )
    cfg = ERA5SelectionConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours=smoke ? 20 : 100,
        min_stations_per_time=12,
        k=5,
        seed=20260816,
        bandwidth_candidates=[30, 50, 80, 120, 160],
        association_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
        vif_threshold=5.0,
        time_offset_hours=9,
    )
    return mger, cfg
end

function main(args=ARGS)
    mode = isempty(args) ? :smoke : Symbol(lowercase(args[1]))
    mger, cfg = experiment_inputs(mode)
    station_meta = load_station_meta(mger.station_meta_path)
    products, ids, product_data = load_global_common_product_data(mger)
    lonlat = build_X_lonlat(station_meta, ids)
    times = product_data[first(products)].times
    Yobs = Matrix{Float64}(product_data[first(products)].Y_obs)
    satellite = Dict(product => Matrix{Float64}(product_data[product].Y_sat) for product in products)
    years = sort(unique(year(time - Hour(cfg.time_offset_hours)) for time in times))
    annual_paths = Dict(year_value => joinpath(
        ERA5_DATA, "era5_land_station_hourly_utc_$(year_value).csv",
    ) for year_value in years)
    println("Loading ERA5-Land with feature_time = time_utc + $(cfg.time_offset_hours) hours")
    loaded = load_era5_panel(annual_paths, ids, times; offset_hours=cfg.time_offset_hours)
    println("Running standalone ERA5 selection: mode=$mode, stations=$(length(ids)), times=$(length(times))")
    result = run_era5_variable_selection(
        cfg, products, ids, times, Yobs, satellite, loaded.values, lonlat;
        data_qc=loaded.qc,
    )
    println("Finished ERA5 variable selection: $(cfg.outdir)")
    println(result.consensus)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
