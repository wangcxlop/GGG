#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using DataFrames, Dates

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules("NDVIVariableSelection")
using Main.NDVIVariableSelection

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const NDVI_PATH = joinpath(
    ROOT, "data", "processed", "covariates", "station_ndvi_16day_2022_2024.csv",
)

function experiment_inputs(mode::Symbol)
    mode in (:smoke, :full) || throw(ArgumentError("mode must be smoke or full"))
    smoke = mode == :smoke
    outdir = joinpath(ROOT, "output", "ndvi_variable_selection", string(mode))
    mkpath(outdir)
    mger = MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA, "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
        ),
        outdir=outdir,
        analysis_start=DateTime(2022, 6, 1, 9),
        analysis_end=DateTime(2024, 10, 1, 8),
        expected_common_time_count=8067,
    )
    cfg = NDVISelectionConfig(
        outdir=outdir,
        wet_threshold=0.1,
        min_wet_hours_per_period=5,
        min_periods_per_station=8,
        min_stations_per_period=12,
        max_age_days=32,
        k=5,
        seed=20260817,
        bandwidth_candidates=[30, 50, 80, 120, 160],
        association_permutations=smoke ? 99 : 999,
        spatial_permutations=smoke ? 99 : 999,
        q_threshold=0.05,
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
    ndvi = load_ndvi_covariates(NDVI_PATH, ids)
    println("Running standalone NDVI selection: mode=$mode, stations=$(length(ids)), times=$(length(times))")
    result = run_ndvi_variable_selection(
        cfg, products, ids, times, Yobs, satellite, ndvi.table, lonlat;
        data_qc=ndvi.data_qc, station_qc=ndvi.station_qc,
    )
    println("Finished NDVI variable selection: $(cfg.outdir)")
    println(result.consensus)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
