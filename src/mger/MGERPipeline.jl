using MixedGWR
using CSV, DataFrames, Dates, Statistics, Random
using Distances


"""
Configuration
- The residual model uses target-centred local east/north coordinates in km
- Both adaptive=true and adaptive=false bandwidths are scanned
- The ridge strength on the local spatial slopes is scanned as well
"""
Base.@kwdef struct MGERConfig
	station_meta_path::String
	obs_hourly_wide_path::String
	sat_paths::Dict{String, String} # "FY4B"=>path, "GPM"=>path, "GSMaP"=>path
	outdir::String = "output/mger_final"

	time_col::Symbol = :time
	station_id_col::Symbol = :station_id
	lon_col::Symbol = :lon
	lat_col::Symbol = :lat

	kernels::Vector{Int} = [GAUSSIAN, EXPONENTIAL, BISQUARE, TRICUBE, BOXCAR]
	bw_adaptive::Vector{Float64} = [30.0, 50.0, 80.0, 120.0]
	bw_fixed_km::Vector{Float64} = [10.0, 20.0, 30.0, 50.0]
	slope_ridge_candidates::Vector{Float64} = [1e-8, 1e-6, 1e-4, 1e-2]
	min_scan_coverage::Float64 = 0.95
	rain_threshold::Float64 = 0.1

	use_loocv_eval::Bool = true
	analysis_start::Union{Nothing, DateTime} = nothing
	analysis_end::Union{Nothing, DateTime} = nothing
	expected_common_time_count::Union{Nothing, Int} = nothing
end


function config_for_kernel(cfg::MGERConfig, kernel::Int, outdir::AbstractString)
	return MGERConfig(
		station_meta_path=cfg.station_meta_path,
		obs_hourly_wide_path=cfg.obs_hourly_wide_path,
		sat_paths=copy(cfg.sat_paths),
		outdir=String(outdir),
		time_col=cfg.time_col,
		station_id_col=cfg.station_id_col,
		lon_col=cfg.lon_col,
		lat_col=cfg.lat_col,
		kernels=[kernel],
		bw_adaptive=copy(cfg.bw_adaptive),
		bw_fixed_km=copy(cfg.bw_fixed_km),
		slope_ridge_candidates=copy(cfg.slope_ridge_candidates),
		min_scan_coverage=cfg.min_scan_coverage,
		rain_threshold=cfg.rain_threshold,
		use_loocv_eval=cfg.use_loocv_eval,
		analysis_start=cfg.analysis_start,
		analysis_end=cfg.analysis_end,
		expected_common_time_count=cfg.expected_common_time_count,
	)
end

include(joinpath(@__DIR__, "pipeline/mger_io.jl"))
include(joinpath(@__DIR__, "pipeline/mger_residual_predict.jl"))
include(joinpath(@__DIR__, "pipeline/mger_writers.jl"))
include(joinpath(@__DIR__, "pipeline/mger_splits.jl"))
include(joinpath(@__DIR__, "pipeline/mger_runners.jl"))

"""
Main audited entry points:

    julia --project=. scripts/run_mger_three_products_timealigned.jl
    julia --project=. scripts/run_mger_three_products_timealigned_loocv.jl
"""
