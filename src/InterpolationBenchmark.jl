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

# The benchmark's implementation is split across concern-specific files (config, fold-splitting,
# DEM, joint-covariates, predictors, tuning, metrics, bootstrap, orchestrator) purely for
# readability; each is a plain top-level fragment, not a module, so include order only needs to
# put shared consts/structs before the files that reference them as call-time names.
include(joinpath(@__DIR__, "InterpolationBenchmarkConfig.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkFolds.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkDEM.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkJoint.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkPredictors.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkTuning.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkMetrics.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkBootstrap.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkRun.jl"))
