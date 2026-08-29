using CSV, DataFrames, Dates, LinearAlgebra, Random, Statistics
using MixedGWR

# Every standalone module is loaded through the shared loader, which puts exactly one copy of
# each into `Main`. Loading a module file twice compiles a second, type-incompatible copy of it -
# `DEMTerrainExperiment.jl` used to be included here and again inside `JointVariableSelection`
# and `JointCovariateModels`. See src/load_modules.jl.
include(joinpath(@__DIR__, "load_modules.jl"))
load_pipeline("MGERPipeline")
load_standalone_modules(
    "TraditionalInterpolation", "DEMTerrainExperiment", "JointCovariateModels",
    "ERA5VariableSelection", "NDVIVariableSelection", "JointVariableSelection",
)

using .TraditionalInterpolation
using .DEMTerrainExperiment: DEMExperimentConfig, terrain_screen, spatial_variability_test
using .DEMTerrainExperiment: mean_wet_residual, monthly_correlation_rows
using .DEMTerrainExperiment: terrain_model_designs, terrain_groups, terrain_columns
using .DEMTerrainExperiment: mixed_gwr_predict, multiscale_gwr_predict
using .DEMTerrainExperiment: select_mixed_bandwidth, select_multiscale_bandwidths
using .JointCovariateModels
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
include(joinpath(@__DIR__, "InterpolationBenchmarkHurdle.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkTuning.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkMetrics.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkBootstrap.jl"))
include(joinpath(@__DIR__, "InterpolationBenchmarkRun.jl"))
