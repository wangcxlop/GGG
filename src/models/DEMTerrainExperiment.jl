module DEMTerrainExperiment

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics

export DEMExperimentConfig, bh_adjust, terrain_screen, spatial_variability_test
export mixed_gwr_predict, multiscale_gwr_predict, select_mixed_bandwidth
export select_multiscale_bandwidths, run_dem_experiment
export mean_wet_residual, monthly_correlation_rows, terrain_model_designs
export terrain_groups, terrain_columns

const TERRAIN_GROUPS = ("elevation", "slope", "aspect")
const TERRAIN_COLUMNS = Dict(
    "elevation" => [:elevation_m],
    "slope" => [:slope_deg],
    "aspect" => [:aspect_sin, :aspect_cos],
)

terrain_groups() = collect(TERRAIN_GROUPS)
terrain_columns(group::String) = copy(TERRAIN_COLUMNS[group])

Base.@kwdef struct DEMExperimentConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    k::Int = 5
    seed::Int = 20260815
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    screen_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
    tolerance::Float64 = 1e-5
    max_iterations::Int = 200
end

# Split by concern, in `dem/`, all `include`d into this module rather than made modules of their
# own - so every name stays where callers already reach it, including the private `_local_hat`,
# `_bisquare_kernel`, `_global_projection` and `_haversine_matrix` that the tests call qualified.
#
# `experiment.jl` is deliberately last and deliberately alone: `run_dem_experiment` is a study
# orchestrator that writes ~15 fixed-name CSVs, and everything before it is the reusable
# screening/GWR library it happens to sit beside. The boundary is now visible in the directory
# listing rather than having to be discovered by reading 1500 lines.
#
# Note that `local_weights.jl` and `backfit.jl` hold two independent takes on "local GWR weights"
# that never meet: the smoother-matrix one, used only by `spatial_variability_test`, and the
# hat-matrix one the mixed/multiscale predictors are built on.

include("dem/screening.jl")
include("dem/local_weights.jl")
include("dem/spatial_variability.jl")
include("dem/backfit.jl")
include("dem/predict.jl")
include("dem/summaries.jl")
include("dem/experiment.jl")

end
