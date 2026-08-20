#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REQUIRED = [
    "prerequisite_full_audit.csv",
    "spatial_folds.csv",
    "joint_independent_candidates.csv",
    "joint_panel_quality_control.csv",
    "joint_scaling.csv",
    "joint_weight_audit.csv",
    "joint_vif.csv",
    "joint_bandwidth_scan.csv",
    "joint_spatial_variability.csv",
    "joint_fold_roles.csv",
    "joint_role_stability.csv",
    "joint_final_full_data_spec.csv",
    "joint_cross_product_consensus.csv",
    "run_status.csv",
]

function verify(mode::Symbol)
    mode in (:smoke, :full) || error("mode must be smoke or full")
    outdir = joinpath(ROOT, "output", "joint_variable_selection", string(mode))
    missing_files = filter(file -> !isfile(joinpath(outdir, file)), REQUIRED)
    isempty(missing_files) || error("Missing joint outputs: $(join(missing_files, ", "))")

    prerequisite = CSV.read(joinpath(outdir, "prerequisite_full_audit.csv"), DataFrame)
    nrow(prerequisite) == 3 || error("Expected three prerequisite families")
    Set(prerequisite.family) == Set(["dem", "era5", "ndvi"]) || error("Prerequisite families differ")
    all(prerequisite.verified) || error("Independent full prerequisite verification failed")
    all(==(8067), prerequisite.time_count) || error("Independent full time count differs")
    all(==(999), prerequisite.association_permutations) || error("Independent full association permutations differ")
    all(==(999), prerequisite.spatial_permutations) || error("Independent full spatial permutations differ")

    candidates = CSV.read(joinpath(outdir, "joint_independent_candidates.csv"), DataFrame)
    nrow(candidates) == 3 * 6 * 9 || error("Expected 162 product/scheme/candidate rows")
    Set(candidates.variable_group) == Set([
        "elevation", "slope", "aspect", "t2m_c", "d2m_c", "u10", "v10", "sp_hpa", "ndvi",
    ]) || error("Joint candidate groups differ")
    all(candidates[candidates.variable_group .== "aspect", :predictor_columns] .==
        "aspect_sin+aspect_cos") || error("Aspect components are not grouped")
    all(candidates[candidates.variable_group .== "ndvi", :predictor_columns] .==
        "ndvi_qc") || error("NDVI source column is not explicit")

    roles = CSV.read(joinpath(outdir, "joint_fold_roles.csv"), DataFrame)
    nrow(roles) == 3 * 6 * 9 || error("Joint role table is incomplete")
    all(in(Set(["not_selected", "local", "global", "uncertain"])), roles.role) ||
        error("Unexpected role value")
    all((.!roles.final_included) .| (roles.independent_selected .& roles.joint_vif_retained .&
        in.(roles.role, Ref(Set(["local", "global"]))))) ||
        error("A final variable bypassed screening, VIF, or role checks")

    final_spec = CSV.read(joinpath(outdir, "joint_final_full_data_spec.csv"), DataFrame)
    nrow(final_spec) == 27 || error("Expected 27 full-data product/group specifications")
    all((.!final_spec.final_included) .| (final_spec.independent_selected .&
        final_spec.joint_vif_retained .& in.(final_spec.role, Ref(Set(["local", "global"]))))) ||
        error("Final full-data configuration is internally inconsistent")

    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    nrow(status) == 18 || error("Expected 18 product/scheme statuses")
    all(==("ok"), status.status) || error("At least one joint product/scheme failed")

    folds = CSV.read(joinpath(outdir, "spatial_folds.csv"), DataFrame)
    nrow(folds) == 237 || error("Expected 237 stations in spatial folds")
    Set(folds.fold) == Set(1:5) || error("Spatial fold IDs differ")

    stability = CSV.read(joinpath(outdir, "joint_role_stability.csv"), DataFrame)
    nrow(stability) == 27 || error("Joint stability table is incomplete")
    all(==(5), stability.fold_count) || error("Each stability row must summarize five folds")

    println("Verified joint variable selection: mode=$mode, output=$outdir")
    println(final_spec[final_spec.final_included, :])
    return true
end

mode = isempty(ARGS) ? :smoke : Symbol(lowercase(ARGS[1]))
verify(mode)
