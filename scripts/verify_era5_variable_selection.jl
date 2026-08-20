#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VARIABLES = Set(["t2m_c", "d2m_c", "u10", "v10", "sp_hpa"])

function main(args=ARGS)
    mode = isempty(args) ? "smoke" : lowercase(args[1])
    mode in ("smoke", "full") || throw(ArgumentError("mode must be smoke or full"))
    outdir = joinpath(ROOT, "output", "era5_variable_selection", mode)
    required = [
        "era5_data_qc.csv", "era5_fold_quality_control.csv", "era5_fold_association.csv",
        "era5_fold_monthly_direction.csv", "era5_fold_vif.csv",
        "era5_fold_bandwidth_scan.csv", "era5_fold_spatial_variability.csv",
        "era5_fold_roles.csv", "era5_role_stability.csv",
        "era5_final_full_data_spec.csv", "era5_cross_product_consensus.csv",
        "spatial_folds.csv", "run_status.csv",
    ]
    for filename in required
        @assert isfile(joinpath(outdir, filename)) "missing output: $filename"
    end
    qc = CSV.read(joinpath(outdir, "era5_data_qc.csv"), DataFrame)
    @assert all(qc.complete)
    @assert all(==(9), qc.feature_time_offset_hours)
    @assert all(==(0), qc.duplicate_keys)
    @assert all(==(0), qc.missing_station_hours)
    @assert all(==(0), qc.nonfinite_values)
    @assert all(qc.annual_complete)
    @assert all(==(0), qc.annual_duplicate_keys)
    @assert all(==(0), qc.annual_missing_station_hours)
    @assert all(==(0), qc.annual_nonfinite_values)
    association = CSV.read(joinpath(outdir, "era5_fold_association.csv"), DataFrame)
    @assert Set(association.variable) == VARIABLES
    @assert all((0 .<= association.pvalue) .& (association.pvalue .<= 1))
    @assert all((0 .<= association.qvalue) .& (association.qvalue .<= 1))
    roles = CSV.read(joinpath(outdir, "era5_fold_roles.csv"), DataFrame)
    @assert Set(roles.variable) == VARIABLES
    @assert all(role -> role in ("local", "global", "uncertain", "not_selected"), roles.role)
    vif = CSV.read(joinpath(outdir, "era5_fold_vif.csv"), DataFrame)
    @assert all(name -> name in names(vif), ["product", "scheme", "fold", "iteration", "variable", "vif", "removed", "reason"])
    bandwidth = CSV.read(joinpath(outdir, "era5_fold_bandwidth_scan.csv"), DataFrame)
    @assert all(name -> name in names(bandwidth), ["product", "scheme", "fold", "neighbors", "loocv_rmse", "available"])
    variability = CSV.read(joinpath(outdir, "era5_fold_spatial_variability.csv"), DataFrame)
    @assert all(name -> name in names(variability), ["product", "scheme", "fold", "variable", "role", "status"])
    spec = CSV.read(joinpath(outdir, "era5_final_full_data_spec.csv"), DataFrame)
    @assert nrow(spec) == 3 * length(VARIABLES)
    @assert Set(spec.variable) == VARIABLES
    folds = CSV.read(joinpath(outdir, "spatial_folds.csv"), DataFrame)
    @assert sort(unique(folds.fold)) == collect(1:5)
    consensus = CSV.read(joinpath(outdir, "era5_cross_product_consensus.csv"), DataFrame)
    @assert Set(consensus.variable) == VARIABLES
    println("ERA5 variable-selection verification passed: $outdir")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
