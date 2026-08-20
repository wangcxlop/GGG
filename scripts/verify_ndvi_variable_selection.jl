#!/usr/bin/env julia

using CSV, DataFrames, Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))

function main(args=ARGS)
    mode = isempty(args) ? "smoke" : lowercase(args[1])
    mode in ("smoke", "full") || throw(ArgumentError("mode must be smoke or full"))
    outdir = joinpath(ROOT, "output", "ndvi_variable_selection", mode)
    required = [
        "ndvi_data_qc.csv", "ndvi_station_qc.csv", "ndvi_alignment_qc.csv",
        "ndvi_alignment_intervals.csv",
        "ndvi_fold_quality_control.csv", "ndvi_station_period_residual.csv",
        "ndvi_fold_association.csv", "ndvi_fold_monthly_direction.csv",
        "ndvi_fold_vif.csv", "ndvi_fold_bandwidth_scan.csv",
        "ndvi_fold_spatial_variability.csv", "ndvi_fold_roles.csv",
        "ndvi_role_stability.csv", "ndvi_final_full_data_spec.csv",
        "ndvi_cross_product_consensus.csv", "spatial_folds.csv", "run_status.csv",
    ]
    for filename in required
        @assert isfile(joinpath(outdir, filename)) "missing output: $filename"
    end
    data_qc = CSV.read(joinpath(outdir, "ndvi_data_qc.csv"), DataFrame)
    @assert data_qc.station_count[1] == 237
    @assert data_qc.all_237_stations[1]
    @assert data_qc.complete[1]
    station_qc = CSV.read(joinpath(outdir, "ndvi_station_qc.csv"), DataFrame)
    @assert nrow(station_qc) == 237
    @assert any(station_qc.nonland_pixel_flag) "non-land pixels must be audited rather than removed"
    alignment = CSV.read(joinpath(outdir, "ndvi_alignment_qc.csv"), DataFrame)
    @assert alignment.max_age_days[1] == 32
    @assert alignment.interpolation[1] == "none"
    @assert alignment.effective_time_rule[1] == "observation_date_plus_1_day_at_00_BJT"
    intervals = CSV.read(joinpath(outdir, "ndvi_alignment_intervals.csv"), DataFrame)
    @assert all(intervals.valid_until_exclusive .<= intervals.effective_time .+ Day(32))
    association = CSV.read(joinpath(outdir, "ndvi_fold_association.csv"), DataFrame)
    @assert Set(association.variable) == Set(["ndvi"])
    @assert all(association.qvalue .== association.pvalue)
    @assert all((0 .<= association.pvalue) .& (association.pvalue .<= 1))
    vif = CSV.read(joinpath(outdir, "ndvi_fold_vif.csv"), DataFrame)
    @assert all(==("not_applicable_single_predictor"), vif.status)
    roles = CSV.read(joinpath(outdir, "ndvi_fold_roles.csv"), DataFrame)
    @assert all(role -> role in ("local", "global", "uncertain", "not_selected"), roles.role)
    spec = CSV.read(joinpath(outdir, "ndvi_final_full_data_spec.csv"), DataFrame)
    @assert nrow(spec) == 3
    @assert Set(spec.product) == Set(["FY4B", "GPM", "GSMaP"])
    folds = CSV.read(joinpath(outdir, "spatial_folds.csv"), DataFrame)
    @assert sort(unique(folds.fold)) == collect(1:5)
    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    @assert nrow(status) == 18
    @assert all(==("ok"), status.status)
    println("NDVI variable-selection verification passed: $outdir")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
