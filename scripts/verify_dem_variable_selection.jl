#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))

function main(args=ARGS)
    mode = isempty(args) ? "smoke" : lowercase(args[1])
    mode in ("smoke", "full") || throw(ArgumentError("mode must be smoke or full"))
    outdir = joinpath(ROOT, "output", "dem_variable_selection", mode)
    required = [
        "data_qc.csv", "station_residual_summary.csv", "correlation_screen.csv",
        "vif.csv", "monthly_correlation_direction.csv", "spatial_variability.csv",
        "spatial_bandwidth_scan.csv", "cross_product_consensus.csv", "spatial_folds.csv",
        "validation_metrics.csv", "bandwidth_scan.csv", "run_status.csv",
    ]
    for filename in required
        @assert isfile(joinpath(outdir, filename)) "missing output: $filename"
    end
    screen = CSV.read(joinpath(outdir, "correlation_screen.csv"), DataFrame)
    @assert Set(screen.variable_group) == Set(["elevation", "slope", "aspect"])
    @assert all((0 .<= screen.pvalue) .& (screen.pvalue .<= 1))
    @assert all((0 .<= screen.qvalue) .& (screen.qvalue .<= 1))
    roles = CSV.read(joinpath(outdir, "spatial_variability.csv"), DataFrame)
    @assert all(role -> role in ("local", "global"), roles.role)
    consensus = CSV.read(joinpath(outdir, "cross_product_consensus.csv"), DataFrame)
    @assert Set(consensus.variable_group) == Set(["elevation", "slope", "aspect"])
    folds = CSV.read(joinpath(outdir, "spatial_folds.csv"), DataFrame)
    @assert sort(unique(folds.fold)) == collect(1:5)
    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    @assert all(method -> method in ("residual_gwr", "mixed_gwr", "multiscale_gwr"),
        filter(!=("all"), status.method))
    metrics = CSV.read(joinpath(outdir, "validation_metrics.csv"), DataFrame)
    overall = filter([:quantity, :stratum] =>
        (quantity, stratum) -> quantity == "corrected_precipitation" && stratum == "all", metrics)
    @assert all(overall.coverage .>= 0.95) "one or more models have insufficient coverage"
    println("DEM variable-selection verification passed: $outdir")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
