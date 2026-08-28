#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]
const REMOVED_METHODS = ["residual_idw", "residual_adw", "residual_tps"]

function main(args=ARGS)
    mode = isempty(args) ? "smoke" : lowercase(args[1])
    mode in ("smoke", "full") || throw(ArgumentError("mode must be smoke or full"))
    outdir = joinpath(ROOT, "output", mode == "smoke" ?
        "interpolation_benchmark_smoke_dem" :
        "interpolation_benchmark_full_dem")

    required = [
        "metrics_stratified.csv", "metrics_folds.csv", "metrics_pooled.csv",
        "parameter_scan.csv", "paired_comparisons.csv", "run_status.csv",
        "claim_assessment.csv", "benchmark_scope.csv", "global_common_time_qc.csv",
        "dem_fold_quality_control.csv", "dem_fold_correlation.csv",
        "dem_fold_aspect_joint_test.csv", "dem_fold_vif.csv",
        "dem_fold_monthly_correlation.csv", "dem_fold_spatial_variability.csv",
        "dem_fold_roles.csv", "dem_role_stability.csv", "dem_final_full_data_spec.csv",
        "dem_bandwidths.csv", "model_status.csv",
    ]
    for filename in required
        @assert isfile(joinpath(outdir, filename)) "missing output: $filename"
    end

    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    @assert nrow(status) > 0 "empty run status"
    @assert !any(status.status .== "failed") "one or more fold/method runs failed"
    partial = filter(:status => ==("partial"), status)
    @assert all(partial.prediction_coverage .> 0) "partial runs must retain valid predictions"
    # `auto` is reported alongside the fitted methods but is never *run* on this legacy DEM
    # path: `inner_selection_prediction` carries no `dem_context`, so it records "skipped" rather
    # than a failure. Asserting the reason keeps a genuine auto breakage from hiding here.
    @assert Set(status.method) == Set(vcat(filter(!=("raw"), METHODS), "auto"))
    auto_status = filter(:method => ==("auto"), status)
    @assert nrow(auto_status) > 0 "auto rows missing from run status"
    @assert all(auto_status.status .== "skipped") "auto should not run on the legacy DEM path"
    @assert all(row -> occursin("legacy DEM path", row.error), eachrow(auto_status))
    @assert !isfile(joinpath(outdir, "auto_selection.csv")) "auto must not select on the DEM path"
    @assert all(name -> name in names(status),
        ["dem_variable_count", "dem_selection_status", "dem_variables", "dem_roles"])
    residual_status = filter(:method => in(["residual_gwr", "mixed_gwr", "mgwr"]), status)
    # `repeat` distinguishes repeated cross-validation partitions, which share these tables.
    @assert all(group -> length(unique(group.dem_variables)) == 1,
        groupby(residual_status, [:repeat, :scheme, :product, :fold]))
    scans = CSV.read(joinpath(outdir, "parameter_scan.csv"), DataFrame)
    @assert all(group -> count(group.selected) == 1,
        groupby(scans, [:repeat, :scheme, :product, :fold, :mode, :method, :group]))
    @assert all(row -> row.mode == "direct" || row.method in ("gwr", "mixed_gwr", "mgwr"),
        eachrow(scans))

    roles = CSV.read(joinpath(outdir, "dem_fold_roles.csv"), DataFrame)
    @assert Set(roles.variable_group) == Set(["elevation", "slope", "aspect"])
    @assert all(role -> role in ("local", "global", "uncertain", "not_selected"), roles.role)
    final_spec = CSV.read(joinpath(outdir, "dem_final_full_data_spec.csv"), DataFrame)
    @assert Set(final_spec.product) == Set(["FY4B", "GPM", "GSMaP"])

    metrics = CSV.read(joinpath(outdir, "metrics_pooled.csv"), DataFrame)
    overall = filter([:group, :level] => (group, level) -> group == "overall" && level == "all", metrics)
    @assert Set(overall.method) == Set(METHODS)
    @assert all(isfinite, overall.RMSE)
    @assert all(overall.coverage .> 0)

    schemes = unique(status.scheme)
    products = sort(unique(status.product))
    for scheme in schemes, product in products
        product_dir = joinpath(outdir, scheme, lowercase(product))
        reference_names = nothing
        reference_rows = nothing
        for method in METHODS
            path = joinpath(product_dir, "oof_$(method).csv")
            @assert isfile(path) "missing OOF output: $path"
            table = CSV.read(path, DataFrame)
            if reference_names === nothing
                reference_names = names(table)
                reference_rows = nrow(table)
            else
                @assert names(table) == reference_names "station order differs across methods"
                @assert nrow(table) == reference_rows "time count differs across methods"
            end
        end
        for method in REMOVED_METHODS
            @assert !isfile(joinpath(product_dir, "oof_$(method).csv")) "obsolete OOF output: $method"
        end
        @assert isfile(joinpath(product_dir, "common_evaluation_mask.csv"))
    end
    println("Interpolation benchmark verification passed: $outdir")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
