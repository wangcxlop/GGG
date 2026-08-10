#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const METHODS = [
    "raw", "idw", "adw", "tps", "gwr",
    "residual_idw", "residual_adw", "residual_tps", "residual_gwr",
]

function main(args=ARGS)
    mode = isempty(args) ? "smoke" : lowercase(args[1])
    mode in ("smoke", "full") || throw(ArgumentError("mode must be smoke or full"))
    outdir = joinpath(ROOT, "output", mode == "smoke" ?
        "interpolation_benchmark_smoke_202206" : "interpolation_benchmark_full")

    required = [
        "metrics_stratified.csv", "metrics_folds.csv", "metrics_pooled.csv",
        "parameter_scan.csv", "paired_comparisons.csv", "run_status.csv",
        "claim_assessment.csv", "benchmark_scope.csv", "global_common_time_qc.csv",
    ]
    for filename in required
        @assert isfile(joinpath(outdir, filename)) "missing output: $filename"
    end

    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    @assert nrow(status) > 0 "empty run status"
    @assert all(status.status .== "success") "one or more fold/method runs failed"
    scans = CSV.read(joinpath(outdir, "parameter_scan.csv"), DataFrame)
    @assert all(group -> count(group.selected) == 1,
        groupby(scans, [:scheme, :product, :fold, :mode, :method]))

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
        @assert isfile(joinpath(product_dir, "common_evaluation_mask.csv"))
    end
    println("Interpolation benchmark verification passed: $outdir")
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

