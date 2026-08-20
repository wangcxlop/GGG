#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]

function verify(mode::Symbol)
    mode in (:smoke, :full) || error("mode must be smoke or full")
    outdir = joinpath(ROOT, "output", mode == :smoke ?
        "interpolation_benchmark_smoke_joint_covariates" :
        "interpolation_benchmark_full_joint_covariates")
    required = [
        "joint_variable_spec_used.csv", "joint_spec_provenance.csv",
        "joint_era5_input_qc.csv", "joint_fold_scaling.csv",
        "joint_fold_quality_control.csv", "joint_bandwidths.csv",
        "covariate_model_status.csv", "metrics_stratified.csv",
        "metrics_folds.csv", "metrics_pooled.csv", "parameter_scan.csv",
        "run_status.csv", "benchmark_scope.csv",
    ]
    missing_files = filter(file -> !isfile(joinpath(outdir, file)), required)
    isempty(missing_files) || error("Missing joint benchmark outputs: $(join(missing_files, ", "))")
    !isfile(joinpath(outdir, "paired_comparisons.csv")) ||
        error("Joint benchmark must not produce paired significance output")
    !isfile(joinpath(outdir, "claim_assessment.csv")) ||
        error("Joint benchmark must not produce a claim assessment")

    spec = CSV.read(joinpath(outdir, "joint_variable_spec_used.csv"), DataFrame)
    included = filter(:final_included => identity, spec)
    nrow(included) == 7 || error("Expected seven included product-variable rows")
    expected = Set([
        ("FY4B", "d2m_c", "local"), ("FY4B", "u10", "local"),
        ("GPM", "elevation", "global"), ("GPM", "u10", "local"),
        ("GPM", "sp_hpa", "local"), ("GSMaP", "u10", "local"),
        ("GSMaP", "sp_hpa", "local"),
    ])
    Set((String(row.product), String(row.variable_group), String(row.role))
        for row in eachrow(included)) == expected || error("Fixed joint specification differs")

    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    residual = filter(:method => in(["residual_gwr", "mixed_gwr", "mgwr"]), status)
    all(==("fixed_full_data"), residual.covariate_selection_mode) ||
        error("Residual models did not use the fixed specification")
    gpm_residual = filter([:product, :method] =>
        (product, method) -> product == "GPM" && method == "residual_gwr", residual)
    all(value -> occursin("elevation=local", value),
        gpm_residual.covariate_effective_roles) ||
        error("GPM elevation is not all-local in residual_gwr")
    gpm_mixed = filter([:product, :method] =>
        (product, method) -> product == "GPM" && method in ("mixed_gwr", "mgwr"), residual)
    all(value -> occursin("elevation=global", value),
        gpm_mixed.covariate_effective_roles) ||
        error("GPM elevation is not global in Mixed GWR/MGWR")

    metrics = CSV.read(joinpath(outdir, "metrics_pooled.csv"), DataFrame)
    Set(metrics.method) == Set(METHODS) || error("Pooled metrics do not contain eight methods")
    schemes = mode == :smoke ? ["balanced_spatial"] : ["balanced_spatial", "random"]
    for scheme in schemes, product in ["FY4B", "GPM", "GSMaP"]
        product_dir = joinpath(outdir, scheme, lowercase(product))
        for method in METHODS
            isfile(joinpath(product_dir, "oof_$(method).csv")) ||
                error("Missing OOF output: $scheme/$product/$method")
        end
        isfile(joinpath(product_dir, "common_evaluation_mask.csv")) ||
            error("Missing common mask: $scheme/$product")
    end

    scope = CSV.read(joinpath(outdir, "benchmark_scope.csv"), DataFrame)
    scope_values = Dict(String(row.key) => String(row.value) for row in eachrow(scope))
    get(scope_values, "covariate_selection_mode", "") == "fixed_full_data" ||
        error("Scope does not disclose fixed full-data selection")
    occursin("non-nested", get(scope_values, "covariate_selection_leakage", "")) ||
        error("Scope does not disclose non-nested selection")
    occursin("no paired significance", get(scope_values, "inference_policy", "")) ||
        error("Scope inference policy differs")
    println("Verified joint-covariate interpolation benchmark: mode=$mode, output=$outdir")
    return true
end

mode = isempty(ARGS) ? :smoke : Symbol(lowercase(ARGS[1]))
verify(mode)
