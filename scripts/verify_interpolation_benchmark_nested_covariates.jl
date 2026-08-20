#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]

function verify(mode::Symbol)
    mode in (:smoke, :full) || error("mode must be smoke or full")
    outdir = joinpath(ROOT, "output", (mode == :smoke ?
        "interpolation_benchmark_smoke_joint_covariates" :
        "interpolation_benchmark_full_joint_covariates") * "_nested")
    required = [
        "joint_spec_provenance.csv", "joint_era5_input_qc.csv", "joint_fold_scaling.csv",
        "joint_fold_quality_control.csv", "joint_bandwidths.csv", "covariate_model_status.csv",
        "joint_fold_candidates.csv", "joint_fold_vif.csv", "joint_fold_spatial_variability.csv",
        "joint_fold_roles.csv", "joint_role_stability.csv",
        "metrics_stratified.csv", "metrics_folds.csv", "metrics_pooled.csv",
        "parameter_scan.csv", "run_status.csv", "benchmark_scope.csv",
        "paired_comparisons.csv", "claim_assessment.csv",
    ]
    missing_files = filter(file -> !isfile(joinpath(outdir, file)), required)
    isempty(missing_files) || error("Missing nested joint benchmark outputs: $(join(missing_files, ", "))")
    # The fixed-mode path's single spec file has no nested-mode equivalent — a fresh spec is
    # computed per fold instead, which is exactly what joint_fold_roles.csv records.
    !isfile(joinpath(outdir, "joint_variable_spec_used.csv")) ||
        error("Nested benchmark should not produce a single fixed variable specification")

    provenance = CSV.read(joinpath(outdir, "joint_spec_provenance.csv"), DataFrame)
    only(provenance.selection_mode) == "nested_per_fold" ||
        error("Provenance does not disclose nested per-fold selection")

    status = CSV.read(joinpath(outdir, "run_status.csv"), DataFrame)
    residual = filter(:method => in(["residual_gwr", "mixed_gwr", "mgwr"]), status)
    all(==("nested_per_fold"), residual.covariate_selection_mode) ||
        error("Residual models did not use nested per-fold selection")

    roles = CSV.read(joinpath(outdir, "joint_fold_roles.csv"), DataFrame)
    nrow(roles) > 0 || error("joint_fold_roles.csv is empty")
    all(in(("scheme", "product", "fold", "repeat", "seed")), string.(names(roles))[1:5]) ||
        error("joint_fold_roles.csv is missing repeat/seed provenance columns")
    all(value -> value in ("local", "global", "not_selected"), roles.role) ||
        error("joint_fold_roles.csv contains an unresolved role")
    # Selection must vary the specification per fold to be a genuine fix for the full-data leak
    # (a nested run that always reselects the identical set would be indistinguishable from the
    # fixed-mode leak it is meant to close).
    per_cell_specs = Dict{Tuple{String,String},Set{String}}()
    for group in groupby(filter(:final_included => identity, roles), [:scheme, :product])
        key = (String(group.scheme[1]), String(group.product[1]))
        per_cell_specs[key] = Set(
            "$(row.fold)|$(row.variable_group)|$(row.role)" for row in eachrow(group)
        )
    end
    stability = CSV.read(joinpath(outdir, "joint_role_stability.csv"), DataFrame)
    nrow(stability) > 0 || error("joint_role_stability.csv is empty")
    any(0 .< stability.selection_frequency .< 1) ||
        @warn "Every covariate was selected in either all folds or none — check this is expected " *
            "for the fixture in use, since real data is not expected to be this stable."

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

    bootstrap = CSV.read(joinpath(outdir, "paired_comparisons.csv"), DataFrame)
    nrow(bootstrap) > 0 ||
        error("Nested benchmark should re-enable the paired significance test")
    claim = CSV.read(joinpath(outdir, "claim_assessment.csv"), DataFrame)
    nrow(claim) > 0 || error("Nested benchmark should produce a claim assessment")

    scope = CSV.read(joinpath(outdir, "benchmark_scope.csv"), DataFrame)
    scope_values = Dict(String(row.key) => String(row.value) for row in eachrow(scope))
    get(scope_values, "covariate_selection_mode", "") == "nested_per_fold" ||
        error("Scope does not disclose nested per-fold selection")
    occursin("re-run inside every training fold", get(scope_values, "covariate_selection_leakage", "")) ||
        error("Scope does not disclose the nested selection leakage control")
    !occursin("no paired significance", get(scope_values, "inference_policy", "")) ||
        error("Scope inference policy still disclaims significance testing")
    println("Verified nested-covariate interpolation benchmark: mode=$mode, output=$outdir")
    return true
end

mode = isempty(ARGS) ? :smoke : Symbol(lowercase(ARGS[1]))
verify(mode)
