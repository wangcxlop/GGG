# ------------------------------------------------- D4: do the covariates buy anything?

"""
Join each fold's selected-covariate count onto that fold's RMSE, so "more covariates" can
be checked against "better fold". Folds that selected nothing run intercept-only and are the
natural control.
"""
function covariate_contribution_table(status::DataFrame, folds::DataFrame; level::String="all")
    metrics = filter(row -> row.group == "overall" && row.level == level, folds)
    joined = innerjoin(
        select(status, [:scheme, :product, :fold, :method, :covariate_variable_count,
            :covariate_variables, :covariate_effective_roles, :prediction_coverage]),
        select(metrics, [:scheme, :product, :fold, :method, :RMSE, :MAE, :Bias, :r, :n]),
        on=[:scheme, :product, :fold, :method],
    )
    baseline = filter(row -> row.method == "adw", metrics)
    rename!(baseline, :RMSE => :RMSE_adw)
    joined = leftjoin(
        joined, select(baseline, [:scheme, :product, :fold, :RMSE_adw]),
        on=[:scheme, :product, :fold],
    )
    joined.rmse_gap_vs_adw = joined.RMSE .- joined.RMSE_adw
    return sort!(joined, [:scheme, :product, :method, :fold])
end
