#!/usr/bin/env julia

"""
Does the permutation role test disagree with the MGWR bandwidth search, and if so, why?

Two mechanisms in this benchmark decide whether a covariate is treated as spatially varying.
`joint_spatial_variability_test` (`src/JointVariableSelection.jl`) assigns `role = "local"` or
`"global"` by permutation with a BH-FDR threshold, on a coefficient surface pooled over every wet
station-hour. `select_mgwr_bandwidths!` (`src/InterpolationBenchmarkTuning.jl`) then picks a
bandwidth per group by held-out RMSE, per hour, and `bw = Inf` makes a group operationally global
(`_local_hat` at `Inf` equals `_global_projection` to machine precision — pinned in
`test/test-dem-terrain-experiment.jl`).

They frequently disagree. This script quantifies the disagreement and then tests the four
mechanical explanations for it, because "the two criteria measure different things" is only
interesting once the boring causes are excluded:

  1. weak signal / FDR thresholding — are the disagreeing cells the marginally significant ones?
  2. tuning noise or ties       — did `Inf` win by a hair over a finite bandwidth?
  3. a flat objective           — is RMSE indistinguishable across the whole candidate set?
  4. the search grid or start   — does the coordinate descent simply never leave its `Inf` start?

    julia --project=. scripts/verify_role_vs_bandwidth.jl [RUN_DIR]

Read-only; RUN_DIR defaults to the canonical full nested run.

On the run this was written against: 64 included cells, 51 local / 13 global; 0 of 13 role-global
cells ever received a bandwidth (they enter the global design and never reach the search); 35 of
51 role-local cells chose `bw = Inf`. All four mechanical explanations were ruled out — 27 of the
35 sit at the permutation floor p = 0.001, no `Inf` win was within 0.01% of its runner-up, the
median objective spread was an order of magnitude above the winning margin, and the descent moved
groups off `Inf` in 16 of 51 cases while converging at iteration 2 in 74 of 81 rows.
"""

using CSV, DataFrames, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_RUN = joinpath(
    ROOT, "output", "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only",
)

"""Read `name` from the run directory, failing with the missing column list rather than a
`KeyError` thrown from somewhere deep in a comprehension."""
function read_output(run_dir::AbstractString, name::AbstractString, required::Vector{Symbol})
    path = joinpath(run_dir, name)
    isfile(path) || error("missing $(path) — is this a completed benchmark run directory?")
    frame = CSV.read(path, DataFrame)
    absent = filter(column -> !hasproperty(frame, column), required)
    isempty(absent) || error("$(path) has no column(s) $(join(absent, ", "))")
    return frame
end

"""The (scheme, product, fold, group) key a row belongs to, as comparable `String`s."""
_key(row, group_column::Symbol) =
    (String(row.scheme), String(row.product), row.fold, String(getproperty(row, group_column)))

function main(arguments::Vector{String})
    run_dir = isempty(arguments) ? DEFAULT_RUN : arguments[1]
    isdir(run_dir) || error("run directory not found: $(run_dir)")
    println("run ", run_dir)

    roles = read_output(run_dir, "joint_fold_roles.csv",
                        [:scheme, :product, :fold, :variable_group, :role, :final_included])
    variability = read_output(run_dir, "joint_fold_spatial_variability.csv",
                              [:scheme, :product, :fold, :variable_group, :pvalue, :qvalue])
    bandwidths = read_output(run_dir, "joint_bandwidths.csv",
                             [:scheme, :product, :fold, :method, :group, :bw, :kernel,
                              :adaptive, :iteration, :selected, :RMSE])
    scan = read_output(run_dir, "parameter_scan.csv",
                       [:scheme, :product, :fold, :method, :group, :bw, :kernel,
                        :adaptive, :iteration, :RMSE, :status])

    selected = filter(r -> r.selected === true && r.method == "mgwr", bandwidths)
    included = filter(:final_included => identity, roles)

    println("\n== What each mechanism decided ==")
    @printf("included covariate cells (scheme x product x fold x covariate): %d\n", nrow(included))
    @printf("  role = local  %d\n  role = global %d\n",
            count(==("local"), included.role), count(==("global"), included.role))

    # A role=global covariate goes into the global design block and is never offered to the
    # bandwidth search, so MGWR cannot promote it. Anything but 0 here means the design changed.
    with_bandwidth = Set(_key(r, :group) for r in eachrow(selected))
    global_cells = filter(:role => ==("global"), included)
    promoted = count(r -> _key(r, :variable_group) in with_bandwidth, eachrow(global_cells))
    @printf("\nrole = global cells MGWR gave a bandwidth to: %d of %d (0 expected — they enter\n",
            promoted, nrow(global_cells))
    println("  the global block, so MGWR can demote local -> global but never the reverse)")

    local_cells = filter(:role => ==("local"), included)
    rows = NamedTuple[]
    for r in eachrow(local_cells)
        key = _key(r, :variable_group)
        match = filter(x -> _key(x, :group) == key, selected)
        nrow(match) == 0 && continue
        test = filter(x -> _key(x, :variable_group) == key, variability)
        push!(rows, (; scheme=key[1], product=key[2], fold=key[3], group=key[4],
                     bw=match.bw[1], is_inf=!isfinite(match.bw[1]),
                     pvalue=nrow(test) == 1 ? test.pvalue[1] : NaN,
                     qvalue=nrow(test) == 1 ? test.qvalue[1] : NaN))
    end
    cross = DataFrame(rows)
    if isempty(cross)
        println("\nno role = local cell reached the bandwidth search; nothing further to report")
        return nothing
    end
    @printf("\nrole = local cells MGWR did give a bandwidth to: %d\n", nrow(cross))
    @printf("  bw = Inf  (test says local, prediction says global): %d (%.0f%%)\n",
            count(cross.is_inf), 100count(cross.is_inf) / nrow(cross))
    @printf("  bw finite (both say local):                         %d (%.0f%%)\n",
            count(.!cross.is_inf), 100count(.!cross.is_inf) / nrow(cross))

    println("\n== 1. Is the disagreement concentrated in marginal significance? ==")
    println("   (if it were, the disagreeing cells would sit near the q-threshold, not the floor)")
    for (label, subset) in (("bw = Inf ", filter(:is_inf => identity, cross)),
                            ("bw finite", filter(:is_inf => !, cross)))
        isempty(subset) && continue
        sorted = sort(subset.pvalue)
        @printf("   %s n=%2d  median p %.3f  max p %.3f  at permutation floor (p=0.001): %d\n",
                label, nrow(subset), sorted[cld(nrow(subset), 2)], maximum(subset.pvalue),
                count(≈(0.001), subset.pvalue))
    end

    println("\n== Per covariate ==")
    for group in sort(unique(cross.group))
        subset = filter(:group => ==(group), cross)
        finite_bw = sort(unique(filter(isfinite, subset.bw)))
        @printf("   %-12s n=%2d  Inf=%2d  finite=%2d  %s\n", group, nrow(subset),
                count(subset.is_inf), count(.!subset.is_inf),
                isempty(finite_bw) ? "" : "finite bws: " * join(finite_bw, ", "))
    end

    # How decisively did each selection win? A bandwidth chosen by a hair, on an objective that is
    # flat anyway, would make the whole comparison noise rather than signal.
    successes = filter(r -> r.method == "mgwr" && r.status == "success", scan)
    margins = NamedTuple[]
    for r in eachrow(filter(x -> x.group != "intercept", selected))
        candidates = filter(x -> _key(x, :group) == _key(r, :group) && x.kernel == r.kernel &&
                                 x.adaptive == r.adaptive && x.iteration == r.iteration, successes)
        nrow(candidates) < 2 && continue
        others = filter(x -> x.bw != r.bw, candidates)
        nrow(others) == 0 && continue
        best_other = minimum(others.RMSE)
        push!(margins, (; group=String(r.group), is_inf=!isfinite(r.bw),
                        margin=(best_other - r.RMSE) / r.RMSE,
                        spread=(maximum(candidates.RMSE) - minimum(candidates.RMSE)) /
                               minimum(candidates.RMSE)))
    end
    println("\n== 2-3. Tuning noise, ties, and how flat the objective is ==")
    if isempty(margins)
        println("   no selection had a scorable alternative in parameter_scan.csv")
    else
        table = DataFrame(margins)
        @printf("   selections with a scorable alternative: %d\n", nrow(table))
        for (label, subset) in (("chose Inf   ", filter(:is_inf => identity, table)),
                                ("chose finite", filter(:is_inf => !, table)))
            isempty(subset) && continue
            @printf("   %s n=%2d  median margin over next-best %.3g%%  median full spread %.3g%%\n",
                    label, nrow(subset), 100median(subset.margin), 100median(subset.spread))
        end
        ties = count(r -> r.is_inf && r.margin < 1e-4, eachrow(table))
        @printf("   Inf chosen with next-best within 0.01%% (a tie): %d of %d Inf selections\n",
                ties, count(table.is_inf))
        println("   (a margin well above 0 with a much larger spread means the winner was picked")
        println("    on a curved objective, not plucked from a flat one)")
    end

    # The descent starts every group at `last(family_candidates)`, which is Inf for the fixed
    # family. If it never moved, "chose Inf" would mean "never searched".
    println("\n== 4. Did the coordinate descent actually move off its Inf start? ==")
    descent = combine(groupby(filter(x -> x.group != "intercept", selected),
                              [:scheme, :product, :fold]),
                      :iteration => maximum => :final_iteration)
    settled = length(unique(collect(zip(String.(selected.scheme), String.(selected.product),
                                        selected.fold))))
    # Not every fold-cell has a covariate group: where the role test called everything global, the
    # only selected row is the intercept. `descent` therefore covers a subset of `settled`, and a
    # gap between the two is expected rather than a sign that a descent failed.
    @printf("   fold-cells that produced a winner: %d (with a covariate group: %d)\n",
            settled, nrow(descent))
    @printf("   final iteration over those: min %d, median %d, max %d\n",
            minimum(descent.final_iteration), round(Int, median(descent.final_iteration)),
            maximum(descent.final_iteration))
    @printf("   groups moved off Inf: %d of %d\n", count(.!cross.is_inf), nrow(cross))
    println("   (converging in 2 iterations is one sweep scoring every candidate plus one")
    println("    confirming sweep — the search ran; hitting mgwr_max_tuning_iterations would")
    println("    instead mean candidates were dropped unscored, which is finding F8)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
