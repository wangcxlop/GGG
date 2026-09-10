# --------------------------------------------------- D2a: is the bandwidth grid binding?

"""
Summarise a `parameter_scan.csv` / `joint_bandwidths.csv` scan: for every selected
candidate, whether it sits on an endpoint of the grid it was chosen from, and whether the
CV curve was monotone up to that endpoint (which is what says the grid, not the data,
picked the value).
"""
function bandwidth_saturation_table(scan::DataFrame)
    rows = NamedTuple[]
    keys_of_interest = [:scheme, :product, :fold, :mode, :method, :group, :kernel, :adaptive]
    for subgroup in groupby(unique(scan), keys_of_interest)
        # `!isnan` rather than `isfinite`: the filter's job is excluding the `NaN` bw of the
        # idw/adw/tps rows, and the GWR family's explicit global candidate is `bw = Inf`, which
        # sorts last and so still reads as "chose the widest candidate offered".
        candidates = filter(row -> row.status == "success" && !isnan(row.bw), subgroup)
        nrow(candidates) < 2 && continue
        selected = filter(row -> row.selected === true, candidates)
        nrow(selected) == 1 || continue
        # MGWR tunes each group by coordinate descent, so a group spans several sweeps over the
        # same grid. Only the sweep that produced the winner is a comparable candidate set.
        candidates = filter(row -> row.iteration == selected.iteration[1], candidates)
        nrow(candidates) < 2 && continue
        order = sortperm(candidates.bw)
        widths = candidates.bw[order]
        errors = candidates.RMSE[order]
        chosen = selected.bw[1]
        at_min = chosen == first(widths)
        at_max = chosen == last(widths)
        # Monotone toward the chosen endpoint means the search was clipped, not resolved.
        increasing = all(diff(errors) .> 0)
        decreasing = all(diff(errors) .< 0)
        # Landing on the widest candidate only counts as clipping if a wider one could exist.
        # The GWR family's top candidate is the explicit global fit (`bw = Inf`), which is the
        # end of the bandwidth continuum, not the end of an arbitrary grid - selecting it is a
        # resolved answer ("this coefficient wants to be global"), so it must not be reported
        # as a grid that needs extending.
        extensible_max = isfinite(last(widths))
        push!(rows, (;
            scheme=subgroup[1, :scheme], product=subgroup[1, :product],
            fold=subgroup[1, :fold], mode=subgroup[1, :mode],
            method=subgroup[1, :method], group=subgroup[1, :group],
            kernel=subgroup[1, :kernel], adaptive=subgroup[1, :adaptive],
            n_candidates=nrow(candidates), bw_min=first(widths), bw_max=last(widths),
            bw_selected=chosen, at_grid_min=at_min, at_grid_max=at_max,
            monotone_increasing=increasing, monotone_decreasing=decreasing,
            clipped_by_grid=(at_min && increasing) || (at_max && decreasing && extensible_max),
            RMSE_selected=selected.RMSE[1],
            RMSE_at_min=first(errors), RMSE_at_max=last(errors),
        ))
    end
    return DataFrame(rows)
end

