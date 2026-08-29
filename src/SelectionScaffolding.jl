"""
Bookkeeping shared by the four variable-selection paths.

`ERA5VariableSelection`, `NDVIVariableSelection` and `JointVariableSelection` each carried their
own byte-identical copy of the two DataFrame accumulation helpers below, and each rebuilt the
`(scheme, fold, train_indices)` list with the same two lines. Nothing here is specific to a
covariate family: it is the scaffolding those searches hang on, so a change to how diagnostic
tables are stamped or how the CV schemes are enumerated is now one edit rather than three.

Deliberately narrow. The searches themselves - screening, VIF pruning, the spatial role tests -
look alike across the families but are not the same computation (weighted vs unweighted VIF,
different fold builders), and are left where they are.
"""
module SelectionScaffolding

using DataFrames

export annotate_selection!, append_selection!, selection_schemes

"""
Stamp `product`, `scheme` and `fold` onto the front of a diagnostic table, in that column order.

Every table written by a variable-selection run carries this triple, so a row can be traced to the
product and fold that produced it.
"""
function annotate_selection!(table::DataFrame, product::String, scheme::String, fold::Int)
    insertcols!(table, 1, :product => fill(product, nrow(table)),
        :scheme => fill(scheme, nrow(table)), :fold => fill(fold, nrow(table)))
    return table
end

"""
Append `source` onto `target`, seeding `target`'s columns from `source` when it is still empty.

The accumulators start as bare `DataFrame()`s with no schema, so the first append has to establish
one; `cols=:union` then tolerates a later table carrying an extra column.
"""
function append_selection!(target::DataFrame, source::DataFrame)
    if ncol(target) == 0
        for column in propertynames(source)
            target[!, column] = similar(source[!, column], 0)
        end
    end
    nrow(source) > 0 && append!(target, source; cols=:union)
    return target
end

"""
The `(scheme, fold, train_indices)` list a variable-selection run walks.

One `"full_data"` pass over all `n` stations, then one `"spatial_cv"` pass per fold holding that
fold out. `folds` is a per-station fold assignment; `k` is the fold count.
"""
function selection_schemes(n::Int, folds::AbstractVector{<:Integer}, k::Int)
    schemes = [("full_data", 0, collect(1:n))]
    append!(schemes, [("spatial_cv", fold, findall(!=(fold), folds)) for fold in 1:k])
    return schemes
end

end # module SelectionScaffolding
