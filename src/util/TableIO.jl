"""
Table plumbing shared across tiers.

Two helpers that had drifted into five and two near-identical copies. Neither is domain logic:
`write_csv_atomic` is the write-temp -> `mv(force=true)` -> clean-up dance the `prepare_*` pipelines
depend on for rewriting `data/processed/` safely, and `column` is the case-insensitive header
lookup every station reader needs.

It has its own tier because its callers span `sources/` and `mger/`, so it belongs to neither. Keep
it that way: this is for plumbing with no opinion about precipitation, not a drawer for anything
that happens to be shared.
"""
module TableIO

using CSV

export write_csv_atomic

"""
Write `table` to `path` through a temporary file, so a reader never sees a partial CSV.

`mkpath`s the parent, writes `<path>.tmp-<pid>`, then moves it into place with `force=true`. The
`finally` removes the temporary when the write or the move threw - which four of the five copies
this replaces did, and `FY4BPreprocessing`'s did not: that one left a `.tmp-<pid>` file behind on
failure. Successful writes behave identically to every copy it replaces.

Returns `path`.
"""
function write_csv_atomic(path::AbstractString, table)
    mkpath(dirname(path))
    temporary = string(path, ".tmp-", getpid())
    try
        CSV.write(temporary, table)
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

"""
The column of `df` matching `requested` ignoring case, returned as the name `df` actually uses.

`label` names the thing in the error message, so a reader with its own wording keeps it. `df` is
untyped so this module needs no DataFrames dependency; anything `names` works on will do.
"""
function column(df, requested::Symbol; label::AbstractString="column")
    mapping = Dict(Symbol(lowercase(String(name))) => name for name in names(df))
    found = get(mapping, Symbol(lowercase(String(requested))), nothing)
    found === nothing && throw(ArgumentError("missing $label: $requested"))
    return found
end

end # module
