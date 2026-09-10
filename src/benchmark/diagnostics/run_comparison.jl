# ------------------------------------------------------- D6: comparing two benchmark runs

"""
Score two runs of the benchmark against each other on cells they both evaluated.

A raw RMSE delta between two run directories is meaningless here. The shared evaluation mask
keeps only cells where every mask-defining method is finite, so any change that alters one
method's coverage resizes the denominator for all of them — and the cells that appear or vanish
are systematically the hard ones. Two masks are therefore intersected before anything is scored,
and each method is additionally restricted to cells where *both* runs produced a finite value,
which makes `delta_paired` a genuine paired difference rather than two numbers from two samples.

`RMSE_before` / `RMSE_after` are each run's own published-style number (its own mask, its own
NaNs dropped) and are reported alongside so the size of the mask effect is visible rather than
hidden. Trust `delta_paired`; read the own-mask columns only to see how much the mask moved.
"""
function run_comparison_table(
    y_obs::Matrix{Float64},
    before::AbstractDict{String,Matrix{Float64}}, before_mask::AbstractMatrix,
    after::AbstractDict{String,Matrix{Float64}}, after_mask::AbstractMatrix;
    scheme::String, product::String,
)
    shared = BitMatrix(before_mask .& after_mask)
    rows = NamedTuple[]
    for method in sort(collect(union(keys(before), keys(after))))
        in_before = haskey(before, method)
        in_after = haskey(after, method)
        own_before = in_before ? _metrics(y_obs, before[method], before_mask) : nothing
        own_after = in_after ? _metrics(y_obs, after[method], after_mask) : nothing
        paired_before = (; n=0, RMSE=NaN)
        paired_after = (; n=0, RMSE=NaN)
        if in_before && in_after
            pair = BitMatrix(shared .& .!isnan.(before[method]) .& .!isnan.(after[method]))
            paired_before = _metrics(y_obs, before[method], pair)
            paired_after = _metrics(y_obs, after[method], pair)
        end
        delta = paired_after.RMSE - paired_before.RMSE
        push!(rows, (;
            scheme, product, method,
            present_in=in_before && in_after ? "both" : (in_before ? "before" : "after"),
            mask_cells_before=count(before_mask), mask_cells_after=count(after_mask),
            mask_cells_shared=count(shared),
            n_before=own_before === nothing ? 0 : own_before.n,
            n_after=own_after === nothing ? 0 : own_after.n,
            RMSE_before=own_before === nothing ? NaN : own_before.RMSE,
            RMSE_after=own_after === nothing ? NaN : own_after.RMSE,
            n_paired=paired_before.n,
            RMSE_paired_before=paired_before.RMSE,
            RMSE_paired_after=paired_after.RMSE,
            delta_paired=delta,
            relative_paired=isnan(delta) || paired_before.RMSE == 0 ? NaN :
                delta / paired_before.RMSE,
        ))
    end
    return DataFrame(rows)
end

