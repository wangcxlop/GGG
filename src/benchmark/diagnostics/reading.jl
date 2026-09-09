# ---------------------------------------------------------------------------- reading

"""
Where a diagnosis of `run_dir` is written: `output/benchmark_diagnostics/<run name>`.

Each diagnostic run gets its own subdirectory named for the benchmark run it read, so diagnosing a
second run never overwrites the first one's numbers. Repeated verbatim in
`run_benchmark_diagnostics.jl`, `run_mgwr_diagnostics.jl` and `run_claim_reassessment.jl`, which
also each checked `run_dir` exists first - done here so the error is worded once.
"""
function benchmark_diagnostics_outdir(root::AbstractString, run_dir::AbstractString)
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    outdir = joinpath(root, "output", "benchmark_diagnostics", basename(run_dir))
    mkpath(outdir)
    return outdir
end


"""Read a wide `time × station` benchmark CSV into a `station × time` matrix ordered by `ids`."""
function read_wide_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = Matrix{Float64}(undef, length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        column = df[!, Symbol(id)]
        out[index, :] = [value === missing ? NaN : Float64(value) for value in column]
    end
    return out
end

"""Read `common_evaluation_mask.csv` into a `station × time` `BitMatrix` ordered by `ids`."""
function read_mask_matrix(path::AbstractString, ids::Vector{String})
    df = CSV.read(path, DataFrame)
    out = falses(length(ids), nrow(df))
    for (index, id) in enumerate(ids)
        out[index, :] = Bool.(df[!, Symbol(id)])
    end
    return out
end

"""Read every `oof_<method>.csv` in one `(scheme, product)` directory into `method => matrix`."""
function load_prediction_matrices(product_dir::AbstractString, ids::Vector{String})
    predictions = Dict{String,Matrix{Float64}}()
    for path in readdir(product_dir; join=true)
        name = basename(path)
        startswith(name, "oof_") && endswith(name, ".csv") || continue
        predictions[name[5:end-4]] = read_wide_matrix(path, ids)
    end
    return predictions
end

"""Read `split_common.csv` into `station_id => fold`."""
function read_fold_map(path::AbstractString)
    df = CSV.read(path, DataFrame; types=Dict(:station_id => String))
    return Dict(String(row.station_id) => Int(row.fold) for row in eachrow(df))
end

