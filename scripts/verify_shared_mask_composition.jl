#!/usr/bin/env julia

"""
What the shared evaluation mask costs, and who is responsible for the cost.

A cell is scored only when every method in `MASK_METHODS` produced a finite value there, so the
mask is an intersection over eight methods. `BenchmarkDiagnostics.mask_cost_table` measures each
method's contribution by leaving it out one at a time — which is structurally blind to a
*correlated* failure: when the three joint methods go NaN on the same cells, dropping any one of
them recovers nothing, and the table reports zero for all three. This script reports the
leave-one-out numbers alongside the leave-the-group-out number, so a shared failure cannot hide.

It also asks whether the lost cells are a random sample. A mask that quietly removes the wettest
hours would bias every RMSE computed on it, so kept and lost cells are compared on wet fraction
and intensity, and the loss is checked for concentration in whole hours or whole stations.

    julia --project=. scripts/verify_shared_mask_composition.jl [RUN_DIR] [OBS_CSV]

Both arguments are optional and default to the canonical full nested run and the hourly gauge
file. Read-only: it reads a completed run directory and writes nothing.

Findings from the run this was written against (full, nested covariates, mgwr :intercept_only):
under `balanced_spatial` the mask keeps ~99.85% of evaluable cells and essentially 100% of the
loss is the joint family failing together — leaving any single method out recovers 0 cells for
GPM and GSMaP, while leaving the group out recovers all ~4000. The loss is concentrated in a few
stations (9 for GPM, 57 for FY4B), the worst losing 19% of its hours, and never removes a whole
hour. Lost cells are slightly *drier* than kept ones, so the mask does not flatter the methods.
Under `random` nothing is lost at all, which is itself the point: the cost is a property of the
spatial fold geometry, not of the methods.
"""

using CSV, DataFrames, Dates, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_RUN = joinpath(
    ROOT, "output", "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only",
)
const DEFAULT_OBS = joinpath(
    ROOT, "data", "processed", "study_area", "hubei_obs_hourly_2022_2025_JunSep.csv",
)

# Must match `MASK_METHODS` in src/InterpolationBenchmarkConfig.jl. Hard-coded rather than loaded
# because this script is a plain CSV reader and pulling in the benchmark would mean loading the
# whole package. Every name is checked against the run directory below, so drift shows up as a
# missing-file error rather than a silently different mask.
const MASK_METHODS = ["raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr"]
# `BenchmarkDiagnostics.JOINT_MASK_METHODS` — the three that share a back-fit and fail together.
const JOINT_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]

"""Normalise a time value from either CSV to a `DateTime`, tolerating a trailing `Z`."""
_timestamp(value) = value isa DateTime ? value : DateTime(rstrip(String(value), 'Z'))

"""Every (scheme, product) pair under `run_dir` that has a written evaluation mask."""
function mask_directories(run_dir::AbstractString)
    pairs = Tuple{String,String}[]
    for scheme in sort(readdir(run_dir))
        scheme_dir = joinpath(run_dir, scheme)
        isdir(scheme_dir) || continue
        for product in sort(readdir(scheme_dir))
            isfile(joinpath(scheme_dir, product, "common_evaluation_mask.csv")) &&
                push!(pairs, (scheme, product))
        end
    end
    return pairs
end

"""The out-of-fold prediction matrix for `method`, as hours × stations with NaN for missing."""
function prediction_matrix(dir::AbstractString, method::AbstractString, stations::Vector{String})
    path = joinpath(dir, "oof_$(method).csv")
    isfile(path) || error("""
        missing $(path).
        The mask is an intersection over MASK_METHODS; without every one of them this script
        would report a different mask than the run actually used. If MASK_METHODS has changed,
        update the copy at the top of this file.""")
    frame = CSV.read(path, DataFrame; select=vcat("time", stations))
    return Matrix{Float64}(coalesce.(Matrix(frame[:, 2:end]), NaN))
end

"""Wet fraction, mean and upper tail of a set of gauge values, for kept-vs-lost comparison."""
function intensity(values::Vector{Float64})
    isempty(values) && return (; n=0, wet=NaN, mean=NaN, p99=NaN, max=NaN)
    return (; n=length(values), wet=100mean(values .>= 0.1), mean=mean(values),
            p99=quantile(values, 0.99), max=maximum(values))
end

function report_cell(run_dir::AbstractString, scheme::AbstractString, product::AbstractString,
                     observations::DataFrame, row_of::Dict{DateTime,Int})
    dir = joinpath(run_dir, scheme, product)
    mask_frame = CSV.read(joinpath(dir, "common_evaluation_mask.csv"), DataFrame)
    stations = String.(names(mask_frame))[2:end]
    rows = [row_of[_timestamp(t)] for t in mask_frame.time]
    n_time, n_station = length(rows), length(stations)

    gauge = Matrix{Float64}([
        coalesce(observations[rows[i], station], NaN) for i in 1:n_time, station in stations
    ])
    finite = Dict(m => .!isnan.(prediction_matrix(dir, m, stations)) for m in MASK_METHODS)

    observed = .!isnan.(gauge)
    # "Evaluable" is the honest denominator: a cell nobody could ever score — no gauge value, or
    # no satellite value to correct — is not a cost the mask imposes.
    evaluable = observed .& finite["raw"]
    shared = copy(observed)
    for m in MASK_METHODS
        shared .&= finite[m]
    end
    n_evaluable, n_shared = count(evaluable), count(shared)

    @printf("\n%s / %s\n", scheme, uppercase(product))
    @printf("  grid %d x %d = %d cells; evaluable (gauge + satellite) %d; shared mask %d (%.5f)\n",
            n_time, n_station, n_time * n_station, n_evaluable, n_shared, n_shared / n_evaluable)

    @printf("\n  Leave one method out — cells the mask would recover:\n")
    for m in MASK_METHODS
        reduced = copy(observed)
        for other in MASK_METHODS
            other == m || (reduced .&= finite[other])
        end
        recovered = count(reduced) - n_shared
        @printf("    without %-13s %7d cells (%.4f%% of evaluable)\n",
                m, recovered, 100recovered / n_evaluable)
    end
    # The number one-at-a-time cannot see.
    grouped = copy(observed)
    for m in MASK_METHODS
        m in JOINT_METHODS || (grouped .&= finite[m])
    end
    group_recovered = count(grouped) - n_shared
    @printf("    without %-13s %7d cells (%.4f%% of evaluable)  <- the correlated failure\n",
            join(JOINT_METHODS, "+"), group_recovered, 100group_recovered / n_evaluable)

    lost = evaluable .& .!shared
    n_lost = count(lost)
    @printf("\n  Evaluable but outside the shared mask: %d (%.4f%%)\n",
            n_lost, 100n_lost / n_evaluable)
    if n_lost == 0
        println("    nothing lost; the remaining diagnostics do not apply")
        return nothing
    end
    for m in MASK_METHODS
        responsible = count(lost .& .!finite[m])
        responsible == 0 || @printf("    %-13s is NaN on %7d of them (%.1f%%)\n",
                                    m, responsible, 100responsible / n_lost)
    end

    kept_values, lost_values = gauge[shared], gauge[lost]
    k, l = intensity(kept_values), intensity(lost_values)
    @printf("\n  Are the lost cells a random sample of the kept ones?\n")
    @printf("    kept %8d cells  wet(>=0.1mm) %6.2f%%  mean %.4f  p99 %6.3f  max %5.1f\n",
            k.n, k.wet, k.mean, k.p99, k.max)
    @printf("    lost %8d cells  wet(>=0.1mm) %6.2f%%  mean %.4f  p99 %6.3f  max %5.1f\n",
            l.n, l.wet, l.mean, l.p99, l.max)

    lost_per_hour = vec(sum(lost, dims=2))
    lost_per_station = vec(sum(lost, dims=1))
    evaluable_per_hour = vec(sum(evaluable, dims=2))
    entire = count(i -> lost_per_hour[i] > 0 && lost_per_hour[i] == evaluable_per_hour[i], 1:n_time)
    @printf("\n  Is the loss concentrated?\n")
    @printf("    hours with any loss    %5d of %5d; hours lost entirely %d\n",
            count(>(0), lost_per_hour), n_time, entire)
    @printf("    stations with any loss %5d of %5d; worst loses %d of %d hours (%.0f%%)\n",
            count(>(0), lost_per_station), n_station, maximum(lost_per_station), n_time,
            100maximum(lost_per_station) / n_time)
    return nothing
end

function main(arguments::Vector{String})
    run_dir = isempty(arguments) ? DEFAULT_RUN : arguments[1]
    obs_path = length(arguments) >= 2 ? arguments[2] : DEFAULT_OBS
    isdir(run_dir) || error("run directory not found: $(run_dir)")
    isfile(obs_path) || error("gauge file not found: $(obs_path)")
    cells = mask_directories(run_dir)
    isempty(cells) && error(
        "no scheme/product directory under $(run_dir) contains common_evaluation_mask.csv"
    )

    observations = CSV.read(obs_path, DataFrame)
    row_of = Dict(_timestamp(observations.time[i]) => i for i in 1:nrow(observations))

    println("run   ", run_dir)
    println("gauge ", obs_path)
    println("mask  ", join(MASK_METHODS, ", "))
    for (scheme, product) in cells
        report_cell(run_dir, scheme, product, observations, row_of)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
