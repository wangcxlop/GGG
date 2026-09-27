#!/usr/bin/env julia

# Can blend_mgwr on MERGED_OLS_LAGNBR be improved by a better combiner over predictions the benchmark
# already made? A screen, not a result: it fits no GWR and reruns nothing.
#
#   julia -t 4 --project=. scripts/screen_blend_mgwr_improvements.jl
#   julia -t 4 --project=. scripts/screen_blend_mgwr_improvements.jl <benchmark_output_dir>
#
# Reads the stored out-of-fold predictions of a finished `--fused-anchor-lagnbr` run and, per scheme,
# refits each candidate combiner leave-one-fold-out: weights are fitted on the other folds' cells and
# scored on the held-out fold. Those other folds' predictions came from models that saw the held-out
# gauges, so the screen leaks slightly - through a handful of weights, never through a prediction -
# and every candidate, including `C0`/`C1` (the benchmark's own blends refitted the same way),
# enjoys the same leak. Compare candidates against `C0`/`C1`, not against the saved `blend_mgwr`,
# when judging what the combiner itself adds.
#
# Writes `candidates.csv` and `weights.csv` to output/blend_mgwr_screen/.

const ROOT = normpath(joinpath(@__DIR__, ".."))

using MixedGWR
using CSV, DataFrames, Dates, LinearAlgebra, Random, Statistics

include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")
load_standalone_modules("BenchmarkDiagnostics", "TableIO")
using Main.BenchmarkDiagnostics
using Main.TableIO: write_csv_atomic

const STUDY_DATA = joinpath(ROOT, "data", "processed", "study_area")
const DEFAULT_RUN = joinpath(ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_" *
    "satwetblend_blendagrenv_fusedanchor_fusedlagnbr")
const OUTDIR = joinpath(ROOT, "output", "blend_mgwr_screen")
const SCHEMES = ["balanced_spatial", "random"]
const PRODUCT = "MERGED_OLS_LAGNBR"
const SOURCE_PRODUCTS = ["FY4B", "GPM", "GSMaP"]
const THRESHOLD = 0.1
const SEED = 20260927
const BOOTSTRAP_REPS = 1000
# A band with fewer training cells than this takes the pooled fit instead of its own.
const MIN_BAND_CELLS = 500
const DISTANCE_EDGES = [20.0, 50.0]
const STRATA = [("all", -Inf, Inf), ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
    ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf)]

"""The `MGERConfig` the full benchmark run used, so the common station/time grid matches."""
function study_config(outdir::AbstractString)
    return MGERConfig(
        station_meta_path=joinpath(STUDY_DATA, "station_meta.csv"),
        obs_hourly_wide_path=joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        sat_paths=Dict(
            "FY4B" => joinpath(
                STUDY_DATA, "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
            ),
            "GPM" => joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2024_full_aligned.csv"),
            "GSMaP" => joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2024_full_aligned.csv"),
        ),
        outdir=outdir,
        rain_threshold=THRESHOLD,
        analysis_start=DateTime(2022, 1, 1, 9),
        analysis_end=DateTime(2025, 1, 1, 8),
        expected_common_time_count=13471,
    )
end

"""
`M` shifted by one hour, `step = -1` for t-1 and `+1` for t+1. A cell whose neighbouring column is
not exactly one hour away, or is NaN, keeps its own value, so the feature degrades to the member.
"""
function hour_shift(M::Matrix{Float64}, times::Vector{DateTime}, step::Int)
    out = copy(M)
    for t in axes(M, 2)
        u = t + step
        1 <= u <= size(M, 2) && abs(times[u] - times[t]) == Hour(1) || continue
        for s in axes(M, 1)
            isnan(M[s, u]) || (out[s, t] = M[s, u])
        end
    end
    return out
end

"""
Every masked, fully predicted cell of one scheme as a column table, plus what the fits band on.
"""
function cell_table(scheme_dir, grid, ids, lonlat, sources)
    product_dir = joinpath(scheme_dir, lowercase(PRODUCT))
    oof(method) = read_wide_matrix(joinpath(product_dir, "oof_$(method).csv"), ids)
    members = Dict(m => oof(m) for m in
        ("mgwr", "adw", "raw", "tps", "idw", "blend_mgwr", "blend_agrenv_mgwr"))
    mask = read_mask_matrix(joinpath(product_dir, "common_evaluation_mask.csv"), ids)
    mask .&= .!isnan.(grid.Y_obs)
    for M in values(members)
        mask .&= .!isnan.(M)
    end
    lags = Dict(
        "mgwr_lag" => hour_shift(members["mgwr"], grid.times, -1),
        "mgwr_lead" => hour_shift(members["mgwr"], grid.times, 1),
        "adw_lag" => hour_shift(members["adw"], grid.times, -1),
        "adw_lead" => hour_shift(members["adw"], grid.times, 1),
    )
    agreement = blend_band_matrix(:agreement_envelope, first(sources), sources, THRESHOLD)

    fold_map = read_fold_map(joinpath(scheme_dir, "split_common.csv"))
    fold_of = [fold_map[id] for id in ids]
    distance = nearest_train_km(fold_of, lonlat)
    distance_band = [1 + count(e -> d >= e, DISTANCE_EDGES) for d in distance]

    cells = findall(mask)
    station = [c[1] for c in cells]
    time = [c[2] for c in cells]
    days = Date.(grid.times)
    unique_days = sort(unique(days))
    day_index = Dict(d => i for (i, d) in enumerate(unique_days))
    source_max = reduce((a, b) -> max.(a, b), (replace(S, NaN => 0.0) for S in sources))

    table = Dict{String,Vector{Float64}}(
        "y" => grid.Y_obs[cells],
        (m => M[cells] for (m, M) in members)...,
        (m => M[cells] for (m, M) in lags)...,
        "gmax" => source_max[cells],
        "one" => ones(length(cells)),
    )
    meta = (;
        station, time, cells, size=size(mask),
        fold=fold_of[station],
        wet=members["raw"][cells] .>= THRESHOLD,
        agreement=agreement[cells],
        distance=distance_band[station],
        day=[day_index[days[t]] for t in time],
        n_days=length(unique_days),
    )
    return table, meta
end

# ---------------------------------------------------------------------------------------------
# Bands. `0` is a band no candidate touches: the cell keeps the source prediction.

band_of(kind::Symbol, meta) =
    kind === :wet ? Int.(meta.wet) :
    kind === :agr ? meta.agreement :
    kind === :wet_dist ? [w ? d : 0 for (w, d) in zip(meta.wet, meta.distance)] :
    kind === :agr_dist ? [a == 0 ? 0 : a + 9 * (d - 1) for (a, d) in zip(meta.agreement, meta.distance)] :
    kind === :wetdry ? [w ? 1 : 2 for w in meta.wet] :
    kind === :agr0 ? meta.agreement .+ 1 :
    kind === :agr0_dist ? [a + 1 + 10 * (d - 1) for (a, d) in zip(meta.agreement, meta.distance)] :
    kind === :top ? [a > 0 && a % 3 == 0 ? a : 0 for a in meta.agreement] :
    throw(ArgumentError("unknown band kind $kind"))

"""Per `(fold, band)` Gram matrix of `[features..., y]`."""
function grams(table, features, band, fold, n_folds, n_bands)
    Z = hcat((table[f] for f in features)..., table["y"])
    k = size(Z, 2)
    G = [zeros(k, k) for _ in 1:n_folds, _ in 0:n_bands]
    n = zeros(Int, n_folds, n_bands + 1)
    for i in axes(Z, 1)
        b = band[i]
        b == 0 && continue
        z = view(Z, i, :)
        BLAS.ger!(1.0, z, z, G[fold[i], b + 1])
        n[fold[i], b + 1] += 1
    end
    return G, n
end

sse(G, w) = (k = size(G, 1); G[k, k] - 2 * dot(w, G[1:k-1, k]) + dot(w, G[1:k-1, 1:k-1] * w))

function ols(G, idx)
    k = size(G, 1)
    A = G[idx, idx]
    A += 1e-9 * tr(A) / length(idx) * I
    w = zeros(k - 1)
    w[idx] = A \ G[idx, k]
    return w
end

"""Non-negative least squares by subset enumeration; fine for the four or so members used here."""
function nnls(G, idx)
    best, best_w = Inf, zeros(size(G, 1) - 1)
    for mask in 1:(2^length(idx) - 1)
        subset = [idx[j] for j in eachindex(idx) if (mask >> (j - 1)) & 1 == 1]
        w = ols(G, subset)
        all(>=(0), w) || continue
        s = sse(G, w)
        s < best && ((best, best_w) = (s, w))
    end
    return best_w
end

"""Convex blend `(1 - λ) source + λ fallback`, λ by closed form clamped to [0, 1]."""
function convex(G, source, fallback)
    k = size(G, 1)
    dd = G[fallback, fallback] - 2 * G[source, fallback] + G[source, source]
    dr = G[fallback, k] - G[source, fallback] - G[source, k] + G[source, source]
    lambda = dd > 0 ? clamp(dr / dd, 0.0, 1.0) : 0.0
    w = zeros(k - 1)
    w[source] = 1 - lambda
    w[fallback] += lambda
    return w
end

"""
One candidate: which features, how the cells are banded, and how each band is fitted.
`fit(G, band)` returns a weight vector over `features`.
"""
struct Candidate
    name::String
    family::String
    features::Vector{String}
    band::Symbol
    fit::Function
end

source_weights(features) = (w = zeros(length(features)); w[findfirst(==("mgwr"), features)] = 1; w)

function candidates()
    mem = ["mgwr", "adw", "raw", "tps"]
    lagged = vcat(mem, ["mgwr_lag", "mgwr_lead", "adw_lag", "adw_lead"])
    blend(src, fb) = (G, b) -> convex(G, src, fb)
    # Band 0 of the agreement axes is the satellite-dry band: toward the raw anchor, not adw.
    dry_raw = (G, b) -> b == 1 ? convex(G, 1, 3) : convex(G, 1, 2)
    return [
        Candidate("C0_const_adw", "reference", ["mgwr", "adw"], :wet, blend(1, 2)),
        Candidate("C1_agrenv_adw", "reference", ["mgwr", "adw"], :agr, blend(1, 2)),
        Candidate("P1a_dist_adw", "pooled", ["mgwr", "adw"], :wet_dist, blend(1, 2)),
        Candidate("P1b_agrenv_dist_adw", "pooled", ["mgwr", "adw"], :agr_dist, blend(1, 2)),
        Candidate("P2a_wet_and_dry_adw", "pooled", ["mgwr", "adw"], :wetdry, blend(1, 2)),
        Candidate("P2b_agrenv_dry_adw", "pooled", ["mgwr", "adw"], :agr0, blend(1, 2)),
        Candidate("P2c_agrenv_dry_raw", "pooled", ["mgwr", "adw", "raw"], :agr0, dry_raw),
        Candidate("P3a_nnls_agr", "pooled", mem, :agr0, (G, b) -> nnls(G, 1:4)),
        Candidate("P3b_nnls_agr_dist", "pooled", mem, :agr0_dist, (G, b) -> nnls(G, 1:4)),
        Candidate("P3c_ols_agr", "pooled", vcat(mem, ["one"]), :agr0, (G, b) -> ols(G, 1:5)),
        Candidate("P3d_ols_agr_dist", "pooled", vcat(mem, ["one"]), :agr0_dist,
            (G, b) -> ols(G, 1:5)),
        Candidate("P4a_ols_lag_agr", "pooled", vcat(lagged, ["one"]), :agr0,
            (G, b) -> ols(G, 1:9)),
        Candidate("P4b_ols_lag_agr_dist", "pooled", vcat(lagged, ["one"]), :agr0_dist,
            (G, b) -> ols(G, 1:9)),
        Candidate("P4c_nnls_lag_agr_dist", "pooled", lagged, :agr0_dist, (G, b) -> nnls(G, 1:8)),
    ]
end

"""Leave-one-fold-out predictions of `c`, and one weights row per `(fold, band)`."""
function fit_candidate(c::Candidate, table, meta, scheme)
    band = band_of(c.band, meta)
    n_bands = maximum(band)
    folds = sort(unique(meta.fold))
    G, n = grams(table, c.features, band, meta.fold, maximum(folds), n_bands)
    identity_w = source_weights(c.features)
    Z = hcat((table[f] for f in c.features)...)
    prediction = copy(table["mgwr"])
    rows = NamedTuple[]
    for fold in folds
        others = [f for f in folds if f != fold]
        train(b) = sum(G[f, b + 1] for f in others)
        n_train(b) = sum(n[f, b + 1] for f in others)
        pooled = c.fit(sum(train(b) for b in 1:n_bands), 0)
        weights = Dict{Int,Vector{Float64}}()
        for b in 1:n_bands
            weights[b] = n_train(b) >= MIN_BAND_CELLS ? c.fit(train(b), b) :
                (n_train(b) == 0 ? identity_w : pooled)
            push!(rows, (; scheme, candidate=c.name, fold, band=b, n_train=n_train(b),
                weights=join(("$(f)=$(round(w; digits=4))" for (f, w) in
                    zip(c.features, weights[b]) if w != 0), "|")))
        end
        for i in eachindex(prediction)
            meta.fold[i] == fold && band[i] != 0 || continue
            prediction[i] = dot(view(Z, i, :), weights[band[i]])
        end
    end
    return prediction, rows
end

"""
Heavy-rain candidates on top of a base prediction. `:top` bands are the three with the envelope at
or above the last edge (8 mm/h).
"""
function heavy_candidates(base, base_name, table, meta, scheme)
    folds = sort(unique(meta.fold))
    top = band_of(:top, meta)
    y = table["y"]
    out = Dict{String,Vector{Float64}}()
    rows = NamedTuple[]

    # H1: least-squares scale per top band; H2: blend toward the max of the raw products there.
    scaled, toward_max = copy(base), copy(base)
    for fold in folds, b in (3, 6, 9)
        train = findall(i -> meta.fold[i] != fold && top[i] == b, eachindex(y))
        test = findall(i -> meta.fold[i] == fold && top[i] == b, eachindex(y))
        isempty(train) && continue
        p, g, yy = base[train], table["gmax"][train], y[train]
        s = dot(p, yy) / dot(p, p)
        d = g .- p
        lambda = clamp(dot(d, yy .- p) / dot(d, d), 0.0, 1.0)
        scaled[test] .= s .* base[test]
        toward_max[test] .= (1 - lambda) .* base[test] .+ lambda .* table["gmax"][test]
        push!(rows, (; scheme, candidate="H1_scale_top($base_name)", fold, band=b,
            n_train=length(train), weights="scale=$(round(s; digits=4))"))
        push!(rows, (; scheme, candidate="H2_to_gmax_top($base_name)", fold, band=b,
            n_train=length(train), weights="lambda=$(round(lambda; digits=4))"))
    end
    out["H1_scale_top($base_name)"] = scaled
    out["H2_to_gmax_top($base_name)"] = toward_max

    # H3: quantile mapping of the prediction onto the training folds' gauge distribution.
    mapped = copy(base)
    for fold in folds
        train = findall(!=(fold), meta.fold)
        test = findall(==(fold), meta.fold)
        p_sorted, y_sorted = sort(base[train]), sort(y[train])
        n = length(p_sorted)
        for i in test
            rank = clamp(searchsortedlast(p_sorted, base[i]), 1, n)
            mapped[i] = y_sorted[rank]
        end
    end
    out["H3_quantile_map($base_name)"] = mapped
    return out, rows
end

"""Daily squared-error sums and counts per stratum, for the day-block bootstrap."""
function daily_sse(prediction, table, meta, mask)
    sse = zeros(meta.n_days)
    n = zeros(Int, meta.n_days)
    y = table["y"]
    for i in eachindex(y)
        mask[i] || continue
        sse[meta.day[i]] += abs2(prediction[i] - y[i])
        n[meta.day[i]] += 1
    end
    return sse, n
end

function score_rows(scheme, name, family, prediction, table, meta, reference, samples)
    y = table["y"]
    strata = [(label, (y .>= low) .& (y .< high)) for (label, low, high) in STRATA]
    for (band, label) in enumerate(("dist_0_20", "dist_20_50", "dist_50_plus"))
        push!(strata, (label, meta.distance .== band))
    end
    rows = NamedTuple[]
    for (label, mask) in strata
        count(mask) == 0 && continue
        e = prediction[mask] .- y[mask]
        rmse = sqrt(mean(abs2, e))
        ref_rmse = sqrt(mean(abs2, reference[mask] .- y[mask]))
        adw_rmse = sqrt(mean(abs2, table["adw"][mask] .- y[mask]))
        ci_low, ci_high = NaN, NaN
        if label in ("all", "heavy")
            s_t, n_t = daily_sse(prediction, table, meta, mask)
            s_r, _ = daily_sse(reference, table, meta, mask)
            gains = [begin
                    nn = sum(n_t[d]); 1 - sqrt(sum(s_t[d]) / nn) / sqrt(sum(s_r[d]) / nn)
                end for d in samples]
            ci_low, ci_high = quantile(gains, 0.025), quantile(gains, 0.975)
        end
        push!(rows, (; scheme, candidate=name, family, stratum=label, n=count(mask),
            RMSE=rmse, Bias=mean(e), gain_vs_blend_mgwr=1 - rmse / ref_rmse,
            ci_low, ci_high, gain_vs_adw=1 - rmse / adw_rmse))
    end
    return rows
end

function screen_scheme(scheme, run_dir, grid, ids, lonlat, sources)
    table, meta = cell_table(joinpath(run_dir, scheme), grid, ids, lonlat, sources)
    println("$scheme: $(length(table["y"])) cells")
    rng = MersenneTwister(SEED)
    days_with_cells = sort(unique(meta.day))
    samples = [rand(rng, days_with_cells, length(days_with_cells)) for _ in 1:BOOTSTRAP_REPS]
    reference = table["blend_mgwr"]

    score_rows_, weight_rows = NamedTuple[], NamedTuple[]
    for (name, family) in (("mgwr", "saved"), ("adw", "saved"), ("raw", "saved"),
            ("blend_mgwr", "saved"), ("blend_agrenv_mgwr", "saved"))
        append!(score_rows_, score_rows(scheme, name, family, table[name], table, meta,
            reference, samples))
    end
    fitted = Dict{String,Vector{Float64}}()
    for c in candidates()
        prediction, rows = fit_candidate(c, table, meta, scheme)
        fitted[c.name] = prediction
        append!(weight_rows, rows)
        append!(score_rows_, score_rows(scheme, c.name, c.family, prediction, table, meta,
            reference, samples))
        println("  $(c.name): RMSE $(round(sqrt(mean(abs2, prediction .- table["y"])); digits=5))")
    end
    # Heavy candidates on the benchmark's own blend and on the best pooled candidate.
    pooled_rmse(name) = sqrt(mean(abs2, fitted[name] .- table["y"]))
    best = argmin(pooled_rmse, [c.name for c in candidates() if c.family == "pooled"])
    for base_name in ("C0_const_adw", best)
        heavy, rows = heavy_candidates(fitted[base_name], base_name, table, meta, scheme)
        append!(weight_rows, rows)
        for (name, prediction) in heavy
            append!(score_rows_, score_rows(scheme, name, "heavy", prediction, table, meta,
                reference, samples))
        end
    end
    return score_rows_, weight_rows
end

function main(args=ARGS)
    run_dir = isempty(args) ? DEFAULT_RUN : abspath(args[1])
    isdir(run_dir) || error("benchmark output directory not found: $run_dir")
    mkpath(OUTDIR)
    write(joinpath(OUTDIR, "source_run.txt"), run_dir * "\n")

    mger = study_config(OUTDIR)
    _, ids, product_data = load_global_common_product_data(mger)
    ids = String.(ids)
    grid = (; times=product_data["GPM"].times, Y_obs=Matrix{Float64}(product_data["GPM"].Y_obs))
    sources = [Matrix{Float64}(product_data[p].Y_sat) for p in SOURCE_PRODUCTS]
    station_meta = load_station_meta(mger.station_meta_path;
        station_id_col=mger.station_id_col, lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    product_data = nothing

    scores, weights = NamedTuple[], NamedTuple[]
    for scheme in SCHEMES
        s, w = screen_scheme(scheme, run_dir, grid, ids, lonlat, sources)
        append!(scores, s)
        append!(weights, w)
        GC.gc()
    end
    write_csv_atomic(joinpath(OUTDIR, "candidates.csv"), DataFrame(scores))
    write_csv_atomic(joinpath(OUTDIR, "weights.csv"), DataFrame(weights))
    println("Wrote the screen to $OUTDIR")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
