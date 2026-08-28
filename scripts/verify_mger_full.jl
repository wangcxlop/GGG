#!/usr/bin/env julia

# Structural verification for `scripts/run_mger_full.jl`, mirroring
# `scripts/verify_mger_smoke_202206.jl` at full scale. Pass `--full-year` to verify a run
# produced with the matching flag.
#
# This checks shape, coverage and cross-kernel comparability -- it deliberately makes no claim
# about whether the corrected fields are any *good*. That question is answered by the
# interpolation benchmark's claim assessment, not here.

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("MGERDataPrep")
using Main.MGERDataPrep
const PROCESSED = joinpath(ROOT, "data", "processed")
const STUDY_DATA = joinpath(PROCESSED, "study_area")
const KERNELS = (
    (0, "gaussian"),
    (1, "exponential"),
    (2, "bisquare"),
    (3, "tricube"),
    (4, "boxcar"),
)



function main(args=ARGS)
    full_year = "--full-year" in args
    result = joinpath(
        ROOT, "output",
        "mger_full_5kernels_5fold" * (full_year ? "_fullyear" : ""),
    )
    # Common-time counts measured with `scripts/audit_mger_inputs.jl`; these are the same
    # constants `run_mger_full.jl` pins via `expected_common_time_count`.
    expected_times = full_year ? 11_426 : 8_067
    fy4b_name = full_year ?
        "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv" :
        "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"
    fy4b_rows = full_year ? 11_426 : 10_770

    station_meta = CSV.read(joinpath(STUDY_DATA, "station_meta.csv"), DataFrame)
    @assert nrow(station_meta) == 237
    @assert all(109.4 .<= station_meta.lon .<= 111.6)
    @assert all(31.2 .<= station_meta.lat .<= 33.4)
    study_ids = string.(station_meta.station_id)
    @assert issorted(study_ids)

    satellite_inputs = full_year ?
        [fy4b_name, "hubei_gpm_hourly_2022_2024_full_aligned.csv",
         "hubei_gsmap_hourly_2022_2024_full_aligned.csv"] :
        [fy4b_name, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv",
         "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"]
    for name in vcat("hubei_obs_hourly_2022_2025_JunSep.csv", satellite_inputs)
        assert_station_columns(joinpath(STUDY_DATA, name), study_ids)
    end
    assert_wide(joinpath(STUDY_DATA, fy4b_name), fy4b_rows, 237)
    for name in satellite_inputs[2:end]
        assert_wide(joinpath(STUDY_DATA, name), full_year ? 26_304 : 11_712, 237)
    end

    # The Jun-Sep audit is the one `scripts/audit_mger_inputs.jl` produces without flags. There
    # is no full-year audit directory, so that window is verified against the run itself only.
    if !full_year
        input_audit = CSV.read(
            joinpath(ROOT, "output", "input_audit", "input_audit.csv"), DataFrame,
        )
        @assert all(input_audit.global_common_station_count .== 237)
        @assert all(input_audit.global_common_timestamp_count .== expected_times)
        @assert all(input_audit.raw_station_count .== 318)
        @assert all(input_audit.study_area_station_count .== 237)
    end

    status = CSV.read(joinpath(result, "kernel_run_status.csv"), DataFrame)
    @assert nrow(status) == 5
    @assert all(status.status .== "success")
    @assert collect(zip(status.kernel, status.kernel_name)) == collect(KERNELS)

    combined = CSV.read(joinpath(result, "summary_five_kernels_pooled.csv"), DataFrame)
    @assert nrow(combined) == 15
    @assert sort(unique(combined.kernel)) == collect(0:4)
    @assert sort(unique(combined.product)) == ["FY4B", "GPM", "GSMaP"]
    @assert nrow(CSV.read(joinpath(result, "summary_five_kernels_folds.csv"), DataFrame)) == 75
    @assert nrow(CSV.read(joinpath(result, "summary_five_kernels_fold_stats.csv"), DataFrame)) == 15

    reference_split = DataFrame()
    for (kernel, kernel_label) in KERNELS
        kernel_result = joinpath(result, kernel_label)
        split = CSV.read(joinpath(kernel_result, "split_common_spatial5fold.csv"), DataFrame)
        @assert nrow(split) == 237
        @assert length(unique(split.station_id)) == 237
        @assert sort(unique(split.fold)) == collect(1:5)
        @assert sort(combine(groupby(split, :fold), nrow => :n).n) == [47, 47, 47, 48, 48]
        # The filename says "spatial5fold" but the run passes `fold_scheme=:random`, so assert the
        # scheme the split was actually built with. Without this a scheme swap - which changes what
        # the numbers mean - passes verification silently.
        @assert all(==("random"), split.fold_scheme)
        scope = CSV.read(joinpath(kernel_result, "validation_scope.csv"), DataFrame)
        @assert only(scope.value[scope.key .== "fold_scheme"]) == "random"
        @assert occursin("NOT supported", only(scope.value[scope.key .== "supported_claim"]))
        # Every kernel must see the identical partition, or the pooled table is not a paired
        # comparison and the kernel ranking in it means nothing.
        if isempty(reference_split)
            reference_split = split
        else
            @assert split == reference_split
        end

        summary = CSV.read(joinpath(kernel_result, "summary_three_products.csv"), DataFrame)
        fold_summary = CSV.read(
            joinpath(kernel_result, "summary_three_products_folds.csv"), DataFrame,
        )
        @assert sort(summary.product) == ["FY4B", "GPM", "GSMaP"]
        @assert all(summary.k .== 5)
        @assert all(summary.n_station .== 237)
        @assert all(summary.coverage .>= 0.95)
        @assert all(fold_summary.scan_coverage .>= 0.95)
        @assert all(isfinite, summary.RMSE_post)

        for product in ("FY4B", "GPM", "GSMaP")
            pooled = assert_wide(
                joinpath(kernel_result, "corr_$(product)_spatial5fold_val.csv"),
                expected_times, 237,
            )
            @assert names(pooled)[2:end] == string.(sort(parse.(Int, names(pooled)[2:end])))
            product_split = CSV.read(
                joinpath(kernel_result, "split_$(product)_spatial5fold.csv"), DataFrame,
            )
            @assert product_split == split
            for fold in 1:5
                fold_dir = joinpath(kernel_result, "fold_$(fold)")
                scan = CSV.read(joinpath(fold_dir, "scan_$(product)_spatialcv.csv"), DataFrame)
                failures = CSV.read(
                    joinpath(fold_dir, "scan_$(product)_spatialcv_failures.csv"), DataFrame,
                )
                # 8 bandwidths (4 adaptive + 4 fixed) x 4 slope_ridge candidates.
                @assert nrow(scan) + nrow(failures) == 32
                @assert all(scan.kernel .== kernel)
                @assert all(failures.kernel .== kernel)
            end
        end
    end

    println(
        "MGER full five-kernel verification passed (",
        full_year ? "full-year" : "Jun-Sep", ", ", expected_times, " common hours)",
    )
    return combined
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
