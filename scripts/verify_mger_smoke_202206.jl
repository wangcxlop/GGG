#!/usr/bin/env julia

using CSV, DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PROCESSED = joinpath(ROOT, "data", "processed")
const STUDY_DATA = joinpath(PROCESSED, "study_area")
const RESULT = joinpath(ROOT, "output", "mger_smoke_202206_5kernels_5fold")
const KERNELS = (
    (0, "gaussian"),
    (1, "exponential"),
    (2, "bisquare"),
    (3, "tricube"),
    (4, "boxcar"),
)

function assert_wide(path::AbstractString, expected_rows::Int, expected_stations::Int)
    table = CSV.read(path, DataFrame)
    @assert nrow(table) == expected_rows "$path has unexpected row count"
    @assert ncol(table) - 1 == expected_stations "$path has unexpected station count"
    return table
end

function assert_station_columns(path::AbstractString, expected_ids::Vector{String})
    header = propertynames(CSV.File(path; limit=1))
    @assert string.(header[2:end]) == expected_ids "$path has unexpected station columns"
end

function main()
    station_meta = CSV.read(joinpath(STUDY_DATA, "station_meta.csv"), DataFrame)
    @assert nrow(station_meta) == 237
    @assert all(109.4 .<= station_meta.lon .<= 111.6)
    @assert all(31.2 .<= station_meta.lat .<= 33.4)
    study_ids = string.(station_meta.station_id)
    @assert issorted(study_ids)
    for path in (
        joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
        joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
        joinpath(STUDY_DATA, "hubei_fy4b_hourly_202206_strict_navcorrected.csv"),
        joinpath(STUDY_DATA, "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"),
    )
        assert_station_columns(path, study_ids)
    end

    assert_wide(
        joinpath(STUDY_DATA, "hubei_gpm_hourly_2022_2025_JunSep_aligned.csv"),
        11_712, 237,
    )
    assert_wide(
        joinpath(STUDY_DATA, "hubei_gsmap_hourly_2022_2025_JunSep_aligned.csv"),
        11_712, 237,
    )
    assert_wide(
        joinpath(STUDY_DATA, "hubei_fy4b_hourly_202206_strict_navcorrected.csv"),
        671, 237,
    )
    assert_wide(
        joinpath(STUDY_DATA, "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"),
        10_770, 237,
    )

    for product in ("gpm", "gsmap")
        duplicate_qc = CSV.read(
            joinpath(ROOT, "output", "input_audit", "$(product)_duplicate_station_times.csv"),
            DataFrame,
        )
        @assert nrow(duplicate_qc) == 11_712
        @assert all(duplicate_qc.source_rows .== 2)
        @assert all(duplicate_qc.unique_value_count .== 1)
    end

    # FY4B source-extraction QC is checked when present. It is not regenerated
    # when the existing statewide FY4B cache is only subset locally.
    qc_path = joinpath(
        ROOT, "output", "input_audit", "fy4b_hourly_qc_202206_strict_navcorrected.csv",
    )
    if isfile(qc_path)
        qc = CSV.read(qc_path, DataFrame)
        @assert nrow(qc) == 720
        @assert count(==("complete"), qc.status) == 671
        @assert count(==("incomplete"), qc.status) == 49
        @assert sum(qc.missing_count) == 82
    end

    full_qc_path = joinpath(
        ROOT, "output", "input_audit",
        "fy4b_hourly_qc_2022_2025_JunSep_strict_navcorrected.csv",
    )
    if isfile(full_qc_path)
        full_qc = CSV.read(full_qc_path, DataFrame)
        @assert nrow(full_qc) == 11_712
        @assert count(==("complete"), full_qc.status) == 10_770
        @assert count(==("incomplete"), full_qc.status) == 942
        @assert sum(full_qc.missing_count) == 1_562
        @assert sum(full_qc.unreadable_count) == 3
        @assert count((full_qc.status .== "incomplete") .& full_qc.native_hourly_available) == 674
    end

    input_audit = CSV.read(
        joinpath(ROOT, "output", "input_audit_smoke_202206", "input_audit.csv"),
        DataFrame,
    )
    @assert all(input_audit.global_common_station_count .== 237)
    @assert all(input_audit.global_common_timestamp_count .== 671)
    @assert all(input_audit.raw_station_count .== 318)
    @assert all(input_audit.study_area_station_count .== 237)
    @assert all(input_audit.study_west .== 109.4)
    @assert all(input_audit.study_east .== 111.6)
    @assert all(input_audit.study_south .== 31.2)
    @assert all(input_audit.study_north .== 33.4)

    full_input_audit = CSV.read(
        joinpath(ROOT, "output", "input_audit", "input_audit.csv"),
        DataFrame,
    )
    @assert all(full_input_audit.global_common_station_count .== 237)
    @assert all(full_input_audit.global_common_timestamp_count .== 8_067)

    status = CSV.read(joinpath(RESULT, "kernel_run_status.csv"), DataFrame)
    @assert nrow(status) == 5
    @assert all(status.status .== "success")
    @assert collect(zip(status.kernel, status.kernel_name)) == collect(KERNELS)

    combined = CSV.read(joinpath(RESULT, "summary_five_kernels_pooled.csv"), DataFrame)
    @assert nrow(combined) == 15
    @assert sort(unique(combined.kernel)) == collect(0:4)
    @assert sort(unique(combined.product)) == ["FY4B", "GPM", "GSMaP"]
    @assert nrow(CSV.read(joinpath(RESULT, "summary_five_kernels_folds.csv"), DataFrame)) == 75
    @assert nrow(CSV.read(joinpath(RESULT, "summary_five_kernels_fold_stats.csv"), DataFrame)) == 15

    reference_split = DataFrame()
    for (kernel, kernel_label) in KERNELS
        kernel_result = joinpath(RESULT, kernel_label)
        split = CSV.read(joinpath(kernel_result, "split_common_spatial5fold.csv"), DataFrame)
        @assert nrow(split) == 237
        @assert length(unique(split.station_id)) == 237
        @assert sort(unique(split.fold)) == collect(1:5)
        @assert sort(combine(groupby(split, :fold), nrow => :n).n) == [47, 47, 47, 48, 48]
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
                671, 237,
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
                @assert nrow(scan) + nrow(failures) == 32
                @assert all(scan.kernel .== kernel)
                @assert all(failures.kernel .== kernel)
            end
        end
    end

    println("MGER full-input and June 2022 five-kernel verification passed")
    return combined
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
