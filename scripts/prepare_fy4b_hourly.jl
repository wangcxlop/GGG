#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "FY4BPreprocessing.jl"))

using .FY4BPreprocessing

function main()
    smoke = "--smoke-202206" in ARGS
    output_name = smoke ?
        "hubei_fy4b_hourly_202206_strict_navcorrected.csv" :
        "hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv"
    qc_name = smoke ?
        "fy4b_hourly_qc_202206_strict_navcorrected.csv" :
        "fy4b_hourly_qc_2022_2025_JunSep_strict_navcorrected.csv"
    return aggregate_fy4b_hourly(
        data_dir=joinpath(ROOT, "data", "FY4B"),
        station_meta_file=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        output_file=joinpath(ROOT, "data", "processed", output_name),
        qc_file=joinpath(ROOT, "output", "input_audit", qc_name),
        years=smoke ? (2022:2022) : (2022:2025),
        start_month=6,
        stop_month=smoke ? 7 : 10,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
