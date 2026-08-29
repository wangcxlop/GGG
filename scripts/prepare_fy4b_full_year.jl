#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("FY4BPreprocessing")

using Main.FY4BPreprocessing

function main()
    return aggregate_fy4b_hourly(
        data_dir=joinpath(ROOT, "data", "FY4B"),
        station_meta_file=joinpath(ROOT, "data", "hubei_station_meta.csv"),
        output_file=joinpath(
            ROOT, "data", "processed", "hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv",
        ),
        qc_file=joinpath(
            ROOT, "output", "input_audit", "fy4b_hourly_qc_2022_2024_full_strict_navcorrected.csv",
        ),
        years=2022:2024,
        months=1:12,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
