#!/usr/bin/env julia

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "FY4BPreprocessing.jl"))

using .FY4BPreprocessing

const SPECS = [
    (
        year=2022,
        months=[10, 11, 12],
        label="202210_202212",
    ),
    (
        year=2023,
        months=[1, 2, 3, 4, 5, 10, 11, 12],
        label="202301_202305_202310_202312",
    ),
    (
        year=2024,
        months=[1, 2, 3, 4, 5, 10, 11, 12],
        label="202401_202405_202410_202412",
    ),
]

function main()
    selected_years = isempty(ARGS) ?
        Set(spec.year for spec in SPECS) :
        Set(parse.(Int, ARGS))
    known_years = Set(spec.year for spec in SPECS)
    issubset(selected_years, known_years) ||
        error("Supported years are $(sort(collect(known_years)))")

    results = NamedTuple[]
    failures = NamedTuple[]
    for spec in SPECS
        spec.year in selected_years || continue
        println("\nProcessing FY4B $(spec.year), months $(join(spec.months, ", "))")
        try
            result = aggregate_fy4b_hourly(
                data_dir=joinpath(ROOT, "data", "FY4B"),
                station_meta_file=joinpath(ROOT, "data", "hubei_station_meta.csv"),
                output_file=joinpath(
                    ROOT,
                    "data",
                    "processed",
                    "hubei_fy4b_hourly_$(spec.label)_strict_navcorrected.csv",
                ),
                qc_file=joinpath(
                    ROOT,
                    "output",
                    "input_audit",
                    "fy4b_hourly_qc_$(spec.label)_strict_navcorrected.csv",
                ),
                years=spec.year:spec.year,
                months=spec.months,
            )
            push!(results, merge((year=spec.year, months=spec.months), result))
        catch error_value
            message = sprint(showerror, error_value)
            @error "FY4B year could not produce an hourly data file" year=spec.year message
            push!(failures, (year=spec.year, message=message))
        end
    end
    return (; results, failures)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
