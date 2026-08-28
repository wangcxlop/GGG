#!/usr/bin/env julia

"""
Golden-diff gate for the conservative refactor.

Compares every file under a baseline snapshot against the corresponding file in a freshly
produced run directory, byte for byte. Byte equality rather than `≈` is deliberate: the refactor
is only allowed to move code, never numbers, so any formatting or value change at all is a
finding. Run it after re-running the smoke pipelines:

    julia --project=. scripts/verify_refactor_invariance.jl \
        output/_refactor_baseline_benchmark \
        output/interpolation_benchmark_smoke_joint_covariates_nested_mgwrintercept_only

Exits non-zero when anything differs, so it can gate a commit.

Delete this script (and the output/_refactor_baseline_* snapshots) once the refactor lands.
"""

const ROOT = normpath(joinpath(@__DIR__, ".."))

"""Every file under `dir`, as paths relative to `dir`, sorted."""
function relative_files(dir::AbstractString)
    files = String[]
    for (root, _, names) in walkdir(dir), name in names
        push!(files, relpath(joinpath(root, name), dir))
    end
    return sort(files)
end

"""First line index at which two files differ, plus both lines, for a diagnosable report."""
function first_difference(baseline_path::AbstractString, current_path::AbstractString)
    baseline_lines = readlines(baseline_path)
    current_lines = readlines(current_path)
    for index in 1:min(length(baseline_lines), length(current_lines))
        if baseline_lines[index] != current_lines[index]
            return (index, baseline_lines[index], current_lines[index])
        end
    end
    return (min(length(baseline_lines), length(current_lines)) + 1,
        "<$(length(baseline_lines)) lines>", "<$(length(current_lines)) lines>")
end

function compare_snapshots(baseline_dir::AbstractString, current_dir::AbstractString)
    isdir(baseline_dir) || error("baseline snapshot does not exist: $baseline_dir")
    isdir(current_dir) || error("run directory does not exist: $current_dir")
    baseline_files = relative_files(baseline_dir)
    current_files = relative_files(current_dir)

    missing_files = setdiff(baseline_files, current_files)
    extra_files = setdiff(current_files, baseline_files)
    changed = String[]
    for name in intersect(baseline_files, current_files)
        read(joinpath(baseline_dir, name)) == read(joinpath(current_dir, name)) ||
            push!(changed, name)
    end

    println("baseline: $baseline_dir")
    println("current:  $current_dir")
    println("$(length(baseline_files)) baseline files, $(length(current_files)) current files")
    for name in missing_files
        println("MISSING  $name")
    end
    for name in extra_files
        println("EXTRA    $name")
    end
    for name in changed
        line, before, after = first_difference(
            joinpath(baseline_dir, name), joinpath(current_dir, name))
        println("CHANGED  $name (first difference at line $line)")
        println("    baseline: $before")
        println("    current:  $after")
    end

    clean = isempty(missing_files) && isempty(extra_files) && isempty(changed)
    println(clean ? "OK: byte-identical" :
        "DRIFT: $(length(changed)) changed, $(length(missing_files)) missing, " *
        "$(length(extra_files)) extra")
    return clean
end

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: verify_refactor_invariance.jl <baseline dir> <run dir>")
    return compare_snapshots(abspath(args[1]), abspath(args[2])) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
