"""
Single-source loader for the standalone `src/` modules.

Each of these files defines its own `module X ... end`, so `include`ing one twice does not
"reload" it - it compiles a second, independent copy. Types are then not interchangeable between
the copies: a `DEMExperimentConfig` built against one is rejected by a function from the other.
Before this loader, `DEMTerrainExperiment.jl` was included into three different namespaces and
`ERA5VariableSelection.jl` into two, and the entry points guarded against it with three competing
idioms (`isdefined(@__MODULE__, ...)`, `isdefined(Main, ...)`, and no guard at all).

Every entry point - scripts, tests, and `InterpolationBenchmark.jl` - now goes through
`load_standalone_modules`, which loads each file into `Main` at most once and pulls in that
module's own sibling dependencies first. Callers then reach the module as `Main.X` (or, once it is
loaded, plain `using .X` from a top-level script).

This file is deliberately safe to `include` more than once: it defines only functions, never
`const`, so re-inclusion is a no-op rather than an "invalid redefinition of constant" error.
"""

"""
Sibling `src/` modules that `name` must be loaded after.

Kept as a function rather than a `Dict` const so this file stays re-includable. The graph is
small and acyclic: everything else in `src/` is a leaf.
"""
function standalone_module_dependencies(name::AbstractString)
    name == "ERA5VariableSelection" && return ["SelectionScaffolding"]
    name == "NDVIVariableSelection" && return ["SelectionScaffolding", "ERA5VariableSelection"]
    name == "JointCovariateModels" && return ["DEMTerrainExperiment", "CovariateGroups"]
    name == "JointVariableSelection" &&
        return ["SelectionScaffolding", "DEMTerrainExperiment", "ERA5VariableSelection",
                "NDVIVariableSelection", "CovariateGroups"]
    # Everything that reads a station table or writes a CSV atomically. Listed as one clause
    # because the dependency is the same for all of them and the list is the interesting part.
    name in ("AppEEARSNDVI", "ERA5LandStations", "FY4BPreprocessing", "MGERDataPrep",
             "MOD13A2NDVIProcessing", "StudyArea", "TerrainFeatures") &&
        return ["TableIO"]
    name == "HeavyRainEvents" && return ["TraditionalInterpolation"]
    return String[]
end

"""
Absolute path of `name.jl`, wherever under `src/` it lives.

`src/` is arranged in tiers rather than flat, so a module name no longer determines a directory.
Searching rather than tabulating keeps the property this loader has always had - a new file in
`src/` is loadable by name with nothing to register anywhere - where a name-to-directory map would
have to enumerate every module with no safe default, unlike `standalone_module_dependencies`
above, which is a sparse exception list. The walk costs nothing: both callers below return early
on `isdefined`, so it runs at most once per name per process.

Two matches is the failure worth naming. It cannot happen while every basename is unique, but a
leftover copy from a half-finished move would otherwise be picked by walk order, silently, and the
symptom would be a module whose contents do not match the file anyone is editing.
"""
function locate_src_file(name::AbstractString)
    matches = String[]
    for (root, _, files) in walkdir(@__DIR__)
        "$(name).jl" in files && push!(matches, joinpath(root, "$(name).jl"))
    end
    isempty(matches) &&
        throw(ArgumentError("no file named $(name).jl anywhere under $(@__DIR__)"))
    length(matches) > 1 && throw(ArgumentError(
        "$(name).jl exists in more than one place under src/, refusing to guess: " *
        join(sort(matches), ", ")))
    return only(matches)
end

"""Load one standalone module into `Main`, dependencies first, skipping anything already there."""
function load_standalone_module(name::AbstractString)
    isdefined(Main, Symbol(name)) && return nothing
    for dependency in standalone_module_dependencies(name)
        load_standalone_module(dependency)
    end
    # Into `Main` regardless of who called us, so there is exactly one copy no matter which
    # script, test, or module triggered the load.
    Base.include(Main, locate_src_file(name))
    return nothing
end

"""
Load the named standalone `src/` modules into `Main`, once each, in dependency order.

    load_standalone_modules("JointCovariateModels", "TraditionalInterpolation")

Names are the module names, which match the file basenames.
"""
function load_standalone_modules(names::AbstractString...)
    for name in names
        load_standalone_module(name)
    end
    return nothing
end

"""
Sentinel name a top-level `src/` pipeline fragment defines once it is loaded.

`MGERPipeline.jl` and `InterpolationBenchmark.jl` are not modules - they are plain top-level
fragments evaluated straight into `Main` - so there is no module name to test for. Each is
detected by a struct it defines instead.
"""
function pipeline_sentinel(name::AbstractString)
    name == "MGERPipeline" && return :MGERConfig
    name == "InterpolationBenchmark" && return :InterpolationBenchmarkConfig
    throw(ArgumentError("unknown pipeline fragment: $name"))
end

"""
Load a top-level `src/` pipeline fragment into `Main` once.

    load_pipeline("InterpolationBenchmark")

`InterpolationBenchmark.jl` pulls in `MGERPipeline.jl` and the standalone modules it needs itself,
so asking for it is enough.
"""
function load_pipeline(name::AbstractString)
    isdefined(Main, pipeline_sentinel(name)) && return nothing
    Base.include(Main, locate_src_file(name))
    return nothing
end
