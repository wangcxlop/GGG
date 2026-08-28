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
    name == "NDVIVariableSelection" && return ["ERA5VariableSelection"]
    name == "JointCovariateModels" && return ["DEMTerrainExperiment"]
    name == "JointVariableSelection" &&
        return ["DEMTerrainExperiment", "ERA5VariableSelection", "NDVIVariableSelection"]
    return String[]
end

"""Load one standalone module into `Main`, dependencies first, skipping anything already there."""
function load_standalone_module(name::AbstractString)
    isdefined(Main, Symbol(name)) && return nothing
    for dependency in standalone_module_dependencies(name)
        load_standalone_module(dependency)
    end
    path = joinpath(@__DIR__, "$(name).jl")
    isfile(path) || throw(ArgumentError("no standalone module at $path"))
    # Into `Main` regardless of who called us, so there is exactly one copy no matter which
    # script, test, or module triggered the load.
    Base.include(Main, path)
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
    Base.include(Main, joinpath(@__DIR__, "$(name).jl"))
    return nothing
end
