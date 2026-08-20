using Test

if !isdefined(@__MODULE__, :TerrainFeatures)
    include(joinpath(@__DIR__, "..", "src", "TerrainFeatures.jl"))
end
using .TerrainFeatures

@testset "terrain feature helpers" begin
    @test collect(aspect_components(90, 10)) ≈ [1.0, 0.0] atol=1e-12
    @test collect(aspect_components(180, 10)) ≈ [0.0, -1.0] atol=1e-12
    @test collect(aspect_components(270, 10)) ≈ [-1.0, 0.0] atol=1e-12
    @test aspect_components(123, 0.01) == (0.0, 0.0)
end
