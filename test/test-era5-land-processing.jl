using Test

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("ERA5LandProcessing")
using Main.ERA5LandProcessing

@testset "ERA5-Land derived variables" begin
    @test relative_humidity(293.15, 293.15) == 100.0
    @test 49 < relative_humidity(293.15, 282.42) < 51
    @test wind_speed(3, 4) == 5.0
    @test wind_direction(0, -1) == 0.0
    @test wind_direction(-1, 0) == 90.0
    @test wind_direction(0, 1) == 180.0
    @test wind_direction(1, 0) == 270.0
end
