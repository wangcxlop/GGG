using Test
using Dates

if !isdefined(@__MODULE__, :ERA5LandCovariates)
    include(joinpath(@__DIR__, "..", "src", "ERA5LandCovariates.jl"))
end
using .ERA5LandCovariates

@testset "ERA5-Land covariate helpers" begin
    @test parse_era5_datetime("2022-01-01T00:00:00.0") == DateTime(2022, 1, 1)
    @test parse_era5_datetime("2024-12-31T23:00:00") == DateTime(2024, 12, 31, 23)
    @test parse_era5_datetime(DateTime(2023, 6, 1, 8)) == DateTime(2023, 6, 1, 8)
end
