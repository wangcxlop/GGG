using Test
using Dates
using DataFrames

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("MOD13A2NDVIProcessing")
using Main.MOD13A2NDVIProcessing

const TEST_NDVI_COLUMN = Symbol("MOD13A2_061__1_km_16_days_NDVI")
const TEST_DOY_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_composite_day_of_the_year")
const TEST_RELIABILITY_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_pixel_reliability")
const TEST_LAND_COLUMN =
    Symbol("MOD13A2_061__1_km_16_days_VI_Quality_Land/Water_Mask_Description")

@testset "MOD13A2 NDVI helpers" begin
    @test parse_mod13a2_date("2022-01-01") == Date(2022, 1, 1)
    @test parse_mod13a2_date(DateTime(2024, 12, 18, 12)) == Date(2024, 12, 18)
    @test composite_observation_date("2024-01-01", 60.0) == Date(2024, 2, 29)
    @test ismissing(composite_observation_date(Date(2022, 1, 1), -1.0))

    @test quality_class(-1.0) == "fill"
    @test quality_class(0.0) == "good"
    @test quality_class(1) == "marginal"
    @test quality_class(2) == "snow_ice"
    @test quality_class(3) == "cloudy"
    @test quality_class(99) == "unknown"

    @test clean_ndvi(-0.2, 0) == -0.2
    @test clean_ndvi(1.0, 1) == 1.0
    @test ismissing(clean_ndvi(-3000.0, -1))
    @test ismissing(clean_ndvi(0.5, 2))
    @test ismissing(clean_ndvi(0.5, 3))
    @test clean_ndvi(0.5, 0; land_class="land", require_land=true) == 0.5
    @test ismissing(clean_ndvi(0.5, 0; land_class="shoreline", require_land=true))
end

@testset "MOD13A2 table processing" begin
    land = "Land (Nothing else but land)"
    water = "Shallow inland water"
    raw = DataFrame(
        ID=["A", "A", "B", "B"],
        Latitude=[32.0, 32.0, 33.0, 33.0],
        Longitude=[110.0, 110.0, 111.0, 111.0],
        Date=["2022-01-01", "2022-01-17", "2022-01-01", "2022-01-17"],
    )
    raw[!, TEST_NDVI_COLUMN] = [0.3, -3000.0, -0.2, 1.0]
    raw[!, TEST_DOY_COLUMN] = [3.0, -1.0, 1.0, 18.0]
    raw[!, TEST_RELIABILITY_COLUMN] = [0.0, -1.0, 1.0, 0.0]
    raw[!, TEST_LAND_COLUMN] = [land, land, water, land]

    result = process_mod13a2_table(
        raw;
        start_date=Date(2022, 1, 1),
        end_date=Date(2022, 1, 31),
        expected_station_count=2,
        expected_period_count=2,
        expected_first_date=Date(2022, 1, 1),
        expected_last_date=Date(2022, 1, 17),
    )
    @test size(result.table) == (4, 11)
    @test result.table.ndvi_qc[1] == 0.3
    @test ismissing(result.table.ndvi_qc[2])
    @test result.table.ndvi_qc[3] == -0.2
    @test ismissing(result.table.ndvi_land_qc[3])
    @test result.table.observation_date[1] == Date(2022, 1, 3)
    @test ismissing(result.table.observation_date[2])
    @test result.audit.periods == [2, 2]
    @test result.audit.has_nonland_pixel == [false, true]
end
