using Test, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("LandformClassification")
const LFC = Main.LandformClassification

@testset "relief_class puts each break in the class above it" begin
    @test LFC.relief_class(0.0) == "plain"
    @test LFC.relief_class(29.9) == "plain"
    @test LFC.relief_class(30.0) == "hills"
    @test LFC.relief_class(199.9) == "hills"
    @test LFC.relief_class(200.0) == "small_relief_mountain"
    @test LFC.relief_class(500.0) == "medium_relief_mountain"
    @test LFC.relief_class(1000.0) == "large_relief_mountain"
    @test_throws ArgumentError LFC.relief_class(-1.0)
    @test_throws ArgumentError LFC.relief_class(NaN)
end

@testset "relief_region pools hills with small-relief mountains below 500 m" begin
    @test LFC.relief_region(10.0) == "relief_lt500"
    @test LFC.relief_region(150.0) == "relief_lt500"
    @test LFC.relief_region(499.9) == "relief_lt500"
    @test LFC.relief_region(500.0) == "relief_500_1000"
    @test LFC.relief_region(1000.0) == "relief_ge1000"
    @test LFC.altitude_class(999.0) == "low_altitude"
    @test LFC.altitude_class(1000.0) == "middle_altitude"
end

@testset "mean_change_point finds the start of the second level" begin
    x = [5.0, 5.1, 4.9, 5.0, 1.0, 1.1, 0.9, 1.0, 1.05]
    change = LFC.mean_change_point(x)
    @test change.index == 5
    @test isnan(change.gain[1])
    @test all(k -> change.gain[5] >= change.gain[k], 2:length(x))
    @test_throws ArgumentError LFC.mean_change_point([1.0, 2.0])
    @test_throws ArgumentError LFC.mean_change_point([1.0, NaN, 2.0])
end

@testset "optimal_relief_window analyses log relief per unit area" begin
    side = collect(1.0:6.0)
    # Relief per unit area of 100 for the first three windows, then 10: the change point is window 4.
    relief = [100.0, 400.0, 900.0, 160.0, 250.0, 360.0]
    window = LFC.optimal_relief_window(side, relief)
    @test window.log_relief_per_area ≈ log.(relief ./ side .^ 2)
    @test window.index == 4
    @test window.side_km == 4.0
    @test_throws ArgumentError LFC.optimal_relief_window([1.0, 1.0, 2.0], [1.0, 2.0, 3.0])
    @test_throws DimensionMismatch LFC.optimal_relief_window([1.0, 2.0, 3.0], [1.0, 2.0])
end

@testset "block_relief is the block's largest maximum minus its smallest minimum" begin
    max_grid = [1.0 2.0 3.0 4.0 9.0
                5.0 6.0 7.0 8.0 9.0
                1.0 1.0 1.0 NaN 9.0]
    min_grid = max_grid .- 1.0
    relief = LFC.block_relief(max_grid, min_grid, 2)
    @test size(relief) == (1, 2)     # the trailing row and column do not fill a block
    @test relief[1, 1] == 6.0 - 0.0
    @test relief[1, 2] == 8.0 - 2.0
    @test LFC.block_relief(max_grid, min_grid, 3) == reshape([7.0 - 0.0], 1, 1)
    max_grid[2, 2] = NaN
    @test all(isnan, LFC.block_relief(max_grid, min_grid, 3))   # now the only 3x3 block holds a NaN
    max_grid[2, 2] = 6.0
    @test isequal(LFC.block_relief(max_grid, min_grid, 1), max_grid - min_grid)
    area = LFC.mean_block_relief(max_grid, min_grid, 1:2)
    @test area.n_blocks == [14, 2]
    @test area.mean_relief == [1.0, 6.0]
end

@testset "_fixed3 writes GMT region bounds without an exponent" begin
    @test LFC._fixed3(3474245.66387) == "3474245.664"      # Julia itself would print 3.47424566387e6
    @test LFC._fixed3(417843.1731) == "417843.173"
    @test LFC._fixed3(-12.3456) == "-12.346"
    @test LFC._fixed3(0.0005) == "0.001"                  # the binary value lies just above the tie
    @test LFC._fixed3(7) == "7.000"
end

@testset "parse_xyz_grid rebuilds a north-up grid and maps nodata to NaN" begin
    lines = ["10 25 1", "20 25 2", "", "10 15 3", "20 15 -9999"]
    parsed = LFC.parse_xyz_grid(lines)
    @test parsed.x == [10.0, 20.0]
    @test parsed.y == [25.0, 15.0]
    @test parsed.grid[1, :] == [1.0, 2.0]
    @test parsed.grid[2, 1] == 3.0 && isnan(parsed.grid[2, 2])
    @test_throws ArgumentError LFC.parse_xyz_grid(["10 25 1", "20 25 2", "10 15 3"])
end
