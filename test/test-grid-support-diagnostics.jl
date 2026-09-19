using Test, DataFrames, Dates, Statistics

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("GridSupportDiagnostics")
const GSD = Main.GridSupportDiagnostics

# The sub-satellite point of the FY4B files this project holds (`N_DISK_1330E`), and a point in
# the middle of the study box, which is 23 degrees of longitude away from it.
const SAT_LON = 133.0
const STUDY_LON, STUDY_LAT = 110.5, 32.3

@testset "cell indices follow the grid each product was sampled on" begin
    @test GSD.latlon_cell_index(110.5, 32.3) == (2905, 1223)
    @test all(GSD.latlon_cell_center(2905, 1223) .≈ (110.55, 32.35))   # spans 110.5-110.6
    # Two points inside one 0.1 degree cell share an index; a point across the edge does not.
    @test GSD.latlon_cell_index(110.51, 32.31) == GSD.latlon_cell_index(110.58, 32.38)
    @test GSD.latlon_cell_index(110.51, 32.31) != GSD.latlon_cell_index(110.61, 32.31)
    @test GSD.latlon_cell_index(110.5, 32.3; step_deg=0.25) !=
        GSD.latlon_cell_index(110.5, 32.3; step_deg=0.1)

    center = (Main.FY4BPreprocessing.GRID_SIZE - 1) / 2
    ix, iy, x, y = GSD.fy4b_pixel_index(SAT_LON, 0.0; sat_lon=SAT_LON)
    @test (x, y) == (center, center)                 # nadir lands on the grid centre
    # An even GRID_SIZE puts that centre exactly on a pixel boundary, so nadir is the one place
    # the rounding is a tie; Julia breaks it to even, giving the pixel above.
    @test (ix, iy) == (round(Int, center) + 1, round(Int, center) + 1) == (1375, 1375)

    # The far hemisphere is NOT rejected - see `fy4b_pixel_index`. Asserted so the day the
    # upstream limb test starts working, this says so instead of silently changing meaning.
    far_side = GSD.fy4b_pixel_index(SAT_LON - 180, 0.0; sat_lon=SAT_LON)
    @test far_side[1] != 0                           # not rejected ...
    @test abs(far_side[1] - ix) <= 1 && abs(far_side[2] - iy) <= 1   # ... it folds onto nadir

    # A gauge-realistic point is far from the centre and comfortably inside the grid.
    sx, sy = GSD.fy4b_pixel_index(STUDY_LON, STUDY_LAT; sat_lon=SAT_LON)[1:2]
    @test 1 <= sx <= Main.FY4BPreprocessing.GRID_SIZE
    @test 1 <= sy <= Main.FY4BPreprocessing.GRID_SIZE
    @test sy < center                                # north of the sub-satellite point
end

@testset "pixel footprint is 4 km at nadir and larger over the study area" begin
    nadir = GSD.pixel_footprint_km(SAT_LON, 0.0; sat_lon=SAT_LON)
    @test nadir.along_scan_km ≈ 4.0 atol = 0.01
    # Cross-scan is slightly wider even at nadir: the projection uses geocentric latitude, so a
    # degree of latitude subtends less scan angle than a degree of longitude does.
    @test 4.0 <= nadir.cross_scan_km <= 4.1
    @test nadir.nadir_ratio ≈ nadir.along_scan_km * nadir.cross_scan_km / 16 atol = 1e-6

    study = GSD.pixel_footprint_km(STUDY_LON, STUDY_LAT; sat_lon=SAT_LON)
    @test study.along_scan_km > nadir.along_scan_km
    @test study.cross_scan_km > nadir.cross_scan_km
    @test study.nadir_ratio > 1.5
    # Symmetric about the sub-satellite longitude and about the equator.
    mirrored = GSD.pixel_footprint_km(2 * SAT_LON - STUDY_LON, -STUDY_LAT; sat_lon=SAT_LON)
    @test mirrored.area_km2 ≈ study.area_km2 rtol = 1e-6

    # The finite-difference step is not load-bearing: halving it moves nothing that is reported.
    @test GSD.pixel_footprint_km(STUDY_LON, STUDY_LAT; sat_lon=SAT_LON, h=5e-4).area_km2 ≈
        study.area_km2 rtol = 1e-6
end

@testset "assignment and multiplicity count the gauges a cell cannot separate" begin
    ids = ["a", "b", "c"]
    lon = [110.51, 110.58, 110.95]
    lat = [32.31, 32.38, 32.05]
    assignment = GSD.pixel_assignment_table(ids, lon, lat; sat_lon=SAT_LON)
    @test nrow(assignment) == 3
    @test assignment.coarse_cell[1] == assignment.coarse_cell[2] != assignment.coarse_cell[3]
    @test all(assignment.fy4b_offset_pixels .<= sqrt(0.5))   # a rounded offset cannot exceed this
    @test all(assignment.fy4b_offset_km .< assignment.fy4b_along_scan_km)

    membership = GSD.cell_membership_table(ids, Dict("coarse" => assignment.coarse_cell))
    @test nrow(membership) == 2
    @test membership[membership.n_stations .== 2, :station_ids] == ["a;b"]

    summary = GSD.cell_multiplicity_summary(membership)
    @test nrow(summary) == 1
    row = only(eachrow(summary))
    @test row.n_stations == 3 && row.n_occupied_cells == 2
    @test row.max_stations_per_cell == 2 && row.n_cells_with_multiple == 1
    @test row.n_stations_sharing == 2 && row.fraction_stations_sharing ≈ 2 / 3
end

@testset "a product sampled through two projections merges only what both merge" begin
    # a and b share a cell on the first grid, b and c on the second: only the pairs that share on
    # both are handed the same value, so nothing is merged here.
    first_grid = ["p1", "p1", "p2"]
    second_grid = ["q1", "q2", "q2"]
    effective = GSD.effective_cell_keys([first_grid, second_grid])
    @test length(unique(effective)) == 3

    # Agreeing grids leave the grouping alone.
    @test length(unique(GSD.effective_cell_keys([first_grid, first_grid]))) == 2
    @test GSD.effective_cell_keys([first_grid]) == first_grid
    @test_throws DimensionMismatch GSD.effective_cell_keys([first_grid, ["q1", "q2"]])
    @test_throws ArgumentError GSD.effective_cell_keys(Vector{String}[])

    # Keys taken from the series grouping: members share one key, everything else is on its own.
    keys = GSD.series_group_keys(["a", "b", "c", "d"], [[1, 3]])
    @test keys[1] == keys[3]
    @test length(unique(keys)) == 3
    @test keys[2] != keys[4]
    @test GSD.series_group_keys(["a", "b"], Vector{Int}[]) == ["solo_a", "solo_b"]
end

@testset "identical series are found without assuming a grid" begin
    #      1  2  3  4  5
    Y = [0.0 1.0 2.0 0.0 3.0        # a
         0.0 1.0 2.0 0.0 3.0        # b: byte-identical to a
         0.0 1.0 2.5 0.0 3.0        # c: differs in one hour only
         NaN 1.0 2.0 NaN 3.0        # d: a with gaps, equal wherever both are finite
         NaN NaN NaN NaN NaN]       # e: no overlap with anything
    ids = ["a", "b", "c", "d", "e"]
    groups, overlaps = GSD.identical_series_groups(Y)
    @test groups == [[1, 2, 4]]
    @test overlaps[1] == 3          # a is joined to d on three shared hours
    @test overlaps[2] == 3
    @test !haskey(overlaps, 3) && !haskey(overlaps, 5)

    # A single differing hour is enough to separate two gauges, and an all-missing series is
    # never called identical to anything.
    @test GSD.identical_series_groups(Y; min_overlap=4) == ([[1, 2]], Dict(1 => 5, 2 => 5))

    table = GSD.identical_series_table(ids, Y, "FY4B")
    @test nrow(table) == 1
    @test only(table.station_ids) == "a;b;d"
    @test only(table.min_overlap_hours) == 3
    @test only(table.n_stations) == 3
    @test only(table.product) == "FY4B"
end

@testset "geometry is checked against what the data can actually distinguish" begin
    Y = [0.0 1.0 2.0 0.0 3.0
         0.0 1.0 2.0 0.0 3.0
         0.0 1.0 2.5 0.0 3.0
         NaN 1.0 2.0 NaN 3.0
         NaN NaN NaN NaN NaN]
    ids = ["a", "b", "c", "d", "e"]

    # The right grid: the one cell holding three gauges is exactly the group the data shows.
    agreed = GSD.cell_grouping_agreement(ids, ["p1", "p1", "p2", "p1", "p3"], Y;
        grid="fy4b_pixel", product="FY4B")
    @test only(agreed.n_cell_groups) == 1
    @test only(agreed.cell_groups_confirmed) == 1
    @test only(agreed.grid_explains_series)
    @test only(agreed.stations_identical_across_cells) == 0
    @test only(agreed.n_stations_in_cell_groups) == 3
    @test only(agreed.n_series_groups) == 1

    # A wrong grid: it puts b with c, which the data separates, and splits b from a, which it
    # does not. Both failures are visible.
    wrong = GSD.cell_grouping_agreement(ids, ["p1", "p2", "p2", "p1", "p3"], Y;
        grid="fy4b_pixel", product="FY4B")
    @test only(wrong.n_cell_groups) == 2
    @test only(wrong.cell_groups_confirmed) == 1
    @test !only(wrong.grid_explains_series)
    @test only(wrong.stations_identical_across_cells) == 1

    # Gauges that are identical only because neither ever reports rain are context, not a failure.
    dry = zeros(2, 4)
    lonely = GSD.cell_grouping_agreement(["a", "b"], ["p1", "p2"], dry;
        grid="coarse", product="GPM")
    @test only(lonely.n_cell_groups) == 0
    @test only(lonely.grid_explains_series)            # nothing to contradict
    @test only(lonely.stations_identical_across_cells) == 1
end

@testset "pair disagreement separates wet hours from the dry majority" begin
    a = [0.0, 1.0, 2.0, NaN]
    b = [0.0, 1.0, 4.0, 5.0]
    stats = GSD.pair_disagreement(a, b)
    @test stats.n == 3
    @test stats.RMSE ≈ sqrt(4 / 3)
    @test stats.MAE ≈ 2 / 3
    @test stats.Bias ≈ 2 / 3
    @test stats.n_wet == 2                       # hours 2 and 3; hour 1 is dry at both gauges
    @test stats.RMSE_wet ≈ sqrt(2)
    @test stats.mean_mm ≈ 1.0

    # A gauge wet at only one of the pair still makes the hour wet: those are the disagreements
    # the number exists to capture.
    @test GSD.pair_disagreement([0.0, 0.0], [0.0, 5.0]).n_wet == 1
    @test GSD.pair_disagreement([NaN, NaN], [1.0, 2.0]).n == 0
end

@testset "gauge pairs carry their cell-sharing flags and bin into a control" begin
    ids = ["a", "b", "c"]
    lonlat = [110.51 32.31; 110.58 32.38; 110.95 32.05]
    Y = [0.0 1.0 2.0 0.0
         0.0 1.0 3.0 0.0
         0.0 2.0 2.0 1.0]
    cells = Dict("coarse" => ["1_1", "1_1", "2_2"], "fine" => ["1_1", "1_2", "2_2"])
    pairs = GSD.gauge_pair_table(ids, lonlat, Y, cells; max_separation_km=100.0)
    @test nrow(pairs) == 3
    @test Set(names(pairs)) ⊇ Set(["shares_coarse", "shares_fine", "separation_km", "RMSE"])
    ab = only(eachrow(pairs[(pairs.station_a .== "a") .& (pairs.station_b .== "b"), :]))
    @test ab.shares_coarse && !ab.shares_fine
    @test ab.RMSE ≈ 0.5                          # one hour differs by 1.0, over four hours
    @test ab.separation_km < 12

    @test nrow(GSD.gauge_pair_table(ids, lonlat, Y, cells; max_separation_km=12.0)) == 1

    summary = GSD.pair_disagreement_summary(pairs)
    @test Set(summary.grid) == Set(["coarse", "fine"])
    shared = summary[(summary.grid .== "coarse") .& summary.shares_cell, :]
    @test nrow(shared) == 1 && only(shared.n_pairs) == 1
    @test sum(summary[summary.grid .== "coarse", :n_pairs]) == 3
end

@testset "bilinear sampling renormalises over the neighbours that passed QC" begin
    indices, weights = GSD.bilinear_neighbours(10.0, 20.0)
    @test indices == ((11, 21), (12, 21), (11, 22), (12, 22))
    @test weights == (1.0, 0.0, 0.0, 0.0)        # a pixel centre reduces to nearest
    @test sum(weights) ≈ 1.0

    _, half = GSD.bilinear_neighbours(10.5, 20.0)
    @test half == (0.5, 0.5, 0.0, 0.0)
    @test GSD.bilinear_combine((1.0, 3.0, NaN, NaN), half) ≈ 2.0
    # One neighbour rejected: the rest are rescaled rather than the cell becoming a hole.
    @test GSD.bilinear_combine((1.0, NaN, NaN, NaN), half) ≈ 1.0
    @test isnan(GSD.bilinear_combine((NaN, NaN, NaN, NaN), half))
    # Finite values that all carry zero weight are not a usable estimate either.
    @test isnan(GSD.bilinear_combine((1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0, 0.0)))

    _, corner = GSD.bilinear_neighbours(10.5, 20.5)
    @test all(corner .≈ 0.25)
    @test GSD.bilinear_combine((1.0, 2.0, 3.0, 4.0), corner) ≈ 2.5

    # Negative fractional coordinates still floor toward the lower-left neighbour.
    @test GSD.bilinear_neighbours(-0.5, 0.25)[1][1] == (0, 1)
end

@testset "extraction sensitivity reports the delta and the wet/dry flips" begin
    ids = ["a", "b"]
    nearest = [0.0 1.0 0.0 2.0
               0.0 0.0 0.0 0.0]
    bilinear = [0.0 1.5 0.2 2.0
                0.0 0.0 0.0 0.0]
    table = GSD.extraction_sensitivity_table(ids, nearest, bilinear)
    @test nrow(table) == 3
    @test table.scope == ["all_stations", "station", "station"]
    overall = only(eachrow(table[table.scope .== "all_stations", :]))
    @test overall.n == 8
    @test overall.mean_abs_delta ≈ 0.7 / 8
    @test overall.rmse_delta ≈ sqrt((0.25 + 0.04) / 8)
    @test overall.max_abs_delta ≈ 0.5
    @test overall.n_wet_flip == 1                # hour 3 at station a: 0.0 vs 0.2
    # Wet-hour views use the union of the two modes' wet calls: hours 2, 3 and 4 at station a.
    @test overall.n_wet_either == 3
    @test overall.rmse_delta_wet ≈ sqrt((0.25 + 0.04) / 3)
    @test overall.mean_abs_delta_wet ≈ 0.7 / 3
    @test overall.fraction_wet_flip_of_wet ≈ 1 / 3
    @test overall.fraction_wet_flip ≈ 1 / 8      # the same flip against every hour, wet or dry
    @test overall.reference_mismatch == -1       # not checked when no reference is given
    @test only(table[table.station_id .== "b", :mean_abs_delta]) ≈ 0.0

    # With the shipped column supplied, an exact reproduction reports zero mismatches.
    @test only(GSD.extraction_sensitivity_table(ids, nearest, bilinear;
        Y_reference=copy(nearest))[1:1, :reference_mismatch]) == 0
    perturbed = copy(nearest); perturbed[1, 2] += 0.1
    @test only(GSD.extraction_sensitivity_table(ids, nearest, bilinear;
        Y_reference=perturbed)[1:1, :reference_mismatch]) == 1
end

@testset "terrain representativeness measures what the gauges never sample" begin
    gauges = [100.0, 200.0]
    domain = [50.0, 100.0, 150.0, 200.0, 300.0]
    table = GSD.terrain_representativeness_table(gauges, domain, "elevation_m";
        probs=[0.0, 0.5, 1.0])
    @test only(table[table.statistic .== "p0", :gauges]) == 100.0
    @test only(table[table.statistic .== "p100", :domain]) == 300.0
    @test only(table[table.statistic .== "mean", :gauges]) == 150.0
    @test only(table[table.statistic .== "n", :domain]) == 5.0
    @test only(table[table.statistic .== "domain_fraction_below_lowest_gauge", :domain]) ≈ 0.2
    @test only(table[table.statistic .== "domain_fraction_above_highest_gauge", :domain]) ≈ 0.2
    @test only(table[table.statistic .== "domain_fraction_outside_gauge_range", :domain]) ≈ 0.4
    @test all(table.variable .== "elevation_m")

    # Non-finite cells are dropped rather than poisoning the quantiles.
    with_gaps = GSD.terrain_representativeness_table(gauges, [domain; NaN], "elevation_m";
        probs=[1.0])
    @test only(with_gaps[with_gaps.statistic .== "p100", :domain]) == 300.0
    @test_throws ErrorException GSD.terrain_representativeness_table(
        Float64[], domain, "elevation_m")
end
