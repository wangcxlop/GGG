using Test
using CSV, DataFrames, Dates, Random, Statistics

if !isdefined(Main, :NDVIVariableSelection)
    include(joinpath(@__DIR__, "..", "src", "NDVIVariableSelection.jl"))
end
using .NDVIVariableSelection

function synthetic_ndvi_inputs(n::Int=40, nperiods::Int=12; seed::Int=71)
    rng = MersenneTwister(seed)
    ids = ["S$(i)" for i in 1:n]
    periods = [Date(2022, 5, 1) + Day(16 * (period - 1)) for period in 1:nperiods]
    rows = NamedTuple[]
    spatial = collect(range(-1.0, 1.0; length=n))
    for station in 1:n, period in 1:nperiods
        value = clamp(0.55 + 0.18 * spatial[station] + 0.08 * sin(period), -0.2, 1.0)
        push!(rows, (
            station_id=ids[station], lon=109.0 + station / n * 3,
            lat=29.5 + 0.4 * sin(station / 5), composite_start=periods[period],
            observation_date=periods[period] + Day(2), ndvi_qc=value,
            pixel_reliability=0, quality_class="good",
            land_water_class=station == n ? "shoreline" : "land",
        ))
    end
    table = DataFrame(rows)
    start_time = DateTime(first(periods))
    end_time = DateTime(last(periods) + Day(20))
    times = collect(start_time:Hour(12):end_time)
    aligned = align_ndvi_asof(table, ids, times)
    residual = fill(0.0, n, length(times))
    for station in 1:n, time in eachindex(times)
        value = aligned.values[station, time]
        residual[station, time] = isfinite(value) ? 3value + 0.02randn(rng) : 0.0
    end
    Yobs = 2 .+ residual .- minimum(residual)
    Ysat = Yobs .- residual
    lonlat = hcat(109 .+ collect(1:n) ./ n .* 3, 29.5 .+ 0.4 .* sin.(collect(1:n) ./ 5))
    return (; ids, periods, table, times, aligned, Yobs, Ysat, lonlat)
end

@testset "NDVI as-of alignment" begin
    table = DataFrame(
        station_id=["A", "A"], lon=[110.0, 110.0], lat=[30.0, 30.0],
        composite_start=[Date(2022, 1, 1), Date(2022, 2, 18)],
        observation_date=[Date(2022, 1, 3), Date(2022, 2, 20)],
        ndvi_qc=Union{Missing,Float64}[0.5, missing], pixel_reliability=[0, 3],
        quality_class=["good", "cloudy"], land_water_class=["shoreline", "shoreline"],
    )
    times = [
        DateTime(2022, 1, 3, 23), DateTime(2022, 1, 4),
        DateTime(2022, 2, 4, 23), DateTime(2022, 2, 5),
    ]
    aligned = align_ndvi_asof(table, ["A"], times; max_age_days=32)
    @test aligned.source_period[1, 1] == 0
    @test aligned.values[1, 2] == 0.5
    @test aligned.values[1, 3] == 0.5
    @test aligned.source_period[1, 4] == 0
    @test aligned.qc.interpolation[1] == "none"
    @test aligned.intervals.effective_time[1] == DateTime(2022, 1, 4)
    @test aligned.intervals.valid_until_exclusive[1] == DateTime(2022, 2, 5)
end

@testset "NDVI quality audit retains non-land stations" begin
    data = synthetic_ndvi_inputs(12, 8)
    mktempdir() do directory
        path = joinpath(directory, "ndvi.csv")
        CSV.write(path, data.table)
        loaded = load_ndvi_covariates(path, data.ids)
        @test loaded.data_qc.station_count[1] == 12
        @test loaded.data_qc.nonland_station_count[1] == 1
        @test nrow(loaded.station_qc) == 12
        @test "S12" in loaded.station_qc.station_id[loaded.station_qc.nonland_pixel_flag]
    end
end

@testset "NDVI station-period aggregation and fold isolation" begin
    data = synthetic_ndvi_inputs()
    cfg = NDVISelectionConfig(
        outdir="unused", min_wet_hours_per_period=5, min_periods_per_station=8,
        min_stations_per_period=10, association_permutations=39,
        spatial_permutations=0, seed=55,
    )
    train = collect(1:30)
    panel = prepare_ndvi_panel(
        data.Yobs, data.Ysat, data.aligned, train, data.ids, cfg,
    )
    @test panel.qc.status[1] == "ok"
    @test nrow(filter(:used => identity, panel.cells)) <= length(train) * length(data.periods)
    @test all(panel.cells.wet_hour_count .>= 1)
    first_screen = screen_ndvi(panel, data.aligned.periods, cfg; rng=MersenneTwister(6))

    changed_obs = copy(data.Yobs); changed_obs[31:40, :] .= 1.0e7
    changed_sat = copy(data.Ysat); changed_sat[31:40, :] .= -1.0e7
    changed_values = copy(data.aligned.values); changed_values[31:40, :] .= 0.99
    changed_aligned = merge(data.aligned, (; values=changed_values))
    changed_panel = prepare_ndvi_panel(
        changed_obs, changed_sat, changed_aligned, train, data.ids, cfg,
    )
    second_screen = screen_ndvi(changed_panel, data.aligned.periods, cfg; rng=MersenneTwister(6))
    @test first_screen.association == second_screen.association
    @test first_screen.vif.status[1] == "not_applicable_single_predictor"
end

@testset "NDVI global and local spatial roles" begin
    rng = MersenneTwister(301)
    n, periods = 45, 14
    lon = collect(range(108.0, 112.0; length=n))
    lat = 30 .+ 0.6 .* sin.(range(0, 3pi; length=n))
    lonlat = hcat(lon, lat)
    x = randn(rng, n, periods)
    cfg = NDVISelectionConfig(
        outdir="unused", min_stations_per_period=12,
        bandwidth_candidates=[30], association_permutations=0,
        spatial_permutations=49, seed=302,
    )
    for (local_effect, expected) in ((false, "global"), (true, "local"))
        coefficient = local_effect ? 0.5 .+ 5 .* (lon .> mean(lon)) : fill(2.0, n)
        y = coefficient .* x .+ 0.03 .* randn(rng, n, periods)
        panel = (; y, x=Dict(:ndvi => x), mask=trues(n, periods), eligible_station=trues(n))
        result = ndvi_spatial_variability_test(
            panel, [:ndvi], lonlat, cfg; rng=MersenneTwister(303),
        )
        @test result.variability.role[1] == expected
    end
end

@testset "NDVI independent output contract" begin
    data = synthetic_ndvi_inputs()
    mktempdir() do outdir
        cfg = NDVISelectionConfig(
            outdir=outdir, min_wet_hours_per_period=5, min_periods_per_station=8,
            min_stations_per_period=10, k=5, bandwidth_candidates=[30],
            association_permutations=19, spatial_permutations=9, seed=404,
        )
        products = ["FY4B", "GPM", "GSMaP"]
        satellite = Dict(product => data.Ysat for product in products)
        data_qc = DataFrame(station_count=[40], all_237_stations=[true], complete=[true])
        station_qc = DataFrame(
            station_id=data.ids, nonland_pixel_flag=[fill(false, 39); true],
        )
        result = run_ndvi_variable_selection(
            cfg, products, data.ids, data.times, data.Yobs, satellite,
            data.table, data.lonlat; data_qc, station_qc,
        )
        required = [
            "ndvi_data_qc.csv", "ndvi_alignment_qc.csv", "ndvi_alignment_intervals.csv",
            "ndvi_fold_association.csv",
            "ndvi_fold_roles.csv", "ndvi_final_full_data_spec.csv",
            "ndvi_cross_product_consensus.csv", "run_status.csv", "spatial_folds.csv",
        ]
        @test all(isfile(joinpath(outdir, filename)) for filename in required)
        @test nrow(result.specifications) == 3
        @test result.consensus.variable[1] == "ndvi"
    end
end
