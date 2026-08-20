using Test
using CSV, DataFrames, Dates, Random, Statistics

if !isdefined(Main, :ERA5VariableSelection)
    include(joinpath(@__DIR__, "..", "src", "ERA5VariableSelection.jl"))
end
using .ERA5VariableSelection

function synthetic_era5_panel(n::Int, nt::Int; local_effect::Bool=false, seed::Int=41)
    rng = MersenneTwister(seed)
    lon = collect(range(108.0, 112.0; length=n))
    lat = 30 .+ 0.8 .* sin.(range(0, 3pi; length=n))
    lonlat = hcat(lon, lat)
    times = [DateTime(2022, 6, 1) + Hour(t - 1) for t in 1:nt]
    spatial = (lon .- mean(lon)) ./ std(lon)
    temporal = sin.(range(0, 6pi; length=nt))
    values = Dict{Symbol,Matrix{Float64}}()
    values[:t2m_c] = spatial .+ temporal' .+ 0.15 .* randn(rng, n, nt)
    values[:d2m_c] = values[:t2m_c] .+ 0.001 .* randn(rng, n, nt)
    values[:u10] = randn(rng, n, nt)
    values[:v10] = randn(rng, n, nt)
    values[:sp_hpa] = randn(rng, n, nt)
    coefficient = local_effect ? 0.5 .+ 4 .* (spatial .> 0) : fill(2.0, n)
    residual = coefficient .* values[:t2m_c] .+ 0.05 .* randn(rng, n, nt)
    Yobs = 1 .+ residual .- minimum(residual)
    Ysat = Yobs .- residual
    return (; Yobs, Ysat, values, lonlat, times)
end

@testset "ERA5 feature-time alignment and integrity" begin
    mktempdir() do directory
        path = joinpath(directory, "era5_2022.csv")
        utc = DateTime(2022, 6, 1, 0)
        frame = DataFrame(
            station_id=["A", "B"], time_utc=[utc, utc],
            # This is deliberately UTC+8; matching must nevertheless use UTC+9.
            time_bjt=[utc + Hour(8), utc + Hour(8)],
            t2m_c=[20.0, 21.0], d2m_c=[18.0, 18.5], u10=[1.0, 2.0],
            v10=[-1.0, -2.0], sp_hpa=[1000.0, 999.0],
        )
        CSV.write(path, frame)
        loaded = load_era5_panel(Dict(2022 => path), ["A", "B"], [utc + Hour(9)];
            validate_annual_complete=false)
        @test loaded.values[:t2m_c][:, 1] == [20.0, 21.0]
        @test loaded.qc.feature_time_offset_hours[1] == 9
        @test loaded.qc.complete[1]

        CSV.write(path, vcat(frame, frame[1:1, :]))
        @test_throws ErrorException load_era5_panel(
            Dict(2022 => path), ["A", "B"], [utc + Hour(9)]; validate_annual_complete=false,
        )
    end
end

@testset "Station blocks and deterministic VIF" begin
    rng = MersenneTwister(7)
    order = station_block_permutation(8, rng)
    panel = reshape(collect(1:32), 8, 4)
    @test sort(order) == collect(1:8)
    @test all(panel[order[i], :] == panel[order, :][i, :] for i in 1:8)

    data = synthetic_era5_panel(40, 48)
    cfg = ERA5SelectionConfig(
        outdir="unused", min_wet_hours=10, min_stations_per_time=12,
        association_permutations=99, spatial_permutations=0,
        bandwidth_candidates=[30], seed=12,
    )
    prepared = prepare_dynamic_panel(
        data.Yobs, data.Ysat, data.values, collect(1:40), data.times, cfg,
    )
    screened = dynamic_panel_screen(prepared, data.times, cfg)
    @test :t2m_c in screened.selected
    @test !(:d2m_c in screened.selected)
    removed = filter(:removed => identity, screened.vif)
    @test "d2m_c" in removed.variable
    # u10 and v10 are separate rows and never treated as a vector group.
    @test count(==("u10"), screened.association.variable) == 1
    @test count(==("v10"), screened.association.variable) == 1

    wind_data = synthetic_era5_panel(40, 48; seed=88)
    wind_residual = 4 .* wind_data.values[:u10] .+ 0.03 .* randn(MersenneTwister(89), 40, 48)
    wind_obs = 1 .+ wind_residual .- minimum(wind_residual)
    wind_sat = wind_obs .- wind_residual
    wind_cfg = ERA5SelectionConfig(
        outdir="unused", min_wet_hours=10, min_stations_per_time=12,
        association_permutations=199, spatial_permutations=0, seed=90,
    )
    wind_panel = prepare_dynamic_panel(wind_obs, wind_sat, wind_data.values,
        collect(1:40), wind_data.times, wind_cfg)
    wind_screen = dynamic_panel_screen(wind_panel, wind_data.times, wind_cfg)
    @test :u10 in wind_screen.selected
    @test !(:v10 in wind_screen.selected)
end

@testset "Training-fold isolation" begin
    data = synthetic_era5_panel(40, 36)
    cfg = ERA5SelectionConfig(
        outdir="unused", min_wet_hours=8, min_stations_per_time=12,
        association_permutations=39, spatial_permutations=0, seed=99,
    )
    train = collect(1:30)
    first_panel = prepare_dynamic_panel(data.Yobs, data.Ysat, data.values, train, data.times, cfg)
    first = dynamic_panel_screen(first_panel, data.times, cfg; rng=MersenneTwister(10))
    changed_obs = copy(data.Yobs); changed_obs[31:40, :] .= 1.0e6
    changed_sat = copy(data.Ysat); changed_sat[31:40, :] .= -1.0e6
    changed_era5 = Dict(key => copy(value) for (key, value) in data.values)
    changed_era5[:u10][31:40, :] .= 1.0e6
    second_panel = prepare_dynamic_panel(changed_obs, changed_sat, changed_era5, train, data.times, cfg)
    second = dynamic_panel_screen(second_panel, data.times, cfg; rng=MersenneTwister(10))
    @test first.association == second.association
    @test first.vif == second.vif
end

@testset "Global and spatially varying panel effects" begin
    cfg = ERA5SelectionConfig(
        outdir="unused", min_wet_hours=10, min_stations_per_time=12,
        association_permutations=0, spatial_permutations=49,
        bandwidth_candidates=[30], seed=123,
    )
    global_data = synthetic_era5_panel(45, 40; local_effect=false)
    global_panel = prepare_dynamic_panel(global_data.Yobs, global_data.Ysat,
        global_data.values, collect(1:45), global_data.times, cfg)
    global_result = panel_spatial_variability_test(
        global_panel, [:t2m_c], global_data.lonlat, cfg; rng=MersenneTwister(22),
    )
    @test global_result.variability.role[1] == "global"

    local_data = synthetic_era5_panel(45, 40; local_effect=true)
    local_panel = prepare_dynamic_panel(local_data.Yobs, local_data.Ysat,
        local_data.values, collect(1:45), local_data.times, cfg)
    local_result = panel_spatial_variability_test(
        local_panel, [:t2m_c], local_data.lonlat, cfg; rng=MersenneTwister(22),
    )
    @test local_result.variability.role[1] == "local"
end

@testset "Independent ERA5 output contract" begin
    mktempdir() do outdir
        data = synthetic_era5_panel(40, 28)
        cfg = ERA5SelectionConfig(
            outdir=outdir, min_wet_hours=5, min_stations_per_time=10, k=5,
            association_permutations=19, spatial_permutations=9,
            bandwidth_candidates=[30], seed=2026,
        )
        products = ["FY4B", "GPM", "GSMaP"]
        satellite = Dict(product => data.Ysat for product in products)
        qc = DataFrame(complete=[true], feature_time_offset_hours=[9], duplicate_keys=[0],
            missing_station_hours=[0], nonfinite_values=[0])
        result = run_era5_variable_selection(cfg, products,
            ["S$(i)" for i in 1:40], data.times, data.Yobs, satellite,
            data.values, data.lonlat; data_qc=qc)
        required = [
            "era5_data_qc.csv", "era5_fold_association.csv", "era5_fold_roles.csv",
            "era5_role_stability.csv", "era5_final_full_data_spec.csv",
            "era5_cross_product_consensus.csv", "run_status.csv", "spatial_folds.csv",
        ]
        @test all(isfile(joinpath(outdir, file)) for file in required)
        @test nrow(result.specifications) == 15
        @test Set(result.specifications.variable) == Set(String.(ERA5_VARIABLES))
    end
end
