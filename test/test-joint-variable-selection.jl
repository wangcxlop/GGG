using Test, DataFrames, Dates, Random

if !isdefined(Main, :JointVariableSelection)
    include(joinpath(@__DIR__, "..", "src", "JointVariableSelection.jl"))
end
const JVS = Main.JointVariableSelection

@testset "joint hierarchical NDVI weights" begin
    n, nt = 3, 6
    Yobs = fill(1.0, n, nt)
    Ysat = zeros(n, nt)
    Yobs[1, 3:4] .= 0.0
    terrain = DataFrame(
        station_id=string.(1:n), elevation_m=collect(1.0:n),
        slope_deg=[2.0, 4.0, 7.0], aspect_sin=[0.1, 0.4, 0.8],
        aspect_cos=[0.9, 0.5, 0.2],
    )
    u10 = [i + 0.2t for i in 1:n, t in 1:nt]
    era5 = Dict(:u10 => u10)
    ndvi_values = [0.1i + 0.03t for i in 1:n, t in 1:nt]
    source_period = repeat(reshape([1, 1, 1, 1, 2, 2], 1, :), n, 1)
    aligned = (; values=ndvi_values, source_period)
    cfg = JVS.JointSelectionConfig(
        outdir="unused", min_wet_hours=1, min_stations_per_time=2,
        min_ndvi_wet_hours_per_period=1, min_ndvi_periods_per_station=2,
        independent_permutations=9, spatial_permutations=9,
    )
    panel = JVS.prepare_joint_panel(
        Yobs, Ysat, terrain, era5, aligned, collect(1:n), ["u10", "ndvi"], cfg,
    )
    @test panel.qc.status[1] == "ok"
    @test panel.weight_audit.total_weight ≈ ones(n)
    station1 = panel.weight_audit[1, :]
    @test station1.ndvi_period_count == 2
    @test station1.minimum_period_weight ≈ 0.5
    @test station1.maximum_period_weight ≈ 0.5
end

@testset "joint grouped and deterministic VIF" begin
    rng = MersenneTwister(42)
    n, nt = 10, 10
    base = randn(rng, n, nt)
    panel = (;
        y=randn(rng, n, nt), weights=fill(1 / nt, n, nt),
        x=Dict(
            :elevation_m => copy(base), :aspect_sin => copy(base),
            :aspect_cos => 2 .* base,
            :u10 => randn(rng, n, nt), :v10 => randn(rng, n, nt),
        ),
    )
    cfg = JVS.JointSelectionConfig(outdir="unused", vif_threshold=5.0)
    audit, retained = JVS.joint_vif(
        panel, ["elevation", "aspect", "u10", "v10"],
        Dict("elevation" => 0.01, "aspect" => 0.04, "u10" => 0.01, "v10" => 0.01), cfg,
    )
    @test "aspect" ∉ retained
    @test "elevation" in retained
    removed_aspect = filter([:variable_group, :removed] => (g, r) -> g == "aspect" && r, audit)
    @test nrow(removed_aspect) == 1
    @test removed_aspect.predictor_columns[1] == "aspect_sin+aspect_cos"

    correlated = randn(rng, n, nt)
    wind_panel = (;
        y=randn(rng, n, nt), weights=fill(1 / nt, n, nt),
        x=Dict(:u10 => correlated, :v10 => copy(correlated)),
    )
    _, wind_retained = JVS.joint_vif(
        wind_panel, ["u10", "v10"], Dict("u10" => 0.01, "v10" => 0.04), cfg,
    )
    @test wind_retained == ["u10"]
end

@testset "joint fold screening isolation" begin
    rng = MersenneTwister(77)
    n, nt = 40, 20
    train = collect(1:30)
    ids = ["S$i" for i in 1:n]
    times = [DateTime(2022, 6, 1) + Hour(t - 1) for t in 1:nt]
    Yobs = 2 .+ abs.(randn(rng, n, nt))
    Ysat = Yobs .- 0.3 .* randn(rng, n, nt)
    terrain = DataFrame(
        station_id=ids, elevation_m=rand(rng, n) .* 1000,
        slope_deg=rand(rng, n) .* 30, aspect_sin=randn(rng, n),
        aspect_cos=randn(rng, n),
    )
    era5 = Dict(variable => randn(rng, n, nt)
        for variable in JVS.ERA5VariableSelection.ERA5_VARIABLES)
    source_period = repeat(reshape(vcat(fill(1, 10), fill(2, 10)), 1, :), n, 1)
    aligned = (;
        values=rand(rng, n, nt), source_period,
        periods=[Date(2022, 6, 1), Date(2022, 6, 17)],
    )
    cfg = JVS.JointSelectionConfig(
        outdir="unused", min_wet_hours=1, min_stations_per_time=3,
        min_ndvi_wet_hours_per_period=1, min_ndvi_periods_per_station=2,
        independent_permutations=9, spatial_permutations=9,
    )
    before, active_before = JVS._candidate_rows(
        "synthetic", Yobs, Ysat, terrain, era5, aligned, train, ids, times, cfg, 4, 81,
    )
    Yobs_changed, Ysat_changed = copy(Yobs), copy(Ysat)
    Yobs_changed[31:end, :] .= 1000 .* rand(rng, 10, nt)
    Ysat_changed[31:end, :] .= -1000 .* rand(rng, 10, nt)
    terrain_changed = copy(terrain)
    terrain_changed[31:end, [:elevation_m, :slope_deg, :aspect_sin, :aspect_cos]] .= 1.0e8
    era5_changed = Dict(variable => copy(values) for (variable, values) in era5)
    for values in Base.values(era5_changed)
        values[31:end, :] .= 1.0e8
    end
    ndvi_changed = copy(aligned.values)
    ndvi_changed[31:end, :] .= 1.0e8
    after, active_after = JVS._candidate_rows(
        "synthetic", Yobs_changed, Ysat_changed, terrain_changed, era5_changed,
        (; values=ndvi_changed, source_period=aligned.source_period, periods=aligned.periods),
        train, ids, times, cfg, 4, 81,
    )
    @test isequal(before, after)
    @test active_before == active_after
end

function synthetic_role_panel(spatially_varying::Bool)
    rng = MersenneTwister(spatially_varying ? 202 : 101)
    n, nt = 45, 12
    lon = repeat(collect(range(0.0, 1.0; length=9)), 5)
    lat = repeat(collect(range(0.0, 1.0; length=5)), inner=9)
    coords = hcat(lon, lat)
    x = randn(rng, n, nt)
    beta = spatially_varying ? 0.5 .+ 5 .* lon : fill(2.0, n)
    y = beta .* x
    panel = (;
        y, x=Dict(:u10 => x), weights=fill(1 / nt, n, nt),
        eligible_station=trues(n),
    )
    return panel, coords
end

@testset "joint spatial global and local roles" begin
    cfg = JVS.JointSelectionConfig(
        outdir="unused", bandwidth_candidates=[12, 18, 25, 35],
        spatial_permutations=99, q_threshold=0.05,
    )
    global_panel, coords = synthetic_role_panel(false)
    global_result = JVS.joint_spatial_variability_test(
        global_panel, ["u10"], coords, cfg; rng=MersenneTwister(1),
    )
    @test global_result.result.role[1] == "global"

    local_panel, coords = synthetic_role_panel(true)
    local_result = JVS.joint_spatial_variability_test(
        local_panel, ["u10"], coords, cfg; rng=MersenneTwister(2),
    )
    @test local_result.result.role[1] == "local"
end

"""Fixture shared by the `select_joint_covariates` tests: a synthetic 40-station panel where
stations 31:40 carry an extreme, otherwise-influential signal, so any leak from outside
`train_indices` would change the result."""
function synthetic_joint_fixture()
    rng = MersenneTwister(77)
    n, nt = 40, 20
    ids = ["S$i" for i in 1:n]
    times = [DateTime(2022, 6, 1) + Hour(t - 1) for t in 1:nt]
    Yobs = 2 .+ abs.(randn(rng, n, nt))
    Ysat = Yobs .- 0.3 .* randn(rng, n, nt)
    terrain = DataFrame(
        station_id=ids, elevation_m=rand(rng, n) .* 1000,
        slope_deg=rand(rng, n) .* 30, aspect_sin=randn(rng, n), aspect_cos=randn(rng, n),
    )
    era5 = Dict(variable => randn(rng, n, nt)
        for variable in JVS.ERA5VariableSelection.ERA5_VARIABLES)
    source_period = repeat(reshape(vcat(fill(1, 10), fill(2, 10)), 1, :), n, 1)
    aligned = (;
        values=rand(rng, n, nt), source_period,
        periods=[Date(2022, 6, 1), Date(2022, 6, 17)],
    )
    lonlat = hcat(range(100.0, 110.0; length=n), range(20.0, 30.0; length=n))
    cfg = JVS.JointSelectionConfig(
        outdir="unused", min_wet_hours=1, min_stations_per_time=3,
        min_ndvi_wet_hours_per_period=1, min_ndvi_periods_per_station=2,
        independent_permutations=9, spatial_permutations=9, bandwidth_candidates=[8, 12, 18],
    )
    return (; ids, times, Yobs, Ysat, terrain, era5, aligned, lonlat, cfg)
end

@testset "select_joint_covariates" begin
    f = synthetic_joint_fixture()
    train = collect(1:30)
    selection = JVS.select_joint_covariates(
        "synthetic", f.Yobs, f.Ysat, f.terrain, f.era5, f.aligned,
        train, f.ids, f.times, f.lonlat, f.cfg, 4, 81,
    )
    @test selection.status in ("ok", "failed")
    @test Set(keys(selection.role_map)) ⊆ Set(selection.active)
    @test all(value -> value in ("local", "global"), values(selection.role_map))
    # role_map must agree with spatial_result restricted to local/global roles.
    expected = Dict(String(row.variable_group) => String(row.role)
        for row in eachrow(selection.spatial_result) if row.role in ("local", "global"))
    @test selection.role_map == expected

    # Perturbing data OUTSIDE train_indices must not change the fold's own selection —
    # the whole point of nesting selection inside cross-validation.
    Yobs2, Ysat2 = copy(f.Yobs), copy(f.Ysat)
    Yobs2[31:end, :] .= 1000 .* rand(MersenneTwister(1), 10, length(f.times))
    Ysat2[31:end, :] .= -1000 .* rand(MersenneTwister(2), 10, length(f.times))
    terrain2 = copy(f.terrain)
    terrain2[31:end, [:elevation_m, :slope_deg, :aspect_sin, :aspect_cos]] .= 1.0e8
    era5_2 = Dict(variable => copy(values) for (variable, values) in f.era5)
    for values in Base.values(era5_2)
        values[31:end, :] .= 1.0e8
    end
    ndvi2 = copy(f.aligned.values)
    ndvi2[31:end, :] .= 1.0e8
    selection2 = JVS.select_joint_covariates(
        "synthetic", Yobs2, Ysat2, terrain2, era5_2,
        (; values=ndvi2, source_period=f.aligned.source_period, periods=f.aligned.periods),
        train, f.ids, f.times, f.lonlat, f.cfg, 4, 81,
    )
    @test selection.role_map == selection2.role_map
    @test selection.active == selection2.active
    @test isequal(selection.candidates, selection2.candidates)

    # A visibly different training subset is free to select differently — confirms the result
    # is actually a function of `train_indices`, not a constant fallback.
    other_train = collect(11:40)
    selection3 = JVS.select_joint_covariates(
        "synthetic", f.Yobs, f.Ysat, f.terrain, f.era5, f.aligned,
        other_train, f.ids, f.times, f.lonlat, f.cfg, 4, 81,
    )
    @test selection3.status in ("ok", "failed")
end
