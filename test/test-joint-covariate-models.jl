using Test, CSV, DataFrames, Random

if !isdefined(Main, :JointCovariateModels)
    include(joinpath(@__DIR__, "..", "src", "JointCovariateModels.jl"))
end
const JCM = Main.JointCovariateModels

function joint_model_fixture(; seed=91)
    rng = MersenneTwister(seed)
    n, nt = 32, 5
    train, target = collect(1:25), collect(26:32)
    lonlat = hcat(
        collect(range(109.0, 113.0; length=n)),
        30 .+ 0.7 .* sin.(range(0, 3pi; length=n)),
    )
    terrain = DataFrame(
        station_id=string.(1:n), elevation_m=100 .+ 1400 .* rand(rng, n),
        slope_deg=rand(rng, n) .* 20, aspect_sin=randn(rng, n),
        aspect_cos=randn(rng, n),
    )
    era5 = Dict(variable => randn(rng, n, nt)
        for variable in [:t2m_c, :d2m_c, :u10, :v10, :sp_hpa])
    era5[:u10][:, 2] .+= collect(range(-2.0, 2.0; length=n))
    residual = 0.4 .* era5[:u10] .+ 0.001 .* terrain.elevation_m
    Yobs = 5 .+ residual .+ 0.01 .* randn(rng, n, nt)
    Ysat = Yobs .- residual
    cfg = JCM.JointCovariateBenchmarkConfig(
        spec_path="unused", terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[10, 15, 20], max_iterations=200, tolerance=1e-5,
    )
    roles = Dict("elevation" => "global", "u10" => "local")
    context = JCM.build_joint_fold_context(
        "GPM", roles, train, target, lonlat, Yobs, Ysat,
        terrain, era5, nothing, cfg,
    )
    return (; context, roles, train, target, lonlat, terrain, era5, Yobs, Ysat, cfg)
end

@testset "joint covariate specification" begin
    mktempdir() do directory
        path = joinpath(directory, "spec.csv")
        CSV.write(path, DataFrame([
            (product="FY4B", family="era5", variable_group="u10",
                predictor_columns="u10", independent_selected=true,
                joint_vif_retained=true, role="local", final_included=true),
            (product="GPM", family="dem", variable_group="elevation",
                predictor_columns="elevation_m", independent_selected=true,
                joint_vif_retained=true, role="global", final_included=true),
            (product="GSMaP", family="era5", variable_group="sp_hpa",
                predictor_columns="sp_hpa", independent_selected=true,
                joint_vif_retained=true, role="local", final_included=true),
        ]))
        loaded = JCM.load_joint_covariate_spec(path, ["FY4B", "GPM", "GSMaP"])
        @test nrow(loaded.included) == 3
        @test loaded.role_maps["GPM"] == Dict("elevation" => "global")
        @test length(JCM.joint_spec_sha256(path)) == 64
    end
end

@testset "dynamic joint designs and roles" begin
    fixture = joint_model_fixture()
    context = fixture.context
    @test context.variables == ["elevation", "u10"]
    @test JCM.joint_effective_roles(context, "residual_gwr") ==
        Dict("elevation" => "local", "u10" => "local")
    @test JCM.joint_effective_roles(context, "mixed_gwr") == fixture.roles
    # The default layout drops the coordinate columns; `:split` is exercised by the grouping
    # testset below.
    @test JCM.joint_group_names(context, "mgwr") == ["intercept", "u10"]
    @test context.predictor_train["elevation"][:, :, 1] ≈
        context.predictor_train["elevation"][:, :, 5]
    @test !(context.predictor_train["u10"][:, :, 1] ≈
        context.predictor_train["u10"][:, :, 2])

    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    prediction, converged = JCM.dynamic_covariate_predict(
        context, residuals, "residual_gwr", [15],
    )
    @test size(prediction) == (length(fixture.target), size(fixture.Yobs, 2))
    @test all(converged)
    @test all(isfinite, prediction)
    mgwr_prediction, mgwr_converged = JCM.dynamic_covariate_predict(
        context, residuals, "mgwr", fill(20, 2),
    )
    @test count(mgwr_converged) >= 4
    @test count(isfinite, mgwr_prediction) >= 4length(fixture.target)
end

"""Rebuild the fixture's context under a different MGWR spatial grouping."""
function regrouped_context(fixture, grouping::Symbol)
    cfg = JCM.JointCovariateBenchmarkConfig(
        spec_path="unused", terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[10, 15, 20], max_iterations=200, tolerance=1e-5,
        mgwr_spatial_grouping=grouping,
    )
    return JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        fixture.Yobs, fixture.Ysat, fixture.terrain, fixture.era5, nothing, cfg,
    )
end

@testset "mgwr spatial grouping controls only the spatial block" begin
    fixture = joint_model_fixture()
    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    expected = Dict(
        :split => ["intercept", "longitude", "latitude", "u10"],
        :shared => ["shared_spatial", "u10"],
        :intercept_only => ["intercept", "u10"],
    )
    for (grouping, names) in expected
        context = regrouped_context(fixture, grouping)
        @test JCM.joint_group_names(context, "mgwr") == names
        # The covariate always keeps its own group, which is what makes the model multiscale.
        @test last(names) == "u10"
        prediction, converged = JCM.dynamic_covariate_predict(
            context, residuals, "mgwr", fill(20, length(names)),
        )
        @test size(prediction) == (length(fixture.target), size(fixture.Yobs, 2))
        @test all(converged)
        @test all(isfinite, prediction)
    end
    # The non-mgwr methods keep one shared spatial group whatever the setting says.
    for grouping in (:split, :shared, :intercept_only)
        context = regrouped_context(fixture, grouping)
        @test JCM.joint_group_names(context, "mixed_gwr") == ["shared_local"]
        @test JCM.joint_group_names(context, "residual_gwr") == ["shared_all_local"]
    end
    @test_throws ArgumentError JCM.joint_group_names(
        regrouped_context(fixture, :nonsense), "mgwr",
    )
end

@testset "back-fitting rejects a bandwidth vector of the wrong length" begin
    fixture = joint_model_fixture()
    n = length(fixture.train)
    lonlat = fixture.lonlat[fixture.train, :]
    groups = [ones(n, 1), reshape(collect(range(-1.0, 1.0; length=n)), :, 1)]
    response = reshape(collect(range(0.0, 1.0; length=n)), :, 1)
    empty_global = zeros(n, 0)
    # Two groups, three bandwidths: the loop used to `zip` these and silently drop the extra,
    # producing a wrong fit rather than an error.
    @test_throws DimensionMismatch JCM._multiscale_predict_damped(
        groups, empty_global, response, lonlat, groups, empty_global, lonlat,
        [10, 10, 10], fixture.cfg,
    )
    prediction, ok = JCM._multiscale_predict_damped(
        groups, empty_global, response, lonlat, groups, empty_global, lonlat,
        [10, 10], fixture.cfg,
    )
    @test only(ok)
    @test all(isfinite, prediction)
end

@testset "relaxation is configurable and reaches the same fixed point" begin
    fixture = joint_model_fixture()
    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    @test fixture.cfg.relaxation == 1.5
    # Over-relaxation changes the convergence rate, not the solution, so a run at a different
    # admissible factor must agree with the default to within the convergence tolerance.
    # Pinned to `:split`, the layout where the relaxation factor actually decides whether the
    # back-fit converges at all.
    relaxed = JCM.JointCovariateBenchmarkConfig(
        spec_path="unused", terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[10, 15, 20], max_iterations=4000, tolerance=1e-10,
        relaxation=1.0, mgwr_spatial_grouping=:split,
    )
    context = JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        fixture.Yobs, fixture.Ysat, fixture.terrain, fixture.era5, nothing, relaxed,
    )
    tight = JCM.JointCovariateBenchmarkConfig(
        spec_path="unused", terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[10, 15, 20], max_iterations=4000, tolerance=1e-10,
        relaxation=1.5, mgwr_spatial_grouping=:split,
    )
    tight_context = JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        fixture.Yobs, fixture.Ysat, fixture.terrain, fixture.era5, nothing, tight,
    )
    slow, slow_ok = JCM.dynamic_covariate_predict(context, residuals, "mgwr", fill(20, 4))
    fast, fast_ok = JCM.dynamic_covariate_predict(tight_context, residuals, "mgwr", fill(20, 4))
    @test all(slow_ok) && all(fast_ok)
    @test slow ≈ fast atol = 1e-4
end

@testset "joint training-fold isolation and missing target" begin
    fixture = joint_model_fixture()
    changed_obs = copy(fixture.Yobs); changed_sat = copy(fixture.Ysat)
    changed_obs[fixture.target, :] .= 1.0e6
    changed_sat[fixture.target, :] .= -1.0e6
    changed_era5 = Dict(variable => copy(values) for (variable, values) in fixture.era5)
    changed_era5[:u10][fixture.target, :] .+= 10.0
    changed = JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        changed_obs, changed_sat, fixture.terrain, changed_era5, nothing, fixture.cfg,
    )
    @test fixture.context.scaling == changed.scaling
    @test fixture.context.predictor_train == changed.predictor_train
    @test fixture.context.predictor_target != changed.predictor_target

    missing_era5 = Dict(variable => copy(values) for (variable, values) in fixture.era5)
    missing_era5[:u10][first(fixture.target), 1] = NaN
    missing_context = JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        fixture.Yobs, fixture.Ysat, fixture.terrain, missing_era5, nothing, fixture.cfg,
    )
    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    prediction, _ = JCM.dynamic_covariate_predict(
        missing_context, residuals, "mixed_gwr", [15],
    )
    @test isnan(prediction[1, 1])
end

@testset "global leave-one-out excludes held response" begin
    rng = MersenneTwister(19)
    X = hcat(randn(rng, 25), randn(rng, 25))
    y = randn(rng, 25)
    first_prediction = JCM._global_predict(X, y, X, 1e-8; leave_one_out=true)
    changed = copy(y); changed[7] += 1.0e6
    second_prediction = JCM._global_predict(X, changed, X, 1e-8; leave_one_out=true)
    @test first_prediction[7] ≈ second_prediction[7]
end
