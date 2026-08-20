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
    @test JCM.joint_group_names(context, "mgwr") ==
        ["intercept", "longitude", "latitude", "u10"]
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
        context, residuals, "mgwr", fill(20, 4),
    )
    @test count(mgwr_converged) >= 4
    @test count(isfinite, mgwr_prediction) >= 4length(fixture.target)
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
