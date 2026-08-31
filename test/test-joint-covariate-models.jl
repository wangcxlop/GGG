using Test, CSV, DataFrames, LinearAlgebra, Random

include(joinpath(@__DIR__, "..", "src", "load_modules.jl"))
load_standalone_modules("JointCovariateModels")
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

    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    prediction, converged = JCM.dynamic_covariate_predict(
        context, residuals, "residual_gwr", [15.0], bisquare,
    )
    @test size(prediction) == (length(fixture.target), size(fixture.Yobs, 2))
    @test all(converged)
    @test all(isfinite, prediction)
    mgwr_prediction, mgwr_converged = JCM.dynamic_covariate_predict(
        context, residuals, "mgwr", fill(20.0, 2), bisquare,
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
    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
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
            context, residuals, "mgwr", fill(20.0, length(names)), bisquare,
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
    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
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
        [10.0, 10.0, 10.0], bisquare, fixture.cfg,
    )
    prediction, ok = JCM._multiscale_predict_damped(
        groups, empty_global, response, lonlat, groups, empty_global, lonlat,
        [10.0, 10.0], bisquare, fixture.cfg,
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
    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
    slow, slow_ok = JCM.dynamic_covariate_predict(
        context, residuals, "mgwr", fill(20.0, 4), bisquare,
    )
    fast, fast_ok = JCM.dynamic_covariate_predict(
        tight_context, residuals, "mgwr", fill(20.0, 4), bisquare,
    )
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
        missing_context, residuals, "mixed_gwr", [15.0], JCM.DEMTerrainExperiment._bisquare_kernel,
    )
    @test isnan(prediction[1, 1])
end

"""Rebuild the fixture's context with a chosen `unsupported_local_target`."""
function unsupported_context(fixture, unsupported::Symbol)
    cfg = JCM.JointCovariateBenchmarkConfig(
        spec_path="unused", terrain_path="unused", era5_paths=Dict{Int,String}(),
        bandwidth_candidates=[10, 15, 20], max_iterations=200, tolerance=1e-5,
        unsupported_local_target=unsupported,
    )
    return JCM.build_joint_fold_context(
        "GPM", fixture.roles, fixture.train, fixture.target, fixture.lonlat,
        fixture.Yobs, fixture.Ysat, fixture.terrain, fixture.era5, nothing, cfg,
    )
end

@testset "an unfittable local target is missing, not a zero correction" begin
    # A target with fewer than `p + 1` positively weighted training stations cannot be fitted
    # locally. `:zero` (the historical encoding) left `_local_hat`'s all-zero row in place, so
    # the target contributed exactly 0 and was reported as a successful prediction of "no local
    # correction" — invisible to every coverage gate the benchmark has. `:missing` marks it NaN.
    #
    # Only a compact kernel at a *fixed* bandwidth can reach that branch: an adaptive bandwidth
    # always keeps `k - 1` stations, and `bw = Inf` gives every kernel weight 1.0. The fixture's
    # seven targets sit east of the training block, so a fixed 60 km reaches the first three and
    # not the last four.
    fixture = joint_model_fixture()
    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
    residuals = fixture.Yobs[fixture.train, :] .- fixture.Ysat[fixture.train, :]
    zero_context = unsupported_context(fixture, :zero)
    missing_context = unsupported_context(fixture, :missing)
    @test zero_context.config.unsupported_local_target === :zero
    # The benchmark's default has to be the honest one; `:zero` is the opt-in escape hatch.
    @test fixture.context.config.unsupported_local_target === :missing

    # MGWR, `:intercept_only`: the intercept group is one column, so it needs two supporting
    # stations. `Inf` on the covariate group keeps it out of the picture.
    bandwidths = [60.0, Inf]
    zero_mgwr, zero_ok = JCM.dynamic_covariate_predict(
        zero_context, residuals, "mgwr", bandwidths, bisquare; adaptive=false,
    )
    missing_mgwr, missing_ok = JCM.dynamic_covariate_predict(
        missing_context, residuals, "mgwr", bandwidths, bisquare; adaptive=false,
    )
    @test zero_ok == missing_ok
    @test any(zero_ok)
    supported, unsupported_rows = 1:3, 4:7
    for time in findall(zero_ok)
        # The bug: an unsupported target used to come back as a real-looking number.
        @test all(isfinite, zero_mgwr[unsupported_rows, time])
        @test all(isnan, missing_mgwr[unsupported_rows, time])
        # The fix must not move a value the local fit could actually support. `_local_hat`
        # fills one row per target and the back-fit hats are untouched, so these are equal
        # bit for bit, not merely close.
        @test isequal(
            collect(zero_mgwr[supported, time]), collect(missing_mgwr[supported, time]),
        )
    end

    # The other joint path, `mixed_gwr_predict`. Its single shared block carries the three
    # spatial columns plus the local covariate, so at 60 km no target clears `p + 1 = 5` and
    # every one of them was previously a fabricated zero correction.
    zero_mixed, _ = JCM.dynamic_covariate_predict(
        zero_context, residuals, "mixed_gwr", [60.0], bisquare; adaptive=false,
    )
    missing_mixed, _ = JCM.dynamic_covariate_predict(
        missing_context, residuals, "mixed_gwr", [60.0], bisquare; adaptive=false,
    )
    @test all(isfinite, zero_mixed)
    @test all(isnan, missing_mixed)

    # An adaptive bandwidth cannot reach the branch, so the two settings agree everywhere.
    for method in ("residual_gwr", "mixed_gwr")
        zero_adaptive, _ = JCM.dynamic_covariate_predict(
            zero_context, residuals, method, [15.0], bisquare,
        )
        missing_adaptive, _ = JCM.dynamic_covariate_predict(
            missing_context, residuals, method, [15.0], bisquare,
        )
        @test all(isfinite, missing_adaptive)
        @test isequal(zero_adaptive, missing_adaptive)
    end
end

@testset "residual_gwr and mixed_gwr coincide when no covariate is global" begin
    # `residual_gwr` forces every covariate local; `mixed_gwr` honours the role map. With no
    # "global" in that map the two designs are the same matrix and only the group name differs,
    # so the benchmark reports one model under two names. The fixture has elevation global, so
    # it is the distinct case; dropping that role gives the degenerate one.
    distinct = joint_model_fixture()
    @test !JCM.joint_models_coincide(distinct.context, "mixed_gwr", "residual_gwr")

    all_local = JCM.build_joint_fold_context(
        "GPM", Dict("elevation" => "local", "u10" => "local"),
        distinct.train, distinct.target, distinct.lonlat, distinct.Yobs, distinct.Ysat,
        distinct.terrain, distinct.era5, nothing, distinct.cfg,
    )
    @test JCM.joint_models_coincide(all_local, "mixed_gwr", "residual_gwr")
    # The claim the predicate is standing in for: same local design, same global design, same
    # predictions — only the group name differs.
    residual_design = JCM._design_at(all_local, "residual_gwr", 1)
    mixed_design = JCM._design_at(all_local, "mixed_gwr", 1)
    @test residual_design.local_groups == mixed_design.local_groups
    @test residual_design.global_design == mixed_design.global_design
    @test residual_design.group_names != mixed_design.group_names

    bisquare = JCM.DEMTerrainExperiment._bisquare_kernel
    residuals = distinct.Yobs[distinct.train, :] .- distinct.Ysat[distinct.train, :]
    same = [JCM.dynamic_covariate_predict(all_local, residuals, method, [15.0], bisquare)[1]
        for method in ("residual_gwr", "mixed_gwr")]
    @test isequal(same[1], same[2])

    # `mgwr` is never a duplicate: per-group bandwidths make it a different model whatever the
    # roles say, and `:intercept_only` also drops the coordinate columns.
    @test !JCM.joint_models_coincide(all_local, "mgwr", "residual_gwr")
    @test !JCM.joint_models_coincide(all_local, "mgwr", "mixed_gwr")
    @test JCM.joint_models_coincide(all_local, "mgwr", "mgwr")
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

@testset "local weighted solve returns the ridge WLS coefficients" begin
    rng = MersenneTwister(37)
    n, m, ridge = 20, 4, 1e-8
    Xtrain = hcat(ones(n), randn(rng, n), randn(rng, n))
    Xtarget = hcat(ones(m), randn(rng, m), randn(rng, m))
    y = randn(rng, n)
    # Every weight is positive, so each target selects the whole training set and the
    # prediction must reproduce the closed-form weighted ridge solution exactly.
    weights = rand(rng, m, n) .+ 0.5
    prediction = JCM._local_predict(Xtrain, y, Xtarget, weights, ridge)
    for target in 1:m
        w = weights[target, :]
        beta = (Xtrain' * (w .* Xtrain) + ridge * I) \ (Xtrain' * (w .* y))
        @test prediction[target] ≈ dot(Xtarget[target, :], beta)
    end
end
