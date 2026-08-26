using Test
using DataFrames
using LinearAlgebra
using Random
using Statistics

if !isdefined(@__MODULE__, :DEMTerrainExperiment)
    include(joinpath(@__DIR__, "..", "src", "DEMTerrainExperiment.jl"))
end
using .DEMTerrainExperiment

@testset "DEM terrain screening helpers" begin
    @test bh_adjust([0.01, 0.04, 0.03]) ≈ [0.03, 0.04, 0.04]

    rng = MersenneTwister(42)
    n = 80
    elevation = collect(range(100.0, 1800.0; length=n))
    slope = rand(rng, n) .* 40
    aspect = collect(range(1.0, 359.0; length=n))
    terrain = DataFrame(
        station_id=string.(1:n), elevation_m=elevation, slope_deg=slope,
        aspect_deg=aspect, aspect_sin=sind.(aspect), aspect_cos=cosd.(aspect),
    )
    response = 0.01 .* elevation .+ 0.05 .* randn(rng, n)
    screen, vif, selected = terrain_screen(
        terrain, response; product="synthetic", permutations=99, seed=7,
    )
    @test "elevation" in selected
    @test screen.selected[screen.variable_group .== "elevation"] == [true]
    @test all(isfinite, vif.VIF)

    constant_screen, _, constant_selected = terrain_screen(
        terrain, ones(n); product="constant", permutations=19, seed=8,
    )
    @test isempty(constant_selected)
    @test !any(constant_screen.selected)

    collinear = copy(terrain)
    collinear.slope_deg = collinear.elevation_m ./ 50
    _, collinear_vif, collinear_selected = terrain_screen(
        collinear, response; product="collinear", permutations=49, seed=9,
    )
    @test length(intersect(collinear_selected, ["elevation", "slope"])) <= 1
    @test any(collinear_vif.VIF .>= 5)

    near_north = [359.0, 1.0]
    @test maximum(abs.(sind.(near_north))) < 0.02
    @test maximum(abs.(cosd.(near_north) .- 1.0)) < 0.001

    target = terrain[1:5, :]
    target_lonlat = hcat(fill(110.0, 5), collect(range(30.0, 31.0; length=5)))
    train_lonlat = hcat(
        collect(range(109.0, 112.0; length=n)),
        collect(range(29.0, 33.0; length=n)),
    )
    mixed_roles = Dict("elevation" => "global", "slope" => "local", "aspect" => "local")
    designs = terrain_model_designs(
        terrain, target, train_lonlat, target_lonlat, mixed_roles,
    )
    @test size(designs.mixed_local_train, 2) == 3 + 1 + 2
    @test size(designs.global_train, 2) == 1
    @test designs.local_group_names ==
        ["intercept", "longitude", "latitude", "slope", "aspect"]
    @test size(designs.multiscale_train[end], 2) == 2

    all_local_roles = Dict(group => "local" for group in keys(mixed_roles))
    all_local = terrain_model_designs(
        terrain, target, train_lonlat, target_lonlat, all_local_roles,
    )
    @test size(all_local.mixed_local_train, 2) == 3 + 1 + 1 + 2
    @test size(all_local.global_train, 2) == 0

    changed_target = copy(target)
    changed_target.elevation_m .+= 10_000
    changed_designs = terrain_model_designs(
        terrain, changed_target, train_lonlat, target_lonlat, mixed_roles,
    )
    @test changed_designs.mixed_local_train == designs.mixed_local_train
    @test changed_designs.global_train == designs.global_train
end

@testset "DEM spatial variability and standalone models" begin
    rng = MersenneTwister(123)
    side = 8
    lon = repeat(collect(range(109.7, 111.2; length=side)), side)
    lat = repeat(collect(range(31.4, 33.0; length=side)), inner=side)
    lonlat = hcat(lon, lat)
    n = length(lon)
    elevation = 500 .+ 200 .* randn(rng, n)
    slope = 15 .+ 5 .* randn(rng, n)
    aspect = rand(rng, n) .* 360
    terrain = DataFrame(
        station_id=string.(1:n), elevation_m=elevation, slope_deg=slope,
        aspect_deg=aspect, aspect_sin=sind.(aspect), aspect_cos=cosd.(aspect),
    )

    global_response = 2 .* ((elevation .- mean(elevation)) ./ std(elevation))
    global_roles, _ = spatial_variability_test(
        terrain, lonlat, global_response, ["elevation"];
        bandwidth_candidates=[20, 30], permutations=49, seed=11,
    )
    @test global_roles.role == ["global"]

    elevation_z = (elevation .- mean(elevation)) ./ std(elevation)
    longitude_z = (lon .- mean(lon)) ./ std(lon)
    local_response = elevation_z .* (1 .+ 3 .* longitude_z)
    local_roles, _ = spatial_variability_test(
        terrain, lonlat, local_response, ["elevation"];
        bandwidth_candidates=[16, 24], permutations=99, seed=12,
    )
    @test local_roles.role == ["local"]

    spatial_z = hcat(
        (lon .- mean(lon)) ./ std(lon),
        (lat .- mean(lat)) ./ std(lat),
    )
    Xlocal = hcat(ones(n), spatial_z)
    Xglobal = reshape(elevation_z, :, 1)
    y = 1.5 .* elevation_z .+ 0.3 .* longitude_z
    Y = reshape(y, :, 1)
    bisquare = DEMTerrainExperiment._bisquare_kernel
    mixed_prediction, mixed_converged = mixed_gwr_predict(
        Xlocal, Xglobal, Y, lonlat, Xlocal, Xglobal, lonlat, 24.0, bisquare;
        max_iterations=100,
    )
    @test all(mixed_converged)
    @test all(isfinite, mixed_prediction)
    @test sqrt(mean(abs2, vec(mixed_prediction) - y)) < 0.2

    empty_global = zeros(Float64, n, 0)
    Y_missing = hcat(Y, Y)
    Y_missing[1, 2] = NaN
    loo_prediction, loo_converged = mixed_gwr_predict(
        Xlocal, empty_global, Y_missing, lonlat,
        Xlocal, empty_global, lonlat, 24.0, bisquare; exclude_self=true,
    )
    @test all(loo_converged)
    @test all(isfinite, loo_prediction[:, 1])
    @test isnan(loo_prediction[1, 2])
    @test all(isfinite, loo_prediction[2:end, 2])

    manual_prediction, manual_converged = mixed_gwr_predict(
        Xlocal[2:end, :], zeros(Float64, n - 1, 0), Y[2:end, :], lonlat[2:end, :],
        Xlocal[1:1, :], zeros(Float64, 1, 0), lonlat[1:1, :], 24.0, bisquare,
    )
    @test all(manual_converged)
    @test loo_prediction[1, 1] ≈ manual_prediction[1, 1]

    local_groups = Matrix{Float64}[ones(n, 1), spatial_z[:, 1:1], spatial_z[:, 2:2]]
    multiscale_prediction, multiscale_converged = multiscale_gwr_predict(
        local_groups, Xglobal, Y, lonlat, local_groups, Xglobal, lonlat,
        [24.0, 24.0, 24.0], bisquare; max_iterations=200,
    )
    @test all(multiscale_converged)
    @test all(isfinite, multiscale_prediction)
    @test sqrt(mean(abs2, vec(multiscale_prediction) - y)) < 0.2

    multiscale_missing = hcat(Y, Y)
    multiscale_missing[1, 2] = NaN
    multiscale_loo, multiscale_loo_converged = multiscale_gwr_predict(
        local_groups, empty_global, multiscale_missing, lonlat,
        local_groups, empty_global, lonlat, [24.0, 24.0, 24.0], bisquare;
        exclude_self=true,
    )
    @test all(multiscale_loo_converged)
    @test all(isfinite, multiscale_loo[:, 1])
    @test isnan(multiscale_loo[1, 2])
    @test all(isfinite, multiscale_loo[2:end, 2])

    single_group = Matrix{Float64}[ones(n, 1)]
    changed_Y = copy(Y)
    changed_Y[1] += 1000
    single_loo, _ = multiscale_gwr_predict(
        single_group, empty_global, Y, lonlat,
        single_group, empty_global, lonlat, [24.0], bisquare; exclude_self=true,
    )
    changed_loo, _ = multiscale_gwr_predict(
        single_group, empty_global, changed_Y, lonlat,
        single_group, empty_global, lonlat, [24.0], bisquare; exclude_self=true,
    )
    @test single_loo[1] ≈ changed_loo[1]
    @test_throws ArgumentError multiscale_gwr_predict(
        local_groups, empty_global, Y, lonlat,
        local_groups, empty_global, reverse(lonlat; dims=1), [24.0, 24.0, 24.0], bisquare;
        exclude_self=true,
    )

    scale_rng = MersenneTwister(3)
    scale_side = 10
    scale_lon = repeat(collect(range(0.0, 1.0; length=scale_side)), scale_side)
    scale_lat = repeat(collect(range(0.0, 1.0; length=scale_side)), inner=scale_side)
    scale_lonlat = hcat(109.7 .+ scale_lon, 31.4 .+ scale_lat)
    scale_n = length(scale_lon)
    broad_x = randn(scale_rng, scale_n)
    narrow_x = randn(scale_rng, scale_n)
    broad_beta = 1 .+ 0.5 .* scale_lon
    narrow_beta = 2 .* exp.(-((scale_lon .- 0.5) .^ 2 + (scale_lat .- 0.5) .^ 2) ./ 0.01)
    scale_response = broad_beta .* broad_x + narrow_beta .* narrow_x +
        0.8 .* randn(scale_rng, scale_n)
    scale_groups = Matrix{Float64}[
        ones(scale_n, 1), reshape(broad_x, :, 1), reshape(narrow_x, :, 1),
    ]
    recovered_bandwidths, _, bandwidth_converged = select_multiscale_bandwidths(
        scale_groups, zeros(scale_n, 0), scale_response, scale_lonlat,
        [12.0, 24.0, 48.0, 80.0], bisquare; max_iterations=100,
    )
    @test bandwidth_converged
    @test recovered_bandwidths[2] > recovered_bandwidths[3]

    failed_prediction, failed_convergence = multiscale_gwr_predict(
        local_groups, Xglobal, Y, lonlat, local_groups, Xglobal, lonlat,
        [24.0, 24.0, 24.0], bisquare; max_iterations=0,
    )
    @test !all(failed_convergence)
    @test all(isnan, failed_prediction)
end

"""
Weighted-least-squares prediction at every target, written straight from the definition.

`_local_hat` returns the linear operator that produces this, so `hat * y` must match it
for any `y`. Deriving the reference from the definition rather than from a copy of the
implementation is what makes these tests bite on a rewrite of the inner loop.
"""
function reference_local_prediction(
    Xtrain, Xtarget, distances, y, neighbors; ridge=1e-8, exclude_self=false,
)
    p = size(Xtrain, 2)
    prediction = zeros(Float64, size(Xtarget, 1))
    for target in axes(Xtarget, 1)
        d = copy(distances[:, target])
        exclude_self && (d[target] = Inf)
        w = DEMTerrainExperiment._adaptive_bisquare(d, neighbors)
        valid = w .> 0
        count(valid) >= p + 1 || continue
        Xv = Xtrain[valid, :]
        wv = w[valid]
        beta = (Xv' * (wv .* Xv) + ridge * I) \ (Xv' * (wv .* y[valid]))
        prediction[target] = dot(Xtarget[target, :], beta)
    end
    return prediction
end

@testset "adaptive bisquare weights" begin
    distances = [5.0, 1.0, 9.0, 3.0, Inf, 7.0]
    bandwidth = sort(filter(isfinite, distances))[3]
    @test bandwidth == 5.0
    weights = DEMTerrainExperiment._adaptive_bisquare(distances, 3)
    for (index, d) in enumerate(distances)
        @test weights[index] ==
            (isfinite(d) && d < bandwidth ? (1 - (d / bandwidth)^2)^2 : 0.0)
    end
    @test all(iszero, weights[distances .>= bandwidth])
    @test count(>(0), weights) == 2

    # Asking for more neighbours than there are finite distances clamps to the largest.
    @test DEMTerrainExperiment._adaptive_bisquare(distances, 99) ==
        DEMTerrainExperiment._adaptive_bisquare(distances, 5)
    # No finite distance at all: no station gets weight.
    @test all(iszero, DEMTerrainExperiment._adaptive_bisquare(fill(Inf, 4), 2))
    # Every finite distance zero: the kernel degenerates to exact coincidence.
    @test DEMTerrainExperiment._adaptive_bisquare([0.0, 0.0, Inf], 2) == [1.0, 1.0, 0.0]
    @test DEMTerrainExperiment._adaptive_bisquare(zeros(4), 2) == ones(4)
end

@testset "local hat operator matches weighted least squares" begin
    bisquare = DEMTerrainExperiment._bisquare_kernel
    rng = MersenneTwister(2026)
    n, m = 40, 11
    lonlat = hcat(109.0 .+ 4 .* rand(rng, n), 30.0 .+ 4 .* rand(rng, n))
    target_lonlat = hcat(109.0 .+ 4 .* rand(rng, m), 30.0 .+ 4 .* rand(rng, m))
    train_distances = DEMTerrainExperiment._haversine_matrix(lonlat, lonlat)
    target_distances = DEMTerrainExperiment._haversine_matrix(lonlat, target_lonlat)
    y = randn(rng, n)
    for p in 1:5
        Xtrain = hcat(ones(n), randn(rng, n, p - 1))
        Xtarget = hcat(ones(m), randn(rng, m, p - 1))
        for bandwidth in (10, 20, n)
            square = DEMTerrainExperiment._local_hat(
                Xtrain, Xtrain, train_distances, Float64(bandwidth), bisquare,
            )
            @test square * y ≈ reference_local_prediction(
                Xtrain, Xtrain, train_distances, y, bandwidth,
            )

            held_out = DEMTerrainExperiment._local_hat(
                Xtrain, Xtrain, train_distances, Float64(bandwidth), bisquare; exclude_self=true,
            )
            @test all(iszero, diag(held_out))
            @test held_out * y ≈ reference_local_prediction(
                Xtrain, Xtrain, train_distances, y, bandwidth; exclude_self=true,
            )

            rectangular = DEMTerrainExperiment._local_hat(
                Xtrain, Xtarget, target_distances, Float64(bandwidth), bisquare,
            )
            @test size(rectangular) == (m, n)
            @test rectangular * y ≈ reference_local_prediction(
                Xtrain, Xtarget, target_distances, y, bandwidth,
            )
        end
    end

    # A bandwidth too small to support the design leaves the row at zero rather than
    # throwing — `_gwr_smoothers` is the one that throws, `_local_hat` is not.
    Xwide = hcat(ones(n), randn(rng, n, 5))
    @test all(iszero, DEMTerrainExperiment._local_hat(
        Xwide, Xwide, train_distances, 2.0, bisquare,
    ))

    # The ridge is applied to the normal equations, so a larger one shrinks the operator.
    Xtrain = hcat(ones(n), randn(rng, n, 2))
    weak = DEMTerrainExperiment._local_hat(
        Xtrain, Xtrain, train_distances, 12.0, bisquare; ridge=1e-8,
    )
    strong = DEMTerrainExperiment._local_hat(
        Xtrain, Xtrain, train_distances, 12.0, bisquare; ridge=1.0,
    )
    @test !(weak ≈ strong)
    @test strong * y ≈ reference_local_prediction(
        Xtrain, Xtrain, train_distances, y, 12; ridge=1.0,
    )

    @test_throws DimensionMismatch DEMTerrainExperiment._local_hat(
        Xtrain, randn(rng, m, 2), target_distances, 12.0, bisquare,
    )
    @test_throws DimensionMismatch DEMTerrainExperiment._local_hat(
        Xtrain, Xtrain, target_distances, 12.0, bisquare,
    )
    @test_throws DimensionMismatch DEMTerrainExperiment._local_hat(
        Xtrain, hcat(ones(m), randn(rng, m, 2)), target_distances, 12.0, bisquare;
        exclude_self=true,
    )
end
