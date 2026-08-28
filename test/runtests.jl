using MixedGWR, RTableTools, Distances, Test
include("main_pkgs.jl");

# These four validate against R's GWmodel on the Shiyan station fixture. `data/` is gitignored,
# so the fixture is absent from a fresh checkout and from CI; they are skipped rather than
# aborting the suite. Everything below them is self-contained and always runs.
if HAS_SHIYAN_DATA
  @testset "GWR" begin
    for adaptive in [true, false]
      for kernel in 0:4
        r = R"GWmodel:::gwr_q($x1, $y, $dMat, 20.0, $kernel + 1, $adaptive)" |> rcopy
        jl = GWR(x1, y, dMat, 20.0; kernel, adaptive)
        @test r ≈ jl
      end
    end
  end

  include("test-solver.jl")
  include("test-ST_GWR.jl")
  include("test-GWR_mixed.jl")
else
  @warn "Skipping the R/GWmodel comparison testsets: $SHIYAN_DATA_PATH is missing"
end
include("test-study-area.jl")
include("test-mger-five-kernels.jl")
include("test-interpolation-benchmark.jl")
include("test-appeears-ndvi.jl")
include("test-era5-land-stations.jl")
include("test-era5-land-processing.jl")
include("test-terrain-features.jl")
include("test-era5-land-covariates.jl")
include("test-mod13a2-ndvi-processing.jl")
include("test-dem-terrain-experiment.jl")
include("test-era5-variable-selection.jl")
include("test-ndvi-variable-selection.jl")
include("test-joint-variable-selection.jl")
include("test-joint-covariate-models.jl")
include("test-benchmark-diagnostics.jl")
