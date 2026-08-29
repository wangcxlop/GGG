using MixedGWR, RTableTools, Distances, Test
using RCall
R"""
library(GWmodel)
"""


function GWR_mixed_r(x1, x2, y, dMat; kernel=BISQUARE)
  R"GWmodel:::gwr_mixed_2($x1, $x2, $y, $dMat, $dMat, 20.0, $kernel+1, TRUE)" |> rcopy
end

function GWR_mixed_trace_r(x1, x2, y, dMat; kernel=BISQUARE)
  R"GWmodel:::gwr_mixed_trace($x1, $x2, $y, $dMat, 20.0, $kernel+1, TRUE)" |> rcopy
end


# load data
#
# `data/` is gitignored, so this fixture is not present in a fresh checkout (nor in CI). The
# R-comparison testsets that need it are skipped rather than allowed to abort the whole suite -
# see the `HAS_SHIYAN_DATA` guard in runtests.jl. Restoring the file re-enables them with no
# other change.
const SHIYAN_DATA_PATH = abspath(joinpath(@__DIR__, "..", "data", "prcp_st174_shiyan.csv"))
const HAS_SHIYAN_DATA = isfile(SHIYAN_DATA_PATH)
if HAS_SHIYAN_DATA
  d = fread(SHIYAN_DATA_PATH)

  coords = Matrix(d[:, [:lon, :lat]])
  points = map(x -> x, eachrow(coords))

  fun_dist = Haversine(6378.388)
  dMat = pairwise(fun_dist, points)

  x1 = d[:, [:lon, :lat]] |> Matrix
  x2 = d[:, [:alt]] .* 1.0 |> Matrix
  y = d[:, :prcp];
end
