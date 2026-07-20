@testset "ST_GWR" begin
  coords = Matrix(d[:, [:lon, :lat]])
  points = map(x -> x, eachrow(coords))

  fun_dist = Haversine(6378.388)
  dMat = pairwise(fun_dist, points)

  x1 = d[:, [:lon, :lat]] |> Matrix
  x2 = d[:, [:alt]] .* 1.0 |> Matrix
  y = d[:, :prcp]
  X = d[:, [:lon, :lat, :alt]] |> Matrix
  Y = repeat(y, outer = (1, 10_0000))

  ## Parameters
  kernel = BISQUARE
  adaptive = true
  bw = 20.0 # 6个站点
  wMat = gw_weight(dMat, bw; kernel, adaptive)

  ## Run
  np = 2 # lon + lat
  Xlocal = X[:, 1:np]
  @time Ypred2 = ST_GWR(Xlocal, Y, wMat; Xpred=Xlocal)
  @time Ypred_fast = ST_GWR_fast(Xlocal, Y, wMat; Xpred=Xlocal, n_max=Int(bw))
  @test Ypred2[:, 1] == Ypred2[:, 2]
  @test Ypred_fast ≈ Ypred2
end

## 
# using BenchmarkTools
# @btime Ypred2 = ST_GWR(Xlocal, Y, wMat; Xpred=Xlocal);
# @btime Ypred_fast = ST_GWR_fast(Xlocal, Y, wMat; 
#   Xpred=Xlocal, n_max=Int(bw));
