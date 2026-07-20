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
  @time Ypred2 = ST_GWR(Xlocal, Y, wMat; Xpred=Xlocal);

  Ypred_fast = similar(Ypred2)
  inds, ws = gwr_neighbors(wMat, Int(bw))
  @time ST_GWR_fast!(Ypred_fast, Xlocal, Y, inds, ws; Xpred=Xlocal);

  @test Ypred2[:, 1] == Ypred2[:, 2]
  @test Ypred_fast ≈ Ypred2
end

# 3.080392 seconds (1.96 k allocations: 3.976 GiB, 10.74% gc time, 27 lock conflicts)
# 0.462219 seconds (1.73 k allocations: 664.364 MiB, 13.87% gc time)

## 
# using BenchmarkTools
# @btime Ypred2 = ST_GWR(Xlocal, Y, wMat; Xpred=Xlocal);
# @btime ST_GWR_fast!($Ypred_fast, $Xlocal, $Y, $inds, $ws;
#   Xpred=$Xlocal);
