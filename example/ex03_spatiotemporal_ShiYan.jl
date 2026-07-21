# ~3 minutes
using MixedGWR, SpatialRasterLite, ArchGDAL, Distances
using Ipaper, RTableTools, NetCDFTools
using JLD2, UnPack, NaNStatistics

include(joinpath(@__DIR__, "kfold.jl"))
using .KfoldGWR

fun_dist = Haversine(6378.388)

# 目标网格
begin
  dem = rast("data/dem_ShiYan_2km.tif", FT=Float64)
  X1 = st_coords(dem)
  Points = map(x -> x, eachrow(X1))

  Xpred = cbind(X1, dem.A[:]) # [lon, lat, alt]
end

f = "/mnt/z/GitHub/jl-pkgs/SpatInterp.jl/Project_十堰/data/ShiYan_Pobs_interpolated_by_IDW.jld2"
l = jldopen(f)
@unpack st, dates, P = l

## 
# INPUT
begin
  coords = Matrix(st[:, [:lon, :lat]])
  points = map(x -> (x[1], x[2]), eachrow(coords))
  dMat = pairwise(fun_dist, points)
  dMat_rp = pairwise(fun_dist, points, Points)

  x1 = st[:, [:lon, :lat]] |> Matrix
  x2 = st_extract(dem, points).value' |> Matrix
  X = cbind(x1, x2)

  ## 只对有降水的日期进行插值
  _prcp = NaNStatistics.nansum(P, dims=1)[:]
  inds = findall(_prcp .>= 1.0) # 所有站点总降水大于1mm日期，才进行插值, ~1/3
  Y = P[:, inds]
end

# 站点空间5折交叉验证：每个站点恰好作为验证样本一次
spatial_params = [
  (; adaptive, bw, n_max)
  for (adaptive, bandwidths) in (
    (true, (6.0, 10.0, 20.0)),
    (false, (10.0, 20.0, 30.0, 50.0)),
  )
  for bw in bandwidths
  for n_max in (8, 12)
]
param_grid = [merge((; np), p) for np in (2, 3) for p in spatial_params]

# V1：原始版本，不标准化、无截距
scores_raw = st_gwr_kfold(
  :gwr, X, Y, dMat, param_grid; k=5, seed=42, standardize=false,
)
fwrite(scores_raw, "OUTPUT/ShiYan_ST_GWR_kfold.csv")

# V2：所有解释变量标准化，并加入截距
scores = st_gwr_kfold(
  :gwr, X, Y, dMat, param_grid; k=5, seed=42, standardize=true,
)
fwrite(scores, "OUTPUT/ShiYan_ST_GWR_kfold_zscore.csv")

# V3：加入截距，仅保留标准化后的经、纬度
scores_noalt = scores[scores.np .== 2, :]
fwrite(scores_noalt, "OUTPUT/ShiYan_ST_GWR_kfold_intercept_noalt.csv")
show(scores_noalt; allrows=true, allcols=true)

# V4：经、纬度为局部变量，高程为全局变量
mixed_scores = st_gwr_kfold(
  :mixed, X, Y, dMat, spatial_params; k=5, seed=42,
)
fwrite(mixed_scores, "OUTPUT/ShiYan_ST_GWR_mixed_kfold.csv")
show(mixed_scores; allrows=true, allcols=true)

## 5 second, (69696, 1000) [n_target, n_time]
adaptive = false
bandwidths = (10.0, 20.0, 30.0, 50.0) # km

outdir = "OUTPUT" |> path_mnt

for bw in bandwidths
  fout = "$outdir/ShiYan_Prcp_Gauged237_201404-202501_2km_GWR3(adaptive=$adaptive,bw=$bw).nc"
  isfile(fout) && continue

  kernel = BISQUARE
  wMat_rp = gw_weight(dMat_rp, bw; kernel, adaptive)

  np = 3 # lon + lat + alt
  @time Ypred = ST_GWR_fast(X[:, 1:np], Y, wMat_rp; Xpred=Xpred[:, 1:np], n_max=8)

  nlon, nlat = size(dem)[1:2]
  R = zeros(Float32, nlon, nlat, length(dates))
  R[:, :, inds] .= reshape(Ypred, nlon, nlat, length(inds))

  lon, lat = st_dims(dem)
  dims = (; lon, lat, time=dates)
  ncsave(fout, true, (; units="mm h-1"); dims, P=R)
end
