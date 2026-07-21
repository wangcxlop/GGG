# 与 ex03_spatiotemporal_ShiYan.R 对齐
# JULIA_NUM_THREADS=32 N_TIME=1000 julia --project=. example/ex03_spatiotemporal_ShiYan_bench.jl
using MixedGWR, SpatialRasterLite, ArchGDAL, Distances
using Ipaper, RTableTools, NaNStatistics

n_time_req = parse(Int, get(ENV, "N_TIME", "100"))
bw = parse(Float64, get(ENV, "BW", "20"))
n_max = parse(Int, get(ENV, "N_MAX", "8"))
fun_dist = Haversine(6378.388)

# 网格
dem = rast("data/dem_ShiYan_2km.tif", FT=Float64)
X1 = st_coords(dem)
Points = map(x -> x, eachrow(X1))
Xpred = cbind(X1, dem.A[:] * 1.0)

# 站点 / 降水
st = fread("/mnt/z/GitHub/jl-pkgs/SpatInterp.jl/Project_十堰/data/十堰_雨量站_sp237_v20250824.csv")
pr = fread("/mnt/z/GitHub/jl-pkgs/SpatInterp.jl/Project_十堰/data/十堰_prcp_sp237_mat_v20250824.csv")
coords = Matrix{Float64}(st[:, [:lon, :lat]])
points = map(x -> (x[1], x[2]), eachrow(coords))
P = coalesce.(Matrix{Union{Missing,Float64}}(pr[:, Symbol.(string.(st.site))]), NaN)
inds = findall(NaNStatistics.nansum(P, dims=2)[:] .>= 1)
n_time_req > 0 && (inds = inds[1:min(n_time_req, length(inds))])
Y = Matrix{Float64}(P[inds, :])'
Y[.!isfinite.(Y)] .= 0
X = cbind(coords, st_extract(dem, points).value' |> Matrix)

# 计时
dMat = pairwise(fun_dist, points, Points)
t_w = @elapsed wMat = gw_weight(dMat, bw; kernel=BISQUARE, adaptive=false)
ST_GWR_fast(X, Y[:, 1:2], wMat; Xpred, n_max)  # warmup

t_gwr = @elapsed Ypred = ST_GWR_fast(X, Y, wMat; Xpred, n_max)

n_st, n_time = size(Y)
@info "n_st=$n_st n_target=$(size(Xpred,1)) n_time=$n_time bw=$bw n_max=$n_max threads=$(Threads.nthreads())"
@info "gw_weight  $(round(t_w;digits=3)) s"
@info "ST_GWR_fast $(round(t_gwr;digits=3)) s  ($(round(t_gwr/n_time;digits=4)) s/day)"
@info "Ypred range=$(extrema(Ypred))"
