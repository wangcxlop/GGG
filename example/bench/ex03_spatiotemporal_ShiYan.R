#!/usr/bin/env Rscript
# R + GWmodel 对照 ex03_spatiotemporal_ShiYan.jl
# N_TIME=1000 BW=20 N_MAX=8 Rscript example/ex03_spatiotemporal_ShiYan.R
pacman::p_load(GWmodel, terra, data.table)

n_time_req <- as.integer(Sys.getenv("N_TIME", "100"))
bw <- as.numeric(Sys.getenv("BW", "20"))
n_max <- as.integer(Sys.getenv("N_MAX", "8"))

## 网格
dem <- rast("data/dem_ShiYan_2km.tif")
xy <- crds(dem, na.rm = FALSE)
alt <- values(dem)[, 1]
alt[!is.finite(alt)] <- mean(alt, na.rm = TRUE)
Xpred <- cbind(xy, alt)

## 站点 / 降水
st <- fread("/mnt/z/GitHub/jl-pkgs/SpatInterp.jl/Project_十堰/data/十堰_雨量站_sp237_v20250824.csv")
st[, alt := terra::extract(dem, cbind(lon, lat))$elevation]
pr <- fread("/mnt/z/GitHub/jl-pkgs/SpatInterp.jl/Project_十堰/data/十堰_prcp_sp237_mat_v20250824.csv")
P <- as.matrix(pr[, as.character(st$site), with = FALSE])
inds <- which(rowSums(P, na.rm = TRUE) >= 1)
if (n_time_req > 0) inds <- head(inds, n_time_req)
Y <- t(P[inds, ])                          # [n_st, n_time]
Y[!is.finite(Y)] <- 0
X <- as.matrix(st[, .(lon, lat, alt)])

## 距离 / 权重（longlat 时 gw.dist 行列颠倒，需转置）
dMat <- {
  d <- gw.dist(dp.locat = X[, 1:2], rp.locat = xy, longlat = TRUE)
  if (nrow(d) == nrow(X)) d else t(d)
}
t_w <- system.time({
  wMat <- apply(dMat, 2, gw.weight, bw = bw, kernel = "bisquare", adaptive = FALSE)
})[3]

## ST-GWR_fast：每格点取 n_max 邻域，多时刻一次求解
st_gwr_fast <- function(X, Y, wMat, Xpred, n_max = 8L, λ = 1e-8) {
  p <- ncol(X)
  n_max <- min(n_max, nrow(X))
  Iλ <- diag(λ, p)
  Ypred <- matrix(NA_real_, ncol(wMat), ncol(Y))
  for (j in seq_len(ncol(wMat))) {
    idx <- order(wMat[, j], decreasing = TRUE)[1:n_max]
    ww <- wMat[idx, j]
    if (sum(ww) <= 0) next
    Xw <- X[idx, , drop = FALSE] * ww
    beta <- tryCatch(solve(crossprod(X[idx, , drop = FALSE], Xw) + Iλ,
                           crossprod(Xw, Y[idx, , drop = FALSE])),
                     error = function(e) NULL)
    if (!is.null(beta)) Ypred[j, ] <- drop(Xpred[j, ] %*% beta)
  }
  Ypred
}

t_gwr <- system.time(Ypred <- st_gwr_fast(X, Y, wMat, Xpred, n_max))[3]

## 原生 gwr.predict 单日（对照）
keep <- is.finite(Y[, 1])
dp <- SpatialPointsDataFrame(X[keep, 1:2],
  data.frame(y = Y[keep, 1], lon = X[keep, 1], lat = X[keep, 2], alt = X[keep, 3]))
rp <- SpatialPointsDataFrame(xy, data.frame(lon = Xpred[, 1], lat = Xpred[, 2], alt = Xpred[, 3]))
t_pred <- system.time(
  gwr.predict(y ~ lon + lat + alt - 1, dp, rp,
    bw = max(n_max, 12), kernel = "bisquare", adaptive = TRUE,
    dMat1 = dMat[keep, ], dMat2 = gw.dist(dp.locat = coordinates(dp), longlat = TRUE))
)[3]

## 汇总
nt <- ncol(Y)
cat(sprintf(
  "n_st=%d  n_target=%d  n_time=%d  bw=%.0f  n_max=%d\n",
  nrow(X), nrow(Xpred), nt, bw, n_max))
cat(sprintf("  gw.weight          %7.3f s\n", t_w))
cat(sprintf("  ST_GWR_fast        %7.3f s  (%.4f s/day)\n", t_gwr, t_gwr / nt))
cat(sprintf("  gwr.predict 1day   %7.3f s\n", t_pred))
cat(sprintf("  外推 gwr.predict×n %7.0f s  → 加速 %.0fx\n", t_pred * nt, t_pred * nt / t_gwr))
cat(sprintf("  Ypred range [%.2f, %.2f]\n", min(Ypred, na.rm = TRUE), max(Ypred, na.rm = TRUE)))
