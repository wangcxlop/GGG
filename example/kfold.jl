module KfoldGWR

using DataFrames: DataFrame
using MixedGWR
using ModelParams: GOF
using Random: MersenneTwister, randperm
using Statistics: mean, std

export kfold_gof, st_gwr_kfold

function kfold_gof(Y, params, predict_fold; k=5, seed=42, metadata=(;))
  n = size(Y, 1)
  @assert 2 <= k <= n

  order = randperm(MersenneTwister(seed), n)
  folds = [order[i:k:end] for i in 1:k]
  valid_time = vec(all(isfinite, Y; dims=1))
  any(valid_time) || error("no complete observations for cross-validation")
  Y = Y[:, valid_time]
  rows = NamedTuple[]

  for p in params
    pred = fill(NaN, size(Y))
    for test in folds
      train = setdiff(1:n, test)
      pred[test, :] = predict_fold(train, test, Y[train, :], p)
    end

    gof = GOF(vec(Y), vec(pred))
    push!(rows, merge(
      metadata, p, gof,
      (; n_time=size(Y, 2), coverage=gof.n_valid / length(Y)),
    ))
  end

  sort!(DataFrame(rows), [:coverage, :RMSE]; rev=[true, false])
end

function zscore_train_test(X_train, X_test)
  μ = mean(X_train; dims=1)
  σ = std(X_train; dims=1)
  (X_train .- μ) ./ σ, (X_test .- μ) ./ σ
end

function predict_fold(::Val{:gwr}, X, dMat, train, test, Y_train, p;
  kernel, standardize)
  X_train = X[train, 1:p.np]
  X_test = X[test, 1:p.np]
  if standardize
    X_train, X_test = zscore_train_test(X_train, X_test)
    X_train = hcat(ones(length(train)), X_train)
    X_test = hcat(ones(length(test)), X_test)
  end

  wMat = gw_weight(dMat[train, test], p.bw;
    kernel, adaptive=p.adaptive)
  ST_GWR_fast(
    X_train, Y_train, wMat; Xpred=X_test,
    n_max=min(p.n_max, length(train)),
  )
end

function st_gwr_mixed_fast(X_local, X_global, Y, wMat, wMat_rp;
  Xlocal_pred, Xglobal_pred, n_max)
  x3 = X_global - ST_GWR_fast(
    X_local, X_global, wMat; Xpred=X_local, n_max,
  )
  y2 = Y - ST_GWR_fast(X_local, Y, wMat; Xpred=X_local, n_max)
  β_global = x3 \ y2
  y_local = Y - X_global * β_global

  ST_GWR_fast(
    X_local, y_local, wMat_rp; Xpred=Xlocal_pred, n_max,
  ) + Xglobal_pred * β_global
end

function predict_fold(::Val{:mixed}, X, dMat, train, test, Y_train, p;
  kernel, standardize)
  standardize || error("Mixed GWR requires standardized predictors")
  Z_train, Z_test = zscore_train_test(X[train, :], X[test, :])
  X_local = hcat(ones(length(train)), Z_train[:, 1:2])
  Xlocal_pred = hcat(ones(length(test)), Z_test[:, 1:2])
  X_global = Z_train[:, 3:3]
  Xglobal_pred = Z_test[:, 3:3]

  wMat = gw_weight(dMat[train, train], p.bw;
    kernel, adaptive=p.adaptive)
  wMat_rp = gw_weight(dMat[train, test], p.bw;
    kernel, adaptive=p.adaptive)
  st_gwr_mixed_fast(
    X_local, X_global, Y_train, wMat, wMat_rp;
    Xlocal_pred, Xglobal_pred,
    n_max=min(p.n_max, length(train)),
  )
end

model_metadata(::Val{:gwr}, standardize) =
  (; model=standardize ? "gwr_zscore_intercept" : "gwr_raw", standardize)

model_metadata(::Val{:mixed}, standardize) =
  (; model="mixed", standardize, local_vars="lon+lat", global_var="alt")

function st_gwr_kfold(model::Symbol, X, Y, dMat, params;
  k=5, seed=42, kernel=BISQUARE, standardize=true)
  method = Val(model)
  predict = (train, test, Y_train, p) -> predict_fold(
    method, X, dMat, train, test, Y_train, p; kernel, standardize,
  )
  kfold_gof(Y, params, predict;
    k, seed, metadata=model_metadata(method, standardize))
end

end
