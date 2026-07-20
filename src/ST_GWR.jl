"""
- `X`    : [n_control, k_local], k_local = 3, lon + lat + alt
- `Y`    : [n_control, n_time]
- `wMat` : [n_control, n_target]

- `Xpred`: [n_target, k_local]
"""
function ST_GWR(X::AbstractMatrix{T}, Y::AbstractMatrix{T},
  wMat::AbstractMatrix{T}; Xpred::AbstractMatrix{T}) where {T<:Real}
  n_target = size(wMat, 2)
  ntime = size(Y, 2)
  # k_local = size(X, 2)
  # β = zeros(T, n_target, k_local, ntime)
  n_target = size(wMat, 2)
  solvers = map(i -> GWRSolver(X, Y), 1:get_nthread())

  Ypred = zeros(T, n_target, ntime)
  p = Progress(n_target)
  @inbounds @threads for i in 1:n_target
    next!(p)
    k = Threads.threadid()
    solver = solvers[k]

    w = @view wMat[:, i]
    _β = solve_chol!(solver, X, Y, w) # [np, ntime]
    # β[i, :, :] .= _β
    _X = Xpred[i, :] # [np]
    _pred = @view Ypred[i, :] # [ntime]
    fitted!(_pred, _X, _β)
  end
  return Ypred
end

function ST_GWR_fast(X::AbstractMatrix{T}, Y::AbstractMatrix{T},
  wMat::AbstractMatrix{T}; Xpred::AbstractMatrix{T}) where {T<:Real}


end
