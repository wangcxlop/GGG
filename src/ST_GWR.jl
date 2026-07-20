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

function gwr_neighbors(wMat::AbstractMatrix{T}, n_max::Int) where {T<:Real}
  1 <= n_max <= size(wMat, 1) || throw(ArgumentError("invalid n_max"))
  _, n_target = size(wMat)
  inds = Matrix{Int}(undef, n_max, n_target)
  ws = Matrix{T}(undef, n_max, n_target)

  @inbounds for i in 1:n_target
    w = @view wMat[:, i]
    idx = partialsortperm(w, 1:n_max; rev=true)
    for r in 1:n_max
      inds[r, i] = idx[r]
      ws[r, i] = w[idx[r]]
    end
  end
  return inds, ws
end

function ST_GWR_fast!(Ypred::AbstractMatrix{T}, X::AbstractMatrix{T},
  Y::AbstractMatrix{T}, inds::Matrix{Int}, ws::Matrix{T};
  Xpred::AbstractMatrix{T}) where {T<:Real}
  size(inds) == size(ws) || throw(DimensionMismatch("inds and ws differ"))

  n_max, n_target = size(inds)
  solvers = map(_ -> _fast_solver(X, Y, n_max), 1:Threads.nthreads(:default))
  thread_offset = Threads.nthreads(:interactive)

  @inbounds @threads :static for i in 1:n_target
    solver = solvers[Threads.threadid() - thread_offset]
    β = _solve_chol_fast!(solver, X, Y, @view(inds[:, i]), @view(ws[:, i]))
    fitted!(@view(Ypred[i, :]), @view(Xpred[i, :]), β)
  end
  return Ypred
end

"""Fast ST-GWR using the `n_max` largest spatial weights per target."""
function ST_GWR_fast(X::AbstractMatrix{T}, Y::AbstractMatrix{T},
  inds::Matrix{Int}, ws::Matrix{T};
  Xpred::AbstractMatrix{T}) where {T<:Real}

  Ypred = zeros(T, size(inds, 2), size(Y, 2))
  ST_GWR_fast!(Ypred, X, Y, inds, ws; Xpred)
end

function ST_GWR_fast(X::AbstractMatrix{T}, Y::AbstractMatrix{T},
  wMat::AbstractMatrix{T}; Xpred::AbstractMatrix{T}, n_max::Int=8) where {T<:Real}

  inds, ws = gwr_neighbors(wMat, n_max)
  ST_GWR_fast(X, Y, inds, ws; Xpred)
end
