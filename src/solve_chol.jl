# https://chatgpt.com/c/68aadc51-a978-8326-8880-9df82cecb413
using LinearAlgebra

export GWRSolver
export solve_chol!, solve_chol


# 这里存储的是临时变量
@with_kw mutable struct GWRSolver{T}
  n_control::Int = 100
  n_max::Int = 8
  p::Int = 2
  n_time::Int = 10
  # X::Matrix{T} = zeros(T, n_control, p)           # k×p
  # Y::Matrix{T} = zeros(T, n_control, n_time)
  # w::Vector{T} = zeros(T, n_control)              # k, 
  WY::Matrix{T} = zeros(T, n_max, n_time)
  XtWX::Matrix{T} = zeros(T, p, p)
  XtWY::Matrix{T} = zeros(T, p, n_time)
  β::Matrix{T} = zeros(T, p, n_time)
end

function GWRSolver(X::AbstractMatrix{T}, Y::AbstractMatrix{T};
  n_max::Int=size(X, 1)) where {T<:Real}
  n_control, p = size(X)
  n_time = size(Y, 2)
  GWRSolver{T}(; n_control, n_max, p, n_time)
end


"""
# Arguments
- `X`: [k, p]
- `w`: [k]
- `Y`: [k, ntime]
"""
function _chol_solve!(β::AbstractMatrix{T}, XtWX::AbstractMatrix{T},
  XtWY::AbstractMatrix{T}; λ::T=T(1e-8)) where {T<:Real}
  @inbounds @simd for i in axes(XtWX, 1)
    XtWX[i, i] += λ
  end
  F = cholesky!(Symmetric(XtWX, :L); check=false)
  ldiv!(β, F, XtWY)
  return β
end

function solve_chol!(
  β::AbstractMatrix{T}, WY::AbstractMatrix{T}, XtWY::AbstractMatrix{T},
  X::AbstractMatrix{T}, Y::AbstractMatrix{T}, w::AbstractVector{T};
  λ::T=T(1e-8)) where {T<:Real}
  n_control = size(X, 1)
  ntime = size(Y, 2)

  Xw = X .* w
  XtWX = transpose(X) * Xw
  @turbo for i in 1:n_control
    wi = w[i]
    for j in 1:ntime
      WY[i, j] = Y[i, j] * wi
    end
  end
  mul!(XtWY, transpose(X), WY)
  _chol_solve!(β, XtWX, XtWY; λ)
end


function solve_chol!(solver::GWRSolver{T},
  X::AbstractMatrix{T}, Y::AbstractMatrix{T}, w::AbstractVector{T};
  λ::T=T(1e-8)) where {T<:Real}
  (; β, WY, XtWY) = solver
  solve_chol!(β, WY, XtWY, X, Y, w; λ)
end

function _solve_chol_fast!(solver::GWRSolver{T},
  X::AbstractMatrix{T}, Y::AbstractMatrix{T}, inds::AbstractVector{Int},
  w::AbstractVector{T}; λ::T=T(1e-8)) where {T<:Real}
  (; β, XtWX, XtWY) = solver
  p = size(X, 2)
  ntime = size(Y, 2)
  fill!(XtWX, zero(T))
  fill!(XtWY, zero(T))

  @inbounds for r in eachindex(inds, w)
    i = inds[r]
    wi = w[r]
    for b in 1:p
      xwb = wi * X[i, b]
      for a in b:p
        XtWX[a, b] += X[i, a] * xwb
      end
    end

    for t in 1:ntime
      ywi = wi * Y[i, t]
      for a in 1:p
        XtWY[a, t] += X[i, a] * ywi
      end
    end
  end
  _chol_solve!(β, XtWX, XtWY; λ)
end


function solve_chol(
  X::AbstractMatrix{T}, Y::AbstractMatrix{T}, w::AbstractVector{T}; λ::T=T(1e-8)) where {T<:Real}
  solver = GWRSolver(X, Y)
  solve_chol!(solver, X, Y, w; λ)
end

solve_chol(X::AbstractMatrix{T}, y::AbstractVector{T}, w::AbstractVector{T}) where {T<:Real} =
  solve_chol(X, Matrix(y), w)
