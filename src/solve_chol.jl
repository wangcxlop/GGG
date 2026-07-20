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
  WY::Matrix{T} = zeros(T, n_control, n_time)
  XW::Matrix{T} = zeros(T, p, n_max)
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

function _fast_solver(X::AbstractMatrix{T}, Y::AbstractMatrix{T},
  n_max::Int) where {T<:Real}
  n_control, p = size(X)
  n_time = size(Y, 2)
  R = zeros(T, p, n_time)
  GWRSolver{T}(;
    n_control, n_max, p, n_time,
    WY=zeros(T, 0, 0), XW=zeros(T, p, n_max), XtWY=R, β=R)
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
  β === XtWY || copyto!(β, XtWY)
  ldiv!(F, β)
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
  (; β, XW, XtWX, XtWY) = solver
  p = size(X, 2)
  ntime = size(Y, 2)
  n = length(inds)

  @inbounds for r in 1:n
    i = inds[r]
    for a in 1:p
      XW[a, r] = X[i, a] * w[r]
    end
  end

  # XtWX是对称矩阵，只计算下三角（a ≥ b), Symmetric(XtWX, :L)补齐
  @inbounds for b in 1:p
    for a in b:p
      s = zero(T)
      for r in 1:n
        s += X[inds[r], a] * XW[b, r]
      end
      XtWX[a, b] = s
    end
  end

  @turbo for t in 1:ntime
    for a in 1:p
      s = zero(T)
      for r in 1:n
        s += XW[a, r] * Y[inds[r], t]
      end
      XtWY[a, t] = s
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
