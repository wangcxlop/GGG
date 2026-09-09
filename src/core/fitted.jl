"""
时间上的处理
- x : [np]
- β : [np, ntime]
- R : [ntime]
"""
function fitted!(R::AbstractVector{T},
  x::AbstractVector{T}, β::AbstractMatrix{T}) where {T<:Real}
  p, ntime = size(β)
  @turbo for i in 1:ntime
    ∑ = zero(T)
    for j in 1:p
      ∑ += (x[j] * β[j, i])
    end
    R[i] = ∑
  end
  return R
end

"""
空间上的处理
- x : [n, np]
- β : [n, np]
- R : [n]
"""
function fitted!(R::AbstractVector{T},
  x::AbstractMatrix{T}, β::AbstractMatrix{T}) where {T<:Real}
  n, p = size(x)
  @turbo for i in 1:n
    ∑ = zero(T)
    for j in 1:p
      ∑ += x[i, j] * β[i, j]
    end
    R[i] = ∑
  end
  return R
end

fitted(x::AbstractMatrix{T}, β::AbstractMatrix{T}) where {T<:Real} =
  fitted!(zeros(T, size(x, 1)), x, β)
