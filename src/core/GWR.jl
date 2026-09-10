"""
GWR with specified distance matrix

# Arguments
- `x`: control variables, [n_control, k_local]
- `y`: response variable
- `dMat`: distance matrix, [n_control, n_target]
- `β`: [n_target, k_local]
提前算好权重，进行加速
"""
function GWR!(β::AbstractMatrix{T}, 
  X::AbstractMatrix{T}, Y::AbstractVector{T},
  wMat::AbstractMatrix{T})::Matrix{T} where {T<:Real}
  # k_local = size(X, 2)
  # n_control = size(x, 1)
  # β = zeros(T, n_target, k_local)
  n_target = size(wMat, 2)
  @inbounds @threads for i in 1:n_target
    w = @view wMat[:, i]
    _β = solve_reg(X, Y, w) # [1, k_local]
    β[i, :] = _β
  end
  return β
end

# AbstractVecOrMat
function GWR(X::AbstractMatrix{T}, Y::AbstractVector{T},
  wMat::AbstractMatrix{T})::Matrix{T} where {T<:Real}
  n_target = size(wMat, 2)
  k_local = size(X, 2)
  β = zeros(T, n_target, k_local)
  GWR!(β, X, Y, wMat)
end





#实现了一个基础的 GWR 求解接口


"
这组函数实现的是：

对每个目标位置 i

取出：

w
i
	​

=wMat[:,i]

然后求：

β
^
	​

i
	​

=(X
T
W
i
	​

X)
−1
X
T
W
i
	​

Y

最后把所有目标位置的局部系数排列成：

β=
	​

β
^
	​

1
T
	​

β
^
	​

2
T
	​

⋮
β
^
	​

n
target
	​

T
	​

	​

	​


即尺寸为：

(n_target, k_local)




"