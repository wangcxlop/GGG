"""
Apply one kernel across a distance vector.

Split out of `gw_weight!` purely as a function barrier. `GWR_KERNELS` is a `Vector{Function}`,
so `GWR_KERNELS[kernel+1]` is only known to be some `Function` and every call through it in the
loop below would otherwise be a dynamic dispatch with a boxed return - once per control point,
per target. Passing the kernel as an argument lets Julia specialise this method on its concrete
type, so the call resolves at compile time. The arithmetic is unchanged.
"""
function _apply_kernel!(wv::AbstractVector{T}, vdist::AbstractVector{T}, bw::T,
  kerf::F) where {T<:Real,F}
  @inbounds for i in 1:length(vdist)
    wv[i] = kerf(vdist[i], bw)
  end
  return wv
end

"""
Calculate weight vector for GWR.

`buffer`, when given, is scratch of at least `length(vdist)` used for the adaptive bandwidth's
order statistic. The bandwidth is the `bw`-th smallest distance, so `partialsort!` on a copy in
`buffer` returns exactly what `sort(vdist)[Int(bw)]` returned, without allocating a fresh sorted
copy for every target. `vdist` itself is never modified either way.
"""
function gw_weight!(wv::AbstractVector{T}, vdist::AbstractVector{T}, bw::T;
  kernel::Int=0, adaptive::Bool=true,
  buffer::Union{Nothing,AbstractVector{T}}=nothing) where {T<:Real}

  kerf = GWR_KERNELS[kernel+1]  # Julia is 1-indexed

  n = length(vdist)
  if adaptive
    dn = bw / n # bw此时理解为n个点
    if dn <= 1
      k = Int(bw)  # Julia is 1-indexed
      if buffer === nothing
        svdist = sort(vdist)
        bw = svdist[k]
      else
        scratch = view(buffer, 1:n)
        copyto!(scratch, vdist)
        bw = partialsort!(scratch, k)
      end
    else
      bw = dn * maximum(vdist)
    end
  end

  return _apply_kernel!(wv, vdist, bw, kerf)
end

# function gw_weight(vdist::AbstractVector{T}, bw::T;
#   kernel::Int=0, adaptive::Bool=true) where {T<:Real}
#   wv = zeros(T, size(vdist))
#   gw_weight!(wv, vdist, bw; kernel, adaptive)
# end


## 二维数据
function gw_weight!(ws::AbstractMatrix{T}, dist::AbstractMatrix{T}, bw::T;
  kernel::Int=0, adaptive::Bool=true) where {T<:Real}
  n_control, n_target = size(dist)

  # One scratch vector for the whole matrix instead of a fresh sorted copy per column.
  buffer = adaptive ? Vector{T}(undef, n_control) : nothing
  @inbounds for j in 1:n_target
    _w = @view ws[:, j]
    _dist = @view dist[:, j]
    gw_weight!(_w, _dist, bw; kernel, adaptive, buffer)
  end
  return ws
end


function gw_weight(dist::AbstractMatrix{T}, bw::T; kernel::Int=0, adaptive::Bool=true) where {T<:Real}
  wMat = zeros(T, size(dist))
  gw_weight!(wMat, dist, bw; kernel, adaptive)
end


export gw_weight, gw_weight!
