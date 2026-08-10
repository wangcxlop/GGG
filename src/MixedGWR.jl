


" 整个项目的模块入口文件。 把依赖包、导出函数和各个源代码文件组织到一个名为 MixedGWR 的模块中
可以被理解为整个项目的“总目录”或“主入口”"


module MixedGWR

using LoopVectorization
using ProgressMeter
using LinearAlgebra, Statistics
using Parameters
using Base.Threads
import Base.summary
# using Polyester: @batch

export MGWR, update_weight!
export GWR
export ST_GWR, ST_GWR_fast, ST_GWR_fast!, gwr_neighbors
export GWR_mixed, GWR_mixed_trace, gwr_q, gw_weight_vec, solver_reg, fitted
export fitted, fitted!, predict, summary

export cor

Base.Matrix(x::Vector) = reshape(x, length(x), 1)


include("MGWR.jl")
include("kernel.jl")
include("gw_weight.jl")
include("PrecipitationCorrection.jl")
include("solve_chol.jl")
include("solve_reg.jl")
include("GWR.jl")
include("GWR_calib.jl")
include("deprecated.jl")

include("ST_GWR.jl")

include("GWR_mixed.jl")
include("GWR_mixed_trace.jl")

get_nthread() = Threads.nthreads(:interactive) + Threads.nthreads(:default)

end # module MGWR
