


" 整个项目的模块入口文件。 把依赖包、导出函数和各个源代码文件组织到一个名为 MixedGWR 的模块中
可以被理解为整个项目的“总目录”或“主入口”"


module MixedGWR

using LoopVectorization
using ProgressMeter
using LinearAlgebra, Statistics
using Parameters
using Base.Threads
# using Polyester: @batch

export GWR
export ST_GWR, ST_GWR_fast, ST_GWR_fast!, gwr_neighbors
export fitted, fitted!

export cor
export metric_continuous, metric_event, common_valid_mask, complete_time_mask

Base.Matrix(x::Vector) = reshape(x, length(x), 1)


include("core/fitted.jl")
include("core/metrics.jl")
include("core/kernel.jl")
include("core/gw_weight.jl")
include("core/PrecipitationCorrection.jl")
include("core/solve_chol.jl")
include("core/solve_reg.jl")
include("core/GWR.jl")
include("core/GWR_calib.jl")
include("core/deprecated.jl")

include("core/ST_GWR.jl")

get_nthread() = Threads.nthreads(:interactive) + Threads.nthreads(:default)

end # module MGWR
