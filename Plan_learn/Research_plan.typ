= 现有结构
计算权重(把距离数据按指定核函数和带宽，转换成 GWR 用的权重。)\
F:\G\GeoWeightedRegression.jl\src\gw_weight.jl  


在已经获得空间权重矩阵 wMat 的前提下，对每个样本位置执行局部加权最小二乘，估计局部系数，并计算拟合误差、系数标准误、帽子矩阵迹和 AIC
F:\G\GeoWeightedRegression.jl\src\GWR_calib.jl


计算混合地理加权回归模型帽子矩阵迹的函数
F:\G\GeoWeightedRegression.jl\src\GWR_mixed_trace.jl


拟合 Mixed GWR(混合地理加权回归) 的主函数\
F:\G\GeoWeightedRegression.jl\src\GWR_mixed.jl


对每个目标位置 i 取出 $w_i = w M a t [:,i]$ 然后求 $hat(beta_i) = (X^upright(T) W_i X)^(-1) X^upright(T) W_i Y$, 最后把所有目标位置的局部系数排列为$ beta = mat(
  delim:"[",
  hat(beta)_1^T;
  hat(beta)_2^T;
  dots.v;
  hat(beta)_(n_"target")^T
) $
得到尺寸为 (n_target, k_local) 的局部系数矩阵\
F:\G\GeoWeightedRegression.jl\src\GWR.jl


几种核函数公式代码(gaussion, exponential, bisquare, tricube, boxcar)
F:\G\GeoWeightedRegression.jl\src\kernel.jl


这段代码定义了一个完整的 MGWR 模型对象，以及权重更新、拟合值计算、预测和模型评价方法。需要先明确：从这份代码的结构看，MGWR 表示 Mixed GWR，而不是 Multiscale GWR。原因是它把变量划分为：
x1：局部变量\
x2：全局变量\
只有一个局部带宽 bw。\ 真正的 Multiscale GWR 通常为每个解释变量设置独立带宽。
F:\G\GeoWeightedRegression.jl\src\MGWR.jl


整个项目的模块入口文件。 把依赖包、导出函数和各个源代码文件组织到一个名为 MixedGWR 的模块中可以被理解为整个项目的“总目录”或“主入口
F:\G\GeoWeightedRegression.jl\src\MixedGWR.jl


实现的是基于 Choleshy 分解的加权最小二乘求解器， 主要用于GWR或ST-GWR中重复计算局部加权回归系数\
F:\G\GeoWeightedRegression.jl\src\solve_chol.jl


用于求解GWR中某一个目标位置的局部加权回归系数(并提供了两个求解函数)
F:\G\GeoWeightedRegression.jl\src\solve_reg.jl



这段代码实现的是一种批量时序 GWR:\
对每个空间目标点单独建立一个局部加权回归；
Y 的每一列代表一个时间点；
同一目标点的所有时间点共用同一个空间权重向量；
一次矩阵求解同时得到该目标点所有时间的局部系数；
最终输出每个目标点、每个时间点的预测值\
F:\G\GeoWeightedRegression.jl\src\ST_GWR.jl




= 输入输出
F:\G\GeoWeightedRegression.jl\src\gw_weight.jl  
一维版本 gw_weight.jl:2：输入是距离向量 vdist、带宽 bw，以及可选参数 kernel 和 adaptive；输出是同一个权重向量 wv，函数会把每个距离通过核函数转换成权重后写入 wv。
二维版本 gw_weight.jl:23：输入是距离矩阵 dist、带宽 bw，以及可选参数 kernel 和 adaptive；输出是权重矩阵 ws，按列对每一组距离计算权重。
非原地版本 gw_weight.jl:33：输入是距离矩阵 dist 和带宽 bw；输出是新创建的权重矩阵 wMat。




F:\G\GeoWeightedRegression.jl\src\GWR_calib.jl
GWR_calib(model::MGWR) 的输入是一个 MGWR 模型对象，它会从里面取出 x1、y 和 wMat，然后调用真正的校准函数。GWR_calib(x, y, wMat) 的输入分别是自变量矩阵 x、因变量向量 y、以及权重矩阵 wMat。\
输出是一个命名元组，里面包含局部回归结果和评估指标：β、se_β、n、r、RSS、RMSE、σ、trace、AIC。


F:\G\GeoWeightedRegression.jl\src\GWR_mixed_trace.jl
算混合地理加权回归的“帽子矩阵迹”\ 
GWR_mixed_trace(x1, x2, dMat, bw; kernel=0, adaptive=false)输入是：
x1：局部变量矩阵
x2：全局变量矩阵
dMat：距离矩阵
bw：带宽
可选参数 kernel、adaptive\
输出是：
一个标量 T，表示混合 GWR 模型帽子矩阵的迹，也就是有效自由度相关的量。\
GWR_mixed_trace(model::MGWR)
输入是：
一个 MGWR 模型对象
输出是：
同样返回帽子矩阵迹的标量。






F:\G\GeoWeightedRegression.jl\src\GWR_mixed.jl
mixed(x1, x2, y, dMat, dMat_rp, bw; kernel=0, adaptive=false)
输入是：

x1：局部变量矩阵
x2：全局变量矩阵
y：因变量向量
dMat：用于局部拟合的距离矩阵
dMat_rp：用于重拟合的距离矩阵
bw：带宽
可选参数 kernel、adaptive
输出是：
一个命名元组，包含 local 和 global 两部分系数结果，也就是局部回归系数和全局回归系数。

WR_mixed(model::MGWR{T})
输入是：
一个 MGWR 模型对象
输出是：
同样返回 local 和 global 两部分系数
另外还会把结果写回到 model.β1 和 model.β2 里，等于同时修改模型内部状态。



F:\G\GeoWeightedRegression.jl\src\GWR.jl
#image("/assets/image.png")



F:\G\GeoWeightedRegression.jl\src\MGWR.jl
#image("/assets/image-1.png")




F:\G\GeoWeightedRegression.jl\src\MixedGWR.jl
#image("/assets/image-2.png")




F:\G\GeoWeightedRegression.jl\src\solve_chol.jl
#image("/assets/image-3.png")


F:\G\GeoWeightedRegression.jl\src\solve_reg.jl
#image("/assets/image-4.png")





F:\G\GeoWeightedRegression.jl\src\ST_GWR.jl
#image("/assets/image-5.png")





