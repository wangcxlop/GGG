#set math.equation(
  numbering:"(1)",
  supplement: [式],
)








GWR认为不同地理位置上的关系可能不同，所以在每个位置 i 都拟合一个局部回归模型（核心假设）
$
y_i = beta_0(u_i, v_i) + beta_1(u_i, v_i)x_(1i) + beta_2(u_i, v_i)x_(2i) + dots.c + epsilon_i
$

变量解释：

 $u_i, v_i$ 为位置 i 的经纬度或空间坐标


$beta_k(u_i, v_i)$ 为第$k$个自变量 $x_k$ 在位置$i$附近的局部影响


与其他模组的区别在于其空间加权，当估计某个位置的权重时，通常会给附近的样本更高的权重，而较远处的样本则权重更低


就这样可以得到每个位置的一组局部系数，可进行下一步的分析计算


= 常见权重函数

== 高斯核 Gaussian

$
w_(i,j) = exp[-0.5(d_(i,j) / b)^2]
$

权重随距离平滑下降，但下降永远不为 0 ， 适合空间关系连续变化，连续衰减，不截断，远处样本仍参与

== 指数核 Exponential

$
w_(i,j) = exp[- d_(i,j) / b]
$

比高斯核下降更慢一些，远处样本仍有一定影响。适合空间影响随距离逐渐衰减，但不想完全截断远处样本的情况


== 双平方核 Bisquare
$
w_(i,j) = cases(
  [1 - (d_(i,j) / b)^2]^2 & d_(i,j) < b,
  0 & "否则",
)
$

截断型，超过带宽直接不参与，常用于GWR, 适合强调局部邻域，减少远处样本干扰的情况

== 三次核 Tricube
$
w_(i, j) = cases([1 - (d_(i, j) / b)^3]^3 & d_(i, j)  < b,
  0 & "否则"


)
$

截断型，更强调近邻，近距离权重高，接近带宽边界快速下降为0。更适合强调邻近样本的局部拟合




== 矩形核 Boxcar

$
w_(i, j) = cases(1 & d_(i, j) < b,
  0 & "否则"



)
$

带宽内所有样本权重一样，带宽外的不用。适合粗略局部平均或对距离衰减要求不高的情况。

= 加权最小二乘公式

用于从空间权重算出局部系数的核心公式

相关公式为:(局部系数($hat(beta)(u_i, v_i)$)

$
hat(beta)(u_i, v_i) = (X^upright(T) W_i X^)^(-1) X^upright(T) W_i y
$ <eq:heta>

筛选标准为:

$
min sum_j w_(i, j) (y_j - hat(y)_j)^2
$ <eq:sa>



$
limits(min)_(beta_i)
sum_(j = 1)^n w_(i,j)
lr([y_i - x_j^upright(T) beta_i]) ^ 2
$



由 @eq:sa 来确定找到的局部系数来让所有样本的 *“加权残差平方和”* 最小



使用加权最小二乘公式:

首先，对某一位置 $i$ ,我们已经有了一组权重:

$
w(i,1),w(i,2),w(i,3),w(i,4),w(i,5),dots w(i,n)
$

令$W_i = op("diag") lr([w(i,1),w(i,2),w(i,3),dots.c w(i,n)])$。若展开成完整的对角矩阵则为:

$
W_i = mat(
  w_(i, 1),  0,        dots.c,        0;
  0,         w_(i, 2), dots.c,        0;
  dots.v,    dots.v,   dots.down,  dots.v;
  0,     0 ,   dots.c  ,  w_(i, n)

)

$




接着把这些权重放入最小二乘公式@eq:heta:






= 局部回归方程

$

hat(Y)_i = hat(beta)_0(u_i, v_i) + hat(beta)_1(u_i, v_i) times x_(i, 1) + hat(beta)_2(u_i, v_i) times x_(i, 2) + dots.c + hat(beta)_p(u_i, v_i) times x_(i, p)
$


注:$x_(i, 1), x_(i, 2), dots.c, x_(i, p)$中，意思是第i个位置上第1/2/.../p个自变量



#text(
  font:(
    "Times New Roman",
    "Microsoft YaHei"),
  weight:900,
)[



= 简单流程


空间距离 $d_(i,j) $

↓ 代入权重函数（Gaussian, Exponential, Bisquare, Tricube, Boxcar）

空间权重 $w_(i,j)$

↓组成目标位置 $i$的空间权重矩阵。对角元素为

$
W_i = op("diag")(w_(i, 1),w_(i, 2),dots.c,w_(i, n))
$

↓ 进入加权回归 (加权最小二乘回归, Weighted Least Squares,WLS)(注:带空间权重的加权最小二乘,X,y已准备)

局部系数 $hat(beta)(u_i, v_i)$

↓使用局部回归方程得到局部预测值($p$ 为自变量个数)

$
hat(Y)_i = hat(beta)_0(u_i, v_i) + hat(beta)_1(u_i, v_i) times x_(i, 1) + hat(beta)_2(u_i, v_i) times x_(i, 2) + dots.c + hat(beta)_p(u_i, v_i) times x_(i, p)
$



↓与真实值比较计算残差

$
e_i = Y_i - hat(Y)_i
$

↓计算模型评价指标

得到 $"RMSE"$、$"MAE"$、$"Bias"$、$R^2$







i 是目前要算的局部系数的目标站点或位置

j 是使用的参与该位置局部回归的样本点编号
]


_eg_ 现在有站点1,站点2,站点3,站点4,站点5。现在要计算站点1的局部系数

则 i = 1

j = 1, 2, 3, 4, 5

由计算得到
d(1,1),d(1,2),d(1,3),d(1,4),d(1,5)

↓ 权重函数

接着得到 w(1,1),w(1,2),w(1,3),w(1,4),w(1,5)


$
W_1 = op("diag") lr([
  w_(1,1),w_(1,2),w_(1,3),w_(1,4),w_(1,5)
])
$


构造$X$:
$
X = mat(
  1, "GPM"_1;
  1, "GPM"_2;
  1, "GPM"_3;
  1, "GPM"_4;
  1, "GPM"_5;
)
$ 

构造$y$:
$
y = mat(
  "Obs"_1;
  "Obs"_2;
  "Obs"_3;
  "Obs"_4;
  "Obs"_5
)

$




↓加权回归，加权最小二乘回归(代入@eq:heta)

得到局部系数: $hat(beta)(u_i, v_i)$

$
hat(beta)(u_i, v_i) = mat(
  delim: "[",
  hat(beta)_0(u_1, v_1);
  hat(beta)_1(u_1, v_1)
)
$

其中$hat(beta)_0(u_1, v_1)$ 为站点1的局部截距
$hat(beta)_1(u_1, v_1)$  为站点1的GPM局部系数

可得站点1的局部回归方程:
$
"Obs" approx hat(beta)_0(u_i, v_i) + hat(beta)_1(u_1, v_1) times "GPM"
$

则站点1的局部预测值:

$
hat(Y)_1 = hat(beta)_0(u_1, v_1) + hat(beta)_1(u_1, v_1) times "GPM"_1
$


最后可得站点1的残差为:

$
e_1 = "Obs"_1 - hat(Y)_1
$


= 评价指标

$"RMSE"$、$"MAE"$、$"Bias"$、$r$、$"POD"$、$"FAR"$、$"CSI"$

#table(
  columns:3,
  [],[],[],
  ["RMSE"],[$"RMSE" = sqrt(1 / n sum_(i = 1)^n (hat(Y)_i - Y_i)^2)$],[0],
  ["MAE"],[$"MAE" = 1 / n sum_(i = 1)^n abs(hat(Y)_i - Y_i) $],[0],
  ["Bias"],[$"Bias" = 1 / n sum_(i = 1)^n (hat(Y)_i - Y_i)$],[0],
  ["r"],[$"r" = (
  sum_(i=1)^n
  (Y_i - overline(Y))
  (hat(Y)_i - overline(hat(Y)))
)
/
sqrt(
  sum_(i=1)^n (Y_i - overline(Y))^2
  times
  sum_(i=1)^n
  (hat(Y)_i - overline(hat(Y)))^2
)$],[1],
  ["POD"],[$"POD" = H / (H + M)$],[1],
  ["FAR"],[$"FAR" = F / (H + F)$],[0],
  [CSI],[$"CSI = H / (H + M + F)"$],[1]




)







= 融合


先估计卫星降水的误差，再把误差加回卫星降水进行校正



= 当前版本





↓ 对每个产品构造残差校正目标
$
R^upright(P) = Y_"obs" - Y_"sat"^upright(P)
$

↓ 构造 $X$


$
X = lr([
  1,
  l o n - overline(l o n ),
  l a t - overline(l a t )
])
$

↓ 五折空间交叉验证

↓ 在训练站点上计算空间距离


$
d_(i, j)
$


↓ 代入权重函数并扫描带宽

↓ 在训练集内进行 LOOCV 调参

得到最优 kernel / bandwidth / adaptive 参数

↓ 用训练站点残差拟合 ST-GWR / GWR 残差模型


$
hat(beta)(u_i, v_i) = (X^upright(T) W_i X)^(-1) X^upright(T) W_i R^upright(P)
$



↓ 在验证站点预测卫星误差

$
hat(R)^upright(P_i) = 
hat(beta)_0(u_i, v_i)
+ hat(beta)_1(u_i, v_i) times (l o n_i - overline(l o n))
+ hat(beta)_2(u_i, v_i) times (l a t_i - overline(l a t))
$

↓ 校正该产品的卫星降水

$
Y_("corr",i)^P = Y_("sat", i)^P + hat(R)_i^P
$

↓ 与验证站点观测值比较

$
e_i^P
= Y_("obs", i) - Y_("corr", i)^P
$

↓ 计算评价指标

$"RMSE"$、$"MAE"$、$"Bias"$、$r$、$"POD"$、$"FAR"$、$"CSI"$

↓ 汇总五折结果


注:P可以依次取 FY4B、GPM、GSMaP


























