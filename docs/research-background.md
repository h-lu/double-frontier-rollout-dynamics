# 研究背景

这份文档用于整理当前双前沿推进模型背后的研究计划。它的定位是“代码仓库的研究背景说明”，不是论文终稿。

## 一、研究目标

这个项目最核心的目标，不是把一个 `4D` 系统堆得很复杂，而是先回答几个更基础、也更关键的问题：

- 降维模型里到底有没有 fold 几何？
- 是否存在 Hopf？
- 有没有稳定极限环？
- 振荡起始是不是 canard 型的窄窗口跳变？
- 更高维的降维系统是否真的支持 MMO 或爆发振荡？

所以工作流必须先从降维模型开始，等低维几何结构清楚以后，再回到完整系统。

## 二、Julia 数值路线

最初的讨论里提到过 SciPy、MatCont 和 AUTO 这类经典工具，它们代表的是一条标准的慢快系统数值路线：刚性直接积分和延拓并用。在当前仓库里，这条路线全部改写成 Julia 版本：

- `OrdinaryDiffEq.jl`：承担刚性初值问题积分
- `BifurcationKit.jl`：承担平衡点和分岔延拓
- `Plots.jl`：负责图像输出

虽然工具换成了 Julia，但数值原则不变：

- 小 `epsilon` 会使系统变得刚性
- `min` 和正部函数最好先光滑化
- 只靠粗网格扫参数很容易漏掉 canard 窗口
- 延拓和直接积分必须结合起来看，而不是只做其中之一

## 三、主模型

当前主模型是一个光滑化的 `4D` 系统，变量为：

- `d`：doer 份额，快变量
- `q`：监督 / 集成负担，中间尺度变量
- `k`：采用资本，慢变量
- `h`：判断资本，慢变量

这个模型保留了研究里最关键的经济含义：

- advisor 前沿
- doer 前沿
- 中间带
- 推进带来的监督负担
- 采用资本的积累
- 判断资本的积累与侵蚀

### 0. 完整模型方程

为了后续写论文、对照代码和扩展降维模型，这里先把完整模型集中写出来。

#### 光滑化函数

$$
\operatorname{smin}_{\sigma}(x,y)
=-\frac{1}{\sigma}\log\left(e^{-\sigma x}+e^{-\sigma y}\right),
\qquad \sigma \gg 1
$$

$$
P_{\sigma}(r)
=\frac{1}{\sigma}\log\left(1+e^{\sigma r}\right),
\qquad \sigma \gg 1
$$

它们分别近似：

- \(\min(x,y)\)
- \(r^+ = \max(r,0)\)

#### advisor 前沿

$$
a(h)=F_A(\alpha_0+\alpha_z z+\alpha_h h+\alpha_G G)
$$

$$
F_A(s)=\frac{1}{1+e^{-\lambda_A(s-c_A)}}
$$

#### 监督覆盖率

$$
\phi(h,q)=\frac{1}{1+\exp[-(\omega_0+\omega_G G+\omega_h h-\omega_q q)]}
$$

#### doer 目标份额

$$
d^\star(d,q,k,h)
=
\operatorname{smin}_{\sigma}\left(
F_D(\theta_0+\theta_z z+\theta_h h+\theta_G G+\theta_E E+\theta_k k+\theta_{\phi}\phi(h,q)-\gamma q+\beta_n d),
\ a(h)
\right)
$$

$$
F_D(s)=\frac{1}{1+e^{-\lambda_D(s-c_D)}}
$$

#### 完整的 4D 动态系统

$$
\varepsilon \dot d = d^\star(d,q,k,h)-d
$$

$$
\dot q
=
\chi_I P_\sigma(d^\star-d)
+\chi_R d^\rho
-(\nu_0+\nu_G G+\nu_h h)q,
\qquad \rho>1
$$

$$
\dot k
=
\eta_k\left[
\psi_M P_\sigma(a(h)-d)
+\psi_D \phi(h,q)d
-\delta_k k
\right]
$$

$$
\dot h
=
\eta_h\left[
(\vartheta_M P_\sigma(a(h)-d)+\vartheta_D \phi(h,q)d+\vartheta_H(1-a(h)))(1-h)
-
(\mu_0+\mu_1(1-\phi(h,q))d)h
\right]
$$

这组方程就是仓库里 `4D` 主模型的数学原型。

### 1. 光滑化

为了避免硬拐角，模型里用两类光滑函数：

- `smin_sigma(x, y)` 近似 `min(x, y)`
- `softplus_sigma(r)` 近似 `max(r, 0)`

它们直接写在
[double_frontier_rollout.jl](/Users/wangxq/Documents/AUTO/scripts/double_frontier_rollout.jl)
里，作用就是让延拓和刚性积分更稳定。

### 2. 前沿与覆盖率

主模型里有三个最核心的非线性对象：

1. advisor 前沿 \(a(h)\)
2. 监督覆盖率 `phi(h, q)`
3. doer 目标份额 \(d^\star(d,q,k,h)\)

经济含义可以概括成：

- `z`、`G`、`h` 越高，advisor 前沿越往右
- `G`、`h` 越高，监督覆盖率越好
- `E`、`k`、短期反馈 `beta_n` 越高，doer 推进压力越强
- 积压变量 `q` 越大，对进一步推进的抑制越强

## 四、三套降维模型

整个项目始终围绕三套嵌套系统展开。

### 1. 2D 模型

固定 `k = kbar`、`h = hbar`，研究 `(d, q)`。

对应的降维系统写成：

$$
\varepsilon \dot d = d^\star(d,q;\bar k,\bar h)-d
$$

$$
\dot q
=
\chi_I P_\sigma(d^\star-d)
+\chi_R d^\rho
-(\nu_0+\nu_G G+\nu_h \bar h)q
$$

这是最先做筛查的层次，用来判断：

- 临界流形是否 S 形
- 是否存在 fold 型结构
- 是否出现 Hopf 候选
- 是否有迟滞和双稳态
- 是否存在 canard 型振幅跳变
- 是否进入松弛振荡

### 2. 3D 模型

固定 `k = kbar`，研究 `(d, q, h)`。

对应的降维系统可以记成：

$$
\varepsilon \dot d = d^\star(d,q,\bar k,h)-d
$$

$$
\dot q = g(d,q,h)
$$

$$
\dot h = \tilde{\eta}_h H(d,q,h)
$$

这是第一次可能出现折叠奇点和 MMO 的层次，也是后面几何分析最重要的过渡模型。

### 3. 4D 模型

完整的 `4D` 系统只在低维结构已经看清楚以后再重点使用。

这一层主要关心：

- 推进波
- `(k, h)` 带来的长记忆效应
- 静默阶段与振荡阶段之间的慢切换
- 爆发振荡

## 五、参数盒

当前项目使用的是“数值勘探参数盒”，不是经验校准。

### 1. 基础值

计划中建议的基础值包括：

- `lambda_A = 8`
- `alpha0 = -1`
- `alpha_z = alpha_h = alpha_G = 1`
- `theta0 = -2`
- `theta_z = 1`
- `theta_h = theta_G = 0.8`
- `theta_E = 1.2`
- `theta_k = theta_phi = 1`
- `nu0 = 0.3`
- `delta_k = 0.2`
- `mu0 = 0.05`
- `mu1 = 1`
- `omega0 = -1`
- `omega_G = 1.2`
- `omega_h = 1`
- `omega_q = 1.5`

### 2. 第一轮重点扫描参数

最重要的扫描参数是：

- `beta_n`
- `lambda_D`
- `E`
- `G`
- `z`
- `gamma`
- `chi_R`
- `epsilon`
- `eta_h`

同时还需要考虑一些离散设定：

- `rho in {1.5, 2, 3}`
- `chi_I in {0, 0.5, 1}`
- `sigma in {20, 50, 100, 200}`

### 3. 默认初值

建议的基础初值是：

- `d(0) = 0.05`
- `q(0) = 0`
- `k(0) = 0.1`
- `h(0) = 0.6`

除此之外，还必须加入高采用和低采用初值，用来检验双稳态。

## 六、实验路线

整个数值路线是分层推进的。

### 实验 0：数值稳定性检查

目的：

- 确认平滑化没有制造虚假的动力学
- 确认结论不是容差幻觉
- 确认刚性求解器工作正常

输出：

- 不同容差下的时序图
- 不同 `sigma` 下的相图
- 求解器诊断信息

### 实验 1：无 fold 的基线区

目的：

- 先看 `beta_n` 很小时模型是否只是普通慢流形调整
- 把真正的非线性几何与简单的滞后区分开

输出：

- `2D` 相图
- 零流线
- 对 `E` 的平衡点延拓

### 实验 2：fold 前置筛查

在忽略 advisor cap 的局部近似下，logistic 情形给出一个非常有用的必要条件：

`beta_n * lambda_D > 4`

目的：

- 在进入大规模数值扫描前，先缩小参数空间

输出：

- 理论上的 fold 可行区
- 数值上的三平衡点区域
- 代表性的临界流形

### 实验 3：advisor cap 是否剪掉 S 形

目的：

- 检查 advisor 前沿是否会直接截断降维模型的右支
- 证明双前沿结构确实改变了系统拓扑，而不仅仅是换了记账方式

输出：

- cap 与 fold 位置的比较
- `(G, h)` 区域图
- 被剪掉和未被剪掉时的相图对比

### 实验 4：双稳态与迟滞

目的：

- 先拿下最稳健、最容易站住脚的动力学结论

输出：

- 单参数延拓分支
- 前向 / 后向 sweep
- 高低初值收敛对比

### 实验 5：Hopf 与极限环

目的：

- 判断降维模型是否真的有周期解

输出：

- Hopf 检测
- 极限环延拓
- 极值包络与周期图

### 实验 6：canard 窗口精扫

目的：

- 判断系统是否存在“极小振幅振荡在极窄窗口中突然放大成大振幅松弛震荡”的现象

输出：

- 局部放大的振幅图
- 局部放大的周期图
- 小振幅 / canard 型 / 大振幅周期轨对比

### 实验 7：松弛振荡

目的：

- 描述 canard 窗口之后的大振幅周期几何

输出：

- 大振幅周期的相图
- 快跳与慢漂移并存的时序图

### 实验 8：3D 折叠奇点与 MMO

目的：

- 判断 `(d, q, h)` 降维系统里是否真的存在折叠节点型的几何组织机制

输出：

- fold curve
- 折叠奇点条件检查
- 长时序图
- SAO 计数图

### 实验 9：4D 爆发振荡

目的：

- 理解 `(k, h)` 的慢漂移怎样把系统带入和带出振荡区

输出：

- 长窗口 `4D` 时序图
- `(k, h)` 慢变量轨道
- 快子系统骨架与完整轨道投影

### 实验 10：中间带的不可替代性

目的：

- 说明相同的最终 doer 份额，不同推进路径也可能积累出完全不同的 `k` 和 `h`

输出：

- 路径比较时序图
- `m-h` 图
- `m-k` 图
- matched-final-`d` 比较图

## 七、图像与现象的对应关系

不同图不是重复，而是分别回答不同问题：

- 时间序列图：看快跳、慢漂移、过冲、爆发振荡
- `(d, q)` 相图：看零流线、fold、Hopf 前后周期轨、松弛回路
- 临界流形图：看 S 形、fold、cap clipping
- 单参数延拓图：看双稳态、迟滞、Hopf、周期包络
- 双参数延拓图：看 fold / Hopf 边界
- 局部放大振幅图：看 canard 窗口
- fold curve 与折叠奇点图：看 MMO 的几何骨架
- 峰值返回图 / SAO 计数图：看 MMO 模式
- 快骨架 + 完整轨道：看爆发振荡的组织机制

## 八、止损规则

这些规则和实验本身同样重要。

1. 如果在合理参数盒里始终找不到 fold，就不要继续追 canard 或 MMO，直接转向慢流形、路径依赖和 tipping。
2. 如果有 fold，但没有可靠的 Hopf 或周期轨，就主打迟滞和多稳态，而不是强卖 canard。
3. 如果 `2D` 里有 canard 型行为，但 `3D` 里没有折叠奇点结构，就把数学主线定位在二维 canard / 松弛振荡论文。
4. 只有在高维证据非常明确时，才把 MMO 或爆发振荡写进题目和摘要。

## 九、可能的论文路线

按结果不同，论文方向也不同：

- 如果只看到慢流形、迟滞或 tipping：
  更适合经济学理论文
- 如果看到 fold + canard 型起始 + 松弛振荡：
  可以写 `2D` 应用数学文
- 如果看到折叠节点 + MMO：
  这是最强的应用数学主线
- 如果明确看到爆发振荡：
  可以作为长文后半部，或者拆成第二篇

## 十、当前仓库所处位置

当前代码库还处在这条路线的前半段：

- `2D` 降维系统已经实现，并可直接做延拓
- `4D` 系统已经实现，可用于代表性仿真
- 部分参数区已经能看到较强的路径依赖
- canard / 松弛振荡的搜索仍然处于探索阶段，还不能当作最终结论

所以现在的仓库更适合被理解为一个认真搭建好的数值研究框架，而不是已经完成证明或论文封装的最终版本。

## 十一、参考说明

最初的规划讨论引用过几类经典资料：

- [SciPy `solve_ivp` 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
- [Guckenheimer / MMO 综述](https://pi.math.cornell.edu/~gucken/PDF/mmo_review.pdf)
- [AUTO-07P 概览](https://archive-dsweb.siam.org/Software/auto-07p.html)
- [canard explosion 参考条目](https://www.researchgate.net/publication/247384460_Relaxation_Oscillation_and_Canard_Explosion?utm_source=chatgpt.com)

这些资料在这里保留为研究背景参考。仓库本身的实现已经全部切换为 Julia。
