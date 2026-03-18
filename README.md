# 双前沿推进动力学

本项目研究一个基于 Julia 的慢快动力系统模型，用来刻画组织在 advisor 区、doer 区与中间带之间的动态推进过程。模型同时包含监督负担、采用资本与判断资本三类组织性状态，因此既适合讨论路径依赖与迟滞，也适合进一步探索 canard、松弛振荡、MMO 与爆发振荡等更复杂的动力学机制。

项目的基本策略是：先在低维系统里确认几何结构，再决定是否值得进入更高维、更昂贵的数值分析。完整的研究背景、模型方程、实验路线与止损规则见
[research-background.md](/Users/wangxq/Documents/AUTO/docs/research-background.md)。

## 研究问题

当前数值工作优先回答以下问题：

- `2D` 降维模型是否存在 fold、Hopf、双稳态与迟滞？
- 振荡起始是否表现为 canard 型窄窗口跳变？
- 大振幅周期是否具有典型的松弛振荡几何？
- 更高维系统是否值得继续追踪折叠奇点、MMO 与爆发振荡？

## 模型结构

主脚本实现的是一个光滑化的 `4D` 系统，状态变量为：

- `d`：doer 份额
- `q`：监督 / 集成负担
- `k`：采用资本
- `h`：判断资本

在完整系统之外，代码还保留了两层降维模型：

- `2D` 模型：固定 `k` 与 `h`，研究 `(d, q)`
- `3D` 模型：固定 `k`，研究 `(d, q, h)`

其中，`2D` 系统是当前最关键的筛查层。它最适合先判断临界流形是否呈现 S 形、是否存在 fold 与 Hopf、是否出现迟滞，以及是否有 canard 型振幅跳变和松弛振荡。

## 仓库内容

- `Project.toml` 与 `Manifest.toml`：Julia 环境
- `docs/research-background.md`：研究背景、完整方程、实验设计与论文路线
- `scripts/install_status.jl`：依赖检查
- `scripts/double_frontier_rollout.jl`：主模型、`2D` 平衡点延拓、迟滞扫描与代表性 `4D` 仿真
- `scripts/double_frontier_canard_search.jl`：围绕 `2D` canard / 松弛振荡区域的搜索脚本
- `results/`：当前图像与文本结果

## 当前结果

当前仓库已保存一组围绕双前沿模型的代表性结果：

- `double_frontier_fold_screen.png`：基于 `beta_n * lambda_D > 4` 的理论前置筛查
- `double_frontier_branch.png`：以推进压力 `E` 为参数的 `2D` 分支延拓结果
- `double_frontier_hysteresis.png`：低 / 高初值下的迟滞与路径依赖
- `double_frontier_phase.png`：代表性的 `2D` 相图与时序
- `double_frontier_full_timeseries.png`：代表性的 `4D` 长时序
- `double_frontier_summary.txt`：当前降维模型结果摘要

这些结果说明：部分参数区已经表现出较强的路径依赖，以及 fold / Hopf 候选结构。当前最重要的后续问题是，振荡起始究竟只是进入一个大振幅稳定极限环，还是确实存在 canard 型窄窗口跳变。

## 运行方式

使用本地 Julia：

```bash
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/install_status.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_rollout.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_canard_search.jl
```

## 数值工具

- `OrdinaryDiffEq.jl`：刚性 ODE 积分，当前主要使用 `Rodas5P()`
- `BifurcationKit.jl`：平衡点延拓与分岔检测
- `Plots.jl`：图像输出

## 当前阶段

当前仓库更接近一个已经搭建好的研究框架，而不是完成版论文代码。现阶段的目标很明确：优先把 `2D` 系统里的几何结构、迟滞区间与振荡起始机制看清楚，再决定是否值得投入 `3D` 折叠奇点、MMO 或 `4D` 爆发振荡。
