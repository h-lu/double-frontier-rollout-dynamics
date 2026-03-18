# 双前沿推进动力学

这个仓库用于研究一个基于 Julia 的慢快动力系统模型。模型的核心对象是：

- advisor 前沿
- doer 前沿
- 监督 / 集成积压
- 采用资本
- 判断资本

当前代码遵循一条很明确的研究顺序：

1. 先使用光滑化、可延拓的主模型；
2. 先在 `2D` 降维模型里判断 fold、Hopf、迟滞与 canard 型跳变是否存在；
3. 只有在低维几何结构比较清楚以后，才回到 `3D` 和 `4D` 系统讨论 MMO 或爆发振荡。

整个项目只使用 Julia，不使用 Python 或 MATLAB。

更完整的研究背景、实验路线和止损规则见
[research-background.md](/Users/wangxq/Documents/AUTO/docs/research-background.md)。

## 模型层次

主脚本实现的是一个光滑化的 `4D` 系统，状态变量为：

- `d`：doer 份额
- `q`：监督 / 集成负担
- `k`：采用资本
- `h`：判断资本

模型的完整方程以及 `2D/3D` 降维系统的显式写法已经补在
[research-background.md](/Users/wangxq/Documents/AUTO/docs/research-background.md) 里。

在完整系统之外，代码还保留了两层降维模型：

- `3D` 模型：固定 `k`
- `2D` 模型：固定 `k` 和 `h`

目前 `2D` 系统是最重要的数值筛查层，因为它最适合先回答下面这些问题：

- 临界流形是否有 S 形几何
- 是否存在 fold 型结构
- 是否出现 Hopf
- 是否存在明显迟滞或双稳态
- 是否存在 canard 型的窄窗口振幅跳变
- 是否进入大振幅松弛振荡

## 仓库结构

- `Project.toml` 与 `Manifest.toml`：Julia 环境
- `docs/research-background.md`：研究背景、实验路线、图像逻辑与论文方向
- `scripts/install_status.jl`：环境与依赖检查
- `scripts/double_frontier_rollout.jl`：主模型、`2D` 平衡点延拓、迟滞扫描与代表性 `4D` 仿真
- `scripts/double_frontier_canard_search.jl`：围绕 `2D` canard / 松弛振荡区域的搜索脚本
- `results/`：当前保存的图与文本结果

## 当前结果文件

当前仓库跟踪的结果都围绕双前沿模型：

- `double_frontier_fold_screen.png`：基于 `beta_n * lambda_D > 4` 的理论前置筛查
- `double_frontier_branch.png`：以推进压力 `E` 为参数的降维分支延拓结果
- `double_frontier_hysteresis.png`：低 / 高初值下的迟滞与路径依赖结果
- `double_frontier_phase.png`：代表性的 `2D` 相图与时序
- `double_frontier_full_timeseries.png`：代表性的 `4D` 长时序
- `double_frontier_summary.txt`：当前降维模型数值结果摘要

## 运行方式

使用本地 Julia：

```bash
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/install_status.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_rollout.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_canard_search.jl
```

## 数值工具栈

- `OrdinaryDiffEq.jl`：刚性 ODE 积分，目前主要使用 `Rodas5P()`
- `BifurcationKit.jl`：平衡点延拓与分岔检测
- `Plots.jl`：图像输出

## 当前研究位置

截至目前，项目已经在部分参数区看到较强的路径依赖，以及降维模型中的 fold / Hopf 候选结构。下一步最关键的问题已经比较明确：

- 振荡起始只是“直接切换到一个大振幅稳定极限环”吗？
- 还是确实存在一个很窄的 canard 型跳变窗口，使系统从极小振幅迅速进入松弛振荡？

当前仓库的设计就是为了优先回答这个 `2D` 问题，再决定是否值得进一步投入 `3D` 折叠奇点、MMO 或 `4D` 爆发振荡。
