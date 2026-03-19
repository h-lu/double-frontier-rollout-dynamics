# 双前沿推进动力学

本项目研究一个基于 Julia 的慢快动力系统模型，用来刻画组织在 advisor 区、doer 区与中间带之间的动态推进过程。模型同时包含监督负担、采用资本与判断资本三类组织性状态，因此既适合讨论路径依赖与迟滞，也适合进一步探索 canard、松弛振荡、MMO 与爆发振荡等更复杂的动力学机制。

项目的基本策略是：先在低维系统里确认几何结构，再决定是否值得进入更高维、更昂贵的数值分析。完整的研究背景、模型方程、实验路线与止损规则见
[docs/research-background.md](/mnt/d/double-frontier-rollout-dynamics/docs/research-background.md)。

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
- `docs/stage-results-2026-03-19.md`：截至 `2026-03-19` 的阶段性结果文档
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
- `double_frontier_hopf_neighborhood.png`：精修 Hopf 点后的局部周期轨延拓
- `double_frontier_summary.txt`：当前降维模型结果摘要
- `double_frontier_research_note_2026-03-19.md`：截至 `2026-03-19` 的阶段研究总结
- `docs/stage-results-2026-03-19.md`：更接近论文结果节写法的正式阶段文档
- `double_frontier_4d_period2_check.txt`、`double_frontier_4d_period2_boundary.txt`：`4D` 交替大周期的长窗与边界诊断
- `double_frontier_4d_period2_micro_windows.txt`、`double_frontier_4d_theta_micro_slices.txt`、`double_frontier_4d_eta_k_micro_slices.txt`、`double_frontier_4d_G_micro_slices.txt`：`4D` 局部微窗与高维碎裂结构搜索结果

截至当前阶段，更稳的结论已经不只是“存在 fold / Hopf 候选”。现有结果支持这样一条层次化叙述：

- `2D` 已经建立出 `fold + subcritical Hopf + large-cycle backbone`
- `3D` 主要表现为 `strongly subcritical Hopf + single relaxation cycle`
- `4D` 当前最强现象是高维碎裂的 `period-2-like relaxation pattern`

也就是说，项目已经从“是否值得研究”推进到“高维复杂结构具体长什么样”，但当前最强证据仍然指向倍周期调制，而不是 classical `MMO` 或 fully developed `bursting`。

## 运行方式

使用本地 Julia。为了让参数扫描真正利用多核，建议显式开启线程：

```bash
julia --threads=auto --project=. scripts/install_status.jl
julia --threads=auto --project=. scripts/double_frontier_rollout.jl
julia --threads=auto --project=. scripts/double_frontier_canard_search.jl
```

如果希望固定线程数，也可以：

```bash
JULIA_NUM_THREADS=32 julia --project=. scripts/double_frontier_canard_search.jl
```

当前仓库里真正适合并行的是参数扫描与独立仿真；`BifurcationKit` 的平衡点延拓仍然按单线程流程运行。

## 数值工具

- `OrdinaryDiffEq.jl`：刚性 ODE 积分，当前主要使用 `Rodas5P()`
- `BifurcationKit.jl`：平衡点延拓与分岔检测
- `Plots.jl`：图像输出

## 当前阶段

当前仓库已经超过“纯框架搭建”阶段，进入了有明确阶段性结论的探索期。当前最合理的项目状态是：

- `2D`：几何骨架已基本建立，可视为低维筛查成功
- `3D`：确认升维后骨架保留，但没有出现可信 `MMO`
- `4D`：已发现强的 `period-2-like relaxation` 与碎裂参数骨架，是当前主研究前沿

如果继续推进，最值得做的不是继续大范围盲扫，而是围绕 `4D` 局部参数骨架做更系统的几何解释与结果整理。
