# 双前沿推进动力学：阶段性结果

日期：2026-03-19

## 摘要

本阶段的数值结果表明，双前沿推进模型在不同维度上呈现出一条清晰的结构链条。

在 `2D` 降维系统中，已经确认 `fold / branch-point`、`subcritical Hopf`、迟滞，以及从 Hopf 邻域延伸到大振幅的周期轨 backbone。`2D` 因而已经足够支撑“低维筛查成功”的结论。

在 `3D` 系统中，`2D` 的核心骨架没有消失，但最主要的动力学仍然是强亚临界 Hopf 与单一松弛振荡分支，而不是可信的 `MMO`。

在 `4D` 系统中，最强现象不是 classical `bursting`，而是一个高维碎裂的 `period-2-like relaxation` 参数骨架。局部参数放大表明，这个大振幅交替周期并不落在平滑 continuation surface 上，而是由薄层、细带、孔洞和微岛共同构成。

## 1. 研究目标

本项目的目标不是一开始就在高维系统中盲目搜索复杂振荡，而是先在低维中确认几何骨架，再逐层判断升维后是否值得继续投入更昂贵的数值分析。

当前阶段优先回答三个问题：

1. `2D` 中是否存在足够清楚的 `fold / Hopf / hysteresis / large-cycle` 骨架。
2. 这些结构在 `3D` 中是否保留，以及是否自然进入 `MMO`。
3. `4D` 中是否出现比 `3D` 更强的多时标调制，并且这种调制在参数空间中的几何形态是什么。

## 2. 2D 结果

### 2.1 主要数值事实

- `fold / branch-point` 候选位于 `E ≈ 1.2116721025`
- 粗 Hopf 位于 `E ≈ 0.9772856335`
- 精修 Hopf 位于 `E ≈ 0.9747525391`
- Hopf normal form 类型为 `SubCritical`
- 周期轨 continuation 已从 Hopf 邻域延伸到
  `E ∈ [0.9744884922, 0.9753187708]`
- 周期轨 `d` 振幅达到 `0.9118387`
- 外侧 branch 上出现 period-doubling

### 2.2 解释

这些结果意味着，`2D` 的大振幅振荡不再只是时域积分里“突然冒出来”的现象。现在可以明确地把它与 Hopf 邻域的小周期 branch 连起来，形成一条具体的周期轨 backbone。

因此，当前对 `2D` 最稳的描述是：

`2D` 已经确认存在 `fold + subcritical Hopf + hysteresis + large relaxation branch`。

这并不等于已经给出严格的 canard 证明，但对于研究路线判断来说，`2D` 已经足够成功。

## 3. 3D 结果

### 3.1 主要数值事实

- `3D fold / bp`：`E ≈ 0.9839955747`
- `3D Hopf`：`E ≈ 1.6729540451`
- Hopf 邻域振荡塌缩窗口：
  `E ≈ 1.6557040451 -> 1.6579540451`
- `3D` 振荡区最大 period CV：`0.012466`
- 中间带占比最大值：`0.001400`

### 3.2 解释

`3D` 确认了一个重要事实：引入动态 `h` 后，`2D` 的骨架并没有被洗掉。但与此同时，`3D` 并没有给出强 `MMO` 证据。

在当前扫描下，`3D` 更像：

- 强亚临界 Hopf
- 单一大振幅松弛振荡 branch
- 很窄的振荡塌缩窗口

因此，`3D` 的主要意义是结构连续性，而不是直接产出 mixed-mode 现象。

## 4. 4D 结果

### 4.1 从调制到稳定交替大周期

`4D` 中最强候选首先表现为很高的 inter-burst 和 burst-width 变异。进一步做事件级诊断后可以确认，这并不是一般意义上的“随机包络调制”，而是非常强的 `period-2-like relaxation pattern`。

代表性的长窗诊断给出：

- odd/even period means = `[5.050674, 28.845005]`
- odd/even width means = `[27.058366, 3.372802]`
- lag-1 相关接近 `-1`
- lag-2 相关接近 `+1`

这说明系统已经进入了稳定的长短大周期交替，而不是弱调制或噪声样波动。

### 4.2 不是 smooth surface，而是 fractured scaffold

更关键的发现来自局部参数放大。

在 `(theta0, E)` 平面里，强 `period-2` 区域初看像一块 tongue-like wedge。但继续放大后发现，这个 wedge 的脆弱角会分裂成很多窄窗；在 `ΔE = 1e-4` 的尺度上，强的大振幅 `period-2` 会反复消失和重现。

继续对 `theta0`、`eta_k` 和 `G` 做微切片后，可以得到更完整的局部几何图景：

- 有些切片更像离散微岛
- 有些切片更像细条纹
- 有些细条纹内部又被许多孔洞打断
- 沿第三、第四参数方向，这种碎裂结构仍然保留

因此，对当前 `4D` 结果最准确的描述不是“找到一张 continuation surface”，而是：

强的大振幅 `period-2 relaxation` 存活在一个高维碎裂参数骨架上。局部上，这个骨架由短 slab、细带、孔洞和微岛共同构成。

## 5. 当前没有得到的结论

到本阶段为止，下面这些更强结论还没有得到支持：

- 还不能说已经找到 classical `MMO`
- 还不能说已经找到 fully developed `bursting`
- 还不能说高阶倍周期级联已经成为当前主要动力学
- 还不能说 `4D` 的碎裂结构已经被严格几何分类

因此，最稳的表述应当保持克制：

当前最强现象是高维碎裂的 large-amplitude `period-2 relaxation scaffold`，而不是 bursting attractor。

## 6. 阶段结论

如果把当前工作压缩成最核心的研究判断，可以概括为四句：

1. `2D` 已经建立了清楚的低维几何骨架。
2. `3D` 说明这条骨架在升维后仍然存在，但没有自然发展成 `MMO`。
3. `4D` 出现了更复杂的慢快调制，但最强证据指向的是大振幅 `period-2`，不是 bursting。
4. 这个 `period-2` 结构本身在参数空间中高度碎裂，说明模型的高维慢快组织远比“单一分岔曲线”更复杂。

## 7. 后续建议

如果继续研究，最值得做的是两类工作。

第一类是整理工作。当前已经有足够多的结果，需要进一步压缩成论文式的图文叙述，而不是继续无边界地扩散参数搜索。

第二类是解释工作。接下来真正重要的问题，不再是“还能不能扫出更细的窗”，而是：

- 为什么 `4D` 会形成碎裂的 `period-2` 参数骨架？
- 哪些慢变量耦合负责把平滑 branch 打碎成 slab、stripe 与 micro-island？

如果未来仍想继续寻找真正的 `bursting`，也应该基于这条几何线索有针对性地改模型耦合，而不是在当前参数块附近继续机械细扫。

## 8. 文件导航

建议按下面顺序阅读：

1. 总摘要：
   [results/double_frontier_summary.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_summary.txt)
2. 阶段研究笔记：
   [results/double_frontier_research_note_2026-03-19.md](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_research_note_2026-03-19.md)
3. `2D` 关键图：
   [results/double_frontier_hopf_neighborhood.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_hopf_neighborhood.png)
   [results/double_frontier_branch.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_branch.png)
4. `3D` 关键图：
   [results/double_frontier_3d_hopf_zoom.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_3d_hopf_zoom.png)
   [results/double_frontier_3d_regime.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_3d_regime.png)
5. `4D` 关键图和微窗分析：
   [results/double_frontier_4d_modulated_candidate.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_modulated_candidate.png)
   [results/double_frontier_4d_alternation.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_alternation.png)
   [results/double_frontier_4d_period2_check.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_period2_check.txt)
   [results/double_frontier_4d_period2_micro_windows.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_period2_micro_windows.txt)
   [results/double_frontier_4d_theta_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_theta_micro_slices.txt)
   [results/double_frontier_4d_eta_k_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_eta_k_micro_slices.txt)
   [results/double_frontier_4d_G_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_G_micro_slices.txt)
