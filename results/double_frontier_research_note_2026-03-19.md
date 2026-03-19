# 双前沿推进动力学阶段总结

日期：2026-03-19

## 一句话结论

当前最强结论不是 `MMO` 或 `bursting`，而是：

- `2D` 已经建立出完整的 `fold + subcritical Hopf + large-cycle backbone`
- `3D` 主要保留为 `strongly subcritical Hopf + single relaxation cycle`
- `4D` 里出现了很强的 `period-2-like relaxation modulation`
- 这条 `period-2` 结构在局部参数空间里不是光滑曲面，而是一个由薄层、细带、孔洞和微岛构成的碎裂参数骨架

## 当前研究状态

- `2D`：基本做完。几何骨架已经清楚，足够支撑“低维筛查成功”的结论。
- `3D`：筛查结论也比较清楚。没有看到可信 `MMO`，但确认了 `2D` 的核心骨架没有在升维后消失。
- `4D`：是当前最值得继续的层，但当前找到的最强现象仍然是高维碎裂的 `period-2 relaxation`，而不是 fully developed bursting。

## 2D 结论

核心数值证据：

- `fold / branch-point` 候选：`E ≈ 1.2116721025`
- 粗 Hopf：`E ≈ 0.9772856335`
- 精修 Hopf：`E ≈ 0.9747525391`
- Hopf normal form：`SubCritical`
- 周期轨 branch 已从 Hopf 邻域继续追出，参数范围
  `E ∈ [0.9744884922, 0.9753187708]`
- 周期轨 `d` 振幅可达 `0.9118387`
- branch 上检测到 period-doubling

解释：

- `2D` 不再只是“看到 Hopf 然后时域里突然跳到大振幅”的黑箱现象。
- 现在可以明确说，大振幅周期轨是从一个真实的 Hopf 邻域周期轨 backbone 连出去的。
- 因此 `2D` 的工作已经足够支持：
  `fold + subcritical Hopf + hysteresis + large relaxation branch`

还没做的事情：

- 没有把 `2D` 写成严格的 canard 证明。
- 但作为研究路线判断，`2D` 已经足够成功。

## 3D 结论

核心数值证据：

- `3D fold / bp`：`E ≈ 0.9839955747`
- `3D Hopf`：`E ≈ 1.6729540451`
- Hopf 邻域振荡塌缩窗口：
  `E ≈ 1.6557040451 -> 1.6579540451`
- `3D` 直接积分中的 low/high mean-d gap 最大只有 `0.015327`
- `3D` 振荡区最大 period CV 只有 `0.012466`
- 中间带占比最大只有 `0.001400`

解释：

- `3D` 更像“强亚临界 Hopf 下挂着单一 relaxation cycle branch”，而不是 `MMO`。
- 没有看到小峰夹大峰，也没有看到明显的 mixed-mode 序列。
- 因此 `3D` 的主要价值不是直接产出 `MMO`，而是说明：
  `2D` 的骨架在升维后仍然成立。

## 4D 结论

### 第一层：不是 regular single cycle

几何参数搜索找到的最强候选在：

- `theta0 = -4.9`
- `beta_n = 1.5`
- `E = 1.07`
- `G = 1.0`
- `epsilon = 0.01`
- `eta_h = 0.002`
- `eta_k = 0.003`

其长窗诊断表现为：

- `inter-burst CV = 0.703592`
- `burst-width CV = 0.782510`
- `crossing ratio = 1.0`

这说明它已经明显不是规则单周期。

### 第二层：最强现象是稳定的 period-2-like relaxation

事件级诊断显示：

- odd/even period means = `[5.050674, 28.845005]`
- odd/even width means = `[27.058366, 3.372802]`
- lag-1 correlations 接近 `-1`
- lag-2 correlations 接近 `+1`

这比“只是有 envelope modulation”更强，最准确的描述是：

`4D` 出现了非常强的长短交替大周期，也就是 `period-2-like relaxation pattern`。

### 第三层：这个 period-2 不是光滑曲面，而是碎裂参数骨架

后续局部放大显示：

- 在 `(theta0, E)` 平面里，强 `period-2` 区域先表现为一块 tongue-like wedge
- 继续放大后发现，这个 wedge 的脆弱角不是平滑边界，而是碎裂成很多窄窗
- 在 `ΔE = 1e-4` 级别上，强窗口会反复出现和消失
- 一些切片更像“离散微岛”
- 一些切片更像“细条纹，但被很多小洞打断”
- 沿 `eta_k` 和 `G` 方向再做微切片后，这种碎裂结构仍然存在

因此当前最准确的高维图景是：

强的大振幅 `period-2 relaxation` 并不是落在一张平滑 continuation surface 上，而是活在一个高维碎裂集合里。局部上更像由短 slab、细带、孔洞和微岛拼起来的参数骨架。

## 当前没有得到的结论

下面这些结论，当前证据还不支持：

- 还不能说已经找到 classical `MMO`
- 还不能说已经找到 fully developed `bursting`
- 还不能说已经观察到更高周期级联是主导机制
- 还不能说 `4D` 的碎裂结构已经被完整几何分类

到目前为止，更稳的说法仍然是：

`4D` 里最强现象是一个高维碎裂的 large-amplitude period-2 relaxation scaffold，而不是 bursting attractor。

## 最可靠的研究叙述

如果要把当前阶段写成论文前的内部结论，最稳的版本是：

1. `2D` 降维系统已经确认存在 `fold`、`subcritical Hopf`、迟滞和大振幅周期轨 backbone。
2. `3D` 说明这些结构在加入动态 `h` 后没有消失，但其主要动力学仍是单一 relaxation cycle，而不是 `MMO`。
3. `4D` 中出现了比 `3D` 更强的调制行为，并最终解析为一个高维碎裂的 `period-2-like relaxation` 参数骨架。
4. 这说明模型确实能产生高复杂度慢快结构，但当前最强证据指向的是复杂的倍周期调制，而不是 classical bursting。

## 建议的下一步

如果继续研究，最值得做的不是再大范围盲扫，而是以下三类工作：

1. 做正式图文整理
   把 `2D -> 3D -> 4D` 的证据链写成一版更像论文导言/结果节的说明。

2. 做局部几何解释
   围绕 `4D period-2` 微窗区，尝试找更稳定的 reduced description，解释为什么会出现“碎裂 slab + micro-islands”。

3. 谨慎继续找 bursting
   如果还想继续找真正 `bursting`，应该有针对性地改几何参数或慢变量耦合，而不是在现有窗口附近机械细扫。

## 结果文件导航

推荐优先看这些文件：

- 原始总摘要：
  [double_frontier_summary.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_summary.txt)
- 本阶段整理版：
  [double_frontier_research_note_2026-03-19.md](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_research_note_2026-03-19.md)
- `2D`：
  [double_frontier_hopf_neighborhood.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_hopf_neighborhood.png)
  [double_frontier_branch.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_branch.png)
- `3D`：
  [double_frontier_3d_hopf_zoom.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_3d_hopf_zoom.png)
  [double_frontier_3d_regime.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_3d_regime.png)
- `4D`：
  [double_frontier_4d_modulated_candidate.png](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_modulated_candidate.png)
  [double_frontier_4d_alternation.png](/mnt/d/double-frontier-rollout-dynamics/results/double-frontier_4d_alternation.png)
  [double_frontier_4d_period2_check.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_period2_check.txt)
  [double_frontier_4d_period2_boundary.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_period2_boundary.txt)
  [double_frontier_4d_period2_micro_windows.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_period2_micro_windows.txt)
  [double_frontier_4d_theta_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_theta_micro_slices.txt)
  [double_frontier_4d_eta_k_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_eta_k_micro_slices.txt)
  [double_frontier_4d_G_micro_slices.txt](/mnt/d/double-frontier-rollout-dynamics/results/double_frontier_4d_G_micro_slices.txt)
