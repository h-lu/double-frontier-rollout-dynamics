# Next-Phase Experiment Report

日期：2026-03-19

本文件是对 `results/next_phase/` 下本轮新实验的统一书面整理。目标不是记录脚本日志，而是用一份人类可读的文档，把实验目的、实验内容、方法、主要结果和结论写清楚。

## 1. 总体目标

本轮工作的核心目标不是继续盲目扩大 `MMO / bursting` 搜索，而是先把“当前结果到底说明了什么”坐实。围绕这个目标，本轮实验集中回答四个问题：

1. 现有 3D / 4D 诊断里使用的 `mid_frac = share of time with 0.2 < d < 0.9` 只是代理量。真实的双前沿诊断会不会改写对当前动力学的解释？
2. 4D 的 `period-2-like` 现象，能不能从长窗时序证据升级成更正式的返回映射 / Poincare 证据？
3. 目前看到的复杂现象，离默认经济学参数盒子有多远？是温和变形就出现，还是必须推到比较特殊的参数角落？
4. 3D 系统到底应不应该继续被当成 `MMO` 搜索器，还是更适合被表述为 “2D backbone 在动态 `h` 下的持续层”？

## 2. 本轮新增实验

本轮新增了四组实验和一组公共工具：

- 真实双前沿诊断：
  `scripts/next_phase_frontier_metrics.jl`
- 4D Poincare / 返回映射诊断：
  `scripts/next_phase_period2_poincare.jl`
- 默认参数盒到 modulation ridge 的单参数 homotopy：
  `scripts/next_phase_homotopy.jl`
- 3D persistence 扫描：
  `scripts/next_phase_3d_persistence.jl`
- 公共函数层：
  `scripts/next_phase_utils.jl`

所有新结果都写入 `results/next_phase/`，未覆盖原有 `results/` 文件。

## 3. 实验 A：真实双前沿诊断

### 3.1 实验目的

这组实验的目标是替换或补充原先的代理诊断，直接计算真实双前沿相关量，回答下面三个问题：

1. 4D ridge 的 `period-2-like` 行为，是否真的发生在 advisor cap 经常绑定的区域？
2. 4D ridge 上真实 middle band 是否仍有非平凡宽度？
3. 当前复杂动力学的主机制，更像是双前沿 cap 几何，还是更像 doer–burden–memory 耦合？

### 3.2 诊断内容

对每个仿真都直接计算：

- advisor frontier `a_t`
- raw doer target `r_t`
- true middle band `m_t = max(a_t - d_t, 0)`
- cap-binding indicator `1{r_t > a_t}`
- cap gap `max(r_t - a_t, 0)`
- `phi_t`

并汇总：

- `m_bar`
- `m_max`
- `capbind_frac`
- `capgap_mean`
- `phi_bar`
- `d_mean, q_mean, k_mean, h_mean`
- `m` 的 `5% / 50% / 95%` 分位数

### 3.3 参数点

共检查 8 个点：

1. 默认 `base_params()`
2. `2D refined Hopf`
3. `2D representative large-cycle point`
4. `3D representative oscillatory point`
5. `4D ridge_peak`
6. `4D ridge_inner`
7. `4D ridge_end`
8. `4D boundary`

其中 2D / 3D 的代表点都按现有脚本逻辑精确复原，完整参数集已经写在：
`results/next_phase/frontier_metrics_summary.md`

### 3.4 主要结果

4D ridge 四个代表点在 full 结果下的核心统计是：

- `ridge_peak`:
  `m_bar = 0.1152026461`,
  `capbind_frac = 0.0748739601`,
  `capgap_mean = 0.0000040860`
- `ridge_inner`:
  `m_bar = 0.1174075874`,
  `capbind_frac = 0.0768878210`,
  `capgap_mean = 0.0000041170`
- `ridge_end`:
  `m_bar = 0.1144731814`,
  `capbind_frac = 0.0743323009`,
  `capgap_mean = 0.0000040718`
- `boundary`:
  `m_bar = 0.1130995123`,
  `capbind_frac = 0.0734156470`,
  `capgap_mean = 0.0000040590`

### 3.5 结论

- 真实 middle band 在 4D ridge 上并不为零，`m_bar` 稳定在大约 `0.113` 到 `0.117`。
- 但 advisor cap 并没有“强绑定”：
  `capbind_frac` 只有约 `7.3%` 到 `7.7%`，
  `capgap_mean` 只有约 `4e-6`。
- 这说明当前 4D ridge 现象并不适合被表述成“系统大部分时间都处在强 advisor-cap 约束下的复杂动力学”。
- 更准确的说法是：
  当前复杂动力学更像 doer–burden–memory 耦合中的大振幅振荡，伴随一个非零 middle band，但不是由强 cap-binding 主导。

## 4. 实验 B：4D Poincare / 返回映射与 period-2 证据

### 4.1 实验目的

这组实验的目标是把原先基于长窗时序的 `period-2-like` 证据，升级成更正式的返回映射证据，并检验以下判断：

1. `ridge_inner` 是否真的呈现稳定的 2-cycle-like 返回结构？
2. `boundary` 点是同一分支的弱化，还是机制已经变化？
3. 如果使用更严格的 section-based 指标，现有结论是否仍然成立？

### 4.2 方法

使用完整 4D 模型，固定：

- `G = 1.0`
- `epsilon = 0.01`
- `eta_h = 0.002`
- `eta_k = 0.003`
- `beta_n = 1.5`

检查四个点：

- `ridge_peak: theta0 = -4.90, E = 1.08`
- `ridge_inner: theta0 = -4.90, E = 1.07`
- `ridge_end: theta0 = -4.86, E = 1.05`
- `boundary: theta0 = -4.84, E = 1.04`

积分设置：

- `tmax = 36000`
- `transient = 18000`
- `saveat = 0.25`

Poincare 截面定义为：

- `d = 0.9`
- 仅保留向上穿越

每次穿越记录：

- `t_n`
- `q_n, k_n, h_n`
- 局部 `d_dot`

并构造：

- `T_n = t_{n+1} - t_n`
- `W_n`（高态宽度）
- `lag1, lag2`
- odd/even means
- odd/even CV
- `R1 = mean(||X_{n+1}-X_n||)`
- `R2 = mean(||X_{n+2}-X_n||)`
- `R2_over_R1`

其中 `X_n = (q_n, k_n, h_n)`。

### 4.3 主要结果

full 结果如下：

- `ridge_peak`:
  `lag1(T) = -0.998741`,
  `lag2(T) = 0.997453`,
  `R2/R1 = 1.284142`
- `ridge_inner`:
  `lag1(T) = -0.999955`,
  `lag2(T) = 0.999894`,
  odd/even means = `[5.044425, 28.851719]`,
  odd/even CVs = `[0.016263, 0.003898]`,
  `R2/R1 = 1.242701`
- `ridge_end`:
  `R2/R1 = 1.274999`
- `boundary`:
  odd/even means = `[15.081641, 20.710030]`,
  `R2/R1 = 0.791489`

### 4.4 结论

- 如果只看 `T_n` 的 alternation，ridge 上的 `period-2-like` 证据仍然很强。
  `ridge_inner` 尤其接近完美长短交替。
- 但如果按更严格的 full-state section 指标来判断，结论明显变弱。
  我们原本希望看到的是 `R2 << R1`，
  也就是 `X_{n+2}` 比 `X_{n+1}` 更接近 `X_n`。
  实际上 ridge interior 三个点的 `R2/R1` 都大于 `1`。
- 因此，当前最诚实的说法是：
  当前证据支持强的 scalar timing alternation，
  但对完整 `(q,k,h)` 返回映射上的干净稳定 2-cycle，证据在 stricter diagnostics 下减弱。
- `boundary` 的 split 更弱，更像同类结构的弱化或局部机制变化，而不是更干净的 2-cycle。

## 5. 实验 C：从默认参数盒子到 modulation ridge 的 homotopy

### 5.1 实验目的

这组实验的目标是回答：
“当前这些复杂现象离默认经济学参数盒子到底有多远？”

也就是区分两种可能：

1. 复杂现象是默认盒子附近的温和变形；
2. 复杂现象只出现在特殊的多参数角落。

### 5.2 方法

只做五条单参数 homotopy，其余参数固定在默认 `base_params()`：

1. `beta_n: 0.7 -> 1.5`
2. `lambda_D: 12.0 -> 20.0`
3. `theta0: -2.0 -> -5.0`
4. `eta_h: 0.01 -> 0.002`
5. `eta_k: 0.03 -> 0.003`

full 模式每条路径用 `61` 个点。

每个点计算：

- `theoretical_fold_viable`
- 2D 低 / 高初值便宜代理下的 bistability 信号
- 4D sustained oscillation 代理
- true `m_bar`
- `capbind_frac`
- `period2_strength / R2_over_R1`
- `d_amp, q_amp, k_amp, h_amp`

### 5.3 主要结果

- `beta_n` 路径：
  没有触发 sustained 4D oscillation，也没有 2D bistability 的明显切换。
- `lambda_D` 路径：
  同样基本平稳，没有把系统推到 ridge-like 状态。
- `eta_h` 和 `eta_k` 路径：
  单独变化几乎不产生显著新结构。
- `theta0` 路径是唯一明显的例外：
  4D sustained oscillation 在 full 网格下大约从 `theta0 ≈ -3.95` 开始出现；
  到 `theta0 = -5.0` 时，`d_amp ≈ 0.966`。

但必须强调：

- 这条 `theta0` 路径虽然能点燃大振幅 4D 振荡，
  却没有产生 ridge 式的 cap-binding 或 clean 2-cycle 证据。
- 在这条路径上：
  `capbind_frac = 0`,
  `capgap_mean ≈ 0`,
  `R2/R1` 也不支持 clean scaffold。

### 5.4 结论

- default box 到 ridge 的距离不是“所有参数都很远”，但也绝不是“随便推一条参数就能到”。
- 最关键的单参数杠杆是 `theta0`。
- 但 modulation ridge 仍然是一个更特殊的多参数组织结果，不是默认盒子附近广泛存在的通用现象。

## 6. 实验 D：3D persistence，而不是继续 MMO 搜索

### 6.1 实验目的

这组实验的目标是重新定位 3D 模型：

- 不再把它当成主要的 `MMO` 搜索器；
- 而是检验它是不是更像 “2D backbone 在动态 `h` 下的持续层”。

### 6.2 方法

研究固定 `k = kbar = 0.3` 的 3D 模型，扫描：

- `E in [1.0, 1.8]`，full 模式 `31` 个点
- `eta_h in [0.001, 0.03]`，full 模式 `21` 个点

每个点输出：

- `d_amp`
- `q_amp`
- `h_amp`
- `return_time_mean`
- `return_time_cv`
- `crossing_ratio`
- true `m_bar`
- `capbind_frac`
- `phi_bar`
- `mid_frac_proxy`

### 6.3 主要结果

full 结果中的几个关键点：

- 最大 `d_amp`：
  `E = 1.000000`, `eta_h = 0.002450`,
  `d_amp = 0.991336`
- 最大 true `m_bar`：
  `E = 1.026667`, `eta_h = 0.001000`,
  `m_bar = 0.839927`
- 最大 `capbind_frac`：
  `E = 1.480000`, `eta_h = 0.030000`,
  `capbind_frac = 0.057064`
- 最大 `return_time_cv`：
  `E = 1.613333`, `eta_h = 0.002450`,
  `return_time_cv = 0.066699`

full 扫描下确实出现了几个局部可疑区，例如：

- `E = 1.613333`, `eta_h = 0.002450`, `return_time_cv = 0.066699`
- `E = 1.533333`, `eta_h = 0.002450`, `return_time_cv = 0.066344`
- `E = 1.506667`, `eta_h = 0.002450`, `return_time_cv = 0.041534`

但这些点同时具有：

- `crossing_ratio = 1`
- `mid_frac_proxy` 仍极小
- 没有形成强 mixed small/large 事件证据

### 6.4 结论

- 3D 系统整体上更像 persistent large-cycle layer，而不是 convincing MMO layer。
- true `m_bar` 在 3D 中可以很大，
  这再次说明旧 proxy 不能替代真实 middle-band 诊断。
- full 扫描确实给出了一小片需要单独保留的可疑区域，
  但这些区域还不足以支持“3D 已经找到强 MMO”的判断。
- 因此，3D 最合适的定位仍然是：
  `2D backbone` 在动态 `h` 下的持续存在，而不是 MMO 搜索器。

## 7. 综合结论

本轮实验最重要的结论可以压缩成五句：

1. 真实双前沿诊断改写了旧解释：
   4D ridge 上 true middle band 非零，但 advisor cap 并不强绑定。
2. 4D ridge 的 `period-2-like` 证据在 timing sequence 上很强，但在严格的 full-state Poincare 指标下减弱。
3. modulation ridge 不是默认参数盒子附近的普遍现象；
   在单参数 homotopy 里，只有 `theta0` 显示出明显的 4D 大振幅振荡开关。
4. 3D 更适合被当成 `2D backbone` 的 persistence layer，
   而不是 MMO 搜索器。
5. 因而，当前最合适的研究收束不是继续追 `MMO / bursting`，
   而是先把现有 4D `period-2-like` 结构做扎实。

## 8. 下一步建议

本轮实验支持的明确建议是：

**C. 暂时不要追 MMO / bursting，先把现有 period-2 结构做扎实。**

原因是：

- `2D backbone` 已经足够稳；
- `3D persistence` 也已经基本成立；
- `4D` 的 timing alternation 是真的；
- 但 full-state 2-cycle 几何证据还不够强；
- 同时，真实双前沿诊断又表明当前 4D 现象并不是一个强 advisor-cap-binding 机制。

所以最合理的顺序是：

1. 先把现有 4D 返回映射和 period-2 结构解释扎实；
2. 再决定是否需要为“advisor cap 真正进入内部机制”而改模型；
3. 在那之前，不宜继续把主要资源投入到 `MMO / bursting` 搜索上。

## 9. 相关结果文件

建议配合阅读：

- `results/next_phase/frontier_metrics_summary.md`
- `results/next_phase/poincare_summary.md`
- `results/next_phase/homotopy_summary.md`
- `results/next_phase/threeD_persistence_summary.md`
- `results/next_phase/next_phase_summary.md`
