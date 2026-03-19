# 4D Poincare / Return-Map Evidence

- Mode: `full`
- Section: `d = 0.9`, upward crossings only
- Integration window: `tmax = 36000`, `transient = 18000`, `saveat = 0.25`

| point | n_crossings | lag1(T) | lag2(T) | odd mean T | even mean T | odd CV | even CV | R1 | R2 | R2/R1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4D ridge_peak | 1034 | -0.998741 | 0.997453 | 5.093169 | 29.731263 | 0.037317 | 0.021137 | 0.033891 | 0.043520 | 1.284142 |
| 4D ridge_inner | 1062 | -0.999955 | 0.999894 | 5.044425 | 28.851719 | 0.016263 | 0.003898 | 0.039959 | 0.049657 | 1.242701 |
| 4D ridge_end | 1024 | -0.999670 | 0.999338 | 5.053684 | 30.079775 | 0.020569 | 0.010685 | 0.037328 | 0.047594 | 1.274999 |
| 4D boundary | 1006 | -0.999618 | 0.999113 | 15.081641 | 20.710030 | 0.006976 | 0.003665 | 0.036012 | 0.028503 | 0.791489 |

## Interpretation

- `4D ridge_peak`: current evidence weakens under the stricter Poincare diagnostic. R2/R1 = 1.284142.
- `4D ridge_inner`: current evidence weakens under the stricter Poincare diagnostic. R2/R1 = 1.242701.
- `4D ridge_end`: current evidence weakens under the stricter Poincare diagnostic. R2/R1 = 1.274999.
- `4D boundary`: current evidence weakens under the stricter Poincare diagnostic. R2/R1 = 0.791489.

- `ridge_inner` still has the tightest scalar alternation in `T_n` among the ridge-interior points: R2/R1 = 1.242701, lag1 = -0.999955, lag2 = 0.999894. The stricter state-space metric therefore weakens, rather than strengthens, the 2-cycle claim.
- The smallest observed R2/R1 is at `4D boundary`, with R2/R1 = 0.791489.
- `boundary` has R2/R1 = 0.791489 and odd/even period means [15.081641, 20.710030]. This indicates a weakened or changed mechanism relative to the ridge interior.

- If the odd/even clouds remain tight and R2 << R1, that supports a stable 2-cycle scaffold rather than diffuse modulation.
