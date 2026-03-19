# Next-Phase Summary

This summary reflects the `--full` reruns of the four next-phase scripts.

## 1. What was added

- New scripts:
  `scripts/next_phase_utils.jl`,
  `scripts/next_phase_frontier_metrics.jl`,
  `scripts/next_phase_period2_poincare.jl`,
  `scripts/next_phase_homotopy.jl`,
  `scripts/next_phase_3d_persistence.jl`.
- New shared functions in `scripts/next_phase_utils.jl`:
  `advisor_frontier_series`,
  `raw_doer_target`,
  `raw_doer_target_series`,
  `frontier_metric_series`,
  `frontier_metric_summary`,
  `representative_parameter_sets`,
  `parameter_point`,
  `simulate_point`,
  `threshold_crossing_data`,
  `poincare_metrics`,
  `period2_strength_from_crossings`,
  `cheap_2d_bistability_proxy`,
  `cheap_4d_summary`,
  plus CSV / parameter-format helpers.
- New result bundles in `results/next_phase/`:
  `frontier_metrics_table.csv`,
  `frontier_metrics_summary.md`,
  `frontier_timeseries_*.png`,
  `frontier_hist_*.png`,
  `poincare_*.csv`,
  `poincare_summary.md`,
  `poincare_*_ridge_inner.png`,
  `poincare_qh_odd_even_*.png`,
  `homotopy_*.csv`,
  `homotopy_regime_map.png`,
  `homotopy_summary.md`,
  `threeD_persistence.csv`,
  `threeD_persistence_map.png`,
  `threeD_persistence_summary.md`.
- Parameter reconstruction is now explicit rather than implicit:
  the full NamedTuple-equivalent parameter sets for the 2D refined Hopf point,
  the 2D representative large-cycle point,
  the 3D representative oscillatory point,
  and the 4D ridge points are written out in
  `results/next_phase/frontier_metrics_summary.md`.

## 2. True double-frontier diagnostics

- The old proxy `mid_frac = share of time with 0.2 < d < 0.9` was not measuring the true middle band.
  The new diagnostics separate:
  advisor frontier `a_t`,
  raw doer target `r_t`,
  true middle band `m_t = max(a_t - d_t, 0)`,
  cap-binding indicator `1{r_t > a_t}`,
  cap gap `max(r_t - a_t, 0)`,
  and `phi_t`.
- For the four 4D ridge points, the true middle band is not negligible:
  `m_bar` is about `0.1131` to `0.1174`,
  with ridge-average `m_bar ≈ 0.115`.
  So the ridge does not collapse to a zero-width middle band.
- But the same four ridge points are not strongly cap-binding in the sense needed for a clean advisor-cap story:
  `capbind_frac` is only about `0.0734` to `0.0769`,
  and `capgap_mean` is only about `4.06e-6` to `4.12e-6`.
  In other words, the raw doer target is only very slightly above the advisor frontier, and only intermittently.
- This materially changes the interpretation.
  The 4D ridge dynamics do not look like a regime that lives mainly in a robust “advisor cap frequently binds hard” region.
  They look more like large doer–burden–memory oscillations that happen while a nonzero middle band is present, with the cap only lightly active.
- The same distinction matters in 3D.
  In the 3D persistence scan, the legacy proxy stays tiny while true `m_bar` can be large.
  So “little time with `0.2 < d < 0.9`” does not imply “little true middle band.”

Working answer for section A:

- 4D ridge period-2-like behavior does **not** appear to live in a strongly cap-binding zone.
- True `m_t` remains nontrivial on the ridge.
- The present complex dynamics look driven more by doer–burden–memory coupling than by a hard advisor-cap geometry.

## 3. Period-2 evidence

- The stricter Poincare section was defined exactly as requested:
  `d = 0.9`,
  upward crossings only,
  `tmax = 36000`,
  `transient = 18000`,
  `saveat = 0.25`.
- On scalar timing diagnostics, the ridge points still look extremely alternating.
  For example at `ridge_inner`:
  `lag1(T) = -0.999955`,
  `lag2(T) = 0.999894`,
  odd/even means are `[5.044425, 28.851719]`,
  and odd/even CVs are `[0.016263, 0.003898]`.
  So the return-time sequence itself is very close to perfect alternation.
- But the stricter state-space test weakens the period-2 claim.
  Using `X_n = (q_n, k_n, h_n)` on the section, the requested metric
  `R2_over_R1 = mean(||X_{n+2}-X_n||) / mean(||X_{n+1}-X_n||)`
  is:
  `1.284` at `ridge_peak`,
  `1.243` at `ridge_inner`,
  `1.275` at `ridge_end`,
  and `0.791` at `boundary`.
- So the result is mixed:
  scalar alternation is very strong,
  but the full-section point cloud does not satisfy the hoped-for `R2 << R1` pattern.
  Under the stricter diagnostic, current evidence for a clean stable 2-cycle in the full `(q,k,h)` return map therefore weakens.
- `boundary` is not “the same branch but more clearly 2-cyclic.”
  Its `R2/R1` is lower than the ridge-interior points, but its odd/even split is much smaller:
  odd/even period means are `[15.081641, 20.710030]`.
  That is more consistent with a weaker or partially changed mechanism than with a stronger 2-cycle.

Working answer for section B:

- `ridge_inner` is still the strongest scalar alternation point, but it is **not** a clean state-space 2-cycle by the requested `R2/R1` test.
- `boundary` looks weaker, not stronger.
- The honest statement is:
  current evidence for “period-2-like” survives in timing space,
  but **current evidence weakens under stricter Poincare diagnostics** when judged by full-section contraction.

## 4. Homotopy conclusion

- The five single-parameter homotopies were run from the default `base_params()` box with all other parameters fixed at default values.
- `beta_n`, `lambda_D`, `eta_h`, and `eta_k` alone do essentially nothing dramatic.
  In quick mode they do not create sustained 4D oscillation, do not create a convincing 2D bistability signal, and do not create ridge-like period-2 structure.
- `theta0` is the one clear exception.
  Along the default-box homotopy `theta0: -2 -> -5`, sustained large-amplitude 4D oscillation turns on around `theta0 ≈ -3.95`.
  By `theta0 = -5.0`, `d_amp ≈ 0.966`.
- But that theta-only deformation still does **not** look like the modulation ridge story.
  Along that path:
  `capbind_frac = 0`,
  cap-gap is numerically zero,
  and `R2/R1` does not indicate a clean 2-cycle scaffold.
- So the homotopy answer is not “the ridge is right next to the default box,” but also not “everything is maximally special.”
  A single parameter, `theta0`, can already turn on large 4D oscillation.
  What remains far away is the specific ridge package:
  negative `theta0` together with the small `epsilon`, slow `eta_h`, slow `eta_k`, and altered high-dimensional geometry near the ridge.

Working answer for section C:

- The modulation ridge is **not** reached by gentle one-parameter deformations of `beta_n`, `lambda_D`, `eta_h`, or `eta_k`.
- `theta0` is the dominant lever for turning oscillations on.
- But the ridge-like alternation / modulation structure still looks like a more special multivariate corner than a generic nearby extension of the default box.

## 5. Recommended next step

Recommendation: **C. 暂时不要追 MMO / bursting，先把现有 period-2 结构做扎实。**

Reason:

- Option A is premature.
  The 4D story is not yet strong enough to headline as a clean “period-2 scaffold,” because the stricter Poincare test weakens the full-section 2-cycle claim.
- Option B is conceptually justified if the future goal is specifically a paper about advisor-cap-driven complex dynamics.
  The new diagnostics show that the current ridge is not strongly cap-binding, so the present 4D phenomenon is not a clean advisor-cap mechanism.
  But changing the model immediately would blur what the current model already does.
- Option C is the right immediate move because it matches the evidence actually in hand:
  2D backbone is solid,
  3D persistence is solid overall,
  4D timing alternation is real,
  but the full-section geometric interpretation is still incomplete.

- The full 3D persistence scan did surface a few localized suspicious patches with return-time CV up to about `0.0667`, mainly near `E ≈ 1.51–1.61` and small `eta_h ≈ 0.00245`.
  Those patches should be listed as caveats, not promoted to MMO evidence:
  crossing ratio stays at `1`,
  true middle-band width can be moderate,
  and the old proxy remains tiny.

Concretely, the next pass should be:

- keep the current model fixed;
- treat the current 4D claim as
  “strong scalar alternation with incomplete full-section 2-cycle confirmation”;
- tighten the 4D return-map story before making a paper-level claim;
- only after that decide whether a model revision is needed to make advisor-cap geometry genuinely internal.

If the eventual research objective is specifically
“double-frontier cap geometry generates high-dimensional complex dynamics,”
then option B is likely the next step **after** this consolidation pass.
