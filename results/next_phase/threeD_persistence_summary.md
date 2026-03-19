# 3D Persistence Scan

- Mode: `full`
- Fixed 3D k = 0.300000
- Base parameter family: `base_params(beta_n=2, lambda_D=20, theta0=-5)` with only `E` and `eta_h` varied
- Legacy proxy retained as `mid_frac_proxy = share of time with 0.2 < d < 0.9`; it is not the true middle band.

## Highlights

- Max d_amp occurs at E = 1.000000, eta_h = 0.002450 with d_amp = 0.991336.
- Max true m_bar occurs at E = 1.026667, eta_h = 0.001000 with m_bar = 0.839927.
- Max capbind_frac occurs at E = 1.480000, eta_h = 0.030000 with capbind_frac = 0.057064.
- Max return_time_cv occurs at E = 1.613333, eta_h = 0.002450 with CV = 0.066699.

## Interpretation

- The main 3D signature should be read against `d_amp`, `return_time_cv`, `true m_bar`, and `capbind_frac`, not against the old mid-band proxy alone.
- Some localized return-time irregularity appears, but it is not enough on its own to establish MMO.
- The true middle band becomes nontrivial in parts of the scan, but that does not automatically imply MMO.
- Local suspicious regions worth keeping separate from a hard MMO claim:
  - E = 1.613333, eta_h = 0.002450, return_time_cv = 0.066699, crossing_ratio = 1.000000, true m_bar = 0.074097, mid_frac_proxy = 0.000311
  - E = 1.533333, eta_h = 0.002450, return_time_cv = 0.066344, crossing_ratio = 1.000000, true m_bar = 0.158693, mid_frac_proxy = 0.000356
  - E = 1.506667, eta_h = 0.002450, return_time_cv = 0.041534, crossing_ratio = 1.000000, true m_bar = 0.200475, mid_frac_proxy = 0.000533
  - E = 1.613333, eta_h = 0.003900, return_time_cv = 0.039963, crossing_ratio = 1.000000, true m_bar = 0.076157, mid_frac_proxy = 0.000711
  - E = 1.400000, eta_h = 0.001000, return_time_cv = 0.038441, crossing_ratio = 1.000000, true m_bar = 0.264037, mid_frac_proxy = 0.000044
