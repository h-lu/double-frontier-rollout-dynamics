# True Double-Frontier Diagnostics

- Mode: `full`
- Script: `scripts/next_phase_frontier_metrics.jl`
- Output table: `results/next_phase/frontier_metrics_table.csv`

## Parameter Points

### default base_params()

- source: base_params() default
- dimension: `full`
- E = 0.2000000000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0300000000
- eta_h = 0.0100000000
- eta_k = 0.0300000000
- beta_n = 0.7000000000
- lambda_D = 12.0000000000
- theta0 = -2.0000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 2D refined Hopf

- source: scripts/double_frontier_rollout.jl: refined Hopf from base_params(beta_n=2, lambda_D=20, theta0=-5)
- dimension: `twoD`
- E = 0.9747525391
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0300000000
- eta_h = 0.0100000000
- eta_k = 0.0300000000
- beta_n = 2.0000000000
- lambda_D = 20.0000000000
- theta0 = -5.0000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 2D representative large-cycle point

- source: scripts/double_frontier_canard_search.jl: representative above-jump simulation at epsilon=0.02
- dimension: `twoD`
- E = 1.2161870000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0200000000
- eta_h = 0.0100000000
- eta_k = 0.0300000000
- beta_n = 2.0000000000
- lambda_D = 20.0000000000
- theta0 = -5.0000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 2.0000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.1000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 3D representative oscillatory point

- source: scripts/double_frontier_rollout.jl: rep3d_idx = argmax(max(amp_low, amp_high) + gap_d + 0.25*gap_h)
- dimension: `threeD`
- E = 1.1666666667
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0300000000
- eta_h = 0.0100000000
- eta_k = 0.0300000000
- beta_n = 2.0000000000
- lambda_D = 20.0000000000
- theta0 = -5.0000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 4D ridge_peak

- source: results/double_frontier_4d_period2_check.txt representative point
- dimension: `full`
- E = 1.0800000000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0100000000
- eta_h = 0.0020000000
- eta_k = 0.0030000000
- beta_n = 1.5000000000
- lambda_D = 20.0000000000
- theta0 = -4.9000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 4D ridge_inner

- source: results/double_frontier_4d_period2_check.txt representative point
- dimension: `full`
- E = 1.0700000000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0100000000
- eta_h = 0.0020000000
- eta_k = 0.0030000000
- beta_n = 1.5000000000
- lambda_D = 20.0000000000
- theta0 = -4.9000000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 4D ridge_end

- source: results/double_frontier_4d_period2_check.txt representative point
- dimension: `full`
- E = 1.0500000000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0100000000
- eta_h = 0.0020000000
- eta_k = 0.0030000000
- beta_n = 1.5000000000
- lambda_D = 20.0000000000
- theta0 = -4.8600000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

### 4D boundary

- source: results/double_frontier_4d_period2_check.txt representative point
- dimension: `full`
- E = 1.0400000000
- G = 1.0000000000
- z = 1.0000000000
- epsilon = 0.0100000000
- eta_h = 0.0020000000
- eta_k = 0.0030000000
- beta_n = 1.5000000000
- lambda_D = 20.0000000000
- theta0 = -4.8400000000
- lambda_A = 8.0000000000
- c_A = 0.0000000000
- alpha0 = -1.0000000000
- alpha_z = 1.0000000000
- alpha_h = 1.0000000000
- alpha_G = 1.0000000000
- c_D = 0.0000000000
- theta_z = 1.0000000000
- theta_h = 0.8000000000
- theta_G = 0.8000000000
- theta_E = 1.2000000000
- theta_k = 1.0000000000
- theta_phi = 1.0000000000
- gamma_q = 1.4000000000
- omega0 = -1.0000000000
- omega_G = 1.2000000000
- omega_h = 1.0000000000
- omega_q = 1.5000000000
- nu0 = 0.3000000000
- nu_G = 0.5000000000
- nu_h = 0.4000000000
- chi_I = 0.5000000000
- chi_R = 0.8000000000
- rho = 2.0000000000
- delta_k = 0.2000000000
- psi_M = 0.7000000000
- psi_D = 0.6000000000
- vartheta_M = 0.7000000000
- vartheta_D = 0.3500000000
- vartheta_H = 0.4500000000
- mu0 = 0.0500000000
- mu1 = 1.0000000000
- sigma = 80.0000000000
- hbar = 0.6500000000
- kbar = 0.3000000000

## Metric Table

| point | dim | m_bar | m_max | capbind_frac | capgap_mean | phi_bar | d_mean | q_mean | k_mean | h_mean | m_q05 | m_q50 | m_q95 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| default base_params() | full | 0.009240 | 0.009240 | 0.000000 | 0.000000 | 0.250739 | 0.990620 | 0.935760 | 0.794579 | 0.108965 | 0.009240 | 0.009240 | 0.009240 |
| 2D refined Hopf | twoD | 0.999966 | 0.999966 | 0.000000 | 0.000000 | 0.699280 | 0.000033 | 0.004087 | 0.300000 | 0.650000 | 0.999966 | 0.999966 | 0.999966 |
| 2D representative large-cycle point | twoD | 0.876296 | 0.999998 | 0.061288 | 0.000000 | 0.663511 | 0.123702 | 0.107923 | 0.300000 | 0.650000 | 0.008667 | 0.982807 | 0.999998 |
| 3D representative oscillatory point | threeD | 0.558290 | 0.999993 | 0.010622 | 0.000000 | 0.567648 | 0.441707 | 0.339175 | 0.300000 | 0.600160 | 0.008747 | 0.988864 | 0.999785 |
| 4D ridge_peak | full | 0.115203 | 0.999939 | 0.074874 | 0.000004 | 0.321269 | 0.884736 | 0.778372 | 1.234015 | 0.212579 | 0.008637 | 0.037288 | 0.999805 |
| 4D ridge_inner | full | 0.117408 | 0.999940 | 0.076888 | 0.000004 | 0.322587 | 0.882532 | 0.775844 | 1.242484 | 0.214731 | 0.008637 | 0.037287 | 0.999836 |
| 4D ridge_end | full | 0.114473 | 0.999939 | 0.074332 | 0.000004 | 0.320812 | 0.885465 | 0.779270 | 1.231207 | 0.211866 | 0.008637 | 0.037296 | 0.999795 |
| 4D boundary | full | 0.113100 | 0.999938 | 0.073416 | 0.000004 | 0.319946 | 0.886838 | 0.780928 | 1.225581 | 0.210440 | 0.008637 | 0.037320 | 0.999777 |

## Interpretation

- Across the four 4D ridge points, average cap-binding frequency is 0.074877 and average true middle-band width is 0.115046.
- The strongest cap-binding point is `ridge_inner` with capbind_frac = 0.076888; the weakest is `boundary` with capbind_frac = 0.073416.
- The largest true middle band is at `ridge_inner` with m_bar = 0.117408; the smallest is `boundary` with m_bar = 0.113100.
- The ridge-inner point has capbind_frac = 0.076888, capgap_mean = 0.000004, and m_bar = 0.117408.
- The boundary point has capbind_frac = 0.073416, capgap_mean = 0.000004, and m_bar = 0.113100.

- Provisional answer: the 4D ridge does not spend much time in a strongly cap-binding regime; the advisor cap is active only intermittently in these samples.
- The true middle band remains nontrivial across the ridge, so the cap geometry does not collapse to zero-width.
- Mean cap-gap is very small, suggesting the raw doer target rarely pushes far past the advisor frontier.
