# Research Background

This note organizes the current research plan behind the double-frontier rollout
model. It is intended as background for the code in this repository rather than
as a polished paper draft.

## Objective

The main goal is to determine which dynamical mechanisms are genuinely present
in the model before committing to a paper narrative. The key question is not
"can we make a complicated `4D` system look interesting?" but:

- does the reduced model have fold geometry?
- is there a Hopf bifurcation?
- are there stable limit cycles?
- is there a narrow canard-style transition?
- do higher-dimensional reductions support MMO or bursting?

That is why the workflow starts from reduced models and only returns to the full
system once the geometry is already understood.

## Julia-Only Numerical Strategy

The original planning discussion referenced SciPy, MatCont, and AUTO as the
standard numerical toolkit for stiff slow-fast systems and continuation. In this
repository we implement the same research logic entirely in Julia:

- `OrdinaryDiffEq.jl` replaces the stiff-IVP role
- `BifurcationKit.jl` replaces the equilibrium / periodic continuation role
- `Plots.jl` handles diagnostic and presentation figures

The practical numerical lessons still carry over:

- small `epsilon` makes the IVP stiff
- smooth approximations are preferable to hard `min` or positive-part kinks
- coarse parameter sweeps can easily miss narrow canard windows
- continuation and direct simulation should be used together rather than as
  substitutes

## Main Model

The baseline object is a smooth `4D` system with:

- `d`: doer share, fast variable
- `q`: supervision / integration backlog, intermediate variable
- `k`: adoption capital, slow variable
- `h`: judgment capital, slow variable

The model preserves the core economic ideas:

- an advisor frontier
- a doer frontier
- a middle band between them
- supervision burden from rollout
- adoption capital accumulated through repeated use
- judgment capital accumulated or eroded depending on oversight quality

### Smooth Components

To avoid hard nonsmooth corners:

- `smin_sigma(x, y)` approximates `min(x, y)`
- `softplus_sigma(r)` approximates `max(r, 0)`

These enter the code in
[double_frontier_rollout.jl](/Users/wangxq/Documents/AUTO/scripts/double_frontier_rollout.jl)
so that continuation and stiff integration remain well behaved.

### Frontiers and Coverage

The model contains three central nonlinear pieces:

1. advisor frontier `a(h)`
2. supervision coverage `phi(h, q)`
3. doer target `d_star(d, q, k, h)`

The economic interpretation is:

- higher `z`, `G`, and `h` expand the advisor frontier
- higher `G` and `h` improve supervision coverage
- higher `E`, `k`, and short-run feedback `beta_n` raise doer pressure
- larger backlog `q` pushes back against further rollout

## Reduced Models

Three nested systems guide the project.

### 2D model

Fix `k = kbar` and `h = hbar`, then study `(d, q)`.

This is the first screening layer for:

- S-shaped critical manifolds
- fold-like structure
- Hopf candidates
- hysteresis
- canard-like sharp amplitude jumps
- relaxation oscillations

### 3D model

Fix `k = kbar`, then study `(d, q, h)`.

This is the first layer where folded singularities and MMO can appear in a way
that matches the standard slow-fast geometry literature.

### 4D model

Use the full system only after the reduced structure is already mapped.

This layer is mainly for:

- rollout waves
- long-memory effects from `(k, h)`
- bursting-style slow passage between quiet and oscillatory regimes

## Parameter Box

The project uses a numerical exploration box rather than an empirical
calibration.

### Baseline values

Suggested baseline values in the planning stage were:

- `lambda_A = 8`
- `alpha0 = -1`
- `alpha_z = alpha_h = alpha_G = 1`
- `theta0 = -2`
- `theta_z = 1`
- `theta_h = theta_G = 0.8`
- `theta_E = 1.2`
- `theta_k = theta_phi = 1`
- `nu0 = 0.3`
- `delta_k = 0.2`
- `mu0 = 0.05`
- `mu1 = 1`
- `omega0 = -1`
- `omega_G = 1.2`
- `omega_h = 1`
- `omega_q = 1.5`

### First-pass scan parameters

The most important early scan parameters are:

- `beta_n`
- `lambda_D`
- `E`
- `G`
- `z`
- `gamma`
- `chi_R`
- `epsilon`
- `eta_h`

Additional discrete choices of interest:

- `rho in {1.5, 2, 3}`
- `chi_I in {0, 0.5, 1}`
- `sigma in {20, 50, 100, 200}`

### Default initial conditions

The baseline initial state is:

- `d(0) = 0.05`
- `q(0) = 0`
- `k(0) = 0.1`
- `h(0) = 0.6`

High-adoption and low-adoption alternatives are also essential for testing
bistability.

## Experimental Roadmap

The numerical plan is staged deliberately.

### Experiment 0: numerical stability check

Purpose:

- make sure smoothing does not create fake dynamics
- make sure tolerance choices are not driving conclusions
- make sure the stiff solver is behaving credibly

Outputs:

- time series under multiple tolerances
- phase portraits under multiple `sigma`
- solver diagnostics

### Experiment 1: baseline region without fold

Purpose:

- confirm what the model looks like when `beta_n` is very small
- separate ordinary slow-manifold adjustment from genuine nonlinear geometry

Outputs:

- `2D` phase portraits
- nullclines
- equilibrium continuation in `E`

### Experiment 2: analytic fold pre-screen

Ignoring the advisor cap locally, the reduced equilibrium relation suggests the
necessary condition

`beta_n * lambda_D > 4`

for fold potential in the logistic case.

Purpose:

- shrink the parameter search space before brute-force scanning

Outputs:

- theoretical fold-viability map
- numerical three-equilibrium regions
- representative critical manifolds

### Experiment 3: advisor cap clipping

Purpose:

- test whether the advisor frontier truncates the right branch of the reduced
  S-shape
- show that the double-frontier structure changes topology, not just accounting

Outputs:

- cap-versus-fold comparisons
- `(G, h)` region maps
- clipped versus unclipped phase portraits

### Experiment 4: bistability and hysteresis

Purpose:

- establish the most robust nonlinear phenomenon first

Outputs:

- continuation branches in `E`
- forward/backward sweeps
- low/high initial condition comparisons

### Experiment 5: Hopf and periodic orbits

Purpose:

- determine whether the reduced model supports genuine oscillations

Outputs:

- Hopf detection
- periodic branch continuation
- extrema and period plots

### Experiment 6: narrow canard window

Purpose:

- test whether the model shows a sharp transition from very small oscillations
  to large relaxation cycles over a narrow interval

Outputs:

- zoomed amplitude plots
- zoomed period plots
- representative small / canard-like / large-cycle trajectories

### Experiment 7: relaxation oscillations

Purpose:

- characterize the geometry after the jump window

Outputs:

- large-cycle phase portraits
- time series with slow drift and fast jumps

### Experiment 8: folded singularities and MMO in 3D

Purpose:

- determine whether the `(d, q, h)` reduction contains folded-node-type
  organization

Outputs:

- fold curves
- folded singularity checks
- long time series
- SAO counting plots

### Experiment 9: bursting in 4D

Purpose:

- understand how slow drift in `(k, h)` moves the system through fast-subsystem
  regimes

Outputs:

- long-run `4D` trajectories
- `(k, h)` slow-drift plots
- fast-subsystem skeleton overlays

### Experiment 10: middle-band irreducibility

Purpose:

- show that two rollout paths with similar final doer share can generate
  different `k` and `h` histories

Outputs:

- path-comparison time series
- `m-h` and `m-k` plots
- matched-final-`d` comparisons

## Figure Logic

Each figure family answers a different question.

- Time series: jumps, drift, overshoot, bursting
- `(d, q)` phase portraits: nullclines, folds, cycles, relaxation loops
- Critical-manifold plots: S-shape and cap clipping
- Single-parameter continuation: hysteresis, Hopf, periodic envelopes
- Two-parameter continuation: fold and Hopf boundaries
- Zoomed amplitude plots: canard-window evidence
- Fold-curve and folded-singularity plots: MMO geometry
- Peak-return or SAO-count plots: MMO patterning
- Fast-skeleton plus full-trajectory overlays: bursting organization

## Stop-Loss Rules

These rules matter as much as the experiments themselves.

1. If no fold appears in a reasonable parameter box, stop chasing canards or
   MMO and pivot to slow-manifold / path-dependence / tipping.
2. If folds appear but no reliable Hopf or periodic branch appears, emphasize
   hysteresis and multistability instead of forcing a canard story.
3. If `2D` canard-like behavior appears but `3D` has no folded singularity
   structure, frame the math contribution as a planar canard / relaxation paper.
4. Only talk about MMO or bursting in titles or abstracts if the higher
   dimensional evidence is unambiguous.

## Recommended Paper Routes

Depending on what the numerics ultimately show:

- only slow manifolds, hysteresis, or tipping:
  economic theory paper
- fold plus canard-style onset plus relaxation oscillations:
  `2D` applied-math paper
- folded node plus MMO:
  strongest applied-math route
- clear bursting:
  long-form paper or second paper

## Current Repository Position

The current codebase is still at the early part of this roadmap.

- the `2D` reduced system is implemented and continuation-ready
- the `4D` system is implemented for representative simulations
- strong path dependence is already visible in selected regions
- canard / relaxation search is still exploratory rather than definitive

That means the present repository should be read as a serious numerical
exploration scaffold, not yet as a final theorem-or-paper package.

## Source Notes

The planning discussion referenced the following sources as background:

- [SciPy `solve_ivp` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
- [Guckenheimer / MMO review](https://pi.math.cornell.edu/~gucken/PDF/mmo_review.pdf)
- [AUTO-07P overview](https://archive-dsweb.siam.org/Software/auto-07p.html)
- [Canard explosion reference note](https://www.researchgate.net/publication/247384460_Relaxation_Oscillation_and_Canard_Explosion?utm_source=chatgpt.com)

These are preserved here as planning references. The repository implementation
itself is Julia-based.
