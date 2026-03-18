# Double-Frontier Rollout Dynamics

This repository contains a Julia workflow for studying a slow-fast rollout model
with advisor and doer frontiers, supervision backlog, adoption capital, and
judgment capital. The code is organized around the research sequence we agreed
on:

1. start with a smooth continuation-friendly model;
2. reduce to `2D` for fold, Hopf, hysteresis, and canard-style screening;
3. return to the `4D` system only after the reduced geometry is understood.

The project intentionally uses Julia only. There is no Python or MATLAB path in
the workflow.

## Model Layers

The main script is a smooth `4D` model with states:

- `d`: doer share
- `q`: supervision / integration backlog
- `k`: adoption capital
- `h`: judgment capital

From that full model, the code also exposes:

- a `3D` reduction that fixes `k`
- a `2D` reduction that fixes `k` and `h`

The `2D` system is the current workhorse for numerical screening because it is
the right place to determine whether the model actually supports:

- fold-like geometry
- Hopf onset
- hysteresis / bistability
- large-amplitude relaxation oscillations
- canard-like sharp amplitude jumps

## Repository Layout

- `Project.toml` and `Manifest.toml`: Julia environment
- `scripts/install_status.jl`: package sanity check
- `scripts/double_frontier_rollout.jl`: full model definition, `2D` continuation,
  hysteresis sweeps, and representative `4D` simulations
- `scripts/double_frontier_canard_search.jl`: targeted `2D` search around
  relaxation / canard-like jump windows
- `results/`: generated plots and text summaries

## Current Outputs

The tracked results currently focus on the double-frontier model:

- `double_frontier_fold_screen.png`: analytic pre-screen based on
  `beta_n * lambda_D > 4`
- `double_frontier_branch.png`: reduced-branch continuation in rollout pressure
  `E`
- `double_frontier_hysteresis.png`: low/high initial condition sweeps
- `double_frontier_phase.png`: representative `2D` phase portrait and time series
- `double_frontier_full_timeseries.png`: representative `4D` trajectories
- `double_frontier_summary.txt`: summary of the reduced-model screening results

## Run

Use the local Julia installation:

```bash
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/install_status.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_rollout.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_canard_search.jl
```

## Numerical Stack

- `OrdinaryDiffEq.jl` is used for stiff integration, currently with `Rodas5P()`
- `BifurcationKit.jl` is used for equilibrium continuation and bifurcation
  detection
- `Plots.jl` is used for all figures

## Research Status

At the current checkpoint, the reduced model already shows strong path
dependence and a fold/Hopf candidate structure in selected parameter regions.
That makes the next numerical question very specific:

- is the oscillation onset merely a switch into a large stable cycle?
- or is there a genuine narrow canard-style transition from tiny oscillations to
  relaxation oscillations?

The repository is set up so that this question can be answered by iterating on
the `2D` scripts first, before investing in more expensive `3D` or `4D`
continuation work.
