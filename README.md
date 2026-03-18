# BifurcationKit Slow-Fast Demo

This workspace contains a Julia environment for studying slow-fast systems with
`BifurcationKit.jl` and `OrdinaryDiffEq.jl`.

The first reproducible example is a classical canard-explosion analysis for the
van der Pol slow-fast system:

```math
\dot{x} = y - x^3/3 + x,\qquad
\dot{y} = \varepsilon(a - x).
```

We continue equilibria in the parameter `a`, detect the Hopf bifurcation near
`a = 1`, switch to the periodic-orbit branch, and compare direct simulations on
both sides of the canard explosion.

The second example is a `3`-variable Hindmarsh-Rose neuron model with one slow
adaptation variable:

```math
\dot{x} = y - ax^3 + bx^2 - z + I,\qquad
\dot{y} = c - dx^2 - y,\qquad
\dot{z} = r(s(x-x_R)-z).
```

For this model we continue the equilibrium branch in the applied current `I`,
detect the Hopf onset of oscillations, and then use direct simulation to map
the transition from rest to bursting and tonic spiking.

The third example is a Julia-only scaffold for the double-frontier rollout
model discussed in the research plan. It keeps the full `4`-variable system
available, but follows the staged workflow numerically:

- use `OrdinaryDiffEq.jl` stiff solvers for direct simulation of the `2D`, `3D`,
  and `4D` systems;
- use `BifurcationKit.jl` continuation on the `2D` reduced model to screen for
  fold and Hopf structure in the pressure parameter `E`;
- use direct low/high-initial-condition sweeps to check hysteresis before moving
  to MMO or bursting work.

## Layout

- `scripts/install_status.jl`: quick environment sanity check
- `scripts/canard_vdp.jl`: continuation + simulation workflow
- `scripts/hindmarsh_rose_bursting.jl`: 3-variable bursting analysis
- `scripts/double_frontier_rollout.jl`: double-frontier rollout model, `2D`
  reduced continuation, hysteresis sweep, and `4D` long-run simulation
- `results/`: generated plots and data

## Run

Use the locally installed Julia runtime:

```bash
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/install_status.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/canard_vdp.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/hindmarsh_rose_bursting.jl
/Users/wangxq/.local/julia/julia-1.12.5/bin/julia --project=. scripts/double_frontier_rollout.jl
```

## Julia Tooling Notes

This workspace intentionally avoids Python and MATLAB.

- `OrdinaryDiffEq.jl` replaces the SciPy stiff-IVP workflow; the current
  scripts use `Rodas5P()` for small-`epsilon` slow-fast integration.
- `BifurcationKit.jl` replaces the MatCont/AUTO branch-continuation role for
  equilibria, Hopf detection, and periodic-orbit continuation.
- The double-frontier script is set up as the first executable checkpoint in the
  larger plan: confirm fold viability, continue the reduced branch, test for
  hysteresis, then decide whether the next stage should target canards, MMOs, or
  only tipping/path-dependence.
