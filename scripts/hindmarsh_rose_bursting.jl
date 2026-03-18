using Printf
using Statistics

import BifurcationKit as BK
import BifurcationKit: @optic
import OrdinaryDiffEq as ODE
using Plots

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results"))
mkpath(RESULTS_DIR)

function hr(u, p)
    (; a, b, c, d, r, s, xR, I) = p
    x, y, z = u
    return [
        y - a * x^3 + b * x^2 - z + I,
        c - d * x^2 - y,
        r * (s * (x - xR) - z),
    ]
end

hr_ode(u, p, t) = hr(u, p)

function base_params(; I = 1.0)
    return (a = 1.0, b = 3.0, c = 1.0, d = 5.0, r = 0.006, s = 4.0, xR = -1.6, I = I)
end

function integrate_to_state(I; tmax = 3000.0)
    p = base_params(I = I)
    prob = ODE.ODEProblem(hr_ode, [-1.2, -8.0, 3.0], (0.0, tmax), p)
    sol = ODE.solve(
        prob,
        ODE.Rodas5P();
        abstol = 1e-10,
        reltol = 1e-9,
        saveat = 1.0,
        maxiters = Int(1e8),
    )
    return sol.u[end], p
end

function make_problem(; I = 1.0)
    u0, params = integrate_to_state(I)
    prob = BK.ODEBifProblem(
        hr,
        collect(u0),
        params,
        (@optic _.I);
        record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], z = x[3]),
    )
    return prob, params
end

function continue_equilibria(prob)
    opts = BK.ContinuationPar(
        p_min = 0.8,
        p_max = 2.2,
        ds = 1e-3,
        dsmin = 1e-6,
        dsmax = 5e-3,
        max_steps = 1400,
        nev = 5,
        detect_bifurcation = 3,
        save_sol_every_step = 1,
        tol_stability = 1e-10,
    )

    return BK.continuation(
        prob,
        BK.PALC(tangent = BK.Bordered()),
        opts;
        plot = false,
        normC = BK.norminf,
    )
end

function hopf_special_point_index(br)
    for (idx, sp) in enumerate(br.specialpoint)
        if sp.type == :hopf
            return idx
        end
    end
    error("No Hopf point found on the Hindmarsh-Rose equilibrium branch.")
end

function peak_times(ts, xs; threshold = 0.0)
    pts = Float64[]
    for i in 2:length(xs)-1
        if xs[i] > threshold && xs[i] > xs[i - 1] && xs[i] >= xs[i + 1]
            push!(pts, ts[i])
        end
    end
    return pts
end

function burst_sizes(pts; burst_gap = 50.0)
    isempty(pts) && return Int[]

    sizes = Int[]
    current = 1
    for dt in diff(pts)
        if dt > burst_gap
            push!(sizes, current)
            current = 1
        else
            current += 1
        end
    end
    push!(sizes, current)
    return sizes
end

function classify_regime(num_peaks, max_gap; burst_gap = 50.0)
    if num_peaks == 0
        return "rest"
    elseif max_gap > burst_gap
        return "bursting"
    else
        return "tonic spiking"
    end
end

function simulate_parameter(I; tmax = 4000.0, transient = 2000.0, saveat = 0.2)
    p = base_params(I = I)
    prob = ODE.ODEProblem(hr_ode, [-1.2, -8.0, 3.0], (0.0, tmax), p)
    sol = ODE.solve(
        prob,
        ODE.Rodas5P();
        abstol = 1e-10,
        reltol = 1e-9,
        saveat = saveat,
        maxiters = Int(1e8),
    )

    mask = sol.t .>= transient
    tvals = sol.t[mask]
    xvals = [u[1] for u in sol.u[mask]]
    yvals = [u[2] for u in sol.u[mask]]
    zvals = [u[3] for u in sol.u[mask]]

    pts = peak_times(tvals, xvals)
    gaps = length(pts) > 1 ? diff(pts) : Float64[]
    max_gap = isempty(gaps) ? 0.0 : maximum(gaps)
    bursts = burst_sizes(pts)

    return (
        I = I,
        t = tvals,
        x = xvals,
        y = yvals,
        z = zvals,
        amplitude = maximum(xvals) - minimum(xvals),
        xmin = minimum(xvals),
        xmax = maximum(xvals),
        peak_count = length(pts),
        max_interspike_gap = max_gap,
        burst_count = length(bursts),
        mean_spikes_per_burst = isempty(bursts) ? 0.0 : mean(bursts),
        regime = classify_regime(length(pts), max_gap),
    )
end

function sweep_regimes(I_values)
    sims = [simulate_parameter(I) for I in I_values]
    return (
        I = collect(I_values),
        amplitude = [sim.amplitude for sim in sims],
        max_gap = [sim.max_interspike_gap for sim in sims],
        spikes_per_burst = [sim.regime == "bursting" ? sim.mean_spikes_per_burst : NaN for sim in sims],
        regimes = [sim.regime for sim in sims],
        sims = sims,
    )
end

function pick_representatives(sweep)
    rest = findfirst(==("rest"), sweep.regimes)
    tonic = findlast(==("tonic spiking"), sweep.regimes)

    burst_candidates = [idx for (idx, regime) in enumerate(sweep.regimes) if regime == "bursting"]
    burst_scores = [sweep.sims[idx].mean_spikes_per_burst for idx in burst_candidates]
    burst = isempty(burst_candidates) ? nothing : burst_candidates[argmax(burst_scores)]

    if isnothing(rest) || isnothing(burst) || isnothing(tonic)
        error("Could not identify rest, bursting, and tonic-spiking examples from the sweep.")
    end

    return sweep.sims[rest], sweep.sims[burst], sweep.sims[tonic]
end

function save_equilibrium_plot(br)
    p = plot(
        br;
        xlabel = "I",
        ylabel = "x equilibrium",
        title = "Hindmarsh-Rose equilibrium continuation",
        legend = :bottomleft,
        markersize = 3,
    )
    savefig(p, joinpath(RESULTS_DIR, "hr_equilibrium_branch.png"))
end

function save_sweep_plot(br, sweep, hopf_I)
    colors = Dict("rest" => :steelblue, "bursting" => :tomato, "tonic spiking" => :seagreen)
    regime_colors = [colors[r] for r in sweep.regimes]

    p1 = scatter(
        sweep.I,
        sweep.amplitude;
        ms = 7,
        markerstrokewidth = 0,
        color = regime_colors,
        xlabel = "I",
        ylabel = "x amplitude",
        title = "Amplitude sweep across Hindmarsh-Rose regimes",
        label = false,
    )
    vline!(p1, [hopf_I]; ls = :dash, label = "Hopf")

    p2 = scatter(
        sweep.I,
        sweep.max_gap;
        ms = 7,
        markerstrokewidth = 0,
        color = regime_colors,
        xlabel = "I",
        ylabel = "max inter-spike gap",
        title = "Bursting has large silent intervals",
        label = false,
    )
    hline!(p2, [50.0]; ls = :dash, label = "burst threshold")

    p3 = plot(
        br;
        xlabel = "I",
        ylabel = "x equilibrium",
        title = "Equilibrium branch and Hopf onset",
        markersize = 3,
        legend = :bottomleft,
    )

    p4 = scatter(
        sweep.I,
        sweep.spikes_per_burst;
        ms = 7,
        markerstrokewidth = 0,
        color = regime_colors,
        xlabel = "I",
        ylabel = "mean spikes per burst",
        title = "Burst complexity inside the bursting regime",
        label = false,
    )

    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(combo, joinpath(RESULTS_DIR, "hr_regime_sweep.png"))
end

function save_example_plot(rest, burst, tonic)
    p1 = plot(
        rest.t,
        rest.x;
        xlabel = "t",
        ylabel = "x(t)",
        title = @sprintf("rest, I = %.2f", rest.I),
        lw = 2,
        label = "x",
    )
    p2 = plot(
        burst.t,
        burst.x;
        xlabel = "t",
        ylabel = "x(t)",
        title = @sprintf("bursting, I = %.2f", burst.I),
        lw = 2,
        label = "x",
    )
    p3 = plot(
        tonic.t,
        tonic.x;
        xlabel = "t",
        ylabel = "x(t)",
        title = @sprintf("tonic spiking, I = %.2f", tonic.I),
        lw = 2,
        label = "x",
    )
    p4 = plot(
        burst.x,
        burst.z;
        xlabel = "x",
        ylabel = "z",
        title = "bursting slow-fast loop in (x, z)",
        lw = 2,
        label = "trajectory",
    )
    p5 = plot(
        tonic.x,
        tonic.y;
        xlabel = "x",
        ylabel = "y",
        title = "tonic-spiking loop in (x, y)",
        lw = 2,
        label = "trajectory",
    )
    p6 = plot(
        burst.t,
        burst.z;
        xlabel = "t",
        ylabel = "z(t)",
        title = "slow adaptation variable during bursting",
        lw = 2,
        label = "z",
    )

    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 1100))
    savefig(combo, joinpath(RESULTS_DIR, "hr_examples.png"))
end

function save_summary(br, sweep, rest, burst, tonic, hopf_idx)
    hopf = br.specialpoint[hopf_idx]
    counts = Dict(
        "rest" => count(==("rest"), sweep.regimes),
        "bursting" => count(==("bursting"), sweep.regimes),
        "tonic spiking" => count(==("tonic spiking"), sweep.regimes),
    )

    open(joinpath(RESULTS_DIR, "hr_summary.txt"), "w") do io
        println(io, "Hindmarsh-Rose 3-variable slow-fast analysis")
        println(io, @sprintf("Hopf point in I ≈ %.10f", hopf.param))
        println(io, @sprintf("Sweep points = %d", length(sweep.I)))
        println(io, @sprintf("Regime counts: rest = %d, bursting = %d, tonic spiking = %d",
            counts["rest"], counts["bursting"], counts["tonic spiking"]))
        println(io)
        println(io, @sprintf("Representative rest: I = %.2f, amplitude = %.6f", rest.I, rest.amplitude))
        println(io, @sprintf("Representative bursting: I = %.2f, amplitude = %.6f, max gap = %.2f, mean spikes/burst = %.2f",
            burst.I, burst.amplitude, burst.max_interspike_gap, burst.mean_spikes_per_burst))
        println(io, @sprintf("Representative tonic spiking: I = %.2f, amplitude = %.6f, max gap = %.2f",
            tonic.I, tonic.amplitude, tonic.max_interspike_gap))
        println(io)
        println(io, "Interpretation:")
        println(io, "- The equilibrium branch loses stability at the Hopf point as I increases.")
        println(io, "- Just beyond that loss of stability the model exhibits pronounced bursting.")
        println(io, "- For larger I the silent intervals collapse and the system settles into tonic spiking.")
    end
end

function main()
    prob, params = make_problem()
    @info "Starting Hindmarsh-Rose equilibrium continuation" params
    br = continue_equilibria(prob)
    hopf_idx = hopf_special_point_index(br)
    hopf = br.specialpoint[hopf_idx]
    @info "Hopf point detected" hopf_idx hopf

    I_values = 1.0:0.2:4.0
    @info "Sweeping applied current values for rest/bursting/spiking regimes" first(I_values) last(I_values)
    sweep = sweep_regimes(I_values)
    rest_coarse, burst_coarse, tonic_coarse = pick_representatives(sweep)
    rest = simulate_parameter(rest_coarse.I; saveat = 0.05)
    burst = simulate_parameter(burst_coarse.I; saveat = 0.05)
    tonic = simulate_parameter(tonic_coarse.I; saveat = 0.05)

    save_equilibrium_plot(br)
    save_sweep_plot(br, sweep, hopf.param)
    save_example_plot(rest, burst, tonic)
    save_summary(br, sweep, rest, burst, tonic, hopf_idx)

    println()
    println("Saved results to: ", RESULTS_DIR)
    println(@sprintf("Hopf point: I ≈ %.10f", hopf.param))
    println(@sprintf("Representative rest: I = %.2f", rest.I))
    println(@sprintf("Representative bursting: I = %.2f, max gap = %.2f, mean spikes/burst = %.2f",
        burst.I, burst.max_interspike_gap, burst.mean_spikes_per_burst))
    println(@sprintf("Representative tonic spiking: I = %.2f, max gap = %.2f",
        tonic.I, tonic.max_interspike_gap))
end

main()
