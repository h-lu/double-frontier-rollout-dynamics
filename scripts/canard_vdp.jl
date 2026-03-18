using Printf
using LinearAlgebra

import BifurcationKit as BK
import BifurcationKit: @optic
import OrdinaryDiffEq as ODE
using Plots

mkpath(joinpath(@__DIR__, "..", "results"))

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results"))

function vdp_slowfast!(du, u, p, t = 0.0)
    (; ε, a) = p
    x, y = u
    du[1] = y - x^3 / 3 + x
    du[2] = ε * (a - x)
    return du
end

function equilibrium_guess(p)
    x = p.a
    y = x^3 / 3 - x
    return [x, y]
end

function make_problem(; ε = 0.05, a = 1.2)
    params = (ε = ε, a = a)
    u0 = equilibrium_guess(params)

    prob = BK.ODEBifProblem(
        vdp_slowfast!,
        u0,
        params,
        (@optic _.a);
        record_from_solution = (x, p; k...) -> (
            x = x[1],
            y = x[2],
            amp_guess = abs(x[1]),
        ),
    )

    return prob, params
end

function continue_equilibria(prob)
    opts = BK.ContinuationPar(
        p_min = 0.6,
        p_max = 1.2,
        ds = -1e-3,
        dsmin = 1e-6,
        dsmax = 5e-3,
        max_steps = 900,
        nev = 4,
        detect_bifurcation = 3,
        save_sol_every_step = 1,
        tol_stability = 1e-10,
    )

    br = BK.continuation(
        prob,
        BK.PALC(tangent = BK.Bordered()),
        opts;
        bothside = false,
        normC = BK.norminf,
        plot = false,
    )

    return br
end

function hopf_special_point_index(br)
    for (idx, sp) in enumerate(br.specialpoint)
        if sp.type == :hopf
            return idx
        end
    end
    error("No Hopf point detected on the equilibrium branch.")
end

function continue_periodic_orbits(br, hopf_idx)
    args_po = (
        record_from_solution = (x, p; k...) -> begin
            orbit = BK.get_periodic_orbit(p.prob, x, p.p)
            xvals = orbit[1, :]
            return (
                xmax = maximum(xvals),
                xmin = minimum(xvals),
                xamp = maximum(xvals) - minimum(xvals),
                period = BK.getperiod(p.prob, x, p.p),
            )
        end,
        normC = BK.norminf,
    )

    opts_po = BK.ContinuationPar(
        p_min = 0.6,
        p_max = 1.2,
        ds = -2e-4,
        dsmin = 1e-7,
        dsmax = 2e-3,
        max_steps = 600,
        save_sol_every_step = 1,
        tol_stability = 1e-8,
    )

    return BK.continuation(
        br,
        hopf_idx,
        opts_po,
        BK.PeriodicOrbitOCollProblem(40, 4; meshadapt = true);
        δp = -2e-4,
        plot = false,
        args_po...,
    )
end

function simulate_parameter(a; ε = 0.05, tmax = 4000.0, transient = 2000.0)
    params = (ε = ε, a = a)
    u0 = equilibrium_guess(params) .+ [0.03, 0.0]
    odeprob = ODE.ODEProblem(vdp_slowfast!, u0, (0.0, tmax), params)
    sol = ODE.solve(
        odeprob,
        ODE.Rodas5P();
        abstol = 1e-11,
        reltol = 1e-9,
        saveat = 0.1,
        maxiters = Int(1e8),
    )

    mask = sol.t .>= transient
    xvals = [u[1] for u in sol.u[mask]]
    yvals = [u[2] for u in sol.u[mask]]
    tvals = sol.t[mask]

    return (
        t = tvals,
        x = xvals,
        y = yvals,
        a = a,
        xamp = maximum(xvals) - minimum(xvals),
        xmin = minimum(xvals),
        xmax = maximum(xvals),
    )
end

function sweep_amplitudes(a_values; ε = 0.05, tmax = 2500.0, transient = 1500.0)
    sims = [simulate_parameter(a; ε = ε, tmax = tmax, transient = transient) for a in a_values]
    amps = [sim.xamp for sim in sims]
    xmins = [sim.xmin for sim in sims]
    xmaxs = [sim.xmax for sim in sims]
    return (a = collect(a_values), amp = amps, xmin = xmins, xmax = xmaxs, sims = sims)
end

function refine_canard_window(ε, hopf_a)
    coarse = sweep_amplitudes(range(hopf_a + 6e-4, hopf_a - 1.4e-2, length = 15); ε = ε)
    jump_idx = argmax(diff(coarse.amp))
    a_left = coarse.a[jump_idx]
    a_right = coarse.a[jump_idx + 1]
    refined = sweep_amplitudes(range(a_left, a_right, length = 11); ε = ε, tmax = 4000.0, transient = 2500.0)
    refined_jump_idx = argmax(diff(refined.amp))
    return coarse, refined, refined_jump_idx
end

function save_equilibrium_plot(br)
    p = plot(
        br;
        xlabel = "a",
        ylabel = "x equilibrium",
        title = "van der Pol equilibrium continuation",
        legend = :bottomleft,
        markersize = 3,
    )
    savefig(p, joinpath(RESULTS_DIR, "vdp_equilibrium_branch.png"))
end

function save_periodic_orbit_plot(br, br_po)
    p = plot(
        br,
        br_po;
        xlabel = "a",
        ylabel = "x",
        title = "Canard branch from Hopf to relaxation oscillations",
        markersize = 3,
        legend = :bottomleft,
    )
    plot!(p, br_po.param, br_po.xmin; label = "periodic xmin", lw = 2)
    plot!(p, br_po.param, br_po.xmax; label = "periodic xmax", lw = 2)
    savefig(p, joinpath(RESULTS_DIR, "vdp_canard_branch.png"))
end

function save_sweep_plot(br, coarse, refined, hopf_a)
    p1 = plot(
        coarse.a,
        coarse.amp;
        marker = :circle,
        lw = 2,
        xlabel = "a",
        ylabel = "post-transient x-amplitude",
        title = "Coarse amplitude sweep across the canard explosion",
        label = "coarse sweep",
    )
    vline!(p1, [hopf_a]; label = "Hopf", ls = :dash)

    p2 = plot(
        refined.a,
        refined.amp;
        marker = :circle,
        lw = 2,
        xlabel = "a",
        ylabel = "post-transient x-amplitude",
        title = "Refined sweep near the jump window",
        label = "refined sweep",
    )
    vline!(p2, [hopf_a]; label = "Hopf", ls = :dash)

    p3 = plot(
        br;
        xlabel = "a",
        ylabel = "x equilibrium",
        title = "Equilibrium continuation from BifurcationKit",
        legend = :bottomleft,
        markersize = 3,
    )

    p4 = plot(
        refined.a,
        refined.xmin;
        lw = 2,
        marker = :circle,
        xlabel = "a",
        ylabel = "x extrema",
        title = "Directly simulated periodic orbit envelope",
        label = "xmin",
    )
    plot!(p4, refined.a, refined.xmax; lw = 2, marker = :circle, label = "xmax")
    vline!(p4, [hopf_a]; label = "Hopf", ls = :dash)

    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(combo, joinpath(RESULTS_DIR, "vdp_canard_sweep.png"))
end

function save_simulation_plot(sim_low, sim_high)
    p1 = plot(
        sim_low.t,
        sim_low.x;
        xlabel = "t",
        ylabel = "x(t)",
        title = @sprintf("small-amplitude oscillation, a = %.6f", sim_low.a),
        label = "x",
        lw = 2,
    )
    p2 = plot(
        sim_high.t,
        sim_high.x;
        xlabel = "t",
        ylabel = "x(t)",
        title = @sprintf("relaxation oscillation, a = %.6f", sim_high.a),
        label = "x",
        lw = 2,
    )
    p3 = plot(
        sim_low.x,
        sim_low.y;
        xlabel = "x",
        ylabel = "y",
        title = "phase portrait below explosion",
        label = "trajectory",
        lw = 2,
    )
    p4 = plot(
        sim_high.x,
        sim_high.y;
        xlabel = "x",
        ylabel = "y",
        title = "phase portrait above explosion",
        label = "trajectory",
        lw = 2,
    )
    layout_plot = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(layout_plot, joinpath(RESULTS_DIR, "vdp_canard_simulations.png"))
end

function save_summary(br, br_po, coarse, refined, refined_jump_idx, sim_low, sim_high, hopf_idx)
    hopf = br.specialpoint[hopf_idx]
    path = joinpath(RESULTS_DIR, "vdp_canard_summary.txt")
    open(path, "w") do io
        println(io, "van der Pol slow-fast canard analysis")
        println(io, @sprintf("epsilon = %.6f", br.prob.params.ε))
        println(io, @sprintf("Hopf parameter a ≈ %.10f", hopf.param))
        println(io, @sprintf("Local periodic-orbit continuation points saved = %d", length(br_po.param)))
        println(io, @sprintf("Refined jump window: [%.10f, %.10f]", refined.a[refined_jump_idx], refined.a[refined_jump_idx + 1]))
        println(io, @sprintf("Refined amplitudes across jump: %.6f -> %.6f", refined.amp[refined_jump_idx], refined.amp[refined_jump_idx + 1]))
        println(io, @sprintf("Small-orbit sample: a = %.6f, amplitude = %.6f", sim_low.a, sim_low.xamp))
        println(io, @sprintf("Large-orbit sample: a = %.6f, amplitude = %.6f", sim_high.a, sim_high.xamp))
        println(io)
        println(io, "Interpretation:")
        println(io, "- BifurcationKit continuation locates the Hopf point and the local periodic branch onset.")
        println(io, "- A direct amplitude sweep shows the sharp jump from small oscillations to O(1) relaxation oscillations.")
        println(io, "- That steep jump over a very narrow parameter interval is the canard explosion signature.")
    end
end

function main()
    prob, params = make_problem()
    @info "Starting equilibrium continuation" params
    br = continue_equilibria(prob)

    hopf_idx = hopf_special_point_index(br)
    hopf = br.specialpoint[hopf_idx]
    @info "Hopf point detected" hopf_idx hopf

    @info "Starting periodic-orbit continuation from Hopf"
    br_po = continue_periodic_orbits(br, hopf_idx)

    coarse, refined, refined_jump_idx = refine_canard_window(params.ε, hopf.param)
    a_low = refined.a[refined_jump_idx]
    a_high = refined.a[refined_jump_idx + 1]

    @info "Simulating two parameters around the canard explosion window" a_low a_high
    sim_low = simulate_parameter(a_low; ε = params.ε)
    sim_high = simulate_parameter(a_high; ε = params.ε)

    save_equilibrium_plot(br)
    save_periodic_orbit_plot(br, br_po)
    save_sweep_plot(br, coarse, refined, hopf.param)
    save_simulation_plot(sim_low, sim_high)
    save_summary(br, br_po, coarse, refined, refined_jump_idx, sim_low, sim_high, hopf_idx)

    println()
    println("Saved results to: ", RESULTS_DIR)
    println(@sprintf("Hopf point: a ≈ %.10f", hopf.param))
    println(@sprintf("Refined jump: a = %.10f -> %.10f, amplitude = %.6f -> %.6f",
        refined.a[refined_jump_idx],
        refined.a[refined_jump_idx + 1],
        refined.amp[refined_jump_idx],
        refined.amp[refined_jump_idx + 1],
    ))
    println(@sprintf("Below explosion: a = %.6f, x-amplitude = %.6f", sim_low.a, sim_low.xamp))
    println(@sprintf("Above explosion: a = %.6f, x-amplitude = %.6f", sim_high.a, sim_high.xamp))
end

main()
