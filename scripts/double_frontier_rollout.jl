using Printf

import BifurcationKit as BK
import BifurcationKit: @optic
import OrdinaryDiffEq as ODE
using Plots

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results"))
mkpath(RESULTS_DIR)

stable_exp(x) = exp(clamp(x, -60.0, 60.0))
logistic(x) = inv(1 + stable_exp(-x))

function softplus_sigma(x, sigma)
    z = sigma * x
    if z > 40
        return x
    elseif z < -40
        return stable_exp(z) / sigma
    else
        return log1p(exp(z)) / sigma
    end
end

function smin_sigma(x, y, sigma)
    ax = -sigma * x
    ay = -sigma * y
    m = max(ax, ay)
    return -(m + log(exp(ax - m) + exp(ay - m))) / sigma
end

function base_params(;
    E = 0.2,
    G = 1.0,
    z = 1.0,
    epsilon = 0.03,
    eta_h = 0.01,
    eta_k = 0.03,
    beta_n = 0.7,
    lambda_D = 12.0,
    theta0 = -2.0,
    sigma = 80.0,
    hbar = 0.65,
    kbar = 0.30,
)
    return (
        lambda_A = 8.0,
        c_A = 0.0,
        alpha0 = -1.0,
        alpha_z = 1.0,
        alpha_h = 1.0,
        alpha_G = 1.0,
        lambda_D = lambda_D,
        c_D = 0.0,
        theta0 = theta0,
        theta_z = 1.0,
        theta_h = 0.8,
        theta_G = 0.8,
        theta_E = 1.2,
        theta_k = 1.0,
        theta_phi = 1.0,
        beta_n = beta_n,
        gamma_q = 1.4,
        omega0 = -1.0,
        omega_G = 1.2,
        omega_h = 1.0,
        omega_q = 1.5,
        nu0 = 0.3,
        nu_G = 0.5,
        nu_h = 0.4,
        chi_I = 0.5,
        chi_R = 0.8,
        rho = 2.0,
        eta_k = eta_k,
        delta_k = 0.2,
        psi_M = 0.7,
        psi_D = 0.6,
        eta_h = eta_h,
        vartheta_M = 0.7,
        vartheta_D = 0.35,
        vartheta_H = 0.45,
        mu0 = 0.05,
        mu1 = 1.0,
        sigma = sigma,
        epsilon = epsilon,
        E = E,
        G = G,
        z = z,
        hbar = hbar,
        kbar = kbar,
    )
end

function advisor_frontier(h, p)
    score = p.alpha0 + p.alpha_z * p.z + p.alpha_h * h + p.alpha_G * p.G
    return logistic(p.lambda_A * (score - p.c_A))
end

function supervision(h, q, p)
    score = p.omega0 + p.omega_G * p.G + p.omega_h * h - p.omega_q * q
    return logistic(score)
end

function doer_target(d, q, k, h, p)
    a = advisor_frontier(h, p)
    phi = supervision(h, q, p)
    score =
        p.theta0 +
        p.theta_z * p.z +
        p.theta_h * h +
        p.theta_G * p.G +
        p.theta_E * p.E +
        p.theta_k * k +
        p.theta_phi * phi -
        p.gamma_q * q +
        p.beta_n * d
    raw = logistic(p.lambda_D * (score - p.c_D))
    return smin_sigma(raw, a, p.sigma)
end

function double_frontier_full!(du, u, p, t = 0.0)
    d, q, k, h = u
    a = advisor_frontier(h, p)
    phi = supervision(h, q, p)
    dstar = doer_target(d, q, k, h, p)

    du[1] = (dstar - d) / p.epsilon
    du[2] = p.chi_I * softplus_sigma(dstar - d, p.sigma) + p.chi_R * d^p.rho - (p.nu0 + p.nu_G * p.G + p.nu_h * h) * q
    du[3] = p.eta_k * (p.psi_M * softplus_sigma(a - d, p.sigma) + p.psi_D * phi * d - p.delta_k * k)
    du[4] = p.eta_h * (
        (p.vartheta_M * softplus_sigma(a - d, p.sigma) + p.vartheta_D * phi * d + p.vartheta_H * (1 - a)) * (1 - h) -
        (p.mu0 + p.mu1 * (1 - phi) * d) * h
    )
    return du
end

function double_frontier_3d!(du, u, p, t = 0.0)
    d, q, h = u
    k = p.kbar
    a = advisor_frontier(h, p)
    phi = supervision(h, q, p)
    dstar = doer_target(d, q, k, h, p)

    du[1] = (dstar - d) / p.epsilon
    du[2] = p.chi_I * softplus_sigma(dstar - d, p.sigma) + p.chi_R * d^p.rho - (p.nu0 + p.nu_G * p.G + p.nu_h * h) * q
    du[3] = p.eta_h * (
        (p.vartheta_M * softplus_sigma(a - d, p.sigma) + p.vartheta_D * phi * d + p.vartheta_H * (1 - a)) * (1 - h) -
        (p.mu0 + p.mu1 * (1 - phi) * d) * h
    )
    return du
end

function double_frontier_2d!(du, u, p, t = 0.0)
    d, q = u
    h = p.hbar
    k = p.kbar
    dstar = doer_target(d, q, k, h, p)

    du[1] = (dstar - d) / p.epsilon
    du[2] = p.chi_I * softplus_sigma(dstar - d, p.sigma) + p.chi_R * d^p.rho - (p.nu0 + p.nu_G * p.G + p.nu_h * h) * q
    return du
end

function rhs_2d(u, p)
    du = zeros(2)
    double_frontier_2d!(du, u, p, 0.0)
    return du
end

function integrate_to_state_2d(p; u0 = [0.05, 0.0], tmax = 1500.0)
    prob = ODE.ODEProblem(double_frontier_2d!, collect(u0), (0.0, tmax), p)
    sol = ODE.solve(
        prob,
        ODE.Rodas5P();
        abstol = 1e-10,
        reltol = 1e-9,
        saveat = tmax,
        maxiters = Int(1e8),
    )
    return copy(sol.u[end])
end

function simulate_2d(p; u0 = [0.05, 0.0], tmax = 3000.0, transient = 1500.0, saveat = 0.5)
    prob = ODE.ODEProblem(double_frontier_2d!, collect(u0), (0.0, tmax), p)
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
    dvals = [u[1] for u in sol.u[mask]]
    qvals = [u[2] for u in sol.u[mask]]

    return (
        E = p.E,
        t = tvals,
        d = dvals,
        q = qvals,
        damp = maximum(dvals) - minimum(dvals),
        qamp = maximum(qvals) - minimum(qvals),
        dmin = minimum(dvals),
        dmax = maximum(dvals),
        qmin = minimum(qvals),
        qmax = maximum(qvals),
    )
end

function simulate_full(p; u0 = [0.05, 0.0, 0.1, 0.6], tmax = 6000.0, transient = 3000.0, saveat = 1.0)
    prob = ODE.ODEProblem(double_frontier_full!, collect(u0), (0.0, tmax), p)
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
    dvals = [u[1] for u in sol.u[mask]]
    qvals = [u[2] for u in sol.u[mask]]
    kvals = [u[3] for u in sol.u[mask]]
    hvals = [u[4] for u in sol.u[mask]]
    mvals = [max(advisor_frontier(hvals[i], p) - dvals[i], 0.0) for i in eachindex(dvals)]

    return (
        E = p.E,
        G = p.G,
        t = tvals,
        d = dvals,
        q = qvals,
        k = kvals,
        h = hvals,
        m = mvals,
        damp = maximum(dvals) - minimum(dvals),
        qamp = maximum(qvals) - minimum(qvals),
    )
end

function make_2d_problem(p)
    u0 = integrate_to_state_2d(p)
    base_p = p
    prob = BK.ODEBifProblem(
        double_frontier_2d!,
        u0,
        p,
        (@optic _.E);
        record_from_solution = (x, E; k...) -> begin
            par = merge(base_p, (E = E,))
            return (
                d = x[1],
                q = x[2],
                a_cap = advisor_frontier(par.hbar, par),
                phi = supervision(par.hbar, x[2], par),
            )
        end,
    )
    return prob
end

function continue_equilibria(prob)
    opts = BK.ContinuationPar(
        p_min = 0.0,
        p_max = 2.0,
        ds = 5e-3,
        dsmin = 1e-6,
        dsmax = 2e-2,
        max_steps = 700,
        nev = 4,
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

function find_special_point_index(br, target)
    for (idx, sp) in enumerate(br.specialpoint)
        if sp.type == target
            return idx
        end
    end
    return nothing
end

function continue_periodic_from_hopf(br, hopf_idx)
    args_po = (
        record_from_solution = (x, par; k...) -> begin
            orbit = BK.get_periodic_orbit(par.prob, x, par.p)
            dvals = orbit[1, :]
            return (
                dmax = maximum(dvals),
                dmin = minimum(dvals),
                damp = maximum(dvals) - minimum(dvals),
                period = BK.getperiod(par.prob, x, par.p),
            )
        end,
        normC = BK.norminf,
    )

    opts_po = BK.ContinuationPar(
        p_min = 0.0,
        p_max = 2.0,
        ds = 1e-4,
        dsmin = 1e-7,
        dsmax = 1e-3,
        max_steps = 300,
        save_sol_every_step = 1,
        tol_stability = 1e-8,
    )

    return BK.continuation(
        br,
        hopf_idx,
        opts_po,
        BK.PeriodicOrbitOCollProblem(40, 4; meshadapt = true);
        δp = 1e-4,
        plot = false,
        args_po...,
    )
end

function sweep_multistability(E_values, p)
    low_branch = Float64[]
    high_branch = Float64[]
    amp_low = Float64[]
    amp_high = Float64[]

    a_cap = advisor_frontier(p.hbar, p)
    high_init = [min(0.95, max(0.15, a_cap - 0.02)), 0.5]

    for E in E_values
        par = merge(p, (E = E,))
        sim_low = simulate_2d(par; u0 = [0.05, 0.0], tmax = 2500.0, transient = 1200.0)
        sim_high = simulate_2d(par; u0 = high_init, tmax = 2500.0, transient = 1200.0)
        push!(low_branch, sim_low.d[end])
        push!(high_branch, sim_high.d[end])
        push!(amp_low, sim_low.damp)
        push!(amp_high, sim_high.damp)
    end

    return (
        E = collect(E_values),
        low_branch = low_branch,
        high_branch = high_branch,
        amp_low = amp_low,
        amp_high = amp_high,
        gap = abs.(high_branch .- low_branch),
    )
end

function fold_indicator(beta_n, lambda_D)
    return beta_n * lambda_D > 4 ? 1.0 : 0.0
end

function save_fold_screen_plot()
    beta_vals = range(0.0, 1.5, length = 200)
    lambda_vals = range(4.0, 30.0, length = 220)
    indicator = [fold_indicator(beta, lambda) for lambda in lambda_vals, beta in beta_vals]

    p = heatmap(
        beta_vals,
        lambda_vals,
        indicator;
        xlabel = "beta_n",
        ylabel = "lambda_D",
        title = "Theoretical fold screen: beta_n * lambda_D > 4",
        color = cgrad([:white, :tomato]),
        colorbar = false,
    )
    plot!(p, beta_vals, 4.0 ./ max.(beta_vals, 1e-3); lw = 2, color = :black, label = "beta_n * lambda_D = 4")
    xlims!(p, extrema(beta_vals))
    ylims!(p, extrema(lambda_vals))
    savefig(p, joinpath(RESULTS_DIR, "double_frontier_fold_screen.png"))
end

function save_branch_plot(br, br_po)
    if br_po === nothing
        p = plot(
            br;
            xlabel = "E",
            ylabel = "d equilibrium",
            title = "2D reduced model equilibrium continuation",
            markersize = 3,
            legend = :bottomleft,
        )
    else
        p = plot(
            br,
            br_po;
            xlabel = "E",
            ylabel = "d",
            title = "Equilibria and periodic orbits in the 2D reduced model",
            markersize = 3,
            legend = :bottomleft,
        )
        plot!(p, br_po.param, br_po.dmin; lw = 2, label = "periodic dmin")
        plot!(p, br_po.param, br_po.dmax; lw = 2, label = "periodic dmax")
    end
    savefig(p, joinpath(RESULTS_DIR, "double_frontier_branch.png"))
end

function save_hysteresis_plot(sweep)
    p1 = plot(
        sweep.E,
        sweep.low_branch;
        lw = 2,
        marker = :circle,
        xlabel = "E",
        ylabel = "post-transient d",
        title = "Low vs high initial conditions",
        label = "low init",
    )
    plot!(p1, sweep.E, sweep.high_branch; lw = 2, marker = :square, label = "high init")

    p2 = plot(
        sweep.E,
        sweep.gap;
        lw = 2,
        marker = :circle,
        xlabel = "E",
        ylabel = "|d_high - d_low|",
        title = "Initial-condition sensitivity",
        label = "branch gap",
    )

    p3 = plot(
        sweep.E,
        sweep.amp_low;
        lw = 2,
        marker = :circle,
        xlabel = "E",
        ylabel = "d amplitude",
        title = "Oscillation amplitude by initial condition",
        label = "low init",
    )
    plot!(p3, sweep.E, sweep.amp_high; lw = 2, marker = :square, label = "high init")

    combo = plot(p1, p2, p3; layout = (3, 1), size = (900, 1100))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_hysteresis.png"))
end

function save_phase_plot(sim, p)
    d_grid = range(0.0, 1.02, length = 180)
    q_grid = range(0.0, 1.8, length = 180)
    f1 = [rhs_2d([d, q], p)[1] for q in q_grid, d in d_grid]
    f2 = [rhs_2d([d, q], p)[2] for q in q_grid, d in d_grid]

    p1 = contour(
        d_grid,
        q_grid,
        f1;
        levels = [0.0],
        linewidth = 2,
        color = :royalblue,
        xlabel = "d",
        ylabel = "q",
        title = @sprintf("2D phase plane at E = %.3f", p.E),
        label = "d-nullcline",
        colorbar = false,
    )
    contour!(p1, d_grid, q_grid, f2; levels = [0.0], linewidth = 2, color = :darkorange, label = "q-nullcline", colorbar = false)
    plot!(p1, sim.d, sim.q; lw = 2, color = :black, label = "trajectory")

    p2 = plot(
        sim.t,
        sim.d;
        lw = 2,
        xlabel = "t",
        ylabel = "d(t)",
        title = "Reduced-model time series",
        label = "d",
    )
    plot!(p2, sim.t, sim.q; lw = 2, label = "q")

    combo = plot(p1, p2; layout = (1, 2), size = (1200, 500))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_phase.png"))
end

function save_full_plot(sim)
    p1 = plot(sim.t, sim.d; lw = 2, xlabel = "t", ylabel = "d(t)", title = "Doer share", label = "d")
    p2 = plot(sim.t, sim.q; lw = 2, xlabel = "t", ylabel = "q(t)", title = "Supervision backlog", label = "q")
    p3 = plot(sim.t, sim.k; lw = 2, xlabel = "t", ylabel = "k(t)", title = "Adoption capital", label = "k")
    p4 = plot(sim.t, sim.h; lw = 2, xlabel = "t", ylabel = "h(t)", title = "Judgment capital", label = "h")
    p5 = plot(sim.t, sim.m; lw = 2, xlabel = "t", ylabel = "m(t)", title = "Middle band", label = "m")
    p6 = plot(sim.k, sim.h; lw = 2, xlabel = "k", ylabel = "h", title = "Slow drift in (k, h)", label = "trajectory")
    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_full_timeseries.png"))
end

function save_summary(br, br_po, sweep, rep_sim, full_sim, p, hopf_idx, fold_idx)
    path = joinpath(RESULTS_DIR, "double_frontier_summary.txt")
    open(path, "w") do io
        println(io, "Double-frontier rollout model: Julia reduced-model workflow")
        println(io, @sprintf("Continuation parameter: E in [%.2f, %.2f]", 0.0, 2.0))
        println(io, @sprintf("Theoretical fold screen: beta_n * lambda_D = %.4f", p.beta_n * p.lambda_D))
        println(io, @sprintf("Advisor cap at hbar = %.3f: a(hbar) = %.6f", p.hbar, advisor_frontier(p.hbar, p)))
        println(io)

        if !isnothing(fold_idx)
            fold = br.specialpoint[fold_idx]
            println(io, @sprintf("Fold/branch-point candidate (%s) detected near E ≈ %.10f", String(fold.type), fold.param))
        else
            println(io, "Fold detected near E ≈ none on the current equilibrium branch")
        end

        if !isnothing(hopf_idx)
            hopf = br.specialpoint[hopf_idx]
            println(io, @sprintf("Hopf detected near E ≈ %.10f", hopf.param))
        else
            println(io, "Hopf detected near E ≈ none on the current equilibrium branch")
        end

        if br_po === nothing
            println(io, "Periodic-orbit continuation: not run because no robust Hopf point was available")
        else
            println(io, @sprintf("Periodic-orbit continuation points saved = %d", length(br_po.param)))
        end

        println(io)
        println(io, @sprintf("Largest low/high branch gap in direct sweeps = %.6f", maximum(sweep.gap)))
        println(io, @sprintf("Representative 2D amplitude at E = %.3f: d_amp = %.6f, q_amp = %.6f", rep_sim.E, rep_sim.damp, rep_sim.qamp))
        println(io, @sprintf("Representative 4D amplitude at E = %.3f, G = %.3f: d_amp = %.6f, q_amp = %.6f", full_sim.E, full_sim.G, full_sim.damp, full_sim.qamp))
        println(io)
        println(io, "Interpretation:")
        println(io, "- The script keeps the full 4D model available, but starts with the 2D reduced system for fold/Hopf screening.")
        println(io, "- Theoretical fold viability is controlled by beta_n * lambda_D > 4, matching the analytic pre-screen in the plan.")
        println(io, "- Direct low/high-initial-condition sweeps provide a Julia-native hysteresis check before moving to deeper MMO or bursting studies.")
    end
end

function main()
    params = base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0)
    prob = make_2d_problem(params)
    @info "Starting 2D equilibrium continuation" params
    br = continue_equilibria(prob)

    fold_idx = find_special_point_index(br, :fold)
    if isnothing(fold_idx)
        fold_idx = find_special_point_index(br, :bp)
    end
    hopf_idx = find_special_point_index(br, :hopf)
    br_po = nothing

    if !isnothing(hopf_idx)
        hopf = br.specialpoint[hopf_idx]
        @info "Hopf point detected on reduced model branch" hopf
        if hopf.precision < 1e-3
            try
                br_po = continue_periodic_from_hopf(br, hopf_idx)
            catch err
                @warn "Periodic-orbit continuation failed; keeping equilibrium-only outputs" err
            end
        else
            @info "Skipping periodic-orbit continuation because the current Hopf localization needs refinement" hopf_precision = hopf.precision
        end
    else
        @info "No Hopf point detected on the current reduced-model branch"
    end

    E_values = range(0.0, 2.0, length = 31)
    sweep = sweep_multistability(E_values, params)

    rep_E = !isnothing(hopf_idx) ? min(br.specialpoint[hopf_idx].param + 0.08, 2.0) : 1.1
    rep_params = merge(params, (E = rep_E,))
    rep_sim = simulate_2d(rep_params; u0 = [0.08, 0.0], tmax = 4000.0, transient = 2200.0, saveat = 0.2)

    full_params = merge(params, (E = 1.05, G = 1.1, eta_k = 0.04, eta_h = 0.015))
    full_sim = simulate_full(full_params; transient = 0.0)

    save_fold_screen_plot()
    save_branch_plot(br, br_po)
    save_hysteresis_plot(sweep)
    save_phase_plot(rep_sim, rep_params)
    save_full_plot(full_sim)
    save_summary(br, br_po, sweep, rep_sim, full_sim, params, hopf_idx, fold_idx)

    println()
    println("Saved results to: ", RESULTS_DIR)
    println(@sprintf("Fold-screen product beta_n * lambda_D = %.4f", params.beta_n * params.lambda_D))
    if !isnothing(fold_idx)
        println(@sprintf("Fold/branch-point candidate (%s): E ≈ %.10f", String(br.specialpoint[fold_idx].type), br.specialpoint[fold_idx].param))
    else
        println("Fold point: none detected on the current continuation branch")
    end
    if !isnothing(hopf_idx)
        println(@sprintf("Hopf point: E ≈ %.10f", br.specialpoint[hopf_idx].param))
    else
        println("Hopf point: none detected on the current continuation branch")
    end
    println(@sprintf("Max low/high branch gap from direct sweeps = %.6f", maximum(sweep.gap)))
    println(@sprintf("Representative 4D amplitudes: d = %.6f, q = %.6f", full_sim.damp, full_sim.qamp))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
