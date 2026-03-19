using Printf
using Base.Threads: @threads, nthreads
using LinearAlgebra
using Statistics

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

function integrate_to_state_3d(p; u0 = [0.05, 0.0, 0.6], tmax = 2500.0)
    prob = ODE.ODEProblem(double_frontier_3d!, collect(u0), (0.0, tmax), p)
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

function simulate_3d(p; u0 = [0.05, 0.0, 0.6], tmax = 5000.0, transient = 2500.0, saveat = 0.5)
    prob = ODE.ODEProblem(double_frontier_3d!, collect(u0), (0.0, tmax), p)
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
    hvals = [u[3] for u in sol.u[mask]]

    return (
        E = p.E,
        t = tvals,
        d = dvals,
        q = qvals,
        h = hvals,
        damp = maximum(dvals) - minimum(dvals),
        qamp = maximum(qvals) - minimum(qvals),
        hamp = maximum(hvals) - minimum(hvals),
        dmean = mean(dvals),
        qmean = mean(qvals),
        hmean = mean(hvals),
        dmin = minimum(dvals),
        dmax = maximum(dvals),
        qmin = minimum(qvals),
        qmax = maximum(qvals),
        hmin = minimum(hvals),
        hmax = maximum(hvals),
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

function make_3d_problem(p)
    u0 = integrate_to_state_3d(p)
    base_p = p
    prob = BK.ODEBifProblem(
        double_frontier_3d!,
        u0,
        p,
        (@optic _.E);
        record_from_solution = (x, E; k...) -> begin
            par = merge(base_p, (E = E,))
            return (
                d = x[1],
                q = x[2],
                h = x[3],
                a_cap = advisor_frontier(x[3], par),
                phi = supervision(x[3], x[2], par),
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

function record_periodic_solution(x, par; k...)
    orbit = BK.get_periodic_orbit(par.prob, x, par.p)
    dvals = orbit[1, :]
    return (
        dmax = maximum(dvals),
        dmin = minimum(dvals),
        damp = maximum(dvals) - minimum(dvals),
        period = BK.getperiod(par.prob, x, par.p),
    )
end

function collect_periodic_branch_data(segments)
    param = Float64[]
    dmin = Float64[]
    dmax = Float64[]
    damp = Float64[]
    period = Float64[]
    stable = Bool[]
    specialpoints = NamedTuple[]
    global_idx = 0

    for (seg_idx, segment) in enumerate(segments)
        branch = segment.branch
        skip = min(segment.skip, length(branch.param))
        start = skip + 1

        append!(param, branch.param[start:end])
        append!(dmin, branch.dmin[start:end])
        append!(dmax, branch.dmax[start:end])
        append!(damp, branch.damp[start:end])
        append!(period, branch.period[start:end])
        append!(stable, branch.stable[start:end])

        for sp in branch.specialpoint
            if sp.type == :endpoint && seg_idx != length(segments)
                continue
            end
            if sp.idx <= skip
                continue
            end
            push!(specialpoints, (type = sp.type, param = sp.param, idx = global_idx + sp.idx - skip, status = sp.status))
        end

        global_idx += length(branch.param) - skip
    end

    return (
        param = param,
        dmin = dmin,
        dmax = dmax,
        damp = damp,
        period = period,
        stable = stable,
        specialpoints = specialpoints,
    )
end

function dedupe_specialpoints(specialpoints; tol = 1e-9)
    ordered = sort!(collect(specialpoints); by = sp -> (string(sp.type), sp.param))
    deduped = NamedTuple[]
    for sp in ordered
        if isempty(deduped) || sp.type != deduped[end].type || abs(sp.param - deduped[end].param) > tol
            push!(deduped, sp)
        end
    end
    return deduped
end

function restart_periodic_branch(branch; x0_idx, x1_idx, theta = 0.95, ds = 5e-4, dsmin = 1e-8, dsmax = 2e-3, max_steps = 60, p_min = 0.90, p_max = 1.30, nev = 15)
    x0 = branch.sol[x0_idx].x.sol
    p0 = branch.sol[x0_idx].p
    x1 = branch.sol[x1_idx].x.sol
    p1 = branch.sol[x1_idx].p
    par0 = merge(BK.getparams(branch.prob), (E = p0,))
    opts = BK.ContinuationPar(
        branch.contparams;
        ds = abs(ds),
        dsmin = dsmin,
        dsmax = dsmax,
        max_steps = max_steps,
        p_min = p_min,
        p_max = p_max,
        nev = nev,
        save_sol_every_step = 1,
        tol_stability = 1e-9,
    )

    return BK.continuation(
        branch.prob,
        x0,
        par0,
        x1,
        p1,
        BK.PALC(θ = theta),
        BK.getlens(branch.prob),
        opts;
        kind = BK.PeriodicOrbitCont(),
        plot = false,
        verbosity = 0,
        callback_newton = BK.cbMaxNorm(1e2),
        normC = BK.norminf,
    )
end

function extend_periodic_branch(branch; theta = 0.95, target_damp = 0.9, max_stages = 8)
    segments = NamedTuple[]
    push!(segments, (branch = branch, skip = 0))
    failures = NamedTuple[]

    if length(branch.param) < 2
        return (segments = segments, data = collect_periodic_branch_data(segments), failures = failures)
    end

    current = branch

    try
        current = restart_periodic_branch(
            branch;
            x0_idx = 2,
            x1_idx = 1,
            theta = theta,
            ds = 5e-4,
            dsmin = 1e-7,
            dsmax = 2e-3,
            max_steps = 40,
            p_min = 0.90,
            p_max = 1.10,
            nev = 15,
        )
        push!(segments, (branch = current, skip = 2))
    catch err
        push!(failures, (stage = "restart-1", error = sprint(showerror, err)))
        return (segments = segments, data = collect_periodic_branch_data(segments), failures = failures)
    end

    for stage in 2:max_stages
        if maximum(current.damp) >= target_damp
            break
        end

        try
            current = restart_periodic_branch(
                current;
                x0_idx = length(current.param) - 1,
                x1_idx = length(current.param),
                theta = theta,
                ds = 5e-4,
                dsmin = 1e-8,
                dsmax = 2e-3,
                max_steps = 60,
                p_min = 0.90,
                p_max = 1.30,
                nev = 15,
            )
            push!(segments, (branch = current, skip = 2))
        catch err
            push!(failures, (stage = @sprintf("restart-%d", stage), error = sprint(showerror, err)))
            break
        end
    end

    return (segments = segments, data = collect_periodic_branch_data(segments), failures = failures)
end

function refine_hopf_normal_form(prob, br, hopf_idx, base_params)
    hopfsol = BK.newton(
        br,
        hopf_idx;
        options = BK.NewtonPar(tol = 1e-12, max_iterations = 25, verbose = false),
        normN = BK.norminf,
    )

    x = Vector{Float64}(hopfsol.u.u)
    E = hopfsol.u.p[1]
    ω_guess = abs(hopfsol.u.p[2])
    par = merge(base_params, (E = E,))

    L = Matrix(BK.jacobian(prob, x, par))
    ev = eigen(L)
    ev_idx = argmin(abs.(real.(ev.values)) .+ abs.(abs.(imag.(ev.values)) .- ω_guess))
    λ = ev.values[ev_idx]
    ζ = ComplexF64.(ev.vectors[:, ev_idx])

    evl = eigen(adjoint(L))
    evl_idx = argmin(abs.(evl.values .- conj(λ)))
    ζstar = ComplexF64.(evl.vectors[:, evl_idx])
    ζstar ./= dot(ζ, ζstar)

    hopf0 = BK.Hopf(
        x,
        br.specialpoint[hopf_idx].τ,
        E,
        imag(λ),
        par,
        (@optic _.E),
        ζ,
        ζstar,
        BK.HopfNormalForm(a = missing, b = missing, Ψ110 = missing, Ψ001 = missing, Ψ200 = missing),
        Symbol("?"),
    )

    hopf = BK.__hopf_normal_form(prob, hopf0, BK.DefaultLS(); verbose = false, L = L, autodiff = true)
    return (solution = hopfsol, point = hopf)
end

function continue_periodic_from_hopf(prob, br, hopf_idx, base_params)
    refined = refine_hopf_normal_form(prob, br, hopf_idx, base_params)
    failures = NamedTuple[]
    best = nothing

    strategies = (
        (
            name = "subcritical-left",
            cont = BK.ContinuationPar(
                p_min = 0.90,
                p_max = 1.05,
                ds = -1e-3,
                dsmin = 1e-5,
                dsmax = 1e-2,
                max_steps = 80,
                nev = 7,
                detect_bifurcation = 3,
                save_sol_every_step = 1,
                tol_stability = 1e-8,
                newton_options = BK.NewtonPar(tol = 1e-10, max_iterations = 20, verbose = false),
            ),
            coll = BK.PeriodicOrbitOCollProblem(20, 5; meshadapt = true, jacobian = BK.AutoDiffDense()),
            δp = -5e-4,
            usedeflation = false,
        ),
        (
            name = "subcritical-left-small-step",
            cont = BK.ContinuationPar(
                p_min = 0.90,
                p_max = 1.05,
                ds = -5e-4,
                dsmin = 1e-6,
                dsmax = 5e-3,
                max_steps = 40,
                nev = 15,
                detect_bifurcation = 3,
                save_sol_every_step = 1,
                tol_stability = 1e-9,
                newton_options = BK.NewtonPar(tol = 1e-8, max_iterations = 25, verbose = false),
            ),
            coll = BK.PeriodicOrbitOCollProblem(30, 5; meshadapt = true, jacobian = BK.AutoDiffDense()),
            δp = -2e-4,
            usedeflation = false,
        ),
    )

    for strategy in strategies
        try
            branch = BK._po_from_hopf(
                prob,
                refined.point,
                strategy.cont,
                strategy.coll;
                verbose = false,
                alg = BK.PALC(),
                δp = strategy.δp,
                usedeflation = strategy.usedeflation,
                eigsolver = BK.FloquetColl(),
                plot = false,
                verbosity = 0,
                callback_newton = BK.cbMaxNorm(1e2),
                record_from_solution = record_periodic_solution,
                normC = BK.norminf,
            )

            candidate = (
                name = strategy.name,
                branch = branch,
                points = length(branch.param),
                max_damp = maximum(branch.damp),
            )

            if isnothing(best) || candidate.points > best.points || (candidate.points == best.points && candidate.max_damp > best.max_damp)
                best = candidate
            end
        catch err
            push!(failures, (name = strategy.name, error = sprint(showerror, err)))
        end
    end

    if !isnothing(best)
        extension = extend_periodic_branch(best.branch.γ)
        best = (
            name = best.name,
            branch = best.branch,
            local_points = best.points,
            local_max_damp = best.max_damp,
            extension = extension,
            data = extension.data,
            points = length(extension.data.param),
            max_damp = maximum(extension.data.damp),
        )
        for failure in extension.failures
            push!(failures, (name = "$(best.name)-$(failure.stage)", error = failure.error))
        end
    end

    return (refined = refined, best = best, failures = failures)
end

function sweep_multistability(E_values, p)
    E = collect(E_values)
    low_branch = Vector{Float64}(undef, length(E))
    high_branch = Vector{Float64}(undef, length(E))
    amp_low = Vector{Float64}(undef, length(E))
    amp_high = Vector{Float64}(undef, length(E))

    a_cap = advisor_frontier(p.hbar, p)
    high_init = [min(0.95, max(0.15, a_cap - 0.02)), 0.5]

    @threads for idx in eachindex(E)
        par = merge(p, (E = E[idx],))
        sim_low = simulate_2d(par; u0 = [0.05, 0.0], tmax = 2500.0, transient = 1200.0)
        sim_high = simulate_2d(par; u0 = high_init, tmax = 2500.0, transient = 1200.0)
        low_branch[idx] = sim_low.d[end]
        high_branch[idx] = sim_high.d[end]
        amp_low[idx] = sim_low.damp
        amp_high[idx] = sim_high.damp
    end

    return (
        E = E,
        low_branch = low_branch,
        high_branch = high_branch,
        amp_low = amp_low,
        amp_high = amp_high,
        gap = abs.(high_branch .- low_branch),
    )
end

function sweep_multistability_3d(E_values, p)
    E = collect(E_values)
    low_d = Vector{Float64}(undef, length(E))
    high_d = Vector{Float64}(undef, length(E))
    low_h = Vector{Float64}(undef, length(E))
    high_h = Vector{Float64}(undef, length(E))
    amp_low = Vector{Float64}(undef, length(E))
    amp_high = Vector{Float64}(undef, length(E))
    hamp_low = Vector{Float64}(undef, length(E))
    hamp_high = Vector{Float64}(undef, length(E))

    high_init = [0.95, 0.5, 0.9]

    @threads for idx in eachindex(E)
        par = merge(p, (E = E[idx],))
        sim_low = simulate_3d(par; u0 = [0.05, 0.0, 0.2], tmax = 4000.0, transient = 1800.0)
        sim_high = simulate_3d(par; u0 = high_init, tmax = 4000.0, transient = 1800.0)
        low_d[idx] = sim_low.dmean
        high_d[idx] = sim_high.dmean
        low_h[idx] = sim_low.hmean
        high_h[idx] = sim_high.hmean
        amp_low[idx] = sim_low.damp
        amp_high[idx] = sim_high.damp
        hamp_low[idx] = sim_low.hamp
        hamp_high[idx] = sim_high.hamp
    end

    return (
        E = E,
        low_d = low_d,
        high_d = high_d,
        low_h = low_h,
        high_h = high_h,
        amp_low = amp_low,
        amp_high = amp_high,
        hamp_low = hamp_low,
        hamp_high = hamp_high,
        gap_d = abs.(high_d .- low_d),
        gap_h = abs.(high_h .- low_h),
    )
end

function scan_3d_hopf_neighborhood(E_values, p)
    E = collect(E_values)
    low_damp = Vector{Float64}(undef, length(E))
    high_damp = Vector{Float64}(undef, length(E))
    low_hamp = Vector{Float64}(undef, length(E))
    high_hamp = Vector{Float64}(undef, length(E))
    low_dmean = Vector{Float64}(undef, length(E))
    high_dmean = Vector{Float64}(undef, length(E))
    low_hmean = Vector{Float64}(undef, length(E))
    high_hmean = Vector{Float64}(undef, length(E))

    @threads for idx in eachindex(E)
        par = merge(p, (E = E[idx],))
        sim_low = simulate_3d(par; u0 = [0.05, 0.0, 0.2], tmax = 6000.0, transient = 3000.0, saveat = 0.2)
        sim_high = simulate_3d(par; u0 = [0.95, 0.5, 0.9], tmax = 6000.0, transient = 3000.0, saveat = 0.2)
        low_damp[idx] = sim_low.damp
        high_damp[idx] = sim_high.damp
        low_hamp[idx] = sim_low.hamp
        high_hamp[idx] = sim_high.hamp
        low_dmean[idx] = sim_low.dmean
        high_dmean[idx] = sim_high.dmean
        low_hmean[idx] = sim_low.hmean
        high_hmean[idx] = sim_high.hmean
    end

    return (
        E = E,
        low_damp = low_damp,
        high_damp = high_damp,
        low_hamp = low_hamp,
        high_hamp = high_hamp,
        low_dmean = low_dmean,
        high_dmean = high_dmean,
        low_hmean = low_hmean,
        high_hmean = high_hmean,
    )
end

function summarize_3d_hopf_scan(scan; amp_threshold = 0.1)
    osc = max.(scan.low_damp, scan.high_damp) .> amp_threshold
    last_osc = findlast(identity, osc)
    first_quiet = isnothing(last_osc) ? nothing : findfirst(.!osc[last_osc+1:end])
    quiet_idx = isnothing(first_quiet) ? nothing : last_osc + first_quiet
    return (
        last_osc_idx = last_osc,
        quiet_idx = quiet_idx,
        last_osc_E = isnothing(last_osc) ? nothing : scan.E[last_osc],
        first_quiet_E = isnothing(quiet_idx) ? nothing : scan.E[quiet_idx],
    )
end

function threshold_crossings(t, x; thresh = 0.9)
    out = Float64[]
    for i in 2:length(x)
        if x[i-1] < thresh <= x[i]
            α = (thresh - x[i-1]) / (x[i] - x[i-1])
            push!(out, t[i-1] + α * (t[i] - t[i-1]))
        end
    end
    return out
end

function threshold_episodes(t, x; thresh = 0.9)
    starts = Float64[]
    stops = Float64[]
    inside = x[1] >= thresh
    if inside
        push!(starts, t[1])
    end

    for i in 2:length(x)
        if !inside && x[i-1] < thresh <= x[i]
            α = (thresh - x[i-1]) / (x[i] - x[i-1])
            push!(starts, t[i-1] + α * (t[i] - t[i-1]))
            inside = true
        elseif inside && x[i-1] >= thresh > x[i]
            α = (x[i-1] - thresh) / (x[i-1] - x[i])
            push!(stops, t[i-1] + α * (t[i] - t[i-1]))
            inside = false
        end
    end

    if inside
        push!(stops, t[end])
    end

    n = min(length(starts), length(stops))
    starts = starts[1:n]
    stops = stops[1:n]
    return (starts = starts, stops = stops, widths = stops .- starts)
end

function safe_cv(values)
    if length(values) <= 1
        return 0.0
    end
    μ = mean(values)
    return iszero(μ) ? 0.0 : std(values) / μ
end

function lag1corr(values)
    if length(values) <= 2
        return 0.0
    end
    x1 = values[1:end-1]
    x2 = values[2:end]
    sx = std(x1)
    sy = std(x2)
    if iszero(sx) || iszero(sy)
        return 0.0
    end
    return cor(x1, x2)
end

function diagnose_3d_pattern(sim)
    tc_hi = threshold_crossings(sim.t, sim.d; thresh = 0.9)
    tc_mid = threshold_crossings(sim.t, sim.d; thresh = 0.5)
    periods_hi = length(tc_hi) > 1 ? diff(tc_hi) : Float64[]
    periods_mid = length(tc_mid) > 1 ? diff(tc_mid) : Float64[]
    high_frac = sum(sim.d .> 0.9) / length(sim.d)
    mid_frac = sum((sim.d .> 0.2) .& (sim.d .< 0.9)) / length(sim.d)
    low_frac = sum(sim.d .< 0.2) / length(sim.d)
    crossing_ratio = length(tc_mid) / max(length(tc_hi), 1)
    period_mean = isempty(periods_hi) ? 0.0 : mean(periods_hi)
    period_cv = safe_cv(periods_hi)
    subcycle_cv = safe_cv(periods_mid)
    # Regular relaxation cycles typically keep crossing_ratio ≈ 1, tiny period CV,
    # and negligible time in the middle band. MMO-like cases should score higher.
    mmo_score =
        sim.damp *
        (
            0.6 * max(crossing_ratio - 1.0, 0.0) +
            2.0 * period_cv +
            1.5 * mid_frac +
            0.5 * subcycle_cv
        )

    return (
        damp = sim.damp,
        hamp = sim.hamp,
        period_mean = period_mean,
        period_cv = period_cv,
        subcycle_cv = subcycle_cv,
        high_frac = high_frac,
        mid_frac = mid_frac,
        low_frac = low_frac,
        n_hi = length(tc_hi),
        n_mid = length(tc_mid),
        crossing_ratio = crossing_ratio,
        mmo_score = mmo_score,
    )
end

function diagnose_3d_relaxation_regime(E_values, p)
    E = collect(E_values)
    damp = Vector{Float64}(undef, length(E))
    hamp = Vector{Float64}(undef, length(E))
    period_mean = Vector{Float64}(undef, length(E))
    period_cv = Vector{Float64}(undef, length(E))
    crossing_ratio = Vector{Float64}(undef, length(E))
    subcycle_cv = Vector{Float64}(undef, length(E))
    high_frac = Vector{Float64}(undef, length(E))
    mid_frac = Vector{Float64}(undef, length(E))
    low_frac = Vector{Float64}(undef, length(E))

    @threads for idx in eachindex(E)
        par = merge(p, (E = E[idx],))
        sim = simulate_3d(par; u0 = [0.95, 0.5, 0.9], tmax = 7000.0, transient = 3000.0, saveat = 0.1)
        diag = diagnose_3d_pattern(sim)
        damp[idx] = diag.damp
        hamp[idx] = diag.hamp
        period_mean[idx] = diag.period_mean
        period_cv[idx] = diag.period_cv
        crossing_ratio[idx] = diag.crossing_ratio
        subcycle_cv[idx] = diag.subcycle_cv
        high_frac[idx] = diag.high_frac
        mid_frac[idx] = diag.mid_frac
        low_frac[idx] = diag.low_frac
    end

    return (
        E = E,
        damp = damp,
        hamp = hamp,
        period_mean = period_mean,
        period_cv = period_cv,
        crossing_ratio = crossing_ratio,
        subcycle_cv = subcycle_cv,
        high_frac = high_frac,
        mid_frac = mid_frac,
        low_frac = low_frac,
    )
end

function summarize_3d_relaxation_regime(diag; amp_threshold = 0.1)
    mask = diag.damp .> amp_threshold
    if !any(mask)
        return (n_osc = 0, max_mid_frac = 0.0, max_period_cv = 0.0, min_period = 0.0, max_period = 0.0)
    end

    return (
        n_osc = count(mask),
        max_mid_frac = maximum(diag.mid_frac[mask]),
        max_period_cv = maximum(diag.period_cv[mask]),
        min_period = minimum(diag.period_mean[mask]),
        max_period = maximum(diag.period_mean[mask]),
    )
end

function targeted_3d_mmo_search(E_values, epsilon_values, eta_h_values, p)
    E = collect(Float64, E_values)
    epsilon = collect(Float64, epsilon_values)
    eta_h = collect(Float64, eta_h_values)
    dims = (length(eta_h), length(epsilon))
    best_score = Matrix{Float64}(undef, dims...)
    best_E = Matrix{Float64}(undef, dims...)
    best_damp = Matrix{Float64}(undef, dims...)
    best_hamp = Matrix{Float64}(undef, dims...)
    best_period_cv = Matrix{Float64}(undef, dims...)
    best_subcycle_cv = Matrix{Float64}(undef, dims...)
    best_mid_frac = Matrix{Float64}(undef, dims...)
    best_crossing_ratio = Matrix{Float64}(undef, dims...)

    cells = CartesianIndices(best_score)
    @threads for linear_idx in eachindex(best_score)
        iη, iϵ = Tuple(cells[linear_idx])
        local_score = -Inf
        local_E = E[1]
        local_diag = nothing

        for Eval in E
            par = merge(p, (E = Eval, epsilon = epsilon[iϵ], eta_h = eta_h[iη]))
            sim = simulate_3d(par; u0 = [0.95, 0.5, 0.9], tmax = 9000.0, transient = 4000.0, saveat = 0.1)
            diag = diagnose_3d_pattern(sim)
            if diag.mmo_score > local_score
                local_score = diag.mmo_score
                local_E = Eval
                local_diag = diag
            end
        end

        best_score[iη, iϵ] = local_score
        best_E[iη, iϵ] = local_E
        best_damp[iη, iϵ] = local_diag.damp
        best_hamp[iη, iϵ] = local_diag.hamp
        best_period_cv[iη, iϵ] = local_diag.period_cv
        best_subcycle_cv[iη, iϵ] = local_diag.subcycle_cv
        best_mid_frac[iη, iϵ] = local_diag.mid_frac
        best_crossing_ratio[iη, iϵ] = local_diag.crossing_ratio
    end

    best_idx = argmax(best_score)
    best_ieta, best_ieps = Tuple(cells[best_idx])
    return (
        E = E,
        epsilon = epsilon,
        eta_h = eta_h,
        score = best_score,
        best_E = best_E,
        best_damp = best_damp,
        best_hamp = best_hamp,
        best_period_cv = best_period_cv,
        best_subcycle_cv = best_subcycle_cv,
        best_mid_frac = best_mid_frac,
        best_crossing_ratio = best_crossing_ratio,
        best = (
            epsilon = epsilon[best_ieps],
            eta_h = eta_h[best_ieta],
            E = best_E[best_ieta, best_ieps],
            score = best_score[best_ieta, best_ieps],
            damp = best_damp[best_ieta, best_ieps],
            hamp = best_hamp[best_ieta, best_ieps],
            period_cv = best_period_cv[best_ieta, best_ieps],
            subcycle_cv = best_subcycle_cv[best_ieta, best_ieps],
            mid_frac = best_mid_frac[best_ieta, best_ieps],
            crossing_ratio = best_crossing_ratio[best_ieta, best_ieps],
        ),
    )
end

function sweep_4d_grid(E_values, G_values, p)
    E = collect(E_values)
    G = collect(G_values)
    damp = Matrix{Float64}(undef, length(G), length(E))
    qamp = Matrix{Float64}(undef, length(G), length(E))
    kamp = Matrix{Float64}(undef, length(G), length(E))
    hamp = Matrix{Float64}(undef, length(G), length(E))

    @threads for linear_idx in eachindex(damp)
        iG, iE = Tuple(CartesianIndices(damp)[linear_idx])
        par = merge(p, (E = E[iE], G = G[iG]))
        sim = simulate_full(par; u0 = [0.05, 0.0, 0.1, 0.6], tmax = 16000.0, transient = 8000.0, saveat = 0.5)
        damp[iG, iE] = sim.damp
        qamp[iG, iE] = sim.qamp
        kamp[iG, iE] = maximum(sim.k) - minimum(sim.k)
        hamp[iG, iE] = maximum(sim.h) - minimum(sim.h)
    end

    return (E = E, G = G, damp = damp, qamp = qamp, kamp = kamp, hamp = hamp)
end

function diagnose_4d_candidate(p; u0 = [0.05, 0.0, 0.1, 0.6])
    sim = simulate_full(p; u0 = u0, tmax = 18000.0, transient = 9000.0, saveat = 0.2)
    diag = diagnose_4d_pattern(sim)
    alt = diagnose_4d_alternation(sim)
    return (
        sim = sim,
        period_mean = diag.period_mean,
        period_cv = diag.period_cv,
        burst_width_mean = diag.burst_width_mean,
        burst_width_cv = diag.burst_width_cv,
        high_frac = diag.high_frac,
        mid_frac = diag.mid_frac,
        low_frac = diag.low_frac,
        kamp = diag.kamp,
        hamp = diag.hamp,
        crossing_ratio = diag.crossing_ratio,
        burst_score = diag.burst_score,
        n_events = alt.n,
        odd_period_mean = alt.odd_period_mean,
        even_period_mean = alt.even_period_mean,
        odd_width_mean = alt.odd_width_mean,
        even_width_mean = alt.even_width_mean,
        period_lag1 = alt.period_lag1,
        width_lag1 = alt.width_lag1,
        period_alt_gap = alt.period_alt_gap,
        width_alt_gap = alt.width_alt_gap,
        periods = alt.periods,
        widths = alt.widths,
    )
end

function diagnose_4d_pattern(sim)
    tc_hi = threshold_crossings(sim.t, sim.d; thresh = 0.9)
    tc_mid = threshold_crossings(sim.t, sim.d; thresh = 0.5)
    periods_hi = length(tc_hi) > 1 ? diff(tc_hi) : Float64[]
    episodes = threshold_episodes(sim.t, sim.d; thresh = 0.9)
    high_frac = sum(sim.d .> 0.9) / length(sim.d)
    mid_frac = sum((sim.d .> 0.2) .& (sim.d .< 0.9)) / length(sim.d)
    low_frac = sum(sim.d .< 0.2) / length(sim.d)
    kamp = maximum(sim.k) - minimum(sim.k)
    hamp = maximum(sim.h) - minimum(sim.h)
    crossing_ratio = length(tc_mid) / max(length(tc_hi), 1)
    period_mean = isempty(periods_hi) ? 0.0 : mean(periods_hi)
    period_cv = safe_cv(periods_hi)
    burst_width_mean = isempty(episodes.widths) ? 0.0 : mean(episodes.widths)
    burst_width_cv = safe_cv(episodes.widths)
    burst_score =
        sim.damp *
        (
            2.0 * period_cv +
            1.5 * burst_width_cv +
            1.2 * max(crossing_ratio - 1.0, 0.0) +
            1.5 * mid_frac +
            0.3 * (kamp + hamp)
        )

    return (
        damp = sim.damp,
        qamp = sim.qamp,
        kamp = kamp,
        hamp = hamp,
        period_mean = period_mean,
        period_cv = period_cv,
        burst_width_mean = burst_width_mean,
        burst_width_cv = burst_width_cv,
        high_frac = high_frac,
        mid_frac = mid_frac,
        low_frac = low_frac,
        crossing_ratio = crossing_ratio,
        burst_score = burst_score,
    )
end

function diagnose_4d_alternation(sim)
    tc_hi = threshold_crossings(sim.t, sim.d; thresh = 0.9)
    periods = length(tc_hi) > 1 ? diff(tc_hi) : Float64[]
    episodes = threshold_episodes(sim.t, sim.d; thresh = 0.9)
    widths = copy(episodes.widths)
    n = min(length(periods), length(widths))
    periods = periods[1:n]
    widths = widths[1:n]

    odd_periods = periods[1:2:end]
    even_periods = periods[2:2:end]
    odd_widths = widths[1:2:end]
    even_widths = widths[2:2:end]
    odd_period_mean = isempty(odd_periods) ? 0.0 : mean(odd_periods)
    even_period_mean = isempty(even_periods) ? 0.0 : mean(even_periods)
    odd_width_mean = isempty(odd_widths) ? 0.0 : mean(odd_widths)
    even_width_mean = isempty(even_widths) ? 0.0 : mean(even_widths)
    period_alt_gap = iszero(mean(periods)) ? 0.0 : abs(even_period_mean - odd_period_mean) / mean(periods)
    width_alt_gap = iszero(mean(widths)) ? 0.0 : abs(even_width_mean - odd_width_mean) / mean(widths)

    return (
        n = n,
        periods = periods,
        widths = widths,
        odd_period_mean = odd_period_mean,
        even_period_mean = even_period_mean,
        odd_width_mean = odd_width_mean,
        even_width_mean = even_width_mean,
        period_lag1 = lag1corr(periods),
        width_lag1 = lag1corr(widths),
        period_alt_gap = period_alt_gap,
        width_alt_gap = width_alt_gap,
    )
end

function targeted_4d_burst_search(E_values, G_values, eta_h_values, eta_k_values, p)
    E = collect(Float64, E_values)
    G = collect(Float64, G_values)
    eta_h = collect(Float64, eta_h_values)
    eta_k = collect(Float64, eta_k_values)
    dims = (length(eta_h), length(eta_k))
    best_score = Matrix{Float64}(undef, dims...)
    best_E = Matrix{Float64}(undef, dims...)
    best_G = Matrix{Float64}(undef, dims...)
    best_damp = Matrix{Float64}(undef, dims...)
    best_kamp = Matrix{Float64}(undef, dims...)
    best_hamp = Matrix{Float64}(undef, dims...)
    best_period_cv = Matrix{Float64}(undef, dims...)
    best_width_cv = Matrix{Float64}(undef, dims...)
    best_mid_frac = Matrix{Float64}(undef, dims...)
    best_crossing_ratio = Matrix{Float64}(undef, dims...)

    cells = CartesianIndices(best_score)
    @threads for linear_idx in eachindex(best_score)
        iηh, iηk = Tuple(cells[linear_idx])
        local_score = -Inf
        local_E = E[1]
        local_G = G[1]
        local_diag = nothing

        for Gval in G, Eval in E
            par = merge(p, (E = Eval, G = Gval, eta_h = eta_h[iηh], eta_k = eta_k[iηk]))
            sim = simulate_full(par; u0 = [0.05, 0.0, 0.1, 0.6], tmax = 18000.0, transient = 9000.0, saveat = 0.2)
            diag = diagnose_4d_pattern(sim)
            if diag.burst_score > local_score
                local_score = diag.burst_score
                local_E = Eval
                local_G = Gval
                local_diag = diag
            end
        end

        best_score[iηh, iηk] = local_score
        best_E[iηh, iηk] = local_E
        best_G[iηh, iηk] = local_G
        best_damp[iηh, iηk] = local_diag.damp
        best_kamp[iηh, iηk] = local_diag.kamp
        best_hamp[iηh, iηk] = local_diag.hamp
        best_period_cv[iηh, iηk] = local_diag.period_cv
        best_width_cv[iηh, iηk] = local_diag.burst_width_cv
        best_mid_frac[iηh, iηk] = local_diag.mid_frac
        best_crossing_ratio[iηh, iηk] = local_diag.crossing_ratio
    end

    best_idx = argmax(best_score)
    best_iηh, best_iηk = Tuple(cells[best_idx])
    return (
        E = E,
        G = G,
        eta_h = eta_h,
        eta_k = eta_k,
        score = best_score,
        best_E = best_E,
        best_G = best_G,
        best_damp = best_damp,
        best_kamp = best_kamp,
        best_hamp = best_hamp,
        best_period_cv = best_period_cv,
        best_width_cv = best_width_cv,
        best_mid_frac = best_mid_frac,
        best_crossing_ratio = best_crossing_ratio,
        best = (
            eta_h = eta_h[best_iηh],
            eta_k = eta_k[best_iηk],
            E = best_E[best_iηh, best_iηk],
            G = best_G[best_iηh, best_iηk],
            score = best_score[best_iηh, best_iηk],
            damp = best_damp[best_iηh, best_iηk],
            kamp = best_kamp[best_iηh, best_iηk],
            hamp = best_hamp[best_iηh, best_iηk],
            period_cv = best_period_cv[best_iηh, best_iηk],
            burst_width_cv = best_width_cv[best_iηh, best_iηk],
            mid_frac = best_mid_frac[best_iηh, best_iηk],
            crossing_ratio = best_crossing_ratio[best_iηh, best_iηk],
        ),
    )
end

function targeted_4d_geometry_search(theta0_values, beta_n_values, E_values, p)
    theta0 = collect(Float64, theta0_values)
    beta_n = collect(Float64, beta_n_values)
    E = collect(Float64, E_values)
    dims = (length(theta0), length(beta_n))
    best_score = Matrix{Float64}(undef, dims...)
    best_E = Matrix{Float64}(undef, dims...)
    best_damp = Matrix{Float64}(undef, dims...)
    best_kamp = Matrix{Float64}(undef, dims...)
    best_hamp = Matrix{Float64}(undef, dims...)
    best_period_cv = Matrix{Float64}(undef, dims...)
    best_width_cv = Matrix{Float64}(undef, dims...)
    best_mid_frac = Matrix{Float64}(undef, dims...)
    best_crossing_ratio = Matrix{Float64}(undef, dims...)

    cells = CartesianIndices(best_score)
    @threads for linear_idx in eachindex(best_score)
        iθ, iβ = Tuple(cells[linear_idx])
        local_score = -Inf
        local_E = E[1]
        local_diag = nothing

        for Eval in E
            par = merge(p, (theta0 = theta0[iθ], beta_n = beta_n[iβ], E = Eval))
            sim = simulate_full(par; u0 = [0.05, 0.0, 0.1, 0.6], tmax = 18000.0, transient = 9000.0, saveat = 0.2)
            diag = diagnose_4d_pattern(sim)
            if diag.burst_score > local_score
                local_score = diag.burst_score
                local_E = Eval
                local_diag = diag
            end
        end

        best_score[iθ, iβ] = local_score
        best_E[iθ, iβ] = local_E
        best_damp[iθ, iβ] = local_diag.damp
        best_kamp[iθ, iβ] = local_diag.kamp
        best_hamp[iθ, iβ] = local_diag.hamp
        best_period_cv[iθ, iβ] = local_diag.period_cv
        best_width_cv[iθ, iβ] = local_diag.burst_width_cv
        best_mid_frac[iθ, iβ] = local_diag.mid_frac
        best_crossing_ratio[iθ, iβ] = local_diag.crossing_ratio
    end

    best_idx = argmax(best_score)
    best_iθ, best_iβ = Tuple(cells[best_idx])
    return (
        theta0 = theta0,
        beta_n = beta_n,
        E = E,
        score = best_score,
        best_E = best_E,
        best_damp = best_damp,
        best_kamp = best_kamp,
        best_hamp = best_hamp,
        best_period_cv = best_period_cv,
        best_width_cv = best_width_cv,
        best_mid_frac = best_mid_frac,
        best_crossing_ratio = best_crossing_ratio,
        best = (
            theta0 = theta0[best_iθ],
            beta_n = beta_n[best_iβ],
            E = best_E[best_iθ, best_iβ],
            score = best_score[best_iθ, best_iβ],
            damp = best_damp[best_iθ, best_iβ],
            kamp = best_kamp[best_iθ, best_iβ],
            hamp = best_hamp[best_iθ, best_iβ],
            period_cv = best_period_cv[best_iθ, best_iβ],
            burst_width_cv = best_width_cv[best_iθ, best_iβ],
            mid_frac = best_mid_frac[best_iθ, best_iβ],
            crossing_ratio = best_crossing_ratio[best_iθ, best_iβ],
        ),
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

function save_branch_plot(br, po_data)
    if isnothing(po_data)
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
            br;
            xlabel = "E",
            ylabel = "d",
            title = "Equilibria and periodic orbits in the 2D reduced model",
            markersize = 3,
            legend = :bottomleft,
        )
        plot!(p, po_data.param, po_data.dmin; lw = 2, label = "periodic dmin")
        plot!(p, po_data.param, po_data.dmax; lw = 2, label = "periodic dmax")
    end
    savefig(p, joinpath(RESULTS_DIR, "double_frontier_branch.png"))
end

function save_hopf_neighborhood_plot(br, hopf_data, po_data)
    if isnothing(hopf_data) || isnothing(po_data)
        return
    end

    hopf_E = hopf_data.point.p
    eq_d = hopf_data.point.x0[1]
    yvals = vcat(po_data.dmin, po_data.dmax, [eq_d])
    ypad = max(5e-4, 0.15 * (maximum(yvals) - minimum(yvals)))

    p = plot(
        br;
        xlabel = "E",
        ylabel = "d",
        title = "Hopf-neighborhood periodic-orbit continuation",
        markersize = 3,
        legend = :bottomleft,
    )
    plot!(p, po_data.param, po_data.dmin; lw = 2, label = "periodic dmin")
    plot!(p, po_data.param, po_data.dmax; lw = 2, label = "periodic dmax")
    vline!(p, [hopf_E]; lw = 2, ls = :dash, color = :black, label = @sprintf("refined Hopf %.6f", hopf_E))
    xlims!(p, (min(minimum(po_data.param), hopf_E) - 5e-4, max(maximum(po_data.param), hopf_E) + 5e-4))
    ylims!(p, (minimum(yvals) - ypad, maximum(yvals) + ypad))
    savefig(p, joinpath(RESULTS_DIR, "double_frontier_hopf_neighborhood.png"))
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

function save_3d_branch_plot(br)
    p1 = plot(br.param, br.d; lw = 2, marker = :circle, xlabel = "E", ylabel = "d", title = "3D equilibrium continuation", label = "d")
    p2 = plot(br.param, br.q; lw = 2, marker = :circle, xlabel = "E", ylabel = "q", title = "3D equilibrium continuation", label = "q")
    p3 = plot(br.param, br.h; lw = 2, marker = :circle, xlabel = "E", ylabel = "h", title = "3D equilibrium continuation", label = "h")
    combo = plot(p1, p2, p3; layout = (3, 1), size = (900, 1100))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_branch.png"))
end

function save_3d_hysteresis_plot(sweep)
    p1 = plot(sweep.E, sweep.low_d; lw = 2, marker = :circle, xlabel = "E", ylabel = "mean d", title = "3D low vs high initial conditions", label = "low init")
    plot!(p1, sweep.E, sweep.high_d; lw = 2, marker = :square, label = "high init")

    p2 = plot(sweep.E, sweep.low_h; lw = 2, marker = :circle, xlabel = "E", ylabel = "mean h", title = "3D judgment-capital response", label = "low init")
    plot!(p2, sweep.E, sweep.high_h; lw = 2, marker = :square, label = "high init")

    p3 = plot(sweep.E, sweep.gap_d; lw = 2, marker = :circle, xlabel = "E", ylabel = "mean-gap", title = "3D path dependence", label = "d gap")
    plot!(p3, sweep.E, sweep.gap_h; lw = 2, marker = :square, label = "h gap")

    p4 = plot(sweep.E, sweep.amp_low; lw = 2, marker = :circle, xlabel = "E", ylabel = "amplitude", title = "3D oscillation amplitudes", label = "d amp low")
    plot!(p4, sweep.E, sweep.amp_high; lw = 2, marker = :square, label = "d amp high")
    plot!(p4, sweep.E, sweep.hamp_low; lw = 2, marker = :diamond, label = "h amp low")
    plot!(p4, sweep.E, sweep.hamp_high; lw = 2, marker = :utriangle, label = "h amp high")

    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_hysteresis.png"))
end

function save_3d_hopf_zoom_plot(scan, hopf_E)
    p1 = plot(scan.E, scan.low_damp; lw = 2, marker = :circle, xlabel = "E", ylabel = "amplitude", title = "3D Hopf neighborhood amplitudes", label = "d amp low")
    plot!(p1, scan.E, scan.high_damp; lw = 2, marker = :square, label = "d amp high")
    plot!(p1, scan.E, scan.low_hamp; lw = 2, marker = :diamond, label = "h amp low")
    plot!(p1, scan.E, scan.high_hamp; lw = 2, marker = :utriangle, label = "h amp high")
    vline!(p1, [hopf_E]; lw = 2, ls = :dash, color = :black, label = @sprintf("3D Hopf %.6f", hopf_E))

    p2 = plot(scan.E, scan.low_dmean; lw = 2, marker = :circle, xlabel = "E", ylabel = "mean state", title = "3D Hopf neighborhood means", label = "mean d low")
    plot!(p2, scan.E, scan.high_dmean; lw = 2, marker = :square, label = "mean d high")
    plot!(p2, scan.E, scan.low_hmean; lw = 2, marker = :diamond, label = "mean h low")
    plot!(p2, scan.E, scan.high_hmean; lw = 2, marker = :utriangle, label = "mean h high")
    vline!(p2, [hopf_E]; lw = 2, ls = :dash, color = :black, label = @sprintf("3D Hopf %.6f", hopf_E))

    combo = plot(p1, p2; layout = (2, 1), size = (1000, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_hopf_zoom.png"))
end

function save_3d_regime_plot(diag)
    p1 = plot(diag.E, diag.period_mean; lw = 2, marker = :circle, xlabel = "E", ylabel = "return time", title = "3D oscillation return time", label = "mean cycle time")
    p2 = plot(diag.E, diag.period_cv; lw = 2, marker = :circle, xlabel = "E", ylabel = "CV", title = "3D return-time variability", label = "period CV")
    p3 = plot(diag.E, diag.high_frac; lw = 2, marker = :circle, xlabel = "E", ylabel = "occupancy fraction", title = "3D occupancy by d-band", label = "high (d > 0.9)")
    plot!(p3, diag.E, diag.mid_frac; lw = 2, marker = :square, label = "mid (0.2 < d < 0.9)")
    plot!(p3, diag.E, diag.low_frac; lw = 2, marker = :diamond, label = "low (d < 0.2)")
    p4 = plot(diag.E, diag.damp; lw = 2, marker = :circle, xlabel = "E", ylabel = "signature", title = "3D oscillation amplitude and subcycle ratio", label = "d amp")
    plot!(p4, diag.E, diag.hamp; lw = 2, marker = :square, label = "h amp")
    plot!(p4, diag.E, diag.crossing_ratio; lw = 2, marker = :diamond, label = "0.5/0.9 crossing ratio")
    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_regime.png"))
end

function save_3d_targeted_search_plot(search)
    eps_perm = sortperm(search.epsilon)
    eta_perm = sortperm(search.eta_h)
    epsilon = search.epsilon[eps_perm]
    eta_h = search.eta_h[eta_perm]
    reorder(mat) = mat[eta_perm, eps_perm]

    p1 = heatmap(epsilon, eta_h, reorder(search.score); xlabel = "epsilon", ylabel = "eta_h", title = "3D targeted MMO score", colorbar_title = "score")
    p2 = heatmap(epsilon, eta_h, reorder(search.best_E); xlabel = "epsilon", ylabel = "eta_h", title = "Best E in targeted search", colorbar_title = "E")
    p3 = heatmap(epsilon, eta_h, reorder(search.best_period_cv); xlabel = "epsilon", ylabel = "eta_h", title = "Best period CV", colorbar_title = "CV")
    p4 = heatmap(epsilon, eta_h, reorder(search.best_mid_frac); xlabel = "epsilon", ylabel = "eta_h", title = "Best mid-band occupancy", colorbar_title = "mid frac")
    p5 = heatmap(epsilon, eta_h, reorder(search.best_crossing_ratio); xlabel = "epsilon", ylabel = "eta_h", title = "Best 0.5/0.9 crossing ratio", colorbar_title = "ratio")
    p6 = heatmap(epsilon, eta_h, reorder(search.best_damp); xlabel = "epsilon", ylabel = "eta_h", title = "Best d amplitude", colorbar_title = "d amp")
    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 1200))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_targeted_search.png"))
end

function save_4d_regime_plot(sweep)
    p1 = heatmap(sweep.E, sweep.G, sweep.damp; xlabel = "E", ylabel = "G", title = "4D sustained d amplitude", colorbar_title = "d amp")
    p2 = heatmap(sweep.E, sweep.G, sweep.qamp; xlabel = "E", ylabel = "G", title = "4D sustained q amplitude", colorbar_title = "q amp")
    p3 = heatmap(sweep.E, sweep.G, sweep.kamp; xlabel = "E", ylabel = "G", title = "4D sustained k amplitude", colorbar_title = "k amp")
    p4 = heatmap(sweep.E, sweep.G, sweep.hamp; xlabel = "E", ylabel = "G", title = "4D sustained h amplitude", colorbar_title = "h amp")
    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_4d_regime.png"))
end

function save_4d_targeted_search_plot(search)
    ηh_perm = sortperm(search.eta_h)
    ηk_perm = sortperm(search.eta_k)
    eta_h = search.eta_h[ηh_perm]
    eta_k = search.eta_k[ηk_perm]
    reorder(mat) = mat[ηh_perm, ηk_perm]

    p1 = heatmap(eta_k, eta_h, reorder(search.score); xlabel = "eta_k", ylabel = "eta_h", title = "4D targeted burst score", colorbar_title = "score")
    p2 = heatmap(eta_k, eta_h, reorder(search.best_E); xlabel = "eta_k", ylabel = "eta_h", title = "Best E in 4D targeted search", colorbar_title = "E")
    p3 = heatmap(eta_k, eta_h, reorder(search.best_G); xlabel = "eta_k", ylabel = "eta_h", title = "Best G in 4D targeted search", colorbar_title = "G")
    p4 = heatmap(eta_k, eta_h, reorder(search.best_period_cv); xlabel = "eta_k", ylabel = "eta_h", title = "Best inter-burst CV", colorbar_title = "CV")
    p5 = heatmap(eta_k, eta_h, reorder(search.best_width_cv); xlabel = "eta_k", ylabel = "eta_h", title = "Best burst-width CV", colorbar_title = "CV")
    p6 = heatmap(eta_k, eta_h, reorder(search.best_mid_frac); xlabel = "eta_k", ylabel = "eta_h", title = "Best mid-band occupancy", colorbar_title = "mid frac")
    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 1200))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_4d_targeted_search.png"))
end

function save_4d_geometry_search_plot(search)
    θ_perm = sortperm(search.theta0)
    β_perm = sortperm(search.beta_n)
    theta0 = search.theta0[θ_perm]
    beta_n = search.beta_n[β_perm]
    reorder(mat) = mat[θ_perm, β_perm]

    p1 = heatmap(beta_n, theta0, reorder(search.score); xlabel = "beta_n", ylabel = "theta0", title = "4D geometry-search burst score", colorbar_title = "score")
    p2 = heatmap(beta_n, theta0, reorder(search.best_E); xlabel = "beta_n", ylabel = "theta0", title = "Best E in geometry search", colorbar_title = "E")
    p3 = heatmap(beta_n, theta0, reorder(search.best_period_cv); xlabel = "beta_n", ylabel = "theta0", title = "Best inter-burst CV", colorbar_title = "CV")
    p4 = heatmap(beta_n, theta0, reorder(search.best_width_cv); xlabel = "beta_n", ylabel = "theta0", title = "Best burst-width CV", colorbar_title = "CV")
    p5 = heatmap(beta_n, theta0, reorder(search.best_mid_frac); xlabel = "beta_n", ylabel = "theta0", title = "Best mid-band occupancy", colorbar_title = "mid frac")
    p6 = heatmap(beta_n, theta0, reorder(search.best_crossing_ratio); xlabel = "beta_n", ylabel = "theta0", title = "Best 0.5/0.9 crossing ratio", colorbar_title = "ratio")
    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 1200))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_4d_geometry_search.png"))
end

function save_4d_modulated_candidate_plot(sim)
    p1 = plot(sim.t, sim.d; lw = 2, xlabel = "t", ylabel = "d(t)", title = "4D modulated candidate: doer share", label = "d")
    p2 = plot(sim.t, sim.q; lw = 2, xlabel = "t", ylabel = "q(t)", title = "4D modulated candidate: backlog", label = "q")
    p3 = plot(sim.t, sim.k; lw = 2, xlabel = "t", ylabel = "k(t)", title = "4D modulated candidate: adoption capital", label = "k")
    p4 = plot(sim.t, sim.h; lw = 2, xlabel = "t", ylabel = "h(t)", title = "4D modulated candidate: judgment capital", label = "h")
    p5 = plot(sim.t, sim.m; lw = 2, xlabel = "t", ylabel = "m(t)", title = "4D modulated candidate: middle band", label = "m")
    p6 = plot(sim.k, sim.h; lw = 2, xlabel = "k", ylabel = "h", title = "4D modulated candidate: slow drift", label = "trajectory")
    combo = plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1200, 900))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_4d_modulated_candidate.png"))
end

function save_4d_alternation_plot(diag)
    nshow = min(40, length(diag.periods))
    idx = collect(1:nshow)
    periods = diag.periods[1:nshow]
    widths = diag.widths[1:nshow]

    p1 = plot(idx, periods; lw = 2, marker = :circle, xlabel = "event index", ylabel = "inter-burst period", title = "4D modulated candidate: period sequence", label = "period")
    p2 = plot(idx, widths; lw = 2, marker = :circle, xlabel = "event index", ylabel = "high-state width", title = "4D modulated candidate: burst-width sequence", label = "width")
    p3 = scatter(periods[1:end-1], periods[2:end]; xlabel = "period_n", ylabel = "period_{n+1}", title = "4D modulated candidate: period return map", label = "return map")
    p4 = scatter(widths[1:end-1], widths[2:end]; xlabel = "width_n", ylabel = "width_{n+1}", title = "4D modulated candidate: width return map", label = "return map")
    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_4d_alternation.png"))
end

function save_3d_timeseries_plot(sim)
    p1 = plot(sim.t, sim.d; lw = 2, xlabel = "t", ylabel = "d(t)", title = "3D doer share", label = "d")
    p2 = plot(sim.t, sim.q; lw = 2, xlabel = "t", ylabel = "q(t)", title = "3D supervision backlog", label = "q")
    p3 = plot(sim.t, sim.h; lw = 2, xlabel = "t", ylabel = "h(t)", title = "3D judgment capital", label = "h")
    p4 = plot(sim.d, sim.h; lw = 2, xlabel = "d", ylabel = "h", title = "3D trajectory projection", label = "trajectory")
    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_3d_timeseries.png"))
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

function save_summary(br, po_data, sweep, rep_sim, br3d, sweep3d, hopf3d_scan, hopf3d_transition, regime3d_diag, regime3d_summary, targeted3d_search, sim3d, sweep4d, targeted4d_search, geometry4d_search, full_diag, modulated4d_diag, p, hopf_idx, fold_idx, hopf_data, po_result, hopf3d_idx, fold3d_idx)
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

        if isnothing(po_data)
            println(io, "Periodic-orbit continuation: no branch was successfully continued from the refined Hopf point")
        else
            println(io, @sprintf("Periodic-orbit continuation points saved = %d", length(po_data.param)))
            println(io, @sprintf("Periodic-orbit parameter span: E in [%.10f, %.10f]", minimum(po_data.param), maximum(po_data.param)))
            println(io, @sprintf("Periodic-orbit d-amplitude span: [%.6e, %.6e]", minimum(po_data.damp), maximum(po_data.damp)))
            pd_points = dedupe_specialpoints(filter(sp -> sp.type == :pd, po_data.specialpoints))
            if !isempty(pd_points)
                println(io, @sprintf("Detected %d unique period-doubling point(s) on the extended periodic branch", length(pd_points)))
                for sp in pd_points
                    println(io, @sprintf("- pd near E ≈ %.10f [idx = %d, status = %s]", sp.param, sp.idx, String(sp.status)))
                end
            end
        end

        if !isnothing(hopf_data)
            hopf = hopf_data.point
            println(io)
            println(io, @sprintf("Refined Hopf point from minimally augmented Newton: E ≈ %.10f", hopf.p))
            println(io, @sprintf("Refined Hopf frequency: |ω| ≈ %.10f, period ≈ %.10f", abs(hopf.ω), 2pi / abs(hopf.ω)))
            println(io, @sprintf("Hopf normal-form type: %s", String(hopf.type)))
            println(io, @sprintf("Hopf normal-form coefficient a ≈ %.6e %+.6ei", real(hopf.nf.a), imag(hopf.nf.a)))
            println(io, @sprintf("Hopf normal-form coefficient b ≈ %.6e %+.6ei", real(hopf.nf.b), imag(hopf.nf.b)))
        end

        if !isnothing(po_result) && !isempty(po_result.failures)
            println(io)
            println(io, "Additional periodic-orbit attempts that failed:")
            for failure in po_result.failures
                println(io, "- ", failure.name, ": ", failure.error)
            end
        end

        println(io)
        println(io, @sprintf("Largest low/high branch gap in direct sweeps = %.6f", maximum(sweep.gap)))
        println(io, @sprintf("Representative 2D amplitude at E = %.3f: d_amp = %.6f, q_amp = %.6f", rep_sim.E, rep_sim.damp, rep_sim.qamp))
        println(io)
        if !isnothing(fold3d_idx)
            println(io, @sprintf("3D fold/branch-point candidate (%s) near E ≈ %.10f", String(br3d.specialpoint[fold3d_idx].type), br3d.specialpoint[fold3d_idx].param))
        else
            println(io, "3D fold detected near E ≈ none on the current equilibrium branch")
        end
        if !isnothing(hopf3d_idx)
            println(io, @sprintf("3D Hopf near E ≈ %.10f", br3d.specialpoint[hopf3d_idx].param))
        else
            println(io, "3D Hopf near E ≈ none on the current equilibrium branch")
        end
        println(io, @sprintf("Largest 3D low/high mean-d gap in direct sweeps = %.6f", maximum(sweep3d.gap_d)))
        println(io, @sprintf("Largest 3D low/high mean-h gap in direct sweeps = %.6f", maximum(sweep3d.gap_h)))
        if !isnothing(hopf3d_transition.last_osc_E) && !isnothing(hopf3d_transition.first_quiet_E)
            println(io, @sprintf("3D Hopf-neighborhood oscillation-collapse window: E ≈ %.10f -> %.10f", hopf3d_transition.last_osc_E, hopf3d_transition.first_quiet_E))
        end
        if regime3d_summary.n_osc > 0
            println(io, @sprintf("3D oscillatory-regime return time span: [%.6f, %.6f]", regime3d_summary.min_period, regime3d_summary.max_period))
            println(io, @sprintf("3D oscillatory-regime max period CV = %.6f", regime3d_summary.max_period_cv))
            println(io, @sprintf("3D oscillatory-regime max mid-band occupancy = %.6f", regime3d_summary.max_mid_frac))
        end
        println(io, @sprintf("Best targeted 3D MMO-search candidate: epsilon = %.4f, eta_h = %.4f, E = %.6f, score = %.6f, d_amp = %.6f, h_amp = %.6f, period CV = %.6f, subcycle CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
            targeted3d_search.best.epsilon,
            targeted3d_search.best.eta_h,
            targeted3d_search.best.E,
            targeted3d_search.best.score,
            targeted3d_search.best.damp,
            targeted3d_search.best.hamp,
            targeted3d_search.best.period_cv,
            targeted3d_search.best.subcycle_cv,
            targeted3d_search.best.mid_frac,
            targeted3d_search.best.crossing_ratio,
        ))
        println(io, @sprintf("Representative 3D amplitude at E = %.3f: d_amp = %.6f, q_amp = %.6f, h_amp = %.6f", sim3d.E, sim3d.damp, sim3d.qamp, sim3d.hamp))
        best_idx = argmax(sweep4d.damp)
        best_G, best_E = Tuple(CartesianIndices(sweep4d.damp)[best_idx])
        println(io, @sprintf("Best sustained 4D candidate at E = %.3f, G = %.3f: d_amp = %.6f, q_amp = %.6f, k_amp = %.6f, h_amp = %.6f",
            sweep4d.E[best_E], sweep4d.G[best_G], full_diag.sim.damp, full_diag.sim.qamp, full_diag.kamp, full_diag.hamp))
        println(io, @sprintf("4D candidate return-time CV = %.6f, middle-band occupancy = %.6f", full_diag.period_cv, full_diag.mid_frac))
        println(io, @sprintf("Best targeted 4D burst-search candidate: eta_h = %.4f, eta_k = %.4f, E = %.3f, G = %.3f, score = %.6f, d_amp = %.6f, k_amp = %.6f, h_amp = %.6f, inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
            targeted4d_search.best.eta_h,
            targeted4d_search.best.eta_k,
            targeted4d_search.best.E,
            targeted4d_search.best.G,
            targeted4d_search.best.score,
            targeted4d_search.best.damp,
            targeted4d_search.best.kamp,
            targeted4d_search.best.hamp,
            targeted4d_search.best.period_cv,
            targeted4d_search.best.burst_width_cv,
            targeted4d_search.best.mid_frac,
            targeted4d_search.best.crossing_ratio,
        ))
        println(io, @sprintf("Best 4D geometry-search candidate: theta0 = %.3f, beta_n = %.3f, E = %.3f, score = %.6f, d_amp = %.6f, k_amp = %.6f, h_amp = %.6f, inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
            geometry4d_search.best.theta0,
            geometry4d_search.best.beta_n,
            geometry4d_search.best.E,
            geometry4d_search.best.score,
            geometry4d_search.best.damp,
            geometry4d_search.best.kamp,
            geometry4d_search.best.hamp,
            geometry4d_search.best.period_cv,
            geometry4d_search.best.burst_width_cv,
            geometry4d_search.best.mid_frac,
            geometry4d_search.best.crossing_ratio,
        ))
        println(io, @sprintf("Long-window diagnostic at the geometry-search candidate: inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
            modulated4d_diag.period_cv,
            modulated4d_diag.burst_width_cv,
            modulated4d_diag.mid_frac,
            modulated4d_diag.crossing_ratio,
        ))
        println(io, @sprintf("4D alternation diagnostic: odd/even period means = [%.6f, %.6f], odd/even width means = [%.6f, %.6f], lag-1 correlations = [period %.6f, width %.6f], normalized alternation gaps = [period %.6f, width %.6f]",
            modulated4d_diag.odd_period_mean,
            modulated4d_diag.even_period_mean,
            modulated4d_diag.odd_width_mean,
            modulated4d_diag.even_width_mean,
            modulated4d_diag.period_lag1,
            modulated4d_diag.width_lag1,
            modulated4d_diag.period_alt_gap,
            modulated4d_diag.width_alt_gap,
        ))
        println(io)
        println(io, "Interpretation:")
        println(io, "- The script keeps the full 4D model available, but starts with the 2D reduced system for fold/Hopf screening.")
        println(io, "- Theoretical fold viability is controlled by beta_n * lambda_D > 4, matching the analytic pre-screen in the plan.")
        println(io, "- Direct low/high-initial-condition sweeps provide a Julia-native hysteresis check before moving to deeper MMO or bursting studies.")
        println(io, "- A refined Hopf solve plus multi-stage two-point periodic-orbit continuation now resolves a long subcritical branch beyond the immediate Hopf neighborhood.")
        println(io, "- The extended branch reaches O(1) d-amplitude and encounters a period-doubling point, so the 2D story is no longer just 'Hopf then black-box jump': the large-cycle regime is connected to a concrete periodic-orbit backbone.")
        println(io, "- The new 3D workflow checks whether that 2D geometry survives once judgment capital h becomes dynamic, using equilibrium continuation plus direct low/high-initial-condition scans.")
        if !isnothing(hopf3d_transition.last_osc_E) && !isnothing(hopf3d_transition.first_quiet_E)
            println(io, "- In the 3D Hopf neighborhood, large-amplitude oscillations persist below the Hopf point and collapse before/around it, indicating a strongly subcritical transition with a coexistence window rather than a gentle small-cycle onset.")
        end
        if regime3d_summary.n_osc > 0
            println(io, "- Across the sampled 3D oscillatory regime, threshold-crossing return times stay regular and the orbit spends almost no time in the middle band 0.2 < d < 0.9, which is more consistent with a single relaxation cycle than with MMO-style mixed small/large oscillations.")
        end
        if targeted3d_search.best.score < 0.15
            println(io, "- A targeted 3D search over smaller epsilon and slower eta_h still fails to produce a strong MMO signature: the best candidate keeps low return-time variability, near-unit threshold-crossing ratio, and little middle-band dwell time.")
        else
            println(io, "- The targeted 3D search identifies a parameter patch with elevated MMO score; that region should be treated as the next high-priority follow-up for longer integrations and continuation.")
        end
        println(io, "- In 4D, sustained oscillations survive after long transients mainly at lower G; the earlier large-amplitude example at G = 1.1 was a transient, not a settled bursting state.")
        println(io, "- The best sustained 4D candidate found here still has low return-time variability and low middle-band occupancy, so it also looks like a regular relaxation cycle rather than bursting.")
        if targeted4d_search.best.score < 0.20
            println(io, "- A broader 4D targeted search over eta_h, eta_k, E, and G also fails to uncover a convincing bursting regime: the strongest candidates still have near-unit crossing ratio, tiny middle-band occupancy, and modest timing variability.")
        else
            println(io, "- The 4D targeted search surfaces a parameter patch with elevated burst score; that patch should be the next priority for longer integrations and slow-fast analysis.")
        end
        if geometry4d_search.best.score > 1.0
            println(io, "- A geometry-focused 4D search over theta0 and beta_n reveals a much stronger burst-like modulation regime than the eta_h/eta_k sweep alone, with large inter-cycle and burst-width variability concentrated near lower beta_n and more negative theta0.")
        end
        if modulated4d_diag.crossing_ratio ≈ 1.0 && modulated4d_diag.mid_frac < 0.01
            println(io, "- Even at that stronger 4D modulation candidate, the orbit still behaves like one large excursion per cycle rather than mixed small/large oscillations, so the current evidence supports irregular envelope modulation rather than fully developed bursting.")
        end
        if modulated4d_diag.period_lag1 < -0.9 && modulated4d_diag.width_lag1 < -0.9
            println(io, "- Event-level diagnostics show an almost perfectly alternating long/short sequence in both inter-burst periods and burst widths, which is more consistent with a period-2-like relaxation pattern than with stochastic or weakly modulated variability.")
        end
        println(io, "- Floquet-based stability labels on the outer branch should still be treated cautiously because the eigensolver emits missing-conjugate-pair warnings during continuation.")
    end
end

function main()
    params = base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0)
    @info "Julia thread count" threads = nthreads()
    prob = make_2d_problem(params)
    prob3d = make_3d_problem(params)
    @info "Starting 2D equilibrium continuation" params
    br = continue_equilibria(prob)
    @info "Starting 3D equilibrium continuation" params
    br3d = continue_equilibria(prob3d)

    fold_idx = find_special_point_index(br, :fold)
    if isnothing(fold_idx)
        fold_idx = find_special_point_index(br, :bp)
    end
    hopf_idx = find_special_point_index(br, :hopf)
    fold3d_idx = find_special_point_index(br3d, :fold)
    if isnothing(fold3d_idx)
        fold3d_idx = find_special_point_index(br3d, :bp)
    end
    hopf3d_idx = find_special_point_index(br3d, :hopf)
    po_data = nothing
    hopf_data = nothing
    po_result = nothing

    if !isnothing(hopf_idx)
        hopf = br.specialpoint[hopf_idx]
        @info "Hopf point detected on reduced model branch" hopf

        try
            po_result = continue_periodic_from_hopf(prob, br, hopf_idx, params)
            hopf_data = po_result.refined
            if !isnothing(po_result.best)
                po_data = po_result.best.data
                @info "Periodic-orbit continuation succeeded from refined Hopf point" strategy = po_result.best.name points = po_result.best.points max_damp = po_result.best.max_damp
            else
                @warn "All Hopf-neighborhood periodic-orbit continuation attempts failed" failures = po_result.failures
            end
        catch err
            @warn "Refined Hopf / periodic-orbit workflow failed; keeping equilibrium-only outputs" err
        end
    else
        @info "No Hopf point detected on the current reduced-model branch"
    end

    E_values = range(0.0, 2.0, length = 31)
    sweep = sweep_multistability(E_values, params)
    E_values_3d = range(0.0, 2.0, length = 25)
    sweep3d = sweep_multistability_3d(E_values_3d, params)
    hopf3d_scan = isnothing(hopf3d_idx) ? nothing : scan_3d_hopf_neighborhood(range(max(0.0, br3d.specialpoint[hopf3d_idx].param - 0.06), min(2.0, br3d.specialpoint[hopf3d_idx].param + 0.03), length = 41), params)
    hopf3d_transition = isnothing(hopf3d_scan) ? (last_osc_idx = nothing, quiet_idx = nothing, last_osc_E = nothing, first_quiet_E = nothing) : summarize_3d_hopf_scan(hopf3d_scan)
    regime3d_E = collect(range(1.15, min(1.65, isnothing(hopf3d_idx) ? 1.65 : br3d.specialpoint[hopf3d_idx].param - 0.02), length = 8))
    regime3d_diag = diagnose_3d_relaxation_regime(regime3d_E, params)
    regime3d_summary = summarize_3d_relaxation_regime(regime3d_diag)
    targeted3d_E = collect(range(max(1.25, isnothing(hopf3d_idx) ? 1.45 : br3d.specialpoint[hopf3d_idx].param - 0.22), min(1.95, isnothing(hopf3d_idx) ? 1.68 : br3d.specialpoint[hopf3d_idx].param - 0.01), length = 8))
    targeted3d_search = targeted_3d_mmo_search(targeted3d_E, (0.03, 0.02, 0.01, 0.005), (0.02, 0.01, 0.005, 0.002, 0.001), params)

    rep_E = !isnothing(hopf_idx) ? min(br.specialpoint[hopf_idx].param + 0.08, 2.0) : 1.1
    rep_params = merge(params, (E = rep_E,))
    rep_sim = simulate_2d(rep_params; u0 = [0.08, 0.0], tmax = 4000.0, transient = 2200.0, saveat = 0.2)

    signal3d = max.(sweep3d.amp_low, sweep3d.amp_high) .+ sweep3d.gap_d .+ 0.25 .* sweep3d.gap_h
    rep3d_idx = argmax(signal3d)
    rep3d_E = sweep3d.E[rep3d_idx]
    rep3d_u0 = sweep3d.amp_high[rep3d_idx] >= sweep3d.amp_low[rep3d_idx] ? [0.95, 0.5, 0.9] : [0.05, 0.0, 0.2]
    rep3d_params = merge(params, (E = rep3d_E,))
    sim3d = simulate_3d(rep3d_params; u0 = rep3d_u0, tmax = 5000.0, transient = 2400.0, saveat = 0.2)

    sweep4d = sweep_4d_grid((0.9, 1.0, 1.05, 1.1, 1.2), (0.9, 1.0, 1.1, 1.2), merge(params, (eta_k = 0.04, eta_h = 0.015)))
    targeted4d_search = targeted_4d_burst_search((0.85, 0.9, 0.95, 1.0, 1.05), (0.9, 1.0, 1.1), (0.015, 0.01, 0.005, 0.002), (0.02, 0.01, 0.005, 0.003), params)
    geometry4d_params = merge(params, (G = 1.0, epsilon = 0.01, eta_h = 0.002, eta_k = 0.003))
    geometry4d_search = targeted_4d_geometry_search((-5.1, -5.0, -4.9), (1.4, 1.5, 1.6), (1.05, 1.07, 1.09), geometry4d_params)
    best4d_idx = argmax(sweep4d.damp)
    best4d_G_idx, best4d_E_idx = Tuple(CartesianIndices(sweep4d.damp)[best4d_idx])
    full_params = merge(params, (E = sweep4d.E[best4d_E_idx], G = sweep4d.G[best4d_G_idx], eta_k = 0.04, eta_h = 0.015))
    full_diag = diagnose_4d_candidate(full_params)
    modulated4d_params = merge(geometry4d_params, (theta0 = geometry4d_search.best.theta0, beta_n = geometry4d_search.best.beta_n, E = geometry4d_search.best.E))
    modulated4d_diag = diagnose_4d_candidate(modulated4d_params)

    save_fold_screen_plot()
    save_branch_plot(br, po_data)
    save_hopf_neighborhood_plot(br, hopf_data, po_data)
    save_hysteresis_plot(sweep)
    save_phase_plot(rep_sim, rep_params)
    save_3d_branch_plot(br3d)
    save_3d_hysteresis_plot(sweep3d)
    if !isnothing(hopf3d_scan)
        save_3d_hopf_zoom_plot(hopf3d_scan, br3d.specialpoint[hopf3d_idx].param)
    end
    save_3d_regime_plot(regime3d_diag)
    save_3d_targeted_search_plot(targeted3d_search)
    save_3d_timeseries_plot(sim3d)
    save_4d_regime_plot(sweep4d)
    save_4d_targeted_search_plot(targeted4d_search)
    save_4d_geometry_search_plot(geometry4d_search)
    save_4d_modulated_candidate_plot(modulated4d_diag.sim)
    save_4d_alternation_plot(modulated4d_diag)
    save_full_plot(full_diag.sim)
    save_summary(br, po_data, sweep, rep_sim, br3d, sweep3d, hopf3d_scan, hopf3d_transition, regime3d_diag, regime3d_summary, targeted3d_search, sim3d, sweep4d, targeted4d_search, geometry4d_search, full_diag, modulated4d_diag, params, hopf_idx, fold_idx, hopf_data, po_result, hopf3d_idx, fold3d_idx)

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
    if !isnothing(hopf_data)
        println(@sprintf("Refined Hopf: E ≈ %.10f, type = %s, period ≈ %.10f", hopf_data.point.p, String(hopf_data.point.type), 2pi / abs(hopf_data.point.ω)))
    end
    if !isnothing(po_data)
        println(@sprintf("Extended periodic branch: %d points, E ∈ [%.10f, %.10f], max d-amplitude = %.6e",
            length(po_data.param), minimum(po_data.param), maximum(po_data.param), maximum(po_data.damp)))
        pd_points = dedupe_specialpoints(filter(sp -> sp.type == :pd, po_data.specialpoints))
        if !isempty(pd_points)
            sp = pd_points[end]
            println(@sprintf("Representative period-doubling point: E ≈ %.10f [idx = %d]", sp.param, sp.idx))
        end
    end
    println(@sprintf("Max low/high branch gap from direct sweeps = %.6f", maximum(sweep.gap)))
    if !isnothing(fold3d_idx)
        println(@sprintf("3D fold/branch-point candidate (%s): E ≈ %.10f", String(br3d.specialpoint[fold3d_idx].type), br3d.specialpoint[fold3d_idx].param))
    else
        println("3D fold point: none detected on the current continuation branch")
    end
    if !isnothing(hopf3d_idx)
        println(@sprintf("3D Hopf point: E ≈ %.10f", br3d.specialpoint[hopf3d_idx].param))
    else
        println("3D Hopf point: none detected on the current continuation branch")
    end
    println(@sprintf("Max 3D low/high mean-d gap from direct sweeps = %.6f", maximum(sweep3d.gap_d)))
    if !isnothing(hopf3d_transition.last_osc_E) && !isnothing(hopf3d_transition.first_quiet_E)
        println(@sprintf("3D Hopf-neighborhood oscillation-collapse window: E ≈ %.10f -> %.10f", hopf3d_transition.last_osc_E, hopf3d_transition.first_quiet_E))
    end
    println(@sprintf("3D oscillatory-regime max period CV = %.6f, max mid-band occupancy = %.6f", regime3d_summary.max_period_cv, regime3d_summary.max_mid_frac))
    println(@sprintf("Best targeted 3D MMO-search candidate: epsilon = %.4f, eta_h = %.4f, E = %.6f, score = %.6f, period CV = %.6f, subcycle CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
        targeted3d_search.best.epsilon,
        targeted3d_search.best.eta_h,
        targeted3d_search.best.E,
        targeted3d_search.best.score,
        targeted3d_search.best.period_cv,
        targeted3d_search.best.subcycle_cv,
        targeted3d_search.best.mid_frac,
        targeted3d_search.best.crossing_ratio,
    ))
    println(@sprintf("Representative 3D amplitudes: d = %.6f, q = %.6f, h = %.6f", sim3d.damp, sim3d.qamp, sim3d.hamp))
    println(@sprintf("Best sustained 4D candidate: E = %.3f, G = %.3f, d = %.6f, q = %.6f, k = %.6f, h = %.6f",
        full_params.E, full_params.G, full_diag.sim.damp, full_diag.sim.qamp, full_diag.kamp, full_diag.hamp))
    println(@sprintf("4D candidate period CV = %.6f, mid-band occupancy = %.6f", full_diag.period_cv, full_diag.mid_frac))
    println(@sprintf("Best targeted 4D burst-search candidate: eta_h = %.4f, eta_k = %.4f, E = %.3f, G = %.3f, score = %.6f, inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
        targeted4d_search.best.eta_h,
        targeted4d_search.best.eta_k,
        targeted4d_search.best.E,
        targeted4d_search.best.G,
        targeted4d_search.best.score,
        targeted4d_search.best.period_cv,
        targeted4d_search.best.burst_width_cv,
        targeted4d_search.best.mid_frac,
        targeted4d_search.best.crossing_ratio,
    ))
    println(@sprintf("Best 4D geometry-search candidate: theta0 = %.3f, beta_n = %.3f, E = %.3f, score = %.6f, inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
        geometry4d_search.best.theta0,
        geometry4d_search.best.beta_n,
        geometry4d_search.best.E,
        geometry4d_search.best.score,
        geometry4d_search.best.period_cv,
        geometry4d_search.best.burst_width_cv,
        geometry4d_search.best.mid_frac,
        geometry4d_search.best.crossing_ratio,
    ))
    println(@sprintf("Long-window geometry candidate diagnostic: inter-burst CV = %.6f, burst-width CV = %.6f, mid-band occupancy = %.6f, crossing ratio = %.6f",
        modulated4d_diag.period_cv,
        modulated4d_diag.burst_width_cv,
        modulated4d_diag.mid_frac,
        modulated4d_diag.crossing_ratio,
    ))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
