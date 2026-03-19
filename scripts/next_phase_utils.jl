using Printf
using Statistics
using LinearAlgebra
using Base.Threads: @threads, nthreads

using Plots

include(joinpath(@__DIR__, "double_frontier_rollout.jl"))

const NEXT_PHASE_RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results", "next_phase"))
mkpath(NEXT_PHASE_RESULTS_DIR)

const DEFAULT_FULL_U0 = [0.05, 0.0, 0.1, 0.6]
const DEFAULT_3D_LOW_U0 = [0.05, 0.0, 0.2]
const DEFAULT_3D_HIGH_U0 = [0.95, 0.5, 0.9]
const DEFAULT_2D_LOW_U0 = [0.05, 0.0]
const DEFAULT_2D_CANARD_U0 = [0.08, 0.0]

function parse_mode(args = ARGS)
    if any(arg -> arg == "--full", args)
        return :full
    elseif any(arg -> arg == "--quick", args)
        return :quick
    else
        return :quick
    end
end

function mode_label(mode::Symbol)
    return mode === :full ? "full" : "quick"
end

function to_float_array(x)
    return collect(Float64, x)
end

function ensure_parent_dir(path::AbstractString)
    mkpath(dirname(path))
    return path
end

function write_csv(path::AbstractString, rows::Vector{<:NamedTuple}, columns::Vector{Symbol})
    ensure_parent_dir(path)
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            vals = String[]
            for col in columns
                val = getproperty(row, col)
                if val isa Bool
                    push!(vals, val ? "1" : "0")
                elseif val isa Integer
                    push!(vals, string(val))
                elseif val isa AbstractFloat
                    if isnan(val)
                        push!(vals, "NaN")
                    else
                        push!(vals, @sprintf("%.10f", val))
                    end
                else
                    sval = replace(string(val), '"' => '\'')
                    if occursin(",", sval)
                        push!(vals, "\"$sval\"")
                    else
                        push!(vals, sval)
                    end
                end
            end
            println(io, join(vals, ","))
        end
    end
    return path
end

function default_param_order()
    return [
        :E, :G, :z, :epsilon, :eta_h, :eta_k, :beta_n, :lambda_D, :theta0,
        :lambda_A, :c_A, :alpha0, :alpha_z, :alpha_h, :alpha_G,
        :c_D, :theta_z, :theta_h, :theta_G, :theta_E, :theta_k, :theta_phi,
        :gamma_q, :omega0, :omega_G, :omega_h, :omega_q,
        :nu0, :nu_G, :nu_h, :chi_I, :chi_R, :rho,
        :delta_k, :psi_M, :psi_D, :vartheta_M, :vartheta_D, :vartheta_H,
        :mu0, :mu1, :sigma, :hbar, :kbar,
    ]
end

function parameter_lines(p)
    lines = String[]
    for key in default_param_order()
        if hasproperty(p, key)
            val = getproperty(p, key)
            if val isa AbstractFloat
                push!(lines, @sprintf("- %s = %.10f", String(key), val))
            else
                push!(lines, "- $(String(key)) = $(val)")
            end
        end
    end
    return lines
end

function parameter_string(p)
    return join(parameter_lines(p), "\n")
end

function ordered_param_row(p)
    pairs = Pair{Symbol, Any}[]
    for key in default_param_order()
        if hasproperty(p, key)
            push!(pairs, key => getproperty(p, key))
        end
    end
    return (; pairs...)
end

function raw_doer_target(d, q, k, h, p)
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
    return logistic(p.lambda_D * (score - p.c_D))
end

function state_series(sim, p)
    d = to_float_array(sim.d)
    q = to_float_array(sim.q)
    k = hasproperty(sim, :k) ? to_float_array(sim.k) : fill(Float64(p.kbar), length(d))
    h = hasproperty(sim, :h) ? to_float_array(sim.h) : fill(Float64(p.hbar), length(d))
    return (d = d, q = q, k = k, h = h)
end

function advisor_frontier_series(sim, p)
    s = state_series(sim, p)
    return [advisor_frontier(h, p) for h in s.h]
end

function raw_doer_target_series(sim, p)
    s = state_series(sim, p)
    return [raw_doer_target(s.d[i], s.q[i], s.k[i], s.h[i], p) for i in eachindex(s.d)]
end

function frontier_metric_series(sim, p)
    s = state_series(sim, p)
    a_t = advisor_frontier_series(sim, p)
    r_t = raw_doer_target_series(sim, p)
    phi_t = [supervision(s.h[i], s.q[i], p) for i in eachindex(s.q)]
    m_t = max.(a_t .- s.d, 0.0)
    capgap_t = max.(r_t .- a_t, 0.0)
    capbind_t = Float64.(r_t .> a_t)
    return (
        a_t = a_t,
        r_t = r_t,
        m_t = m_t,
        capbind_t = capbind_t,
        capgap_t = capgap_t,
        phi_t = phi_t,
    )
end

function frontier_metric_summary(sim, p)
    s = state_series(sim, p)
    fm = frontier_metric_series(sim, p)
    return (
        m_bar = mean(fm.m_t),
        m_max = maximum(fm.m_t),
        m_q05 = quantile(fm.m_t, 0.05),
        m_q50 = quantile(fm.m_t, 0.50),
        m_q95 = quantile(fm.m_t, 0.95),
        capbind_frac = mean(fm.capbind_t),
        capgap_mean = mean(fm.capgap_t),
        phi_bar = mean(fm.phi_t),
        d_mean = mean(s.d),
        q_mean = mean(s.q),
        k_mean = mean(s.k),
        h_mean = mean(s.h),
    )
end

function representative_parameter_sets()
    rollout_base = base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0)
    canard_base = merge(
        base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0, epsilon = 0.02),
        (chi_R = 0.8, gamma_q = 2.0, nu0 = 0.1),
    )
    ridge_base = merge(
        base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0),
        (G = 1.0, epsilon = 0.01, eta_h = 0.002, eta_k = 0.003),
    )

    return [
        (
            name = "base_default",
            label = "default base_params()",
            dim = :full,
            source = "base_params() default",
            p = base_params(),
            u0 = copy(DEFAULT_FULL_U0),
        ),
        (
            name = "twoD_refined_hopf",
            label = "2D refined Hopf",
            dim = :twoD,
            source = "scripts/double_frontier_rollout.jl: refined Hopf from base_params(beta_n=2, lambda_D=20, theta0=-5)",
            p = merge(rollout_base, (E = 0.974752539138242,)),
            u0 = copy(DEFAULT_2D_LOW_U0),
        ),
        (
            name = "twoD_large_cycle",
            label = "2D representative large-cycle point",
            dim = :twoD,
            source = "scripts/double_frontier_canard_search.jl: representative above-jump simulation at epsilon=0.02",
            p = merge(canard_base, (E = 1.216187,)),
            u0 = copy(DEFAULT_2D_CANARD_U0),
        ),
        (
            name = "threeD_oscillatory",
            label = "3D representative oscillatory point",
            dim = :threeD,
            source = "scripts/double_frontier_rollout.jl: rep3d_idx = argmax(max(amp_low, amp_high) + gap_d + 0.25*gap_h)",
            p = merge(rollout_base, (E = 1.1666666666666667,)),
            u0 = copy(DEFAULT_3D_LOW_U0),
        ),
        (
            name = "ridge_peak",
            label = "4D ridge_peak",
            dim = :full,
            source = "results/double_frontier_4d_period2_check.txt representative point",
            p = merge(ridge_base, (theta0 = -4.90, E = 1.08, beta_n = 1.5)),
            u0 = copy(DEFAULT_FULL_U0),
        ),
        (
            name = "ridge_inner",
            label = "4D ridge_inner",
            dim = :full,
            source = "results/double_frontier_4d_period2_check.txt representative point",
            p = merge(ridge_base, (theta0 = -4.90, E = 1.07, beta_n = 1.5)),
            u0 = copy(DEFAULT_FULL_U0),
        ),
        (
            name = "ridge_end",
            label = "4D ridge_end",
            dim = :full,
            source = "results/double_frontier_4d_period2_check.txt representative point",
            p = merge(ridge_base, (theta0 = -4.86, E = 1.05, beta_n = 1.5)),
            u0 = copy(DEFAULT_FULL_U0),
        ),
        (
            name = "boundary",
            label = "4D boundary",
            dim = :full,
            source = "results/double_frontier_4d_period2_check.txt representative point",
            p = merge(ridge_base, (theta0 = -4.84, E = 1.04, beta_n = 1.5)),
            u0 = copy(DEFAULT_FULL_U0),
        ),
    ]
end

function parameter_point(name::AbstractString)
    for point in representative_parameter_sets()
        if point.name == name
            return point
        end
    end
    error("Unknown parameter point: $name")
end

function simulate_point(point; tmax, transient, saveat)
    if point.dim === :twoD
        return simulate_2d(point.p; u0 = point.u0, tmax = tmax, transient = transient, saveat = saveat)
    elseif point.dim === :threeD
        return simulate_3d(point.p; u0 = point.u0, tmax = tmax, transient = transient, saveat = saveat)
    else
        return simulate_full(point.p; u0 = point.u0, tmax = tmax, transient = transient, saveat = saveat)
    end
end

function threshold_crossing_data(sim, p; thresh = 0.9)
    s = state_series(sim, p)
    t = to_float_array(sim.t)
    times = Float64[]
    qvals = Float64[]
    kvals = Float64[]
    hvals = Float64[]
    dvals = Float64[]
    d_dots = Float64[]

    for i in 2:length(t)
        d0 = s.d[i - 1]
        d1 = s.d[i]
        if d0 < thresh <= d1
            α = (thresh - d0) / (d1 - d0)
            tc = t[i - 1] + α * (t[i] - t[i - 1])
            qc = s.q[i - 1] + α * (s.q[i] - s.q[i - 1])
            kc = s.k[i - 1] + α * (s.k[i] - s.k[i - 1])
            hc = s.h[i - 1] + α * (s.h[i] - s.h[i - 1])
            dc = thresh
            du = zeros(4)
            double_frontier_full!(du, [dc, qc, kc, hc], p, tc)
            if du[1] > 0
                push!(times, tc)
                push!(qvals, qc)
                push!(kvals, kc)
                push!(hvals, hc)
                push!(dvals, dc)
                push!(d_dots, du[1])
            end
        end
    end

    return (
        t = times,
        q = qvals,
        k = kvals,
        h = hvals,
        d = dvals,
        d_dot = d_dots,
    )
end

function width_sequence(sim; thresh = 0.9)
    return threshold_episodes(to_float_array(sim.t), to_float_array(sim.d); thresh = thresh).widths
end

function lag2corr(values)
    if length(values) <= 3
        return 0.0
    end
    x1 = values[1:end-2]
    x2 = values[3:end]
    sx = std(x1)
    sy = std(x2)
    if iszero(sx) || iszero(sy)
        return 0.0
    end
    return cor(x1, x2)
end

function odd_even_stats(values)
    odd = values[1:2:end]
    even = values[2:2:end]
    return (
        odd_mean = isempty(odd) ? NaN : mean(odd),
        even_mean = isempty(even) ? NaN : mean(even),
        odd_cv = isempty(odd) ? NaN : safe_cv(odd),
        even_cv = isempty(even) ? NaN : safe_cv(even),
    )
end

function mean_pair_distance(x, lag)
    n = length(x)
    if n <= lag
        return NaN
    end
    vals = Float64[]
    for i in 1:(n - lag)
        push!(vals, norm(x[i + lag] - x[i]))
    end
    return mean(vals)
end

function poincare_metrics(crossings, widths = Float64[])
    times = crossings.t
    q = crossings.q
    k = crossings.k
    h = crossings.h
    T = length(times) > 1 ? diff(times) : Float64[]
    n = min(length(times), length(widths))
    W = n > 0 ? widths[1:n] : Float64[]
    X = [[q[i], k[i], h[i]] for i in eachindex(q)]
    Tstats = odd_even_stats(T)
    Wstats = odd_even_stats(W)
    qstats = odd_even_stats(q)
    hstats = odd_even_stats(h)
    return (
        T = T,
        W = W,
        lag1 = lag1corr(T),
        lag2 = lag2corr(T),
        width_lag1 = lag1corr(W),
        width_lag2 = lag2corr(W),
        T_odd_mean = Tstats.odd_mean,
        T_even_mean = Tstats.even_mean,
        T_odd_cv = Tstats.odd_cv,
        T_even_cv = Tstats.even_cv,
        W_odd_mean = Wstats.odd_mean,
        W_even_mean = Wstats.even_mean,
        W_odd_cv = Wstats.odd_cv,
        W_even_cv = Wstats.even_cv,
        q_odd_mean = qstats.odd_mean,
        q_even_mean = qstats.even_mean,
        h_odd_mean = hstats.odd_mean,
        h_even_mean = hstats.even_mean,
        R1 = mean_pair_distance(X, 1),
        R2 = mean_pair_distance(X, 2),
        R2_over_R1 = isnan(mean_pair_distance(X, 1)) || iszero(mean_pair_distance(X, 1)) ? NaN : mean_pair_distance(X, 2) / mean_pair_distance(X, 1),
    )
end

function period2_strength_from_crossings(crossings, widths = Float64[])
    pm = poincare_metrics(crossings, widths)
    if isempty(pm.T)
        return NaN
    end
    Tmean = mean(pm.T)
    Wmean = isempty(pm.W) ? NaN : mean(pm.W)
    Tgap = (isnan(pm.T_odd_mean) || isnan(pm.T_even_mean) || iszero(Tmean)) ? 0.0 : abs(pm.T_even_mean - pm.T_odd_mean) / Tmean
    Wgap = (isempty(pm.W) || isnan(pm.W_odd_mean) || isnan(pm.W_even_mean) || iszero(Wmean)) ? 0.0 : abs(pm.W_even_mean - pm.W_odd_mean) / Wmean
    return 0.5 * (Tgap + Wgap) - 0.5 * (pm.lag1 + pm.width_lag1) + 0.5 * (pm.lag2 + pm.width_lag2)
end

function cheap_2d_bistability_proxy(p; tmax = 2500.0, transient = 1200.0, saveat = 0.5, gap_threshold = 0.05)
    a_cap = advisor_frontier(p.hbar, p)
    high_init = [min(0.95, max(0.15, a_cap - 0.02)), 0.5]
    sim_low = simulate_2d(p; u0 = copy(DEFAULT_2D_LOW_U0), tmax = tmax, transient = transient, saveat = saveat)
    sim_high = simulate_2d(p; u0 = high_init, tmax = tmax, transient = transient, saveat = saveat)
    gap = abs(sim_high.d[end] - sim_low.d[end])
    return (
        gap = gap,
        has_bistability = gap > gap_threshold,
        low_d_end = sim_low.d[end],
        high_d_end = sim_high.d[end],
        low_amp = sim_low.damp,
        high_amp = sim_high.damp,
    )
end

function sustained_oscillation_proxy(sim; amp_threshold = 0.1)
    return (maximum(sim.d) - minimum(sim.d)) > amp_threshold
end

function cheap_4d_summary(p; tmax = 8000.0, transient = 4000.0, saveat = 0.5, thresh = 0.9)
    sim = simulate_full(p; u0 = copy(DEFAULT_FULL_U0), tmax = tmax, transient = transient, saveat = saveat)
    fm = frontier_metric_summary(sim, p)
    crossings = threshold_crossing_data(sim, p; thresh = thresh)
    widths = width_sequence(sim; thresh = thresh)
    pm = poincare_metrics(crossings, widths)
    return (
        sim = sim,
        frontier = fm,
        crossings = crossings,
        widths = widths,
        poincare = pm,
        sustained_oscillation = sustained_oscillation_proxy(sim),
        d_amp = sim.damp,
        q_amp = sim.qamp,
        k_amp = maximum(sim.k) - minimum(sim.k),
        h_amp = maximum(sim.h) - minimum(sim.h),
        period2_strength = period2_strength_from_crossings(crossings, widths),
    )
end

function heatmap_matrix(xs, ys, rows, field)
    mat = Matrix{Float64}(undef, length(ys), length(xs))
    lookup = Dict((row.eta_h, row.E) => row for row in rows)
    for (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        mat[iy, ix] = getproperty(lookup[(y, x)], field)
    end
    return mat
end

function print_generated_files(paths::Vector{String})
    println("Generated files:")
    for path in paths
        println(path)
    end
end
