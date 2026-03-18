using Printf

using Plots

include(joinpath(@__DIR__, "double_frontier_rollout.jl"))

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results"))
mkpath(RESULTS_DIR)

function scan_relaxation_candidates()
    base = base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0, epsilon = 0.01)
    candidates = NamedTuple[]

    for chi_R in (0.8, 1.2)
        for gamma_q in (2.0, 2.8)
            for nu0 in (0.1, 0.2, 0.3)
                best_amp = -Inf
                best_E = NaN
                best_sim = nothing

                for E in range(1.10, 1.35, length = 11)
                    p = merge(base, (chi_R = chi_R, gamma_q = gamma_q, nu0 = nu0, E = E))
                    sim = simulate_2d(p; u0 = [0.08, 0.0], tmax = 2600.0, transient = 1400.0, saveat = 0.2)
                    if sim.damp > best_amp
                        best_amp = sim.damp
                        best_E = E
                        best_sim = sim
                    end
                end

                push!(candidates, (
                    chi_R = chi_R,
                    gamma_q = gamma_q,
                    nu0 = nu0,
                    best_E = best_E,
                    best_amp = best_amp,
                    dmin = isnothing(best_sim) ? NaN : best_sim.dmin,
                    dmax = isnothing(best_sim) ? NaN : best_sim.dmax,
                ))
            end
        end
    end

    sort!(candidates; by = x -> x.best_amp, rev = true)
    return candidates
end

function refine_jump_scan(base, epsilons)
    scans = Dict{Float64, NamedTuple}()

    for epsilon in epsilons
        E_vals = collect(range(1.205, 1.223, length = 10))
        amps = Float64[]
        dmins = Float64[]
        dmaxs = Float64[]

        for E in E_vals
            p = merge(base, (epsilon = epsilon, E = E))
            sim = simulate_2d(p; u0 = [0.08, 0.0], tmax = 4200.0, transient = 2200.0, saveat = 0.2)
            push!(amps, sim.damp)
            push!(dmins, sim.dmin)
            push!(dmaxs, sim.dmax)
        end

        diffs = diff(amps)
        jump_idx = argmax(diffs)
        scans[epsilon] = (
            E = E_vals,
            amp = amps,
            dmin = dmins,
            dmax = dmaxs,
            jump_idx = jump_idx,
            jump_left = E_vals[jump_idx],
            jump_right = E_vals[jump_idx + 1],
            jump_size = diffs[jump_idx],
        )
    end

    return scans
end

function representative_simulations(base, epsilon)
    below = simulate_2d(merge(base, (epsilon = epsilon, E = 1.215)); u0 = [0.08, 0.0], tmax = 6000.0, transient = 2500.0, saveat = 0.1)
    above = simulate_2d(merge(base, (epsilon = epsilon, E = 1.217)); u0 = [0.08, 0.0], tmax = 6000.0, transient = 2500.0, saveat = 0.1)
    return below, above
end

function save_candidate_plot(candidates)
    top = candidates[1:min(end, 10)]
    labels = [@sprintf("chi_R=%.1f, gamma=%.1f, nu0=%.1f", c.chi_R, c.gamma_q, c.nu0) for c in top]
    amps = [c.best_amp for c in top]

    p = bar(
        labels,
        amps;
        xlabel = "candidate",
        ylabel = "max 2D d-amplitude",
        title = "Top 2D relaxation candidates from coarse screening",
        legend = false,
        xrotation = 35,
        size = (1200, 500),
    )
    savefig(p, joinpath(RESULTS_DIR, "double_frontier_2d_candidate_scan.png"))
end

function save_jump_plot(scans)
    p1 = plot(size = (1200, 800), layout = (2, 1))
    p2 = plot(size = (1200, 800), layout = (2, 1))

    for epsilon in sort(collect(keys(scans)); rev = true)
        scan = scans[epsilon]
        plot!(p1[1], scan.E, scan.amp; lw = 2, marker = :circle, label = @sprintf("epsilon = %.3f", epsilon))
        plot!(p1[2], scan.E, scan.dmin; lw = 2, marker = :circle, label = @sprintf("dmin, epsilon = %.3f", epsilon))
        plot!(p2[1], scan.E, scan.dmax; lw = 2, marker = :circle, label = @sprintf("dmax, epsilon = %.3f", epsilon))
    end

    xlabel!(p1[1], "E")
    ylabel!(p1[1], "d amplitude")
    title!(p1[1], "2D direct-simulation jump scan near the oscillation onset")
    xlabel!(p1[2], "E")
    ylabel!(p1[2], "dmin")
    title!(p1[2], "Lower envelope collapses toward 0 on the large cycle")
    xlabel!(p2[1], "E")
    ylabel!(p2[1], "dmax")
    title!(p2[1], "Upper envelope stays near the doer cap")

    combo = plot(p1[1], p1[2], p2[1]; layout = (3, 1), size = (1100, 1100))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_2d_jump_scan.png"))
end

function save_representative_plot(below, above)
    p1 = plot(
        below.t,
        below.d;
        lw = 2,
        xlabel = "t",
        ylabel = "d(t)",
        title = "Below the jump window: converges to equilibrium",
        label = @sprintf("E = %.3f", below.E),
    )
    plot!(p1, below.t, below.q; lw = 2, label = "q")

    p2 = plot(
        above.t,
        above.d;
        lw = 2,
        xlabel = "t",
        ylabel = "d(t)",
        title = "Above the jump window: large relaxation cycle",
        label = @sprintf("E = %.3f", above.E),
    )
    plot!(p2, above.t, above.q; lw = 2, label = "q")

    p3 = plot(
        below.d,
        below.q;
        lw = 2,
        xlabel = "d",
        ylabel = "q",
        title = "Phase portrait below the jump",
        label = "trajectory",
    )

    p4 = plot(
        above.d,
        above.q;
        lw = 2,
        xlabel = "d",
        ylabel = "q",
        title = "Phase portrait above the jump",
        label = "trajectory",
    )

    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1200, 800))
    savefig(combo, joinpath(RESULTS_DIR, "double_frontier_2d_relaxation_examples.png"))
end

function save_summary(candidates, scans, below, above, base)
    path = joinpath(RESULTS_DIR, "double_frontier_2d_canard_search.txt")
    open(path, "w") do io
        println(io, "2D canard / relaxation search for the double-frontier model")
        println(io, @sprintf("Base parameters: beta_n = %.2f, lambda_D = %.2f, theta0 = %.2f", base.beta_n, base.lambda_D, base.theta0))
        println(io, @sprintf("Selected feedback candidate: chi_R = %.2f, gamma_q = %.2f, nu0 = %.2f", base.chi_R, base.gamma_q, base.nu0))
        println(io)
        println(io, "Top coarse-scan candidates by maximum post-transient d-amplitude:")
        for c in candidates[1:min(end, 5)]
            println(io, @sprintf("- chi_R = %.2f, gamma_q = %.2f, nu0 = %.2f, best E = %.3f, max amp = %.6f, d-range = [%.6f, %.6f]",
                c.chi_R, c.gamma_q, c.nu0, c.best_E, c.best_amp, c.dmin, c.dmax))
        end
        println(io)
        println(io, "Refined jump windows:")
        for epsilon in sort(collect(keys(scans)); rev = true)
            scan = scans[epsilon]
            println(io, @sprintf("- epsilon = %.3f: jump between E = %.6f and E = %.6f with amplitude increase %.6f",
                epsilon, scan.jump_left, scan.jump_right, scan.jump_size))
        end
        println(io)
        println(io, @sprintf("Representative below-jump simulation: E = %.3f, d_amp = %.6f", below.E, below.damp))
        println(io, @sprintf("Representative above-jump simulation: E = %.3f, d_amp = %.6f, d-range = [%.6f, %.6f]",
            above.E, above.damp, above.dmin, above.dmax))
        println(io)
        println(io, "Interpretation:")
        println(io, "- The coarse screen finds genuine stable large-amplitude 2D oscillations once q-feedback is strengthened.")
        println(io, "- In the selected parameter set the post-transient amplitude jumps from numerical zero to O(1) over a very narrow E interval.")
        println(io, "- This is strong evidence for a relaxation-oscillation regime and a canard-like onset, but not yet a full proof of a classical canard explosion.")
        println(io, "- To upgrade the claim from canard-like to canard, the next step is higher-accuracy periodic-orbit continuation anchored to the Hopf neighborhood.")
    end
end

function main()
    candidates = scan_relaxation_candidates()
    chosen = candidates[1]
    chosen_base = merge(
        base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0),
        (chi_R = chosen.chi_R, gamma_q = chosen.gamma_q, nu0 = chosen.nu0),
    )

    scans = refine_jump_scan(chosen_base, (0.05, 0.02, 0.01))
    below, above = representative_simulations(chosen_base, 0.02)

    save_candidate_plot(candidates)
    save_jump_plot(scans)
    save_representative_plot(below, above)
    save_summary(candidates, scans, below, above, chosen_base)

    println()
    println("Saved canard/relaxation search results to: ", RESULTS_DIR)
    println(@sprintf("Selected candidate: chi_R = %.2f, gamma_q = %.2f, nu0 = %.2f", chosen.chi_R, chosen.gamma_q, chosen.nu0))
    for epsilon in sort(collect(keys(scans)); rev = true)
        scan = scans[epsilon]
        println(@sprintf("epsilon = %.3f jump window: E = %.6f -> %.6f, amplitude increase = %.6f",
            epsilon, scan.jump_left, scan.jump_right, scan.jump_size))
    end
    println(@sprintf("Below jump: E = %.3f, d_amp = %.6f", below.E, below.damp))
    println(@sprintf("Above jump: E = %.3f, d_amp = %.6f", above.E, above.damp))
end

main()
