using Printf
using Statistics

include(joinpath(@__DIR__, "next_phase_utils.jl"))

function persistence_grid(mode::Symbol)
    if mode === :full
        return (
            E = collect(range(1.0, 1.8, length = 31)),
            eta_h = collect(range(0.001, 0.03, length = 21)),
            sim = (tmax = 9000.0, transient = 4500.0, saveat = 0.2),
        )
    else
        return (
            E = collect(range(1.0, 1.8, length = 16)),
            eta_h = collect(range(0.001, 0.03, length = 11)),
            sim = (tmax = 6000.0, transient = 3000.0, saveat = 0.2),
        )
    end
end

function persistence_row(E, eta_h, sim, p)
    fm = frontier_metric_summary(sim, p)
    diag = diagnose_3d_pattern(sim)
    return (
        E = E,
        eta_h = eta_h,
        d_amp = sim.damp,
        q_amp = sim.qamp,
        h_amp = sim.hamp,
        return_time_mean = diag.period_mean,
        return_time_cv = diag.period_cv,
        crossing_ratio = diag.crossing_ratio,
        m_bar = fm.m_bar,
        capbind_frac = fm.capbind_frac,
        phi_bar = fm.phi_bar,
        mid_frac_proxy = diag.mid_frac,
    )
end

function run_persistence(mode)
    grid = persistence_grid(mode)
    base = base_params(beta_n = 2.0, lambda_D = 20.0, theta0 = -5.0)
    rows = Vector{NamedTuple}(undef, length(grid.E) * length(grid.eta_h))
    cells = CartesianIndices((length(grid.eta_h), length(grid.E)))

    @threads for linear_idx in eachindex(rows)
        iη, iE = Tuple(cells[linear_idx])
        p = merge(base, (E = grid.E[iE], eta_h = grid.eta_h[iη]))
        sim = simulate_3d(p; u0 = copy(DEFAULT_3D_HIGH_U0), grid.sim...)
        rows[linear_idx] = persistence_row(grid.E[iE], grid.eta_h[iη], sim, p)
    end

    return (rows = rows, E = grid.E, eta_h = grid.eta_h, k_fixed = base.kbar)
end

function save_persistence_plot(path, data)
    rows = data.rows
    d_amp = heatmap_matrix(data.E, data.eta_h, rows, :d_amp)
    rt_cv = heatmap_matrix(data.E, data.eta_h, rows, :return_time_cv)
    m_bar = heatmap_matrix(data.E, data.eta_h, rows, :m_bar)
    capbind = heatmap_matrix(data.E, data.eta_h, rows, :capbind_frac)

    p1 = heatmap(data.E, data.eta_h, d_amp; xlabel = "E", ylabel = "eta_h", title = "3D d_amp", colorbar_title = "d amp")
    p2 = heatmap(data.E, data.eta_h, rt_cv; xlabel = "E", ylabel = "eta_h", title = "3D return_time_cv", colorbar_title = "CV")
    p3 = heatmap(data.E, data.eta_h, m_bar; xlabel = "E", ylabel = "eta_h", title = "3D true m_bar", colorbar_title = "m_bar")
    p4 = heatmap(data.E, data.eta_h, capbind; xlabel = "E", ylabel = "eta_h", title = "3D capbind_frac", colorbar_title = "frac")
    combo = plot(p1, p2, p3, p4; layout = (2, 2), size = (1300, 900))
    savefig(combo, path)
    return path
end

function suspicious_regions(rows)
    suspects = filter(row -> row.return_time_cv > 0.03 || row.crossing_ratio > 1.05 || row.mid_frac_proxy > 0.01, rows)
    sort!(suspects; by = row -> (row.return_time_cv + abs(row.crossing_ratio - 1.0) + row.mid_frac_proxy), rev = true)
    return suspects[1:min(end, 5)]
end

function write_persistence_summary(path, mode, data)
    rows = data.rows
    max_d_amp_row = rows[argmax([row.d_amp for row in rows])]
    max_m_row = rows[argmax([row.m_bar for row in rows])]
    max_capbind_row = rows[argmax([row.capbind_frac for row in rows])]
    max_cv_row = rows[argmax([row.return_time_cv for row in rows])]
    suspects = suspicious_regions(rows)

    open(path, "w") do io
        println(io, "# 3D Persistence Scan")
        println(io)
        println(io, "- Mode: `$(mode_label(mode))`")
        println(io, @sprintf("- Fixed 3D k = %.6f", data.k_fixed))
        println(io, "- Base parameter family: `base_params(beta_n=2, lambda_D=20, theta0=-5)` with only `E` and `eta_h` varied")
        println(io, "- Legacy proxy retained as `mid_frac_proxy = share of time with 0.2 < d < 0.9`; it is not the true middle band.")
        println(io)
        println(io, "## Highlights")
        println(io)
        println(io, @sprintf("- Max d_amp occurs at E = %.6f, eta_h = %.6f with d_amp = %.6f.", max_d_amp_row.E, max_d_amp_row.eta_h, max_d_amp_row.d_amp))
        println(io, @sprintf("- Max true m_bar occurs at E = %.6f, eta_h = %.6f with m_bar = %.6f.", max_m_row.E, max_m_row.eta_h, max_m_row.m_bar))
        println(io, @sprintf("- Max capbind_frac occurs at E = %.6f, eta_h = %.6f with capbind_frac = %.6f.", max_capbind_row.E, max_capbind_row.eta_h, max_capbind_row.capbind_frac))
        println(io, @sprintf("- Max return_time_cv occurs at E = %.6f, eta_h = %.6f with CV = %.6f.", max_cv_row.E, max_cv_row.eta_h, max_cv_row.return_time_cv))
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "- The main 3D signature should be read against `d_amp`, `return_time_cv`, `true m_bar`, and `capbind_frac`, not against the old mid-band proxy alone.")
        if max_cv_row.return_time_cv < 0.05
            println(io, "- Return-time variation stays low across the scan, which supports the reading of 3D as a persistent large cycle rather than a convincing MMO regime.")
        else
            println(io, "- Some localized return-time irregularity appears, but it is not enough on its own to establish MMO.")
        end
        if max_m_row.m_bar < 0.01
            println(io, "- The true middle band remains small throughout most of the scan, so the oscillation is better interpreted as a near-cap relaxation cycle than as wide mixed-mode roaming through a middle band.")
        else
            println(io, "- The true middle band becomes nontrivial in parts of the scan, but that does not automatically imply MMO.")
        end
        if isempty(suspects)
            println(io, "- No local region crosses even a weak irregularity threshold; current evidence against strong MMO remains intact.")
        else
            println(io, "- Local suspicious regions worth keeping separate from a hard MMO claim:")
            for row in suspects
                println(io, @sprintf("  - E = %.6f, eta_h = %.6f, return_time_cv = %.6f, crossing_ratio = %.6f, true m_bar = %.6f, mid_frac_proxy = %.6f",
                    row.E, row.eta_h, row.return_time_cv, row.crossing_ratio, row.m_bar, row.mid_frac_proxy))
            end
        end
    end
    return path
end

function main()
    mode = parse_mode()
    data = run_persistence(mode)
    csv_path = joinpath(NEXT_PHASE_RESULTS_DIR, "threeD_persistence.csv")
    write_csv(csv_path, data.rows, [
        :E, :eta_h,
        :d_amp, :q_amp, :h_amp,
        :return_time_mean, :return_time_cv, :crossing_ratio,
        :m_bar, :capbind_frac, :phi_bar, :mid_frac_proxy,
    ])

    plot_path = joinpath(NEXT_PHASE_RESULTS_DIR, "threeD_persistence_map.png")
    save_persistence_plot(plot_path, data)

    summary_path = joinpath(NEXT_PHASE_RESULTS_DIR, "threeD_persistence_summary.md")
    write_persistence_summary(summary_path, mode, data)

    print_generated_files([csv_path, plot_path, summary_path])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
