using Printf
using Statistics

include(joinpath(@__DIR__, "next_phase_utils.jl"))

function frontier_metrics_config(mode::Symbol, dim::Symbol)
    if dim === :twoD
        return mode === :full ? (tmax = 9000.0, transient = 4500.0, saveat = 0.1) : (tmax = 6000.0, transient = 3000.0, saveat = 0.1)
    elseif dim === :threeD
        return mode === :full ? (tmax = 9000.0, transient = 4500.0, saveat = 0.2) : (tmax = 6000.0, transient = 3000.0, saveat = 0.2)
    else
        return mode === :full ? (tmax = 36000.0, transient = 18000.0, saveat = 0.25) : (tmax = 18000.0, transient = 9000.0, saveat = 0.25)
    end
end

function tail_indices(t; tail_width = 600.0)
    tmax = t[end]
    mask = t .>= max(t[1], tmax - tail_width)
    return findall(identity, mask)
end

function save_frontier_timeseries(path, point, sim)
    fm = frontier_metric_series(sim, point.p)
    idx = tail_indices(sim.t; tail_width = point.dim === :full ? 800.0 : 200.0)
    t = sim.t[idx]
    d = sim.d[idx]
    a = fm.a_t[idx]
    r = fm.r_t[idx]
    m = fm.m_t[idx]
    capgap = fm.capgap_t[idx]
    phi = fm.phi_t[idx]

    p1 = plot(t, d; lw = 2, xlabel = "t", ylabel = "frontier level", title = "$(point.label): d(t), a(t), r(t)", label = "d")
    plot!(p1, t, a; lw = 2, label = "a")
    plot!(p1, t, r; lw = 2, label = "r")

    p2 = plot(t, m; lw = 2, xlabel = "t", ylabel = "diagnostic", title = "$(point.label): m(t), capgap(t), phi(t)", label = "m")
    plot!(p2, t, capgap; lw = 2, label = "capgap")
    plot!(p2, t, phi; lw = 2, label = "phi")

    combo = plot(p1, p2; layout = (2, 1), size = (1200, 900))
    savefig(combo, path)
    return path
end

function save_frontier_hist(path, point, sim)
    fm = frontier_metric_series(sim, point.p)
    p1 = histogram(fm.m_t; bins = 50, xlabel = "m_t", ylabel = "count", title = "$(point.label): true middle band", label = false)
    p2 = histogram(fm.capgap_t; bins = 50, xlabel = "capgap_t", ylabel = "count", title = "$(point.label): advisor-cap gap", label = false)
    combo = plot(p1, p2; layout = (1, 2), size = (1200, 450))
    savefig(combo, path)
    return path
end

function summary_row(point, sim)
    fm = frontier_metric_summary(sim, point.p)
    return merge(
        (
            name = point.name,
            label = point.label,
            dim = String(point.dim),
            source = point.source,
        ),
        ordered_param_row(point.p),
        fm,
    )
end

function write_frontier_summary(path, mode, rows, points)
    lookup = Dict(row.name => row for row in rows)
    ridge_names = ["ridge_peak", "ridge_inner", "ridge_end", "boundary"]

    capbind_vals = [lookup[name].capbind_frac for name in ridge_names]
    mbar_vals = [lookup[name].m_bar for name in ridge_names]
    capgap_vals = [lookup[name].capgap_mean for name in ridge_names]
    ridge_capbinding = mean(capbind_vals)
    ridge_mbar = mean(mbar_vals)

    strongest_capbind = ridge_names[argmax(capbind_vals)]
    weakest_capbind = ridge_names[argmin(capbind_vals)]
    strongest_m = ridge_names[argmax(mbar_vals)]
    weakest_m = ridge_names[argmin(mbar_vals)]

    open(path, "w") do io
        println(io, "# True Double-Frontier Diagnostics")
        println(io)
        println(io, "- Mode: `$(mode_label(mode))`")
        println(io, "- Script: `scripts/next_phase_frontier_metrics.jl`")
        println(io, "- Output table: `results/next_phase/frontier_metrics_table.csv`")
        println(io)
        println(io, "## Parameter Points")
        println(io)
        for point in points
            println(io, "### $(point.label)")
            println(io)
            println(io, "- source: $(point.source)")
            println(io, "- dimension: `$(point.dim)`")
            for line in parameter_lines(point.p)
                println(io, line)
            end
            println(io)
        end
        println(io, "## Metric Table")
        println(io)
        println(io, "| point | dim | m_bar | m_max | capbind_frac | capgap_mean | phi_bar | d_mean | q_mean | k_mean | h_mean | m_q05 | m_q50 | m_q95 |")
        println(io, "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for row in rows
            println(io, @sprintf("| %s | %s | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f |",
                row.label, row.dim, row.m_bar, row.m_max, row.capbind_frac, row.capgap_mean, row.phi_bar,
                row.d_mean, row.q_mean, row.k_mean, row.h_mean, row.m_q05, row.m_q50, row.m_q95))
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, @sprintf("- Across the four 4D ridge points, average cap-binding frequency is %.6f and average true middle-band width is %.6f.", ridge_capbinding, ridge_mbar))
        println(io, @sprintf("- The strongest cap-binding point is `%s` with capbind_frac = %.6f; the weakest is `%s` with capbind_frac = %.6f.",
            strongest_capbind, lookup[strongest_capbind].capbind_frac, weakest_capbind, lookup[weakest_capbind].capbind_frac))
        println(io, @sprintf("- The largest true middle band is at `%s` with m_bar = %.6f; the smallest is `%s` with m_bar = %.6f.",
            strongest_m, lookup[strongest_m].m_bar, weakest_m, lookup[weakest_m].m_bar))
        println(io, @sprintf("- The ridge-inner point has capbind_frac = %.6f, capgap_mean = %.6f, and m_bar = %.6f.", lookup["ridge_inner"].capbind_frac, lookup["ridge_inner"].capgap_mean, lookup["ridge_inner"].m_bar))
        println(io, @sprintf("- The boundary point has capbind_frac = %.6f, capgap_mean = %.6f, and m_bar = %.6f.", lookup["boundary"].capbind_frac, lookup["boundary"].capgap_mean, lookup["boundary"].m_bar))
        println(io)
        if ridge_capbinding < 0.1
            println(io, "- Provisional answer: the 4D ridge does not spend much time in a strongly cap-binding regime; the advisor cap is active only intermittently in these samples.")
        else
            println(io, "- Provisional answer: the 4D ridge does spend substantial time in a cap-binding regime, so advisor-cap geometry remains dynamically relevant.")
        end
        if ridge_mbar < 1e-3
            println(io, "- The true middle band is numerically tiny across the ridge, which weakens a literal double-frontier interpretation.")
        else
            println(io, "- The true middle band remains nontrivial across the ridge, so the cap geometry does not collapse to zero-width.")
        end
        if mean(capgap_vals) < 1e-3
            println(io, "- Mean cap-gap is very small, suggesting the raw doer target rarely pushes far past the advisor frontier.")
        else
            println(io, "- Mean cap-gap stays visibly positive, so the raw doer target frequently wants to exceed the advisor frontier before capping.")
        end
    end
    return path
end

function main()
    mode = parse_mode()
    points = representative_parameter_sets()
    rows = NamedTuple[]
    sims = Dict{String, Any}()

    for point in points
        cfg = frontier_metrics_config(mode, point.dim)
        sim = simulate_point(point; cfg...)
        sims[point.name] = sim
        push!(rows, summary_row(point, sim))
    end

    csv_path = joinpath(NEXT_PHASE_RESULTS_DIR, "frontier_metrics_table.csv")
    columns = [
        :name, :label, :dim, :source,
        default_param_order()...,
        :m_bar, :m_max, :m_q05, :m_q50, :m_q95,
        :capbind_frac, :capgap_mean, :phi_bar,
        :d_mean, :q_mean, :k_mean, :h_mean,
    ]
    write_csv(csv_path, rows, Symbol.(columns))

    generated = String[csv_path]

    timeseries_specs = [
        ("base_default", "frontier_timeseries_baseline.png"),
        ("ridge_peak", "frontier_timeseries_ridge_peak.png"),
        ("ridge_inner", "frontier_timeseries_ridge_inner.png"),
        ("boundary", "frontier_timeseries_boundary.png"),
    ]
    for (name, file) in timeseries_specs
        path = joinpath(NEXT_PHASE_RESULTS_DIR, file)
        save_frontier_timeseries(path, parameter_point(name), sims[name])
        push!(generated, path)
    end

    hist_specs = [
        ("ridge_peak", "frontier_hist_ridge_peak.png"),
        ("ridge_inner", "frontier_hist_ridge_inner.png"),
        ("boundary", "frontier_hist_boundary.png"),
    ]
    for (name, file) in hist_specs
        path = joinpath(NEXT_PHASE_RESULTS_DIR, file)
        save_frontier_hist(path, parameter_point(name), sims[name])
        push!(generated, path)
    end

    summary_path = joinpath(NEXT_PHASE_RESULTS_DIR, "frontier_metrics_summary.md")
    write_frontier_summary(summary_path, mode, rows, points)
    push!(generated, summary_path)

    print_generated_files(generated)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
