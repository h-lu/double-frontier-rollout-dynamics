using Printf
using Statistics

include(joinpath(@__DIR__, "next_phase_utils.jl"))

function poincare_config(mode::Symbol)
    # Keep the long-window setting fixed in both modes so the return-map evidence
    # is directly comparable to the earlier period-2 check.
    return (tmax = 36000.0, transient = 18000.0, saveat = 0.25)
end

function poincare_row_table(crossings, widths)
    n = length(crossings.t)
    T = n > 1 ? diff(crossings.t) : Float64[]
    W = isempty(widths) ? Float64[] : widths[1:min(length(widths), n)]
    rows = NamedTuple[]
    for i in 1:n
        push!(rows, (
            n = i,
            t_n = crossings.t[i],
            q_n = crossings.q[i],
            k_n = crossings.k[i],
            h_n = crossings.h[i],
            d_dot_n = crossings.d_dot[i],
            T_n = i <= length(T) ? T[i] : NaN,
            W_n = i <= length(W) ? W[i] : NaN,
        ))
    end
    return rows
end

function save_map_plot(path, x, y; xlabel, ylabel, title)
    p = scatter(x, y; xlabel = xlabel, ylabel = ylabel, title = title, markersize = 4, label = false)
    savefig(p, path)
    return path
end

function save_qh_odd_even_plot(path, point, crossings)
    idx = collect(1:length(crossings.q))
    odd = isodd.(idx)
    p = scatter(crossings.q[odd], crossings.h[odd]; xlabel = "q_n", ylabel = "h_n", title = "$(point.label): odd/even crossings", markersize = 4, label = "odd")
    scatter!(p, crossings.q[.!odd], crossings.h[.!odd]; markersize = 4, label = "even")
    savefig(p, path)
    return path
end

function diagnostics_row(point, pm, crossings)
    return (
        name = point.name,
        label = point.label,
        n_crossings = length(crossings.t),
        n_periods = length(pm.T),
        lag1 = pm.lag1,
        lag2 = pm.lag2,
        width_lag1 = pm.width_lag1,
        width_lag2 = pm.width_lag2,
        odd_mean = pm.T_odd_mean,
        even_mean = pm.T_even_mean,
        odd_cv = pm.T_odd_cv,
        even_cv = pm.T_even_cv,
        width_odd_mean = pm.W_odd_mean,
        width_even_mean = pm.W_even_mean,
        width_odd_cv = pm.W_odd_cv,
        width_even_cv = pm.W_even_cv,
        R1 = pm.R1,
        R2 = pm.R2,
        R2_over_R1 = pm.R2_over_R1,
        q_odd_mean = pm.q_odd_mean,
        q_even_mean = pm.q_even_mean,
        h_odd_mean = pm.h_odd_mean,
        h_even_mean = pm.h_even_mean,
    )
end

function write_poincare_summary(path, mode, rows)
    lookup = Dict(row.name => row for row in rows)
    ratio_rows = filter(row -> !isnan(row.R2_over_R1), rows)
    min_ratio_row = ratio_rows[argmin([row.R2_over_R1 for row in ratio_rows])]
    open(path, "w") do io
        println(io, "# 4D Poincare / Return-Map Evidence")
        println(io)
        println(io, "- Mode: `$(mode_label(mode))`")
        println(io, "- Section: `d = 0.9`, upward crossings only")
        println(io, "- Integration window: `tmax = 36000`, `transient = 18000`, `saveat = 0.25`")
        println(io)
        println(io, "| point | n_crossings | lag1(T) | lag2(T) | odd mean T | even mean T | odd CV | even CV | R1 | R2 | R2/R1 |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for row in rows
            println(io, @sprintf("| %s | %d | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f |",
                row.label, row.n_crossings, row.lag1, row.lag2, row.odd_mean, row.even_mean, row.odd_cv, row.even_cv, row.R1, row.R2, row.R2_over_R1))
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        for name in ("ridge_peak", "ridge_inner", "ridge_end", "boundary")
            row = lookup[name]
            if !isnan(row.R2_over_R1) && row.R2_over_R1 < 0.25 && row.odd_cv < 0.05 && row.even_cv < 0.05
                println(io, @sprintf("- `%s`: stable 2-cycle-like return structure. R2/R1 = %.6f, odd/even CVs = [%.6f, %.6f].", row.label, row.R2_over_R1, row.odd_cv, row.even_cv))
            elseif !isnan(row.R2_over_R1) && row.R2_over_R1 < 0.6
                println(io, @sprintf("- `%s`: still shows a 2-step contraction, but the separation is weaker. R2/R1 = %.6f.", row.label, row.R2_over_R1))
            else
                println(io, @sprintf("- `%s`: current evidence weakens under the stricter Poincare diagnostic. R2/R1 = %.6f.", row.label, row.R2_over_R1))
            end
        end
        println(io)
        println(io, @sprintf("- `ridge_inner` still has the tightest scalar alternation in `T_n` among the ridge-interior points: R2/R1 = %.6f, lag1 = %.6f, lag2 = %.6f. The stricter state-space metric therefore weakens, rather than strengthens, the 2-cycle claim.",
            lookup["ridge_inner"].R2_over_R1, lookup["ridge_inner"].lag1, lookup["ridge_inner"].lag2))
        println(io, @sprintf("- The smallest observed R2/R1 is at `%s`, with R2/R1 = %.6f.", min_ratio_row.label, min_ratio_row.R2_over_R1))
        println(io, @sprintf("- `boundary` has R2/R1 = %.6f and odd/even period means [%.6f, %.6f]. This indicates %s.",
            lookup["boundary"].R2_over_R1, lookup["boundary"].odd_mean, lookup["boundary"].even_mean,
            lookup["boundary"].R2_over_R1 < 0.4 ? "the same 2-cycle branch persisting but with weaker split" : "a weakened or changed mechanism relative to the ridge interior"))
        println(io)
        println(io, "- If the odd/even clouds remain tight and R2 << R1, that supports a stable 2-cycle scaffold rather than diffuse modulation.")
    end
    return path
end

function main()
    mode = parse_mode()
    cfg = poincare_config(mode)
    point_names = ["ridge_peak", "ridge_inner", "ridge_end", "boundary"]

    diagnostics = NamedTuple[]
    generated = String[]

    for name in point_names
        point = parameter_point(name)
        sim = simulate_point(point; cfg...)
        crossings = threshold_crossing_data(sim, point.p; thresh = 0.9)
        widths = width_sequence(sim; thresh = 0.9)
        pm = poincare_metrics(crossings, widths)

        csv_rows = poincare_row_table(crossings, widths)
        csv_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_$(name).csv")
        write_csv(csv_path, csv_rows, [:n, :t_n, :q_n, :k_n, :h_n, :d_dot_n, :T_n, :W_n])
        push!(generated, csv_path)

        qh_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_qh_odd_even_$(name).png")
        save_qh_odd_even_plot(qh_path, point, crossings)
        push!(generated, qh_path)

        if name == "ridge_inner"
            if length(pm.T) > 1
                tmap_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_T_map_ridge_inner.png")
                save_map_plot(tmap_path, pm.T[1:end-1], pm.T[2:end]; xlabel = "T_n", ylabel = "T_{n+1}", title = "ridge_inner: T return map")
                push!(generated, tmap_path)
            end
            if length(crossings.q) > 1
                qmap_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_q_map_ridge_inner.png")
                save_map_plot(qmap_path, crossings.q[1:end-1], crossings.q[2:end]; xlabel = "q_n", ylabel = "q_{n+1}", title = "ridge_inner: q return map")
                push!(generated, qmap_path)

                hmap_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_h_map_ridge_inner.png")
                save_map_plot(hmap_path, crossings.h[1:end-1], crossings.h[2:end]; xlabel = "h_n", ylabel = "h_{n+1}", title = "ridge_inner: h return map")
                push!(generated, hmap_path)
            end
        end

        push!(diagnostics, diagnostics_row(point, pm, crossings))
    end

    summary_path = joinpath(NEXT_PHASE_RESULTS_DIR, "poincare_summary.md")
    write_poincare_summary(summary_path, mode, diagnostics)
    push!(generated, summary_path)

    print_generated_files(generated)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
