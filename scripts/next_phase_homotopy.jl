using Printf
using Statistics

include(joinpath(@__DIR__, "next_phase_utils.jl"))

function homotopy_specs(mode::Symbol)
    npts = mode === :full ? 61 : 31
    return [
        (name = "beta_n", label = "beta_n", values = collect(range(0.7, 1.5, length = npts))),
        (name = "lambda_D", label = "lambda_D", values = collect(range(12.0, 20.0, length = npts))),
        (name = "theta0", label = "theta0", values = collect(range(-2.0, -5.0, length = npts))),
        (name = "eta_h", label = "eta_h", values = collect(range(0.01, 0.002, length = npts))),
        (name = "eta_k", label = "eta_k", values = collect(range(0.03, 0.003, length = npts))),
    ]
end

function homotopy_sim_config(mode::Symbol)
    if mode === :full
        return (
            bistability = (tmax = 3500.0, transient = 1800.0, saveat = 0.5),
            sim4d = (tmax = 12000.0, transient = 6000.0, saveat = 0.5),
        )
    else
        return (
            bistability = (tmax = 2500.0, transient = 1200.0, saveat = 0.5),
            sim4d = (tmax = 8000.0, transient = 4000.0, saveat = 0.5),
        )
    end
end

function apply_single_parameter(base, name::String, value)
    return merge(base, (Symbol(name) => value,))
end

function homotopy_row(spec, value, p, cfg)
    bistable = cheap_2d_bistability_proxy(p; cfg.bistability...)
    diag4d = cheap_4d_summary(p; cfg.sim4d...)
    return (
        homotopy = spec.name,
        parameter = spec.label,
        value = value,
        theoretical_fold_viable = p.beta_n * p.lambda_D > 4 ? 1 : 0,
        twoD_bistable_proxy = bistable.has_bistability ? 1 : 0,
        twoD_gap = bistable.gap,
        sustained_oscillation = diag4d.sustained_oscillation ? 1 : 0,
        m_bar = diag4d.frontier.m_bar,
        capbind_frac = diag4d.frontier.capbind_frac,
        capgap_mean = diag4d.frontier.capgap_mean,
        phi_bar = diag4d.frontier.phi_bar,
        period2_strength = diag4d.period2_strength,
        R2_over_R1 = diag4d.poincare.R2_over_R1,
        d_amp = diag4d.d_amp,
        q_amp = diag4d.q_amp,
        k_amp = diag4d.k_amp,
        h_amp = diag4d.h_amp,
        E = p.E,
        G = p.G,
        epsilon = p.epsilon,
        eta_h = p.eta_h,
        eta_k = p.eta_k,
        beta_n = p.beta_n,
        lambda_D = p.lambda_D,
        theta0 = p.theta0,
    )
end

function run_single_homotopy(spec, mode)
    base = base_params()
    cfg = homotopy_sim_config(mode)
    rows = Vector{NamedTuple}(undef, length(spec.values))
    @threads for idx in eachindex(spec.values)
        value = spec.values[idx]
        p = apply_single_parameter(base, spec.name, value)
        rows[idx] = homotopy_row(spec, value, p, cfg)
    end
    return rows
end

function switch_points(rows, field::Symbol)
    vals = [getproperty(row, field) for row in rows]
    out = String[]
    for i in 2:length(rows)
        if vals[i] != vals[i - 1]
            push!(out, @sprintf("%.6f -> %.6f", rows[i - 1].value, rows[i].value))
        end
    end
    return out
end

function save_homotopy_plot(path, specs, all_rows)
    panels = Any[]
    for spec in specs
        rows = all_rows[spec.name]
        x = [row.value for row in rows]
        push!(panels, plot(x, [row.d_amp for row in rows]; lw = 2, xlabel = spec.label, ylabel = "d_amp", title = "$(spec.label): d_amp", label = false))
        push!(panels, plot(x, [row.capbind_frac for row in rows]; lw = 2, xlabel = spec.label, ylabel = "capbind_frac", title = "$(spec.label): capbind_frac", label = false))
        push!(panels, plot(x, [row.m_bar for row in rows]; lw = 2, xlabel = spec.label, ylabel = "m_bar", title = "$(spec.label): m_bar", label = false))
        push!(panels, plot(x, [row.period2_strength for row in rows]; lw = 2, xlabel = spec.label, ylabel = "period2_strength", title = "$(spec.label): period2_strength", label = false))
    end
    combo = plot(panels...; layout = (length(specs), 4), size = (1600, 2200))
    savefig(combo, path)
    return path
end

function write_homotopy_summary(path, mode, specs, all_rows)
    open(path, "w") do io
        println(io, "# Default-Box to Ridge Homotopy")
        println(io)
        println(io, "- Mode: `$(mode_label(mode))`")
        println(io, "- Base point: `base_params()` with only one parameter varied at a time")
        println(io, "- Output figure: `results/next_phase/homotopy_regime_map.png`")
        println(io)
        println(io, "## Per-Homotopy Findings")
        println(io)
        for spec in specs
            rows = all_rows[spec.name]
            fold_switch = switch_points(rows, :theoretical_fold_viable)
            bistable_switch = switch_points(rows, :twoD_bistable_proxy)
            osc_switch = switch_points(rows, :sustained_oscillation)
            max_d_amp = maximum(row.d_amp for row in rows)
            max_capbind = maximum(row.capbind_frac for row in rows)
            max_mbar = maximum(row.m_bar for row in rows)
            max_p2 = maximum(filter(!isnan, [row.period2_strength for row in rows]); init = NaN)
            println(io, "### $(spec.label)")
            println(io)
            println(io, @sprintf("- max d_amp = %.6f", max_d_amp))
            println(io, @sprintf("- max capbind_frac = %.6f", max_capbind))
            println(io, @sprintf("- max true m_bar = %.6f", max_mbar))
            println(io, @sprintf("- max period2_strength = %.6f", max_p2))
            println(io, "- theoretical fold viability switches: $(isempty(fold_switch) ? "none" : join(fold_switch, ", "))")
            println(io, "- 2D bistability proxy switches: $(isempty(bistable_switch) ? "none" : join(bistable_switch, ", "))")
            println(io, "- 4D sustained-oscillation proxy switches: $(isempty(osc_switch) ? "none" : join(osc_switch, ", "))")
            println(io)
        end
        println(io, "## Overall Conclusion")
        println(io)
        any_bistable = any(any(row.twoD_bistable_proxy == 1 for row in rows) for rows in values(all_rows))
        any_osc = any(any(row.sustained_oscillation == 1 for row in rows) for rows in values(all_rows))
        peak_capbind = maximum(maximum(row.capbind_frac for row in rows) for rows in values(all_rows))
        peak_mbar = maximum(maximum(row.m_bar for row in rows) for rows in values(all_rows))
        println(io, any_bistable ?
            "- A single-parameter deformation can already trigger some low-dimensional regime changes." :
            "- None of the five single-parameter homotopies produce a strong 2D bistability signal inside the default box.")
        println(io, any_osc ?
            "- At least one single-parameter homotopy generates sustained 4D oscillation before the full ridge parameter package is assembled." :
            "- Sustained 4D oscillation does not appear under these single-parameter deformations, which means the modulation ridge lies well outside the default box in a multivariate sense.")
        println(io, @sprintf("- Peak capbind_frac across all homotopies is %.6f and peak true m_bar is %.6f.", peak_capbind, peak_mbar))
        println(io, "- This homotopy result should be read as a distance test, not as a proof of impossibility.")
    end
    return path
end

function main()
    mode = parse_mode()
    specs = homotopy_specs(mode)
    all_rows = Dict{String, Vector{NamedTuple}}()
    generated = String[]

    for spec in specs
        rows = run_single_homotopy(spec, mode)
        all_rows[spec.name] = rows
        csv_path = joinpath(NEXT_PHASE_RESULTS_DIR, "homotopy_$(spec.name).csv")
        write_csv(csv_path, rows, [
            :homotopy, :parameter, :value, :theoretical_fold_viable,
            :twoD_bistable_proxy, :twoD_gap, :sustained_oscillation,
            :m_bar, :capbind_frac, :capgap_mean, :phi_bar,
            :period2_strength, :R2_over_R1,
            :d_amp, :q_amp, :k_amp, :h_amp,
            :E, :G, :epsilon, :eta_h, :eta_k, :beta_n, :lambda_D, :theta0,
        ])
        push!(generated, csv_path)
    end

    plot_path = joinpath(NEXT_PHASE_RESULTS_DIR, "homotopy_regime_map.png")
    save_homotopy_plot(plot_path, specs, all_rows)
    push!(generated, plot_path)

    summary_path = joinpath(NEXT_PHASE_RESULTS_DIR, "homotopy_summary.md")
    write_homotopy_summary(summary_path, mode, specs, all_rows)
    push!(generated, summary_path)

    print_generated_files(generated)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
