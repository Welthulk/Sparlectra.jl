# file: examples/others/apslf_vs_nr_timing.jl
# purpose: one table "case, buses, PV buses, solver, order, outcome, iterations,
#          time" for APSLF against the rectangular Newton solver on the shipped
#          MATPOWER cases and on synthetic tiled grids of growing size.
#          Timing is the median of three warm runs; the first run of every
#          (case, solver) pair is discarded as compile/warm-up.
#
# Aufruf:  julia --startup-file=no --project=. examples/others/apslf_vs_nr_timing.jl
# Ausgabe: Tabelle auf der Konsole, apslf_vs_nr_timing.csv und apslf_vs_nr_timing.svg
#          (Zeit über Knotenzahl, log-log, eine Linie je Solver/Ordnung) im
#          aktuellen Verzeichnis. Nur konvergierte Läufe landen in der Kurve.
# Optional: SPARLECTRA_TIMING_SIZES="500,1000,2000,5000" for the tiled grids,
#           SPARLECTRA_TIMING_ORDERS="24,40,60" for the APSLF orders on the
#           larger cases (every case gets order 24; from 300 buses on, the
#           further orders are added).

using Sparlectra
using Printf
using Statistics

const MPOWER = joinpath(pkgdir(Sparlectra), "data", "mpower")
const SIZES = parse.(Int, split(get(ENV, "SPARLECTRA_TIMING_SIZES", "500,1000,2000,5000"), ","))
const ORDERS = parse.(Int, split(get(ENV, "SPARLECTRA_TIMING_ORDERS", "24,40,60"), ","))
const REPEATS = 3

quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
cfg_nr = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false), output=quiet)
cfg_ap(order) = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=order, convergence_radius=false)), output=quiet)

# the MATPOWER import options of the workshops
import_m(path) = createNetFromMatPowerFile(filename=path, flatstart=false, enable_pq_gen_controllers=true, bus_shunt_model=:admittance, matpower_shift_sign=1.0, matpower_shift_unit=:deg, matpower_ratio=:normal, tap_changer_model=:ideal)

# every entry is (label, builder) where the builder returns a FRESH net; the
# solve mutates the net, so each timed run starts from the same start state
cases = Tuple{String,Function}[]
for name in ("sp_case9", "sp_case118", "sp_case300", "sp_case1354")
    path = joinpath(MPOWER, name * ".m")
    isfile(path) || (println("skip ", name, ": ", path, " not found"); continue)
    push!(cases, (name, () -> import_m(path)))
end
push!(cases, ("sp_case5", () -> importSCF(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case5.scf.json"))))
for n in SIZES
    push!(cases, ("tiled_$(n)", () -> first(build_synthetic_tiled_grid_net(n))))
end

function timed_run(build, cfg)
    r = run_sparlectra(net=build(), config=cfg)          # warm-up, discarded
    ts = Float64[]
    for _ in 1:REPEATS
        net = build()
        push!(ts, @elapsed r = run_sparlectra(net=net, config=cfg))
    end
    return r, median(ts)
end

rows = NamedTuple[]
for (label, build) in cases
    net0 = build()
    nb = length(net0.nodeVec)
    npv = count(isPVNode, net0.nodeVec)
    local runs = [("NR", 0, cfg_nr)]
    orders = nb >= 300 ? ORDERS : (first(ORDERS),)
    for o in orders
        push!(runs, ("APSLF", o, cfg_ap(o)))
    end
    for (solver, order, cfg) in runs
        r, t = try
            timed_run(build, cfg)
        catch err
            println("  ", label, " ", solver, " order ", order, ": ", sprint(showerror, err)[1:min(end, 120)])
            (nothing, NaN)
        end
        push!(rows, (case=label, buses=nb, pv=npv, solver=solver, order=order,
            outcome=r === nothing ? "error" : string(r.outcome),
            iterations=r === nothing ? 0 : r.iterations,
            mismatch=r === nothing ? NaN : r.final_mismatch, ms=1000 * t))
        @printf("  %-12s %5d buses  %-5s order %2d  %-14s %3d it  %.1e  %8.1f ms\n", label, nb, solver, order,
            rows[end].outcome, rows[end].iterations, rows[end].mismatch, rows[end].ms)
    end
end

# ---- console table, fixed width ------------------------------------------
println()
hdr = @sprintf("%-12s %6s %4s  %-6s %5s  %-14s %4s  %-8s  %10s", "case", "buses", "PV", "solver", "order", "outcome", "it", "mismatch", "time")
println(hdr)
println("-" ^ length(hdr))
for r in rows
    @printf("%-12s %6d %4d  %-6s %5s  %-14s %4d  %-8.1e  %8.1f ms\n", r.case, r.buses, r.pv, r.solver,
        r.order == 0 ? "-" : string(r.order), r.outcome, r.iterations, r.mismatch, r.ms)
end
println()
println("threads = ", Threads.nthreads(), ", Julia ", VERSION, ", median of ", REPEATS, " warm runs, APSLF without convergence-radius evaluation")

# ---- CSV ------------------------------------------------------------------
open("apslf_vs_nr_timing.csv", "w") do io
    println(io, "case,buses,pv_buses,solver,order,outcome,iterations,mismatch,time_ms")
    for r in rows
        println(io, join((r.case, r.buses, r.pv, r.solver, r.order, r.outcome, r.iterations, r.mismatch, round(r.ms; digits=2)), ","))
    end
end

# ---- SVG: time over bus count, log-log, converged runs only -----------------
function write_svg(path, rows)
    series = Dict{String,Vector{Tuple{Float64,Float64}}}()
    for r in rows
        r.outcome == "converged" || continue
        key = r.solver == "NR" ? "NR" : "APSLF $(r.order)"
        push!(get!(series, key, Tuple{Float64,Float64}[]), (Float64(r.buses), r.ms))
    end
    isempty(series) && return
    pts = vcat(values(series)...)
    xmin, xmax = extrema(first.(pts))
    ymin, ymax = extrema(last.(pts))
    lx(v) = log10(v)
    W, H, L, B = 720.0, 420.0, 70.0, 50.0
    sx(v) = L + (lx(v) - lx(xmin)) / max(lx(xmax) - lx(xmin), 1e-9) * (W - L - 20)
    sy(v) = H - B - (lx(v) - lx(ymin)) / max(lx(ymax) - lx(ymin), 1e-9) * (H - B - 30)
    colors = ["#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b"]
    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(W)\" height=\"$(H)\" font-family=\"sans-serif\" font-size=\"12\">")
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
        println(io, "<text x=\"$(W/2)\" y=\"18\" text-anchor=\"middle\" font-size=\"14\">Solve time over bus count (log-log), converged runs, Julia $(VERSION)</text>")
        # axes with decade ticks
        for d in floor(Int, lx(xmin)):ceil(Int, lx(xmax))
            x = sx(10.0^d)
            xmin <= 10.0^d <= xmax || continue
            println(io, "<line x1=\"$(x)\" y1=\"30\" x2=\"$(x)\" y2=\"$(H-B)\" stroke=\"#ddd\"/>")
            println(io, "<text x=\"$(x)\" y=\"$(H-B+16)\" text-anchor=\"middle\">$(10^d)</text>")
        end
        for d in floor(Int, lx(ymin)):ceil(Int, lx(ymax))
            y = sy(10.0^d)
            ymin <= 10.0^d <= ymax || continue
            println(io, "<line x1=\"$(L)\" y1=\"$(y)\" x2=\"$(W-20)\" y2=\"$(y)\" stroke=\"#ddd\"/>")
            println(io, "<text x=\"$(L-6)\" y=\"$(y+4)\" text-anchor=\"end\">$(10^d)</text>")
        end
        println(io, "<line x1=\"$(L)\" y1=\"$(H-B)\" x2=\"$(W-20)\" y2=\"$(H-B)\" stroke=\"#333\"/>")
        println(io, "<line x1=\"$(L)\" y1=\"30\" x2=\"$(L)\" y2=\"$(H-B)\" stroke=\"#333\"/>")
        println(io, "<text x=\"$(W/2)\" y=\"$(H-8)\" text-anchor=\"middle\">buses</text>")
        println(io, "<text x=\"14\" y=\"$(H/2)\" text-anchor=\"middle\" transform=\"rotate(-90 14 $(H/2))\">ms</text>")
        for (k, name) in enumerate(sort(collect(keys(series))))
            p = sort(series[name]; by=first)
            c = colors[mod1(k, length(colors))]
            d = join(("$(i == 1 ? "M" : "L")$(round(sx(x); digits=1)),$(round(sy(y); digits=1))" for (i, (x, y)) in enumerate(p)), " ")
            println(io, "<path d=\"$(d)\" fill=\"none\" stroke=\"$(c)\" stroke-width=\"2\"/>")
            for (x, y) in p
                println(io, "<circle cx=\"$(round(sx(x); digits=1))\" cy=\"$(round(sy(y); digits=1))\" r=\"3\" fill=\"$(c)\"/>")
            end
            println(io, "<rect x=\"$(L+10)\" y=\"$(36+16*k)\" width=\"12\" height=\"12\" fill=\"$(c)\"/>")
            println(io, "<text x=\"$(L+28)\" y=\"$(46+16*k)\">$(name)</text>")
        end
        println(io, "</svg>")
    end
end
write_svg("apslf_vs_nr_timing.svg", rows)
println("written apslf_vs_nr_timing.csv and apslf_vs_nr_timing.svg")