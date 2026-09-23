# Copyright 2023-2026 Udo Schmitz                                             #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     http://www.apache.org/licenses/LICENSE-2.0                              #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and        #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_apslf.jl                                            #src
# purpose: Literate.jl source of the APSLF workshop notebook and its          #src
#          Documenter page: the analytic power-series solver used through     #src
#          Sparlectra's model and configuration, series order, the           #src
#          convergence radius as a loadability margin, Q-limits, the hybrid   #src
#          start. Regenerate the committed outputs with                       #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # APSLF in Sparlectra: series solver, convergence radius, hybrid start
#
# > **Level: Intermediate**, after the [basic tour](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour/). The theory of the analytic power-series load flow (APSLF) and the raw solver have their own workshop in [AnalyticLoadFlow.jl](https://github.com/Welthulk/AnalyticLoadFlow.jl); this chapter is the complement from the network side: the same solver used on a Sparlectra model, through the configuration, with the diagnostics a run reports.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Welthulk/Sparlectra.jl/blob/main/notebooks/workshop_apslf.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# Newton-Raphson answers the power-flow equations by iterating from a
# start value. APSLF answers them without one: the bus voltages are
# written as a power series $V(s)$ in an embedding parameter $s$, the
# coefficients follow from a recursion on the admittance matrix, and the
# solution is the series evaluated at $s = 1$ (through Padé approximants
# rather than a plain Taylor sum). Three things follow for daily work,
# and this chapter shows each of them on a network you can inspect:
#
# 1. **No start value, no divergence.** The series either converges at
#    $s = 1$ or it does not; there is no "wrong basin" and no damping to
#    tune (Part 1).
# 2. **A margin comes for free.** The distance of the nearest Padé pole
#    to $s = 1$ (`dmin`, the APSLF convergence radius) shrinks as the
#    network approaches its loadability limit, so one solve carries its
#    own stress indicator (Part 2).
# 3. **A start value for Newton.** Where Newton is still wanted (with
#    controllers, or for a polished result), the series solution is the
#    best start value there is (Part 3).
#
# Nothing here needs an extra installation: AnalyticLoadFlow.jl is a
# dependency of Sparlectra, and `power_flow.solver = apslf` in the
# configuration switches a run over to it.

#nb # > **Colab is slow here.** The install cell below fetches Sparlectra and
#nb # > precompiles it on Colab's two cores; that takes about ten minutes and
#nb # > prints nothing while it runs ("Precompiling packages..." is the last
#nb # > line you see). Wait for the cell to finish before running the next one;
#nb # > a second click on it starts the whole install again.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb ENV["SPARLECTRA_PRECOMPILE_WORKLOAD"] = "minimal"   # notebook: skip the service-layer warm-up, about 40 s less precompile
#nb Pkg.add(url = "https://github.com/Welthulk/Sparlectra.jl", rev = "main")
#nb ## The notebook installs the development version from GitHub.
#nb ## For the latest registered release use: Pkg.add("Sparlectra")
using Sparlectra
using Printf

# ## Part 1: the same network, two solvers
#
# The study network is the 7-bus ring of the basic tour (B1 carries the
# grid connection, B3 a generator, the diagonals are the cross-ties B2-B5
# and B3-B6). It has no controllers, which matters: APSLF runs the whole
# solve inside the series, so the outer control loop (taps, Q(U)) has no
# place in it; Sparlectra rejects such a network for the APSLF solver
# instead of silently ignoring the controllers (Part 3 shows that).

function build_ring7(name::String; lambda::Float64=1.0)
    net = Net(name=name, baseMVA=100.0)
    addBus!(net=net, busName="B1", vn_kV=110.0, vm_pu=1.02, va_deg=0.0)
    for i in 2:7
        addBus!(net=net, busName="B$(i)", vn_kV=110.0, vm_pu=1.0, va_deg=0.0)
    end
    addPIModelACLine!(net=net, fromBus="B1", toBus="B2", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B2", toBus="B3", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B3", toBus="B4", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B4", toBus="B5", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B5", toBus="B6", r_pu=0.011, x_pu=0.085, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B6", toBus="B7", r_pu=0.012, x_pu=0.090, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B7", toBus="B1", r_pu=0.010, x_pu=0.080, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B2", toBus="B5", r_pu=0.009, x_pu=0.070, b_pu=0.0, status=1)
    addPIModelACLine!(net=net, fromBus="B3", toBus="B6", r_pu=0.009, x_pu=0.070, b_pu=0.0, status=1)
    addProsumer!(net=net, busName="B1", type="EXTERNALNETWORKINJECTION", referencePri="B1", vm_pu=1.02, va_deg=0.0)
    addProsumer!(net=net, busName="B3", type="GENERATOR", p=60.0, q=10.0)
    addProsumer!(net=net, busName="B2", type="LOAD", p=lambda * 35.0, q=lambda * 10.0)
    addProsumer!(net=net, busName="B4", type="LOAD", p=lambda * 45.0, q=lambda * 15.0)
    addProsumer!(net=net, busName="B5", type="LOAD", p=lambda * 25.0, q=lambda * 8.0)
    addProsumer!(net=net, busName="B6", type="LOAD", p=lambda * 30.0, q=lambda * 10.0)
    addProsumer!(net=net, busName="B7", type="LOAD", p=lambda * 20.0, q=lambda * 6.0)
    ok, msg = validate!(net=net)
    ok || error("Network validation failed: $msg")
    return net
end

# Two configurations, identical except for the solver. `run_sparlectra`
# takes the whole run from one `SparlectraConfig`: solver choice, control
# loop, output. The console summary is switched off here because the
# chapter prints its own comparisons.

quiet = OutputConfig(logfile_results=:off, console_summary=false, startup_latency_hint=false)
cfg_nr = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, rescue=false), output=quiet)
cfg_apslf = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf), output=quiet)

res_nr = run_sparlectra(net=build_ring7("ring7 NR"), config=cfg_nr)
res_ap = run_sparlectra(net=build_ring7("ring7 APSLF"), config=cfg_apslf)
(res_nr.outcome, res_ap.outcome)

# Both runs report the same outcome symbol. The result object carries the
# solved net; the bus voltages of the two solves agree to solver
# tolerance, which is the first thing to check whenever a second solver
# enters a workflow:

bus_vm(net) = [n._vm_pu for n in net.nodeVec]
bus_va(net) = [n._va_deg for n in net.nodeVec]
max_dvm = maximum(abs.(bus_vm(res_nr.net) .- bus_vm(res_ap.net)))
max_dva = maximum(abs.(bus_va(res_nr.net) .- bus_va(res_ap.net)))
@printf("max |ΔVm| = %.2e pu, max |ΔVa| = %.2e deg\n", max_dvm, max_dva)
@assert max_dvm < 1e-6 && max_dva < 1e-4                                   #src

# What differs is how the two got there. `iterations` counts Newton steps
# for the rectangular solver; for APSLF it counts the outer passes of the
# solver (one for a plain solve, more when reactive limits switch a
# machine, see Part 2), because the series itself has no iterations. The
# APSLF run additionally reports its convergence radius:

apslf_status = Sparlectra.rectangular_pf_status(res_ap.net)
println("NR    : ", res_nr.iterations, " iterations, mismatch ", res_nr.final_mismatch)
println("APSLF : ", res_ap.iterations, " pass(es),   mismatch ", res_ap.final_mismatch)
println("APSLF : ", apslf_status.apslf_convergence_line)

# Reading aid: `dmin` is the distance of the nearest Padé pole to the
# evaluation point $s = 1$. A pole AT $s = 1$ would mean the series cannot
# represent the solution there (the network has no solution, or the
# approximant is defective); the further away the pole, the more
# comfortably the solution sits inside the convergence region. The level
# (`GRN`/`YEL`/`RED`) is AnalyticLoadFlow's traffic light on that distance.
# Part 2 shows what moves it.
#
# ### The series order
#
# The one numerical knob of the solver is the number of series
# coefficients (`power_flow.apslf.order`, default 24). Too few and the
# Padé approximant has not settled; more than needed costs time and
# nothing else. The mismatch after the solve shows where the plateau is:

order_table = Tuple{Int,Float64,Bool}[]
for order in (4, 6, 8, 12, 16, 24, 40)
    cfg = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=order)), output=quiet)
    r = run_sparlectra(net=build_ring7("ring7 order $(order)"), config=cfg)
    push!(order_table, (order, r.final_mismatch, r.final_converged))
end
println(" order   max mismatch [pu]   converged")
for (order, mismatch, ok) in order_table
    @printf("  %3d    %.3e          %s\n", order, mismatch, ok)
end
@assert order_table[end][3]                                                 #src
@assert order_table[end][2] <= order_table[1][2]                            #src

# On a lightly loaded ring the plateau is reached early. Stressed
# networks (large angles, low voltages) need more coefficients, which is
# exactly the situation where `dmin` is small; the two knobs are read
# together in Part 2.

# ## Part 2: the convergence radius as a loadability margin
#
# Scaling every load of the ring by a factor $\lambda$ (the `lambda`
# keyword of `build_ring7`; loads are aggregated onto the buses when they
# are added, so the scaling happens at build time) walks the network
# towards its loadability limit (the "nose" of the PV curve). Newton
# reports that limit only by failing; APSLF reports it in advance
# through `dmin`, because the nearest Padé pole moves towards $s = 1$ as
# the solution branch approaches the fold.

scaled_ring7(lambda::Float64) = build_ring7(@sprintf("ring7 λ=%.2f", lambda); lambda=lambda)

margin_table = NamedTuple[]
for lambda in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
    r_ap = run_sparlectra(net=scaled_ring7(lambda), config=cfg_apslf)
    r_nr = run_sparlectra(net=scaled_ring7(lambda), config=cfg_nr)
    st = Sparlectra.rectangular_pf_status(r_ap.net)
    vmin = r_ap.final_converged ? minimum(bus_vm(r_ap.net)) : NaN
    push!(margin_table, (lambda=lambda, dmin=st.apslf_convergence_radius, level=st.apslf_convergence_level, apslf=r_ap.final_converged, nr=r_nr.final_converged, vmin=vmin))
end
println("   λ      dmin    level   APSLF   NR     min Vm")
for row in margin_table
    @printf("  %.1f   %7.3f   %-5s   %-5s   %-5s  %.3f\n", row.lambda, row.dmin, row.level, row.apslf, row.nr, row.vmin)
end
@assert margin_table[1].apslf && margin_table[1].nr                         #src
@assert margin_table[1].dmin > margin_table[5].dmin                         #src
@assert !margin_table[end].apslf                                            #src

# Reading aid: as $\lambda$ grows, the lowest voltage falls, `dmin`
# shrinks and the level moves from `GRN` towards `RED`, all on a network
# that still solves. Past the fold neither solver converges; the
# difference is that APSLF told you it was coming. In a run log or on the
# Web UI runs page this is the line `APSLF radius`, reported next to the
# Jacobian condition number of a Newton run.
#
# The row before the fold shows the other use of the margin. At that
# load the default order no longer reaches the tolerance while Newton
# still converges: the pole has come close enough to $s = 1$ that 24
# coefficients do not settle the Padé approximant. The `YEL` level with
# a non-converged run is the signal to raise the order before anything
# else is changed:

lambda_yel = margin_table[end-1].lambda
cfg_order40 = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=40)), output=quiet)
r_yel = run_sparlectra(net=scaled_ring7(lambda_yel), config=cfg_order40)
@printf("λ = %.1f with order 40: %s, mismatch %.1e, %s\n", lambda_yel, r_yel.outcome, r_yel.final_mismatch, Sparlectra.rectangular_pf_status(r_yel.net).apslf_convergence_line)
@assert r_yel.final_converged                                               #src

# The radius evaluation costs about as much as the solve itself: for
# every non-slack bus it forms the Padé denominator and finds its roots,
# one small eigenvalue problem per bus. On a large network it can be
# switched off with `power_flow.apslf.convergence_radius: false`, and
# the status line then says so instead of showing a number.
#
# Two reading rules for `dmin`. It is a heuristic continuation margin,
# not a stability certificate: Padé approximants can produce spurious
# poles, so a `dmin` you want to rely on should be checked against a
# second order (Part 1's table shows how cheap that is). And a
# non-converged run at `YEL` is a truncation signal, not a loadability
# verdict; raise the order first, as shown above.

cfg_noradius = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(convergence_radius=false)), output=quiet)
r_noradius = run_sparlectra(net=build_ring7("ring7 no radius"), config=cfg_noradius)
println(Sparlectra.rectangular_pf_status(r_noradius.net).apslf_convergence_line)
@assert !isfinite(Sparlectra.rectangular_pf_status(r_noradius.net).apslf_convergence_radius)   #src

# ### Reactive limits under APSLF
#
# A generator that would leave its reactive band cannot stay a PV bus.
# The rectangular solver switches such a machine to PQ at the binding
# limit inside its active-set loop; APSLF does the same in its own outer
# passes (the `iterations` count of Part 1). Sparlectra registers the
# clamps the solver reports in the net's Q-limit log, so the result table
# and the Q-V check judge an APSLF run exactly like a Newton run. The
# shipped `sp_case118` (a synthetic case with the cardinalities of the
# IEEE 118-bus system: 118 buses, 186 branches, 54 generators, of which
# a good third are synchronous condensers with tight reactive bands) has
# several machines at their limits in the base case:

case118 = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "mpower", "sp_case118.m")
r118_nr = run_sparlectra(casefile=basename(case118), path=dirname(case118), config=cfg_nr)
r118_ap = run_sparlectra(casefile=basename(case118), path=dirname(case118), config=cfg_apslf)
clamped(net) = sort!(collect(keys(net.qLimitEvents)))
println("NR    : ", r118_nr.outcome, ", ", length(clamped(r118_nr.net)), " machine(s) at a Q-limit, ", r118_nr.iterations, " iterations")
println("APSLF : ", r118_ap.outcome, ", ", length(clamped(r118_ap.net)), " machine(s) at a Q-limit, ", r118_ap.iterations, " pass(es)")
println("APSLF : ", Sparlectra.rectangular_pf_status(r118_ap.net).apslf_convergence_line)
@assert r118_nr.final_converged && r118_ap.final_converged                  #src
@assert !isempty(clamped(r118_ap.net))                                      #src

# The two clamp sets need not be identical, and neither the voltages:
# which machines end at a limit depends on the path taken, and a case
# with many machines near their limits can have more than one valid
# limit configuration (the Q-limit strategy
# question of the [control chapter](https://welthulk.github.io/Sparlectra.jl/generated/workshop_tour_control/); the rectangular
# solver's active-set loop and the series solver's outer passes are two
# such paths). This is the one place where "same equations, same answer"
# does not hold, and the run log's Q-limit block is where to look when
# two solvers disagree on a case with many machines at their limits:

common = intersect(clamped(r118_nr.net), clamped(r118_ap.net))
println(length(common), " machine(s) clamped by both solvers, ", length(symdiff(clamped(r118_nr.net), clamped(r118_ap.net))), " differ")
@printf("max |ΔVm| over all buses: %.2e pu\n", maximum(abs.(bus_vm(r118_nr.net) .- bus_vm(r118_ap.net))))

# ## Part 3: the series as a start value for Newton
#
# Where a Newton solve is wanted after all, the series solution is the
# start value with the shortest way to go. `power_flow.apslf_start.enabled`
# keeps the rectangular solver as the executing solver and runs one APSLF
# solve (order `apslf_start.order`, default 40) ahead of it as the start
# value generator. On the stressed ring from Part 2 that shortens the
# Newton run from the flat start visibly:

lambda_hard = 3.5
cfg_hybrid = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, apslf_start=Sparlectra.ApslfStartConfig(enabled=true)), output=quiet)
r_flat = run_sparlectra(net=scaled_ring7(lambda_hard), config=cfg_nr)
r_hyb = run_sparlectra(net=scaled_ring7(lambda_hard), config=cfg_hybrid)
println("flat start  : ", r_flat.outcome, ", ", r_flat.iterations, " Newton iterations")
println("APSLF start : ", r_hyb.outcome, ", ", r_hyb.iterations, " Newton iterations")
@assert r_hyb.final_converged                                               #src
@assert r_hyb.iterations <= r_flat.iterations                               #src

# The hybrid start is also what the automatic rescue ladder of a
# contingency batch reaches for when a warm start fails
# (`contingency.rescue_ladder`), and what the Web UI offers as
# "APSLF start" next to the solver choice.
#
# ### Where APSLF is not the right tool
#
# The series solves the algebraic power-flow equations and nothing else.
# Outer-loop controllers (tap changers with a voltage target, Q(U)
# characteristics, remote voltage control) change the model between
# solves, and Sparlectra refuses the combination rather than running the
# solver on a model whose controllers would stay silent. The shipped
# `sp_case14` carries such a tap controller, so the run below is expected
# to be refused; the printed line is the refusal, not a defect:

case14 = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case14.scf.json")
rejected = try
    run_sparlectra(net=importSCF(case14), config=cfg_apslf)
    ""
catch err
    sprint(showerror, err)
end
println("refused as intended: ", first(rejected, 120), " ...")
@assert !isempty(rejected)                                                  #src

# For such a network the hybrid start of Part 3 is the way to use the
# series: the controllers run in the rectangular outer loop, the series
# only supplies the start value. The same network solves that way, tap
# controller included:

r14 = run_sparlectra(net=importSCF(case14), config=cfg_hybrid)
println("hybrid start on sp_case14: ", r14.outcome, ", ", r14.iterations, " Newton iterations")
@assert r14.final_converged                                                 #src

# ## Summary
#
# - `power_flow.solver: apslf` solves a Sparlectra network by the
#   analytic power series; the result carries the same status contract as
#   a Newton run (outcome, mismatch, Q-limit log, result table).
# - `power_flow.apslf.order` is the one numerical knob; the mismatch
#   plateau shows when it is high enough.
# - The convergence radius `dmin` is a loadability margin that comes with
#   every solve; read it before raising the order or blaming the data.
# - `power_flow.apslf_start.enabled` uses the series as the start value
#   of the rectangular solver, which is the way to combine it with
#   controllers.
#
# The equivalent YAML for a configuration file:
#
# ```yaml
# power_flow:
#   solver: apslf          # rectangular | apslf | dc
#   apslf:
#     order: 24
#     convergence_radius: true
#   apslf_start:
#     enabled: false       # true: series start for the rectangular solver
#     order: 40
# ```
#
# See [Power Flow Configuration](https://welthulk.github.io/Sparlectra.jl/powerflow_configuration/) for every key and the
# [AnalyticLoadFlow.jl workshop](https://github.com/Welthulk/AnalyticLoadFlow.jl) for the theory, the recursion by hand, and the raw solver API.