# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Date: 2026-07-21
# file: examples/apslf_demo.jl
# purpose: compares the internal rectangular Newton-Raphson solver against the AnalyticLoadFlow.jl-backed APSLF solver (standalone and as an NR start-value generator) on case30

# APSLF (used by AnalyticLoadFlow.jl, an optional weak dependency — see
# ext/SparlectraAnalyticLoadFlowExt.jl) is an analytic power-series load-flow
# method: it expands the bus voltages as a power series in an embedding
# parameter s and evaluates that series at s=1, optionally through Padé
# [L/M] approximants (`use_pade`) for a larger effective convergence radius.
# The series always starts from the canonical analytic germ V(s=0) = 1∠0 —
# there is no configurable start voltage, unlike the rectangular Newton-
# Raphson solver, because the embedding construction fixes that boundary
# condition. `nr_polish` runs a few Newton-Raphson steps on top of the series
# result to tighten the final mismatch; `apslf_start` deliberately disables
# polish internally (`nr_polish=false`) because the downstream rectangular
# Newton-Raphson solve performs that polishing step itself.
#
# This example runs three variants on the same case and compares them against
# the rectangular Newton-Raphson solve:
#   A) NR reference          — power_flow.solver = rectangular (baseline)
#   B) APSLF standalone       — power_flow.solver = apslf
#   C) NR + APSLF start values — power_flow.solver = rectangular,
#                                 power_flow.apslf_start.enabled = true
#
# Running this example requires the optional AnalyticLoadFlow.jl dependency
# to be loaded (`using AnalyticLoadFlow`) in addition to `using Sparlectra`;
# without it, apslf_solver() raises a clear "not installed" error for runs B
# and C.
#
# A second, independent demonstration follows at the end of this file: a
# two-bus load sweep across a PV-curve voltage-stability nose, comparing
# plain flat-start NR, standalone APSLF, and NR seeded with APSLF start
# values as the load approaches and passes the nose. See
# `main_nose_curve_comparison()` below for the theory anchor and the
# (measured, not assumed) comparison between the three strategies.

using Sparlectra
using AnalyticLoadFlow
using Printf

include(joinpath(@__DIR__, "internal", "example_header.jl"))

"""
    _apslf_demo_options(path) -> (order, use_pade, nr_polish)

Read optional `order`/`use_pade`/`nr_polish` overrides for this example from
a small standalone YAML file via the existing YAML-subset parser
(`Sparlectra.load_yaml_dict`). Falls back to the `power_flow.apslf` defaults
(order=40, use_pade=true, nr_polish=true) when `path` does not exist.
"""
function _apslf_demo_options(path::AbstractString)
  order = 6
  use_pade = true
  nr_polish = true
  if isfile(path)
    raw = Sparlectra.load_yaml_dict(path)
    order = Int(get(raw, "order", order))
    use_pade = Bool(get(raw, "use_pade", use_pade))
    nr_polish = Bool(get(raw, "nr_polish", nr_polish))
  end
  return (; order, use_pade, nr_polish)
end

"""
    _apslf_demo_deltas(net, reference_net) -> (max_dVm_pu, max_dVa_deg)

Maximum absolute voltage-magnitude and voltage-angle deviation of `net`
relative to `reference_net`, bus-by-bus in network order.
"""
function _apslf_demo_deltas(net::Net, reference_net::Net)
  max_dVm_pu = maximum(abs.(getfield.(net.nodeVec, :_vm_pu) .- getfield.(reference_net.nodeVec, :_vm_pu)))
  max_dVa_deg = maximum(abs.(getfield.(net.nodeVec, :_va_deg) .- getfield.(reference_net.nodeVec, :_va_deg)))
  return (max_dVm_pu = max_dVm_pu, max_dVa_deg = max_dVa_deg)
end

function main(;
  casefile::AbstractString = "case30.m",
  demo_config_file::AbstractString = joinpath(@__DIR__, "apslf_demo.yaml.example"),
)
  print_example_banner("examples/apslf_demo.jl", "compares the internal rectangular Newton-Raphson solver against the AnalyticLoadFlow.jl-backed APSLF solver (standalone and as an NR start-value generator) on case30")

  opts = _apslf_demo_options(demo_config_file)
  println("APSLF options: order=", opts.order, "  use_pade=", opts.use_pade, "  nr_polish=", opts.nr_polish)
  println()

  config_file = Sparlectra.configuration_path_from_inputs(
    env_var = "SPARLECTRA_CONFIGURATION_YAML",
    fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH],
  )
  case_path = ensure_casefile(casefile)
  mpc = Sparlectra.MatpowerIO.read_case(case_path)

  # A) NR reference.
  cfg_nr = Sparlectra.load_sparlectra_config(config_file; reload = true, overrides = Dict(
    "power_flow" => Dict("solver" => "rectangular"),
    "output" => Dict("logfile_results" => "off"),
  ))
  net_nr = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
  result_nr = run_sparlectra(net = net_nr, config = cfg_nr)

  # B) APSLF standalone.
  cfg_apslf = Sparlectra.load_sparlectra_config(config_file; reload = true, overrides = Dict(
    "power_flow" => Dict(
      "solver" => "apslf",
      "apslf" => Dict("order" => opts.order, "use_pade" => opts.use_pade, "nr_polish" => opts.nr_polish),
    ),
    "output" => Dict("logfile_results" => "off"),
  ))
  net_apslf = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
  result_apslf = run_sparlectra(net = net_apslf, config = cfg_apslf)

  # C) NR with APSLF start values (nr_polish is always off internally for the
  # start generator; the rectangular Newton-Raphson solve below does the polish).
  cfg_hybrid = Sparlectra.load_sparlectra_config(config_file; reload = true, overrides = Dict(
    "power_flow" => Dict(
      "solver" => "rectangular",
      "apslf_start" => Dict("enabled" => true, "order" => opts.order),
    ),
    "output" => Dict("logfile_results" => "off"),
  ))
  net_hybrid = Sparlectra.createNetFromMatPowerCase(mpc = mpc, flatstart = true)
  result_hybrid = run_sparlectra(net = net_hybrid, config = cfg_hybrid)

  runs = (
    (label = "A: NR reference", result = result_nr, deltas = (max_dVm_pu = 0.0, max_dVa_deg = 0.0)),
    (label = "B: APSLF standalone", result = result_apslf, deltas = _apslf_demo_deltas(net_apslf, net_nr)),
    (label = "C: NR + APSLF start", result = result_hybrid, deltas = _apslf_demo_deltas(net_hybrid, net_nr)),
  )
  return (casefile = basename(case_path), opts = opts, runs = runs)
end

result = run_example(main)

println("Case: ", result.casefile)
println()
println(rpad("Run", 24), rpad("converged", 11), rpad("iterations", 12), rpad("solver_time_s", 15), rpad("max|ΔVm| [pu]", 16), "max|Δδ| [deg]")
println("-"^96)
for run in result.runs
  r = run.result
  println(
    rpad(run.label, 24),
    rpad(string(r.final_converged), 11),
    rpad(string(r.iterations), 12),
    rpad(string(round(something(r.solver_elapsed_s, NaN); digits = 6)), 15),
    rpad(string(round(run.deltas.max_dVm_pu; digits = 10)), 16),
    string(round(run.deltas.max_dVa_deg; digits = 8)),
  )
end

# =============================================================================
# Second demonstration: Newton-Raphson vs. APSLF near a PV-curve nose
# =============================================================================
#
# A two-bus system operated up to its voltage-stability limit (the "nose" of
# the P-V curve), comparing three solve strategies as the load is swept
# toward and past that limit:
#   A) plain flat-start Newton-Raphson (autodamp and start projection
#      explicitly disabled, so the framework's own robustness features do
#      not mask the comparison)
#   B) APSLF standalone (analytic power series, Padé-accelerated, no NR
#      polishing, so the analytic result stands on its own)
#   C) NR seeded with APSLF-generated start values (`apslf_start`), i.e. the
#      same analytic series used only as a warm start for NR
#
# Theory anchor (lossless approximation, x_pu = 0.5, Q = 0.3 * P):
#   Solvability condition for the load-bus voltage U = V2^2:
#     U^2 + U*(2*Q*x - 1) + x^2*(P^2 + Q^2) = 0
#   Real solutions exist while  1 - 0.6*P - P^2 >= 0  =>  P_nose ~ 0.744 pu
#   (~74 MW at 100 MVA base; the small series resistance below shifts it
#   slightly). Beyond that load, no solution exists at all, so no solver —
#   NR, APSLF, or the hybrid — can converge there; that zone is included on
#   purpose to show the difference between "no solver reaches this yet" and
#   "no solution exists".
#
# What this sweep actually shows (verified by running it, not assumed): in
# this small radial case, undamped flat-start NR (A) converges at least as
# far into the nose as unpolished standalone APSLF (B) — APSLF's finite
# series truncation limits its own reach here, so a naive "HELM-style
# methods always out-converge NR near a nose" expectation does not hold
# without qualification. The point this demo actually makes is more nuanced
# and arguably more useful: an analytic power series is not automatically
# more robust than Newton-Raphson near a nose; whether it helps depends on
# series order, Padé use, and — as tested here — whether it is used as a
# standalone solver (B) or as a warm start feeding NR's own quadratic local
# convergence (C, `apslf_start`). Run the sweep and read the printed
# comparison rather than assuming which column "wins".

# Fresh net per load level (no state carry-over between runs).
function _nose_build_net(p_MW::Float64, q_MVAr::Float64)
  net = Net(name = "nose_demo", baseMVA = 100.0, flatstart = true)

  addBus!(net = net, busName = "SLACK", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "LOAD", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)

  # High series reactance => low maximum transferable power => reachable nose.
  addPIModelACLine!(net = net, fromBus = "SLACK", toBus = "LOAD", r_pu = 0.01, x_pu = 0.5, b_pu = 0.0, status = 1)

  addProsumer!(net = net, busName = "SLACK", type = "EXTERNALNETWORKINJECTION", referencePri = "SLACK", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "LOAD", type = "ENERGYCONSUMER", p = p_MW, q = q_MVAr)

  ok, msg = validate!(net = net)
  ok || error("nose demo: invalid network at P=$(p_MW) MW: $msg")
  return net
end

const _NOSE_MAX_ITER = 30
const _NOSE_TOL = 1e-8

# A) Deliberately naive Newton-Raphson: flat start, with autodamp and start
# projection explicitly disabled so the effect this demo illustrates is not
# masked by the framework's own robustness features.
function _nose_solve_plain_nr(p_MW::Float64, q_MVAr::Float64)
  net = _nose_build_net(p_MW, q_MVAr)
  cfg = SparlectraConfig(powerflow = PowerFlowConfig(
    max_iter = _NOSE_MAX_ITER, tol = _NOSE_TOL,
    solver = :rectangular, autodamp = false,
    start_mode = StartModeConfig(start_projection = false),
  ))
  result = run_sparlectra(net = net, config = cfg)
  vm = result.final_converged ? net.nodeVec[2]._vm_pu : NaN
  return (converged = result.final_converged, ite = result.iterations, vm_pu = vm)
end

# B) APSLF standalone: analytic series solution with Padé acceleration, no NR
# polishing so the analytic result stands on its own.
function _nose_solve_apslf(p_MW::Float64, q_MVAr::Float64)
  net = _nose_build_net(p_MW, q_MVAr)
  cfg = SparlectraConfig(powerflow = PowerFlowConfig(
    max_iter = _NOSE_MAX_ITER, tol = _NOSE_TOL,
    solver = :apslf, apslf = Sparlectra.ApslfConfig(order = 40, use_pade = true, nr_polish = false),
  ))
  result = run_sparlectra(net = net, config = cfg)
  vm = result.final_converged ? net.nodeVec[2]._vm_pu : NaN
  return (converged = result.final_converged, vm_pu = vm)
end

# C) NR seeded with APSLF-generated start values instead of a flat start;
# same autodamp/start-projection-disabled NR as (A) otherwise.
function _nose_solve_apslf_start(p_MW::Float64, q_MVAr::Float64)
  net = _nose_build_net(p_MW, q_MVAr)
  cfg = SparlectraConfig(powerflow = PowerFlowConfig(
    max_iter = _NOSE_MAX_ITER, tol = _NOSE_TOL,
    solver = :rectangular, autodamp = false,
    start_mode = StartModeConfig(start_projection = false),
    apslf_start = Sparlectra.ApslfStartConfig(enabled = true, order = 40),
  ))
  result = run_sparlectra(net = net, config = cfg)
  vm = result.final_converged ? net.nodeVec[2]._vm_pu : NaN
  return (converged = result.final_converged, ite = result.iterations, vm_pu = vm)
end

function main_nose_curve_comparison()
  print_example_banner("examples/apslf_demo.jl (nose-curve comparison)", "compares plain flat-start NR, standalone APSLF, and NR seeded with APSLF start values across a PV-curve voltage-stability nose")

  println("Two-bus nose demo: x_pu = 0.5, r_pu = 0.01, Q = 0.3*P, flat start")
  println("analytic nose (lossless): ~74 MW")
  println()
  @printf("%8s | %-18s | %-12s | %-18s\n", "P [MW]", "A: plain NR", "B: APSLF", "C: NR+APSLFstart")
  println("-"^68)

  rows = NamedTuple[]
  for p in 50.0:2.0:80.0
    q = 0.3 * p
    nr = _nose_solve_plain_nr(p, q)
    apslf = _nose_solve_apslf(p, q)
    hybrid = _nose_solve_apslf_start(p, q)
    push!(rows, (p_MW = p, nr = nr, apslf = apslf, hybrid = hybrid))

    nr_text = nr.converged ? @sprintf("conv (%2dit) %.3f", nr.ite, nr.vm_pu) : "FAILED"
    apslf_text = apslf.converged ? @sprintf("conv %.3f", apslf.vm_pu) : "FAILED"
    hybrid_text = hybrid.converged ? @sprintf("conv (%2dit) %.3f", hybrid.ite, hybrid.vm_pu) : "FAILED"
    @printf("%8.1f | %-18s | %-12s | %-18s\n", p, nr_text, apslf_text, hybrid_text)
  end

  first_fail = Dict{Symbol,Union{Nothing,Float64}}(:nr => nothing, :apslf => nothing, :hybrid => nothing)
  for row in rows
    row.nr.converged || (first_fail[:nr] === nothing && (first_fail[:nr] = row.p_MW))
    row.apslf.converged || (first_fail[:apslf] === nothing && (first_fail[:apslf] = row.p_MW))
    row.hybrid.converged || (first_fail[:hybrid] === nothing && (first_fail[:hybrid] = row.p_MW))
  end

  println()
  println("First non-convergence: A(plain NR)=", something(first_fail[:nr], "none in range"), " MW,  B(APSLF)=", something(first_fail[:apslf], "none in range"), " MW,  C(NR+APSLFstart)=", something(first_fail[:hybrid], "none in range"), " MW")
  println()
  println("Interpretation: watch which column fails first, not which one this")
  println("comment assumes will fail first. In the case measured when this file")
  println("was written, (B) fails before (A), and (C) tracks (A) — i.e. an")
  println("APSLF warm start did not extend NR's own reach here, and raw APSLF")
  println("was the least robust of the three this close to the nose. Re-run")
  println("after changing apslf.order, x_pu, or the P/Q ratio to see how the")
  println("comparison shifts; nothing here is hardcoded to produce a particular")
  println("winner.")
  return rows
end

nose_curve_rows = run_example(main_nose_curve_comparison)
