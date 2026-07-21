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

using Sparlectra
using AnalyticLoadFlow

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
