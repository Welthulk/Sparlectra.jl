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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# DC power-flow reporting and the public rundcpf! entry point.

# file: src/powerflow_dc/dc_report.jl

"""
    DcPowerFlowReport

Structured container for DC power-flow results (see [`rundcpf!`](@ref)).

Deliberately not [`ACPFlowReport`](@ref): the DC model has no voltage
magnitude solution (always 1.0 pu by definition), no reactive power, no
losses (exactly zero, lossless model), and no iteration/tolerance concept
(a direct linear solve, not a Newton iteration) — so those columns are
simply absent rather than filled with placeholder values.

# Fields
- `metadata`: Run/case metadata (`case`, `baseMVA`, `converged`, `elapsed_s`,
  `iterations`, `solver`, `timestamp`, plus `seed_ac_start`/`ac_converged`/
  `ac_iterations`/`ac_elapsed_s` when `rundcpf!` was called with
  `seed_ac_start=true`).
- `nodes`: Per-bus angle and specified generation/load.
- `branches`: Per-branch directional MW flow (`p_to_MW == -p_from_MW`, lossless).
"""
struct DcPowerFlowReport
  metadata::NamedTuple
  nodes::Vector{NamedTuple}
  branches::Vector{NamedTuple}
end

function Base.show(io::IO, report::DcPowerFlowReport)
  print(io, "DcPowerFlowReport(", "case=", report.metadata.case, ", converged=", report.metadata.converged, ", nodes=", length(report.nodes), ", branches=", length(report.branches), ")")
end

"""
    buildDcPowerFlowReport(net::Net; ct=0.0, converged=true) -> DcPowerFlowReport

Builds a structured report from a `net` already solved by [`rundcpf!`](@ref)
(or `_run_dc_powerflow!`). Reuses the same bus-name/power-component helpers
as [`buildACPFlowReport`](@ref) so DC and AC reports stay visually
consistent, but only includes the columns the DC model actually defines.
"""
function buildDcPowerFlowReport(net::Net; ct::Float64 = 0.0, converged::Bool = true)::DcPowerFlowReport
  busNameByIdx = _bus_name_by_idx(net)
  nodes_sorted = sort(net.nodeVec, by = x -> x.busIdx)
  power_components = _bus_power_component_cache(net)

  node_rows = NamedTuple[]
  for n in nodes_sorted
    p_gen, _, p_load, _ = _bus_power_components(power_components, n.busIdx)
    push!(node_rows, (bus = n.busIdx, bus_name = _effective_bus_name(busNameByIdx, net, n.busIdx), type = toString(n._nodeType), va_deg = n._va_deg, p_gen_MW = p_gen, p_load_MW = p_load, original_bus_name = _original_bus_name(busNameByIdx, net, n.busIdx)))
  end

  branch_rows = NamedTuple[]
  for br in net.branchVec
    f = br.fBranchFlow
    t = br.tBranchFlow
    p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
    p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
    rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
    overload = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated
    from_name = _effective_bus_name(busNameByIdx, net, br.fromBus)
    to_name = _effective_bus_name(busNameByIdx, net, br.toBus)
    push!(branch_rows, (branch = br.comp.cName, branch_index = br.branchIdx, from_bus = br.fromBus, to_bus = br.toBus, status = br.status, p_from_MW = p_from, p_to_MW = p_to, rated_MVA = rated, overloaded = overload, from_bus_name = from_name, to_bus_name = to_name, branch_kind = _branch_kind_name(net, br)))
  end

  metadata = (case = net.name, baseMVA = net.baseMVA, converged = converged, elapsed_s = ct, iterations = 1, solver = :dc, timestamp = Dates.now())
  return DcPowerFlowReport(metadata, node_rows, branch_rows)
end

"""
    printDcPowerFlowResults(net::Net, ct::Float64; toFile=false, path="")

Compact text summary of a solved DC power flow, printed to `stdout` (or a
`result_<net.name>.txt` file when `toFile=true`, mirroring
`printACPFlowResults`'s file-output convention). Omits the Q/tap-
controller/wrong-branch sections that don't apply to the DC model.
"""
function printDcPowerFlowResults(net::Net, ct::Float64; toFile::Bool = false, path::String = "")
  if toFile
    filename = strip("result_$(net.name).txt")
    io = open(joinpath(path, filename), "w")
    @info "Results are written to $(joinpath(path, filename))"
  else
    io = Base.stdout
  end

  report = buildDcPowerFlowReport(net; ct = ct, converged = true)
  println(io, "=== DC power flow results for \"$(net.name)\" ===")
  @printf(io, "solver = dc, elapsed = %.4f s\n\n", ct)
  println(io, "Bus results:")
  @printf(io, "%6s  %-20s  %10s  %10s  %10s\n", "Bus", "Name", "Va[deg]", "Pgen[MW]", "Pload[MW]")
  for row in report.nodes
    @printf(io, "%6d  %-20s  %10.4f  %10.4f  %10.4f\n", row.bus, row.bus_name, _default0(row.va_deg), row.p_gen_MW, row.p_load_MW)
  end
  println(io, "\nBranch flows:")
  @printf(io, "%6s  %-12s  %6s  %6s  %10s  %10s\n", "Branch", "Name", "From", "To", "Pfrom[MW]", "Pto[MW]")
  for row in report.branches
    @printf(io, "%6d  %-12s  %6d  %6d  %10.4f  %10.4f\n", row.branch_index, row.branch, row.from_bus, row.to_bus, row.p_from_MW, row.p_to_MW)
  end

  toFile && close(io)
  return nothing
end

"""
    _dc_seed_rectangular_angles!(net, pf_cfg; verbose=0, performance_profile=nothing)

Run the standalone DC power flow and write its bus angles into `net` as a
Newton-Raphson start seed, without disturbing Slack/PV regulated voltage
magnitudes. The DC model has no voltage-magnitude solution — it writes
`vm_pu=1.0` for every bus — which is correct for a standalone DC report but
would erase Slack/PV setpoints ahead of an AC seed (`initialVrect`'s
non-flatstart path reads `node._vm_pu` directly, with no separate setpoint
fallback). Snapshots those setpoints before the DC solve and restores them
immediately after, so only the angles end up DC-seeded. Shared by
[`rundcpf!`](@ref)'s `seed_ac_start=true` path and the config-driven
`power_flow.start_mode.dc_seed_unconditional` path in `execution.jl`.
"""
function _dc_seed_rectangular_angles!(net::Net, pf_cfg::PowerFlowConfig; verbose::Int = 0, performance_profile = nothing)
  vm_regulated = [getNodeType(n) in (Slack, PV) ? n._vm_pu : nothing for n in net.nodeVec]
  _run_dc_powerflow!(net, pf_cfg; verbose, performance_profile)
  @inbounds for (k, n) in enumerate(net.nodeVec)
    vm = vm_regulated[k]
    (vm !== nothing && vm > 0.0) && (n._vm_pu = vm)
  end
  return net
end

"""
    rundcpf!(net::Net; angle_reference_deg=0.0, verbose=0, performance_profile=nothing,
             seed_ac_start=false, ac_kwargs=NamedTuple()) -> DcPowerFlowReport

Solve the standalone, MATPOWER-`rundcpf`-equivalent DC power flow and write
the bus angles (voltage magnitude implicitly 1.0 pu everywhere, by DC-model
definition) and lossless branch flows back into `net`. The run is flagged
via [`dc_pf_status`](@ref) — a status registry kept separate from
`rectangular_pf_status`, so AC reporting/status code never mistakes a DC
solution for an AC one. Respects `power_flow.islands` (each AC island is
solved independently, `matpower_like` reference bus). Ported from
MATPOWER's `rundcpf`/`makeBdc` (BSD-3-Clause; algorithm only, no code
copied).

# Keyword arguments
- `angle_reference_deg`: uniform angle offset added to every bus after the
  slack-referenced solve (see [`solve_dc_powerflow`](@ref) for why this is
  exact, not approximate). The slack bus itself is fixed at this reference.
- `seed_ac_start`: if `true`, after the DC solve succeeds, immediately runs
  the AC rectangular Newton-Raphson solve ([`runpf_rectangular!`](@ref))
  seeded from the just-computed DC angles (`opt_flatstart=false`, so the
  solver picks up `net`'s current `_va_deg`/`_vm_pu` as its starting point —
  the DC step already wrote these). Slack/PV voltage-magnitude setpoints
  (which the DC write would otherwise flatten to 1.0 pu along with every
  other bus) are snapshotted before the DC solve and restored onto `net`
  right before the AC call, so only the angles are DC-seeded. The AC call
  also defaults to `damp=1.0` (full Newton steps, matching what
  `run_sparlectra`'s config-driven path uses — `runpf_rectangular!`'s own
  bare keyword default, `damp=0.2`, is deliberately conservative for
  unmediated direct calls and would waste a good start on needlessly small
  steps). On return, `net` holds the **AC-converged** solution (not the DC
  one) when this option is used; the returned `DcPowerFlowReport` still
  reflects the DC step's own result (a useful before/after comparison
  point), with the AC outcome recorded in
  `report.metadata.ac_converged`/`ac_iterations`/`ac_elapsed_s`.
- `ac_kwargs`: extra keyword arguments forwarded to `runpf_rectangular!`
  when `seed_ac_start=true` (e.g. `(tol = 1e-9, maxiter = 40)`).
  `opt_flatstart` is always forced to `false` and cannot be overridden this
  way, since that is the entire point of seeding; `damp` defaults to `1.0`
  as described above but can be overridden here.
- `verbose`, `performance_profile`: forwarded to both the DC solve and (when
  `seed_ac_start=true`) the AC solve.

Does not support outer-loop controllers (tap/PST/shunt/... control):
the *static* effect of a phase-shifting transformer's current
`phase_shift_deg` is correctly represented in the DC model (it enters the
B′/injection math), but no outer loop adjusts it — mirrors the existing
`power_flow.solver=:apslf` restriction.
"""
function rundcpf!(net::Net; angle_reference_deg::Float64 = 0.0, verbose::Int = 0, performance_profile = nothing, seed_ac_start::Bool = false, ac_kwargs::NamedTuple = NamedTuple())::DcPowerFlowReport
  t0 = time()
  pf_cfg = PowerFlowConfig(dc = DcPowerFlowConfig(angle_reference_deg = angle_reference_deg))
  if seed_ac_start
    _dc_seed_rectangular_angles!(net, pf_cfg; verbose, performance_profile)
  else
    _run_dc_powerflow!(net, pf_cfg; verbose, performance_profile)
  end
  ct = time() - t0
  report = buildDcPowerFlowReport(net; ct = ct, converged = true)

  if seed_ac_start
    t1 = time()
    # runpf_rectangular!'s own bare keyword default (damp=0.2) is deliberately
    # conservative for unmediated direct calls; the config-driven run_sparlectra
    # path always uses damp=1.0 (full Newton steps, PowerFlowConfig has no damp
    # field of its own — the outer active-set/framework loop supplies it as a
    # runtime kwarg). Match that here so a good DC-angle start isn't wasted on
    # needlessly damped steps; still overridable via ac_kwargs.
    ac_damp = get(ac_kwargs, :damp, 1.0)
    ac_kwargs_clean = Base.structdiff(ac_kwargs, NamedTuple{(:opt_flatstart, :damp)})
    ac_iters, ac_erg = Base.invokelatest(runpf_rectangular!, net; opt_flatstart = false, damp = ac_damp, verbose = verbose, performance_profile = performance_profile, ac_kwargs_clean...)
    ac_ct = time() - t1
    verbose > 0 && println("DC-seeded AC power flow: ", ac_erg == 0 ? "converged" : "did NOT converge", " in ", ac_iters, " iteration(s), ", round(ac_ct; digits = 4), " s")
    metadata = merge(report.metadata, (seed_ac_start = true, ac_converged = ac_erg == 0, ac_iterations = ac_iters, ac_elapsed_s = ac_ct))
    report = DcPowerFlowReport(metadata, report.nodes, report.branches)
  end

  return report
end
