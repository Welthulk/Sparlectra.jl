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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 07.09.2023
# file: src/results.jl
# Purpose: functions for formatting and printing results of power flow calculations
function format_version(version::VersionNumber)
  major = lpad(version.major, 2, '0')
  minor = lpad(version.minor, 1, '0')
  patch = lpad(version.patch, 2, '0')
  return "$major.$minor.$patch"
end

"""
    ACPFlowReport

Structured container for AC power flow results.

The vectors (`nodes`, `branches`, `links`, `q_limit_events`) are table-like and can
be converted directly to `DataFrame`s if `DataFrames.jl` is available, e.g.
`DataFrame(report.nodes)`.
"""
struct ACPFlowReport
  metadata::NamedTuple
  nodes::Vector{NamedTuple}
  branches::Vector{NamedTuple}
  links::Vector{NamedTuple}
  q_limit_events::Vector{NamedTuple}
end

function Base.show(io::IO, report::ACPFlowReport)
  print(io,
    "ACPFlowReport(",
    "case=", report.metadata.case,
    ", converged=", report.metadata.converged,
    ", nodes=", length(report.nodes),
    ", branches=", length(report.branches),
    ", links=", length(report.links),
    ", q_limit_events=", length(report.q_limit_events),
    ")")
end

_default0(x) = isnothing(x) ? 0.0 : x

function _bus_name_by_idx(net::Net)
  busNameByIdx = Dict{Int,String}()
  for (name, idx) in net.busDict
    busNameByIdx[idx] = name
  end
  return busNameByIdx
end

function _effective_pf_node_count(net::Net)::Int
  reps = _active_link_representative_map(net)
  return length(unique(reps))
end

function _bus_control_flags(net::Net, bus_idx::Int)
  has_qu = false
  has_pu = false
  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == bus_idx || continue
    has_qu |= has_qu_controller(ps)
    has_pu |= has_pu_controller(ps)
  end
  return has_qu, has_pu
end

function _control_label(net::Net, bus_idx::Int)::String
  has_qu, has_pu = _bus_control_flags(net, bus_idx)
  if has_qu && has_pu
    return "Q(U), P(U)"
  elseif has_qu
    return "Q(U)"
  elseif has_pu
    return "P(U)"
  else
    return "-"
  end
end

function _effective_bus_power_components(net::Net, bus_idx::Int)
  base = net.baseMVA
  vm = net.nodeVec[bus_idx]._vm_pu
  vm_safe = vm > 1e-9 ? vm : 1e-9

  p_gen = 0.0
  q_gen = 0.0
  p_load = 0.0
  q_load = 0.0

  for ps in net.prosumpsVec
    getPosumerBusIndex(ps) == bus_idx || continue

    p_mw = isnothing(ps.pVal) ? 0.0 : ps.pVal
    q_mvar = isnothing(ps.qVal) ? 0.0 : ps.qVal

    if has_pu_controller(ps)
      p_pu, _ = evaluate_controller(ps.puController, vm_safe)
      p_mw = p_pu * base
    end
    if has_qu_controller(ps)
      q_pu, _ = evaluate_controller(ps.quController, vm_safe)
      q_mvar = q_pu * base
    end

    if isGenerator(ps)
      p_gen += p_mw
      q_gen += q_mvar
    else
      p_load += p_mw
      q_load += q_mvar
    end
  end

  return p_gen, q_gen, p_load, q_load
end

function _branch_kind_label(br::Branch)::String
  name = br.comp.cName
  if occursin("_ACL_", name)
    return "Line"
  elseif occursin("_2WT_", name)
    return "Trafo"
  elseif occursin("_PI_", name)
    return "PI"
  else
    return "Branch"
  end
end

"""
    buildACPFlowReport(net::Net; ...)

Builds a structured report object from solved network data.
This provides a machine-readable alternative to `printACPFlowResults`.
"""
function buildACPFlowReport(net::Net;
  ct::Float64 = 0.0,
  ite::Int = 0,
  tol::Float64 = 1e-6,
  converged::Bool = true,
  solver::Symbol = :NR,
)::ACPFlowReport
  busNameByIdx = _bus_name_by_idx(net)
  nodes_sorted = sort(net.nodeVec, by = x -> x.busIdx)

  node_rows = NamedTuple[]
  for n in nodes_sorted
    p_gen, q_gen, p_load, q_load = _effective_bus_power_components(net, n.busIdx)
    push!(node_rows, (
      bus = n.busIdx,
      bus_name = get(busNameByIdx, n.busIdx, n.comp.cName),
      type = toString(n._nodeType),
      vm_pu = n._vm_pu,
      va_deg = n._va_deg,
      vn_kV = n.comp.cVN,
      v_kV = n.comp.cVN * n._vm_pu,
      p_gen_MW = p_gen,
      q_gen_MVar = q_gen,
      p_load_MW = p_load,
      q_load_MVar = q_load,
      p_shunt_MW = _default0(n._pShunt),
      q_shunt_MVar = _default0(n._qShunt),
      is_isolated = isIsolated(n),
      q_limit_hit = haskey(net.qLimitEvents, n.busIdx),
      control = _control_label(net, n.busIdx),
    ))
  end

  branch_rows = NamedTuple[]
  for br in net.branchVec
    f = br.fBranchFlow
    t = br.tBranchFlow
    p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
    q_from = isnothing(f) || isnothing(f.qFlow) ? 0.0 : f.qFlow
    p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
    q_to = isnothing(t) || isnothing(t.qFlow) ? 0.0 : t.qFlow
    p_loss = _default0(br.pLosses)
    q_loss = _default0(br.qLosses)
    rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
    overload = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated

    push!(branch_rows, (
      branch = br.comp.cName,
      branch_index = br.branchIdx,
      from_bus = br.fromBus,
      to_bus = br.toBus,
      status = br.status,
      p_from_MW = p_from,
      q_from_MVar = q_from,
      p_to_MW = p_to,
      q_to_MVar = q_to,
      p_loss_MW = p_loss,
      q_loss_MVar = q_loss,
      rated_MVA = rated,
      overloaded = overload,
    ))
  end

  link_rows = NamedTuple[]
  for l in net.linkVec
    push!(link_rows, (
      link = l.cName,
      link_index = l.linkIdx,
      from_bus = l.fromBus,
      to_bus = l.toBus,
      status = l.status,
      p_MW = _default0(l.pFlow_MW),
      q_MVar = _default0(l.qFlow_MVar),
      ifrom_kA = _default0(l.iFrom_kA),
      ito_kA = _default0(l.iTo_kA),
    ))
  end

  q_events = NamedTuple[]
  for (bus, hit) in sort(collect(net.qLimitEvents), by = x -> x[1])
    push!(q_events, (bus = bus, hit = hit))
  end

  p_loss_total, q_loss_total = getTotalLosses(net = net)
  metadata = (
    case = net.name,
    baseMVA = net.baseMVA,
    converged = converged,
    elapsed_s = ct,
    iterations = ite,
    tolerance = tol,
    solver = solver,
    timestamp = Dates.now(),
    total_p_loss_MW = p_loss_total,
    total_q_loss_MVar = q_loss_total,
  )

  return ACPFlowReport(metadata, node_rows, branch_rows, link_rows, q_events)
end

function formatBranchResults(net::Net)
  busNameByIdx = _bus_name_by_idx(net)
  #! format: off
  formatted_results = @sprintf("\n========================================================================================================================================================\n")
  formatted_results *= @sprintf("| %-25s | %-6s | %-25s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n", "Branch", "Type", "Connection", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]")
  formatted_results *= @sprintf("========================================================================================================================================================\n")
  #! format: on
  for br in net.branchVec
    from = br.fromBus
    to = br.toBus
    bName = br.comp.cName
    branchKind = _branch_kind_label(br)
    fromName = get(busNameByIdx, Int(from), string(from))
    toName = get(busNameByIdx, Int(to), string(to))
    connection = string(fromName, " -> ", toName)
    if br.status == 0 || (isnothing(br.fBranchFlow)) || (isnothing(br.tBranchFlow))
      pfromVal = qfromVal = ptoVal = qtoVal = pLossval = qLossval = 0.0
    else
      pfromVal = (br.fBranchFlow.pFlow === nothing) ? NaN : br.fBranchFlow.pFlow
      qfromVal = (br.fBranchFlow.qFlow === nothing) ? NaN : br.fBranchFlow.qFlow

      ptoVal = (br.tBranchFlow.pFlow === nothing) ? NaN : br.tBranchFlow.pFlow
      qtoVal = (br.tBranchFlow.qFlow === nothing) ? NaN : br.tBranchFlow.qFlow

      pLossval = (br.pLosses === nothing) ? NaN : br.pLosses
      qLossval = (br.qLosses === nothing) ? NaN : br.qLosses
      ratedS = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA

      check = false
      if ratedS > 0.0
        if max(abs(pfromVal), abs(ptoVal)) > ratedS
          check = true
        end
      end

      if check
        bName *= " !"
      end
    end
    #! format: off
    formatted_results *= @sprintf("| %-25s | %-6s | %-25s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f |  %-10.3f|\n", bName, branchKind, connection, pfromVal, qfromVal, ptoVal, qtoVal, pLossval, qLossval)
    #! format: on
  end
  formatted_results *= @sprintf("--------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  (∑pv, ∑qv) = getTotalLosses(net = net)
  total_losses = @sprintf("total network power balance (Σ S_branch): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

  return formatted_results, total_losses
end

function printACPFlowResults(net::Net, ct::Float64, ite::Int, tol::Float64, toFile::Bool = false, path::String = ""; converged::Bool = true, solver::Symbol = :NR)
  if toFile
    filename = strip("result_$(net.name).txt")
    io = open(joinpath(path, filename), "w")
    @info "Results are written to $(joinpath(path, filename))"
  else
    io = Base.stdout
  end

  vers = Sparlectra.SparlectraVersion
  current_date = Dates.format(Dates.now(), "d-u-yy H:M:S")

  formatted_version = format_version(vers)
  flowResults, totalLosses = formatBranchResults(net)

  @printf(io, "================================================================================\n")
  @printf(io, "| SPARLECTRA Version %-10s - AC Power Flow Results                        |\n", formatted_version)
  @printf(io, "================================================================================\n")

  busses = length(net.nodeVec)
  pf_nodes = _effective_pf_node_count(net)
  branches = length(net.branchVec)
  links = length(net.linkVec)
  lines = length(net.linesAC)
  trafos = length(net.trafos)
  gens = 0
  loads = 0
  shunts = 0
  auxb = 0
  busNameByIdx = _bus_name_by_idx(net)

  nodes = sort(net.nodeVec, by = x -> x.busIdx)

  npv = 0
  npq = 0
  niso = 0
  for n in nodes
    npv += n._nodeType == Sparlectra.PV ? 1 : 0
    npq += isPQNode(n) ? 1 : 0
    niso += isIsolated(n) ? 1 : 0
    if occursin("_Aux_", n.comp.cName)
      auxb += 1
    end
  end
  for ps in net.prosumpsVec
    loads += ps.proSumptionType == Sparlectra.Consumption ? 1 : 0
    gens += ps.proSumptionType == Sparlectra.Injection ? 1 : 0
  end
  shunts = length(net.shuntVec)

  @printf(io, "Date           :%20s\n", current_date)
  @printf(io, "Iterations     :%10d\n", ite)
  @printf(io, "Flatstart      :%10s\n", net.flatstart ? "Yes" : "No")
  @printf(io, "Tolerance      : %.1e\n", tol)
  @printf(io, "Solver         :%15s\n", string(solver))
  if converged
    @printf(io, "Converged in   :%10f seconds\n", ct)
  else
    @printf(io, "Status         :%10s\n", "Not Converged")
  end
  @printf(io, "Case           :%15s\n", net.name)  
  @printf(io, "Cooldown iters :%10d\n", net.cooldown_iters)
  @printf(io, "Q-hysteresis   :%10.4f pu\n", net.q_hyst_pu)
  
  @printf(io, "BaseMVA        :%10d\n", net.baseMVA)
  if auxb > 0 && niso > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Iso: %d Slack: %d\n", busses, npv, npq, auxb, niso, 1)
  elseif auxb > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Slack: %d\n", busses, npv, npq, auxb, 1)
  else
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d Slack: %d)\n", busses, npv, npq, 1)
  end
  if pf_nodes != busses
    @printf(io, "PF Nodes       :%10d (after active-link merge)\n", pf_nodes)
  end
  @printf(io, "Branches       :%10d\n", branches)
  @printf(io, "Links          :%10d\n", links)
  @printf(io, "Lines          :%10d\n", lines)
  @printf(io, "Trafos         :%10d\n", trafos)
  @printf(io, "Generators     :%10d\n", gens)
  @printf(io, "Loads          :%10d\n", loads)
  @printf(io, "Shunts         :%10d\n", shunts)

  num_q_limit = length(net.qLimitEvents)
  @printf(io, "PV→PQ (Q-Limit):%10d\n", num_q_limit)

  println(io, "\n", totalLosses)

  @printf(io, "====================================================================================================================================================================================================\n")
  @printf(io, "| %-5s | %-20s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-12s |\n", "Nr", "Bus", "Vn [kV]", "V [kV]", "V [pu]", "phi [deg]", "Pg [MW]", "Qg [MVar]", "Pl [MW]", "Ql [MVar]", "Ps [MW]", "Qs [MVar]", "Type", "Control")
  @printf(io, "====================================================================================================================================================================================================\n")

  pGS = qGS = pLS = qLS = ""
  tpGS = tqGS = tpLS = tqLS = 0.0
  pShunt_str = qShunt_str = ""
  tpShunt = tqShunt = 0.0
  for n in nodes
    p_gen, q_gen, p_load, q_load = _effective_bus_power_components(net, n.busIdx)
    if abs(p_gen) > 1e-6
      pGS = @sprintf("%10.3f", p_gen)
      tpGS += p_gen
    else
      pGS = ""
    end
    if abs(q_gen) > 1e-6
      qGS = @sprintf("%10.3f", q_gen)
      tqGS += q_gen
    else
      qGS = ""
    end

    if abs(p_load) > 1e-6
      pLS = @sprintf("%10.3f", p_load)
      tpLS += p_load
    else
      pLS = ""
    end
    if abs(q_load) > 1e-6
      qLS = @sprintf("%10.3f", q_load)
      tqLS += q_load
    else
      qLS = ""
    end
    if !isnothing(n._pShunt) && abs(n._pShunt) > 1e-6
      pShunt_str = @sprintf("%10.3f", n._pShunt)
      tpShunt += n._pShunt
    else
      pShunt_str = ""
    end
    if !isnothing(n._qShunt) && abs(n._qShunt) > 1e-6
      qShunt_str = @sprintf("%10.3f", n._qShunt)
      tqShunt += n._qShunt
    else
      qShunt_str = ""
    end
    typeStr = toString(n._nodeType)
    controlStr = _control_label(net, n.busIdx)

    # Mark PV→PQ buses (hit Q-limit) with a star in the Type column
    if haskey(net.qLimitEvents, n.busIdx)
      typeStr *= "*"
    end

    v = n.comp.cVN * n._vm_pu
    nodeName = get(busNameByIdx, n.busIdx, n.comp.cName)
    if !isnothing(n._vmin_pu) && !isnothing(n._vmax_pu)
      if !isIsolated(n) && (n._vm_pu < n._vmin_pu || n._vm_pu > n._vmax_pu)
        nodeName *= " !"
      end
    end

    @printf(io, "| %-5d | %-20s | %-10.1f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-12s |\n", n.busIdx, nodeName, n.comp.cVN, v, n._vm_pu, n._va_deg, pGS, qGS, pLS, qLS, pShunt_str, qShunt_str, typeStr, controlStr)
  end

  @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  println(io, flowResults)


  if !isempty(net.linkVec)
    @printf(io, "\n----------------------------------- Link Flows (KCL) -----------------------------------\n")
    @printf(io, "| %-5s | %-8s | %-8s | %-6s | %-12s | %-12s | %-12s | %-12s |\n", "Nr", "From", "To", "Stat", "P [MW]", "Q [MVar]", "Ifrom [kA]", "Ito [kA]")
    @printf(io, "----------------------------------------------------------------------------------------------------------------\n")
    for l in net.linkVec
      p = isnothing(l.pFlow_MW) ? NaN : l.pFlow_MW
      q = isnothing(l.qFlow_MVar) ? NaN : l.qFlow_MVar
      ifrom = isnothing(l.iFrom_kA) ? NaN : l.iFrom_kA
      ito = isnothing(l.iTo_kA) ? NaN : l.iTo_kA
      fromName = get(busNameByIdx, Int(l.fromBus), string(l.fromBus))
      toName = get(busNameByIdx, Int(l.toBus), string(l.toBus))
      @printf(io, "| %-5d | %-8s | %-8s | %-6d | %-12.3f | %-12.3f | %-12.4f | %-12.4f |\n", l.linkIdx, fromName, toName, l.status, p, q, ifrom, ito)
    end
  end
  if toFile
    close(io)
    #println("Results have been written to $(joinpath(path, filename))")
  end
end

function formatProsumerResults(net::Net)
  buf = IOBuffer()

  # Rebuild mapping: bus index -> bus name
  busNameByIdx = _bus_name_by_idx(net)

  # Collect indices per bus, separated into generators and loads
  gens_by_bus  = Dict{Int,Vector{Int}}()
  loads_by_bus = Dict{Int,Vector{Int}}()

  for (idx, ps) in enumerate(net.prosumpsVec)
    bus = getPosumerBusIndex(ps)
    if isGenerator(ps)
      push!(get!(gens_by_bus, bus, Int[]), idx)
    else
      push!(get!(loads_by_bus, bus, Int[]), idx)
    end
  end

  # Helper: compute generator status from Q and minQ/maxQ
  status_from_Q = function (ps, qres)
    # default status
    status = "ok"
    isnothing(qres) && return status

    # tolerance for limit detection
    tol = 1e-6

    if !isnothing(ps.maxQ) && qres >= ps.maxQ - tol
      status = "Q-max limit"
    elseif !isnothing(ps.minQ) && qres <= ps.minQ + tol
      status = "Q-min limit"
    end

    return status
  end

  # =========================
  # Generator section
  # =========================
  println(buf, "\nGenerator results:")
  println(buf, "────────────────────────────────────────────────────────")
  @printf(buf, "%-8s %4s %12s %12s   %-14s\n", "Bus", "Gen#", "P [MW]", "Q [MVar]", "Status")
  println(buf, "────────────────────────────────────────────────────────")

  for bus in sort(collect(keys(gens_by_bus)))
    gens_idx = gens_by_bus[bus]

    # optional: sort generators at a bus by component name
    sort!(gens_idx, by = i -> net.prosumpsVec[i].comp.cName)

    busName = get(busNameByIdx, bus, "Bus_$bus")

    for (k, idx) in enumerate(gens_idx)
      ps = net.prosumpsVec[idx]

      # Prefer results (pRes/qRes); fall back to spec values
      p = ps.pRes === nothing ? ps.pVal : ps.pRes
      q = ps.qRes === nothing ? ps.qVal : ps.qRes

      status = status_from_Q(ps, q)

      @printf(buf, "%-8s %4d %12.3f %12.3f   %-14s\n", busName, k, p === nothing ? 0.0 : p, q === nothing ? 0.0 : q, status)
    end
  end

  println(buf, "────────────────────────────────────────────────────────")

  # =========================
  # Load section
  # =========================
  println(buf, "\nLoad results:")
  println(buf, "────────────────────────────────────────")
  @printf(buf, "%-8s %5s %12s %12s\n", "Bus", "Load#", "P [MW]", "Q [MVar]")
  println(buf, "────────────────────────────────────────")

  if !isnothing(loads_by_bus)
    for bus in sort(collect(keys(loads_by_bus)))
      try
        loads_idx = loads_by_bus[bus]

        # optional: sort loads at a bus by component name
        sort!(loads_idx, by = i -> net.prosumpsVec[i].comp.cName)

        busName = get(busNameByIdx, bus, "Bus_$bus")

        for (k, idx) in enumerate(loads_idx)
          ps = net.prosumpsVec[idx]

          # Prefer results (pRes/qRes); fall back to spec values
          p = ps.pRes === nothing ? ps.pVal : ps.pRes
          q = ps.qRes === nothing ? ps.qVal : ps.qRes

          @printf(buf, "%-8s %5d %12.3f %12.3f\n", busName, k, p === nothing ? 0.0 : p, q === nothing ? 0.0 : q)
        end
      catch e
        @warn "Error formatting load results for bus index $bus: $e"
      end
    end
  end

  println(buf, "────────────────────────────────────────────────────────")

  return String(take!(buf))
end

function printProsumerResults(net::Net)
  prosText = formatProsumerResults(net)
  println(prosText)
end
