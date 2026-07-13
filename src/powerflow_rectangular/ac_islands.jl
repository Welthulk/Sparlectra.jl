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

const AC_ISLAND_DISABLED_MESSAGE = "The AC network contains multiple disconnected islands after MATPOWER dcline import. DC lines are modeled as fixed terminal injections and do not connect the Ybus. Enable island-wise power-flow solving or provide one valid reference bus per island."

function _bus_display_name(net::Net, bus::Int)::String
  return get(net.busOriginalNameDict, bus, getCompName(net.nodeVec[bus].comp))
end

function _active_ac_island_components(net::Net)
  n = length(net.nodeVec)
  parent = collect(1:n)
  find_root(i::Int) = begin
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    i
  end
  function union_set(a::Int, b::Int)
    ra = find_root(a)
    rb = find_root(b)
    ra == rb && return
    ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)
  end
  for br in net.branchVec
    br.status == 1 || continue
    f = Int(br.fromBus)
    t = Int(br.toBus)
    1 <= f <= n && 1 <= t <= n || continue
    getNodeType(net.nodeVec[f]) == Isolated && continue
    getNodeType(net.nodeVec[t]) == Isolated && continue
    union_set(f, t)
  end
  groups = Dict{Int,Vector{Int}}()
  for bus in 1:n
    getNodeType(net.nodeVec[bus]) == Isolated && continue
    push!(get!(groups, find_root(bus), Int[]), bus)
  end
  return sort!(collect(values(groups)); by = xs -> minimum(xs))
end

function _dcline_terminal_rows_by_bus(net::Net)
  by_bus = Dict{Int,Vector{NamedTuple}}()
  bus_by_orig = Dict(orig => bus for (bus, orig) in net.busOrigIdxDict)
  for row in net.matpowerDclineMetadata
    fb = get(bus_by_orig, Int(row.from_bus), 0)
    tb = get(bus_by_orig, Int(row.to_bus), 0)
    fb > 0 && push!(get!(by_bus, fb, NamedTuple[]), merge(row, (terminal = :from, p_mw = -Float64(row.pf_mw), q_mvar = Float64(row.qf_mvar))))
    tb > 0 && push!(get!(by_bus, tb, NamedTuple[]), merge(row, (terminal = :to, p_mw = Float64(row.effective_pt_mw), q_mvar = Float64(row.qt_mvar))))
  end
  return by_bus
end

function detect_ac_islands(net::Net)
  refreshBusTypesFromProsumers!(net)
  components = _active_ac_island_components(net)
  terminal_by_bus = _dcline_terminal_rows_by_bus(net)
  rows = NamedTuple[]
  bus_to_island = Dict{Int,Int}()
  for (island_id, buses) in enumerate(components)
    busset = Set(buses)
    for bus in buses
      bus_to_island[bus] = island_id
    end
    branches = [i for (i, br) in enumerate(net.branchVec) if br.status == 1 && Int(br.fromBus) in busset && Int(br.toBus) in busset]
    prosumers = [ps for ps in net.prosumpsVec if Int(ps.comp.cFrom_bus) in busset]
    generators = [ps for ps in prosumers if isGenerator(ps)]
    loads = [ps for ps in prosumers if !isGenerator(ps)]
    refs = [bus for bus in buses if getNodeType(net.nodeVec[bus]) == Slack]
    pvs = [bus for bus in buses if getNodeType(net.nodeVec[bus]) == PV]
    dc_terms = reduce(vcat, [get(terminal_by_bus, bus, NamedTuple[]) for bus in buses]; init = NamedTuple[])
    total_load_p = sum(something(net.nodeVec[bus]._pƩLoad, 0.0) for bus in buses)
    total_load_q = sum(something(net.nodeVec[bus]._qƩLoad, 0.0) for bus in buses)
    total_gen_p = sum(something(net.nodeVec[bus]._pƩGen, 0.0) for bus in buses)
    total_gen_q = sum(something(net.nodeVec[bus]._qƩGen, 0.0) for bus in buses)
    total_dc_p = sum((Float64(t.p_mw) for t in dc_terms); init = 0.0)
    total_dc_q = sum((Float64(t.q_mvar) for t in dc_terms); init = 0.0)
    chosen = !isempty(refs) ? minimum(refs) : (!isempty(pvs) ? minimum(pvs) : 0)
    status = !isempty(refs) ? "has_ref" : (!isempty(pvs) ? "promote_pv_ref" : "missing_ref")
    note = !isempty(refs) ? "" : (!isempty(pvs) ? "matpower_like will promote PV bus $(chosen) as island angle reference" : "no REF/Slack or PV bus available")
    push!(rows, (
      island_id = island_id,
      buses = buses,
      branches = branches,
      n_bus = length(buses),
      n_branch = length(branches),
      n_pq = count(bus -> getNodeType(net.nodeVec[bus]) == PQ, buses),
      n_pv = length(pvs),
      n_ref = length(refs),
      n_generator = length(generators),
      n_load = length(loads),
      dc_terminal_count = length(dc_terms),
      dc_terminal_buses = [t.terminal === :from ? Int(t.from_bus) : Int(t.to_bus) for t in dc_terms],
      total_load_p_mw = total_load_p,
      total_gen_p_mw = total_gen_p,
      total_dcline_p_mw = total_dc_p,
      total_load_q_mvar = total_load_q,
      total_gen_q_mvar = total_gen_q,
      total_dcline_q_mvar = total_dc_q,
      imbalance_p_mw = total_gen_p - total_load_p,
      has_ref = !isempty(refs),
      chosen_ref_bus = chosen,
      status = status,
      note = note,
    ))
  end
  return (rows = rows, bus_to_island = bus_to_island, terminal_by_bus = terminal_by_bus)
end

function write_ac_island_report(path::AbstractString, report)
  rows = report.rows
  cols = (:island_id, :n_bus, :n_branch, :n_pq, :n_pv, :n_ref, :n_generator, :n_load, :dc_terminal_count, :dc_terminal_buses, :total_load_p_mw, :total_gen_p_mw, :total_dcline_p_mw, :total_load_q_mvar, :total_gen_q_mvar, :total_dcline_q_mvar, :imbalance_p_mw, :has_ref, :chosen_ref_bus, :status, :note)
  _write_namedtuple_csv(path, rows, cols; format = "technical")
  return path
end

function _print_ac_island_summary(report)
  println("AC island detection:")
  println("  islands: ", length(report.rows))
  for row in report.rows
    println("  island ", row.island_id, ": buses=", row.n_bus, ", branches=", row.n_branch, ", ref=", row.chosen_ref_bus == 0 ? "none" : row.chosen_ref_bus, ", status=", row.status)
  end
end

function _validate_island_references!(report)
  bad = [row for row in report.rows if row.chosen_ref_bus == 0]
  isempty(bad) && return
  ids = join(getfield.(bad, :island_id), ", ")
  error("AC island reference validation failed: island(s) $(ids) have no REF/Slack bus and no PV/voltage-controlled generator available for matpower_like reference selection.")
end

function _prepare_island_net(net::Net, row)
  inet = deepcopy(net)
  busmap = Dict(old => new for (new, old) in enumerate(row.buses))
  busset = Set(row.buses)
  inet.nodeVec = [deepcopy(net.nodeVec[old]) for old in row.buses]
  empty!(inet.busDict)
  empty!(inet.busOrigIdxDict)
  empty!(inet.busOriginalNameDict)
  for (new, old) in enumerate(row.buses)
    node = inet.nodeVec[new]
    node.busIdx = new
    inet.busDict[getCompName(node.comp)] = new
    inet.busOrigIdxDict[new] = get(net.busOrigIdxDict, old, old)
    if haskey(net.busOriginalNameDict, old)
      inet.busOriginalNameDict[new] = net.busOriginalNameDict[old]
    end
  end
  inet.branchVec = [deepcopy(net.branchVec[i]) for i in row.branches]
  for (newidx, br) in enumerate(inet.branchVec)
    br.branchIdx = newidx
    br.fromBus = busmap[Int(br.fromBus)]
    br.toBus = busmap[Int(br.toBus)]
    hasproperty(br.comp, :cFrom_bus) && (br.comp.cFrom_bus = Int(br.fromBus))
    hasproperty(br.comp, :cTo_bus) && (br.comp.cTo_bus = Int(br.toBus))
  end
  inet.prosumpsVec = [deepcopy(ps) for ps in net.prosumpsVec if Int(ps.comp.cFrom_bus) in busset]
  for ps in inet.prosumpsVec
    newbus = busmap[Int(ps.comp.cFrom_bus)]
    ps.comp.cFrom_bus = newbus
    ps.comp.cTo_bus = newbus
  end
  inet.shuntVec = [deepcopy(sh) for sh in net.shuntVec if Int(sh.busIdx) in busset]
  for sh in inet.shuntVec
    sh.busIdx = busmap[Int(sh.busIdx)]
    hasproperty(sh.comp, :cFrom_bus) && (sh.comp.cFrom_bus = sh.busIdx)
    hasproperty(sh.comp, :cTo_bus) && (sh.comp.cTo_bus = sh.busIdx)
  end
  inet.slackVec = [busmap[b] for b in net.slackVec if b in busset]
  empty!(inet.isoNodes)
  if row.n_ref == 0 && row.chosen_ref_bus > 0
    chosen = busmap[row.chosen_ref_bus]
    inet.nodeVec[chosen]._nodeType = Slack
    push!(inet.slackVec, chosen)
    for ps in inet.prosumpsVec
      if Int(ps.comp.cFrom_bus) == chosen && isGenerator(ps) && isRegulating(ps)
        ps.referencePri = chosen
        break
      end
    end
  end
  return inet
end

function _sync_island_solution!(net::Net, inet::Net, row)
  for (new, bus) in enumerate(row.buses)
    net.nodeVec[bus]._vm_pu = inet.nodeVec[new]._vm_pu
    net.nodeVec[bus]._va_deg = inet.nodeVec[new]._va_deg
  end
  for (new, bridx) in enumerate(row.branches)
    net.branchVec[bridx].fBranchFlow = inet.branchVec[new].fBranchFlow
    net.branchVec[bridx].tBranchFlow = inet.branchVec[new].tBranchFlow
    net.branchVec[bridx].pLosses = inet.branchVec[new].pLosses
    net.branchVec[bridx].qLosses = inet.branchVec[new].qLosses
  end
end
