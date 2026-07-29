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

# =============================================================================
# AC island handling.
#
# Two different notions of "island" live in this file, and mixing them up has
# already produced one real bug:
#
#   * _active_ac_island_components -- connectivity over in-service BRANCHES
#     only. This is the solver's view: by the time the solver runs, link
#     clusters have been contracted by _merged_pf_net, so switches no longer
#     exist as separate objects.
#
#   * electricalIslandComponents -- connectivity over branches AND closed
#     links. This is the physical view on a raw, not-yet-merged net, e.g.
#     straight out of the CGMES importer where retained switches are links.
#
# On a link-merged net both agree. On a raw net they do not: a CGMES delivery
# with retained switches shows more branch-only "islands" than it physically
# has (Svedala: 8 vs 6). Always ask which net you are holding before picking
# one of the two.
# =============================================================================

const AC_ISLAND_DISABLED_MESSAGE = "The AC network contains multiple disconnected islands after MATPOWER dcline import. DC lines are modeled as fixed terminal injections and do not connect the Ybus. Enable island-wise power-flow solving or provide one valid reference bus per island."

"""Human-readable bus label: the imported original name if we kept one, else the component name."""
function _bus_display_name(net::Net, bus::Int)::String
  return get(net.busOriginalNameDict, bus, getCompName(net.nodeVec[bus].comp))
end

# Connected components over in-service branches only -- the solver's view of
# the network. See the file header for why this differs from the electrical
# view on an unmerged net.
function _active_ac_island_components(net::Net)
  n = length(net.nodeVec)
  # Union-find over bus indices; parent[i] == i marks a root.
  parent = collect(1:n)
  find_root(i::Int) = begin
    # Path halving: every second node on the way up is relinked to its
    # grandparent, which keeps the trees flat without a second pass.
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
    # Always hang the larger root under the smaller one. Union-by-size would
    # be marginally faster, but this makes the representative deterministic
    # (always the smallest bus index of the component), which several callers
    # and the test suite rely on.
    ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)
  end
  for br in net.branchVec
    br.status == 1 || continue
    f = Int(br.fromBus)
    t = Int(br.toBus)
    # Defensive: branches may carry stale indices after topology edits.
    1 <= f <= n && 1 <= t <= n || continue
    # De-energized buses are not part of any island -- a branch that ends in
    # one must not drag its other end into a component.
    getNodeType(net.nodeVec[f]) == Isolated && continue
    getNodeType(net.nodeVec[t]) == Isolated && continue
    union_set(f, t)
  end
  groups = Dict{Int,Vector{Int}}()
  for bus in 1:n
    getNodeType(net.nodeVec[bus]) == Isolated && continue
    push!(get!(groups, find_root(bus), Int[]), bus)
  end
  # Sorted by smallest member so island numbering is stable across runs --
  # island IDs end up in diagnostic artifacts and error messages.
  return sort!(collect(values(groups)); by = xs -> minimum(xs))
end

"""
    electricalIslandComponents(net) -> Vector{Vector{Int}}

Connected bus components of `net` counting **both** in-service branches and
closed links as connections — the electrical view.

`_active_ac_island_components` deliberately walks `branchVec` only, because
the solver runs on a net whose link clusters are already contracted. On a raw
imported net that is misleading: parts joined by a closed switch appear as
separate islands (a CGMES delivery with retained switches can show several
"islands" that are one electrical island). Use this helper whenever islands
are judged on a net that has not been link-merged — diagnostics, and the
importer's slack decisions.
"""
function electricalIslandComponents(net::Net)::Vector{Vector{Int}}
  n = length(net.nodeVec)
  parent = collect(1:n)
  find_root(i::Int) = begin
    while parent[i] != i
      parent[i] = parent[parent[i]]
      i = parent[i]
    end
    i
  end
  # Same union-find as _active_ac_island_components, but the bounds and
  # Isolated guards sit inside union_set: this variant is fed from two vectors
  # (branches and links), so filtering once here beats duplicating it twice.
  function union_set(a::Int, b::Int)
    (1 <= a <= n && 1 <= b <= n) || return
    getNodeType(net.nodeVec[a]) == Isolated && return
    getNodeType(net.nodeVec[b]) == Isolated && return
    ra = find_root(a)
    rb = find_root(b)
    ra == rb && return
    ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)
  end
  for br in net.branchVec
    br.status == 1 && union_set(Int(br.fromBus), Int(br.toBus))
  end
  # The one line that separates this from the branch-only view: a closed link
  # (busbar coupler, retained CGMES switch) is a galvanic connection, so its
  # two buses are electrically one node.
  for l in net.linkVec
    l.status == 1 && union_set(Int(l.fromBus), Int(l.toBus))
  end
  groups = Dict{Int,Vector{Int}}()
  for bus = 1:n
    getNodeType(net.nodeVec[bus]) == Isolated && continue
    push!(get!(groups, find_root(bus), Int[]), bus)
  end
  return sort!(collect(values(groups)); by = xs -> minimum(xs))
end

"""Bus index → electrical island number (see [`electricalIslandComponents`](@ref))."""
function electricalIslandOfBus(net::Net)::Dict{Int,Int}
  out = Dict{Int,Int}()
  for (i, comp) in enumerate(electricalIslandComponents(net))
    for b in comp
      out[b] = i
    end
  end
  return out
end

# MATPOWER DC lines are not part of the Ybus -- they are modelled as a fixed
# injection at each of their two terminals. That is precisely why they can
# leave the AC network split into islands while the case still looks connected
# to a user. Group those injections per bus so the island report can show how
# much power enters or leaves an island through DC.
function _dcline_terminal_rows_by_bus(net::Net)
  by_bus = Dict{Int,Vector{NamedTuple}}()
  # MATPOWER metadata refers to original case bus numbers, not our indices.
  bus_by_orig = Dict(orig => bus for (bus, orig) in net.busOrigIdxDict)
  for row in net.matpowerDclineMetadata
    fb = get(bus_by_orig, Int(row.from_bus), 0)
    tb = get(bus_by_orig, Int(row.to_bus), 0)
    # Sign convention: pf_mw is the power *leaving* the from-terminal, so from
    # the bus's perspective it is a negative injection. The to-terminal uses
    # effective_pt_mw, i.e. pf minus the line losses.
    fb > 0 && push!(get!(by_bus, fb, NamedTuple[]), merge(row, (terminal = :from, p_mw = -Float64(row.pf_mw), q_mvar = Float64(row.qf_mvar))))
    tb > 0 && push!(get!(by_bus, tb, NamedTuple[]), merge(row, (terminal = :to, p_mw = Float64(row.effective_pt_mw), q_mvar = Float64(row.qt_mvar))))
  end
  return by_bus
end

"""
Build the per-island report the solver, the diagnostics writer and the CSV
artifact all share. Purely descriptive -- it decides nothing, it only records
what each island contains and which reference it *would* get.
"""
function detect_ac_islands(net::Net)
  # Bus types are derived from prosumers, and the caller may have edited the
  # net (merges, de-energization) since they were last computed. Refresh first,
  # otherwise refs/pvs below are read from stale types.
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
    # Reference selection, MATPOWER-like: an existing Slack wins; otherwise the
    # island borrows a PV bus as its angle reference; otherwise it has none and
    # _validate_island_references! will reject the whole run. `minimum` keeps
    # the choice deterministic. chosen == 0 encodes "no candidate".
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

# Fail fast and for ALL bad islands at once: an island without any angle
# reference has an underdetermined system, and finding that out one island at a
# time is needlessly slow for the user.
function _validate_island_references!(report)
  bad = [row for row in report.rows if row.chosen_ref_bus == 0]
  isempty(bad) && return
  ids = join(getfield.(bad, :island_id), ", ")
  error("AC island reference validation failed: island(s) $(ids) have no REF/Slack bus and no PV/voltage-controlled generator available for matpower_like reference selection.")
end

"""
Cut one island out of `net` as a standalone `Net` that the solver can run on
unchanged. Every bus index is renumbered to 1..n_island, so all index-carrying
containers have to be rebuilt consistently -- missing one of them yields a net
that looks fine and solves the wrong topology.
"""
function _prepare_island_net(net::Net, row)
  # Start from a full copy to inherit scalar settings (baseMVA, tolerances,
  # flags); every index-bearing vector below is then replaced wholesale.
  inet = deepcopy(net)
  # old (net-wide) bus index -> new (island-local) bus index
  busmap = Dict(old => new for (new, old) in enumerate(row.buses))
  busset = Set(row.buses)
  inet.nodeVec = [deepcopy(net.nodeVec[old]) for old in row.buses]
  # The name/index dictionaries were copied from the full net and are now
  # wrong in both directions -- rebuild them from scratch.
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
  # row.branches already holds only in-service branches with both ends inside
  # the island, so no further filtering is needed here. The component objects
  # carry their own bus references and must be updated alongside the branch.
  inet.branchVec = [deepcopy(net.branchVec[i]) for i in row.branches]
  for (newidx, br) in enumerate(inet.branchVec)
    br.branchIdx = newidx
    br.fromBus = busmap[Int(br.fromBus)]
    br.toBus = busmap[Int(br.toBus)]
    hasproperty(br.comp, :cFrom_bus) && (br.comp.cFrom_bus = Int(br.fromBus))
    hasproperty(br.comp, :cTo_bus) && (br.comp.cTo_bus = Int(br.toBus))
  end
  # Prosumers and shunts are bus-local, so island membership follows from their
  # bus alone. cTo_bus is set to the same bus: for a bus-local element both
  # terminals are the same node.
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
  # Only slacks of *this* island survive, renumbered. Slacks of other islands
  # are dropped, not moved.
  inet.slackVec = [busmap[b] for b in net.slackVec if b in busset]
  # Isolated buses were filtered out when the components were built, so the
  # island net starts without any.
  empty!(inet.isoNodes)

  # MATPOWER-like reference promotion: an island that brought no Slack of its
  # own gets one PV bus turned into the angle reference. detect_ac_islands has
  # already picked the bus; here it is actually applied.
  if row.n_ref == 0 && row.chosen_ref_bus > 0
    chosen = busmap[row.chosen_ref_bus]
    inet.nodeVec[chosen]._nodeType = Slack
    push!(inet.slackVec, chosen)
    # Bus typing is re-derived from prosumers downstream, so the node type
    # alone would not survive: mark the regulating generator at that bus as the
    # reference too. First match wins -- one reference per bus is enough.
    for ps in inet.prosumpsVec
      if Int(ps.comp.cFrom_bus) == chosen && isGenerator(ps) && isRegulating(ps)
        ps.referencePri = chosen
        break
      end
    end
  end
  return inet
end

# Write one island's solution back into the full net. The inverse of the
# renumbering done in _prepare_island_net: row.buses/row.branches map the
# island-local position back to the net-wide index. Only results are copied --
# the island net's topology is discarded.
#
# Note the angles are NOT shifted to a common reference: each island has its
# own angle reference and its angles are only meaningful within that island.
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
