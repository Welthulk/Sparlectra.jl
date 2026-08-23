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
# purpose: functions for formatting and printing results of power flow
#          calculations, including the structured ACPFlowReport container
"""
    _fitColumn(text, width) -> String

Trim `text` to `width` characters, marking the cut with `…`. Fixed-width
`@sprintf` fields pad but never truncate, so long names — CGMES bus and
branch identifiers routinely exceed 25 characters — would otherwise push
every following column out of alignment.
"""
function _fitColumn(text, width::Int)::String
  str = string(text)
  length(str) <= width && return str
  width <= 1 && return first(str, width)
  return first(str, width - 1) * "…"
end

# Column layout of the HVDC Link Flows table (parallel Phase 6): widths are
# DERIVED from the canonical mode/status sets in hvdc_pair_control.jl, and
# header, separator, and rows are generated from the same tuple, so a new
# mode or status value can never misalign the table again.
const _HVDC_COL_MODE = max(maximum(length(string(m)) for m in HVDC_PAIR_MODES), length("fixed"), length("Mode"))
const _HVDC_COL_STATUS = max(maximum(length(string(s)) for s in HVDC_PAIR_STATUS_VALUES), length("Status"))
const _HVDC_COL_BUS = 12
const _HVDC_TABLE_COLUMNS = (("Nr", 5), ("Name", 14), ("From", _HVDC_COL_BUS), ("To", _HVDC_COL_BUS), ("Mode", _HVDC_COL_MODE), ("P_from [MW]", 11), ("P_to [MW]", 9), ("Loss [MW]", 9), ("Q_from [MVar]", 13), ("Q_to [MVar]", 11), ("Rating", 7), ("Status", _HVDC_COL_STATUS))

# every cell passes through _fitColumn with its column width; no raw %-Ns on
# free-form or enum strings
_hvdc_table_line(cells) = string("| ", join((rpad(_fitColumn(string(c), w), w) for (c, (_, w)) in zip(cells, _HVDC_TABLE_COLUMNS)), " | "), " |")
_hvdc_table_header() = _hvdc_table_line(Tuple(name for (name, _) in _HVDC_TABLE_COLUMNS))

function format_version(version::VersionNumber)
  major = lpad(version.major, 1, '0')
  minor = lpad(version.minor, 1, '0')
  patch = lpad(version.patch, 0, '0')
  return "$major.$minor.$patch"
end

"""
    ACPFlowReport

Structured container for AC power flow results.

The vectors (`nodes`, `branches`, `links`, `transformer_controls`, `q_limit_events`)
are table-like and can
be converted directly to `DataFrame`s if `DataFrames.jl` is available, e.g.
`DataFrame(report.nodes)`.

# Fields
- `metadata`: Global run/case metadata (solver, tolerance, elapsed time, losses, ...).
- `nodes`: Per-bus electrical state and power balance values.
- `branches`: Per-branch directional flows and losses.
- `links`: Link-flow values from KCL post-processing.
- `transformer_controls`: Tap-controller state rows with typed `missing` for
  non-applicable engineering values.
- `q_limit_events`: PV→PQ limit-hit markers.
- `hvdc_links`: One row per HVDC link (`net.hvdcLinks`): terminal flows,
  loss, mode, rating, controller status (see `_hvdc_link_flow_rows`).
"""
struct ACPFlowReport
  metadata::NamedTuple
  nodes::Vector{NamedTuple}
  branches::Vector{NamedTuple}
  links::Vector{NamedTuple}
  transformer_controls::Vector{NamedTuple}
  q_limit_events::Vector{NamedTuple}
  hvdc_links::Vector{NamedTuple}
end

"""
    _build_transformer_control_rows(net::Net)

Internal helper that mirrors transformer control state into table-like rows for
`ACPFlowReport`.

Notes:
- Rows are intentionally typed for DataFrame conversion.
- Non-applicable controller fields remain `missing` (not placeholder strings).
"""
function _build_transformer_control_rows(net::Net)::Vector{NamedTuple}
  # Reuse tap-control module helper to keep classic and structured reports aligned.
  return buildTapControllerReportRows(net)
end

function Base.show(io::IO, report::ACPFlowReport)
  print(
    io,
    "ACPFlowReport(",
    "case=",
    report.metadata.case,
    ", converged=",
    report.metadata.converged,
    ", nodes=",
    length(report.nodes),
    ", branches=",
    length(report.branches),
    ", links=",
    length(report.links),
    ", transformer_controls=",
    length(report.transformer_controls),
    ", q_limit_events=",
    length(report.q_limit_events),
    ", hvdc_links=",
    length(report.hvdc_links),
    ")",
  )
end

_default0(x) = isnothing(x) ? 0.0 : x

function _bus_name_by_idx(net::Net)
  busNameByIdx = Dict{Int,String}()
  for (name, idx) in net.busDict
    busNameByIdx[idx] = name
  end
  return busNameByIdx
end

_effective_bus_name(busNameByIdx::AbstractDict, net::Net, bus_idx::Int)::String = get(busNameByIdx, bus_idx, string(bus_idx))

# External-grid source presentation (issue #299): the hidden internal
# reference bus "<bus>__extgrid_int" reports type "SOURCE" instead of
# "SLACK" — the grid connection is a source, the hidden bus only its
# internal anchor. The solver-internal node type stays Slack; this is
# presentation only (result table, structured report, CSV/JSON rows).
_is_external_grid_internal_bus(name::AbstractString)::Bool = endswith(name, "__extgrid_int")

function _bus_type_label(nodeType, name::AbstractString)::String
  nodeType == Sparlectra.Slack && _is_external_grid_internal_bus(name) && return "SOURCE"
  return toString(nodeType)
end

# One-line "slack or source" statement for the result header: names every
# reference of the net, and for an external-grid source also its feeder data
# (Sk''/R-X recovered from the native feeder record) and the internal slack
# bus, so run.log states the chosen connection model explicitly.
function _grid_connection_summary(net::Net, busNameByIdx::AbstractDict)::String
  parts = String[]
  for busIdx in net.slackVec
    name = _effective_bus_name(busNameByIdx, net, busIdx)
    if _is_external_grid_internal_bus(name)
      terminal = name[1:(end-length("__extgrid_int"))]
      vn = getNodeVn(net.nodeVec[busIdx])
      sk = NaN
      rx = NaN
      for e in net.sc_sources.external_network_injections
        (e.bus == terminal && e.maxInitialSymShCCurrent_A !== nothing) || continue
        sk = sqrt(3.0) * vn * Float64(e.maxInitialSymShCCurrent_A) / 1000.0
        rx = e.maxR1ToX1Ratio === nothing ? NaN : Float64(e.maxR1ToX1Ratio)
        break
      end
      data = isfinite(sk) ? @sprintf("Sk'' = %.1f MVA, R/X = %s; ", sk, isfinite(rx) ? @sprintf("%.3g", rx) : "?") : ""
      push!(parts, "external-grid source at $(terminal) ($(data)internal slack: $(name))")
    else
      push!(parts, "slack bus $(name)")
    end
  end
  isempty(parts) && return "no reference registered"
  return join(parts, "; ")
end

function _original_bus_name(busNameByIdx::AbstractDict, net::Net, bus_idx::Int)::String
  return get(net.busOriginalNameDict, bus_idx, _effective_bus_name(busNameByIdx, net, bus_idx))
end

function _original_branch_name(net::Net, br::Branch)::String
  meta = get(net.matpower_branch_metadata, br.branchIdx, nothing)
  if meta !== nothing && hasproperty(meta, :orig_name)
    name = getproperty(meta, :orig_name)
    name === nothing || isempty(String(name)) || return String(name)
  end
  return br.comp.cName
end

function _branch_kind_name(net::Net, br::Branch)::String
  meta = get(net.matpower_branch_metadata, br.branchIdx, nothing)
  if meta !== nothing
    hasproperty(meta, :orig_kind) && return String(getproperty(meta, :orig_kind))
    hasproperty(meta, :dtf_kind) && return String(getproperty(meta, :dtf_kind))
  end
  return string(br.comp.cTyp)
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

struct BusControlFlags
  has_qu::Bool
  has_pu::Bool
  # HVDC pair terminal role: :none, :terminal (steered converter injection),
  # :gridforming (island_feed reference at the converter PCC)
  hvdc::Symbol
end

const _NO_BUS_CONTROL_FLAGS = BusControlFlags(false, false, :none)

function _bus_control_flag_cache(net::Net)::Dict{Int,BusControlFlags}
  cache = Dict{Int,BusControlFlags}()
  for ps in net.prosumpsVec
    bus_idx = getPosumerBusIndex(ps)
    old = get(cache, bus_idx, _NO_BUS_CONTROL_FLAGS)
    cache[bus_idx] = BusControlFlags(old.has_qu | has_qu_controller(ps), old.has_pu | has_pu_controller(ps), old.hvdc)
  end
  # HVDC pair terminals: steered injections read "B2B", a grid-forming
  # reference (island_feed to side) reads "B2B src" — the bus table then
  # explains WHY a SLACK row is a converter and not a power plant
  for c in _hvdc_pair_controllers(net)
    c.enabled || continue
    fidx = geNetBusIdx(net = net, busName = c.from_bus)
    tidx = geNetBusIdx(net = net, busName = c.to_bus)
    fold = get(cache, fidx, _NO_BUS_CONTROL_FLAGS)
    cache[fidx] = BusControlFlags(fold.has_qu, fold.has_pu, :terminal)
    told = get(cache, tidx, _NO_BUS_CONTROL_FLAGS)
    cache[tidx] = BusControlFlags(told.has_qu, told.has_pu, c.mode == :island_feed ? :gridforming : :terminal)
  end
  return cache
end

@inline function _control_label(flags::BusControlFlags)::String
  parts = String[]
  flags.has_qu && push!(parts, "Q(U)")
  flags.has_pu && push!(parts, "P(U)")
  flags.hvdc == :terminal && push!(parts, "B2B")
  flags.hvdc == :gridforming && push!(parts, "B2B src")
  return isempty(parts) ? "-" : join(parts, ", ")
end

@inline _cached_control_label(cache::AbstractDict{Int,BusControlFlags}, bus_idx::Int)::String = _control_label(get(cache, bus_idx, _NO_BUS_CONTROL_FLAGS))

function _control_label(net::Net, bus_idx::Int)::String
  has_qu, has_pu = _bus_control_flags(net, bus_idx)
  return _control_label(BusControlFlags(has_qu, has_pu, :none))
end

function _tap_voltage_target_by_bus(net::Net)::Dict{Int,Float64}
  targets = Dict{Int,Float64}()
  for ctrl in _tap_controllers(net)
    ctrl.enabled || continue
    ctrl.mode in (:voltage, :voltage_and_branch_active_power) || continue
    isnothing(ctrl.target_bus) && continue
    isnothing(ctrl.target_vm_pu) && continue
    bus_idx = geNetBusIdx(net = net, busName = ctrl.target_bus)
    targets[bus_idx] = ctrl.target_vm_pu
  end
  return targets
end

function _controller_counts(net::Net)
  tap = count(c -> c.enabled, _tap_controllers(net))
  qu = count(ps -> has_qu_controller(ps), net.prosumpsVec)
  pu = count(ps -> has_pu_controller(ps), net.prosumpsVec)
  return tap, qu, pu
end

function _effective_bus_power_components(net::Net, bus_idx::Int)
  base = net.baseMVA
  node = net.nodeVec[bus_idx]

  vm = node._vm_pu
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

  if node._nodeType == Sparlectra.Slack
    if !isnothing(node._pƩGen)
      p_gen = node._pƩGen
    end
    if !isnothing(node._qƩGen)
      q_gen = node._qƩGen
    end
  elseif node._nodeType == Sparlectra.PV
    if abs(q_gen) <= 1e-9 && !isnothing(node._qƩGen)
      q_gen = node._qƩGen
    end
    if abs(p_gen) <= 1e-9 && !isnothing(node._pƩGen)
      p_gen = node._pƩGen
    end
  end

  return p_gen, q_gen, p_load, q_load
end

struct BusPowerComponents
  p_gen::Float64
  q_gen::Float64
  p_load::Float64
  q_load::Float64
end

Base.Tuple(components::BusPowerComponents) = (components.p_gen, components.q_gen, components.p_load, components.q_load)

function _bus_power_component_cache(net::Net)::Dict{Int,BusPowerComponents}
  base = net.baseMVA
  totals = Dict{Int,NTuple{4,Float64}}()
  for ps in net.prosumpsVec
    bus_idx = getPosumerBusIndex(ps)
    node = net.nodeVec[bus_idx]
    vm = node._vm_pu
    vm_safe = vm > 1e-9 ? vm : 1e-9
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
    old = get(totals, bus_idx, (0.0, 0.0, 0.0, 0.0))
    totals[bus_idx] = isGenerator(ps) ? (old[1] + p_mw, old[2] + q_mvar, old[3], old[4]) : (old[1], old[2], old[3] + p_mw, old[4] + q_mvar)
  end

  cache = Dict{Int,BusPowerComponents}()
  for node in net.nodeVec
    p_gen, q_gen, p_load, q_load = get(totals, node.busIdx, (0.0, 0.0, 0.0, 0.0))
    if node._nodeType == Sparlectra.Slack
      !isnothing(node._pƩGen) && (p_gen = node._pƩGen)
      !isnothing(node._qƩGen) && (q_gen = node._qƩGen)
    elseif node._nodeType == Sparlectra.PV
      abs(q_gen) <= 1e-9 && !isnothing(node._qƩGen) && (q_gen = node._qƩGen)
      abs(p_gen) <= 1e-9 && !isnothing(node._pƩGen) && (p_gen = node._pƩGen)
    end
    cache[node.busIdx] = BusPowerComponents(p_gen, q_gen, p_load, q_load)
  end
  return cache
end

function _bus_power_components(cache::AbstractDict{Int,BusPowerComponents}, bus_idx::Int)
  return Tuple(get(cache, bus_idx, BusPowerComponents(0.0, 0.0, 0.0, 0.0)))
end

# Printed Type column of the branch table. Consult the import metadata
# first so print and report rows (_branch_kind_name) agree for the same
# branch, then fall back to the name marker. The PST label is a heuristic
# (phase tap present, nonzero shift, no ratio tap) pending the typed
# BranchKind enum (tracked separately).
function _branch_kind_label(net::Net, br::Branch)::String
  meta = get(net.matpower_branch_metadata, br.branchIdx, nothing)
  kind = if meta !== nothing && hasproperty(meta, :orig_kind)
    Symbol(getproperty(meta, :orig_kind))
  elseif meta !== nothing && hasproperty(meta, :dtf_kind)
    Symbol(getproperty(meta, :dtf_kind))
  else
    nothing
  end
  label = if kind === :transformer
    "Trafo"
  elseif kind === :line
    "Line"
  else
    name = br.comp.cName
    if occursin("_ACL_", name)
      "Line"
    elseif occursin("_2WT_", name)
      "Trafo"
    elseif occursin("_PI_", name)
      "PI"
    else
      "Branch"
    end
  end
  if label == "Trafo" && br.has_phase_tap && !br.has_ratio_tap && abs(br.phase_shift_deg) > 0.0
    return "PST"
  end
  return label
end

"""
    buildACPFlowReport(net::Net; ...)

Builds a structured report object from solved network data.
This provides a machine-readable alternative to `printACPFlowResults`.
"""
function buildACPFlowReport(net::Net; ct::Float64 = 0.0, ite::Int = 0, tol::Float64 = 1e-6, converged::Bool = true, solver::Symbol = :NR)::ACPFlowReport
  rect_status = rectangular_pf_status(net)
  wrong_branch_status = _rect_status_get(rect_status, :wrong_branch_status, :not_checked)
  wrong_branch_reason = _rect_status_get(rect_status, :wrong_branch_reason, :not_checked)
  busNameByIdx = _bus_name_by_idx(net)
  nodes_sorted = sort(net.nodeVec, by = x -> x.busIdx)
  power_components = _bus_power_component_cache(net)
  control_labels = _bus_control_flag_cache(net)

  node_rows = NamedTuple[]
  for n in nodes_sorted
    p_gen, q_gen, p_load, q_load = _bus_power_components(power_components, n.busIdx)
    push!(
      node_rows,
      (
        bus = n.busIdx,
        bus_name = _effective_bus_name(busNameByIdx, net, n.busIdx),
        type = _bus_type_label(n._nodeType, _effective_bus_name(busNameByIdx, net, n.busIdx)),
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
        control = _cached_control_label(control_labels, n.busIdx),
        original_bus_name = _original_bus_name(busNameByIdx, net, n.busIdx),
      ),
    )
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

    from_name = _effective_bus_name(busNameByIdx, net, br.fromBus)
    to_name = _effective_bus_name(busNameByIdx, net, br.toBus)
    push!(
      branch_rows,
      (
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
        branch_name = br.comp.cName,
        original_branch_name = _original_branch_name(net, br),
        from_bus_name = from_name,
        to_bus_name = to_name,
        original_from_bus_name = _original_bus_name(busNameByIdx, net, br.fromBus),
        original_to_bus_name = _original_bus_name(busNameByIdx, net, br.toBus),
        branch_kind = _branch_kind_name(net, br),
        # r0.9.10 per-terminal state plus the open-end voltage result of a
        # one-sided open branch (missing when closed or fully open)
        terminal_state = _branch_terminal_state(br),
        open_end_vm_pu = br.open_end_vm_pu === nothing ? missing : br.open_end_vm_pu,
        open_end_va_deg = br.open_end_va_deg === nothing ? missing : br.open_end_va_deg,
      ),
    )
  end

  link_rows = NamedTuple[]
  for l in net.linkVec
    push!(link_rows, (link = l.cName, link_index = l.linkIdx, from_bus = l.fromBus, to_bus = l.toBus, status = l.status, p_MW = _default0(l.pFlow_MW), q_MVar = _default0(l.qFlow_MVar), ifrom_kA = _default0(l.iFrom_kA), ito_kA = _default0(l.iTo_kA)))
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
    wrong_branch_status = wrong_branch_status,
    wrong_branch_reason = wrong_branch_reason,
  )

  # Keep structured tap-controller data as a dedicated relation in the report
  # (preferred over widening the branch table with sparse controller columns).
  tap_control_rows = _build_transformer_control_rows(net)
  return ACPFlowReport(metadata, node_rows, branch_rows, link_rows, tap_control_rows, q_events, _hvdc_link_flow_rows(net))
end

function formatBranchResults(net::Net; max_rows::Union{Nothing,Int} = nothing)
  busNameByIdx = _bus_name_by_idx(net)
  ctrl_by_branch = Dict{Int,NamedTuple}()
  for row in buildTapControllerReportRows(net)
    ctrl_by_branch[row.transformer_branch_index] = row
  end
  # series-reactance (TCSC) controllers reuse the same Ctrl/P_tgt/Ctrl-status
  # columns; TapPos stays "-" (no discrete positions), the moved reactance is
  # in the TCSC summary block of the Control footer
  for c in _series_reactance_controllers(net)
    haskey(ctrl_by_branch, c.branch_idx) && continue
    ctrl_by_branch[c.branch_idx] = (control_type = :TCSC, p_target_mw = c.p_target_mw, ratio_tap_position = missing, phase_tap_position = missing, converged = c.converged, at_limit = c.at_limit, status = c.status)
  end
  # full UPFC (#326) marks its line branch with the same Ctrl/P_tgt/Ctrl-status
  # columns; P_tgt shows the active target (the Q target, series voltage, and
  # DC-link balance are in the UPFC summary block of the Control footer)
  for c in _upfc_full_controllers(net)
    haskey(ctrl_by_branch, c.branch_idx) && continue
    ctrl_by_branch[c.branch_idx] = (control_type = :UPFC, p_target_mw = c.p_target_mw, ratio_tap_position = missing, phase_tap_position = missing, converged = c.converged, at_limit = c.at_limit, status = c.status)
  end

  ctrl_status = function (row)
    if row.converged
      return "converged"
    elseif row.at_limit
      return "at_limit_not_converged"
    else
      return row.status
    end
  end
  #! format: off
  fr_io = IOBuffer()
  @printf(fr_io, "\n==========================================================================================================================================================================================================================\n")
  @printf(fr_io, "| %-25s | %-6s | %-25s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-9s | %-22s |\n", "Branch", "Type", "Connection", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]", "Ctrl", "P_tgt", "TapPos", "Ctrl status")
  @printf(fr_io, "==========================================================================================================================================================================================================================\n")
  #! format: on
  shown_rows = 0
  for br in net.branchVec
    if !isnothing(max_rows) && shown_rows >= max_rows
      break
    end
    from = br.fromBus
    to = br.toBus
    bName = br.comp.cName
    branchKind = _branch_kind_label(net, br)
    fromName = get(busNameByIdx, Int(from), string(from))
    toName = get(busNameByIdx, Int(to), string(to))
    connection = string(fromName, " -> ", toName)
    # only fully open branches print zeros; a one-sided open branch
    # (r0.9.10) carries its stored one-sided flows (closed side = charging
    # draw, open side = 0)
    terminal_state = _branch_terminal_state(br)
    if terminal_state == :open || (isnothing(br.fBranchFlow)) || (isnothing(br.tBranchFlow))
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
    ctrl = get(ctrl_by_branch, br.branchIdx, nothing)
    ctrl_type = isnothing(ctrl) ? "-" : String(ctrl.control_type)
    p_target = isnothing(ctrl) || ismissing(ctrl.p_target_mw) ? "-" : @sprintf("%.3f", ctrl.p_target_mw)
    tap_pos = "-"
    if !isnothing(ctrl)
      if !ismissing(ctrl.ratio_tap_position)
        tap_pos = @sprintf("%+d", ctrl.ratio_tap_position)
      elseif !ismissing(ctrl.phase_tap_position)
        tap_pos = @sprintf("%+d", ctrl.phase_tap_position)
      end
    end
    status = isnothing(ctrl) ? "-" : ctrl_status(ctrl)
    # partial state marker in the status column: open@to / open@from, with
    # the fictitious open-terminal voltage when the solve produced one
    if terminal_state == :open_to || terminal_state == :open_from
      marker = terminal_state == :open_to ? "open@to" : "open@from"
      if br.open_end_vm_pu !== nothing
        marker = string(marker, @sprintf(" %.4fpu", br.open_end_vm_pu))
      end
      status = status == "-" ? marker : string(marker, ", ", status)
    end

    #! format: off
    @printf(fr_io, "| %-25s | %-6s | %-25s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-9s | %-22s |\n", _fitColumn(bName, 25), branchKind, _fitColumn(connection, 25), pfromVal, qfromVal, ptoVal, qtoVal, pLossval, qLossval, ctrl_type, p_target, tap_pos, status)
    #! format: on
    shown_rows += 1
  end
  # HVDC links look like branches (they connect two buses and move power)
  # but are NOT Y-bus elements. They are appended as `Link` rows so the
  # topology reads in one table, marked "not a branch" per row; the sign
  # convention matches the branch rows (from side positive into the link,
  # to side negative = receiving), Q is 0 (a DC link transfers no reactive
  # power), Pv is the converter loss. Details in the HVDC Link Flows table.
  for r in _hvdc_link_flow_rows(net)
    connection = string(r.from_bus_name, " -> ", r.to_bus_name)
    # P_tgt is the controller SETPOINT when one exists, "-" otherwise; the
    # widths of the Ctrl (10) and status (22) columns are fixed by the
    # branch header, so mode and note pass through _fitColumn with them
    p_target = r.p_target_MW === missing ? "-" : @sprintf("%.3f", r.p_target_MW)
    link_mode = r.mode === :fixed ? "-" : _fitColumn(string(r.mode), 10)
    #! format: off
    @printf(fr_io, "| %-25s | %-6s | %-25s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-9s | %-22s |\n", _fitColumn(r.name, 25), "Link", _fitColumn(connection, 25), r.p_from_MW, 0.0, -r.p_to_MW, 0.0, r.loss_MW, 0.0, link_mode, _fitColumn(p_target, 10), "-", _fitColumn("HVDC, not a branch", 22))
    #! format: on
  end
  @printf(fr_io, "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  (∑pv, ∑qv) = getTotalLosses(net = net)
  total_losses = @sprintf("total network power balance (Σ S_branch): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

  if !isnothing(max_rows) && length(net.branchVec) > shown_rows
    @printf(fr_io, "Branch results shown: %d / %d\n", shown_rows, length(net.branchVec))
  end
  return String(take!(fr_io)), total_losses
end

"""
    _print_wrong_branch_summary_line(io, net)

Prints a single console/log line summarizing the wrong-branch detection
outcome when it is suspect or invalid (`status` is neither `:ok` nor
`:not_checked`). Clean runs and disabled detection print nothing, keeping
logs stable for the common case.
"""
function _print_wrong_branch_summary_line(io::IO, net::Net)
  rect_status = rectangular_pf_status(net)
  status = _rect_status_get(rect_status, :wrong_branch_status, :not_checked)
  status === :ok && return nothing
  status === :not_checked && return nothing

  reason = _rect_status_get(rect_status, :wrong_branch_reason, :unknown)
  detail = if reason === :low_voltage_magnitude
    string(_rect_status_get(rect_status, :wrong_branch_low_vm_count, 0), " low-Vm bus(es)")
  elseif reason === :high_voltage_magnitude
    string(_rect_status_get(rect_status, :wrong_branch_high_vm_count, 0), " high-Vm bus(es)")
  elseif reason === :angle_spread_exceeded
    @sprintf("angle spread %.1f°", _rect_status_get(rect_status, :wrong_branch_angle_spread_deg, NaN))
  elseif reason === :branch_angle_exceeded
    string(_rect_status_get(rect_status, :wrong_branch_branch_angle_violation_count, 0), " branch-angle violation(s)")
  else
    string(reason)
  end
  label = status === :fail ? "FAIL" : "SUSPECT"
  # same 15-character label rule as every other classical header line
  println(io, "Wrong-branch   : ", label, " (", reason, ", ", detail, ")")
  return nothing
end

"""
    _print_distributed_slack_summary_line(io, net)

One line inside the classical result header when the last solve ran with
the distributed slack active: mode, the solved `lambda_P`, and the
participant count. The per-bus participation itself lives in the bus table
(columns `dSl alpha` and `Pg eff MW`, see
[`_distributed_slack_bus_shares`](@ref)). Prints nothing when the feature
was off or the net has no rectangular solver status.
"""
function _print_distributed_slack_summary_line(io::IO, net::Net)
  rect_status = rectangular_pf_status(net)
  _rect_status_get(rect_status, :distributed_slack_active, false) || return nothing
  mode = _rect_status_get(rect_status, :distributed_slack_mode, :none)
  lambda_mw = _rect_status_get(rect_status, :distributed_slack_lambda_p_mw, 0.0)
  rows = _rect_status_get(rect_status, :distributed_slack_participation, nothing)
  # label padded to 15 characters: the classical result header aligns every
  # colon at column 16
  @printf(io, "Dist. slack    : mode %s, lambda_P = %+.3f MW (imbalance + losses picked up by %d participant(s), see the dSl alpha column)\n", String(mode), lambda_mw, rows === nothing ? 0 : length(rows))
  return nothing
end

"""
    _distributed_slack_bus_shares(net) -> (active::Bool, shares::Dict{Int,NTuple{2,Float64}})

Per-bus distributed-slack participation for the bus table of the classical
result print: bus index to `(alpha, dp_mw)`, aggregated over the
generators of a bus (the persisted participation table is per generator).
`active` is false when the last solve ran without the distributed slack or
no solver status exists; the table then omits the participation columns.
"""
function _distributed_slack_bus_shares(net::Net)
  rect_status = rectangular_pf_status(net)
  active = _rect_status_get(rect_status, :distributed_slack_active, false)
  shares = Dict{Int,NTuple{2,Float64}}()
  active || return (false, shares)
  rows = _rect_status_get(rect_status, :distributed_slack_participation, nothing)
  rows === nothing && return (true, shares)
  for r in rows
    hasproperty(r, :bus_idx) || continue
    a, d = get(shares, r.bus_idx, (0.0, 0.0))
    shares[r.bus_idx] = (a + r.alpha, d + r.dp_mw)
  end
  return (true, shares)
end

function printACPFlowResults(
  net::Net,
  ct::Float64,
  ite::Int,
  tol::Float64,
  toFile::Bool = false,
  path::String = "";
  converged::Bool = true,
  solver::Symbol = :NR,
  solver_time_s::Union{Nothing,Float64} = nothing,
  result_mode::Symbol = :classic,
  max_rows::Union{Nothing,Int} = nothing,
)
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
  totalLosses = let (∑pv, ∑qv) = getTotalLosses(net = net)
    @sprintf("total network power balance (Σ S_branch): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)
  end
  # converter losses are NOT part of getTotalLosses (no Y-bus branch exists
  # for an HVDC link); the extra line below the balance keeps them visible
  if !isempty(net.hvdcLinks)
    hvdc_loss = sum(r.loss_MW for r in _hvdc_link_flow_rows(net) if r.status == 1; init = 0.0)
    totalLosses *= @sprintf("converter losses (HVDC): P = %10.3f [MW]\n", hvdc_loss)
  end

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
  control_labels = _bus_control_flag_cache(net)

  nodes = sort(net.nodeVec, by = x -> x.busIdx)

  npv = 0
  npq = 0
  niso = 0
  # multi-island runs carry one reference per island, so the slack count is
  # real data, not always 1. References modeled as external-grid sources are
  # counted separately (nsource), matching the SOURCE label in the bus table:
  # total references = nslack + nsource.
  nslack = 0
  nsource = 0
  for n in nodes
    npv += n._nodeType == Sparlectra.PV ? 1 : 0
    npq += isPQNode(n) ? 1 : 0
    niso += isIsolated(n) ? 1 : 0
    if n._nodeType == Sparlectra.Slack
      # same name source as the bus-table label (_bus_type_label), so the
      # header count and the SOURCE rows can never disagree
      if _is_external_grid_internal_bus(_effective_bus_name(busNameByIdx, net, n.busIdx))
        nsource += 1
      else
        nslack += 1
      end
    end
    if occursin("_Aux_", n.comp.cName)
      auxb += 1
    end
  end
  for ps in net.prosumpsVec
    loads += ps.proSumptionType == Sparlectra.Consumption ? 1 : 0
    gens += ps.proSumptionType == Sparlectra.Injection ? 1 : 0
  end
  shunts = length(net.shuntVec)
  tap_ctrl_count, qu_ctrl_count, pu_ctrl_count = _controller_counts(net)
  series_ctrl_count = count(c -> c.enabled, _series_reactance_controllers(net))
  upfc_ctrl_count = count(c -> c.enabled, _upfc_full_controllers(net))
  # the outer-loop FACTS controllers (machine remote-voltage/STATCOM, SVC
  # shunt, HVDC pair) used to be printed in the Control footer but never
  # counted here, so a net with e.g. a UPFC quadrature composite (SSSC +
  # STATCOM) or a lone STATCOM reported "Controllers: 0". Count them too.
  machine_ctrl_count = count(c -> c.enabled, _machine_controllers(net))
  shunt_ctrl_count = count(c -> c.enabled, _shunt_controllers(net))
  hvdc_ctrl_count = count(c -> c.enabled, _hvdc_pair_controllers(net))
  total_ctrl_count = tap_ctrl_count + qu_ctrl_count + pu_ctrl_count + series_ctrl_count + upfc_ctrl_count + machine_ctrl_count + shunt_ctrl_count + hvdc_ctrl_count

  @printf(io, "Date           :%20s\n", current_date)
  @printf(io, "Iterations     :%10d\n", ite)
  @printf(io, "Flatstart      :%10s\n", net.flatstart ? "Yes" : "No")
  @printf(io, "Tolerance      : %.1e\n", tol)
  @printf(io, "Solver         :%15s\n", string(solver))
  if converged
    isnothing(solver_time_s) || println(io, "Solver time    : ", @sprintf("%.6f s", solver_time_s))
    println(io, "Total time     : ", @sprintf("%.6f s", ct))
  else
    @printf(io, "Status         :%10s\n", "Not Converged")
  end
  @printf(io, "Case           :%15s\n", net.name)
  @printf(io, "Cooldown iters :%10d\n", net.cooldown_iters)
  @printf(io, "Q-hysteresis   :%10.4f pu\n", net.q_hyst_pu)
  # Always on: prefers the lazy estimate the solver stored for the exact
  # system it factored (post-merge topology, final active set); falls back
  # to a standalone reconstruction for non-NR runs. A failed estimate never
  # breaks the report, but the reason must be visible (see
  # _jacobian_condest, condition_number.jl).
  kappa = _jacobian_condest(net; context = "result output")
  kappa === nothing || @printf(io, "Jacobian cond. : %s\n", _condition_report_line(kappa; tol = tol))

  @printf(io, "BaseMVA        :%10d\n", net.baseMVA)
  # sources appear in the count only when present, keeping the common
  # no-source header unchanged
  slack_part = nsource > 0 ? @sprintf("Slack: %d Source: %d", nslack, nsource) : @sprintf("Slack: %d", nslack)
  if auxb > 0 && niso > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Iso: %d %s\n", busses, npv, npq, auxb, niso, slack_part)
  elseif auxb > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) %s\n", busses, npv, npq, auxb, slack_part)
  else
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d %s)\n", busses, npv, npq, slack_part)
  end
  @printf(io, "Grid connection: %s\n", _grid_connection_summary(net, busNameByIdx))
  if pf_nodes != busses
    @printf(io, "PF Nodes       :%10d (after active-link merge)\n", pf_nodes)
  end
  @printf(io, "Branches       :%10d\n", branches)
  # partially open branches are IN service (they draw their charging) and
  # are counted separately from out-of-service ones; line only when present
  n_partial = count(br -> _branch_terminal_state(br) in (:open_from, :open_to), net.branchVec)
  n_partial > 0 && @printf(io, "Open terminals :%10d (branches open at one terminal)\n", n_partial)
  # one status line per one-sided open branch: the fictitious open-terminal
  # voltage (Ferranti rise) with its overhang relative to the feeding bus.
  # Printed only when such branches exist, so fully closed nets keep the
  # byte-stable summary.
  for br in net.branchVec
    st = _branch_terminal_state(br)
    st in (:open_from, :open_to) || continue
    br.open_end_vm_pu === nothing && continue
    open_bus = st == :open_to ? Int(br.toBus) : Int(br.fromBus)
    closed_bus = st == :open_to ? Int(br.fromBus) : Int(br.toBus)
    open_name = get(busNameByIdx, open_bus, string(open_bus))
    closed_name = get(busNameByIdx, closed_bus, string(closed_bus))
    closed_vm = net.nodeVec[closed_bus]._vm_pu
    if closed_vm !== nothing && isfinite(closed_vm) && closed_vm > 0.0
      rise_pct = (br.open_end_vm_pu / closed_vm - 1.0) * 100.0
      @printf(io, "  open end     : %s at %s: Vm %.4f pu, Va %.3f deg (%+.2f %% vs feeding bus %s at %.4f pu)\n", getCompName(br.comp), open_name, br.open_end_vm_pu, something(br.open_end_va_deg, NaN), rise_pct, closed_name, closed_vm)
    else
      @printf(io, "  open end     : %s at %s: Vm %.4f pu, Va %.3f deg\n", getCompName(br.comp), open_name, br.open_end_vm_pu, something(br.open_end_va_deg, NaN))
    end
  end
  @printf(io, "Links          :%10d\n", links)
  # always printed, including 0: stable parser anchor (same rationale as
  # "Transformer controls: none")
  @printf(io, "HVDC links     :%10d\n", length(net.hvdcLinks))
  @printf(io, "Lines          :%10d\n", lines)
  @printf(io, "Trafos         :%10d\n", trafos)
  @printf(io, "Generators     :%10d\n", gens)
  @printf(io, "Loads          :%10d\n", loads)
  @printf(io, "Shunts         :%10d\n", shunts)
  # the FACTS sub-counts appear only when present, so nets without such a
  # controller keep the byte-stable line for tools and parsers
  ctrl_detail = @sprintf("Tap: %d, Q(U): %d, P(U): %d", tap_ctrl_count, qu_ctrl_count, pu_ctrl_count)
  machine_ctrl_count > 0 && (ctrl_detail *= @sprintf(", MachV: %d", machine_ctrl_count))
  shunt_ctrl_count > 0 && (ctrl_detail *= @sprintf(", SVC: %d", shunt_ctrl_count))
  series_ctrl_count > 0 && (ctrl_detail *= @sprintf(", TCSC: %d", series_ctrl_count))
  hvdc_ctrl_count > 0 && (ctrl_detail *= @sprintf(", HVDC: %d", hvdc_ctrl_count))
  upfc_ctrl_count > 0 && (ctrl_detail *= @sprintf(", UPFC: %d", upfc_ctrl_count))
  @printf(io, "Controllers    :%10d (%s)\n", total_ctrl_count, ctrl_detail)

  num_guarded_locks = length(net.qLimitEvents)
  num_iterative_events = length(net.qLimitLog)
  @printf(io, "PV→PQ locks    :%10d\n", num_guarded_locks)
  @printf(io, "PV→PQ events   :%10d\n", num_iterative_events)

  _print_wrong_branch_summary_line(io, net)
  _print_distributed_slack_summary_line(io, net)

  println(io, "\n", totalLosses)
  if result_mode === :summary
    if toFile
      close(io)
    end
    return nothing
  end

  flowResults, _ = formatBranchResults(net; max_rows = max_rows)
  tap_target_vm = _tap_voltage_target_by_bus(net)
  # distributed-slack participation columns appear only when the last solve
  # ran with the feature active; plain runs keep the byte-stable table
  ds_active, ds_shares = _distributed_slack_bus_shares(net)
  ds_hdr = ds_active ? @sprintf(" %-10s | %-10s |", "dSl alpha", "Pg eff MW") : ""
  ds_pad_eq = ds_active ? repeat("=", 26) : ""
  ds_pad_dash = ds_active ? repeat("-", 26) : ""
  @printf(io, "==========================================================================================================================================================================================================================%s\n", ds_pad_eq)
  @printf(
    io,
    "| %-5s | %-20s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-12s | %-12s |%s\n",
    "Nr",
    "Bus",
    "Vn [kV]",
    "V [kV]",
    "V [pu]",
    "phi [deg]",
    "Pg [MW]",
    "Qg [MVar]",
    "Pl [MW]",
    "Ql [MVar]",
    "Ps [MW]",
    "Qs [MVar]",
    "Type",
    "Control",
    "Tap Vm tgt",
    ds_hdr
  )
  @printf(io, "==========================================================================================================================================================================================================================%s\n", ds_pad_eq)

  # buses sitting at the OPEN terminal of a one-sided open branch get an
  # "open-end" marker in the Control column. For an ISOLATED such bus the
  # solved state carries no voltage of its own (dead busbar behind the open
  # breaker), so the row shows the branch's open-end (Ferranti) voltage
  # instead of the meaningless start value; the substitution needs exactly
  # one open branch end at the bus to be unambiguous, and an energized bus
  # (fed from elsewhere) always keeps its real solved voltage.
  open_end_buses = Set{Int}()
  open_end_v = Dict{Int,Vector{Tuple{Float64,Float64}}}()
  for br in net.branchVec
    st = _branch_terminal_state(br)
    bus = st == :open_to ? Int(br.toBus) : (st == :open_from ? Int(br.fromBus) : 0)
    bus == 0 && continue
    push!(open_end_buses, bus)
    br.open_end_vm_pu === nothing && continue
    push!(get!(open_end_v, bus, Tuple{Float64,Float64}[]), (br.open_end_vm_pu, something(br.open_end_va_deg, NaN)))
  end

  # buses whose voltage an outer-loop FACTS controller holds get a device
  # marker in the Control column (appended at render time, so the Q(U)/P(U)
  # flag cache and its tests are untouched): the machine controller regulates
  # its target bus (STATCOM in current-limit mode, else RVC), the shunt
  # controller its own bus (SVC continuous, MSC discrete bank).
  facts_bus_labels = Dict{Int,String}()
  for c in _machine_controllers(net)
    c.enabled || continue
    facts_bus_labels[geNetBusIdx(net = net, busName = c.target_bus)] = c.limit_mode === :current ? "STATCOM" : "RVC"
  end
  for c in _shunt_controllers(net)
    c.enabled || continue
    facts_bus_labels[geNetBusIdx(net = net, busName = c.bus)] = c.step_mvar === nothing ? "SVC" : "MSC"
  end

  pGS = qGS = pLS = qLS = ""
  tpGS = tqGS = tpLS = tqLS = 0.0
  pShunt_str = qShunt_str = ""
  tpShunt = tqShunt = 0.0
  shown_bus_rows = 0
  for n in nodes
    if !isnothing(max_rows) && shown_bus_rows >= max_rows
      break
    end
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
    nodeName = get(busNameByIdx, n.busIdx, n.comp.cName)
    typeStr = _bus_type_label(n._nodeType, nodeName)
    controlStr = _cached_control_label(control_labels, n.busIdx)
    if n.busIdx in open_end_buses
      controlStr = (isempty(controlStr) || controlStr == "-") ? "open-end" : string(controlStr, ",open-end")
    end
    facts_label = get(facts_bus_labels, n.busIdx, nothing)
    if facts_label !== nothing
      controlStr = (isempty(controlStr) || controlStr == "-") ? facts_label : string(controlStr, ",", facts_label)
    end

    # Mark PV→PQ buses (hit Q-limit) with a star in the Type column
    if haskey(net.qLimitEvents, n.busIdx)
      typeStr *= "*"
    end

    # isolated open-end bus: show the open-end (Ferranti) voltage in the
    # V/phi columns (see the open_end_v comment above). Isolation is read
    # from net.isoNodes: the solver resets the node TYPE of excluded buses
    # to PQ with a 1.0-pu start value before solving, so isIsolated(n) is
    # false again at print time.
    vm_show = n._vm_pu
    va_show = n._va_deg
    if (n.busIdx in net.isoNodes || isIsolated(n)) && length(get(open_end_v, n.busIdx, Tuple{Float64,Float64}[])) == 1
      vm_show, va_show = only(open_end_v[n.busIdx])
    end
    v = n.comp.cVN * vm_show
    if !isnothing(n._vmin_pu) && !isnothing(n._vmax_pu)
      if !isIsolated(n) && (n._vm_pu < n._vmin_pu || n._vm_pu > n._vmax_pu)
        nodeName *= " !"
      end
    end

    tap_target_str = haskey(tap_target_vm, n.busIdx) ? @sprintf("%.4f", tap_target_vm[n.busIdx]) : ""
    ds_cells = ""
    if ds_active
      share = get(ds_shares, n.busIdx, nothing)
      # Pg eff = scheduled bus generation plus the alpha share of lambda_P;
      # the Pg column keeps the schedule, this column shows what the
      # participant actually delivers
      ds_cells = share === nothing ? @sprintf(" %-10s | %-10s |", "", "") : @sprintf(" %-10s | %-10s |", @sprintf("%.4f", share[1]), @sprintf("%.3f", p_gen + share[2]))
    end
    @printf(
      io,
      "| %-5d | %-20s | %-10.1f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-12s | %-12s |%s\n",
      n.busIdx,
      _fitColumn(nodeName, 20),
      n.comp.cVN,
      v,
      vm_show,
      va_show,
      pGS,
      qGS,
      pLS,
      qLS,
      pShunt_str,
      qShunt_str,
      typeStr,
      controlStr,
      tap_target_str,
      ds_cells
    )
    shown_bus_rows += 1
  end
  if !isnothing(max_rows) && length(nodes) > shown_bus_rows
    @printf(io, "Bus results shown: %d / %d\n", shown_bus_rows, length(nodes))
  end

  @printf(io, "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%s\n", ds_pad_dash)
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
      @printf(io, "| %-5d | %-8s | %-8s | %-6d | %-12.3f | %-12.3f | %-12.4f | %-12.4f |\n", l.linkIdx, _fitColumn(fromName, 8), _fitColumn(toName, 8), l.status, p, q, ifrom, ito)
    end
  end
  if !isempty(net.hvdcLinks)
    # P_from = power leaving the from bus into the link (positive for
    # export), P_to = power delivered into the to bus; Status is the
    # controller status, "-" for Stage-0 fixed injections. Header,
    # separator, and rows are generated from _HVDC_TABLE_COLUMNS, whose
    # widths derive from the canonical mode/status sets (Phase 6).
    header = _hvdc_table_header()
    title = " HVDC Link Flows "
    pad = max(0, length(header) - length(title))
    println(io, "\n", "-"^(pad ÷ 2), title, "-"^(pad - pad ÷ 2))
    println(io, header)
    println(io, "-"^length(header))
    for r in _hvdc_link_flow_rows(net)
      rating = r.p_rating_MW === missing ? "-" : @sprintf("%.1f", r.p_rating_MW)
      println(io, _hvdc_table_line((r.nr, r.name, r.from_bus_name, r.to_bus_name, r.mode, @sprintf("%.3f", r.p_from_MW), @sprintf("%.3f", r.p_to_MW), @sprintf("%.3f", r.loss_MW), @sprintf("%.3f", r.q_from_MVar), @sprintf("%.3f", r.q_to_MVar), rating, r.ctrl_status)))
    end
  end
  println(io, "\nControl")
  if !isnothing(max_rows) && (length(nodes) > shown_bus_rows || length(net.branchVec) > max_rows)
    println(io, "Large-case result output is row-limited: $(shown_bus_rows) bus rows and $(min(max_rows, length(net.branchVec))) branch rows shown.")
    println(io, "For timing-focused runs use output.logfile_results=summary.")
  end
  println(io, "-------")
  if isempty(_tap_controllers(net))
    # Explicit statement keeps report output deterministic for tools/parsers.
    println(io, "Transformer controls: none")
  else
    # Detailed engineering-style control section (OLTC/PST/combined).
    printTapControllerSummary(io, net)
  end
  # Machine remote voltage controllers print only when present — the "none"
  # marker above stays the deterministic anchor for parsers.
  printMachineControllerSummary(io, net)
  # Shunt voltage controllers (SVC / MSC-MSR bank) print only when present.
  printShuntVoltageControllerSummary(io, net)
  # Series-reactance (TCSC) controllers print only when present as well.
  printSeriesReactanceControllerSummary(io, net)
  # HVDC pair controllers print only when present as well.
  printHvdcPairControllerSummary(io, net)
  # full UPFC controllers (#326) print only when present as well.
  printUpfcFullControllerSummary(io, net)
  if toFile
    close(io)
    #println("Results have been written to $(joinpath(path, filename))")
  end
end

# prosumer result accessors shared by the generator table and the HVDC rows:
# prefer solved results, fall back to the specification values
_prosumer_p_result(ps) = something(ps.pRes === nothing ? ps.pVal : ps.pRes, 0.0)
_prosumer_q_result(ps) = something(ps.qRes === nothing ? ps.qVal : ps.qRes, 0.0)

"""
    _hvdc_link_flow_rows(net) -> Vector{NamedTuple}

One row per `HvdcLink` for the `HVDC Link Flows` table, `ACPFlowReport`, and
the CSV export. Sign convention: `p_from_MW` is the power leaving the from
bus into the link (positive for export), `p_to_MW` the power delivered into
the to bus. With an attached controller the values come from its setpoints
and live terminal state (`mode` is the controller mode); without one they
come from the terminal prosumers (`mode = :fixed`,
`loss = -(P_from_injection + P_to_injection)` in the MATPOWER convention
where the from injection is negative).
"""
function _hvdc_link_flow_rows(net::Net)::Vector{NamedTuple}
  ctrl_by_name = Dict{String,Any}(c.name => c for c in _hvdc_pair_controllers(net))
  rows = NamedTuple[]
  busNameByIdx = _bus_name_by_idx(net)
  for (nr, l) in enumerate(net.hvdcLinks)
    from_name = _effective_bus_name(busNameByIdx, net, l.from_bus)
    to_name = _effective_bus_name(busNameByIdx, net, l.to_bus)
    ctrl = l.controller_name === nothing ? nothing : get(ctrl_by_name, l.controller_name, nothing)
    if ctrl === nothing
      fp = _prosumer_p_result(net.prosumpsVec[l.from_prosumer])
      tp = _prosumer_p_result(net.prosumpsVec[l.to_prosumer])
      push!(rows, (
        nr = nr,
        name = l.name,
        from_bus_name = from_name,
        to_bus_name = to_name,
        mode = :fixed,
        p_from_MW = -fp,
        p_to_MW = tp,
        loss_MW = -(fp + tp),
        q_from_MVar = _prosumer_q_result(net.prosumpsVec[l.from_prosumer]),
        q_to_MVar = _prosumer_q_result(net.prosumpsVec[l.to_prosumer]),
        p_rating_MW = missing,
        p_target_MW = missing,
        status = l.status,
        ctrl_status = "-",
      ))
    else
      sp = _hvdc_pair_setpoints(ctrl)
      push!(rows, (
        nr = nr,
        name = l.name,
        from_bus_name = from_name,
        to_bus_name = to_name,
        mode = ctrl.mode,
        p_from_MW = sp.transfer,
        p_to_MW = sp.to_p,
        loss_MW = sp.loss,
        q_from_MVar = ctrl.from_q_now,
        q_to_MVar = ctrl.to_q_now,
        p_rating_MW = ctrl.p_rating_mw === nothing ? missing : ctrl.p_rating_mw,
        # the controller SETPOINT (pre-clamp target), for the P_tgt column
        p_target_MW = ctrl.p_transfer_mw,
        status = l.status,
        ctrl_status = string(ctrl.status),
      ))
    end
  end
  return rows
end

# generator-table terminal marker per prosumer index: "B2B from"/"B2B to"
# for back-to-back links, "HVDC from"/"HVDC to" for point-to-point
function _hvdc_terminal_labels(net::Net)::Dict{Int,String}
  labels = Dict{Int,String}()
  for l in net.hvdcLinks
    prefix = l.kind == :p2p ? "HVDC" : "B2B"
    labels[l.from_prosumer] = string(prefix, " from")
    labels[l.to_prosumer] = string(prefix, " to")
  end
  return labels
end

function formatProsumerResults(net::Net)
  buf = IOBuffer()

  # Rebuild mapping: bus index -> bus name
  busNameByIdx = _bus_name_by_idx(net)
  hvdc_labels = _hvdc_terminal_labels(net)

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

      # HVDC terminals are converters, not machines: the marker replaces the
      # Q-limit status text (P/Q columns stay untouched)
      status = haskey(hvdc_labels, idx) ? hvdc_labels[idx] : status_from_Q(ps, q)

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
