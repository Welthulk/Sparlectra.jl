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
"""
struct ACPFlowReport
  metadata::NamedTuple
  nodes::Vector{NamedTuple}
  branches::Vector{NamedTuple}
  links::Vector{NamedTuple}
  transformer_controls::Vector{NamedTuple}
  q_limit_events::Vector{NamedTuple}
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
      terminal = name[1:(end - length("__extgrid_int"))]
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
end

const _NO_BUS_CONTROL_FLAGS = BusControlFlags(false, false)

function _bus_control_flag_cache(net::Net)::Dict{Int,BusControlFlags}
  cache = Dict{Int,BusControlFlags}()
  for ps in net.prosumpsVec
    bus_idx = getPosumerBusIndex(ps)
    old = get(cache, bus_idx, _NO_BUS_CONTROL_FLAGS)
    cache[bus_idx] = BusControlFlags(old.has_qu | has_qu_controller(ps), old.has_pu | has_pu_controller(ps))
  end
  return cache
end

@inline function _control_label(flags::BusControlFlags)::String
  if flags.has_qu && flags.has_pu
    return "Q(U), P(U)"
  elseif flags.has_qu
    return "Q(U)"
  elseif flags.has_pu
    return "P(U)"
  else
    return "-"
  end
end

@inline _cached_control_label(cache::AbstractDict{Int,BusControlFlags}, bus_idx::Int)::String = _control_label(get(cache, bus_idx, _NO_BUS_CONTROL_FLAGS))

function _control_label(net::Net, bus_idx::Int)::String
  has_qu, has_pu = _bus_control_flags(net, bus_idx)
  return _control_label(BusControlFlags(has_qu, has_pu))
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
    push!(branch_rows, (branch = br.comp.cName, branch_index = br.branchIdx, from_bus = br.fromBus, to_bus = br.toBus, status = br.status, p_from_MW = p_from, q_from_MVar = q_from, p_to_MW = p_to, q_to_MVar = q_to, p_loss_MW = p_loss, q_loss_MVar = q_loss, rated_MVA = rated, overloaded = overload, branch_name = br.comp.cName, original_branch_name = _original_branch_name(net, br), from_bus_name = from_name, to_bus_name = to_name, original_from_bus_name = _original_bus_name(busNameByIdx, net, br.fromBus), original_to_bus_name = _original_bus_name(busNameByIdx, net, br.toBus), branch_kind = _branch_kind_name(net, br)))
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
  metadata = (case = net.name, baseMVA = net.baseMVA, converged = converged, elapsed_s = ct, iterations = ite, tolerance = tol, solver = solver, timestamp = Dates.now(), total_p_loss_MW = p_loss_total, total_q_loss_MVar = q_loss_total, wrong_branch_status = wrong_branch_status, wrong_branch_reason = wrong_branch_reason)

  # Keep structured tap-controller data as a dedicated relation in the report
  # (preferred over widening the branch table with sparse controller columns).
  tap_control_rows = _build_transformer_control_rows(net)
  return ACPFlowReport(metadata, node_rows, branch_rows, link_rows, tap_control_rows, q_events)
end

function formatBranchResults(net::Net; max_rows::Union{Nothing,Int} = nothing)
  busNameByIdx = _bus_name_by_idx(net)
  ctrl_by_branch = Dict{Int,NamedTuple}()
  for row in buildTapControllerReportRows(net)
    ctrl_by_branch[row.transformer_branch_index] = row
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
  formatted_results = @sprintf("\n==========================================================================================================================================================================================================================\n")
  formatted_results *= @sprintf("| %-25s | %-6s | %-25s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-9s | %-22s |\n", "Branch", "Type", "Connection", "P [MW]", "Q [MVar]", "P [MW]", "Q [MVar]", "Pv [MW]", "Qv [MVar]", "Ctrl", "P_tgt", "TapPos", "Ctrl status")
  formatted_results *= @sprintf("==========================================================================================================================================================================================================================\n")
  #! format: on
  shown_rows = 0
  for br in net.branchVec
    if !isnothing(max_rows) && shown_rows >= max_rows
      break
    end
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

    #! format: off
    formatted_results *= @sprintf("| %-25s | %-6s | %-25s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-9s | %-22s |\n", _fitColumn(bName, 25), branchKind, _fitColumn(connection, 25), pfromVal, qfromVal, ptoVal, qtoVal, pLossval, qLossval, ctrl_type, p_target, tap_pos, status)
    #! format: on
    shown_rows += 1
  end
  formatted_results *= @sprintf("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  (∑pv, ∑qv) = getTotalLosses(net = net)
  total_losses = @sprintf("total network power balance (Σ S_branch): P = %10.3f [MW], Q = %10.3f [MVar]\n", ∑pv, ∑qv)

  if !isnothing(max_rows) && length(net.branchVec) > shown_rows
    formatted_results *= @sprintf("Branch results shown: %d / %d\n", shown_rows, length(net.branchVec))
  end
  return formatted_results, total_losses
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
  println(io, "wrong-branch check: ", label, " (", reason, ", ", detail, ")")
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
  @printf(io, "Distributed slack: mode %s, lambda_P = %+.3f MW (imbalance + losses picked up by %d participant(s), see the dSl alpha column)\n", String(mode), lambda_mw, rows === nothing ? 0 : length(rows))
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

function printACPFlowResults(net::Net, ct::Float64, ite::Int, tol::Float64, toFile::Bool = false, path::String = ""; converged::Bool = true, solver::Symbol = :NR, solver_time_s::Union{Nothing,Float64} = nothing, result_mode::Symbol = :classic, max_rows::Union{Nothing,Int} = nothing)
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
  tap_ctrl_count, qu_ctrl_count, pu_ctrl_count = _controller_counts(net)
  series_ctrl_count = count(c -> c.enabled, _series_reactance_controllers(net))
  total_ctrl_count = tap_ctrl_count + qu_ctrl_count + pu_ctrl_count + series_ctrl_count

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

  @printf(io, "BaseMVA        :%10d\n", net.baseMVA)
  if auxb > 0 && niso > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Iso: %d Slack: %d\n", busses, npv, npq, auxb, niso, 1)
  elseif auxb > 0
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d (Aux: %d) Slack: %d\n", busses, npv, npq, auxb, 1)
  else
    @printf(io, "Nodes          :%10d (PV: %d PQ: %d Slack: %d)\n", busses, npv, npq, 1)
  end
  @printf(io, "Grid connection: %s\n", _grid_connection_summary(net, busNameByIdx))
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
  # the TCSC count appears only when present, so nets without a series
  # controller keep the byte-stable line for tools and parsers
  if series_ctrl_count > 0
    @printf(io, "Controllers    :%10d (Tap: %d, Q(U): %d, P(U): %d, TCSC: %d)\n", total_ctrl_count, tap_ctrl_count, qu_ctrl_count, pu_ctrl_count, series_ctrl_count)
  else
    @printf(io, "Controllers    :%10d (Tap: %d, Q(U): %d, P(U): %d)\n", total_ctrl_count, tap_ctrl_count, qu_ctrl_count, pu_ctrl_count)
  end

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

    # Mark PV→PQ buses (hit Q-limit) with a star in the Type column
    if haskey(net.qLimitEvents, n.busIdx)
      typeStr *= "*"
    end

    v = n.comp.cVN * n._vm_pu
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
    @printf(io, "| %-5d | %-20s | %-10.1f | %-10.3f | %-10.3f | %-10.3f | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-12s | %-12s |%s\n", n.busIdx, _fitColumn(nodeName, 20), n.comp.cVN, v, n._vm_pu, n._va_deg, pGS, qGS, pLS, qLS, pShunt_str, qShunt_str, typeStr, controlStr, tap_target_str, ds_cells)
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
  # Series-reactance (TCSC) controllers print only when present as well.
  printSeriesReactanceControllerSummary(io, net)
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
