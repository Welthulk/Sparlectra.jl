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

# file: src/api/run_bus_powers_export.jl
# purpose: bus_powers.csv (issue #376): one row per bus with the solved
#          generation/load/shunt power and the bus's Q-limit/control state,
#          built from the same data as the console Q-V characteristic check
#          and the controller summary so all three agree on one run.

const _BUS_POWERS_CSV_COLUMNS = (
  :bus_index, :bus_id, :bus_name, :vn_kV,
  :bus_type_start, :bus_type_end,
  :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar, :p_shunt_MW, :q_shunt_MVar,
  :q_limit_side, :q_min_MVar, :q_max_MVar,
  :control, :control_target_bus, :control_status,
  :vm_setpoint_pu, :non_physical,
)

# `control_status` vocabulary is the outer-loop three-way state
# (converged/at_limit/not_converged) issue #376 asks for. Every outer
# controller struct (MachineVoltageControl, ShuntVoltageControl,
# PowerTransformerControl, ...) carries plain `converged`/`at_limit` fields
# for exactly this; "not_converged" covers both an outer loop that ran out of
# iterations and one that is silently still active, since the API only
# reports it once the run has finished either way (see
# `_controller_status_label` in controller/tap_control.jl for the closely
# related, but five-valued, print-table version of this mapping). A parked
# discrete shunt bank (`ShuntVoltageControl.parked`, issue #324) blocks like
# `at_limit` and is folded into it here.
_bus_powers_control_status(ctrl)::String = ctrl.converged ? "converged" : (ctrl.at_limit || (hasproperty(ctrl, :parked) && ctrl.parked)) ? "at_limit" : "not_converged"

"""
    _bus_powers_rows(net::Net) -> Vector{NamedTuple}

Build the per-bus rows for `bus_powers.csv`. Reuses the same data sources as
the console result table and the Q-V characteristic check so the values agree
on one run (issue #376):

- generation/load: `_bus_power_component_cache` (same cache the detailed bus
  CSV and the console table read);
- shunt: the node's own `_pShunt`/`_qShunt` (not carried by
  `BusPowerComponents`, read directly as the console table does);
- `bus_type_start`/`bus_type_end`/`q_limit_side`: a bus can only ever move
  PV -> PQ under Q-limit enforcement (Slack and originally-PQ/Isolated buses
  never change type), so `bus_type_start` is `PV` exactly for the buses in
  the pre-solve snapshot (`net.qLimitInitialPVRows`, the same snapshot
  `q_limit_initial_limits.csv` uses) and the CURRENT type otherwise;
  `bus_type_end` is the current type with a `*` suffix while the bus sits in
  `net.qLimitEvents`, matching the star the console table prints;
- `non_physical`: a bus is flagged exactly when [`qvCharacteristicViolations`](@ref)
  reports it (clamped at a limit with the voltage on the wrong side of its
  setpoint) - the same rows `printQVCharacteristicCheck` prints;
- `control`/`control_target_bus`/`control_status`/`vm_setpoint_pu`: combines
  the Q(U)/P(U)/HVDC flags from `_bus_control_flag_cache` with the
  RVC/STATCOM (machine) and SVC/MSC (shunt) outer controllers and the
  voltage-mode tap-changer target (`_tap_voltage_target_by_bus`),
  i.e. exactly the labels the console result table's Control column shows
  (`results.jl`, the `facts_bus_labels`/`_cached_control_label` block) plus
  `OLTC_target` for a tap-controller's regulated bus, which the console
  table does not currently label at all. `control_target_bus` is populated
  for the machine's own bus row of an RVC/STATCOM controller (the bus it
  regulates); for the `OLTC_target` row itself there is no further "other"
  bus to name, so it stays empty there - a deliberate asymmetry, not a gap.
"""
function _bus_powers_rows(net::Net)::Vector{NamedTuple}
  busNameByIdx = _bus_name_by_idx(net)
  power_components = _bus_power_component_cache(net)
  control_flags = _bus_control_flag_cache(net)
  vset = _bus_voltage_setpoints_from_prosumers(net)
  qmin_pu, qmax_pu = getQLimits_pu(net)
  qv_rows = Dict(r.bus => r for r in qvCharacteristicViolations(net))
  initial_pv_rows = isempty(net.qLimitInitialPVRows) ? snapshotPVQLimits!(net) : net.qLimitInitialPVRows
  initial_pv_buses = Set(row.bus for row in initial_pv_rows)

  # RVC/STATCOM label sits on the MACHINE's own bus (the actuator), per issue
  # text ("remote voltage controller at this machine"); its control_target_bus
  # is the far bus it regulates.
  machine_by_own_bus = Dict{Int,MachineVoltageControl}()
  for c in _machine_controllers(net)
    c.enabled || continue
    machine_by_own_bus[geNetBusIdx(net = net, busName = c.bus)] = c
  end
  shunt_by_own_bus = Dict{Int,ShuntVoltageControl}()
  for c in _shunt_controllers(net)
    c.enabled || continue
    shunt_by_own_bus[geNetBusIdx(net = net, busName = c.bus)] = c
  end
  oltc_target_vm = _tap_voltage_target_by_bus(net)
  oltc_ctrl_by_target_bus = Dict{Int,PowerTransformerControl}()
  for c in _tap_controllers(net)
    c.enabled || continue
    c.mode in (:voltage, :voltage_and_branch_active_power) || continue
    isnothing(c.target_bus) && continue
    oltc_ctrl_by_target_bus[geNetBusIdx(net = net, busName = c.target_bus)] = c
  end

  rows = NamedTuple[]
  for n in sort(net.nodeVec, by = x -> x.busIdx)
    bus = n.busIdx
    p_gen, q_gen, p_load, q_load = _bus_power_components(power_components, bus)
    clamped = haskey(net.qLimitEvents, bus)
    type_end_base = toString(n._nodeType)
    bus_type_start = bus in initial_pv_buses ? "PV" : type_end_base
    bus_type_end = clamped ? string(type_end_base, "*") : type_end_base
    q_limit_side = clamped ? String(net.qLimitEvents[bus]) : ""
    q_min_MVar = (bus <= length(qmin_pu) && isfinite(qmin_pu[bus])) ? qmin_pu[bus] * net.baseMVA : ""
    q_max_MVar = (bus <= length(qmax_pu) && isfinite(qmax_pu[bus])) ? qmax_pu[bus] * net.baseMVA : ""

    flags = get(control_flags, bus, _NO_BUS_CONTROL_FLAGS)
    control_parts = String[]
    flags.has_qu && push!(control_parts, "Q(U)")
    flags.has_pu && push!(control_parts, "P(U)")
    flags.hvdc == :terminal && push!(control_parts, "B2B")
    flags.hvdc == :gridforming && push!(control_parts, "B2B src")
    control_target_bus = ""
    control_status = ""
    vm_setpoint_pu = ""
    if haskey(machine_by_own_bus, bus)
      c = machine_by_own_bus[bus]
      push!(control_parts, c.limit_mode === :current ? "STATCOM" : "RVC")
      control_target_bus = c.target_bus
      control_status = _bus_powers_control_status(c)
      vm_setpoint_pu = c.target_vm_pu
    end
    if haskey(shunt_by_own_bus, bus)
      c = shunt_by_own_bus[bus]
      push!(control_parts, c.step_mvar === nothing ? "SVC" : "MSC")
      control_status = _bus_powers_control_status(c)
      vm_setpoint_pu = c.target_vm_pu
    end
    if haskey(oltc_ctrl_by_target_bus, bus)
      c = oltc_ctrl_by_target_bus[bus]
      push!(control_parts, "OLTC_target")
      control_status = _bus_powers_control_status(c)
      vm_setpoint_pu = get(oltc_target_vm, bus, "")
    end
    # PV and clamped (PQ*) buses hold a voltage setpoint even without an
    # outer-loop controller (the classic PV/Q-limit switching setpoint); an
    # outer controller above already set vm_setpoint_pu, so only fill the gap.
    if vm_setpoint_pu === "" && (bus_type_start == "PV" || clamped) && bus <= length(vset) && isfinite(vset[bus])
      vm_setpoint_pu = vset[bus]
    end
    control = isempty(control_parts) ? "none" : join(control_parts, ", ")

    push!(
      rows,
      (
        bus_index = bus,
        bus_id = n.comp.cName,
        bus_name = get(busNameByIdx, bus, n.comp.cName),
        vn_kV = n.comp.cVN,
        bus_type_start = bus_type_start,
        bus_type_end = bus_type_end,
        p_gen_MW = p_gen,
        q_gen_MVar = q_gen,
        p_load_MW = p_load,
        q_load_MVar = q_load,
        p_shunt_MW = _default0(n._pShunt),
        q_shunt_MVar = _default0(n._qShunt),
        q_limit_side = q_limit_side,
        q_min_MVar = q_min_MVar,
        q_max_MVar = q_max_MVar,
        control = control,
        control_target_bus = control_target_bus,
        control_status = control_status,
        vm_setpoint_pu = vm_setpoint_pu,
        non_physical = haskey(qv_rows, bus),
      ),
    )
  end
  return rows
end

"""
    _write_bus_powers_artifact(output_dir, net; format = "technical") -> String

Write `bus_powers.csv` next to `bus_voltages_complex.csv` (issue #376). See
[`_bus_powers_rows`](@ref) for the per-bus values and how they are kept
consistent with the console Q-V characteristic check and controller summary.
`format` selects the delimiter/decimal separator like every other CSV
artifact of the run (`_resolve_detailed_csv_format`).
"""
function _write_bus_powers_artifact(output_dir::AbstractString, net::Net; format = "technical")::String
  path = joinpath(output_dir, "bus_powers.csv")
  _write_namedtuple_csv(path, _bus_powers_rows(net), _BUS_POWERS_CSV_COLUMNS; format = format)
  return "bus_powers.csv"
end
