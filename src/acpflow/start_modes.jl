# Copyright 2023-2026 Udo Schmitz
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
# file: src/acpflow/start_modes.jl
# purpose: applies the configured start modes to the imported net ahead of
#          the rectangular solve: voltage_mode and angle_mode handling,
#          profile blending, and DC angle starts. Format independent since
#          stage 3a (D10): the raw source start values come from the typed
#          case's start_state, not from a MATPOWER container.

"""
    _apply_start_modes!(net, case, start_cfg; performance_profile) -> Nothing

Apply the configured start voltage and angle modes with the RAW source
start values of `case.sparlectra.start_state` as the reference profile
(design decision D10). PV and slack classification and the resolved
voltage setpoints come from the built net itself: the adapter resolved
the source system's setpoint choice at conversion time, so the historical
distinction between the raw generator setpoint and the imported setpoint
collapses onto the value the machines actually regulate to.
"""
# task_import_direct: the projection core is import-path neutral. The raw
# start voltages arrive index-aligned with net.nodeVec (nothing = no raw
# value for that node); each import path builds that vector from ITS
# source of truth, the SCF start_state or the parsed MATPOWER bus rows.
function _apply_start_modes!(net::Net, case::SCFCase, start_cfg::StartModeConfig; performance_profile = nothing)
  spar = case.sparlectra
  start_nodes = spar === nothing ? Dict{Int,SCFStartNode}() : spar.start_state.nodes
  ids = _perf_profile_time!(performance_profile, :start_projection_source_bus_map) do
    scf_node_ids_in_build_order(case)
  end
  raw_starts = Union{Nothing,NamedTuple}[k <= length(ids) ? (raw = get(start_nodes, ids[k], nothing); raw === nothing ? nothing : (vm_pu = raw.vm_pu, va_deg = raw.va_deg)) : nothing for k in eachindex(net.nodeVec)]
  return _apply_start_modes!(net, raw_starts, start_cfg; performance_profile = performance_profile)
end

"""
    _matpower_raw_starts(mpc) -> Vector

The raw per-bus start voltages of a parsed MATPOWER case, index-aligned
with the direct import's node order (which is the bus row order).
"""
function _matpower_raw_starts(mpc)
  busDict, _, _ = _createDict()
  VM = busDict["Vm"]
  VA = busDict["Va"]
  busData = convert(Matrix{Float64}, getproperty(mpc, :bus))
  return Union{Nothing,NamedTuple}[(vm_pu = row[VM], va_deg = row[VA]) for row in eachrow(busData)]
end

function _apply_start_modes!(net::Net, raw_starts::AbstractVector, start_cfg::StartModeConfig; performance_profile = nothing)
  vmode = start_cfg.voltage_mode
  amode = start_cfg.angle_mode
  psource = start_cfg.profile_source
  vmode in (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :profile_blend) || error("power_flow.start_mode.voltage_mode must be classic, pv_gen_vg, pv_bus_vm, all_bus_vm, or profile_blend")
  amode in (:classic, :dc, :bus_va_blend, :matpower_va) || error("power_flow.start_mode.angle_mode must be classic, dc, bus_va_blend, or matpower_va")

  for k in eachindex(net.nodeVec)
    node = net.nodeVec[k]
    k <= length(raw_starts) || continue
    raw = raw_starts[k]
    raw === nothing && continue
    bus_vm = raw.vm_pu === nothing ? something(node._vm_pu, 1.0) : raw.vm_pu
    bus_va = raw.va_deg === nothing ? something(node._va_deg, 0.0) : raw.va_deg
    ntype = getNodeType(node)
    is_setpoint_bus = ntype == PV || ntype == Slack
    # the resolved setpoint the build laid down (regulator u_ref); for a
    # setpoint bus without one the raw start value stands in, exactly like
    # the historical BUS.VM fallback of gen-less PV buses
    setpoint = something(node._vm_pu, bus_vm)
    if vmode == :pv_gen_vg && is_setpoint_bus
      setVmVa!(node = node, vm_pu = setpoint)
    elseif vmode == :pv_bus_vm && is_setpoint_bus
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :all_bus_vm
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :profile_blend
      if psource == :state_estimation || psource == :se_snapshot
        # SE chain start (phase 5): the start voltages come from the
        # estimated state that runpf_from_se!/the SE chain runner already
        # wrote into the net; blending them with the source reference
        # would corrupt the estimate, so the blend is a deliberate no-op.
      else
        psource == :matpower_reference || error("profile_blend currently requires profile_source=:matpower_reference for imported cases (or the SE chain sources :state_estimation/:se_snapshot, which keep the estimated voltages untouched).")
        setVmVa!(node = node, vm_pu = 0.5 * (something(node._vm_pu, 1.0) + bus_vm))
      end
    elseif vmode == :classic && is_setpoint_bus
      setVmVa!(node = node, vm_pu = setpoint)
    end
    if amode == :matpower_va
      setVmVa!(node = node, vm_pu = something(node._vm_pu, bus_vm), va_deg = bus_va)
    elseif amode == :bus_va_blend
      setVmVa!(node = node, vm_pu = something(node._vm_pu, bus_vm), va_deg = 0.5 * (something(node._va_deg, 0.0) + bus_va))
    end
  end
  return nothing
end
