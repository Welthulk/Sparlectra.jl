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
# file: src/acpflow/start_modes.jl

function _apply_matpower_start_modes!(net::Net, mpc, start_cfg::StartModeConfig, mat_cfg::MatpowerImportConfig; performance_profile = nothing)
  vmode = start_cfg.voltage_mode
  amode = start_cfg.angle_mode
  psource = start_cfg.profile_source
  vmode in (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :profile_blend) || error("power_flow.start_mode.voltage_mode must be classic, pv_gen_vg, pv_bus_vm, all_bus_vm, or profile_blend")
  amode in (:classic, :dc, :bus_va_blend, :matpower_va) || error("power_flow.start_mode.angle_mode must be classic, dc, bus_va_blend, or matpower_va")
  busrow = _perf_profile_time!(performance_profile, :start_projection_matpower_bus_map) do
    MatpowerIO.bus_row_index(mpc)
  end
  _perf_profile_time!(performance_profile, :start_projection_matpower_branch_map) do
    count(e -> size(mpc.branch, 2) < 11 || mpc.branch[e, 11] != 0.0, axes(mpc.branch, 1))
  end
  pv_rows = _perf_profile_time!(performance_profile, :start_projection_matpower_reference_lookup) do
    MatpowerIO.pv_voltage_reference_rows(mpc; matpower_pv_voltage_source = mat_cfg.pv_voltage_source, tol = mat_cfg.pv_voltage_mismatch_tol_pu, warn = false)
  end
  imported_vset = Dict(row.busI => row.imported_vset for row in pv_rows)
  gen_vset = Dict(row.busI => (isempty(row.gen_vgs) ? row.bus_vm : row.gen_vgs[1]) for row in pv_rows)

  for k in eachindex(net.nodeVec)
    node = net.nodeVec[k]
    busI = get(net.busOrigIdxDict, k, k)
    haskey(busrow, busI) || continue
    r = busrow[busI]
    btype = Int(mpc.bus[r, 2])
    bus_vm = Float64(mpc.bus[r, 8])
    bus_va = Float64(mpc.bus[r, 9])
    if vmode == :pv_gen_vg && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = get(gen_vset, busI, bus_vm))
    elseif vmode == :pv_bus_vm && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :all_bus_vm
      setVmVa!(node = node, vm_pu = bus_vm)
    elseif vmode == :profile_blend
      psource == :matpower_reference || error("profile_blend currently requires profile_source=:matpower_reference for MATPOWER imports.")
      setVmVa!(node = node, vm_pu = 0.5 * (something(node._vm_pu, 1.0) + bus_vm))
    elseif vmode == :classic && (btype == 2 || btype == 3)
      setVmVa!(node = node, vm_pu = get(imported_vset, busI, bus_vm))
    end
    if amode == :matpower_va
      setVmVa!(node = node, vm_pu = something(node._vm_pu, bus_vm), va_deg = bus_va)
    elseif amode == :bus_va_blend
      setVmVa!(node = node, vm_pu = something(node._vm_pu, bus_vm), va_deg = 0.5 * (something(node._va_deg, 0.0) + bus_va))
    end
  end
  return nothing
end
