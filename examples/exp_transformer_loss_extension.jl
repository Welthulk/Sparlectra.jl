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

using Sparlectra

"""
    main(; output_dir = joinpath(@__DIR__, "_out", "transformer_loss_extension"))

Create a minimal transformer case with active no-load conductance metadata,
export it as MATPOWER, and reimport it to demonstrate Sparlectra's proprietary
transformer-loss extension round trip.
"""
function main(; output_dir = joinpath(@__DIR__, "_out", "transformer_loss_extension"))
  mkpath(output_dir)
  net = Sparlectra.Net(name = "transformer_loss_extension", baseMVA = 100.0)
  Sparlectra.addBus!(net = net, busName = "FROM", vn_kV = 110.0, oBusIdx = 1)
  Sparlectra.addBus!(net = net, busName = "TO", vn_kV = 110.0, oBusIdx = 2)
  Sparlectra.addProsumer!(net = net, busName = "FROM", type = "GENERATOR", p = 0.0, q = 0.0, referencePri = "FROM", vm_pu = 1.0)
  Sparlectra._addPIModelTrafo_by_idx!(net = net, from = 1, to = 2, r_pu = 0.01, x_pu = 0.05, b_pu = -0.002, g_pu = 0.005, status = 1, ratedS = 100.0, ratio = 1.0, shift_deg = 0.0)

  # FOR/DTF transformer G is represented by the native branch PI conductance.
  # The proprietary MATPOWER metadata preserves it for Sparlectra round trips.
  net.matpower_branch_metadata[1] = (
    orig_name = "T1A FROM -> TO",
    source_label = "T1A",
    orig_kind = :transformer,
    orig_index = 1,
    dtf_kind = 'T',
    u_ref_kV = 110.0,
    from_bus_vn_kV = 110.0,
    to_bus_vn_kV = 110.0,
    r_ohm = 1.21,
    x_ohm = 6.05,
    g_s = 4.132231404958678e-5,
    b_s = -1.652892561983471e-5,
    g_pu = 0.005,
    b_pu = -0.002,
    active_no_load_g_pu = 0.005,
    transformer_loss_allocation = :native_branch_pi,
    tap_ratio = 1.0,
    actual_tap_step = 0,
    max_tap_step = 9,
    phase_shift_deg = 0.0,
  )

  path = joinpath(output_dir, "transformer_loss_extension.m")
  Sparlectra.writeMatpowerCasefile(net, path)
  mpc = Sparlectra.MatpowerIO.read_case_m(path; legacy_compat = false)
  roundtrip = Sparlectra.createNetFromMatPowerCase(mpc = mpc, apply_bus_names = true, apply_branch_names = true, apply_branch_kind = true)
  return (path = path, transformer_loss_records = length(mpc.sparlectra.transformer_losses), total_g_pu = roundtrip.branchVec[1].g_pu)
end

if abspath(PROGRAM_FILE) == @__FILE__
  result = Base.invokelatest(main)
  println("Wrote ", result.path)
  println("Transformer-loss records: ", result.transformer_loss_records)
  println("Round-trip total G pu: ", result.total_g_pu)
end
