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

# file: docs/lit/warmup/workshop_tour_control.jl
# purpose: compile warm-up of the parallel-transformer control workshop
#          (docs/lit/workshop_tour_control.jl): runs every path the
#          examples use once on throwaway substations, so the examples do
#          not stall on first-call compilation. Included by the first code
#          cell of the notebook, which defines the helpers used here
#          (build_substation, solve!, trafo_ids). Not part of the library.

"""
    warmup()

Run every compile-heavy path of the parallel-transformer control workshop
once and print how long the first power flow and the tap-group control
loop took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  wnet = build_substation("warmup")
  t_pf = @elapsed solve!(wnet)
  t_ctrl = @elapsed redirect_stdout(devnull) do
    wg = build_substation("warmup_group")
    addPowerTransformerControl!(wg; trafo = trafo_ids[1], followers = [trafo_ids[2]], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.03, deadband_vm_pu = 5e-3)
    wres = run_control!(wg)
    println(wres.status, get_bus_vm_pu(wg, "Feeder"), Sparlectra._find_trafo_branch(wg, trafo_ids[1]).tap_ratio)
    calcNetLosses!(wg)
    for e in controllableElements(wg)
      println(e.element, e.quantity, e.target, e.status)
    end
    wi = build_substation("warmup_single")
    addPowerTransformerControl!(wi; trafo = trafo_ids[1], mode = :voltage, target_bus = "Feeder", target_vm_pu = 1.0)
  end
  println("warm: power flow ", round(t_pf; digits = 2), " s, tap-group control loop ", round(t_ctrl; digits = 2), " s (first calls compile); everything warm")
  return nothing
end
