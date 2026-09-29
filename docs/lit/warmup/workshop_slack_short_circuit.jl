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

# file: docs/lit/warmup/workshop_slack_short_circuit.jl
# purpose: compile warm-up of the slack and short-circuit workshop
#          (docs/lit/workshop_slack_short_circuit.jl): runs every path the
#          scenarios use once on throwaway nets, so the scenarios do not
#          stall on first-call compilation. Included by the first code cell
#          of the notebook, which defines the helpers used here (solve!,
#          build_grid). Not part of the library.

"""
    warmup()

Run every compile-heavy path of the slack and short-circuit workshop once
and print how long the first solve and the rest took. The output of the
paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## tiny warm-up net: a grid connection WITH declared short-circuit data
  wnet = Net(name = "warmup", baseMVA = 100.0)
  addBus!(net = wnet, busName = "A", vn_kV = 110.0)
  addBus!(net = wnet, busName = "B", vn_kV = 110.0)
  addExternalGrid!(net = wnet, busName = "A", vm_pu = 1.0, sk_max_MVA = 2000.0, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  t_pf = @elapsed solve!(wnet)
  t_sc = @elapsed runShortCircuit!(wnet; case = :max)
  println("warm: power flow ", round(t_pf; digits = 2), " s, short circuit ", round(t_sc; digits = 2), " s (first calls compile)")

  ## every further path of the notebook once, output discarded: the three
  ## slack models with their result tables, the loss summary, and both
  ## short-circuit cases with their printout
  t_paths = @elapsed redirect_stdout(devnull) do
    for (mode, kw) in ((:slack, (;)), (:source, (;)), (:slack, (distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)))
      wg = build_grid(mode)
      etime, ite = solve!(wg; kw...)
      printACPFlowResults(wg, etime, ite, 1e-8)
      getTotalLosses(net = wg)
      get_bus_vm_pu(wg, "B1")
    end
    wg = build_grid(:slack)
    solve!(wg)
    printShortCircuitResult(runShortCircuit!(wg; case = :max))
    printShortCircuitResult(runShortCircuit!(wg; case = :min))
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (result tables, three slack models, short circuit max/min); everything warm")
  return nothing
end
