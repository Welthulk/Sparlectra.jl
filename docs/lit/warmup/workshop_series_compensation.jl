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

# file: docs/lit/warmup/workshop_series_compensation.jl
# purpose: compile warm-up of the series compensation workshop
#          (docs/lit/workshop_series_compensation.jl): runs every path the
#          examples use once on throwaway nets, so the examples do not
#          stall on first-call compilation. Included by the first code cell
#          of the notebook, which defines the builders used here
#          (build_loop, build_upfc_mesh). Not part of the library.

"""
    warmup()

Run every compile-heavy path of the series compensation workshop once and
print how long the first controlled power flow and the rest took. The
output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  wnet = Net(name = "warmup", baseMVA = 100.0)
  for b in ("A", "M", "B")
    addBus!(net = wnet, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "M", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = wnet, fromBus = "M", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
  addSeriesReactanceControl!(wnet; fromBus = "A", toBus = "M", p_target_mw = 6.0, x_min_pu = 0.05, x_max_pu = 0.2)
  t_ctrl = @elapsed run_control!(wnet; controllers = collect_outer_controllers(wnet), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 15, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 4, trace = false))
  println("warm: power flow plus series-reactance control ", round(t_ctrl; digits = 2), " s (first calls compile)")

  ## every further path of the examples once, output discarded: the
  ## run_sparlectra control loop, the control result and element views, the
  ## result tables, the quadrature and the full UPFC
  t_paths = @elapsed redirect_stdout(devnull) do
    wl = build_loop()
    run_sparlectra(net = wl)
    addSeriesReactanceControl!(wl; fromBus = "A", toBus = "M2", p_target_mw = 30.0, x_min_pu = 0.02, x_max_pu = 0.30)
    run_sparlectra(net = wl)
    println(latest_control_result(wl).status, get_branch_p_from_to_mw(wl, "A", "M2"))
    foreach(println, controllableElements(wl))
    calcNetLosses!(wl)
    printACPFlowResults(wl, 0.0, 1, 1e-8)
    wu = build_loop()
    addProsumer!(net = wu, busName = "M2", type = "GENERATOR", p = 0.0, q = 0.0)
    addUpfcControl!(wu; fromBus = "A", toBus = "M2", shunt_bus = "M2", target_bus = "B", target_vm_pu = 0.99, p_target_mw = 35.0, v_inj_max_pu = 0.08, s_max_mva = 40.0)
    run_control!(wu)
    calcNetLosses!(wu)
    printACPFlowResults(wu, 0.0, 1, 1e-8)
    wf = build_upfc_mesh()
    addUpfcControl!(wf; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I", p_target_mw = 40.0, q_target_mvar = 10.0, q_shunt_mvar = 0.0, v_inj_max_pu = 0.30, s_max_mva = 120.0, deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
    run_control!(wf; control_config = ControlConfig(max_outer_iterations = 80))
    calcNetLosses!(wf)
    printACPFlowResults(wf, 0.0, 1, 1e-8)
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (control runs, result tables, both UPFC models); everything warm")
  return nothing
end
