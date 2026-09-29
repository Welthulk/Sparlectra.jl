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

# file: docs/lit/warmup/workshop_transformers.jl
# purpose: compile warm-up of the transformer workshop
#          (docs/lit/workshop_transformers.jl): runs every path the chapters
#          use once on throwaway nets, so the chapters do not stall on
#          first-call compilation. Included by the first code cell of the
#          notebook, which defines the helpers used here (solve!, bus_vm,
#          feeder, pst_loop, build_3wt). Not part of the library.

"""
    warmup()

Run every compile-heavy path of the transformer workshop once and print how
long the first power flow, the control loop and the rest took. The output
of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## warm both paths on a throwaway 2-bus net with a transformer
  wnet = Net(name = "warmup", baseMVA = 100.0)
  addBus!(net = wnet, busName = "A", vn_kV = 110.0)
  addBus!(net = wnet, busName = "B", vn_kV = 20.0)
  addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelTrafo!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  t_pf = @elapsed solve!(wnet)
  addPowerTransformerControl!(wnet; trafo = "1", mode = :branch_active_power, target_branch = ("A", "B"), p_target_mw = 10.0, control_ratio = false, control_phase = true, deadband_p_mw = 5.0)
  wnet.branchVec[1].has_phase_tap = true
  wnet.branchVec[1].phase_min_deg = -10.0
  wnet.branchVec[1].phase_max_deg = 10.0
  wnet.branchVec[1].phase_step_deg = 0.5
  t_ctrl = @elapsed run_sparlectra(net = wnet)
  println("warm: power flow ", round(t_pf; digits = 2), " s, control loop ", round(t_ctrl; digits = 2), " s (first calls compile)")

  ## every further path of the chapters once, output discarded: the tap
  ## device formulas, the builders, the phase-tap loop with a reactance
  ## characteristic on the winding, and the 3WT star equivalent
  t_paths = @elapsed redirect_stdout(devnull) do
    taps = PowerTransformerTaps(Vn_kV = 380.0, step = 1, lowStep = -9, highStep = 9, neutralStep = 0, voltageIncrement_kV = 3.8)
    wf = feeder(calcRatioTapCorrection(taps))
    println(bus_vm(wf, "B"))
    wtap = calcPhaseTapAngleRatio(PhaseTapChangerModel(kind = :asymmetrical, step = 1, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.0125, winding_connection_angle_deg = 60.0))
    wp = pst_loop(wtap.effective_ratio, wtap.effective_shift_deg)
    println(get_branch_p_from_to_mw(wp, "S", "M"))
    wc = pst_loop(1.0, 0.0)
    wbr = getNetBranch(net = wc, fromBus = "S", toBus = "M")
    wbr.has_phase_tap = true
    wbr.phase_min_deg = -10.0
    wbr.phase_max_deg = 10.0
    wbr.phase_step_deg = 0.5
    wc.trafos[1].side1.phase_taps = Sparlectra.PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, x_min = 0.08, x_max = 0.16)
    run_sparlectra(net = wc)
    addPowerTransformerControl!(wc; trafo = string(wbr.branchIdx), mode = :branch_active_power, target_branch = ("S", "M"), p_target_mw = get_branch_p_from_to_mw(wc, "S", "M") - 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 4.0)
    run_sparlectra(net = wc)
    println(bus_vm(build_3wt(oltc_step = 1), "B3"))
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (tap formulas, builders, X(alpha) control loop, 3WT); everything warm")
  return nothing
end
