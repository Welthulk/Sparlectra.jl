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

# file: examples/others/exp_pst_reactance_coupling.jl
# purpose: phase-shifting transformer whose series reactance follows the
#          device characteristic X(α) while the outer control loop moves the
#          tap — compared against the same control run with a static
#          reactance.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

Build a small loop network (transformer + parallel line) with a PST that
regulates the active power through the transformer. The PST winding carries a
symmetrical `PhaseTapChangerModel` with `x_min`/`x_max`: every accepted tap
move of the outer control loop then updates the branch series reactance to
X(α) at the new angle before the next solve. The same run without the typed
model keeps the static import-time reactance — the comparison shows both the
different converged reactance and the different flow picture.
"""
function main()
  print_example_banner("examples/others/exp_pst_reactance_coupling.jl", "PST outer-loop control with tap-dependent series reactance X(α) vs. a static-reactance run")

  build = function (; with_model::Bool)
    net = Net(name = "pst_xalpha_demo", baseMVA = 100.0)
    for bus in ("Slack", "Mid", "Load")
      addBus!(net = net, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
    addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
    tbr.has_phase_tap = true
    tbr.phase_min_deg = -10.0
    tbr.phase_max_deg = 10.0
    tbr.phase_step_deg = 0.5
    if with_model
      # symmetrical PST: X grows from x_min at neutral to x_max at the
      # range end (X(α) = x_min + (x_max - x_min)·(sin(α/2)/sin(αmax/2))²)
      net.trafos[1].side1.phase_taps = Sparlectra.PhaseTapChangerModel(kind = :symmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, x_min = 0.08, x_max = 0.16)
    end
    return net, tbr
  end

  run_variant = function (; with_model::Bool)
    net, tbr = build(with_model = with_model)
    baseline = run_sparlectra(net = net)
    p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")
    addPowerTransformerControl!(net; trafo = string(tbr.branchIdx), mode = :branch_active_power, target_branch = ("Slack", "Mid"), p_target_mw = p0 - 8.0, control_ratio = false, control_phase = true, deadband_p_mw = 4.0)
    result = run_sparlectra(net = net)
    return (
      converged = baseline.numerical_converged && result.numerical_converged,
      p_before = p0,
      p_after = get_branch_p_from_to_mw(net, "Slack", "Mid"),
      phase_deg = tbr.phase_shift_deg,
      x_pu = tbr.x_pu,
    )
  end

  tracked = run_variant(with_model = true)
  static = run_variant(with_model = false)
  return (tracked = tracked, static = static)
end

result = run_example(main)
println("PST control with X(α) tracking: converged=", result.tracked.converged, "  P ", round(result.tracked.p_before; digits = 2), " → ", round(result.tracked.p_after; digits = 2), " MW  at α=", result.tracked.phase_deg, "°  x_pu=", round(result.tracked.x_pu; digits = 5))
println("PST control static reactance:   converged=", result.static.converged, "  P ", round(result.static.p_before; digits = 2), " → ", round(result.static.p_after; digits = 2), " MW  at α=", result.static.phase_deg, "°  x_pu=", round(result.static.x_pu; digits = 5))
println("The tracked run ends on the device characteristic (x_pu > 0.08); the static run keeps the import-time reactance.")
