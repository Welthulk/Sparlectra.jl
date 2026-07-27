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
#
# Example: split Schrägregelung — TWO independent controllers on ONE transformer
#
# A Schrägregler combines in-phase (Längs-) and quadrature (Quer-) regulation
# on a single unit. In real installations the two actuators typically have
# SEPARATE controllers: a voltage controller driving the ratio tap and an
# active-power controller driving the phase tap. This demo models exactly
# that split (contrast: tap_control_demo_grid.jl uses ONE combined controller
# in mode :voltage_and_branch_active_power for its T_SCHRAEG unit):
#
# 1) Uncontrolled power flow — record voltage and through-flow.
# 2) Attach two PowerTransformerControl instances to the same transformer:
#    - mode = :voltage,             control_ratio = true  (OLTC channel)
#    - mode = :branch_active_power, control_phase = true  (PST channel)
#    Each has its own target, deadband, and convergence status; the outer
#    control loop alternates them until both are inside their deadbands.
# 3) Show the per-actuator exclusivity rule: a second controller on an
#    already-driven actuator is rejected.
# 4) Controlled run — report per-controller rows and achieved values.

# Date: 2026-07-27
# file: examples/others/tap_control_schraeg_two_controllers.jl
# purpose: split Schrägregelung — independent voltage (ratio tap) and active-power (phase tap) controllers on one transformer, with per-actuator exclusivity

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

function build_net()
  net = Net(name = "schraeg_split_demo", baseMVA = 100.0)
  for bus in ("Slack", "Mid", "Load")
    addBus!(net = net, busName = bus, vn_kV = 110.0)
  end

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)

  # The Schrägregler unit: one transformer carrying BOTH a ratio tap and a
  # phase tap. A parallel line gives the phase tap a loop to shift flow into.
  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.03, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

  tbr = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
  tbr.comp.cName = "T_SCHRAEG"
  tbr.has_ratio_tap = true
  tbr.tap_min = 0.95
  tbr.tap_max = 1.08
  tbr.tap_step = 0.0125
  tbr.has_phase_tap = true
  tbr.phase_min_deg = -10.0
  tbr.phase_max_deg = 10.0
  tbr.phase_step_deg = 0.5

  return net
end

function main()
  print_example_banner(
    "examples/others/tap_control_schraeg_two_controllers.jl",
    "split Schrägregelung: independent voltage (ratio tap) and active-power (phase tap) controllers on one transformer, with per-actuator exclusivity",
  )

  cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-9), output = OutputConfig(logfile_results = :off))

  net = build_net()
  r0 = run_sparlectra(net = net, config = cfg)
  r0.numerical_converged || error("Uncontrolled PF failed: $(r0.reason_text)")
  vm0 = get_bus_vm_pu(net, "Load")
  p0 = get_branch_p_from_to_mw(net, "Slack", "Mid")

  println("Uncontrolled state:")
  @printf("  Vm(Load)        = %.4f pu\n", vm0)
  @printf("  P(Slack -> Mid) = %.2f MW (through the Schrägregler unit)\n\n", p0)

  # Controller 1 — OLTC channel: hold the Load-bus voltage via the ratio tap.
  vm_target = 1.03
  addTapController!(net;
    trafo = "T_SCHRAEG",
    mode = :voltage,
    target_bus = "Load",
    target_vm_pu = vm_target,
    control_ratio = true,
    control_phase = false,
    deadband_vm_pu = 5e-3,
  )

  # Controller 2 — PST channel: shift 8 MW off the transformer path via the
  # phase tap. Same unit, disjoint actuator, own deadband and status.
  # Discrete-step sizing rule: one 0.5° phase step moves roughly 2.5 MW here,
  # so the deadband must cover at least half a step's effect — otherwise the
  # discrete controller hunts around the target without ever converging.
  p_target = round(p0) + 8.0
  addTapController!(net;
    trafo = "T_SCHRAEG",
    mode = :branch_active_power,
    target_branch = ("Slack", "Mid"),
    p_target_mw = p_target,
    control_ratio = false,
    control_phase = true,
    deadband_p_mw = 2.0,
  )

  # Per-actuator exclusivity: the ratio tap is already driven by controller 1,
  # so a third controller claiming it is rejected.
  println("Per-actuator exclusivity check:")
  try
    addTapController!(net; trafo = "T_SCHRAEG", mode = :voltage, target_bus = "Load", target_vm_pu = 1.0)
    error("expected the duplicate-actuator controller to be rejected")
  catch err
    println("  rejected as expected: ", first(split(sprint(showerror, err), ';')))
  end
  println()

  r1 = run_sparlectra(net = net, config = cfg)
  r1.numerical_converged || error("Controlled PF failed: $(r1.reason_text)")
  result = latest_control_result(net)
  result === nothing && error("Expected control result on net.control_result")

  println("Controlled state (targets: Vm = $(vm_target) pu, P = $(p_target) MW):")
  @printf("  Vm(Load)        = %.4f pu\n", get_bus_vm_pu(net, "Load"))
  @printf("  P(Slack -> Mid) = %.2f MW\n", get_branch_p_from_to_mw(net, "Slack", "Mid"))
  println("  Control status: ", result.status, ", outer iterations: ", result.outer_iterations, ", PF solves: ", result.powerflow_solves)
  println()

  println("Per-controller result rows (one row per channel):")
  for row in result.controllers
    @printf(
      "  %-22s mode=%-20s status=%-10s tap_ratio=%.4f phase=%.1f°  target/achieved: %s\n",
      row.controller_name,
      row.mode,
      row.status,
      row.tap_ratio,
      row.phase_shift_deg,
      row.mode == "voltage" ? @sprintf("%.4f / %.4f pu", row.target_vm_pu, row.achieved_vm_pu) : @sprintf("%.1f / %.2f MW", row.p_target_mw, row.achieved_p_mw),
    )
  end
  println()
  println("Interpretation:")
  println("  Both actuators sit on the SAME transformer, yet each channel keeps its")
  println("  own target, deadband, and convergence flag — the outer loop alternates")
  println("  them. Voltage couples mainly to the ratio tap and active power to the")
  println("  phase tap, so the alternating scheme settles despite the shared unit.")
  return nothing
end

run_example(main)
