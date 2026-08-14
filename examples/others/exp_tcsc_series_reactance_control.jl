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

# file: examples/others/exp_tcsc_series_reactance_control.jl
# purpose: TCSC-like series-reactance controller (issue #297) steering the
#          flow split of a two-corridor loop network onto a branch target
#          (in-range secant regulation and the honest at_limit clamping),
#          plus the controller row and the generic controllable-element view.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

Two parallel corridors carry 80 MW from a source to a sink area; the
corridor reactance ratio 1:2 puts only one third of the transfer on the
upper corridor. A TCSC on that corridor (`addSeriesReactanceControl!`)
lowers the branch series reactance within `[0.02, 0.30]` p.u. until the
corridor carries the 35 MW target: the outer control loop measures the
branch flow, steps `x_pu` via secant iteration, and re-solves; every
accepted step re-stamps one branch entry of the Y-bus. A second variant
asks for 70 MW: the target is unreachable inside the range, the reactance
clamps at the capacitive bound and the controller reports `at_limit`; the
branch then behaves as a fixed compensated line. The generic
controllable-element view (`controllableElements`) describes both runs in
the shared actuator/target/limits vocabulary.
"""
function main()
  print_example_banner("examples/others/exp_tcsc_series_reactance_control.jl", "TCSC series-reactance control on a loop network: in-range flow steering, honest at_limit clamping, and the generic controllable-element view")

  build = function ()
    net = Net(name = "tcsc_demo", baseMVA = 100.0)
    for bus in ("A", "M1", "M2", "B")
      addBus!(net = net, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 80.0, q = 20.0)
    # lower corridor: A - M1 - B, x = 0.10 per line
    addPIModelACLine!(net = net, fromBus = "A", toBus = "M1", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "M1", toBus = "B", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    # upper corridor: A - M2 - B, x = 0.20 per line, the TCSC sits on A - M2
    addPIModelACLine!(net = net, fromBus = "A", toBus = "M2", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "M2", toBus = "B", r_pu = 0.02, x_pu = 0.20, b_pu = 0.0, status = 1)
    return net
  end

  run_variant = function (; p_target_mw::Float64, print_classic::Bool = false)
    net = build()
    run_sparlectra(net = net)
    p0 = get_branch_p_from_to_mw(net, "A", "M2")
    ctrl = addSeriesReactanceControl!(net; fromBus = "A", toBus = "M2", p_target_mw = p_target_mw, x_min_pu = 0.02, x_max_pu = 0.30)
    result = run_sparlectra(net = net)
    if print_classic
      printACPFlowResults(net, result.elapsed_s, result.iterations, 1e-8)
    end
    cr = latest_control_result(net)
    return (
      p_before = p0,
      p_after = get_branch_p_from_to_mw(net, "A", "M2"),
      x_pu = ctrl.x_pu,
      status = ctrl.status,
      at_limit = ctrl.at_limit,
      controller_rows = cr.controllers,
      loop_status = cr.status,
      outer = cr.outer_iterations,
      elements = controllableElements(net),
    )
  end

  in_range = run_variant(p_target_mw = 35.0, print_classic = true)
  limited = run_variant(p_target_mw = 70.0)
  return (in_range = in_range, limited = limited)
end

result = run_example(main)
println("TCSC target 35 MW: P(A->M2) ", round(result.in_range.p_before; digits = 2), " to ", round(result.in_range.p_after; digits = 2), " MW at x_pu = ", round(result.in_range.x_pu; digits = 4), "  status=", result.in_range.status, " (outer loop: ", result.in_range.loop_status, ", ", result.in_range.outer, " iterations)")
println("TCSC target 70 MW: P(A->M2) ", round(result.limited.p_before; digits = 2), " to ", round(result.limited.p_after; digits = 2), " MW at x_pu = ", round(result.limited.x_pu; digits = 4), "  status=", result.limited.status, " (at_limit=", result.limited.at_limit, ": fixed compensated line at the range end)")
println()
println("Controller rows (latest_control_result):")
for row in result.in_range.controller_rows
  println("  ", row.controller_name, ": achieved ", ismissing(row.achieved_p_mw) ? "-" : round(row.achieved_p_mw; digits = 2), " MW of ", row.p_target_mw, " MW target, x_pu = ", round(row.x_pu; digits = 4), " in [", row.x_min_pu, ", ", row.x_max_pu, "]")
end
println()
println("Controllable elements (generic view):")
for e in result.in_range.elements
  println("  ", e.name, ": ", e.device, ", actuator ", e.actuator, " in [", e.actuator_min, ", ", e.actuator_max, "], ", e.quantity, " @ ", e.target, " -> ", e.target_value, "  status=", e.status)
end
