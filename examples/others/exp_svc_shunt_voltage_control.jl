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

# file: examples/others/exp_svc_shunt_voltage_control.jl
# purpose: SVC-style variable-shunt voltage controller holding a bus voltage
#          (in-range secant regulation and the honest constant-B limit
#          region), plus the generic controllable-element view.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

A lightly loaded feeder lifts the far bus above 1.1 p.u. An SVC at that bus
(`addShuntVoltageControl!`) moves its shunt susceptance within
±60 MVAr to hold 1.0 p.u. — the outer control loop measures, steps the
susceptance via secant iteration, and re-solves. A second variant gets only
±10 MVAr: the target is unreachable, the susceptance clamps at the inductive
bound and the controller reports `at_limit` — the reactive output then
follows V² through the Y-bus (the constant-B region of a real SVC), it does
not pretend to hold a constant Q. The generic controllable-element view
(`controllableElements`) describes both runs in the shared
actuator/target/limits vocabulary.
"""
function main()
  print_example_banner("examples/others/exp_svc_shunt_voltage_control.jl", "SVC-style variable-shunt voltage control: in-range regulation, honest at_limit clamping, and the generic controllable-element view")

  build = function ()
    svcnet = Net(name = "svc_demo", baseMVA = 100.0)
    for bus in ("Slack", "Mid", "Load")
      addBus!(net = svcnet, busName = bus, vn_kV = 110.0)
    end
    addProsumer!(net = svcnet, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
    addProsumer!(net = svcnet, busName = "Load", type = "ENERGYCONSUMER", p = -80.0, q = -30.0)
    addPIModelACLine!(net = svcnet, fromBus = "Slack", toBus = "Mid", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    addPIModelACLine!(net = svcnet, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    return svcnet
  end

  run_variant = function (; bs_range_mvar::Float64)
    net = build()
    run_sparlectra(net = net)
    vm0 = get_bus_vm_pu(net, "Load")
    addShuntVoltageControl!(net; bus = "Load", target_vm_pu = 1.0, bs_min_mvar = -bs_range_mvar, bs_max_mvar = bs_range_mvar)
    run_sparlectra(net = net)
    ctrl = only(Sparlectra._shunt_controllers(net))
    return (vm_before = vm0, vm_after = get_bus_vm_pu(net, "Load"), bs_mvar = ctrl.bs_mvar, status = ctrl.status, at_limit = ctrl.at_limit, elements = controllableElements(net))
  end

  in_range = run_variant(bs_range_mvar = 60.0)
  limited = run_variant(bs_range_mvar = 10.0)
  return (in_range = in_range, limited = limited)
end

result = run_example(main)
println("SVC ±60 MVAr: Vm(Load) ", round(result.in_range.vm_before; digits = 4), " → ", round(result.in_range.vm_after; digits = 4), " pu at Bs = ", round(result.in_range.bs_mvar; digits = 2), " MVAr  status=", result.in_range.status)
println("SVC ±10 MVAr: Vm(Load) ", round(result.limited.vm_before; digits = 4), " → ", round(result.limited.vm_after; digits = 4), " pu at Bs = ", round(result.limited.bs_mvar; digits = 2), " MVAr  status=", result.limited.status, " (at_limit=", result.limited.at_limit, " — constant-B region, Q follows V²)")
println()
println("Controllable elements (generic view):")
for e in result.in_range.elements
  println("  ", e.name, ": ", e.device, " — actuator ", e.actuator, " ∈ [", e.actuator_min, ", ", e.actuator_max, "], ", e.quantity, " @ ", e.target, " → ", e.target_value, "  status=", e.status)
end
