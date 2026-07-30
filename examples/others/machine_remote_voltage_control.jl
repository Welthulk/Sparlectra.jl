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

# file: examples/others/machine_remote_voltage_control.jl
# purpose: demo of MachineVoltageControl — a PQ machine whose reactive output
#          is driven by the outer control loop until the voltage at a REMOTE
#          bus reaches its target, including the honest at_limit outcome when
#          the reactive range cannot deliver the target.

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Three-bus chain: Slack — GenBus (machine) — Load. The machine's regulated
# bus is NOT its own connection point, which is exactly the situation a CGMES
# RegulatingControl at a foreign terminal describes (#294 point 3).
function build_rvc_net(; qmin_mvar::Float64 = -50.0, qmax_mvar::Float64 = 50.0)
  net = Net(name = "machine_rvc_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "GenBus", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "GenBus", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 0.0, qMin = qmin_mvar, qMax = qmax_mvar)
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "GenBus", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "GenBus", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
  return net
end

function run_case(; title::String, target_vm_pu::Float64, qmin_mvar::Float64, qmax_mvar::Float64)
  println("\n== ", title, " ==")
  net = build_rvc_net(qmin_mvar = qmin_mvar, qmax_mvar = qmax_mvar)

  pf = PowerFlowConfig(max_iter = 30, tol = 1e-9)
  _, erg = runpf!(net; config = pf, verbose = 0)
  erg == 0 || error("Uncontrolled PF failed")
  @printf("uncontrolled : Vm(Load) = %.5f pu, machine Q = %.2f MVAr\n", get_bus_vm_pu(net, "Load"), net.prosumpsVec[2].qVal)

  addMachineVoltageControl!(net; bus = "GenBus", target_bus = "Load", target_vm_pu = target_vm_pu, deadband_vm_pu = 5e-4)
  result = run_control!(net; controllers = collect_outer_controllers(net), pf_config = pf, control_config = ControlConfig(max_outer_iterations = 15), verbose = 0)
  @printf("controlled   : Vm(Load) = %.5f pu (target %.4f), machine Q = %.2f MVAr\n", get_bus_vm_pu(net, "Load"), target_vm_pu, net.prosumpsVec[2].qVal)
  @printf("loop         : status=%s, outer=%d, pf_solves=%d\n", result.status, result.outer_iterations, result.powerflow_solves)
  printMachineControllerSummary(stdout, net)
  return net
end

function main(args::AbstractVector{String} = ARGS)
  print_example_banner("examples/others/machine_remote_voltage_control.jl", "remote voltage control via machine reactive power (MachineVoltageControl + run_control!)")

  # Reachable target: the secant loop settles inside the deadband.
  run_case(title = "reachable target", target_vm_pu = 1.05, qmin_mvar = -50.0, qmax_mvar = 50.0)

  # Unreachable target: with the reactive range cut to ±2 MVAr the controller
  # parks honestly at its limit instead of pretending convergence — the same
  # physics as a PV machine switching to PQ under Q limits.
  run_case(title = "unreachable target (Q range ±2 MVAr)", target_vm_pu = 1.00, qmin_mvar = -2.0, qmax_mvar = 2.0)

  println("\nTip: on CGMES deliveries the same controllers are attached by importCGMES(machine_control = true)")
  println("     (config key cgmes_import.machine_control) for machines whose RegulatingControl points at a foreign bus.")
end

run_example(main, ARGS)
