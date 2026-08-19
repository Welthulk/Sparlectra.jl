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

# Date: 2026-08-20
# file: examples/others/exp_hvdc_b2b_pairing.jl
# purpose: back-to-back HVDC pairing controller demo (#297 Draft B): two AC
#          areas coupled only by a converter pair, solved as Stage-0 fixed
#          injections first, then with addHvdcPairControl! steering the
#          transfer to a new target while the pairing invariant holds, and
#          finally the grid-forming island_feed mode where the receiving
#          converter IS the island reference and the sending side mirrors
#          the island draw

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Two areas, each with its own reference, joined ONLY by the converter pair
# at B2 (from side, exports) and B4 (to side, receives). No AC tie exists,
# so the areas stay separate electrical islands in every mode.
function build_two_area_net(name::String; transfer::Float64, loss::Float64)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B3", type = "EXTERNALNETWORKINJECTION", referencePri = "B3", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  # the converter injections (Stage-0 snapshot: from side exports transfer,
  # to side receives transfer minus loss)
  addProsumer!(net = net, busName = "B2", type = "GENERATOR", p = -transfer, q = 0.0)
  addProsumer!(net = net, busName = "B4", type = "GENERATOR", p = transfer - loss, q = 0.0)
  return net
end

function area_balance(net::Net)
  slack1 = net.prosumpsVec[1]
  slack2 = net.prosumpsVec[2]
  return (b1 = something(slack1.pVal, 0.0), b3 = something(slack2.pVal, 0.0))
end

function main()
  print_example_banner("examples/others/exp_hvdc_b2b_pairing.jl", "back-to-back HVDC pairing controller: Stage-0 snapshot vs steerable transfer")

  println("=== Stage 0: fixed injections reproduce the snapshot ===")
  net0 = build_two_area_net("b2b_stage0"; transfer = 80.0, loss = 4.0)
  ite, erg = runpf!(net0, 30, 1e-8, 0; method = :rectangular, islands_enabled = true)
  erg == 0 || error("Stage-0 power flow did not converge")
  println("two islands solved in ", ite, " iteration(s); the link carries the fixed 80 MW snapshot")

  println()
  println("=== Paired control: steer the same link to a new transfer ===")
  net1 = build_two_area_net("b2b_paired"; transfer = 80.0, loss = 4.0)
  addHvdcPairControl!(net1; from_bus = "B2", to_bus = "B4", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
  result = run_control!(net1; controllers = collect_outer_controllers(net1), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 10, trace = true))
  println("outer loop: ", result.status, " after ", result.outer_iterations, " iteration(s), ", result.powerflow_solves, " solve(s)")
  printHvdcPairControllerSummary(net1)

  println()
  println("The pairing invariant after the retarget (injection convention):")
  m_from = net1.prosumpsVec[end-1]
  m_to = net1.prosumpsVec[end]
  println("  from side (B2): ", round(something(m_from.pVal, 0.0); digits = 3), " MW")
  println("  to side   (B4): ", round(something(m_to.pVal, 0.0); digits = 3), " MW = 120 - 4 loss")

  println()
  println("=== Element vocabulary ===")
  for e in controllableElements(net1)
    println("  ", e.name, ": actuator ", e.actuator, " target ", e.target, " = ", e.target_value, " (status ", e.status, ")")
  end

  println()
  println("=== Grid-forming: the link feeds an island without any other source ===")
  # island C has no classical slack: the receiving converter at C2 is the
  # island reference (grid-forming Vf), the load hangs at the far end C1,
  # and the controller mirrors the island draw plus loss onto the sender
  net2 = Net(name = "b2b_island_feed", baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net2, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net2, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net2, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addProsumer!(net = net2, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net2, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net2, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net2, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net2, busName = "A2", type = "GENERATOR", p = 0.0, q = 0.0)
  ctrl2 = addHvdcPairControl!(net2; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0, p_rating_mw = 150.0)
  result2 = run_control!(net2; controllers = collect_outer_controllers(net2), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 10, trace = true))
  println("outer loop: ", result2.status, " after ", result2.outer_iterations, " iteration(s)")
  printHvdcPairControllerSummary(net2)
  println("the island decided: draw ", round(ctrl2.p_transfer_mw - 4.0; digits = 3), " MW, mirrored as ", round(-ctrl2.p_transfer_mw; digits = 3), " MW on the sending side")

  return (stage0 = net0, paired = net1, result = result, island_feed = net2)
end

run_example(main)
