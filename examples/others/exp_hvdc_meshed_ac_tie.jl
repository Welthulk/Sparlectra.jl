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

# Date: 2026-08-19
# file: examples/others/exp_hvdc_meshed_ac_tie.jl
# purpose: HVDC pair in parallel to an AC tie (r0.9.9): one angle reference
#          per synchronous island (the multi-reference error, caught and
#          shown), the setpoint pair as a parallel PQ path next to the tie,
#          a retarget shifting the exchange between link and tie, and the
#          honest invalid_topology outcome for island_feed in a synchronous
#          grid. Theory: docs/src/hvdc_back_to_back.md, "Meshed operation".

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# The workshop-tour B2B topology (A1, A2, C1, C2) plus the AC tie A1 -> C1.
# c1_model picks how the former island-C reference is modeled once the tie
# closes: :reference keeps it (scene 1, must fail), :pv demotes it to a
# voltage-regulated generator (scenes 2 to 4).
function build_meshed(name::String; c1_model::Symbol, transfer::Float64 = 80.0, loss::Float64 = 4.0)
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "C1", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  if c1_model === :reference
    addProsumer!(net = net, busName = "C1", type = "EXTERNALNETWORKINJECTION", referencePri = "C1", vm_pu = 1.0, va_deg = 0.0)
  else
    addProsumer!(net = net, busName = "C1", type = "GENERATOR", p = 20.0, q = 0.0, vm_pu = 1.0, isRegulated = true)
  end
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "C2", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net, busName = "A2", type = "GENERATOR", p = -transfer, q = 0.0)  ## converter, exports
  addProsumer!(net = net, busName = "C2", type = "GENERATOR", p = transfer - loss, q = 0.0)  ## converter, receives
  return net
end

function print_link_and_tie(net::Net, label::AbstractString)
  tie = get_branch_p_from_to_mw(net, "A1", "C1")
  rows = Sparlectra._hvdc_link_flow_rows(net)
  link = isempty(rows) ? nothing : rows[1]
  println(label)
  link === nothing || println("  HVDC link : ", round(link.p_from_MW; digits = 3), " MW ordered transfer, ", round(link.loss_MW; digits = 3), " MW loss (", link.mode, ")")
  println("  AC tie    : ", round(tie; digits = 3), " MW on A1 -> C1")
end

function main()
  print_example_banner("examples/others/exp_hvdc_meshed_ac_tie.jl", "HVDC pair in parallel to an AC tie: one reference per synchronous island")

  println("=== Scene 1: both references kept, AC tie closed -> rejected ===")
  net1 = build_meshed("meshed_two_refs"; c1_model = :reference)
  try
    runpf!(net1, 30, 1e-8, 0; method = :rectangular, islands_enabled = true)
    error("scene 1 unexpectedly solved")
  catch err
    println(sprint(showerror, err))
  end

  println()
  println("=== Scene 2: C1 demoted to PV, Stage-0 injections in parallel to the tie ===")
  net2 = build_meshed("meshed_stage0"; c1_model = :pv)
  # register the hand-built Stage-0 link so the HVDC Link Flows table and
  # the report carry it (importers do this automatically)
  addHvdcLink!(net2; from_bus = "A2", to_bus = "C2")
  ite, erg = runpf!(net2, 30, 1e-8, 0; method = :rectangular, islands_enabled = true)
  erg == 0 || error("Stage-0 power flow did not converge")
  calcNetLosses!(net2)
  printACPFlowResults(net2, 0.0, ite, 1e-8)
  print_link_and_tie(net2, "Stage-0 snapshot (fixed 80 MW):")

  println()
  println("=== Scene 3: setpoint pair (120 MW) as a parallel PQ path ===")
  net3 = build_meshed("meshed_paired"; c1_model = :pv)
  ctrl = addHvdcPairControl!(net3; from_bus = "A2", to_bus = "C2", p_transfer_mw = 120.0, loss_mw = 4.0, p_rating_mw = 150.0)
  result = run_control!(net3; controllers = collect_outer_controllers(net3), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 10, trace = false))
  println("outer loop: ", result.status, " after ", result.outer_iterations, " iteration(s)")
  calcNetLosses!(net3)
  printACPFlowResults(net3, 0.0, result.last_pf_iterations, 1e-8)
  print_link_and_tie(net3, "Retargeted to 120 MW (invariant P_to = 120 - 4):")

  println()
  println("=== Scene 4: retarget to 40 MW, the tie takes over ===")
  ctrl.p_transfer_mw = 40.0
  ctrl.p_applied = false
  result2 = run_control!(net3; controllers = collect_outer_controllers(net3), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 10, trace = false))
  println("outer loop: ", result2.status, " after ", result2.outer_iterations, " iteration(s)")
  calcNetLosses!(net3)
  print_link_and_tie(net3, "Retargeted to 40 MW (the A/C exchange moves onto the tie):")

  println()
  println("=== Scene 5: island_feed in a synchronous grid -> rejected, then invalid_topology ===")
  # grid-forming reference at C2 plus the AC tie: two references in one
  # synchronous island, same early error as scene 1
  net5 = Net(name = "meshed_island_feed", baseMVA = 100.0)
  for b in ("A1", "A2", "C1", "C2")
    addBus!(net = net5, busName = b, vn_kV = 380.0)
  end
  addPIModelACLine!(net = net5, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net5, fromBus = "C2", toBus = "C1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net5, fromBus = "A1", toBus = "C1", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
  addProsumer!(net = net5, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net5, busName = "A2", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net5, busName = "C2", type = "EXTERNALNETWORKINJECTION", referencePri = "C2", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net5, busName = "C1", type = "ENERGYCONSUMER", p = 50.0, q = 12.0)
  addProsumer!(net = net5, busName = "A2", type = "GENERATOR", p = 0.0, q = 0.0)
  ctrl5 = addHvdcPairControl!(net5; from_bus = "A2", to_bus = "C2", mode = :island_feed, loss_mw = 4.0)
  try
    run_control!(net5; controllers = collect_outer_controllers(net5), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
    error("scene 5 unexpectedly solved with two references")
  catch err
    println(sprint(showerror, err))
  end
  # demoting the grid-forming reference resolves the reference conflict but
  # invalidates the island_feed semantics: the controller says so honestly
  for ps in net5.prosumpsVec
    if Sparlectra.getPosumerBusIndex(ps) == net5.busDict["C2"] && ps.referencePri !== nothing
      ps.referencePri = nothing
    end
  end
  refreshBusTypesFromProsumers!(net5)
  run_control!(net5; controllers = collect_outer_controllers(net5), pf_config = PowerFlowConfig(method = :rectangular, max_iter = 30, tol = 1e-8), control_config = ControlConfig(max_outer_iterations = 8, trace = false))
  println("controller status after demotion: ", ctrl5.status, " (grid-forming in a synchronous grid is a different model)")
  printHvdcPairControllerSummary(net5)

  return (stage0 = net2, paired = net3, island_feed = net5)
end

run_example(main)
