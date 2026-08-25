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

# Date: 2026-08-25
# file: examples/others/exp_facts_base_impedance.jl
# purpose: equipment impedance vs FACTS operating point (issue #329). A full
#          UPFC (model = :full) stamps a compensated series impedance onto the
#          live branch fields during the control run; its resistance part can go
#          negative. Short circuit and the MATPOWER/CGMES exports must read the
#          PHYSICAL base impedance instead, so the compensated net matches the
#          equipment net. Also shows restoreBaseImpedances! and
#          clearUpfcFullControllers!. Theory: docs/src/facts.md.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# One IEC 60909-0 network feeder record (max/min initial symmetrical currents).
_feeder(bus) = (mrid = bus, name = "F_" * bus, bus = bus,
  maxInitialSymShCCurrent_A = 20_000.0, minInitialSymShCCurrent_A = 16_000.0,
  maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1, maxR0ToX0Ratio = nothing,
  maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)

# A small meshed corridor: the S->L parallel path gives the UPFC on I->J real
# flow-steering freedom, so the series converter injects an active component and
# the live branch resistance is driven negative (the effect #329 is about).
function build_mesh()
  net = Net(name = "upfc_base_impedance", baseMVA = 100.0)
  for b in ("S", "I", "J", "L")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "I", type = "GENERATOR", p = 0.0, q = 0.0)          # shunt converter node
  addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 90.0, q = 30.0)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "I", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "I", toBus = "J", r_pu = 0.02, x_pu = 0.18, b_pu = 0.0, status = 1)   # UPFC line
  addPIModelACLine!(net = net, fromBus = "J", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "L", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)   # parallel path
  ok, msg = validate!(net = net)
  ok || error("mesh net invalid: $msg")
  return net
end

function main()
  print_example_banner("examples/others/exp_facts_base_impedance.jl", "equipment impedance vs FACTS operating point (#329): SC and export read the base")

  scd = Sparlectra.CGMESImporter.CGMESShortCircuitData([_feeder("S")], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])

  # 1) the equipment network (no UPFC): the short-circuit reference
  base = build_mesh()
  base_sc = runShortCircuit!(base, scd; case = :max)

  # 2) the same network with a full UPFC steering P and Q on I->J
  net = build_mesh()
  addUpfcControl!(net; model = :full, fromBus = "I", toBus = "J", shunt_bus = "I",
                  p_target_mw = 48.0, q_target_mvar = 8.0, q_shunt_mvar = 0.0,
                  v_inj_max_pu = 0.30, s_max_mva = 120.0,
                  deadband_p_mw = 1e-2, deadband_q_mvar = 1e-2, max_outer_iters = 80)
  run_control!(net; control_config = ControlConfig(max_outer_iterations = 80))
  brij = getNetBranch(net = net, fromBus = "I", toBus = "J")

  println("I->J branch after the UPFC control run:")
  println("  live   r/x = ", round(brij.r_pu; digits = 5), " / ", round(brij.x_pu; digits = 5),
          " pu   <- compensated operating point (r can be < 0)")
  println("  base   r/x = ", round(brij.r_base_pu; digits = 5), " / ", round(brij.x_base_pu; digits = 5),
          " pu   <- physical equipment, unchanged")
  println("  the power flow solves with the LIVE value; short circuit and export read the BASE.")
  println()

  # 3) short circuit reads the base, so the compensated net matches the
  # equipment net bus for bus (the fault bypasses the FACTS converter)
  comp_sc = runShortCircuit!(net, scd; case = :max)
  worst = 0.0
  for (a, b) in zip(sort(base_sc.rows; by = r -> r.bus_idx), sort(comp_sc.rows; by = r -> r.bus_idx))
    (isnan(a.ik_kA) || isnan(b.ik_kA)) && continue
    worst = max(worst, abs(a.ik_kA - b.ik_kA))
  end
  println("short circuit (Ik'' max), equipment net vs UPFC net: worst per-bus deviation ",
          round(worst; digits = 9), " kA  -> identical (base impedance used)")
  println()

  # 4) MATPOWER export writes the base impedance, not the negative-r operating
  # point; both nets export without the #326 guard firing
  mktempdir() do d
    runpf!(base, 30, 1e-8, 0); calcNetLosses!(base)
    calcNetLosses!(net)
    writeMatpowerCasefile(base, joinpath(d, "equipment.m"))
    writeMatpowerCasefile(net, joinpath(d, "with_upfc.m"))
    println("MATPOWER export: both cases written; the I->J row carries r = ",
            round(brij.r_base_pu; digits = 5), " pu (base), not the negative live value.")
  end
  println()

  # 5) reset the net to its equipment model, then drop the controller
  restoreBaseImpedances!(net)
  println("restoreBaseImpedances!: live r/x now ", round(brij.r_pu; digits = 5), " / ",
          round(brij.x_pu; digits = 5), " pu (back to base)")
  clearUpfcFullControllers!(net)
  println("clearUpfcFullControllers!: ", length(Sparlectra._upfc_full_controllers(net)),
          " full-UPFC controllers left on the net.")
  return nothing
end

run_example(main)
