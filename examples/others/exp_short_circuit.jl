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

# Date: 2026-07-30
# file: examples/others/exp_short_circuit.jl
# purpose: demonstrates runShortCircuit! (IEC 60909-0, issue #277) — feeder/machine sources, Ik'' max vs. min, and the safety flag on defaulted data

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Hand-built two-bus network mirroring the CGMES harvest shape: a 110 kV
# feeder (declared 10 kA / R/X 0.1) feeding a line to a second bus. The same
# data arrives from a real delivery via importCGMES(...).shortcircuit — see
# docs/src/short_circuit.md.
function _demo_net_and_data()
  net = Net(name = "sc_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "STATION", vn_kV = 110.0)
  addBus!(net = net, busName = "FEEDER_POINT", vn_kV = 110.0)
  addPIModelACLine!(net = net, fromBus = "FEEDER_POINT", toBus = "STATION", r_pu = 0.05, x_pu = 0.5, b_pu = 0.1, status = 1)
  feeder = (mrid = "F1", name = "Feeder", bus = "FEEDER_POINT", maxInitialSymShCCurrent_A = 10_000.0, minInitialSymShCCurrent_A = 8_000.0, maxR1ToX1Ratio = 0.1, minR1ToX1Ratio = 0.1, maxR0ToX0Ratio = nothing, maxZ0ToZ1Ratio = nothing, ikSecond = nothing, governorSCD = nothing)
  # A machine WITHOUT x''_d: the documented default (0.2 pu) is substituted
  # and — because Ik'' is safety-relevant — flagged on the result rows.
  machine = (mrid = "G1", name = "G1", bus = "STATION", satDirectSubtransX_pu = nothing, satDirectTransX_pu = nothing, r0_pu = nothing, x0_pu = nothing, r2_pu = nothing, x2_pu = nothing, earthing = nothing, ratedS_MVA = 50.0, ratedU_kV = 110.0)
  sc = Sparlectra.CGMESImporter.CGMESShortCircuitData([feeder], [machine], NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[])
  return net, sc
end

function main()
  print_example_banner("examples/others/exp_short_circuit.jl", "demonstrates runShortCircuit! (IEC 60909-0) — feeder/machine sources, Ik'' max vs. min, and the safety flag on defaulted data")
  net, sc = _demo_net_and_data()

  result_max = runShortCircuit!(net, sc; case = :max)
  result_min = runShortCircuit!(net, sc; case = :min)

  println("Maximum case (equipment rating, c_max per IEC 60909-0 Table 1):")
  printShortCircuitResult(result_max)
  println()
  println("Minimum case (protection sensitivity, c_min):")
  printShortCircuitResult(result_min)
  println()
  println("Note the flag column: the machine ships without x''_d, so the")
  println("documented default was substituted and every affected row says so —")
  println("for the :max case such a result is a lower bound.")
  println()
  println("On a CGMES delivery the same call is:")
  println("  result = importCGMES(path = [\"grid.zip\", \"boundary.zip\"])")
  println("  runShortCircuit!(result; case = :max)")
  return result_max
end

run_example(main)
