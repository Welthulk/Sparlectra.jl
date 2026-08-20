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
# file: examples/others/exp_open_terminal_line.jl
# purpose: one-sided open branches (r0.10.0): a long 380 kV line is solved
#          closed, then opened at the to end with setBranchTerminalStatus!.
#          The open line still draws its FULL charging at the closed bus
#          (plus the small r loss of the charging current), and the open-end
#          voltage shows the Ferranti rise. Theory: branchmodel.md,
#          "One-sided open branches".

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# a long 380 kV line (several hundred km equivalent): substantial charging
# so the effect is not lost in rounding
function build_two_bus(name::String)
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "A", vn_kV = 380.0)
  addBus!(net = net, busName = "B", vn_kV = 380.0)
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 120.0, q = 30.0)
  addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.16, b_pu = 0.9, g_pu = 0.004, status = 1)
  ok, msg = validate!(net = net)
  ok || error("validation failed: $msg")
  return net
end

function main()
  print_example_banner("examples/others/exp_open_terminal_line.jl", "one-sided open line: full charging draw and the Ferranti rise")

  println("=== Both ends closed: the line carries the load ===")
  net = build_two_bus("open_terminal_closed")
  ite, erg = runpf!(net, 30, 1e-10, 0)
  erg == 0 || error("closed solve did not converge")
  calcNetLosses!(net)
  printACPFlowResults(net, 0.0, ite, 1e-10)

  println()
  println("=== Open the to end: no through flow, but the charging stays ===")
  net2 = build_two_bus("open_terminal_open_to")
  setBranchTerminalStatus!(net2.branchVec[1]; to = false)
  markIsolatedBuses!(net = net2, log = false)
  ite2, erg2 = runpf!(net2, 30, 1e-10, 0)
  erg2 == 0 || error("open-terminal solve did not converge")
  calcNetLosses!(net2)
  printACPFlowResults(net2, 0.0, ite2, 1e-10)

  br = net2.branchVec[1]
  println()
  println("Reading the numbers:")
  println("  charging draw at A : ", round(br.fBranchFlow.qFlow; digits = 3), " MVAr (the FULL b of the line, not half)")
  println("  active loss        : ", round(br.fBranchFlow.pFlow; digits = 4), " MW (g losses plus the r loss of the charging current)")
  println("  open-end voltage   : ", round(br.open_end_vm_pu; digits = 5), " pu at ", round(br.open_end_va_deg; digits = 3), " deg, HIGHER than the ", round(get_bus_vm_pu(net2, "A"); digits = 5), " pu at the feeding bus A")
  println("                       (Ferranti effect: the charging current through the line reactance lifts the voltage toward the open end)")
  println("  open terminal      : S = 0 by definition; bus B is reported isolated")

  return (closed = net, open_to = net2)
end

run_example(main)
