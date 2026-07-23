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

# Date: 2026-07-23
# file: examples/exp_dc_powerflow.jl
# purpose: solves a standalone DC power flow (MATPOWER rundcpf equivalent) and optionally seeds the AC Newton-Raphson solve from it

using Sparlectra

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function main(; casefile::AbstractString = "case9.m")
  print_example_banner("examples/exp_dc_powerflow.jl", "solves a standalone DC power flow (MATPOWER rundcpf equivalent) and optionally seeds the AC Newton-Raphson solve from it")
  case_path = ensure_casefile(casefile)

  println("=== Standalone DC power flow ===")
  net = createNetFromMatPowerFile(filename = case_path)
  report = rundcpf!(net; verbose = 1)
  println()
  println("Bus angles (deg):")
  for row in report.nodes
    println("  ", row.bus, "  ", rpad(row.bus_name, 12), "  va=", round(row.va_deg; digits = 4), "  Pgen=", round(row.p_gen_MW; digits = 2), " MW  Pload=", round(row.p_load_MW; digits = 2), " MW")
  end
  println()
  println("Branch flows (lossless: p_to_MW == -p_from_MW):")
  for row in report.branches
    println("  ", row.branch, "  ", row.from_bus, "->", row.to_bus, "  Pf=", round(row.p_from_MW; digits = 2), " MW")
  end

  println()
  println("=== Same case, DC-seeded AC Newton-Raphson (seed_ac_start=true) ===")
  net2 = createNetFromMatPowerFile(filename = case_path)
  seeded_report = rundcpf!(net2; verbose = 1, seed_ac_start = true)
  println()
  println("DC step converged   : ", seeded_report.metadata.converged)
  println("AC step converged   : ", seeded_report.metadata.ac_converged)
  println("AC iterations       : ", seeded_report.metadata.ac_iterations)
  println("net now holds the AC-converged solution (voltage magnitudes are no longer flattened to 1.0 pu):")
  for n in sort(net2.nodeVec, by = x -> x.busIdx)
    println("  bus ", n.busIdx, "  vm=", round(n._vm_pu; digits = 4), " pu  va=", round(n._va_deg; digits = 4), " deg")
  end

  return (dc_report = report, seeded_report = seeded_report)
end

run_example(main)
