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
# file: examples/others/exp_short_circuit_cgmes.jl
# purpose: runShortCircuit! on a real CGMES delivery (ENTSO-E MicroGrid BE) — Ik'' max/min per bus straight from the import's short-circuit harvest

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Offline-safe: the ENTSO-E MicroGrid fixture is used from the local test-set
# cache only. Without the cache the example explains how to fetch it and
# exits cleanly instead of downloading on its own.
function _microgrid_paths()
  cache = Sparlectra.CGMESImporter.cgmesTestSetCacheDir()
  bc = joinpath(cache, "extracted", "MicroGrid", "BaseCase_BC")
  be = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BC_BE_v2")
  bd = joinpath(bc, "CGMES_v2.4.15_MicroGridTestConfiguration_BD_v2")
  return isdir(be) && isdir(bd) ? (be, bd) : nothing
end

function main()
  print_example_banner("examples/others/exp_short_circuit_cgmes.jl", "runShortCircuit! on the ENTSO-E MicroGrid BE delivery — Ik'' max/min per bus from the import's short-circuit harvest")
  paths = _microgrid_paths()
  if paths === nothing
    println("MicroGrid fixture not cached — run examples/experimental/cgmes_fetch_testsets.jl once, then retry.")
    return nothing
  end
  be, bd = paths

  result = importCGMES(path = [be, bd], name = "MicroGrid_BE")
  println("Harvested short-circuit source data (read on every CGMES import):")
  printShortCircuitCoverage(stdout, result.shortcircuit)
  println()

  sc_max = runShortCircuit!(result; case = :max)
  sc_min = runShortCircuit!(result; case = :min)
  println("Maximum case (equipment rating):")
  printShortCircuitResult(sc_max)
  println()
  println("Minimum case (protection sensitivity):")
  printShortCircuitResult(sc_min)
  println()

  # Single-substation question ("this busbar"): the per-bus API is the
  # primary contract; the sweep above is just a loop over it.
  bus = sc_max.rows[1].bus
  single = runShortCircuit!(result.net, result.shortcircuit; buses = [bus], case = :max)
  row = single.rows[1]
  println("Single-bus query for $(bus): Ik'' = $(round(row.ik_kA; sigdigits = 5)) kA, Sk'' = $(round(row.sk_MVA; sigdigits = 5)) MVA, ip = $(round(row.ip_kA; sigdigits = 5)) kA")
  return sc_max
end

run_example(main)
