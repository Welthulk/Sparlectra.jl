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

# file: examples/others/exp_cgmes_topology_processor.jl
# purpose: node-breaker import WITHOUT a TP profile (#314): the topology
#          processor derives the bus partition from ConnectivityNodes and
#          switch states; compared bus for bus against the shipped TP.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Offline-safe: the ENTSO-E MiniGrid node-breaker fixture is used from the
# local test-set cache only. Without the cache the example explains how to
# fetch it and exits cleanly instead of downloading on its own.
function _minigrid_nb_paths()
  cache = Sparlectra.CGMESImporter.cgmesTestSetCacheDir()
  nb = joinpath(cache, "extracted", "MiniGrid", "NodeBreaker")
  base = joinpath(nb, "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3")
  bd = joinpath(nb, "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
  return isdir(base) && isdir(bd) ? (base, bd) : nothing
end

"""
    main()

Import the MiniGrid node-breaker delivery twice: once complete (the shipped
TP assigns the buses) and once with the TP and SV files withheld, the
EQ+SSH-only shape real EMS exports ship. The second import triggers the
topology processor, which aggregates connectivity nodes across closed
non-retained switches into derived buses. Both nets are solved from a flat
start and compared: same partition size, same sorted voltage magnitudes.
"""
function main()
  print_example_banner("examples/others/exp_cgmes_topology_processor.jl", "CGMES node-breaker import without a TP profile: derived buses vs. the shipped TP")
  paths = _minigrid_nb_paths()
  if paths === nothing
    println("MiniGrid node-breaker fixture not cached, run examples/experimental/cgmes_fetch_testsets.jl once, then retry.")
    return nothing
  end
  base, bd = paths

  # both sides load the same profile subset apart from TP: no SV, so tap
  # positions and operating point come from SSH in both imports
  no_sv = [f for f in readdir(base; join = true) if endswith(f, ".xml") && !occursin("_SV", basename(f))]
  no_tp = [f for f in no_sv if !occursin("_TP", basename(f))]

  ref = importCGMES(path = vcat(no_sv, bd), name = "minigrid_tp")
  derived = importCGMES(path = vcat(no_tp, bd), name = "minigrid_no_tp")
  for m in derived.messages
    occursin("topology processor", m) && println(m)
  end
  println("with shipped TP: ", length(ref.net.nodeVec), " buses, ", length(ref.net.branchVec), " branches")
  println("derived (no TP): ", length(derived.net.nodeVec), " buses, ", length(derived.net.branchVec), " branches")

  ref.net.flatstart = true
  derived.net.flatstart = true
  _, e1 = runpf!(ref.net, 60, 1e-8, 0; islands_enabled = true)
  _, e2 = runpf!(derived.net, 60, 1e-8, 0; islands_enabled = true)
  vm_ref = sort([n._vm_pu for n in ref.net.nodeVec])
  vm_der = sort([n._vm_pu for n in derived.net.nodeVec])
  dvm = maximum(abs.(vm_ref .- vm_der))
  println("both solve (status ", e1, "/", e2, "); max |Vm| difference across the sorted buses: ", dvm)
  return (buses_equal = length(ref.net.nodeVec) == length(derived.net.nodeVec), solved = e1 == 0 && e2 == 0, max_dvm = dvm)
end

result = run_example(main)
if result !== nothing
  println()
  println("partition identical: ", result.buses_equal, ", both solved: ", result.solved, ", max dVm: ", result.max_dvm)
end
