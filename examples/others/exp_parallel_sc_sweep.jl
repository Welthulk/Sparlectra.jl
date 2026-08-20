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
# file: examples/others/exp_parallel_sc_sweep.jl
# purpose: feel-the-progress demo for the parallel short-circuit all-bus
#          sweep (Phase 3 of the multi-core work): the same IEC 60909-0
#          sweep run serially and in parallel, wall clock side by side plus
#          a row-identity check. Run with threads to see the effect:
#          julia --threads=auto --project=. examples/others/exp_parallel_sc_sweep.jl

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# feeder-fed rings with declared short-circuit power: every island has one
# ExternalNetworkInjection with Sk'' data, so the sweep has real sources
function build_sweep_net(n::Int, m::Int)
  net = Net(name = "sc_sweep_$(n)x$(m)", baseMVA = 100.0)
  for k in 1:n
    name = i -> "S$(k)_B$(i)"
    for i in 1:m
      addBus!(net = net, busName = name(i), vn_kV = 110.0)
    end
    addExternalGrid!(net = net, busName = name(1), vm_pu = 1.0, sk_max_MVA = 2000.0 + 100.0 * k, sk_min_MVA = 1500.0, rx_max = 0.1, internal_impedance = false)
    for i in 1:m
      addPIModelACLine!(net = net, fromBus = name(i), toBus = name(i == m ? 1 : i + 1), r_pu = 0.001, x_pu = 0.004, b_pu = 0.0, status = 1)
    end
  end
  ok, msg = validate!(net = net)
  ok || error("validation failed: $msg")
  return net
end

function main()
  print_example_banner("examples/others/exp_parallel_sc_sweep.jl", "serial vs parallel IEC 60909-0 all-bus sweep (Phase 3 of the multi-core work)")

  n_islands = 8
  buses_per_island = 1000
  println("Julia threads : ", Threads.nthreads(), Threads.nthreads() == 1 ? "  <- start with julia --threads=auto to see the parallel effect" : "")
  println("net           : ", n_islands, " islands, ", buses_per_island, " buses each = ", n_islands * buses_per_island, " fault locations")
  println()

  net = build_sweep_net(n_islands, buses_per_island)

  # warm both paths (compilation), then measure
  runShortCircuit!(net; case = :max, parallel_enabled = false)
  runShortCircuit!(net; case = :max, parallel_enabled = true, parallel_min_work_items = 2)

  t_serial = @elapsed serial = runShortCircuit!(net; case = :max, parallel_enabled = false)
  t_parallel = @elapsed parallel = runShortCircuit!(net; case = :max, parallel_enabled = true, parallel_min_work_items = 2)

  identical = isequal(serial.rows, parallel.rows)
  println("serial sweep   : ", round(t_serial * 1000; digits = 1), " ms")
  println("parallel sweep : ", round(t_parallel * 1000; digits = 1), " ms")
  println("speedup        : ", round(t_serial / t_parallel; digits = 2), "x")
  println("rows identical : ", identical, " (", length(serial.rows), " rows, NaN-aware)")
  identical || error("serial and parallel sweep results differ")
  return nothing
end

run_example(main)
