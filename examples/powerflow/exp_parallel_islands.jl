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
# file: examples/powerflow/exp_parallel_islands.jl
# purpose: feel-the-progress demo for the parallel island path (Phase 2 of
#          the multi-core work): the same multi-island net solved with
#          power_flow.islands.mode = solve_independent (serial) and
#          solve_parallel, wall clock side by side plus a voltage-identity
#          check. Run with threads to see the effect:
#          julia --threads=auto --project=. examples/powerflow/exp_parallel_islands.jl

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# One island = a meshed ring of `m` buses fed by one grid connection, with
# distributed loads. `n` disconnected copies form the multi-island net; the
# islands differ slightly in load so their solves are not trivially equal.
function build_island_net(n::Int, m::Int)
  net = Net(name = "parallel_islands_$(n)x$(m)", baseMVA = 100.0)
  for k in 1:n
    name = i -> "I$(k)_B$(i)"
    for i in 1:m
      addBus!(net = net, busName = name(i), vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = name(1), type = "EXTERNALNETWORKINJECTION", referencePri = name(1), vm_pu = 1.0, va_deg = 0.0)
    for i in 1:m
      addPIModelACLine!(net = net, fromBus = name(i), toBus = name(i == m ? 1 : i + 1), r_pu = 0.0005, x_pu = 0.002, b_pu = 0.0, status = 1)
      if i % 4 == 0
        addProsumer!(net = net, busName = name(i), type = "ENERGYCONSUMER", p = 0.2 + 0.01 * k, q = 0.05)
      end
    end
  end
  ok, msg = validate!(net = net)
  ok || error("validation failed: $msg")
  return net
end

function solve_islands(mode::Symbol, n::Int, m::Int)
  net = build_island_net(n, m)
  profile = Dict{Symbol,Any}(:enabled => true)
  etime = @elapsed begin
    ite, erg = runpf!(net, 30, 1e-8, 0; islands_enabled = true, islands_mode = mode, islands_parallel_min_work_items = 2, performance_profile = profile)
    erg == 0 || error("$(mode) solve did not converge (erg=$(erg))")
  end
  return (; net, profile, etime)
end

function main()
  print_example_banner("examples/powerflow/exp_parallel_islands.jl", "serial vs parallel island solving on one multi-island net (Phase 2 of the multi-core work)")

  n_islands = 8
  buses_per_island = 600
  println("Julia threads : ", Threads.nthreads(), Threads.nthreads() == 1 ? "  <- start with julia --threads=auto to see the parallel effect" : "")
  println("net           : ", n_islands, " disconnected islands, ", buses_per_island, " buses each")
  println()

  # warm both paths on a tiny copy so compilation stays out of the timings
  solve_islands(:solve_independent, 2, 12)
  solve_islands(:solve_parallel, 2, 12)

  serial = solve_islands(:solve_independent, n_islands, buses_per_island)
  parallel = solve_islands(:solve_parallel, n_islands, buses_per_island)

  identical = all(serial.net.nodeVec[i]._vm_pu == parallel.net.nodeVec[i]._vm_pu && serial.net.nodeVec[i]._va_deg == parallel.net.nodeVec[i]._va_deg for i in eachindex(serial.net.nodeVec))

  println("serial   (solve_independent): ", round(serial.etime; digits = 3), " s")
  println("parallel (solve_parallel)   : ", round(parallel.etime; digits = 3), " s")
  if haskey(parallel.profile[:timings], :parallel_wall_time)
    wall = parallel.profile[:timings][:parallel_wall_time].elapsed_s
    println("parallel fan-out wall clock : ", round(wall; digits = 3), " s (island solves only, without island detection and merge)")
  end
  println("speedup                     : ", round(serial.etime / parallel.etime; digits = 2), "x")
  println("voltages bitwise identical  : ", identical)
  identical || error("serial and parallel results differ")
  return nothing
end

run_example(main)
