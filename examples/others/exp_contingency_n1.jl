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
# file: examples/others/exp_contingency_n1.jl
# purpose: N-1 showcase on case1354pegase: case screening (voltage/rating
#          filters) and per-case weights (#331 Phase 2), generator outages with
#          auto/distributed slack (#331 Phase 3), then full branch N-1 serial vs
#          parallel timing side by side, top 10 worst contingencies, plus the
#          deepcopy cost and warm-start iteration numbers it prints. Needs the
#          case1354pegase.m from the Web UI case cache; skips
#          with a message when it is not cached. Run with threads:
#          julia --threads=auto --project=. examples/others/exp_contingency_n1.jl

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

const CASE_PATH = joinpath(homedir(), ".local", "state", "sparlectra", "webui", "data", "mpower", "case1354pegase.m")

function main()
  print_example_banner("examples/others/exp_contingency_n1.jl", "full branch N-1 on case1354pegase, serial vs parallel (Phase 4 of the multi-core work)")

  if !isfile(CASE_PATH)
    println("SKIPPED: case1354pegase.m not found under the Web UI case cache (", CASE_PATH, ").")
    println("Fetch it once through the Web UI or the service path, then rerun.")
    return nothing
  end

  # PEGASE convention (from the case sidecar): branch angles are radians with
  # inverted sign, tap ratio normal (the default; the convention with the
  # fewest NR iterations is the correct one, which the auto-profile residual
  # scan confirms). qlimits stay on: case1354 converges with the active set.
  net = createNetFromMatPowerFile(filename = CASE_PATH, matpower_shift_unit = :rad, matpower_shift_sign = -1.0, matpower_ratio = :normal)
  cases = generateN1Branches(net)
  println("Julia threads : ", Threads.nthreads(), Threads.nthreads() == 1 ? "  <- start with julia --threads=auto to see the parallel effect" : "")
  println("case          : case1354pegase, ", length(net.nodeVec), " buses, ", length(net.branchVec), " branches -> ", length(cases), " contingencies")

  # --- case screening and weights (#331 Phase 2) ---
  # On a large grid you rarely simulate every outage: screen the list down to
  # the top voltage level (and, if you like, to rated branches only).
  ehv = generateN1Branches(net; min_vn_kV = 300.0)
  rated = generateN1Branches(net; min_sn_MVA = 1.0)
  println("screening     : ", length(ehv), " outages at >= 300 kV, ", length(rated), " with a branch rating (of ", length(cases), ")")
  # Attach per-branch weights (here the EHV outages count double); in practice
  # read a name->rate map from a CSV with readContingencyWeightsCSV. The weight
  # rides through to the result and shows up in the table's weight column.
  weights = Dict(c.name => 2.0 for c in ehv)
  ehv_weighted = applyContingencyWeights(ehv, weights)
  ehv_res = runContingencies!(net, ehv_weighted[1:min(8, end)]; parallel_enabled = false)
  println("screened EHV subset (weight column carries the per-case weight):")
  printContingencyResults(stdout, ehv_res; max_rows = 8)
  println()

  # --- generator outages (#331 Phase 3) ---
  # One case per generator; auto_slack lets the solver promote a survivor when
  # the outage removes the system slack itself (otherwise that one case reports
  # "no slack bus"). distributed_slack shares the lost output over the rest.
  gcases = generateN1Generators(net)
  gres = runContingencies!(net, gcases; parallel_enabled = true, parallel_min_work_items = 2, auto_slack = true)
  gres_ds = runContingencies!(net, gcases; parallel_enabled = true, parallel_min_work_items = 2, auto_slack = true, distributed_slack_enabled = true)
  println("generator N-1   : ", length(gcases), " units ; ", count(r -> r.converged, gres),
          " converge (auto_slack), ", count(r -> r.converged, gres_ds), " with distributed slack")
  println()

  deepcopy(net)   # first call compiles deepcopy for the Net type tree (~0.5 s once per process)
  t_copy = @elapsed deepcopy(net)
  println("deepcopy(net) : ", round(t_copy * 1000; digits = 1), " ms per contingency worker (warm; about a quarter of one case's cost)")
  println()

  # warm start: the template is the solved base case, so contingency solves
  # start from the base voltages. The flat
  # variant forces flatstart on the template; its base solve then fails (the
  # printed warning is expected) and most flat contingency solves diverge:
  # on this case the warm start is not an iteration saver but the
  # feasibility condition for plain (non-projected) solves.
  warm = runContingencies!(net, cases[1:24]; parallel_enabled = false)
  flat_net = deepcopy(net)
  flat_net.flatstart = true
  flat = runContingencies!(flat_net, cases[1:24]; parallel_enabled = false)
  warm_it = sum(r.iterations for r in warm if r.converged; init = 0)
  flat_it = sum(r.iterations for r in flat if r.converged; init = 0)
  println("warm start    : ", count(r -> r.converged, warm), "/24 converged, ", warm_it, " NR iterations")
  println("flat start    : ", count(r -> r.converged, flat), "/24 converged, ", flat_it, " NR iterations")
  println()

  println("running full N-1 serially ...")
  t_serial = @elapsed serial = runContingencies!(net, cases; parallel_enabled = false)
  println("running full N-1 in parallel ...")
  t_parallel = @elapsed parallel = runContingencies!(net, cases; parallel_enabled = true, parallel_min_work_items = 2)

  identical = all(all(isequal(getfield(serial[i], f), getfield(parallel[i], f)) for f in fieldnames(ContingencyResult)) for i in eachindex(serial))
  println()
  println("serial   : ", round(t_serial; digits = 1), " s")
  println("parallel : ", round(t_parallel; digits = 1), " s")
  println("speedup  : ", round(t_serial / t_parallel; digits = 2), "x on ", Threads.nthreads(), " threads")
  println("results identical: ", identical)
  identical || error("serial and parallel contingency results differ")
  println("converged: ", count(r -> r.converged, serial), " of ", length(serial), ", islanding cases: ", count(r -> r.island_count > 1, serial))
  println()

  # the table ranks by severity (failures first, then the heaviest weighted
  # overloads), so max_rows = 10 already surfaces the worst contingencies
  println("worst 10 contingencies (severity-ranked, failures first):")
  printContingencyResults(stdout, serial; max_rows = 10)
  println()

  # aggregate report (#331 Phase 4): outcome counts, load shed, worst loading,
  # worst severity, and the branches overloaded by the most contingencies
  printContingencyReport(buildContingencyReport(serial))
  return nothing
end

run_example(main)
