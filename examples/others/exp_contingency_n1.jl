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
# file: examples/others/exp_contingency_n1.jl
# purpose: N-1 showcase (Phase 4 of the multi-core work): full branch N-1 on
#          case1354pegase, serial vs parallel timing side by side, top 10
#          worst contingencies, plus the deepcopy cost and warm-start
#          iteration numbers the analysis report records. Needs the
#          case1354pegase.m from the Web UI case cache; skips with a message
#          when it is not cached. Run with threads to see the effect:
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

  # PEGASE convention: branch angles are radians with inverted sign (the
  # auto-profile residual scan identifies this; recorded in the project
  # notes). qlimits stay on: case1354 converges with the active set.
  net = createNetFromMatPowerFile(filename = CASE_PATH, matpower_shift_unit = :rad, matpower_shift_sign = -1.0)
  cases = generateN1Branches(net)
  println("Julia threads : ", Threads.nthreads(), Threads.nthreads() == 1 ? "  <- start with julia --threads=auto to see the parallel effect" : "")
  println("case          : case1354pegase, ", length(net.nodeVec), " buses, ", length(net.branchVec), " branches -> ", length(cases), " contingencies")

  t_copy = @elapsed deepcopy(net)
  println("deepcopy(net) : ", round(t_copy * 1000; digits = 1), " ms per contingency worker (structural-sharing copy is the known follow-up)")
  println()

  # warm start evidence for the analysis report: the template is the solved
  # base case, so contingency solves start from the base voltages. The flat
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

  # top 10 worst: by loading when ratings exist, else by violation count
  has_loading = any(r -> !isnan(r.max_branch_loading_pct), serial)
  worst = sort([r for r in serial if r.converged]; by = r -> has_loading ? -(isnan(r.max_branch_loading_pct) ? -Inf : r.max_branch_loading_pct) : -(length(r.voltage_violations)))
  println("top 10 worst contingencies (by ", has_loading ? "max branch loading" : "voltage-violation count", "):")
  printContingencyResults(stdout, worst[1:min(10, end)]; max_rows = 10)
  return nothing
end

run_example(main)
