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
#
# file: examples/others/exp_process_batch.jl
# purpose: batches in worker processes on the shipped sp_case1354: an
#          injection batch of 32 load scenarios and the branch N-1 batch,
#          each on threads and in worker processes, with the results
#          compared bit for bit, and the GC pause status of a single run.
#          The pool of worker processes is started ahead and stopped at
#          the end.
#          julia -t 4 --project=. examples/others/exp_process_batch.jl

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const CASE_PATH = joinpath(Sparlectra.MPOWER_DIR, "sp_case1354.m")
const WORKERS = 4

function main()
  print_example_banner("examples/others/exp_process_batch.jl", "injection and N-1 batches on threads and in worker processes, sp_case1354")
  net = createNetFromMatPowerFile(filename = CASE_PATH)
  println("network: ", length(net.nodeVec), " buses, Julia threads: ", Threads.nthreads())

  # the pool start (a few seconds) and the compilation on the workers stay
  # out of the timed batches: start the pool, then one batch as warm-up
  t = @elapsed start_batch_processes!(WORKERS)
  @printf("pool of %d worker processes started in %.1f s: %s\n", WORKERS, t, join(batch_process_ids(), ", "))

  sc = random_load_scenarios(net, 32; band = 0.2, seed = 1)
  cases = generateN1Branches(net)
  try
    for ex in (:threads, :processes)
      runInjectionBatch(net, sc; executor = ex, max_processes = WORKERS, parallel_max_tasks = WORKERS, parallel_min_work_items = 1)
      runContingencies!(net, cases; executor = ex, max_processes = WORKERS, parallel_max_tasks = WORKERS, parallel_min_work_items = 1)
    end
    inj = Dict{Symbol,Any}()
    n1 = Dict{Symbol,Any}()
    for ex in (:threads, :processes)
      t_inj = @elapsed inj[ex] = runInjectionBatch(net, sc; executor = ex, max_processes = WORKERS, parallel_max_tasks = WORKERS, parallel_min_work_items = 1)
      t_n1 = @elapsed n1[ex] = runContingencies!(net, cases; executor = ex, max_processes = WORKERS, parallel_max_tasks = WORKERS, parallel_min_work_items = 1)
      r = inj[ex]
      @printf("%-9s injection batch: ran as %s on %d, %d of %d converged, %.2f s; N-1 of %d branches: %.2f s\n",
        String(ex), String(r.stats[:executor]), r.stats[:workers], count(r.converged), length(r.converged), t_inj, length(cases), t_n1)
    end
    a, b = inj[:threads], inj[:processes]
    println("injection batch, threads against processes: same voltages bit for bit ", a.vm == b.vm && a.va_deg == b.va_deg, ", same iterations ", a.iterations == b.iterations)
    same_n1 = all(i -> n1[:threads][i].converged == n1[:processes][i].converged && n1[:threads][i].iterations == n1[:processes][i].iterations, eachindex(cases))
    println("N-1, threads against processes: same convergence and iterations per case ", same_n1)
  finally
    stop_batch_processes!()
  end

  # every power flow runs with the garbage collector paused; with_gc_paused
  # does the same for any function. Inside, the status counts the run and
  # shows the limit of the memory guard; afterwards the GC is on again.
  w = deepcopy(net)
  inside = with_gc_paused(() -> (runpf!(w); gc_pause_status()))
  after = gc_pause_status()
  @printf("GC pause: inside the run paused %s, %d active run(s), guard limit %.1f GB; afterwards paused %s, %d active\n",
    inside.paused, inside.active_runs, inside.limit_bytes / 2^30, after.paused, after.active_runs)
  return nothing
end

run_example(main)
