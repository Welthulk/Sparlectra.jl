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
# file: benchmark/gpu/cpu_batch_scaling.jl
# purpose: GPU study, stage 2: how the CPU injection batch (C1 or C2)
#          scales with the Julia thread count. One process per thread
#          count (julia -t <k>); N random scenarios, BLAS at 1 thread
#          inside the batch; per run the wall time, the GC time, the time
#          of the serial worker creation and the median runpf! time of one
#          scenario. Median of REPEATS runs after one warm-up.
#
# Run: julia -t <k> --project=benchmark/gpu benchmark/gpu/cpu_batch_scaling.jl [case]
#   default case_ACTIVSg25k; N (default 48), REPEATS (default 3),
#   C2=1 for power mode with power_mode_lu auto. One CSV line on stdout.

using Sparlectra, KLU, LinearAlgebra, Printf, Statistics

const QUIET = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)

function main()
  name = isempty(ARGS) ? "case_ACTIVSg25k" : ARGS[1]
  N = parse(Int, get(ENV, "N", "48"))
  repeats = parse(Int, get(ENV, "REPEATS", "3"))
  c2 = get(ENV, "C2", "0") == "1"
  p = Sparlectra.ensure_casefile(name * ".m"; to_jl = false)
  net = Sparlectra.import_case(p, SparlectraConfig(output = QUIET)).net
  eff = Sparlectra._runpf_effective_config(net, nothing)
  pf = eff === nothing ? Sparlectra.powerflow_config() : (eff isa SparlectraConfig ? eff.powerflow : eff)
  c2 && (pf = Sparlectra._copy_powerflow_with(pf; power_mode = true, power_mode_lu = :auto))
  sc = random_load_scenarios(net, N; seed = 1)
  rows = NamedTuple[]
  for k = 0:repeats
    GC.gc()
    t = @timed runInjectionBatch(net, sc; config = pf, blas_threads = 1, executor = :threads, parallel_enabled = true, parallel_max_tasks = Threads.nthreads(), parallel_min_work_items = 1)
    st = t.value.stats
    k > 0 && push!(rows, (time = t.time, gc = t.gctime, copy = st[:worker_copy_s], runpf = median(st[:cpu_runpf_s])))
  end
  med(f) = median([getfield(r, f) for r in rows])
  println("case,variant,threads,N,repeats,median_s,gc_s,worker_copy_s,runpf_median_s,scen_per_s,scen_per_s_without_copy")
  @printf("%s,%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f\n", name, c2 ? "C2" : "C1", Threads.nthreads(), N, repeats, med(:time), med(:gc), med(:copy), med(:runpf), N / med(:time), N / (med(:time) - med(:copy)))
end

Base.invokelatest(main)
