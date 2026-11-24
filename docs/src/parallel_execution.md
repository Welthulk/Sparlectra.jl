# Parallel Execution

Sparlectra runs independent work items on Julia threads at four sites,
the batches of two of them also in worker processes; everything else is
serial on purpose.

## What runs in parallel

* **Island solves**: a net that splits into several synchronous AC islands
  solves every island of at least `power_flow.islands.parallel_min_buses`
  (default 200) buses on its own task, largest first, whenever two or more
  such islands exist; smaller islands follow serially. There is no mode to
  choose: with one thread, or below the threshold, the solve is serial.
* **State-estimation islands**: the per-island estimation follows the same
  threshold, every island on its own subnet copy; the results are merged in
  island order.
* **Short-circuit sweeps** (`runShortCircuit!` over many fault buses) run
  the per-bus evaluations in chunks. The selected-inverse sweep
  (`sweep_method = :takahashi`, [Short Circuit](short_circuit.md)) is an
  algorithmic speedup on top and itself serial.
* **N-1 contingency and scenario batches** (`runContingencies!`,
  `runScenarios!`): the outage cases are independent power flows and run in
  chunks on threads or in worker processes (below).
* **Injection batches** (`runInjectionBatch`, [Injection batch](@ref contingency_injection_batch)):
  the load scenarios of one network run on reused working copies, in worker
  processes or in chunks on threads, with BLAS at one thread inside the
  batch.

The Web UI's background task per solver job is no computation fan-out. The
Newton iteration of one island does not parallelize (one sparse
factorization per iteration): a single large connected case does not get
faster with more threads, many independent solves do.

## How the CPUs are used

```bash
julia -t auto --project=. yourscript.jl        # or JULIA_NUM_THREADS=8
```

`runtime.parallel.enabled` is the master switch, `max_tasks` caps the
concurrent tasks per site by chunking (100 outages on 8 tasks run as 8
chunks), and work lists shorter than `min_work_items` run serially; keys
and defaults in [Runtime](@ref perf-runtime). With one Julia thread
everything runs serially regardless of the switches; the startup summary
(`runtime.print_thread_config`) prints the resolved
`parallel: enabled=... max_tasks=...` line. BLAS keeps its own thread pool;
the sparse factorizations do not profit from raising it.

```yaml
runtime:
  parallel:
    enabled: true
    max_tasks: auto
    min_work_items: 4
power_flow:
  islands:
    parallel_min_buses: 200
```

## [Batches in worker processes](@id parallel_execution_processes)

Threads of one process do not scale for large networks: the sparse LU of
every solve takes a lock of the BLAS library many times per factorization,
and every garbage collection waits for all threads to leave their current
factorization. Separate processes have their own of both. The batches of the
scenario engine (N-1, scenarios) and the injection batch therefore run in
worker processes on large networks.

`runtime.parallel.executor` (or the keyword `executor` of
`runContingencies!`, `runScenarios!` and `runInjectionBatch`) picks:

| Value | Batch runs |
|---|---|
| `auto` (default) | in processes for networks from 1000 buses on when the worker pool already holds the workers or the batch is large enough to pay the pool start (items times buses at least 3e6); on threads otherwise |
| `processes` | in the worker processes of the batch pool |
| `threads` | on the Julia threads of this process |

`runtime.parallel.enabled = false`, or fewer items than `min_work_items`,
run serially whatever the executor.

The pool holds Julia worker processes of the `Distributed` standard
library, one thread each, started with the project, load path and system
image of the main process. A batch that runs in processes starts the pool
when it is not there and keeps it for the next batch;
`start_batch_processes!(n)` starts it ahead (the first batch then does not
pay for the start and the compilation on the workers),
`stop_batch_processes!()` ends it. `runtime.parallel.max_processes` caps the
workers (`auto`: `Sys.CPU_THREADS`, fewer when the free memory holds
fewer workers of the size of the network). A batch sends the network once per
worker; the scenarios go to the next free worker one by one. A worker keeps
the network of its last batch, so the next batch on the same network (same
content) sends nothing (`stats[:workers_reused]` of an injection batch).
The abort of a Web UI job reaches the scenarios that run on the workers:
each checks it once per Newton iteration, as on threads, and the batch
throws `PowerFlowAborted`. The progress counter runs in the main process.

```julia
start_batch_processes!(8)
r = runInjectionBatch(net, random_load_scenarios(net, 64); executor = :processes)
r.stats[:executor], r.stats[:workers]     # (:processes, 8)
```

Each worker needs its own memory: the loaded package, the network twice
(template and working copy) and the factorizations of one solve; cap the
workers with `max_processes` on a large network. Each worker also compiles
the code paths that are not in the package image on their first use, so
the first batches on a new pool are slower than the later ones. Worker
processes start without the power-mode cache and the solver status of the
network; with `power_flow.power_mode` on, a worker of an injection batch
solves the base case once when it loads the batch, as a chunk on threads
does before its first scenario. When the main process has KLU loaded, the
workers load it too, so power mode factorizes with the same solver there.

## [GC pause](@id parallel_execution_gc_pause)

The garbage collector is off while a power flow runs (`runpf!`, and every
solve of a batch) and on again when the last run of the process has ended,
also after an error. Runs that overlap (threads, Web UI jobs) count
against one process-wide counter. A run allocates several times the size of
its network without collecting; the memory guard checks the resident size
once per Newton iteration and collects once when it passes
`runtime.gc_pause.max_resident_gb` (`auto`: half of the physical memory;
every worker process of the pool gets its share). `with_gc_paused(f)` runs
any function as such a run, `gc_pause_status()` shows the counter, the limit
and the collections of the guard.

```yaml
runtime:
  gc_pause:
    enabled: true
    max_resident_gb: auto
```

## Guarantees

Every parallel site gives bitwise the same result as its serial path (the
same function, not a copy), tested in repeated identity runs; no
factorization is shared across tasks. The process batch gives the result of
the serial batch bit for bit as well. The one difference is the failure
semantics: in the parallel island mode every island reports its real status
and a failure is raised after all islands finish, the serial mode stops at
the first failure ([PowerFlow Configuration](powerflow_configuration.md)).

## Measuring the effect

Phase times are summed CPU-side work; `parallel_wall_time` records the
wall-clock time of the fan-out, and the speedup shows there, not in the
phase sums ([Performance and Profiling](performance_profiling.md)).
