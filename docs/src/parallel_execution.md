# Parallel Execution

Sparlectra runs independent work items on Julia threads. This page says
what actually runs in parallel, how the CPUs are used, and which
guarantees hold, so you can decide how many threads to give the process
and what to expect from them.

## What runs in parallel

Three computation sites fan out over threads; everything else is serial on
purpose.

* **Island solves** (`power_flow.islands.mode = solve_parallel`): a net
  that splits into several synchronous AC islands solves them
  concurrently, largest island first. Multi-area MATPOWER cases and CGMES
  deliveries with stub islands profit directly.
* **Short-circuit sweeps** (`runShortCircuit!` over many fault buses): the
  per-bus evaluations are independent and run in chunks. The optional
  selected-inverse sweep (`sweep_method = :takahashi`, see
  [Short Circuit](short_circuit.md)) is an ALGORITHMIC speedup on top and
  itself serial.
* **N-1 contingency batches** (`runContingencies!`): the outage cases are
  independent power flows and run in chunks.

The Web UI additionally runs every solver job on a background task so the
pages stay responsive; that is concurrency for responsiveness, not a
computation fan-out.

What does NOT parallelize: the Newton iteration of a single island. Its
dominant cost is one sparse factorization per iteration, and that
factorization is serial per island. A single large connected case
therefore does not get faster with more threads; the payoff of threads is
in MANY independent solves (islands, fault buses, outages).

## How the CPUs are used

Start Julia with threads, for example:

```bash
julia -t auto --project=. yourscript.jl        # or JULIA_NUM_THREADS=8
```

* `runtime.parallel.enabled` (default `true`) is the master switch. Off
  means every site takes its serial path, which is the SAME function, not
  a copy.
* `runtime.parallel.max_tasks` (default `auto` = `Threads.nthreads()`)
  caps the concurrent tasks per site. The cap works by chunking the work
  list: 100 outages on 8 tasks run as 8 chunks, not 100 tasks.
* `runtime.parallel.min_work_items` (default `4`): shorter work lists run
  serially, the task overhead would exceed the gain.
* With one Julia thread everything runs serially regardless of the
  switches; the startup summary (`runtime.print_thread_config`) prints the
  resolved `parallel: enabled=... max_tasks=...` line so a run's report
  says what it actually used.
* BLAS keeps its own thread pool; the startup summary prints it. The
  sparse factorizations dominate and do not profit from raising it.

Configuration example:

```yaml
runtime:
  parallel:
    enabled: true
    max_tasks: auto
    min_work_items: 4
power_flow:
  islands:
    mode: solve_parallel
```

## Guarantees

* **Bitwise identical results.** Every parallel site is tested to produce
  exactly the serial result, and the identity tests run REPEATED (a
  wrong sharing can pass a single comparison by luck).
* **No shared factorizations across tasks.** Each chunk works on its own
  copy of the linear-solver factorization; the copies are created
  serially before the tasks start. A shared factorization would either
  serialize the solves (no speedup) or corrupt them silently, so it is
  never used concurrently.
* **Worker state is indexed by chunk, never by thread id.** Thread ids
  can exceed the thread count under an interactive threadpool; chunk
  indexing is stable.
* **Failure semantics.** In the parallel island mode every island reports
  its real status and a failure is raised after all islands finish
  (in-flight islands cannot be skipped); the serial mode keeps the
  immediate stop. This is the one documented semantic difference, see
  [PowerFlow Configuration](powerflow_configuration.md).

## Measuring the effect

The per-phase timings keep their meaning in parallel runs: phase times are
summed CPU-side work, and the extra `parallel_wall_time` entry records the
wall-clock time of the fan-out. Speedup shows in the wall clock, not in
the phase sums. The option reference lives in
[Performance and Profiling](performance_profiling.md).
