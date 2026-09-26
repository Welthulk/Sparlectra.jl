# Parallel Execution

Sparlectra runs independent work items on Julia threads. Three computation
sites fan out; everything else is serial on purpose.

## What runs in parallel

* **Island solves** (`power_flow.islands.mode = solve_parallel`): a net
  that splits into several synchronous AC islands solves them
  concurrently, largest island first.
* **Short-circuit sweeps** (`runShortCircuit!` over many fault buses): the
  per-bus evaluations run in chunks. The selected-inverse sweep
  (`sweep_method = :takahashi`, see [Short Circuit](short_circuit.md)) is
  an algorithmic speedup on top and itself serial.
* **N-1 contingency batches** (`runContingencies!`): the outage cases are
  independent power flows and run in chunks.

The Web UI runs every solver job on a background task for responsiveness;
that is not a computation fan-out.

The Newton iteration of a single island does not parallelize (one sparse
factorization per iteration). A single large connected case does not get
faster with more threads; the payoff is in many independent solves.

## How the CPUs are used

```bash
julia -t auto --project=. yourscript.jl        # or JULIA_NUM_THREADS=8
```

* `runtime.parallel.enabled` (default `true`) is the master switch. Off
  means every site takes its serial path (the same function, not a copy).
* `runtime.parallel.max_tasks` (default `auto` = `Threads.nthreads()`)
  caps the concurrent tasks per site by chunking: 100 outages on 8 tasks
  run as 8 chunks.
* `runtime.parallel.min_work_items` (default `4`): shorter work lists run
  serially.
* With one Julia thread everything runs serially regardless of the
  switches; the startup summary (`runtime.print_thread_config`) prints
  the resolved `parallel: enabled=... max_tasks=...` line.
* BLAS keeps its own thread pool (printed in the startup summary); the
  sparse factorizations do not profit from raising it.

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

* **Bitwise identical results.** Every parallel site is tested against
  the serial result in repeated identity runs.
* **No shared factorizations across tasks.** Each chunk works on its own
  copy of the linear-solver factorization, created serially before the
  tasks start; a shared one would serialize or silently corrupt the
  solves.
* **Worker state is indexed by chunk, never by thread id.**
* **Failure semantics.** In the parallel island mode every island reports
  its real status and a failure is raised after all islands finish; the
  serial mode stops immediately. This is the one semantic difference, see
  [PowerFlow Configuration](powerflow_configuration.md).

## Measuring the effect

Phase times are summed CPU-side work; the extra `parallel_wall_time` entry
records the wall-clock time of the fan-out. Speedup shows in the wall
clock, not in the phase sums. Option reference:
[Performance and Profiling](performance_profiling.md).
