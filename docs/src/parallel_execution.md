# Parallel Execution

Sparlectra runs independent work items on Julia threads at four sites;
everything else is serial on purpose.

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
* **N-1 contingency batches** (`runContingencies!`): the outage cases are
  independent power flows and run in chunks.

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

## Guarantees

Every parallel site gives bitwise the same result as its serial path (the
same function, not a copy), tested in repeated identity runs; no
factorization is shared across tasks. The one difference is the failure
semantics: in the parallel island mode every island reports its real status
and a failure is raised after all islands finish, the serial mode stops at
the first failure ([PowerFlow Configuration](powerflow_configuration.md)).

## Measuring the effect

Phase times are summed CPU-side work; `parallel_wall_time` records the
wall-clock time of the fan-out, and the speedup shows there, not in the
phase sums ([Performance and Profiling](performance_profiling.md)).
