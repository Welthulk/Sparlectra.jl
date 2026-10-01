# Performance and Profiling Configuration

## [Runtime](@id perf-runtime)

Process-level settings read once at startup: the thread report, the Julia and BLAS thread counts, and the parallel work split of the sweeps.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `runtime.print_thread_config` | Bool | `true` | `true`, `false` | Print Julia/BLAS thread summary at startup, including the `parallel: enabled=... max_tasks=...` line. |
| `runtime.julia_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | Julia thread policy for runner setup. Requires process startup (`--threads`) or script re-exec. |
| `runtime.blas_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | BLAS thread policy for runner setup and can be applied at runtime. |
| `runtime.parallel.enabled` | Bool | `true` | `true`, `false` | Master switch for in-process parallel execution of independent work items (island solves, short-circuit sweeps, contingency batches). `false` forces every parallel site onto the serial path (the same functions, not copies). |
| `runtime.parallel.max_tasks` | String | `auto` | `auto` or positive integer string | Task cap for parallel sites. `auto` resolves to `Threads.nthreads()`; the cap is applied via chunking, so it also bounds `@threads` sites. |
| `runtime.parallel.min_work_items` | Int | `4` | integer >= 1 | Work lists shorter than this run serially (avoids task overhead on tiny cases). |

In a single-threaded Julia process the parallel sites run serially even
with `runtime.parallel.enabled: true`; the startup summary prints a hint
to start with `julia --threads=auto`.

Per-phase timings in the performance profile (and `performance.log`) are
CPU-time sums across all islands/workers; the elapsed real time of a
fan-out is accounted separately under `parallel_wall_time`, and the ratio
is the achieved speedup.

Julia thread priority for `examples/powerflow/matpower_import.jl`:

1. CLI override: `--julia-threads=<N|auto|keep>`
2. Environment override: `SPARLECTRA_JULIA_THREADS`
3. YAML `runtime.julia_threads`
4. keep current process setting

```bash
julia --threads=8 --project=. examples/powerflow/matpower_import.jl
julia --project=. examples/powerflow/matpower_import.jl --julia-threads=8
```

```powershell
$env:JULIA_NUM_THREADS = "8"
julia --project=. examples/powerflow/matpower_import.jl
```

## APSLF against Newton: solve times

`examples/others/apslf_vs_nr_timing.jl` solves the shipped `sp_` cases and
synthetic tiled grids (500 to 5000 buses, `SPARLECTRA_TIMING_SIZES`) with
the rectangular Newton solver and with APSLF at orders 24, 40 and 60
(higher orders from 300 buses on, `SPARLECTRA_TIMING_ORDERS`), median of
three warm runs, convergence-radius evaluation off. It writes
`apslf_vs_nr_timing.csv` and `apslf_vs_nr_timing.svg` into
`results/apslf_vs_nr_timing/`; below: Linux, Julia 1.13.0, 16 threads.

On pure-PQ grids (the tiled grids) the series solve is faster than Newton
from 500 buses on, with cost about linear in the order; cases with PV
buses and reactive limits (the `sp_` cases) need outer passes and are
Newton's ground. Networks with controllers do not run under APSLF (use
the hybrid start, see the APSLF workshop).

![APSLF against Newton, solve time over bus count](assets/apslf_vs_nr_timing.svg)

| case | buses | PV buses | solver | order | outcome | it / passes | time |
|---|---:|---:|---|---:|---|---:|---:|
| `sp_case9` | 9 | 3 | NR | - | converged | 5 | 0.3 ms |
| `sp_case9` | 9 | 3 | APSLF | 24 | converged | 1 | 0.2 ms |
| `sp_case118` | 118 | 54 | NR | - | converged | 6 | 2.3 ms |
| `sp_case118` | 118 | 54 | APSLF | 24 | converged | 3 | 4.5 ms |
| `sp_case300` | 300 | 69 | NR | - | converged | 7 | 7.4 ms |
| `sp_case300` | 300 | 69 | APSLF | 24 | converged | 4 | 14.9 ms |
| `sp_case300` | 300 | 69 | APSLF | 40 | converged | 4 | 21.9 ms |
| `sp_case300` | 300 | 69 | APSLF | 60 | converged | 4 | 35.7 ms |
| `sp_case1354` | 1354 | 260 | NR | - | converged | 8 | 34.7 ms |
| `sp_case1354` | 1354 | 260 | APSLF | 24 | converged | 3 | 61.8 ms |
| `sp_case1354` | 1354 | 260 | APSLF | 40 | converged | 3 | 88.3 ms |
| `sp_case1354` | 1354 | 260 | APSLF | 60 | converged | 3 | 129.4 ms |
| `sp_case5` | 5 | 2 | NR | - | converged | 1 | 0.2 ms |
| `sp_case5` | 5 | 2 | APSLF | 24 | converged | 1 | 0.2 ms |
| `tiled_500` | 500 | 1 | NR | - | converged | 5 | 7.9 ms |
| `tiled_500` | 500 | 1 | APSLF | 24 | converged | 1 | 6.1 ms |
| `tiled_500` | 500 | 1 | APSLF | 40 | converged | 1 | 8.4 ms |
| `tiled_500` | 500 | 1 | APSLF | 60 | converged | 1 | 14.3 ms |
| `tiled_1000` | 1000 | 1 | NR | - | converged | 5 | 19.4 ms |
| `tiled_1000` | 1000 | 1 | APSLF | 24 | converged | 1 | 13.3 ms |
| `tiled_1000` | 1000 | 1 | APSLF | 40 | converged | 1 | 20.9 ms |
| `tiled_1000` | 1000 | 1 | APSLF | 60 | converged | 1 | 40.4 ms |
| `tiled_2000` | 2000 | 1 | NR | - | converged | 5 | 63.8 ms |
| `tiled_2000` | 2000 | 1 | APSLF | 24 | converged | 1 | 56.9 ms |
| `tiled_2000` | 2000 | 1 | APSLF | 40 | converged | 1 | 83.3 ms |
| `tiled_2000` | 2000 | 1 | APSLF | 60 | converged | 1 | 110.9 ms |
| `tiled_5000` | 5000 | 1 | NR | - | converged | 5 | 221.2 ms |
| `tiled_5000` | 5000 | 1 | APSLF | 24 | converged | 1 | 162.7 ms |
| `tiled_5000` | 5000 | 1 | APSLF | 40 | converged | 1 | 189.3 ms |
| `tiled_5000` | 5000 | 1 | APSLF | 60 | converged | 1 | 252.3 ms |

## [Output configuration](@id perf-output)

What a run writes besides its result: the console summary and its diagnostics, the log files, and the detailed CSV export with its writer settings.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `output.console_summary` | Bool | `true` | `true`, `false` | Print compact run summary to console. |
| `output.console_auto_profile` | Symbol/String | `compact` | `off`, `compact`, `full` | MATPOWER auto-profile console detail. |
| `output.console_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `summary`, `full` | Diagnostic detail on console. |
| `output.console_q_limit_events` | Symbol/String | `summary` | `off`, `summary`, `full` | Q-limit/PV→PQ event console detail. |
| `output.console_max_rows` | Int | `100` | non-negative integer | Max rows in compact console tables. |
| `output.logfile_results` | Symbol/String | `off` | `off`, `compact`, `classic`, `full` | Solved result table detail in logfile. |
| `output.detailed_result_csv_write_mode` | Symbol/String | `auto` | `auto`, `buffered`, `streaming` | Detailed CSV artifact write strategy; `auto` streams very large outputs. |
| `output.detailed_result_csv_exporter` | Symbol/String | `auto` | `auto`, `report`, `direct` | Detailed CSV row-generation path; `auto` uses the direct streaming exporter for large bus counts. |
| `output.detailed_result_csv_direct_threshold_buses` | Int | `10000` | positive integer | Bus-count threshold where `auto` switches detailed CSV export from report generation to direct streaming. |
| `output.detailed_result_csv_buffer_initial_bytes` | Int | `8388608` | non-negative integer | Initial size hint for buffered detailed CSV artifact writing. |
| `output.detailed_result_csv_buffer_max_bytes` | Int | `67108864` | positive integer | Cheap estimated-size limit above which `auto` prefers streaming. |
| `output.detailed_result_csv_streaming_threshold_rows` | Int | `100000` | positive integer | Row-count threshold above which `auto` prefers streaming. |
| `output.logfile_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `full` | Diagnostic logfile detail. |
| `output.logfile_performance` | Symbol/String | `compact` | `off`, `compact`, `full` | Performance profile logfile detail. |
| `output.logfile_warnings` | Symbol/String | `table` | `off`, `summary`, `table`, `full` | Warning representation in logfile. |

The `Jacobian cond.` line (estimate plus verdict) is always part of the
classic result output; a leftover `output.condition_number` key in a YAML
file is ignored. See [Solver](solver.md).

For API and Web UI runs, `classic` writes the result output and a compact
summary (`solver_time`, `representative_time`, iterations, final
mismatch, outcome; benchmark median and samples only in benchmark mode).
`full` adds a **Full run details** section with the effective typed
configuration, artifact options and status diagnostics.

## Diagnostics configuration

| YAML path | Type | Default | Allowed values | Meaning | Cost notes |
|---|---:|---:|---|---|---|
| `diagnostics.log_effective_config` | Bool | `false` | `true`, `false` | Log merged effective configuration. | true: low |
| `diagnostics.console_summary` | Bool | `true` | `true`, `false` | Emit compact run summary to console. | true: low |
| `diagnostics.console_auto_profile` | Symbol/String | `compact` | `off`, `compact`, `full` | Auto-profile detail on console. | `full`: medium |
| `diagnostics.console_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `summary`, `full` | Solver diagnostics detail on console. | `full`: medium/high |
| `diagnostics.console_q_limit_events` | Symbol/String | `summary` | `off`, `summary`, `full` | PV→PQ event verbosity on console. | `full`: medium |
| `diagnostics.console_max_rows` | Int | `100` | non-negative integer | Row cap for console diagnostics tables. | low |
| `diagnostics.logfile_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `full` | Diagnostic logfile detail level. | `full`: medium/high |

## Performance configuration

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `performance.enabled` | Bool | `true` | `true`, `false` | Enable performance instrumentation. |
| `performance.level` | Symbol/String | `iteration` | `off`, `summary`, `iteration`, `full` | Instrumentation detail level. |
| `performance.print_to_console` | Bool | `true` | `true`, `false` | Emit performance output to console. |
| `performance.write_to_logfile` | Bool | `true` | `true`, `false` | Emit performance output to logfile. |
| `performance.show_allocations` | Bool | `false` | `true`, `false` | Include allocation stats. |
| `performance.show_iteration_table` | Bool | `true` | `true`, `false` | Show iteration-level timing table. |
| `performance.compact_logging` | Bool | `true` | `true`, `false` | Compact performance logging format. |
| `performance.skip_reference_comparison` | Bool | `false` | `true`, `false` | Skip voltage/reference comparisons for speed. |
| `performance.skip_expensive_diagnostics` | Bool | `true` | `true`, `false` | Skip high-cost diagnostics. |
| `performance.skip_branch_neighborhood_report` | Bool | `true` | `true`, `false` | Skip branch neighborhood report. |
| `performance.max_diagnostic_rows` | Int | `25` | non-negative integer | Row cap for diagnostics tables. |

## [Power mode for repeated solves](@id power-mode)

Power mode is for the cases where one network is solved many times: a
benchmark on a persistent model, the scenario engine, N-1 contingency
sweeps, Monte Carlo, outer control loops. It is off by default; a single
run gains nothing from it.

Switching it on:

- Configuration file: `power_flow.power_mode: true` (every run of the
  service, the Web UI and `run_sparlectra` reads it).
- Library call: the keyword `power_mode = true` on `runpf!` or
  `runpf_rectangular!`, for example
  `runpf!(net, 30, 1e-8, 0; power_mode = true)`; the cache lives on the
  network, so every further solve of that `net` reuses it.
- Web UI: Settings page, open **Advanced options**, fieldset **Solver
  backend**, tick **Power mode**; the fieldset is shown for the Newton
  solver only (`power_flow.solver: rectangular`). Save the settings
  (configuration file or this case) and the N-1 and scenario runs of
  that case use it.
- KLU: load the extension with `using KLU` next to `using Sparlectra`
  (the application package loads it, so the service and the Web UI have
  it). Without it power mode factorises with UMFPACK and gains less.

Between solves the network keeps its Ybus
(reused while the fingerprint of the branch admittances, terminals and
shunts is unchanged), the symbolic analysis of the sparse LU with its
factorization object and Jacobian buffers (the context re-analyses by
itself when the Jacobian pattern changes, for example after a bus-type or
a topology change), and the Newton work arrays; the ranked mismatch
diagnostics of the final status are skipped (the maxima and the worst
row stay). Nothing else changes: same Jacobian, same polar update, same
Q-limit handling, same tolerance, the same solution to 1e-12 pu.

The sparse LU of power mode is KLU when the KLU package extension is
loaded (`using KLU` next to `using Sparlectra`; the application package
loads it, so the service and the Web UI have it), UMFPACK with the kept
analysis otherwise. KLU's numeric refactorization is 6 to 20 times faster
than UMFPACK's on power-flow Jacobians of 5k to 60k unknowns; on very
large Jacobians with heavy fill-in (an 82000-bus synthetic case) UMFPACK
is faster, which is one reason power mode is a switch and not the default.

Measured warm solve (second and later solves on the same imported
network, median of 20, single thread, flat start, Q limits off, the
grid-bench adapter's keywords), without and with the switch, KLU loaded:

| case | without | with power mode | factor | allocations per warm solve |
|---|---|---|---|---|
| case2869pegase | 47.3 ms | 7.8 ms | 6.0 | 35.7 MB to 5.7 MB |
| case9241pegase | 267.9 ms | 55.3 ms | 4.9 | 141 MB to 21 MB |
| mvlv29840 | 477.4 ms | 44.7 ms | 10.7 | 407 MB to 60 MB |

Without the KLU extension the same three items (persistent analysis,
Ybus, diagnostics) give a factor of 1.2 on case2869pegase. Through the
service API every call imports the case anew, so only the diagnostics
and the faster LU apply: the solver phase gains a factor of 1.1 to 1.3
there and the whole call 3 to 7 percent.

The rule: a single run leaves it off (the first solve pays the analysis
either way, and the ranked diagnostics are what a single run wants to
read); a loop over one network switches it on. `reset_power_mode!(net)`
drops the kept state when the memory is wanted back.

### [Dishonest Newton](@id dishonest_newton)

`power_flow.jacobian_reuse: true` (keyword `jacobian_reuse = true`, Web
UI: the checkbox "Dishonest Newton" under Solver backend, directly below
power mode) keeps the factorization of the Jacobian for the next Newton
step while the mismatch falls by at least
`jacobian_reuse_min_reduction` (default 10) per step, at most
`jacobian_reuse_max_steps` (default 3) steps in a row; a reused step that
misses the factor is replaced by the honest step. The rules are on the
[Solver Guide](@ref dishonest_newton_solver). Off by default; it changes
the iteration count, not the solution (same tolerance, final state equal
to honest Newton within it), and works with and without power mode.

Does it pay off? In short:

| use | effect of dishonest Newton |
|---|---|
| repeated solves of one large transmission case, power mode with KLU | about 7 percent faster |
| the same with UMFPACK (no KLU extension) | 12 to 15 percent faster |
| radial distribution grids | none or slightly slower |
| N-1 and scenario runs (warm starts, few steps) | slower, up to twice the time |
| a single run | no point: the factorization is not the bottleneck |

A reused step converges linearly instead of quadratically: on the
measured cases three reused steps replace two honest ones, so a run takes
two more steps and one factorization less.
That pays off only where a factorization costs much more than a step's
other work: large meshed cases without KLU. With KLU the refactorization
is cheap, and most of the gain is already in power mode.

Measured warm solve as above (median of 20, single thread, one session,
the grid-bench adapter's keywords), without and with the switch:

| case | power mode, KLU | with dishonest Newton | power mode, UMFPACK | with dishonest Newton | Newton steps (refactorisations) |
|---|---|---|---|---|---|
| case2869pegase | 7.5 ms | 7.0 ms | 30.1 ms | 25.4 ms | 5 (5) to 7 (4) |
| case9241pegase | 54.4 ms | 50.7 ms | 131.4 ms | 115.5 ms | 6 (6) to 8 (5) |
| mvlv29840 | 44.1 ms | 45.4 ms | 203.2 ms | 173.5 ms | 5 (5) to 7 (4) |

With KLU a refactorization costs little, and the gain (about 7 percent
on the two transmission cases) is close to the cost of the extra
steps; on the radial mvlv29840 the extra steps cost more than the saved
factorization. With UMFPACK the factorization is most of a step and the
switch saves 12 to 15 percent. The rule: a loop over one network switches
power mode on, and on large transmission cases dishonest Newton on top.

Not for N-1 and scenario runs. A warm-started outage solve converges in a
few fast steps, the switch adds about two steps per outage and saves a
factorization that KLU makes cheap. Full branch N-1 of case_ACTIVSg2000
(3206 outages, 16 threads):

| variant | wall time |
|---|---|
| serial | 504 s |
| parallel | 128 s |
| parallel, power mode | 69 s |
| parallel, power mode, dishonest Newton | 139 s |

Power mode is the lever for N-1 here; dishonest Newton doubles the time
(serially it costs 7 percent, 8.0 to 9.8 iterations per outage).

## [Benchmark configuration](@id perf-benchmark)

The Web UI's `performance_timing=off|compact|full` option writes
`performance.log` with the phases of one service/API request (request
parsing, case resolution, configuration, case loading/network
construction/solve, postprocessing, artifact writing, total time; `full`
adds the internal profile entries). `benchmark.enabled` instead makes
`run_matpower_case` perform repeated solves and report representative and
median timing. The repeated timing comes from the package BenchmarkTools,
which Sparlectra does not install: load it next to Sparlectra
(`using BenchmarkTools`, after `Pkg.add("BenchmarkTools")` where it is
missing). Without it a run with `benchmark.enabled: true` stops before it
reads the case and names the package. The benchmark mode is a scripting
feature: a power-flow run through the service or the Web UI solves once,
and the Web UI has no benchmark option.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `benchmark.enabled` | Bool | `false` | `true`, `false` | Enable benchmark mode of `run_matpower_case` (needs `using BenchmarkTools`). |
| `benchmark.methods` | Vector{Symbol/String} | `[rectangular]` | `rectangular` (current PF core) | Methods benchmarked. |
| `benchmark.seconds` | Float64 | `2.0` | positive real | Benchmark max. time budget. This is not a minimum runtime, solver timeout, or iteration limit; a running sample is not interrupted. |
| `benchmark.samples` | Int | `50` | positive integer | Max benchmark samples per method. The benchmark may finish earlier when this count is reached before the time budget, or collect fewer samples when the time budget is reached first. |
| `benchmark.show_once` | Bool | `false` | `true`, `false` | Run one full visible solve before timing loop. |
| `benchmark.show_once_output` | Symbol/String | `classic` | `classic`, `dataframe`, `compact` | Output format for `show_once`. |
| `benchmark.show_once_max_nodes` | Int | `0` | non-negative integer | Row cap for one-shot output. |

## Solver workspace and warmup keys

Performance-relevant keys outside the `performance.*` section:

| Key | Default | Meaning |
|---|---|---|
| `power_flow.rectangular_workspace_reuse` | `true` | Reuse the rectangular solver's workspace between solves of one session instead of reallocating. |
| `power_flow.rectangular_preallocate_workspace` | `auto` | Preallocate the workspace up front (`off`, `on`, `auto`; `auto` decides by case size). |
| `power_flow.rectangular_workspace_min_buses` | `1000` | Case size from which `auto` preallocates. |
| `performance.representative_warmup_runs` | `0` | Untimed warmup solves before a timed representative run. |
| `performance.compare_cold_warm` | `false` | Report the cold (first) and warm (subsequent) timings side by side. |
