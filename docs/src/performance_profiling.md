# Performance and Profiling Configuration

## [Runtime](@id perf-runtime)

Process-level settings read once at startup: the thread report, the Julia and BLAS thread counts, and the parallel work split of the sweeps.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `runtime.print_thread_config` | Bool | `true` | `true`, `false` | Print the Julia/BLAS thread summary at startup, including the `parallel: enabled=... max_tasks=...` line. |
| `runtime.julia_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | Julia thread policy for runner setup. Requires process startup (`--threads`) or script re-exec. |
| `runtime.blas_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | BLAS thread policy for runner setup; can be applied at runtime. |
| `runtime.parallel.enabled` | Bool | `true` | `true`, `false` | Master switch for in-process parallel execution of independent work items (island solves, short-circuit sweeps, contingency batches). `false` forces every parallel site onto its serial path (the same functions, not copies). |
| `runtime.parallel.max_tasks` | String | `auto` | `auto` or positive integer string | Task cap for parallel sites. `auto` resolves to `Threads.nthreads()`; the cap is applied via chunking, so it also bounds `@threads` sites. |
| `runtime.parallel.min_work_items` | Int | `4` | integer >= 1 | Work lists shorter than this run serially. |

In a single-threaded Julia process the parallel sites run serially even
with `runtime.parallel.enabled: true`; the startup summary then hints at
`julia --threads=auto`. Per-phase timings in the performance profile and
`performance.log` are CPU-time sums over all islands and workers; the
elapsed time of a fan-out is accounted separately as `parallel_wall_time`,
and the ratio is the achieved speedup. What runs in parallel:
[Parallel Execution](parallel_execution.md).

Julia thread priority for `examples/powerflow/matpower_import.jl`: the CLI
override `--julia-threads=<N|auto|keep>`, then the environment variable
`SPARLECTRA_JULIA_THREADS`, then YAML `runtime.julia_threads`, else the
current process setting.

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
from 500 buses on, with cost about linear in the order; cases with PV buses
and reactive limits (the `sp_` cases) need outer passes and are Newton's
ground. A network with controllers solves under APSLF with the controllers
left out (one warning names them); for controllers use the APSLF start
values ahead of the Newton solve
([External Solver Interface](external_solvers.md)).

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
| `output.console_summary` | Bool | `true` | `true`, `false` | Print the compact run summary to the console. |
| `output.console_auto_profile` | Symbol/String | `compact` | `off`, `compact`, `full` | MATPOWER auto-profile console detail. |
| `output.console_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `summary`, `full` | Diagnostic detail on the console. |
| `output.console_q_limit_events` | Symbol/String | `summary` | `off`, `summary`, `full` | Q-limit/PV→PQ event console detail. |
| `output.console_max_rows` | Int | `100` | non-negative integer | Max rows in compact console tables. |
| `output.logfile_results` | Symbol/String | `off` (`full` in the packaged template) | `off`, `compact`, `classic`, `full` | Solved result table detail in the logfile. A library run without a configuration file stays quiet; the template, which the Web UI and the services start from, logs the full table. |
| `output.detailed_result_csv_write_mode` | Symbol/String | `auto` | `auto`, `buffered`, `streaming` | Detailed CSV artifact write strategy; `auto` streams very large outputs. |
| `output.detailed_result_csv_exporter` | Symbol/String | `auto` | `auto`, `report`, `direct` | Detailed CSV row-generation path; `auto` uses the direct streaming exporter for large bus counts. |
| `output.detailed_result_csv_direct_threshold_buses` | Int | `10000` | positive integer | Bus count from which `auto` switches the detailed CSV export from report generation to direct streaming. |
| `output.detailed_result_csv_buffer_initial_bytes` | Int | `8388608` | non-negative integer | Initial size hint for buffered detailed CSV writing. |
| `output.detailed_result_csv_buffer_max_bytes` | Int | `67108864` | positive integer | Estimated size above which `auto` prefers streaming. |
| `output.detailed_result_csv_streaming_threshold_rows` | Int | `100000` | positive integer | Row count above which `auto` prefers streaming. |
| `output.logfile_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `full` | Diagnostic logfile detail. |
| `output.logfile_performance` | Symbol/String | `compact` (`full` in the packaged template) | `off`, `compact`, `full` | Performance profile logfile detail. |
| `output.logfile_warnings` | Symbol/String | `table` (`full` in the packaged template) | `off`, `summary`, `table`, `full` | Warning representation in the logfile. |

The `Jacobian cond.` line (estimate plus verdict) is always part of the
classic result output; a leftover `output.condition_number` key is ignored
([Solver](solver.md)). For API and Web UI runs, `classic` writes the result
output and a compact summary (`solver_time`, `representative_time`,
iterations, final mismatch, outcome; benchmark median and samples only in
benchmark mode); `full` adds a **Full run details** section with the
effective typed configuration, artifact options and status diagnostics.

## Diagnostics configuration

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `diagnostics.log_effective_config` | Bool | `false` | `true`, `false` | Log the merged effective configuration. |

The former `diagnostics.console_*` and `diagnostics.logfile_diagnostics`
keys are deprecated aliases of the `output.*` keys above and are ignored
with a warning ([Configuration](configuration.md)).

## Performance configuration

The defaults are those of the packaged template
(`src/config/configuration.yaml.example`); a file that omits a key, and a
library run without a configuration file, fall back to `enabled: false`,
`level: summary` and `show_iteration_table: false`.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `performance.enabled` | Bool | `true` | `true`, `false` | Enable performance instrumentation. |
| `performance.level` | Symbol/String | `iteration` | `off`, `summary`, `iteration`, `full` | Instrumentation detail level. |
| `performance.print_to_console` | Bool | `true` | `true`, `false` | Emit performance output to the console. |
| `performance.write_to_logfile` | Bool | `true` | `true`, `false` | Emit performance output to the logfile. |
| `performance.show_allocations` | Bool | `false` | `true`, `false` | Include allocation stats. |
| `performance.show_iteration_table` | Bool | `true` | `true`, `false` | Show the iteration-level timing table. |
| `performance.compact_logging` | Bool | `true` | `true`, `false` | Compact performance logging format. |
| `performance.skip_reference_comparison` | Bool | `false` | `true`, `false` | Skip voltage/reference comparisons for speed. |
| `performance.skip_expensive_diagnostics` | Bool | `true` | `true`, `false` | Skip high-cost diagnostics. |
| `performance.skip_branch_neighborhood_report` | Bool | `true` | `true`, `false` | Skip the branch neighborhood report. |
| `performance.max_diagnostic_rows` | Int | `25` | non-negative integer | Row cap for diagnostics tables. |

## [Power mode for repeated solves](@id power-mode)

Power mode is for the cases where one network is solved many times: a
benchmark on a persistent model, the scenario engine, N-1 contingency
sweeps, Monte Carlo, outer control loops. It is off by default; a single
run gains nothing from it.

Switch it on with `power_flow.power_mode: true` in the configuration file
(read by the service, the Web UI and `run_sparlectra`), or with the keyword
`power_mode = true` on `runpf!` or `runpf_rectangular!`
(`runpf!(net, 30, 1e-8, 0; power_mode = true)`; the cache lives on the
network, so every further solve of that `net` reuses it). Web UI: Settings
page, **Advanced options**, fieldset **Solver backend** (shown for the
Newton solver only, `power_flow.solver: rectangular`), tick **Power mode**
and save; the N-1 and scenario runs of that case then use it. With the KLU
extension (`using KLU` next to `using Sparlectra`; the application package
loads it, so the service and the Web UI have it) power mode factorises
with KLU, without it with UMFPACK and gains less; `power_flow.power_mode_lu`
chooses between the two, see [below](@ref power_mode_lu).

Between solves the network keeps its Ybus (reused while the fingerprint of
the branch admittances, terminals and shunts is unchanged), the symbolic
analysis of the sparse LU with its factorization object and Jacobian buffers
(re-analysed by itself when the Jacobian pattern changes, for example after
a bus-type or a topology change) and the Newton work arrays; the ranked
mismatch diagnostics of the final status are skipped (the maxima and the
worst row stay). Same Jacobian, same polar update, same Q-limit handling,
same tolerance, the same solution to 1e-12 pu. `reset_power_mode!(net)`
drops the kept state.

Measured warm solve (second and later solves on the same imported network,
median of 20, single thread, flat start, Q limits off, the grid-bench
adapter's keywords), without and with the switch, KLU loaded:

| case | without | with power mode | factor | allocations per warm solve |
|---|---|---|---|---|
| case2869pegase | 47.3 ms | 7.8 ms | 6.0 | 35.7 MB to 5.7 MB |
| case9241pegase | 267.9 ms | 55.3 ms | 4.9 | 141 MB to 21 MB |
| mvlv29840 | 477.4 ms | 44.7 ms | 10.7 | 407 MB to 60 MB |

Without the KLU extension the same three items (persistent analysis, Ybus,
diagnostics) give a factor of 1.2 on case2869pegase. Through the service API
every call imports the case anew, so only the diagnostics and the faster LU
apply: the solver phase gains a factor of 1.1 to 1.3 and the whole call 3 to
7 percent.

### [Power-mode LU: KLU or UMFPACK](@id power_mode_lu)

KLU's numeric refactorization is 6 to 20 times faster than UMFPACK's on
power-flow Jacobians of 5k to 60k unknowns. On very large Jacobians it is
the other way round: KLU's partial pivoting takes thousands of
off-diagonal pivots there, its factor grows well beyond the estimate of
its analysis, and on the 82000-bus case_SyntheticUSA a power-mode run
with KLU took more than twice as long as one without power mode. Neither
backend is right for every network, so `auto` chooses per network:

| `power_flow.power_mode_lu` | sparse LU of power mode |
|---|---|
| `auto` (default) | KLU's symbolic analysis of the first power-mode Jacobian of a network (`klu_analyze`: block triangular form and AMD ordering, no numeric work) estimates the factorization flops; up to `Sparlectra.POWER_MODE_LU_KLU_MAX_FLOPS` (5e7) power mode factorizes with KLU, above it with UMFPACK; without the KLU extension always UMFPACK |
| `klu` | KLU (UMFPACK with a warning when the extension is not loaded) |
| `umfpack` | UMFPACK |

The threshold is fitted on measured Jacobians: KLU refactorized faster on
case2869pegase, case9241pegase, mvlv29840, case_ACTIVSg2000,
case13659pegase and the two small islands of case_SyntheticUSA (estimates
1.8e6 to 2.7e7 flops), UMFPACK on case_ACTIVSg25k and the 70000-bus island
of case_SyntheticUSA (1.3e8 and 5.6e8). The estimate cannot see KLU's
numeric pivoting; it stands in for it through the size of the problem, so a
network far from the measured ones may get the slower backend; `klu` or
`umfpack` then fix the choice.

The choice is made once per network, at the first factorization of its
first power-mode solve, and kept for every later solve of that network
(outages, other bus types, other Ybus patterns, N-1 and scenario workers,
which are copies of the solved base case). An island solved on its own
decides for itself; an island an outage splits off takes the choice of its
network. A newly imported network decides again, as does a network after
`reset_power_mode!(net)`. The decision costs one symbolic analysis, which
KLU's factorization continues from when KLU wins. The solution is the same
with every value.

The run log names the choice in the result header with the estimate, for
example

```text
Power-mode LU  : klu (KLU symbolic estimate 1.95e+06 flops, at most the threshold 5e+07; decided in this solve)
```

and the run metadata carries it as `power_mode_lu_choice`,
`power_mode_lu_source`, `power_mode_lu_klu_est_flops` and
`power_mode_lu_threshold_flops` (one row per island in
`power_mode_lu_islands`). Web UI: the select **Power-mode LU** directly
under the power-mode checkbox, greyed while power mode is off.

### [Dishonest Newton](@id dishonest_newton)

`power_flow.jacobian_reuse: true` (keyword `jacobian_reuse = true`; Web UI:
the checkbox "Dishonest Newton" under Solver backend, below power mode)
keeps the Jacobian factorization for the next Newton step while the
mismatch falls by at least `jacobian_reuse_min_reduction` (default 10) per
step, at most `jacobian_reuse_max_steps` (default 3) steps in a row. Off by
default; it changes the iteration count, not the solution, and works with
and without power mode. The rules are on the
[Solver Guide](@ref dishonest_newton_solver).

| use | effect of dishonest Newton |
|---|---|
| repeated solves of one large transmission case, power mode with KLU | about 7 percent faster |
| the same with UMFPACK (no KLU extension) | 12 to 15 percent faster |
| radial distribution grids | none or slightly slower |
| N-1 and scenario runs (warm starts, few steps) | slower (6 percent in the N-1 measurement below) |
| a single run | no point: the factorization is not the bottleneck |

A reused step converges linearly instead of quadratically: on the measured
cases three reused steps replace two honest ones, so a run takes two more
steps and one factorization less. That pays only where a factorization
costs much more than the rest of a step: large meshed cases without KLU.

Measured warm solve as above (median of 20, single thread, one session, the
grid-bench adapter's keywords), without and with the switch:

| case | power mode, KLU | with dishonest Newton | power mode, UMFPACK | with dishonest Newton | Newton steps (refactorisations) |
|---|---|---|---|---|---|
| case2869pegase | 7.5 ms | 7.0 ms | 30.1 ms | 25.4 ms | 5 (5) to 7 (4) |
| case9241pegase | 54.4 ms | 50.7 ms | 131.4 ms | 115.5 ms | 6 (6) to 8 (5) |
| mvlv29840 | 44.1 ms | 45.4 ms | 203.2 ms | 173.5 ms | 5 (5) to 7 (4) |

With KLU the gain is close to the cost of the extra steps, and on the
radial mvlv29840 the extra steps cost more than the saved factorization;
with UMFPACK the factorization is most of a step.

Not for N-1 and scenario runs: a warm-started outage solve converges in a
few fast steps, and the switch adds steps per outage to save a
factorization that KLU makes cheap. Full branch N-1 of case_ACTIVSg2000
(3206 outages):

| variant | wall time |
|---|---|
| serial | 343 s |
| parallel | 83 s |
| parallel, power mode | 70 s |
| parallel, power mode, dishonest Newton | 74 s |

Measured with Sparlectra 0.30.2: `runContingencies!` on the outages of
`generateN1Branches`, 16 threads for the parallel rows, KLU loaded, the
four variants interleaved, median of three series each, after a warm-up of
every variant (compile time not included); every variant gives the same
per-outage results as the serial one.

## [Benchmark configuration](@id perf-benchmark)

The Web UI's `performance_timing=off|compact|full` option writes
`performance.log` with the phases of one service/API request (request
parsing, case resolution, configuration, case loading, network
construction, solve, postprocessing, artifact writing, total; `full` adds
the internal profile entries). `benchmark.enabled` instead makes
`run_matpower_case` repeat the solve and report representative and median
timing; it is a scripting feature, a run through the service or the Web UI
solves once, and the Web UI has no benchmark option. The repeated timing
needs the package BenchmarkTools, which Sparlectra does not install: load
it next to Sparlectra (`using BenchmarkTools`, after
`Pkg.add("BenchmarkTools")` where it is missing); without it a run with
`benchmark.enabled: true` stops before it reads the case and names the
package.

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `benchmark.enabled` | Bool | `false` | `true`, `false` | Enable the benchmark mode of `run_matpower_case` (needs `using BenchmarkTools`). |
| `benchmark.methods` | Vector{Symbol/String} | `[rectangular]` | `rectangular` (current PF core) | Methods benchmarked. |
| `benchmark.seconds` | Float64 | `2.0` | positive real | Benchmark time budget. Not a minimum runtime, solver timeout or iteration limit; a running sample is not interrupted. |
| `benchmark.samples` | Int | `30` | positive integer | Max benchmark samples per method; the run ends at this count or at the time budget, whichever comes first. |
| `benchmark.show_once` | Bool | `false` | `true`, `false` | Run one full visible solve before the timing loop. |
| `benchmark.show_once_output` | Symbol/String | `classic` | `classic`, `dataframe`, `compact` | Output format for `show_once`. |
| `benchmark.show_once_max_nodes` | Int | `0` | non-negative integer | Row cap for the one-shot output. |

## Solver workspace and warmup keys

Performance-relevant keys outside the `performance.*` section:

| Key | Default | Meaning |
|---|---|---|
| `power_flow.rectangular_workspace_reuse` | `true` | Reuse the rectangular solver's workspace between solves of one session instead of reallocating. |
| `power_flow.rectangular_preallocate_workspace` | `auto` | Preallocate the workspace up front (`off`, `on`, `auto`; `auto` decides by case size). |
| `power_flow.rectangular_workspace_min_buses` | `1000` | Case size from which `auto` preallocates. |
| `performance.representative_warmup_runs` | `0` | Untimed warmup solves before a timed representative run. |
| `performance.compare_cold_warm` | `false` | Report the cold (first) and warm (subsequent) timings side by side. |
