# Performance and Profiling Configuration

## Runtime

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `runtime.print_thread_config` | Bool | `true` | `true`, `false` | Print Julia/BLAS thread summary at startup. |
| `runtime.julia_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | Julia thread policy for runner setup. Requires process startup (`--threads`) or script re-exec. |
| `runtime.blas_threads` | String | `keep` | `keep`, `default`, `off`, `auto`, or integer-like string | BLAS thread policy for runner setup and can be applied at runtime. |

For `examples/matpower_import.jl`, Julia thread priority is:

1. CLI override: `--julia-threads=<N|auto|keep>`
2. Environment override: `SPARLECTRA_JULIA_THREADS`
3. YAML `runtime.julia_threads`
4. keep current process setting

Example startup commands:

```bash
julia --threads=8 --project=. examples/matpower_import.jl
julia --project=. examples/matpower_import.jl --julia-threads=8
```

```powershell
$env:JULIA_NUM_THREADS = "8"
julia --project=. examples/matpower_import.jl
```

## Output configuration

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `output.console_summary` | Bool | `true` | `true`, `false` | Print compact run summary to console. |
| `output.console_auto_profile` | Symbol/String | `compact` | `off`, `compact`, `full` | MATPOWER auto-profile console detail. |
| `output.console_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `summary`, `full` | Diagnostic detail on console. |
| `output.console_q_limit_events` | Symbol/String | `summary` | `off`, `summary`, `full` | Q-limit/PV→PQ event console detail. |
| `output.console_max_rows` | Int | `5` | non-negative integer | Max rows in compact console tables. |
| `output.logfile_results` | Symbol/String | `off` | `off`, `compact`, `classic`, `full` | Solved result table detail in logfile. |
| `output.logfile_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `full` | Diagnostic logfile detail. |
| `output.logfile_performance` | Symbol/String | `compact` | `off`, `compact`, `full` | Performance profile logfile detail. |
| `output.logfile_warnings` | Symbol/String | `table` | `off`, `summary`, `table`, `full` | Warning representation in logfile. |

## Diagnostics configuration

| YAML path | Type | Default | Allowed values | Meaning | Cost notes |
|---|---:|---:|---|---|---|
| `diagnostics.log_effective_config` | Bool | `false` | `true`, `false` | Log merged effective configuration. | true: low |
| `diagnostics.console_summary` | Bool | `true` | `true`, `false` | Emit compact run summary to console. | true: low |
| `diagnostics.console_auto_profile` | Symbol/String | `compact` | `off`, `compact`, `full` | Auto-profile detail on console. | `full`: medium |
| `diagnostics.console_diagnostics` | Symbol/String | `compact` | `off`, `compact`, `summary`, `full` | Solver diagnostics detail on console. | `full`: medium/high |
| `diagnostics.console_q_limit_events` | Symbol/String | `summary` | `off`, `summary`, `full` | PV→PQ event verbosity on console. | `full`: medium |
| `diagnostics.console_max_rows` | Int | `20` | non-negative integer | Row cap for console diagnostics tables. | low |
| `diagnostics.logfile_diagnostics` | Symbol/String | `full` | `off`, `compact`, `full` | Diagnostic logfile detail level. | `full`: medium/high |

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

## Benchmark configuration

| YAML path | Type | Default | Allowed values | Meaning |
|---|---:|---:|---|---|
| `benchmark.enabled` | Bool | `true` | `true`, `false` | Enable benchmark mode. |
| `benchmark.methods` | Vector{Symbol/String} | `[rectangular]` | `rectangular` (current PF core) | Methods benchmarked. |
| `benchmark.seconds` | Float64 | `2.0` | positive real | Benchmark target duration. |
| `benchmark.samples` | Int | `50` | positive integer | Max benchmark samples. |
| `benchmark.show_once` | Bool | `false` | `true`, `false` | Run one full visible solve before timing loop. |
| `benchmark.show_once_output` | Symbol/String | `classic` | `classic`, `dataframe`, `compact` | Output format for `show_once`. |
| `benchmark.show_once_max_nodes` | Int | `0` | non-negative integer | Row cap for one-shot output. |
