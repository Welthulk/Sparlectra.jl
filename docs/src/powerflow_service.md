# Local PowerFlow Service

The local PowerFlow service is a synchronous call around
[`run_sparlectra_api`](@ref): one request dictionary in, one run directory
with a fixed set of artifacts and a serialized result out. It starts no
HTTP server; the [Web UI](webui.md) is a layer above it. The result
contract and the override rules are on
[Programmatic API](programmatic_api.md);
`examples/powerflow/exp_powerflow_service.jl` is the runnable example.

```julia
using Sparlectra
using SparlectraApp   # service API: environment app/, or app/ on the load path

request = Dict(
    "casefile" => "data/mpower/case5.m",
    "config_file" => "examples/configuration.yaml",
    "output_root" => "results/powerflow_service",
    "config_overrides" => Dict("power_flow.tol" => 1e-8),
)
result = start_powerflow_run(request)
run_id = result["run_id"]

refresh_powerflow_run_registry!("results/powerflow_service")   # after a Julia restart
runs = list_powerflow_runs("results/powerflow_service")
stored_result = get_powerflow_result(run_id)
artifacts = list_powerflow_artifacts(run_id)
result_json = resolve_powerflow_artifact(run_id, "result.json")
```

| Call | Does |
|---|---|
| [`start_powerflow_run`](@ref) | validates the request, chooses run id and directory before anything is written, runs the API, registers the result and updates the index. The case directory comes from the function keyword, never from the request; bare `.m` names are resolved with `ensure_casefile` (`to_jl = false`), a `.jl` request resolves to its `.m` source; missing path-like inputs are not downloaded, URLs are rejected |
| [`load_powerflow_run_index`](@ref) | reads `powerflow_runs_index.json`; a missing index is an empty index |
| [`list_powerflow_runs`](@ref) | the indexed run summaries, with `available` and a structured `reason` when a run directory or `result.json` is unavailable or unsafe |
| [`refresh_powerflow_run_registry!`](@ref) | rebuilds the in-process registry from the valid `result.json` files after a restart; one corrupt run does not block the others |
| [`get_powerflow_result`](@ref) | the serialized run metadata by id; `casefile` is the effective local `.m` or `.jl` path passed to the API |
| [`list_powerflow_artifacts`](@ref) | the artifact metadata discovered in the run directory |
| [`resolve_powerflow_artifact`](@ref) | exactly one named artifact of the selected run |

**Artifacts** of a run under `output_root/<run_id>/`. The index
`powerflow_runs_index.json` sits in `output_root` itself and holds run
metadata only; failed runs are indexed once they produced `result.json`.

| File | Content |
|---|---|
| `result.json` | the machine-readable result: status, metadata, the raw phase sequence |
| `run.log` | the run narrative: solver time, iterations, final mismatch, outcome, phase timings per step, benchmark median and sample count when benchmarking is on; `output.logfile_results=full` adds run parameters, artifact choices and status diagnostics. Console output is captured here; `output.console_live: true` mirrors it live |
| `effective_config.yaml` | the resolved configuration |
| `run_metadata.yaml` | request and lifecycle metadata |
| `performance.log` | with `performance_timing` (`off`, `compact`, `full`): Y-bus, Newton-iteration, Q-limit and linear-solve aggregates; `solver_elapsed_s` is the pure solver time |
| `diagnose.log` | with `run_diagnostics`: the PowerFlow and Q-limit printers; a diagnostic failure never replaces the primary result |
| `cgmes.log` | the full CGMES import report; only its `warning:` lines are mirrored into `run.log` |
| `bus_voltages_complex.csv`, `branch_flows.csv`, `bus_powers.csv` | with `detailed_result_csv` (default `false`), after a successful solve, from `buildACPFlowReport(raw_result.net)`: voltages per bus (polar and rectangular); flows and losses per branch; per bus the solved generation, load and shunt power, `bus_type_start`/`bus_type_end`, the binding Q-limit, the `control`/`control_status` summary and the `non_physical` flag |
| `q_limit_events.csv`, `q_limit_initial_limits.csv` | Q-limit switching events and the initial limits, with `run_diagnostics` or `detailed_result_csv`; a run with more than five switching events writes `q_limit_events.csv` on its own and `q_limit.log` keeps a five-row preview that names the file |

**CSV request fields**

| Field | Values | Effect |
|---|---|---|
| `config_overrides["output.csv_format"]` | `auto` (default), `technical`, `excel_de`, `excel_us`; delimiters and the `auto` rule in [Configuration](configuration.md) | every CSV a run writes goes through the one writer `write_result_csv`: the files above, the short-circuit, contingency and scenario tables, the state-estimation CSVs, the AC island report, the SV comparison and the DTF outage metrics |
| `detailed_result_csv_format`, `detailed_result_csv_semicolon=true` (meaning `excel_de`) | deprecated, still accepted | become the `output.csv_format` override; an explicit override wins |
| `detailed_result_csv_write_mode` | `buffered`, `streaming`, `auto` (default: streams above the configured thresholds) | how the detailed CSVs are written |
| `detailed_result_csv_exporter` | `auto` | report-based path for small cases, direct streaming from `detailed_result_csv_direct_threshold_buses` buses; thresholds in [Output configuration](performance_profiling.md#perf-output) |

**Notes**

- Public failures are dictionaries with `status`, `success`, `reason` and
  `message`; recovery lists invalid entries in `unavailable_runs` and
  continues. Indexed paths are constrained to `output_root/<run_id>`;
  artifact resolution rejects traversal, missing artifacts and paths
  outside the run directory.
- The measurement CSV (`# sparlectra-measurements v1`) keeps its fixed
  layout; `readSEStateCSV!` accepts the SE state CSV in any of the three
  formats.
- The Web UI wraps the call in `start_webui_powerflow_run`, one active job
  at a time with cooperative abort; an aborted run gets a normal run
  directory, `result.json`, index entry and `run.log` status marker.
- `loading_julia_case`, the evaluation of a generated `.jl` case, can take
  minutes on a cold run of a very large literal case.
