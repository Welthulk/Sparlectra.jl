# Local PowerFlow Service

The local PowerFlow service is a synchronous call around
[`run_sparlectra_api`](@ref), the service layer over `runpf!`: one request
dictionary in, one run directory with a fixed set of artifacts and a
serialized result out. It starts no HTTP server; the [Web UI](webui.md) is
a layer above it. The call surface and result contract of the underlying
API are on [Programmatic API](programmatic_api.md);
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

**Use**

| Call | Does |
|---|---|
| [`start_powerflow_run`](@ref) | validates the request, chooses run id and directory, runs the API, registers the result and updates the index. A trusted caller such as the Web UI passes the server-owned case directory through the function keyword; bare `.m` names are resolved with `ensure_casefile` (`to_jl=false`) and stay the executed source, an explicit `.jl` request resolves to its `.m` source and is rejected when that is missing. Browser form values never control the directory, missing path-like inputs are not downloaded, URLs are rejected |
| [`load_powerflow_run_index`](@ref) | reads the transport-safe `powerflow_runs_index.json`; a missing index is an empty index |
| [`list_powerflow_runs`](@ref) | the indexed run summaries, with `available` and a structured `reason` when a run directory or `result.json` is unavailable or unsafe |
| [`refresh_powerflow_run_registry!`](@ref) | rebuilds the in-process registry from the valid `result.json` files after a restart; one missing or corrupt run does not block the others |
| [`get_powerflow_result`](@ref) | the serialized run metadata by id; `casefile` is the effective local `.m` or `.jl` path passed to the API |
| [`list_powerflow_artifacts`](@ref) | the artifact metadata discovered in the run directory |
| [`resolve_powerflow_artifact`](@ref) | exactly one named artifact of the selected run |

**Artifacts** of a run under `output_root/<run_id>/` (the index
`powerflow_runs_index.json` sits in `output_root` itself and holds run
metadata only, failed runs included once they produced `result.json`):

| File | Content |
|---|---|
| `result.json` | the machine-readable result with status, metadata and the raw phase sequence |
| `run.log` | the run's narrative: solver time, iterations, final mismatch, outcome, case file extension and size, phase timings per step, the large-case timing summary, benchmark median and sample count when benchmarking is on. `output.logfile_results=full` adds run parameters, artifact choices and status diagnostics beyond the `classic` report. Console output is captured here; `output.console_live: true` mirrors it live, the file is identical |
| `effective_config.yaml` | the resolved configuration |
| `run_metadata.yaml` | request and lifecycle metadata |
| `performance.log` | with `performance_timing` (`off`, `compact`, `full`): Y-bus, Newton-iteration, Q-limit and linear-solve aggregates; `solver_elapsed_s` is the pure solver time |
| `diagnose.log` | with `run_diagnostics`: the PowerFlow and Q-limit printers. A diagnostic failure never replaces the primary result; run directories with the legacy `diagnose.txt` stay discoverable |
| `cgmes.log` | the full CGMES import report; only its `warning:` lines are mirrored into `run.log` |
| `bus_voltages_complex.csv`, `branch_flows.csv`, `bus_powers.csv` | with `detailed_result_csv` (default `false`), after a successful solve, from `buildACPFlowReport(raw_result.net)`: polar and rectangular voltages per bus (for Excel), flows at both ends and losses per branch, and one row per bus with solved generation, load and shunt power, `bus_type_start` and `bus_type_end`, the binding Q-limit side and band, the `control`/`control_status` summary (Q(U), P(U), RVC, STATCOM, SVC, MSC, OLTC_target) and the `non_physical` Q-V flag, the data of the console result table and `printQVCharacteristicCheck` |
| `q_limit_events.csv`, `q_limit_initial_limits.csv` | Q-limit switching events and the initial limits |

**CSV request fields**

| Field | Values | Effect |
|---|---|---|
| `config_overrides["output.csv_format"]` | `technical` (comma, decimal point, no grouping), `excel_de` (semicolon, decimal comma, thousands dot), `excel_us` (comma, decimal point, thousands comma, grouped numbers quoted) | every CSV a run writes goes through the one writer `write_result_csv`: the files above, the short-circuit, contingency and scenario tables, the state-estimation state and diagnostic CSVs, the AC island report, the SV comparison and the DTF outage metrics. Key reference: [Configuration](configuration.md) |
| `detailed_result_csv_format`, `detailed_result_csv_semicolon=true` (meaning `excel_de`) | deprecated, still accepted | become the `output.csv_format` override; an explicit override wins |
| `detailed_result_csv_write_mode` | `buffered` (one write from an IOBuffer), `streaming` (rows straight to the file), `auto` (default: streams above the configured row or estimated-byte thresholds) | how the detailed CSVs are written |
| `detailed_result_csv_exporter` | `auto` | keeps the report-based path for small cases and switches to direct streaming when the bus count reaches `detailed_result_csv_direct_threshold_buses` |

**Notes**

- Public failures are dictionaries with `status`, `success`, `reason` and
  `message`; recovery lists invalid entries in `unavailable_runs` and
  continues.
- Indexed paths are normalized and constrained to `output_root/<run_id>`,
  the result file must be that directory's `result.json`, and paths are
  resolved before recovery to reject symlink or absolute-path escapes.
  Artifact resolution rejects traversal components, Windows-style arbitrary
  paths, missing artifacts and paths outside the run directory.
- The measurement CSV (`# sparlectra-measurements v1`) is a data format
  Sparlectra reads back and keeps its fixed layout; `readSEStateCSV!`
  accepts the SE state CSV in any of the three formats. Excel formats avoid
  exponent notation where practical; identifiers that look like scientific
  notation can still trigger Excel's auto-conversion warning.
- The Web UI wraps the synchronous call in `start_webui_powerflow_run`, one
  active job at a time with cooperative abort: an aborted run gets a normal
  run directory, `result.json`, index entry and `run.log` status marker,
  and partial artifacts are never taken for a solve ([Web UI](webui.md)).
  The job states, phase names and abort checkpoints are developer
  material in `DEVELOPER.md` (Web UI job lifecycle).
- Phases are reported per step in `run.log` (repeated per-iteration phases
  aggregated into one line with a count) and in `performance.log`; the raw
  sequence stays in `result.json`. `loading_julia_case`, the evaluation of
  a generated `.jl` case, can take minutes on a cold run of a very large
  literal case.

!!! details "Why it is built this way"
    Each run receives its id before execution, so its directory is known
    before anything is written and an abort or crash leaves a run directory
    that history recovery can still classify. The in-process registry gives
    fast lookup while Julia runs; the persistent index makes completed runs
    discoverable after a restart and stays small because it carries
    metadata only, the details come from each run's `result.json`. Failed
    runs are indexed for the same reason: their status, log, effective
    configuration and runtime metadata stay available for diagnosis.
    `run.log` references content that has its own artifact instead of
    repeating it, so the narrative stays readable on large cases.
