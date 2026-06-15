# Local PowerFlow Service

Web UI jobs track `current_phase`, `phase_started_at`, `last_progress_at`, and
`abort_requested_at`. Cancellation boundaries cover case resolution,
configuration, case loading, Y-bus/start construction, Newton and Q-limit
iteration boundaries, sparse linear solves, diagnostics, artifact writing, and
final success persistence. A late abort is checked before success is
published. If an in-process numerical call cannot return within the Web UI
timeout, the explicit hard-reset path records an invalid `aborted_unknown`
result and shuts down the local Web UI instead of unsafely killing a Julia
task.

The local PowerFlow service is a small, testable boundary above
[`run_sparlectra_api`](@ref). It is intended for a future local Genie.jl web GUI,
but it does not start an HTTP server and has no Genie dependency.

```julia
using Sparlectra

request = Dict(
    "casefile" => "data/mpower/case5.m",
    "config_file" => "examples/configuration.yaml",
    "output_root" => "results/powerflow_service",
    "config_overrides" => Dict(
        "power_flow.tol" => 1e-8,
    ),
)

result = start_powerflow_run(request)
run_id = result["run_id"]

# Later, after restarting Julia:
refresh_powerflow_run_registry!("results/powerflow_service")

runs = list_powerflow_runs("results/powerflow_service")
stored_result = get_powerflow_result(run_id)
artifacts = list_powerflow_artifacts(run_id)
result_json = resolve_powerflow_artifact(run_id, "result.json")
```

Each run receives a unique ID before execution and writes to a GUI-friendly,
deterministic location. Lightweight metadata is also written to the persistent
index:

```text
output_root/
  powerflow_runs_index.json
  <run_id>/
    run.log
    result.json
    effective_config.yaml
    ...
```

The in-process registry provides fast lookup while Julia is running. The
`powerflow_runs_index.json` file makes completed runs discoverable after a
process restart. It stores only lightweight run metadata; full run details are
read from each run's `result.json`. Failed API runs are indexed when they
produce `result.json`, so their status, log, and effective configuration remain
available for diagnosis.

## Service boundary

The Web UI uses `start_webui_powerflow_run` as a small asynchronous lifecycle
layer around the synchronous service call. It keeps one active local Web UI job,
tracks queued/running/success/failed/aborted states, and exposes cooperative
abort requests. Cancellation is checked around case resolution, configuration
and import work, solving, diagnostics, and artifact writing. The rectangular
solver checks before and after Y-bus construction and start projection, at
every Newton iteration boundary, after Q-limit active-set work, and after each
Newton step. A currently executing sparse linear solve remains
non-interruptible, but the following check terminates the run. No unsafe task
interruption is used.

The Web UI supplies a small lifecycle callback to this asynchronous boundary so
`powerflow_started`, completion/failure, and final abort events can be appended
to the user Web UI `logs/webui_operations.jsonl` support log. This does not
change the PowerFlow result schema or per-run artifact contract. Every event
includes the running Sparlectra version and a millisecond-precision UTC
timestamp ending in `Z`. Marked `autorefresh=1` status requests are omitted from
user-action logging. Above 10,000 valid entries, the JSONL log atomically keeps
its newest 1,000 valid entries.

Aborted runs receive a normal run directory, `result.json`, an index entry, and
a `run.log` status marker. This keeps history recovery and artifact path safety
identical to completed runs while making it clear that any partial artifacts do
not represent a successful solve. The worker reaches terminal `aborted` only
after those records and the `powerflow_aborted` operation event are written;
its `finally` cleanup prevents queued, running, or aborting state from leaking
if controlled cancellation exits the worker.

- [`start_powerflow_run`](@ref) validates the dictionary-like request, chooses
  the run ID and directory, invokes the programmatic API, registers the result,
  and updates the persistent index. A trusted caller such as the Web UI can
  provide the server-owned MATPOWER case directory through the function
  keyword; bare `.m` or `.jl` names are then resolved with `ensure_casefile`
  using `to_jl=true`. Existing or generated `.jl` cases are preferred, while a
  readable `.m` case remains the fallback if conversion fails. Browser form
  values never control this directory, missing path-like inputs are not
  downloaded, and URLs are rejected.
  Optional `performance_timing` (`off`, `compact`, or `full`),
  `run_diagnostics`, `detailed_result_csv`, and
  `detailed_result_csv_format` request fields are forwarded to the API artifact
  writer. They produce `performance.log`, `diagnose.log`, and the detailed CSV
  artifacts, respectively.
  `detailed_result_csv` defaults to `false`; when enabled after a successful
  solve it writes `bus_voltages_complex.csv` and `branch_flows.csv` from
  `buildACPFlowReport(raw_result.net)`. The bus file includes polar and numeric
  rectangular voltage columns for Excel, while the branch file includes
  active/reactive flows at both ends and losses.
  `detailed_result_csv_format` accepts `technical` (default comma delimiter,
  decimal point, no grouping), `excel_de` (semicolon delimiter, decimal comma,
  thousands dot), or `excel_us` (comma delimiter, decimal point, thousands
  comma). Numeric fields containing US thousands commas are quoted. The legacy
  `detailed_result_csv_semicolon=true` request remains accepted and maps to
  `excel_de` when no explicit format is supplied. Legacy run directories
  containing `diagnose.txt` remain discoverable and downloadable.
  Diagnostic generation reuses existing PowerFlow/Q-limit printers and a
  diagnostic failure does not replace the primary PowerFlow result.
- [`load_powerflow_run_index`](@ref) reads the transport-safe index structure. A
  missing index returns an empty index.
- [`list_powerflow_runs`](@ref) returns indexed run summaries for a future run
  history table. Entries include `available` and a structured `reason` when the
  run directory or `result.json` is unavailable or unsafe.
- [`refresh_powerflow_run_registry!`](@ref) reloads valid `result.json` files
  and reconstructs the in-process registry after restart. One missing or corrupt
  run does not prevent other runs from loading.
- [`get_powerflow_result`](@ref) returns serialized run metadata by run ID.
  Its `casefile` field records the effective local `.m` or `.jl` path passed to
  the programmatic PowerFlow API.
- [`list_powerflow_artifacts`](@ref) returns metadata discovered inside the
  registered run directory.
- [`resolve_powerflow_artifact`](@ref) resolves only an exact artifact metadata
  name belonging to the selected run.

Indexed paths are normalized and constrained to `output_root/<run_id>`. The
result file must be that run directory's `result.json`, and existing paths are
resolved before recovery to reject symlink or absolute-path escapes. Artifact
resolution separately rejects traversal components, Windows-style arbitrary
paths, missing artifacts, and paths that resolve outside the run directory.

Public service failures are dictionaries containing `status`, `success`,
`reason`, and `message`. Recovery reports invalid entries in
`unavailable_runs` and continues loading valid runs.

Every completed API `run.log` includes solver time where available,
representative time, iterations, final mismatch, and final outcome. Benchmark
median and configured samples are included when benchmarking is enabled.
`output.logfile_results=full` adds effective configuration, artifact choices,
and available status diagnostics beyond the `classic` report.

This layer intentionally introduces no HTTP routes, Genie.jl server, browser
GUI, authentication, or database. See `examples/exp_powerflow_service.jl` for a
runnable local example.

## Web UI runtime paths and cancellation

The browser submit route schedules work and redirects to the run status page before case loading, solving, diagnostics, or artifact writing. Active status pages refresh every two seconds until terminal state while retaining the manual refresh link. Active status and history/form banners expose a POST-only Abort control while the job is queued or running. Cancellation is cooperative: the UI changes to `aborting` immediately, repeated requests are logged as already requested, and the worker records terminal `aborted` plus `powerflow_aborted` at the next safe phase or rectangular-iteration check. Terminal abort releases the single-active-job guard and stops automatic refresh. Delete requests for queued, running, or aborting jobs are rejected with a controlled explanation; terminal aborted jobs can be deleted normally. The Web UI provisions configuration, case cache, run output, and operation-log paths beneath the user-writable Web UI application directories; explicit startup paths and explicit local case paths remain supported.
