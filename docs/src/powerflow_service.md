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
| `cgmes_import_report.md` | the full CGMES import report as Markdown (files, classes, network built, DY profile when the delivery carries one, short-circuit data, importer messages); only its `warning:` lines are mirrored into `run.log` |
| `bus_voltages_complex.csv`, `branch_flows.csv`, `bus_powers.csv` | with `detailed_result_csv` (default `false`), after a successful solve, from `buildACPFlowReport(raw_result.net)`: voltages per bus (polar and rectangular); flows and losses per branch; per bus the solved generation, load and shunt power, `bus_type_start`/`bus_type_end`, the binding Q-limit, the `control`/`control_status` summary and the `non_physical` flag |
| `q_limit_events.csv`, `q_limit_initial_limits.csv` | Q-limit switching events and the initial limits, with `run_diagnostics` or `detailed_result_csv`; a run with more than five switching events writes `q_limit_events.csv` on its own and `q_limit.log` keeps a five-row preview that names the file |

**CSV request fields**

| Field | Values | Effect |
|---|---|---|
| `config_overrides["output.csv_format"]` | `auto` (default), `technical`, `excel_de`, `excel_us`; delimiters and the `auto` rule in [Configuration](configuration.md) | every CSV a run writes goes through the one writer `write_result_csv`: the files above, the short-circuit, contingency and scenario tables, the state-estimation CSVs, the AC island report, the SV comparison and the DTF outage metrics |
| `detailed_result_csv_format`, `detailed_result_csv_semicolon=true` (meaning `excel_de`) | deprecated, still accepted | become the `output.csv_format` override; an explicit override wins |
| `detailed_result_csv_write_mode` | `buffered`, `streaming`, `auto` (default: streams above the configured thresholds) | how the detailed CSVs are written |
| `detailed_result_csv_exporter` | `auto` | report-based path for small cases, direct streaming from `detailed_result_csv_direct_threshold_buses` buses; thresholds in [Output configuration](performance_profiling.md#perf-output) |

**Q-limit stability check.** After every AC run with a converged numerical
solution `run_sparlectra` checks each voltage-regulated generator bus (PV at
the start, switched or locked; reference buses, buses joined to one by a
closed bus link and the reference of an island without a slack excluded)
against the
stability conditions of Zhao, Chiang, Li, Ju (PSCC 2008): a machine at
$Q_{max}$ must not sit above its set point (condition (7)), one at
$Q_{min}$ not below it (condition (6)), a PQ machine inside its band must
hold the set point, and a PV machine is ok inside its band. Tolerances are
fixed at 1e-3 MVAr and 1e-6 pu. The check never changes the run status and
has no switch; a run without Q-limit handling is not checked. A run stopped
at `power_flow.max_iter` (one of the three reasons below) gets the rows of
its last iterate with the status `last_iterate` and no warning; that state is
no solution, the rows show where the bus types stood. After
`power_flow.dc.fallback` the net holds the DC state instead of that iterate,
and the run is `not_checked`.

| Key | Content |
|---|---|
| `diagnostics.qlimit_stability_violations` | number of buses that fail (`Int`) |
| `diagnostics.qlimit_stability_report` | one row per regulated bus from `qlimitStabilityReport`: `bus`, `bus_name`, `end_type`, `q_MVar`, `qmin_MVar`, `qmax_MVar`, `u_pu`, `uset_pu`, `reached_limit` (`max`, `min`, `none`), `check_6`, `check_7`, `in_band`, `ok`, `locked_by`, `q_source`, `q_report_MVar` |
| `diagnostics.qlimit_stability_status` | `ok`, `violations`, `last_iterate` (run stopped at `max_iter`), `not_checked` (no converged solution otherwise, or DC fallback), `qlimits_disabled` |
| `locked_by` | why the bus could not switch freely: `none`, `guard_lock_pq` (narrow or zero band under the guard), `max_switches` (kept PV by `guard.max_switches`), `lock_pv_to_pq_buses`, `oscillation` (the active-set lock after the second hit) |
| metadata `qlimit_stability_status`, `qlimit_stability_violations`, `qlimit_stability_buses`, `qlimit_stability_line` | the same verdict in the result metadata |

Q is the bus generator Q of the result report; on a network with a closed
bus link and at a bus that also carries a Q(U) generator it is computed from
the solved voltages (`q_source = voltages`).
The line `Q-limit stability violations: n`, with the failing buses when
`n > 0`, is printed with the results, written to `run.log` and shown on the
Web UI result page; a run with a violation also logs one warning that names
each failing bus. The Q-V characteristic line (`qv_characteristic_line`) is
the part of this check for the buses with a Q-limit event: on their limit
on the wrong side of the set point. In auto mode the check uses the Q-limit
settings of the attempt that produced the result.

**Reason of a run at `max_iter`.** A rectangular Newton run that uses up
`power_flow.max_iter` keeps the status `not_converged`; its `reason` says
which of three cases the last three iterations show (Zhao, Chiang, Li, Ju,
PSCC 2008 tell the numerical divergence from the type identification
divergence, see [Q-limit switching](q_limit_switching_strategy.md)):

| `reason` | Last three iterations |
|---|---|
| `type_identification_divergence` | a bus type changed between PV and PQ, and the smallest mismatch measured before the switch is at most 10 times `power_flow.tol` |
| `iteration_limit` | no type change, the mismatch fell from each iteration to the next, and the last mismatch is below the first mismatch of the run |
| `numerical_divergence` | everything else |

The status keeps `mismatch_history` (the mismatch of the equations in force
after the switching of each iteration), `mismatch_history_pre_switch` (the
same before the switch) and `type_switch_iterations`. A run with outer
controllers reports outcome `pf_failed` and keeps the solver's reason;
before 0.60 the reason was `pf_failed` as well, and the former reason
`nr_mismatch_not_converged_active_set_unstable` (three or more switching
iterations over the whole run) is gone. The Q-limit log (`net.qLimitLog`,
`q_limit_events.csv` with the column `kind`) holds the PQ->PV releases
(`release`) beside the PV->PQ switches (`hit`, `group_revert`).

**Notes**

- Public failures are dictionaries with `status`, `success`, `reason` and
  `message`; recovery lists invalid entries in `unavailable_runs` and
  continues. Indexed paths are constrained to `output_root/<run_id>`;
  artifact resolution rejects traversal, missing artifacts and paths
  outside the run directory.
- The measurement table (`# sparlectra-measurement-table v1`) keeps its
  fixed layout (comma, decimal point), also as the `measurements.csv`
  artifact of an SE run; its verdict file `measurements_verdict.csv` and
  the SE state CSV follow the format, and `readSEStateCSV!` accepts the SE
  state CSV in any of the three formats.
- The Web UI wraps the call in `start_webui_powerflow_run`, one active job
  at a time with cooperative abort; an aborted run gets a normal run
  directory, `result.json`, index entry and `run.log` status marker.
- `loading_julia_case`, the evaluation of a generated `.jl` case, can take
  minutes on a cold run of a very large literal case.
