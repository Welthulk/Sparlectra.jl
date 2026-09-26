# Developer notes

Internal architecture notes for Sparlectra.jl: the names, states and
contracts a change must keep consistent. The user-facing behaviour is in
the documentation (`docs/src`); the contribution process is in
[CONTRIBUTING.md](CONTRIBUTING.md).

## Documentation anchors are Web UI contracts

The Web UI help (`WEBUI_HELP_TOPICS` in `app/src/webui/docs.jl`) names a
documentation section per topic by the Documenter label of its heading:
`## [Text](@id page-slug)`, as `doc = "page/#page-slug"`. The label is the
contract. Rename the heading text freely, never change or drop the `@id`;
when a section moves to another page, keep the label and update the page in
the topic.

The Web UI reads nothing from `docs/src` at run time. The sections the help
pages (`/help/<topic>`) show are generated into
`app/src/webui/help_excerpts.jl` by `tools/generate_webui_help_excerpts.jl`:
the leading blocks of the section up to a character budget, without the
`!!! details` reasoning blocks, Documenter references turned into text and
relative page links into links to the published site. After editing a
referenced section, run the generator and commit the result:

```sh
julia --project=app tools/generate_webui_help_excerpts.jl
```

The docs gate runs `tools/check_webui_doc_links.jl` (every target resolves,
every hover hint fits its limit) and the generator's `--check` (the
committed excerpts match the sources). A renamed anchor or a section that
was edited without regenerating fails there.

## Web UI job lifecycle

Developer notes on the boundary between the synchronous PowerFlow service
(`app/src/api/powerflow_service.jl`, `run_api.jl`) and the asynchronous Web
UI job layer (`app/src/api/webui_jobs.jl`). The user-facing behaviour is in
`docs/src/powerflow_service.md` and `docs/src/webui.md`; this section holds
the names a change must keep consistent.

**Job states.** `start_webui_powerflow_run` wraps `start_powerflow_run` in
a background task. One job is active at a time: `queued` and `running`
block a new submission, `aborting` does not (an abandoned worker runs to
its next cancellation point and discards its result). Terminal states are
`success`, `not_converged`, `failed`, `aborted` and `aborted_unknown` (after
a hard reset, offered once a run has stayed in `aborting` for
`WEBUI_ABORT_HARD_RESET_AFTER_SECONDS` = 60 s). The job snapshot also
carries `solver_status`, `artifact_status` and `run_status`, which the
phase updates set (`solving_powerflow` -> solver running,
`postprocessing_result` -> solver completed, any `writing_*` phase ->
artifacts running and `run_status = finalizing`).

**Cooperative abort.** The service reads three optional request fields:
`cancellation_token` (a `Threads.Atomic{Bool}`), `phase_callback` and
`operation_callback`. `_check_powerflow_cancelled!` throws `PowerFlowAborted`
(`src/session.jl`) when the token is set; the solver loops read the same
token from task storage (`sparlectra_arm_abort_token!` /
`sparlectra_disarm_abort_token!` around the worker). Checkpoints: every
phase change, after case resolution, after the configuration build, before
diagnostics and artifact writing, before finalizing, after `performance.log`,
and every 1000 rows inside the detailed CSV exporter. No task is ever
interrupted from outside, so a non-interruptible call (a sparse
factorization, a large import) finishes before the abort is observed; the
status page names the phase active at the abort request. An aborted run is
persisted like a completed one: `result.json` with `status = aborted` and
`reason = aborted_by_user`, a `run.log` marker and an index entry, so history
recovery and artifact path safety are the same code path. Runs still active
when the Web UI stopped are marked as stale aborted results by
`refresh_powerflow_run_registry!`.

**Phase names.** The Web UI accepts the phases in `_WEBUI_OPERATION_LOG_PHASES`:
`resolving_case`, `checking_case_cache`, `preparing_case_file`,
`preparing_configuration`, `reading_matpower_case`, `loading_julia_case`
(`MatpowerIO.read_case` evaluating a generated `.jl` case),
`converting_matpower_case`, `parsing_matpower_file`,
`building_sparlectra_net` (`createNetFromMatPowerCase`),
`applying_import_options`, `preparing_start_values`, `solving_powerflow`,
the state-estimation phases `importing_case`, `topology_precheck`,
`state_estimation`, `estimation_diagnostics`, then
`postprocessing_result`, `writing_artifacts`, `finalizing_success`,
`finalizing_failed`, `finalizing_aborted`. A new phase must be added to that
set, or the status page will not show it. The service-side timing keys in
`performance.log` are separate (`phases[:request_parse]`,
`:case_resolution`, `:api_config_build`, `:case_loading_network_solver`,
`:solver`, `:postprocessing`, `:artifact_writing`, `:total`).

**Operation log.** `webui_operations.jsonl` under the Web UI output root
(`app/src/webui/operations.jl`): locked appends, rotated to `.1` above
1 MiB, pruned to the newest 1000 entries above 5000, entries older than
`webui.operation_log_retention_days` (default 10) dropped; the environment
variables `SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS`, `_MAX_BYTES`,
`_MAX_ENTRIES` and `_KEEP_ENTRIES` override these. Job events:
`powerflow_submitted`, `powerflow_started`, `powerflow_phase_started`,
`powerflow_lifecycle_status`, `powerflow_completed`, `powerflow_failed`,
`powerflow_aborted`, plus the `detailed_result_csv_export_*` events. Status
page polls (`autorefresh=1`) and the sysimage page poll are not logged.

## Outer-loop controller hook interface

A controller for `run_control!` subtypes `AbstractOuterController`, with a
state subtype of `AbstractControlState` and an update subtype of
`AbstractControlUpdate` (the three abstract types are exported). The
outer loop calls these methods, none of which is exported:

- `control_initialize!`
- `control_evaluate!`
- `control_propose_update!`
- `control_apply_update!`
- `control_is_converged`
- `control_is_blocked`
- `control_status`
- `control_report_rows`
- `control_trace_rows`
- `control_max_outer_iterations`

`control_max_outer_iterations` provides a controller-specific outer-loop
limit, combined with the global budget `control.max_outer_iterations`;
extend it as `Sparlectra.control_max_outer_iterations(::MyController)`.
A user-defined controller runs via `run_control!(net; controllers = [...])`.

FACTS regression coverage:

- `test/test_tap_controller.jl`: STATCOM registration, limit tracking
  (Q = V * S_max), in-range equivalence with the constant-Q mode, the
  SVC-versus-STATCOM contrast (V^2 versus V collapse).
- `test/test_series_reactance_control.jl`: SSSC registration, converged
  operation inside the live window, pinned operation at V_inj,max,
  TCSC-mode regression.
- `test/test_upfc_control.jl`: the quadrature composite equals the manual
  SSSC+STATCOM pair to machine precision, both limits at their clamps,
  all-or-nothing registration, YAML type `upfc`; the full model reaches
  independent P and Q with the DC-link balance closed, reduces to the SSSC
  in forced quadrature, round-trips through `model: full`, and short
  circuit and the exports read the base impedance afterwards.
