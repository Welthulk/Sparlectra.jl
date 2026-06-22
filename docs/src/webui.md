# Local PowerFlow Web UI

## Cooperative abort and hard reset

Active runs expose their current phase, phase start, last progress, and
abort-request time. Cancellation is cooperative. If abort is requested during
`linear_solve`, the current sparse solve must return before the cancellation
check can run.

After 60 seconds in `aborting`, the status page offers **Hard reset Web UI**.
This does not inject an exception into numerical code: it marks the run
`aborted_unknown`, records that the result is invalid, and requests a clean
server shutdown. Restart with `julia --project=. start_webui.jl`.
Cooperative cancellation logs `powerflow_aborted`; the fallback logs
`webui_hard_reset_requested` and `webui_shutdown_requested`.

Sparlectra includes a small, local-first browser interface above the existing
PowerFlow service. The Web UI contains presentation, form parsing, and route
handling only; numerical execution continues through `start_powerflow_run` and
`run_sparlectra_api`.

## Start after package installation

```julia
using Sparlectra

server = Sparlectra.start_sparlectra_webui(open_browser = true)
wait(server.task)
```

The package installation directory does not need to be known. By default,
results are written beneath `%LOCALAPPDATA%\Sparlectra\WebUI\runs` on Windows,
`$XDG_STATE_HOME/sparlectra/webui/runs` (or
`~/.local/state/sparlectra/webui/runs`) on Linux, and
`~/Library/Application Support/Sparlectra/WebUI/runs` on macOS. Directories
are created automatically. The operation log is in the sibling user Web UI `logs` directory, and downloaded/generated MATPOWER cases are
cached in the sibling user Web UI `data/mpower` directory. On first start, `warmup_case3.jl` is copied there as a small selectable demo case.

On first startup, the Web UI copies the package configuration template to its user-writable `config/configuration.yaml`. Pass `output_root="my_sparlectra_runs"` or `config_file="my_configuration.yaml"` to override these defaults; an explicit configuration file is never overwritten.
The effective configuration, output root, MATPOWER cache, and operation log are displayed by
the Web UI; the browser cannot change the output root.

### Repository developer launcher

From a repository checkout, run:

```sh
julia --project=. start_webui.jl
```

`start_webui.jl` is the single maintained developer launcher and delegates
startup and default-path behavior to `start_sparlectra_webui`.

The call returns a `SparlectraWebUIServer` handle immediately. Stop it with
`close(server)`, `Ctrl+C`, or the **Stop Web UI** button in the shared page
header. The button sends `POST /webui/shutdown`, closes the listening socket,
and allows `wait(server.task)` to return. Pass `open_browser=true` to open the
Web UI in a standalone app-style window without normal browser tabs or controls:

```julia
server = start_sparlectra_webui(open_browser = true)
```

The app-window launcher supports Microsoft Edge, Google Chrome, Chromium, and
Brave. If none is installed, Sparlectra logs the local URL instead of falling
back to a regular tab; open `http://127.0.0.1:8080/powerflow` manually if
needed.

For safety, the first prototype accepts loopback hosts only: `127.0.0.1`,
`localhost`, or `::1`. It is not intended for public or multi-user deployment.

The shared page header uses the existing Sparlectra documentation logo from
`docs/src/assets/logo.png`. The Web UI serves that single PNG through its local
asset route, and no additional branding configuration is required.
The shared header and footer display the running package version from
`Sparlectra.version()`, for example `Sparlectra.jl v0.8.4`.

### Startup warm-up

`start_sparlectra_webui` accepts `warmup`, `warmup_casefile`, and
`warmup_store_result`. The library default is `warmup=false`, while
`start_webui.jl` enables it for the repository developer workflow.
Warm-up runs asynchronously and fails softly. With
`warmup_store_result=false`, it uses a temporary directory and creates no
normal run-history entry. The bundled `data/webui/warmup_case3.jl` file is an
original Sparlectra-owned synthetic three-bus MATPOWER-compatible Julia case,
not a derivative of an external MATPOWER case.

## Starting a PowerFlow run

Web UI submissions start in a background worker task and redirect immediately to
a run-status page. Near the top, a highlighted clock card shows elapsed time in
`HH:MM:SS` format beside the current status. The technical details table retains
the raw `elapsed_seconds` value along with the requested and resolved case paths,
start time, and a manual refresh link. Queued, running, and
aborting pages also refresh every two seconds through an HTML refresh directive.
The marked `autorefresh=1` requests are not recorded as user actions. Terminal
success, failure, and abort pages stop refreshing automatically. While a job is queued or
running, it also shows an **Abort run** form that sends
`POST /powerflow/abort/<run-id>`; JavaScript is not required.

The start form and run-history page also show a prominent active-run banner
with **Open status** and POST-only **Abort** controls. This keeps abort
discoverable even when the user navigates away from the status page. The
controls disappear as soon as the run reaches a completed, failed, or aborted
state.

Abort is cooperative and never kills a Julia task unsafely. The abort request
changes the visible state to `aborting` immediately, and the Web UI PowerFlow
path checks cancellation before and after major service phases and inside each
rectangular Newton iteration. The rectangular path also checks immediately
before and after Y-bus construction and start projection, after Q-limit active
set work, and after each Newton step. A sparse factorization or other
non-interruptible operation may still finish before the next check, so the
status page shows the current phase and explains that the phase may need to
finish before cancellation is observed.
Repeated requests are idempotent. Once cancellation
is observed, the terminal state becomes `aborted`, `powerflow_aborted` is
written to the operation log, the active-run guard is released, and a new
submission is accepted. Aborted runs retain a normal run directory,
`result.json`, and `run.log` status marker, but are never reported as success.

Large MATPOWER cases can spend substantial time before Newton iterations begin,
especially while reading `.m` files, evaluating large `.jl` literal cases,
building the Sparlectra network, assembling Y-bus data, or preparing start
values. The result page, operation log, `run.log`, `result.json`, and
`performance.log` (when timing is enabled) now expose service phase timings so
users can distinguish reader, converter/cache, network-builder, solver, and
artifact costs. These diagnostics do not guarantee that very large cases finish
quickly through the local Web UI; they identify where follow-up optimization
should focus.
The operation log intentionally records only high-level phase starts; detailed
Y-bus, Newton-iteration, Q-limit, and linear-solve timings belong to each run's
`performance.log` and are summarized there instead of being repeated in the Web
UI support log.

Deletion of a queued, running, or aborting run is rejected with an explanation.
After the run reaches terminal `aborted` status, normal deletion is available.

The start page accepts:

- a typed MATPOWER case name or an existing local case path;
- a Sparlectra configuration template file;
- read-only information showing the server's configured output-root directory;
- PowerFlow tolerance and maximum iterations;
- autodamping and its minimum factor;
- Q-limit handling;
- wrong-branch detection mode;
- angle and voltage start modes;
- a visible MATPOWER import conventions section with auto-profile mode
  (`off`, `recommend`, or `apply`) plus manual transformer-ratio,
  phase-shift, bus-shunt, PV-voltage-source, and comparison-reference
  overrides;
- logfile result mode;
- single-run performance timing detail;
- optional post-run diagnostics; and
- benchmark enablement, samples, and seconds.

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are submitted. The output root is
chosen only through `start_sparlectra_webui(; output_root=...)`, is displayed as
read-only information, and cannot be overridden by a submitted browser field.
The page does not offer a generic YAML editor and never modifies the selected
template. The service creates an `effective_config.yaml` artifact for each run.
The header and `webui_operations.jsonl` include the Sparlectra version, package
path, and local Git commit when available so users can confirm which checkout
is serving the browser page.

### Case-specific settings profiles

Case-specific Web UI settings are optional. When a terminal result page shows a
successful or converged run, the compact **Case settings** section offers
**Save settings for this case**. This writes only the Web UI form options that
were used for that completed run, plus traceability metadata, into a sanitized
YAML profile below the Web UI output root. It does not save the run's full
`effective_config.yaml`, solver internals, artifact paths, or transient
convergence diagnostics.

If the run did not converge, the result page does not show the normal save
action. It instead labels the action **Save these settings anyway** and records
that the user explicitly overrode the non-successful-run warning. No profile is
saved automatically.

When the same MATPOWER case is opened again with a saved profile, the form is
prefilled with the profile values and displays a small notice. Precedence stays
conservative: built-in defaults are loaded first, global configuration remains
unchanged, the case-specific Web UI profile only prefills editable form fields,
and any manual browser edit wins for the submitted run.

The existing-case selector lists canonical MATPOWER `.m` files from the user Web UI case
cache. Generated MATPOWER `.jl` cache artifacts are internal and are not user-selectable.
First startup still provisions the small `warmup_case3.jl` demo case for warmup,
but generated MATPOWER cache files are hidden from normal case selection. A
separate manual field accepts a bare case name such as `case14.m`, `case118.m`,
or `case9241pegase.m`; a nonempty manual value overrides the selected cached
case.

The landing page includes a compact, collapsible MATPOWER acknowledgement beside
the case inputs. It distinguishes Sparlectra from MATPOWER, provides links to
the MATPOWER project, its citation guidance, and the standard 2011 paper DOI,
and notes that ACTIVSg, PEGASE, RTE, and other case files may request additional
case-specific citations in their file headers.

A missing bare `.m` case name is resolved in the user Web UI `data/mpower` cache
through the standard MATPOWER download helper and remains the executed source.
Web UI/service PowerFlow runs do not automatically replace selected `.m` files
with generated `.jl` cache files. If a user manually submits a generated `.jl`
file from the Web UI MATPOWER cache and a matching `.m` source exists, the
service resolves back to the `.m` file and records that the generated cache was
bypassed. If the matching `.m` source is missing, the request is rejected with a
clear validation error. Large generated `.jl` MATPOWER cases can fail while
Julia and SparseArrays load literal data, before Sparlectra network construction
or Newton iterations begin, so `.m` remains the canonical Web UI execution
source.

Existing absolute and relative `.m` or `.jl` paths remain supported. A missing
input containing a path separator is rejected instead of downloaded, and URL
input is not accepted. The browser cannot select runtime directories. The
read-only configuration path is the provisioned user file or the explicit file
passed at startup.

## Run artifacts and output modes

The **Logfile output mode** is forwarded through the form, service request, API
configuration override, and `run.log` writer. `classic` keeps the standard
result report plus a compact API timing/status summary. `full` adds a marked
**Full run details** section with the effective typed configuration, artifact
choices, and available status diagnostics. The summary records
`solver_time`, `representative_time`, iterations, final mismatch, and final
outcome where available. Benchmark median and sample count appear when
benchmark mode is enabled.

The **Performance timing** control accepts `off`, `compact`, or `full` and
writes `performance.log` when enabled. It describes phases of one request,
unlike `benchmark.enabled`, which measures repeated solves and median timing.
Service runs can include request parsing and case resolution; API phases include
configuration, case loading/network construction/solve, postprocessing when
separately available, artifact writing, solver time, and total time. `full`
also includes available internal profile entries.

The **Run diagnostics** checkbox writes `diagnose.log` after the PowerFlow
result is available. It reuses the existing Q-limit event, PV-limit, and final
limit-validation printers. A diagnostic exception is contained and recorded in
that file without changing a successful PowerFlow result. Older run directories
can still contain `diagnose.txt`; the artifact viewer continues to list and
download that legacy filename.

The **Export detailed result CSV files** checkbox is off by default because
large networks can produce large files. When enabled for a successful run, it
writes Excel-friendly UTF-8 artifacts:

- `bus_voltages_complex.csv` contains one row per bus, including `vm_pu`,
  `va_deg`, numeric rectangular components `v_re` and `v_im`, a readable
  `v_complex` value, nominal/actual voltage, generation, load, Q-limit, and
  control columns.
- `branch_flows.csv` contains physical branch rows with active/reactive power
  at both ends, losses, rating, status, and overload information.

The CSV files reuse the structured `ACPFlowReport` node and branch rows and,
like the log artifacts, are viewable and downloadable through the normal
artifact list.

The indented **CSV format** selector is subordinate to the opt-in export:

- `technical` (default) uses a comma delimiter, decimal point, and no
  thousands grouping.
- `excel_de` uses a semicolon delimiter, decimal comma, and thousands dot.
- `excel_us` uses a comma delimiter, decimal point, and thousands comma.

US-formatted numbers containing a thousands comma are quoted because comma is
also the field delimiter. Empty values and fields containing delimiters,
quotes, carriage returns, or line feeds follow the same CSV quoting rules in
all formats. The Excel-oriented formats write numeric fields in decimal
notation instead of exponent notation where practical. The readable `v_complex`
column follows the selected decimal notation, while `v_re` and `v_im` remain
separate numeric columns.

Excel may still warn about automatic conversions when opening CSV files
directly, especially if textual identifiers resemble scientific notation such
as `1E5`. Use Excel's **Data > From Text/CSV** import flow and select text
types for exact textual identifiers when that distinction matters. The
`technical` format remains the clean machine-readable default and does not add
Excel-specific text hints.

## PowerFlow input paths

| Help topic | Input | Guidance |
|---|---|---|
| `webui.casefile` | MATPOWER case file | Choose an available canonical `.m` file from the existing-case selector, or type a bare case name or existing local path in the separate manual field. A nonempty manual value takes precedence. Missing bare `.m` names may be downloaded into the server-owned `data/mpower` directory; generated MATPOWER `.jl` cache files in that directory are not user-selectable and are not preferred for execution. Missing path-like inputs and URLs are rejected. |
| `webui.config_file` | Configuration template file | Select a YAML configuration or `*.yaml.example` template discovered in `examples`. Form values create allowlisted per-run overrides, while the selected template remains unchanged. |
| `webui.output_root` | Output root directory | Configure this path when calling `start_sparlectra_webui`; the browser displays it read-only. The service creates its persistent run index and one subdirectory per run beneath this root. |

## Contextual help and documentation

Every editable PowerFlow form option includes a contextual help link. Help pages
include a **Back** button that uses the browser history to return to the existing
PowerFlow form, preserving the values entered before opening help. If no local
history entry is available, the button falls back to `/powerflow`. Help pages
load the matching section or option row from repository Markdown at request
time. Solver options use
[`powerflow_configuration.md`](powerflow_configuration.md), output and benchmark
options use [`performance_profiling.md`](performance_profiling.md), and path
fields use the table above. Explanatory option text is not copied into Julia
views or HTML templates; repository Markdown remains the single source of
truth.

The **Documentation** navigation link opens `/docs`, which lists selected
allowlisted pages under `docs/src`. Each `/docs/<page>` request resolves only a
registered page name, so arbitrary paths and traversal requests are rejected.
Markdown links to another allowlisted page are rewritten to local `/docs/...`
routes, including section fragments such as the
[start-mode options](powerflow_configuration.md#start-mode-options). External
HTTP and HTTPS links remain external; unknown or unsafe local paths are made
inert. This is a lightweight reader for local reference material, not a
replacement for the Documenter.jl site.

## Results and artifacts

Successful and failed runs both have a result page. It shows the run ID, schema
version, status, convergence and solution flags, iteration count, final
mismatch, reason/message fields, input paths, and output directory.
Aborted runs are listed distinctly in history and are never rendered as
successful. After an abort, the user can return to the form and submit new
inputs without restarting the server.

Artifact lists come from `list_powerflow_artifacts`. Artifact requests are
resolved by exact metadata name through `resolve_powerflow_artifact`; browser
input is never joined directly to a filesystem path. JSON, YAML, logs, CSV,
HTML, Markdown, and other text artifacts are displayed as escaped text in a
large, scrollable panel that preserves long lines for horizontal scrolling.
Help excerpts and full documentation pages also use wider content panels and
readable line spacing. Other files are
downloaded, and every artifact page also offers an explicit download response.
The text-artifact viewer uses 75–85 percent of the viewport height and a wider
page layout for practical inspection of long logs and configuration files.

## Persistent operation log

The Web UI appends support-oriented JSON Lines events to its user-writable `logs/webui_operations.jsonl`. This file is independent of individual
run directories and survives Web UI restarts that reuse the same output root.
It records key page opens, submissions and validation failures, asynchronous
run lifecycle changes, artifact views/downloads, abort requests, history
refreshes, deletions, shutdown requests, and enabled diagnostics or timing
modes. Static CSS and image requests are not user-action events.

Use the shared **Operation Log** navigation link to open the escaped text viewer
at `/webui/operation-log`, or download the JSONL file from that page for an
error report. Entries contain concise route, method, status, run/case/artifact,
message, and timing fields when available. Every event also records
`sparlectra_version` and a millisecond-precision UTC timestamp using
`yyyy-mm-ddTHH:MM:SS.sssZ`. They never include artifact
contents, local file contents, or complete configuration bodies. Logging is
best effort and cannot fail a normal Web UI request. After an append takes the
file above 10,000 valid JSONL entries, compaction atomically keeps the newest
1,000 valid entries and drops empty or malformed lines encountered during
compaction. When the current file reaches 10 MiB, it is replaced after being
retained as `webui_operations.jsonl.1`; this byte-size guard remains independent
of entry-count compaction. The viewer and download read the current file.

Development reports distinguish external **Verification limitations**, such as
an unavailable browser executable or an HTTP 403/proxy failure while installing
documentation dependencies, from genuinely unfinished implementation work.

## Persistent run history and run management

`start_sparlectra_webui` refreshes `powerflow_runs_index.json` beneath the
configured output root before it begins serving pages, so valid runs from an
earlier Julia process appear immediately. The **Refresh registry** button
remains available for a manual reload. Missing, corrupt, or unsafe entries are
reported or skipped without preventing valid runs from loading.

History is ordered newest first and shows a readable local date/time, run ID,
status text, status badge, solver summary fields, and actions. Green, yellow,
red, gray, and blue badges distinguish successful, warning/partial, failed,
unknown, and running states while retaining visible text for accessibility.
Older indexes without timestamps use the `result.json` modification time.

Each registered run has a **Delete** action, and **Delete all runs** removes all
safely registered runs for the current configured root. These actions update
the in-memory registry and persistent index as well as deleting the matching
run directories. Run IDs are validated, indexed paths must remain beneath the
configured root, and browser-submitted output roots are ignored; unrelated
files and directories are never deletion targets.

## Browser and application shutdown

All normal pages send a small local heartbeat after they load. With
`auto_shutdown_on_browser_close=true`, heartbeat expiry after
`browser_heartbeat_timeout_seconds` triggers best-effort shutdown after at
least one heartbeat has been received. A server started with
`open_browser=false` therefore remains running until a browser actually
connects, the **Stop Web UI** button is used, `close(server)` is called, or the
terminal receives `Ctrl+C`.

Browser-close detection is necessarily best effort: browser crashes, forced
process termination, or operating-system shutdown may prevent a final clean
lifecycle. `Ctrl+C` remains the fallback. If a port is still occupied, stop the
old Julia process or start the Web UI with a different `port` value.

## Current limitations

This first prototype is intentionally synchronous and local. It has no State
Estimation UI, authentication, public-server mode, background queue, live
progress stream, WebSockets, database, topology view, or advanced plotting.
It uses a compact Julia `Sockets` HTTP layer to avoid a heavy web-framework
dependency.

### Configuration check and refresh

The PowerFlow page includes explicit **Check configuration** and **Refresh configuration** actions for user YAML files. They exist because the package configuration template can gain new options over time while existing local files remain in place. A check performs a dry run only: it compares the selected configuration against `src/configuration.yaml.example`, reports missing keys, known deprecated aliases, duplicate YAML keys, and shows a refreshed YAML preview without writing.

Refresh is conservative and user-initiated. It preserves existing user values, adds missing keys with current template defaults, and may normalize known deprecated aliases in user files to canonical start-mode settings or legacy `matpower_*` Q-limit mode names. It never rewrites YAML during Web UI startup. When writing a server-local configuration file, Sparlectra first creates a timestamped backup next to the original file. If duplicate YAML keys are detected, refresh refuses to write so the file can be reviewed manually. Browser-uploaded or pasted YAML is not rewritten in place; the refreshed YAML is offered as a download instead.
