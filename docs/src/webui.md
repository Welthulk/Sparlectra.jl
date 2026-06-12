# Local PowerFlow Web UI

Sparlectra includes a small, local-first browser interface above the existing
PowerFlow service. The Web UI contains presentation, form parsing, and route
handling only; numerical execution continues through `start_powerflow_run` and
`run_sparlectra_api`.

## Start the server

```julia
using Sparlectra

server = start_sparlectra_webui(
    host = "127.0.0.1",
    port = 8080,
    output_root = "results/powerflow_service",
    open_browser = true,
    auto_shutdown_on_browser_close = true,
    browser_heartbeat_timeout_seconds = 15.0,
    warmup = true,
    warmup_store_result = false,
)
```

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
`examples/exp_webui_powerflow.jl` enables it for the local browser workflow.
Warm-up runs asynchronously and fails softly. With
`warmup_store_result=false`, it uses a temporary directory and creates no
normal run-history entry. The bundled `data/webui/warmup_case3.jl` file is an
original Sparlectra-owned synthetic three-bus MATPOWER-compatible Julia case,
not a derivative of an external MATPOWER case.

## Starting a PowerFlow run

The start page accepts:

- a typed MATPOWER case name or an existing local case path;
- a Sparlectra configuration template file;
- read-only information showing the server's configured output-root directory;
- PowerFlow tolerance and maximum iterations;
- autodamping and its minimum factor;
- Q-limit handling;
- wrong-branch detection mode;
- angle and voltage start modes;
- logfile result mode;
- single-run performance timing detail;
- optional post-run diagnostics; and
- benchmark enablement, samples, and seconds.

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are submitted. The output root is
chosen only through `start_sparlectra_webui(; output_root=...)`, is displayed as
read-only information, and cannot be overridden by a submitted browser field.
The page does not offer a generic YAML editor and never modifies the selected
template. The service creates an `effective_config.yaml` artifact for each run.

The Web UI resolves the Sparlectra application directory from the process start
directory. It supports starting directly in the repository or from a parent
directory containing `Sparlectra` or `Sparlectra.jl` (for example
`C:\Users\scud\.julia\dev`). MATPOWER `.m` and generated `.jl` files found in
`Sparlectra/data/mpower` are shown in an explicit existing-case selector. A
separate manual field accepts a bare case name such as `case14.m`, `case118.m`,
or `case9241pegase.m`; a nonempty manual value overrides the selected local
case. YAML files and `*.yaml.example` templates found in `Sparlectra/examples`
are offered in a separate dropdown.

A missing bare `.m` or `.jl` case name is resolved in the server-owned
`data/mpower` directory through the standard MATPOWER download helper. For an
`.m` case, Sparlectra generates a Julia `.jl` representation and uses it for the
run when generation succeeds. Existing generated `.jl` files are preferred on
later runs to avoid repeatedly parsing the `.m` source. If generation fails but
the `.m` file remains readable, the run falls back to that `.m` file.

Existing absolute and relative `.m` or `.jl` paths remain supported. A missing
input containing a path separator is rejected instead of downloaded, and URL
input is not accepted. The browser cannot select the case download directory or
the run output root. If the `examples` directory contains no supported
configuration file, the configuration field falls back to the package's
default `src/configuration.yaml.example`.

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

The **Run diagnostics** checkbox writes `diagnose.txt` after the PowerFlow
result is available. It reuses the existing Q-limit event, PV-limit, and final
limit-validation printers. A diagnostic exception is contained and recorded in
that file without changing a successful PowerFlow result. Both new files are
viewable and downloadable through the normal artifact list.

## PowerFlow input paths

| Help topic | Input | Guidance |
|---|---|---|
| `webui.casefile` | MATPOWER case file | Choose an available `.m` or `.jl` file from the existing-case selector, or type a bare case name or existing local path in the separate manual field. A nonempty manual value takes precedence. Missing bare names may be downloaded into the server-owned `data/mpower` directory; missing path-like inputs and URLs are rejected. Generated `.jl` cases are preferred for execution. |
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

Artifact lists come from `list_powerflow_artifacts`. Artifact requests are
resolved by exact metadata name through `resolve_powerflow_artifact`; browser
input is never joined directly to a filesystem path. JSON, YAML, logs, CSV,
HTML, Markdown, and other text artifacts are displayed as escaped text in a
large, scrollable, pre-wrapped panel. Help excerpts and full documentation pages
also use wider content panels and readable line spacing. Other files are
downloaded, and every artifact page also offers an explicit download response.
The text-artifact viewer uses 75–85 percent of the viewport height and a wider
page layout for practical inspection of long logs and configuration files.

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
