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
)
```

Open `http://127.0.0.1:8080/powerflow` in a browser. The call returns a
`SparlectraWebUIServer` handle immediately. Stop it with `close(server)`.
Passing `open_browser=true` asks the operating system to open the URL, but is
optional.

For safety, the first prototype accepts loopback hosts only: `127.0.0.1`,
`localhost`, or `::1`. It is not intended for public or multi-user deployment.

The shared page header uses the existing Sparlectra documentation logo from
`docs/src/assets/logo.png`. The Web UI serves that single PNG through its local
asset route, and no additional branding configuration is required.

## Starting a PowerFlow run

The start page accepts:

- a local MATPOWER case file;
- a Sparlectra configuration template file;
- the output-root directory;
- PowerFlow tolerance and maximum iterations;
- autodamping and its minimum factor;
- Q-limit handling;
- wrong-branch detection mode;
- angle and voltage start modes;
- logfile result mode; and
- benchmark enablement, samples, and seconds.

Only keys in `GUI_EDITABLE_CONFIG_KEYS` are submitted. The page does not offer
a generic YAML editor and never modifies the selected template. The service
creates an `effective_config.yaml` artifact for each run.

The case and configuration fields are browser text fields rather than native
file uploads. Enter paths that are readable by the local Julia process.

## PowerFlow input paths

| Help topic | Input | Guidance |
|---|---|---|
| `webui.casefile` | MATPOWER case file | Enter a local `.m` case path readable by the Julia process. The Web UI passes the path to the existing PowerFlow service; it does not upload or modify the file. |
| `webui.config_file` | Configuration template file | Enter a readable Sparlectra YAML template. Form values create allowlisted per-run overrides, while the selected template remains unchanged. |
| `webui.output_root` | Output root directory | Enter the directory beneath which the service creates its persistent run index and one subdirectory per run. |

## Contextual help and documentation

Every visible PowerFlow form option includes a contextual help link. Help pages
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

## Persistent run history

The history page reads `powerflow_runs_index.json` beneath the selected output
root. Use **Refresh registry** to call `refresh_powerflow_run_registry!` and
make valid runs from an earlier Julia process available by run ID. Available
runs can then be opened through the normal result and artifact views.

## Current limitations

This first prototype is intentionally synchronous and local. It has no State
Estimation UI, authentication, public-server mode, background queue, live
progress stream, WebSockets, database, topology view, or advanced plotting.
It uses a compact Julia `Sockets` HTTP layer to avoid a heavy web-framework
dependency.
