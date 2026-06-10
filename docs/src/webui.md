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

## Results and artifacts

Successful and failed runs both have a result page. It shows the run ID, schema
version, status, convergence and solution flags, iteration count, final
mismatch, reason/message fields, input paths, and output directory.

Artifact lists come from `list_powerflow_artifacts`. Artifact requests are
resolved by exact metadata name through `resolve_powerflow_artifact`; browser
input is never joined directly to a filesystem path. JSON, YAML, logs, CSV,
HTML, Markdown, and other text artifacts are displayed as escaped text. Other
files are downloaded, and every artifact page also offers an explicit download
response.

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
