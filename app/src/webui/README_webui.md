# Web UI (developer notes)

## Purpose

The local Web UI server: a dependency-free TCP/HTTP loop with routing, HTML
rendering, form handling, uploads, run management, in-app documentation, an
operation log, and sysimage status. It talks to the service layer in
`src/api/`; it never runs solvers itself.

## Include order

`webui.jl` is the module root and includes:

```julia
include("sysimage.jl")
include("options.jl")
include("forms.jl")
include("docs.jl")
include("operations.jl")
include("views.jl")
include("handlers.jl")
include("routes.jl")
```

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `webui.jl` | Server loop, request parsing, response writing, defaults | `start_sparlectra_webui`, `default_webui_output_root` |
| `routes.jl` | Request routing with URL decoding | `route_sparlectra_webui` |
| `handlers.jl` | Action handlers: imports, runs, downloads, settings | `handle_*` family |
| `views.jl` | HTML rendering for all pages | `_webui_escape` and the page renderers |
| `forms.jl` | Input resolution: case options, upload classification | `_webui_casefile_options` |
| `options.jl` | Central run-option metadata (`WEBUI_OPTION_SPECS`) | (data table) |
| `docs.jl` | In-app documentation from the markdown pages | `load_webui_help_excerpt` |
| `operations.jl` | Operation log: locked JSONL appends with retention | `webui_operation_log_path` |
| `sysimage.jl` | Sysimage location/metadata contract and build state | `sysimage_status`, `start_sysimage_build!` |

## Conventions

- Browser form semantics matter: disabled controls (including the submitting
  button itself) are dropped from the form entry list, and `onsubmit` runs
  before serialization. When a value "does not arrive", check the rendered
  HTML and the actual POST body before the server side.
- Escape everything rendered (`_webui_escape`); artifact paths are validated
  against the run directory before download.
- Keep the operation log high-level (no per-iteration entries); detailed
  timings belong in `performance.log`/`run.log` or structured metadata.
- The Web UI owns MATPOWER case registry, resolution, and caching; reuse it
  through the service path rather than building a second cache.
- These files have no reference page by design; the exported entry points
  render on `docs/src/reference_api.md` (Application entry points), and the
  file list is excused in `test/test_repository_hygiene.jl`.

## Reference

User-facing documentation: `docs/src/powerflow_service.md` and the in-app
help pages.
