# Local PowerFlow Service

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

- [`start_powerflow_run`](@ref) validates the dictionary-like request, chooses
  the run ID and directory, invokes the programmatic API, registers the result,
  and updates the persistent index. A trusted caller such as the Web UI can
  provide the server-owned MATPOWER case directory through the function
  keyword; bare `.m` or `.jl` names are then resolved with `ensure_casefile`
  using `to_jl=true`. Existing or generated `.jl` cases are preferred, while a
  readable `.m` case remains the fallback if conversion fails. Browser form
  values never control this directory, missing path-like inputs are not
  downloaded, and URLs are rejected.
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

This layer intentionally introduces no HTTP routes, Genie.jl server, browser
GUI, authentication, or database. See `examples/exp_powerflow_service.jl` for a
runnable local example.
