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
stored_result = get_powerflow_result(run_id)
artifacts = list_powerflow_artifacts(run_id)
result_json = resolve_powerflow_artifact(run_id, "result.json")
```

Each run receives a unique ID before execution and writes to a GUI-friendly
location:

```text
output_root/
  <run_id>/
    run.log
    result.json
    effective_config.yaml
    ...
```

The service registry is intentionally local and in-process. Generated files are
stored on disk, while `run_id` lookup state lasts for the current Julia process.
This keeps route handlers thin and the service testable without a browser or a
long-running server.

## Service boundary

- [`start_powerflow_run`](@ref) validates the dictionary-like request, chooses
  the run ID and directory, invokes the programmatic API, and returns a
  transport-safe result dictionary.
- [`get_powerflow_result`](@ref) returns serialized run metadata by run ID.
- [`list_powerflow_artifacts`](@ref) returns metadata from the API artifact
  model rather than guessing generated filenames.
- [`resolve_powerflow_artifact`](@ref) resolves only an exact artifact metadata
  name belonging to the selected run. It rejects absolute paths, Windows-style
  arbitrary paths, traversal components, missing artifacts, and paths that
  resolve outside the run directory.

Public service failures are dictionaries containing `status`, `success`,
`reason`, and `message`. Lookup functions use reasons such as `run_not_found`,
`artifact_not_found`, and `unsafe_artifact_name`, which can be mapped directly
to future HTTP responses.

See `examples/exp_powerflow_service.jl` for a runnable local example.
