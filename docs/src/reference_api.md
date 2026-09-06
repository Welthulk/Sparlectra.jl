# API and Services Reference

## Public API

```@autodocs
Modules = [Sparlectra]
Pages = [
  "src/api/api_types.jl",
  "src/api/artifact_registry.jl",
  "src/api/artifacts.jl",
  "src/api/powerflow_service.jl",
  "src/api/run_api.jl",
  "src/api/run_contingency_service.jl",
  "src/api/run_csv_exports.jl",
  "src/api/run_diagnostic_artifacts.jl",
  "src/api/run_failures.jl",
  "src/api/run_finalization.jl",
  "src/api/run_import_analysis_service.jl",
  "src/api/run_index.jl",
  "src/api/run_lifecycle_metadata.jl",
  "src/api/run_matpower_artifacts.jl",
  "src/api/run_metadata.jl",
  "src/api/run_self_check.jl",
  "src/api/run_short_circuit_service.jl",
  "src/api/run_state_estimation_service.jl",
  "src/api/serialization.jl",
  "src/api/service_json.jl",
  "src/api/webui_jobs.jl",
]
Public = true
Private = false
```

## Application entry points

The Web UI server and the build tools live outside the reference page set
(`src/webui/`, `src/build/`); their exported entry points render here.

```@docs
start_sparlectra_webui
default_webui_output_root
default_webui_config_path
default_webui_case_cache_dir
default_webui_operation_log_path
buildSysimage
```

## Internals

```@autodocs
Modules = [Sparlectra]
Pages = [
  "src/api/api_types.jl",
  "src/api/artifact_registry.jl",
  "src/api/artifacts.jl",
  "src/api/powerflow_service.jl",
  "src/api/run_api.jl",
  "src/api/run_contingency_service.jl",
  "src/api/run_csv_exports.jl",
  "src/api/run_diagnostic_artifacts.jl",
  "src/api/run_failures.jl",
  "src/api/run_finalization.jl",
  "src/api/run_import_analysis_service.jl",
  "src/api/run_index.jl",
  "src/api/run_lifecycle_metadata.jl",
  "src/api/run_matpower_artifacts.jl",
  "src/api/run_metadata.jl",
  "src/api/run_self_check.jl",
  "src/api/run_short_circuit_service.jl",
  "src/api/run_state_estimation_service.jl",
  "src/api/serialization.jl",
  "src/api/service_json.jl",
  "src/api/webui_jobs.jl",
]
Public = false
Private = true
```
