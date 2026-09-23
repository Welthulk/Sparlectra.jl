```@meta
CurrentModule = SparlectraApp
```

# API and Services Reference

## Public API

```@autodocs
Modules = [SparlectraApp]
Pages = [
  "app/src/api/api_types.jl",
  "app/src/api/artifact_registry.jl",
  "app/src/api/artifacts.jl",
  "app/src/api/powerflow_service.jl",
  "app/src/api/run_api.jl",
  "app/src/api/run_bus_powers_export.jl",
  "app/src/api/run_contingency_service.jl",
  "app/src/api/run_csv_exports.jl",
  "app/src/api/run_diagnostic_artifacts.jl",
  "app/src/api/run_failures.jl",
  "app/src/api/run_finalization.jl",
  "app/src/api/run_import_analysis_service.jl",
  "app/src/api/run_index.jl",
  "app/src/api/run_lifecycle_metadata.jl",
  "app/src/api/run_matpower_artifacts.jl",
  "app/src/api/run_metadata.jl",
  "app/src/api/run_self_check.jl",
  "app/src/api/run_short_circuit_service.jl",
  "app/src/api/run_state_estimation_service.jl",
  "app/src/api/se_measurement_generator.jl",
  "app/src/api/serialization.jl",
  "app/src/api/service_json.jl",
  "app/src/api/webui_jobs.jl",
]
Public = true
Private = false
```

## Application entry points

The Web UI server and the build tools live outside the reference page set
(`app/src/webui/`, `app/src/sysimage_builder.jl`); their exported entry points render here.

```@docs
SparlectraApp.start_sparlectra_webui
SparlectraApp.default_webui_output_root
SparlectraApp.default_webui_config_path
SparlectraApp.default_webui_case_cache_dir
SparlectraApp.default_webui_operation_log_path
SparlectraApp.buildSysimage
```

## Internals

```@autodocs
Modules = [SparlectraApp]
Pages = [
  "app/src/api/api_types.jl",
  "app/src/api/artifact_registry.jl",
  "app/src/api/artifacts.jl",
  "app/src/api/powerflow_service.jl",
  "app/src/api/run_api.jl",
  "app/src/api/run_bus_powers_export.jl",
  "app/src/api/run_contingency_service.jl",
  "app/src/api/run_csv_exports.jl",
  "app/src/api/run_diagnostic_artifacts.jl",
  "app/src/api/run_failures.jl",
  "app/src/api/run_finalization.jl",
  "app/src/api/run_import_analysis_service.jl",
  "app/src/api/run_index.jl",
  "app/src/api/run_lifecycle_metadata.jl",
  "app/src/api/run_matpower_artifacts.jl",
  "app/src/api/run_metadata.jl",
  "app/src/api/run_self_check.jl",
  "app/src/api/run_short_circuit_service.jl",
  "app/src/api/run_state_estimation_service.jl",
  "app/src/api/se_measurement_generator.jl",
  "app/src/api/serialization.jl",
  "app/src/api/service_json.jl",
  "app/src/api/webui_jobs.jl",
]
Public = false
Private = true
```
