# API and service layer (developer notes)

## Purpose

This directory contains the non-interactive service surface on top of the
framework: `run_sparlectra_api` and its result contract, the PowerFlow run
registry and artifact handling used by the Web UI, and the per-domain service
runs (contingency, short circuit, state estimation, import analysis). The
numerical work always goes through `run_sparlectra`; nothing here duplicates
solver logic.

## File responsibilities

| File | Responsibility | Key functions |
|---|---|---|
| `api_types.jl` | Typed result surface | `SparlectraApiResult`, `SparlectraApiArtifact` |
| `run_api.jl` | The orchestrator: one API run from request to result file | `run_sparlectra_api` |
| `config_overrides.jl` | GUI/API override allowlist and validation | `validate_gui_config_overrides`, `GUI_EDITABLE_CONFIG_KEYS` |
| `run_metadata.jl` | Result and metadata assembly | `_api_result`, `_finalize_api_result` |
| `run_failures.jl` | Failed-result construction with stable reasons | `_api_failure`, `_api_execution_failure` |
| `run_finalization.jl` | Output capture, phase summaries, resolved runtime options | `_capture_run_output`, `_finalize_service_timings!` |
| `artifacts.jl` | Artifact discovery and classification | `collect_sparlectra_api_artifacts` |
| `artifact_registry.jl` | Run-id validation and path-contained artifact access | `list_powerflow_artifacts`, `resolve_powerflow_artifact` |
| `run_index.jl` | Persistent run index (`powerflow_runs_index.json`) | `load_powerflow_run_index`, `list_powerflow_runs` |
| `run_lifecycle_metadata.jl` | Stable Web UI-visible lifecycle keys | `_build_success_lifecycle_metadata` |
| `run_csv_exports.jl` | Buffered and streaming detailed CSV export | `_write_namedtuple_csv` |
| `run_diagnostic_artifacts.jl` | Diagnostic artifacts (mismatch, narrative, Q-limit validation) | `_write_powerflow_diagnostics` |
| `run_matpower_artifacts.jl` | MATPOWER import artifacts (effective config, auto profile) | `_write_matpower_auto_profile_artifact` |
| `run_self_check.jl` | Fixed-reference self check | `run_fixed_reference_self_check` |
| `powerflow_service.jl` | In-memory run registry and case resolution for the Web UI | `start_powerflow_run`, `get_powerflow_result` |
| `webui_jobs.jl` | Web UI job registry, phases, cancellation | `get_active_webui_powerflow_job` |
| `run_contingency_service.jl` | Service N-1 contingency run | `_run_contingency_service` |
| `run_short_circuit_service.jl` | Service short-circuit run (CGMES and SCF) | `_run_short_circuit_service` |
| `run_state_estimation_service.jl` | Service state-estimation run | `_run_state_estimation_service` |
| `run_import_analysis_service.jl` | CGMES delivery analysis without a solve | `_run_import_analysis_service` |
| `serialization.jl` | Dependency-free to_dict/to_json/to_yaml transport | `to_dict`, `to_json`, `to_yaml` |
| `service_json.jl` | Minimal JSON parser and service failure helpers | `_parse_service_json!` family |

## Conventions

- Case input goes through `import_case` (`src/import/case_import.jl`); a
  service that refuses a format words its own refusal.
- Artifact paths are validated to stay inside the run directory
  (`_path_is_within`); never bypass the registry when handing paths out.
- Failure reasons are part of the contract (`"casefile_not_found"`,
  `"invalid_configuration"`, ...); add new ones deliberately and document them
  in `docs/src/programmatic_api.md`.

## Reference

Rendered API documentation: `docs/src/reference_api.md`. The user-facing
contract is described in `docs/src/programmatic_api.md`.
