# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# file: app/src/SparlectraApp.jl
# purpose: the application layer of Sparlectra: the GUI-ready service API
#          (run_sparlectra_api, the local PowerFlow service with its run
#          index and artifacts), the loopback Web UI, the operation log and
#          the sysimage build. It sits on the Sparlectra library and is what
#          start_webui.jl, the launcher and a GUI integration load; a
#          library session (scripts, notebooks, the solver tests) never
#          loads it, which keeps the library's precompile to the solver.
module SparlectraApp

using Sparlectra
using AnalyticLoadFlow
using Dates
using Logging
using Markdown
using Printf
using Random
using SHA
using Sockets
using TOML
using UUIDs
using ZipArchives
using LinearAlgebra
using SparseArrays

# the application package directory (app/) and the checkout that holds the
# library next to it; the sysimage tooling and the launcher agree on both
const SPARLECTRA_APP_ROOT = normpath(joinpath(@__DIR__, ".."))

# Internals of the library the application layer reads or extends. They are
# imported by name so that the boundary stays visible: a name missing here
# fails at load time, never at the first click.
import Sparlectra:
  _flatten_config_values!,
  EnergyConsumer,
  ExternalNetworkInjection,
  IaMeas,
  ImagMeas,
  Isolated,
  LinearShuntCompensator,
  Load,
  MatpowerIO,
  PQ,
  PV,
  PflowMeas,
  PinjMeas,
  QflowMeas,
  QinjMeas,
  Slack,
  SynchronousMachine,
  _cached_control_label,
  _se_start_state,
  FetchMatpowerCase,
  _set_measurement_active,
  ACTIVE_SPARLECTRA_CONFIG,
  BusPowerComponents,
  CGMES_HVDC_MODE_VALUES,
  CGMES_START_VALUES_VALUES,
  CONTINGENCY_SCREENING_MODE_VALUES,
  CsvFormatRuntime,
  DEFAULT_SPARLECTRA_CONFIG_PATH,
  DISTRIBUTED_SLACK_P_MODE_VALUES,
  EXTERNAL_GRID_SOURCE_VALUES,
  ImportedCase,
  MATPOWER_AUTO_PROFILE_VALUES,
  MATPOWER_BUS_SHUNT_MODEL_VALUES,
  MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES,
  MATPOWER_DCLINE_MODE_VALUES,
  MATPOWER_PV_VOLTAGE_SOURCE_VALUES,
  MATPOWER_RATIO_VALUES,
  MATPOWER_SHIFT_UNIT_VALUES,
  OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES,
  OUTPUT_LOGFILE_RESULTS_VALUES,
  POWERFLOW_LINEAR_SOLVER_VALUES,
  POWERFLOW_SOLVER_VALUES,
  POWERFLOW_START_ANGLE_MODE_VALUES,
  POWERFLOW_START_VOLTAGE_MODE_VALUES,
  PowerFlowAborted,
  SPARLECTRA_ROOT,
  SparlectraVersion,
  TRANSFORMER_TAP_CHANGER_MODEL_VALUES,
  TRUST_REGION_STEP_MODE_VALUES,
  WRONG_BRANCH_DETECTION_VALUES,
  ZERO_INJECTION_SIGMA,
  _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT,
  _NO_BUS_CONTROL_FLAGS,
  _SE_DENSE_LINALG_MAX_N,
  _apply_config_net_parameters!,
  _apply_external_grid_config!,
  _branch_anomaly_diagnostics,
  _branch_kind_name,
  _branch_terminal_state,
  _bus_control_flag_cache,
  _bus_mrid,
  _bus_name_by_idx,
  _bus_power_component_cache,
  _bus_power_components,
  _bus_voltage_setpoints_from_prosumers,
  _capture_initial_residual_rows!,
  _cgmes_delivery_paths,
  _column_normalized,
  _condition_report_line,
  _condition_verdict,
  _config_file_hash,
  _config_resolve_reason,
  _config_source_report,
  _config_with_request_csv_format,
  _controller_counts,
  _copy_powerflow_with,
  _criticality_wii_tolerance,
  _csv_field,
  _declared_slack_feeder,
  _default0,
  _case_format_conflict,
  _detect_case_format,
  _detect_yaml_duplicate_keys,
  _dotted_config_set!,
  _dotted_config_value,
  _duplicate_measured_quantities,
  _effective_bus_name,
  _format_branch_anomaly_rows,
  _format_csv_number,
  _format_top_mismatch_rows,
  _hvdc_link_flow_rows,
  _import_sparlectra_net,
  _is_machine_transformer,
  _islandwise_failure_message,
  _jacobian_condest,
  _load_and_validate_config,
  _load_api_config,
  _machine_controllers,
  _machine_side,
  _merge_config_overrides,
  _original_branch_name,
  _original_bus_name,
  _print_qlimit_active_set_summary,
  _rectangular_mismatch_diagnostics,
  _residual_diagnostics,
  _resolve_detailed_csv_format,
  _run_sparlectra,
  _scf_int,
  _scf_measurement_provenance,
  _scf_source_reference,
  _se_config_with,
  _self_check_forced_overrides,
  _shunt_controllers,
  _sigma_max,
  _solver_elapsed_from_profile,
  _sparlectra_config_with,
  _tap_controllers,
  _tap_voltage_target_by_bus,
  _transformer_mrids,
  _write_auto_mode_decision_log,
  _write_namedtuple_csv,
  _yaml_dict_text,
  auto_pf_record,
  case_config_path,
  comment,
  control_status,
  format_tolerance_physical,
  import_case,
  load_case_config,
  matpower_import_auto_profile,
  options_type,
  rectangular_pf_status,
  resolve_config,
  result_csv_format,
  scf_json_parse,
  snapshotPVQLimits!,
  sparlectra_arm_abort_token!,
  sparlectra_disarm_abort_token!,
  write_case_config,
  write_csv_row_direct!,
  write_matpower_import_auto_profile,
  write_result_csv

# --- service layer: run services and their metadata and artifact plumbing
include("api/api_types.jl")
include("api/serialization.jl")
include("api/artifacts.jl")
include("api/run_metadata.jl")
include("api/run_api.jl")
include("api/powerflow_service.jl")
include("api/run_short_circuit_service.jl")
include("api/run_import_analysis_service.jl")
include("api/run_contingency_service.jl")
include("api/run_state_estimation_service.jl")
include("api/se_measurement_generator.jl")

# --- the loopback Web UI and the sysimage build
include("webui/webui.jl")
include("sysimage_builder.jl")

# one export per line: the list stays readable, and every parser (the
# editor's included) agrees on where a statement ends
export run_sparlectra_api  # Stable non-interactive backend contract for GUI/API integrations.
export run_fixed_reference_self_check  # Evaluate mismatch at a case's own stored VM/VA, no corrective Newton step.
export SparlectraApiResult  # Structured API run status, numerical metadata, and artifacts.
export SparlectraApiArtifact  # Explicit metadata for generated API artifacts.
export collect_sparlectra_api_artifacts  # Discover generated files without filename assumptions.
export POWERFLOW_RUN_INDEX_FILENAME  # Persistent local PowerFlow run-index filename.
export start_powerflow_run  # Start and persist a local PowerFlow service run.
export load_powerflow_run_index  # Load the persistent run index from an output root.
export list_powerflow_runs  # List indexed runs and their disk availability.
export refresh_powerflow_run_registry!  # Recover the in-process registry from disk.
export delete_powerflow_run  # Safely delete one registered run beneath an output root.
export delete_all_powerflow_runs  # Safely delete all registered runs beneath an output root.
export get_powerflow_result  # Look up serialized run metadata by run ID.
export list_powerflow_artifacts  # List run artifacts by run ID.
export resolve_powerflow_artifact  # Safely resolve a run artifact by metadata name.
export default_webui_output_root  # Return the user-writable default Web UI output directory.
export default_webui_config_path  # Return the provisioned Web UI configuration path.
export default_webui_case_cache_dir  # Return the user-writable Web UI case cache.
export default_webui_operation_log_path  # Return the user-writable Web UI operation-log path.
export start_sparlectra_webui  # Start the loopback-only local PowerFlow Web UI.
export buildSysimage  # One-call sysimage build (10-20 min, see docstring).
export to_dict  # Convert API results and artifacts to dictionaries.
export to_namedtuple  # Convert API results to named tuples.
export to_json  # Serialize API results as JSON.
export to_yaml  # Serialize API results as YAML.

include("precompile.jl")

end # module SparlectraApp
