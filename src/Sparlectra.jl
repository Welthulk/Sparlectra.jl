# Copyright 2023–2026 Udo Schmitz
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

# Author: Udo Schmitz (https://github.com/Welthulk)

# Naming Conventions:
# The project follows the Julia Naming Conventions for the most part,
# but it's important to note that the naming convention for functions might deviate.
# In this module, functions are written in CamelCase with a lowercase initial letter.

# file: src/Sparlectra.jl
# purpose: main package module: dependency loading, package-wide constants,
#          include order of all source files, and the public export surface
#! format: off

module Sparlectra

using AnalyticLoadFlow
using BenchmarkTools
using Dates
using LinearAlgebra
using Logging
using Printf
using SparseArrays
using TOML

const SPARLECTRA_ROOT = normpath(joinpath(@__DIR__, ".."))
const MPOWER_DIR = normpath(joinpath(SPARLECTRA_ROOT, "data", "mpower"))

"""
The square root of three, the line-to-line factor of three-phase power.
"""
const Wurzel3 = 1.7320508075688772

function _read_project_version()::VersionNumber
    project_file = joinpath(SPARLECTRA_ROOT, "Project.toml")

    isfile(project_file) ||
        error("Project.toml wurde nicht gefunden: $project_file")

    project_data = TOML.parsefile(project_file)

    haskey(project_data, "version") ||
        error("In $project_file fehlt der Eintrag 'version'.")

    return VersionNumber(project_data["version"])
end

# SparlectraVersion is baked in at precompile time; declare Project.toml as a
# precompile dependency so a version bump alone invalidates the cache and
# version() cannot report a stale number.
Base.include_dependency(joinpath(SPARLECTRA_ROOT, "Project.toml"))
const SparlectraVersion = _read_project_version()

"""
    version() -> VersionNumber

The version of the loaded Sparlectra package.
"""
version() = SparlectraVersion

"""
Supertype of every branch-like element (lines, transformers).
"""
abstract type AbstractBranch end

const _MODULE_DOC = """
    Sparlectra $(SparlectraVersion)

Sparlectra is a Julia package for the calculation of electrical networks.
It is designed to be used in the context of power system analysis and optimization.

- GitHub Repository: https://github.com/welthulk/Sparlectra.jl
- Website: https://welthulk.github.io/Sparlectra.jl

"""
@doc _MODULE_DOC Sparlectra
export
  version,                                # Return the package version.

  # constants
  Wurzel3,
  ComponentTyp,

  # classes

  # Component
  AbstractComponent,
  Component,
  ImpPGMComp,
  ImpPGMComp3WT,

  # Node
  Node,

  # Line
  ACLineSegment,

  # Trafo
  TrafoTyp,
  AbstractTapChangerModel,
  PowerTransformerTaps,
  PhaseTapChangerModel,
  TapTablePoint,
  PowerTransformerWinding,
  PowerTransformer,
  TransformerModelParameters,

  # ProSumer
  ProSumer,
  AbstractVoltageDependentController,
  PiecewiseLinearCharacteristic,
  QUController,
  PUController,
  VoltageAdjustConfig,

  # Tap controller
  AbstractOuterController,
  AbstractControlState,
  AbstractControlUpdate,
  PowerTransformerControl,
  MachineVoltageControl,                  # Remote voltage control via machine Q.
  ShuntVoltageControl,                    # SVC-style variable-shunt voltage control.

  # Branch
  AbstractBranch,
  Branch,
  BranchModel,
  BranchFlow,
  getBranchFlow,
  setBranchFlow!,
  getBranchLosses,
  setBranchLosses!,
  setBranchStatus!,
  setBranchTerminalStatus!,               # Open or close individual branch terminals (r0.9.10, one-sided open branches).

  # Link
  BusLink,
  setLinkStatus!,
  setLinkFlow!,
  setLinkCurrent!,

  # Shunt
  Shunt,

  # Net
  Net,

  # functions

  # performance_profile.jl

  # yamlparams.jl
  parse_yaml_scalar,
  load_yaml_dict,
  merge_yaml_dict!,
  as_bool,
  as_int_vector,

  # configuration.jl
  StartModeConfig,
  QLimitConfig,
  PowerFlowConfig,
  ControlConfig,
  ObservabilityConfig,
  StateEstimationConfig,
  MatpowerImportConfig,
  MatpowerExportConfig,                   # MATPOWER export settings (write_solution).
  ModelConfig,                            # Model-construction settings (shunt model, tap-changer model, auto profile, cache, preallocation).
  PerformanceConfig,
  BenchmarkConfig,
  RuntimeConfig,
  ParallelRuntimeConfig,                  # runtime.parallel.* settings (enabled, max_tasks, min_work_items).
  parallel_max_tasks,                     # Resolve max_tasks (auto = Threads.nthreads()) to a task count.
  DiagnosticsConfig,
  OutputConfig,
  SparlectraConfig,
  load_sparlectra_config,                 # Load a configuration file without activating it.
  load_sparlectra_config!,                # Load and activate a configuration file.
  refresh_sparlectra_config_file,        # Check or explicitly refresh a user YAML configuration file.
  refresh_sparlectra_config_text,        # Dry-run refresh for uploaded/browser YAML text.
  set_sparlectra_config!,                 # Replace the active global configuration.
  active_sparlectra_config,               # Return the currently active configuration.
  powerflow_config,                       # Access the active power-flow settings.
  matpower_import_config,                 # Access MATPOWER import settings.
  matpower_export_config,                 # Access MATPOWER export settings.
  model_config,                           # Access model-construction settings.
  configured_matpower_cases,              # Resolve configured MATPOWER batch order.
  state_estimation_config,
  diagnostics_config,
  output_config,
  performance_config,
  benchmark_config,
  runtime_config,
  print_effective_config,                 # Print the resolved effective configuration.
  configuration_path_from_inputs,         # Resolve the configuration file path.
  run_matpower_case,                      # Run a MATPOWER case through the high-level workflow.
  run_synthetic_tiled_grid_pf_perf,       # Run synthetic-grid PF performance tests.
  run_voltage_dependent_control_demo,

  # synthetic_grids.jl
  synthetic_tiled_grid_bus_index,
  build_synthetic_tiled_grid_net,         # Build a synthetic benchmark network.
  build_tiled_grid_net,

  # BusData
  BusData,
  getBusData,
  getBusTypeVec,
  countNodes,
  map_NR_voltage_to_net!,                 # Copy NR voltage results back into the network.
  buildVoltageVector_from_busVec,         # Build a voltage vector from bus data.

  # Compomnent
  toComponentTyp,
  getCompName,
  getCompID,

  # Transformers
  getSideNumber2WT,
  getWinding2WT,
  calcTransformerRatio,                   # Compute transformer off-nominal ratio.
  create2WTRatioTransformerNoTaps,
  create3WTWindings!,
  getTrafoImpPGMComp,
  getWT3AuxBusID,
  isPerUnit_RXGB,
  getWindingRatedS,
  getTrafoRXBG,
  getTrafoRXBG_pu,

  # Nodes
  setVmVa!,                               # Set voltage magnitude and angle on a node.
  addShuntPower!,
  addLoadPower!,
  addGenPower!,
  getNodeVn,
  isSlack,
  isPVNode,
  isPQNode,
  isIsolated,
  toNodeType,
  setNodeType!,
  getNodeType,
  busComparison,
  toString,

  # Branch
  getBranchIdx,
  calcBranchYser,
  calcBranchYserBase,                     # Series admittance from the physical base impedance (#329).
  restoreBaseImpedances!,                 # Reset live branch impedance to the physical base, discarding a series-FACTS operating point (#329).
  calcBranchYshunt,
  calcBranchRatio,
  calcAdmittance,

  # Link
  BusLink,
  setLinkStatus!,
  setLinkFlow!,
  setLinkCurrent!,

  # Shunt
  getGBShunt,
  getPQShunt,
  updatePQShunt!,

  # ACLineSegment
  get_line_parameters,
  isLinePIModel,
  getLineRXBG,
  getLineRXBG_pu,

  # ProSumer
  isSlack,
  isGenerator,
  isAPUNode,
  isRegulating,
  setQGenReplacement!,
  getQGenReplacement,
  toProSumptionType,
  updatePQ!,
  getPosumerBusIndex,
  setPQResult!,
  evaluate_characteristic,                # Evaluate a voltage-dependent characteristic.
  evaluate_controller,                    # Evaluate a voltage-dependent controller.
  make_characteristic,                    # Create a controller characteristic.
  has_qu_controller,
  has_pu_controller,

  # Network
  addBus!,                                # Add a bus to the network.
  addShunt!,                              # Add a shunt element.
  addACLine!,                             # Add an AC line.
  addPIModelACLine!,                      # Add an AC line as a PI model.
  add2WTrafo!,                            # Add a two-winding transformer.
  addPIModelTrafo!,                       # Add a transformer PI model.
  applyTapNameplate!,                     # Declare a transformer's tap-changer nameplate (step size, position band, live position) - the one path importers and hand-built nets share.
  addProsumer!,                           # Add generator/load prosumer data.
  addExternalGrid!,                       # Add an external grid: slack-or-source PF side + native SC feeder data (#299).
  convertSlackToExternalGrid!,            # Replace the marked slack by a non-ideal external-grid source (#299).
  lockNet!,                               # Finalize network topology before solving.
  validate!,                              # Validate network consistency.
  hasBusInNet,
  addBusGenPower!,
  addBusLoadPower!,
  setNodeVoltage!,
  setNodeAngle!,
  getNetOrigBusIdx,
  geNetBusIdx,
  setNetBranchStatus!,
  getNetBranch,
  getNetBranchNumberVec,
  setTotalLosses!,
  getTotalLosses,
  getBusType,
  getEffectiveBusType,
  getBusProsumers,
  refreshBusTypesFromProsumers!,          # Rebuild effective bus types from prosumers.
  get_bus_vn_kV,
  get_vn_kV,
  updateBranchParameters!,
  hasShunt!,
  getShunt!,
  markIsolatedBuses!,
  ensureSlack!,                           # Auto-promote a reference when no slack is registered.
  setTotalBusPower!,
  setPVBusVset!,                          # Set PV voltage target.
  setQLimits!,                            # Set generator reactive-power limits.
  getNodeVm,
  distributeBusResults!,
  getTotalBusPower,
  getTotalLosses,
  buildVoltageVector,                     # Build complex bus-voltage vector.
  initialVrect,                           # Build the rectangular initial voltage vector.
  buildComplexSVec,                       # Build specified complex power vector.
  buildControlledSVec,                    # Build controlled specified power vector.
  has_voltage_dependent_control,
  addShuntMatpower!,
  normalize_bus_shunt_model,              # Normalize bus-shunt modeling mode.
  bus_shunt_totals_pu,
  log_bus_shunt_model,
  add2WTPIModelTrafo!,
  add3WTPiModelTrafo!,
  showNet,                                # Print a network summary.
  buildQLimits!,                          # Build Q-limit arrays for the solver.
  updateShuntPowers!,
  addLink!,
  setNetLinkStatus!,
  getNetLinks,
  calcLinkFlowsKCL!,                      # Compute bus-link flows from KCL.
  calcLinkFlowsSE!,                       # W2 link-flow allocation with link measurements (SE phase 3).
  se_view,                                # Static SE view: frozen controllers, shunt releases, link clusters, excluded measurements.
  print_se_view,                          # Formatted se_view report.
  collect_outer_controllers,
  run_control!,                           # Run outer-loop control workflow.
  latest_control_result,
  ControlRunResult,
  addPowerTransformerControl!,            # Add transformer control data.
  addTapController!,                      # Add a tap-controller wrapper.
  clearTapControllers!,
  get_bus_vm_pu,
  get_branch_p_from_to_mw,
  get_branch_q_from_to_mvar,
  buildTapControllerReportRows,
  printTapControllerSummary,
  addMachineVoltageControl!,              # Add a machine remote voltage controller.
  clearMachineControllers!,
  addShuntVoltageControl!,                # Add an SVC-style shunt voltage controller.
  addSeriesReactanceControl!,             # Add a TCSC-like series-reactance controller on a line branch (#297).
  clearSeriesReactanceControllers!,       # Remove all TCSC/SSSC controllers, restoring each branch to its base impedance (#329).
  addUpfcControl!,                        # Register a UPFC: model=:quadrature SSSC+STATCOM composite (#325) or model=:full DC-link-coupled independent-P/Q model (#326).
  clearUpfcFullControllers!,              # Remove all full-UPFC controllers, restoring each branch to its base impedance (#329).
  addHvdcPairControl!,                    # Add a back-to-back HVDC pairing controller on two converter injections (#297).
  addHvdcLink!,                           # Register a Stage-0 HVDC link record on a hand-built net (r0.9.9 result reporting).
  HvdcPairControl,                        # Controller type of the HVDC pair (element rows, isa checks).
  clearHvdcPairControllers!,              # Remove all HVDC pair controllers from a net.
  printHvdcPairControllerSummary,         # Engineering summary of the registered HVDC pairs.
  applyConfiguredControllers!,            # Instantiate controllers declared under control.controllers (#305).
  printSeriesReactanceControllerSummary,  # Engineering-style summary of the registered TCSC controllers.
  printUpfcFullControllerSummary,         # Engineering-style summary of the registered full UPFC controllers (#326).
  printShuntVoltageControllerSummary,     # Engineering-style summary of the registered SVC/MSC shunt controllers.
  clearShuntControllers!,
  controllableElements,                   # Generic controllable-element records of all registered controllers.
  buildMachineControllerReportRows,
  printMachineControllerSummary,

  # remove_functions.jl
  removeBus!,                             # Remove a bus and dependent data.
  removeBranch!,                          # Remove a branch.
  removeACLine!,                          # Remove an AC line.
  removeTrafo!,                           # Remove a transformer.
  removeShunt!,                           # Remove a shunt.
  removeProsumer!,                        # Remove prosumer data.
  clearIsolatedBuses!,                    # Remove isolated imported buses.

  # import.jl
  createNetFromMatPowerFile,              # Import a MATPOWER case file as Net.
  DTFImporter,                            # Native DTF parser and Net builder.
  createNetFromDTFFile,                   # Import a legacy DTF file as Net.
  CGMESImporter,                          # Lean CGMES reader and diagnostics (Stage 0).
  summarizeCGMES,                         # Diagnose a CGMES delivery without building a Net.
  analyzeCGMES,                           # Explain why a delivery cannot import (missing dependencies, unresolved refs).
  createNetFromCGMES,                     # Import a CGMES delivery as Net (Stage 1).
  importCGMES,                            # CGMES import with full result (Net + SC data + report).
  compareWithSV,                          # Compare a solved CGMES import against its SV profile.
  shortCircuitCoverage,                   # Completeness of harvested short-circuit data (#277 input).
  writeCGMESFiles,                        # Export a Net as CGMES EQ+TP+SSH profile files.
  CGMESLineShortCircuit,                  # Optional per-line zero-sequence data for the CGMES export.
  cgmesLineShortCircuitData,              # Harvested zero-sequence line data as writeCGMESFiles input.
  printShortCircuitCoverage,              # Readable rendering of the coverage rows.
  runShortCircuit!,                       # IEC 60909 balanced Ik'' max/min per fault bus (#277).
  exportSCF,                              # Write a Net as a Sparlectra Case Format file (PGM-compatible, #342).
  net_to_scf,                             # Build the SCF root object of a Net without writing it (#342).
  importSCF,                              # Read a Sparlectra Case Format file into a Net (#342).
  scf_to_net,                             # Build a Net from a parsed SCF root object (#342).
  SCFCase,                                # Typed in-memory form of one SCF document (adapter stage 2, D1).
  read_scf_json,                          # Read a case file into its typed SCFCase form.
  write_scf_json,                         # Serialize an SCFCase to its canonical file bytes.
  build_net,                              # The one network constructor of the run path over an SCFCase (D2/D12).
  net_to_scfcase,                         # The typed case a Net exports as (same keywords as exportSCF).
  FormatAdapter,                          # Abstract input-format adapter contract (stage 3, D3).
  MatpowerAdapter,                        # MATPOWER adapter: convert_case(mpc) -> SCFCase (stage 3a).
  MatpowerAdapterOptions,                 # Adapter-scope options of the MATPOWER conversion.
  convert_case,                           # Adapter contract: source -> SCFCase.
  matpower_adapter_options,               # MatpowerAdapterOptions from an effective run configuration.
  DTFAdapter,                             # DTF adapter: convert_case(DTFCase) -> SCFCase (stage 3b).
  DTFAdapterOptions,                      # Adapter-scope options of the DTF conversion.
  dtf_adapter_options,                    # DTFAdapterOptions from an effective run configuration.
  CGMESAdapter,                           # CGMES adapter: configured import captured as SCFCase (stage 3c).
  CGMESAdapterOptions,                    # Adapter-scope options of the CGMES conversion.
  cgmes_adapter_options,                  # CGMESAdapterOptions from an effective run configuration.
  cgmes_enrich_case!,                     # Attach the mRID registry and per-component mRIDs to a typed case.
  PGMAdapter,                             # power-grid-model adapter: plain dataset through the shared pipeline (stage 3d).
  PGMAdapterOptions,                      # Contract options of the PGM conversion (none).
  PatchOp,                                # One scenario patch operation on an SCF component id (scenario task D1).
  Scenario,                               # Named, weighted, ordered patch list (D2).
  ScenarioSet,                            # Scenario vector plus N-1 mode and exclusions.
  ScenarioIndex,                          # SCF id to component class and internal index, from the typed case.
  validate_scenarios,                     # Load-time validation with scenario name and op index in every error.
  expand_scenarios,                       # Expand the N-1 modes through the existing generators.
  scf_case_scenarios,                     # The scenario set a typed case carries (scenarios or mapped contingencies).
  scenario_set_dict,                      # Document form of a scenario set.
  scenario_set_from_dict,                 # Read a sparlectra.scenarios block.
  scenario_set_from_contingencies,        # Map the legacy contingencies block onto the scenario model (D3).
  UndoLog,                                # Reversible record of one apply! (scenario task D4).
  apply!,                                 # Apply patch operations to a working copy with an undo log.
  restore!,                               # Replay the undo log in reverse; the copy returns to the base bitwise.
  runScenarios!,                          # Evaluate a scenario set on the engine (D9); N-1 status ops match runContingencies! exactly.
  ScenarioResult,                         # ContingencyResult surface plus screened flag and screening estimate (D8).
  scf_case_config,                        # The dotted configuration overrides a case file carries (#342).
  scf_is_case_config_key,                 # Whether a configuration key describes the case or the installation (#342).
  scf_case_studies,                       # The contingency and short-circuit study definitions of a case file (#342).
  scf_extra_names,                        # Component id to reference name, so study blocks can be resolved against a network (#342).
  scf_validate_dataset,                   # File-level validation of a parsed case file, without building the network (#342).
  scf_fault_nodes,                        # The buses a case file's PGM fault rows point at (#342).
  ShortCircuitResult,                     # Result type of runShortCircuit! (safety-flagged rows).

  # contingency.jl (N-1 batch, multi-core Phase 4)
  ContingencyCase,                        # One outage case (name, kind = :branch/:gen, element, weight).
  ContingencyResult,                      # Outcome of one case: convergence, envelope, overloads, shed load, weight.
  OverloadRecord,                         # One overloaded branch under a case (name + loading %).
  ContingencyReport,                      # Aggregate over a batch (counts, load shed, worst loading, top overloads).
  runContingencies!,                      # Evaluate a case batch on template copies, optionally threaded.
  generateN1Branches,                     # One case per in-service branch (transformer/voltage/rating/name filters).
  generateN1Generators,                   # One case per in-service generator (kind = :gen; Pg/name filters).
  generateContingenciesFromFOR001,        # Cases from imported MATPOWER FOR001 metadata.
  applyContingencyWeights,                # Attach per-case weights (by name) to a case list.
  readContingencyWeightsCSV,              # Read case name -> outage weight from a two-column CSV.
  printContingencyResults,                # Fixed-width result table (style of printShortCircuitResult).
  writeContingencyResultsCSV,             # Semicolon CSV writer for the result batch.
  buildContingencyReport,                 # Aggregate a result batch into a ContingencyReport.
  printContingencyReport,                 # Print the aggregate report summary block.
  NativeShortCircuitData,                 # Native SC source container on Net (filled by addExternalGrid!, #299).
  printShortCircuitResult,                # Pretty-printer for ShortCircuitResult (CSV via the service run).
  _createDict,

  # exportMatPower.jl
  writeMatpowerCasefile,                  # Export a network as MATPOWER case file.

  # equicircuit.jl
  calcComplexRatio,
  calcSkewAngleTap,
  calcTapImpedanceCorrectionFactor,       # Tap-changer impedance-correction factor (central model).
  calcTapCorrectedRX,                     # Apply the tap-changer impedance correction to R/X.
  calcRatioTapCorrection,                 # Ratio-tap correction factor for the :neutral_relative convention.
  calcRatioTapRange,                      # Ratio-terms tap range (tap_min, tap_max, tap_step).
  calcPhaseTapFraction,                   # CGMES phase-tap-changer n-n0 tap fraction.
  calcPhaseTapAngleRatio,                 # CGMES phase-tap-changer effective ratio/shift/regulating vector.
  calcPhaseTapReactance,                  # CGMES phase-tap-changer series-reactance dependence on tap angle.
  calcPhaseTapTable,                      # Tabular phase-tap-changer exact lookup (overrides formulas).
  calcNeutralU,
  electricalIslandComponents,             # Island decomposition counting closed links as connections.
  electricalIslandOfBus,                  # Bus -> electrical island number.
  createYBUS,                             # Build the network admittance matrix.
  adjacentBranches,
  toPU_RXBG,
  fromPU_RXBG,
  branchFlow_pu,                          # Compute branch flow in per unit.

  # nbi.jl

  # jacobian.jl
  runpf!,                                 # Run the classic Newton-Raphson PF solver.

  # jacobian_full.jl

  # condition_number.jl
  condestJacobian,                        # Hager 1-norm condition estimate for NR Jacobians.
  reportCondition,                        # Print condition estimate with digits-lost verdict.

  # limits.jl
  printQLimitLog,
  printPVQLimitsTable,                    # Print PV Q-limit diagnostics.
  printFinalLimitValidation,              # Print final Q-limit validation.
  qvCharacteristicViolations,             # Generators at a Q-limit on the wrong voltage side (non-physical solution).
  printQVCharacteristicCheck,             # Print that check; part of the result print and the run report.
  validate_q_limit_signs!,                # Validate Q-limit sign conventions.
  logQLimitHit!,
  lastQLimitIter,
  getQLimits_pu,
  logQLimitHit!,
  lastQLimitIter,
  resetQLimitLog!,                        # Clear Q-limit event log.
  pv_hit_q_limit,
  has_q_limits,
  active_set_q_limits!,                   # Apply active-set PV/PQ Q-limit logic.

  # losses.jl
  calcNetLosses!,                         # Compute total network losses.

  # results.jl
  ACPFlowReport,
  buildACPFlowReport,                     # Build structured AC PF report.
  printACPFlowResults,                    # Print AC PF result report.
  printProsumerResults,

  # powerflow_dc/ — standalone DC power flow (MATPOWER rundcpf/makeBdc equivalent).
  rundcpf!,                               # Solve standalone DC power flow; optional AC-NR seeding via seed_ac_start.
  DcPowerFlowConfig,                      # power_flow.dc.* configuration.
  DcPowerFlowReport,                      # Structured DC PF report (angles + lossless branch flows).
  buildDcPowerFlowReport,                 # Build structured DC PF report.
  printDcPowerFlowResults,                # Print DC PF result report.
  dc_pf_status,                           # Most recent DC PF status for a network.
  solve_dc_powerflow,                     # Core linear DC PF solve (theta, Pf, Pt, slack injection).
  assemble_dc_bbus,                       # MATPOWER makeBdc-equivalent B' + phase-shift injection assembly.

  # Configuration-driven framework runner.
  run_sparlectra,                         # Preferred public import/control/solve/output workflow.
  run_sparlectra_cases,                   # Run configured MATPOWER cases sequentially.
  run_acpflow,                            # Thin AC power-flow alias for run_sparlectra.
  SparlectraRunResult,                    # Stable typed framework-run result.
  run_sparlectra_api,                     # Stable non-interactive backend contract for GUI/API integrations.
  run_fixed_reference_self_check,         # Evaluate mismatch at a case's own stored VM/VA, no corrective Newton step.
  SparlectraApiResult,                    # Structured API run status, numerical metadata, and artifacts.
  SparlectraApiArtifact,                  # Explicit metadata for generated API artifacts.
  GUI_EDITABLE_CONFIG_KEYS,               # Controlled allowlist for GUI configuration overrides.
  validate_gui_config_overrides,          # Validate and nest dotted GUI configuration overrides.
  collect_sparlectra_api_artifacts,        # Discover generated files without filename assumptions.
  POWERFLOW_RUN_INDEX_FILENAME,           # Persistent local PowerFlow run-index filename.
  start_powerflow_run,                    # Start and persist a local PowerFlow service run.
  load_powerflow_run_index,               # Load the persistent run index from an output root.
  list_powerflow_runs,                    # List indexed runs and their disk availability.
  refresh_powerflow_run_registry!,        # Recover the in-process registry from disk.
  delete_powerflow_run,                   # Safely delete one registered run beneath an output root.
  delete_all_powerflow_runs,              # Safely delete all registered runs beneath an output root.
  get_powerflow_result,                   # Look up serialized run metadata by run ID.
  list_powerflow_artifacts,               # List run artifacts by run ID.
  resolve_powerflow_artifact,             # Safely resolve a run artifact by metadata name.
  default_webui_output_root,              # Return the user-writable default Web UI output directory.
  default_webui_config_path,              # Return the provisioned Web UI configuration path.
  default_webui_case_cache_dir,           # Return the user-writable Web UI case cache.
  default_webui_operation_log_path,       # Return the user-writable Web UI operation-log path.
  start_sparlectra_webui,                 # Start the loopback-only local PowerFlow Web UI.
  buildSysimage,                          # One-call sysimage build (10-20 min, see docstring).
  to_dict,                                # Convert API results and artifacts to dictionaries.
  to_namedtuple,                          # Convert API results to named tuples.
  to_json,                                # Serialize API results as JSON.
  to_yaml,                                # Serialize API results as YAML.
  run_matpower_case,                      # Run a MATPOWER case through the high-level workflow.

  # solver_core.jl
  calc_injections,                        # Compute complex bus injections.
  calc_currents,                          # Compute bus currents.
  solve_linear,
  solve_sparse_system,                    # Solve a sparse linear system.
  UmfpackReuseNewtonContext,              # Reusable UMFPACK lu! factorization context for the rectangular Newton step.
  solve_newton_factorized!,               # Analyze/refactor solve (UMFPACK lu! reuse) with shared fallback chain.
  build_pos_map,
  slack_elimination_indices,
  extract_bus_types_and_vset,             # Extract solver bus types and V targets.
  build_qload_pu,
  build_voltage_vector,                   # Build solver voltage vector.
  compute_sbus_and_totals,                # Build Sbus and aggregate P/Q totals.

  # measurements.jl / state_estimation.jl
  MeasurementType,
  Measurement,
  measurementStdDevs,
  measurementSigmaFloors,                 # Per-type sigma floors for relative-sigma measurement generation.
  generateMeasurementsFromPF,             # Generate synthetic SE measurements from PF.
  setMeasurementsFromPF!,                 # Replace measurements from PF results.
  addMeasurementNoise!,                   # Perturb the measurement values a net already holds, without a power flow (an ideal set becomes a realistic one).
  addMeasurement!,
  addVmMeasurement!,                      # Add voltage-magnitude measurement.
  addVaMeasurement!,                      # Add PMU voltage-angle measurement.
  addPmuPhasorMeasurement!,               # Add PMU voltage phasor (Vm + Va pair).
  addPinjMeasurement!,                    # Add active-power injection measurement.
  addQinjMeasurement!,                    # Add reactive-power injection measurement.
  addPflowMeasurement!,                   # Add active branch-flow measurement.
  addQflowMeasurement!,                   # Add reactive branch-flow measurement.
  addImagMeasurement!,                    # Add current-magnitude measurement (ampere, auxiliary): branch or shunt-bay variant (SE phases 1/2).
  addIaMeasurement!,                      # Add PMU current-phasor angle measurement (degrees, alpha-referenced like VaMeas).
  addCurrentPhasorMeasurement!,           # Add a full PMU current phasor as the ImagMeas + IaMeas pair.
  setTapEstimation!,                      # Release transformer tap positions as SE states (cascade model).
  calcMachineTrafoTapFromSE,              # Back-calculate a machine (GSU) transformer tap after an SE (mutually exclusive with a release).
  validate_topology,                      # Stage-1 topology pre-checks (advisory: findings only, never a mutation).
  test_topology_hypotheses,               # Stage-3 status-toggle hypothesis test on working copies (recommendations only).
  addShuntQMeasurement!,                  # Add shunt-bay reactive-power measurement (MVar; SE phase 2 case A).
  deriveShuntPseudoMeasurements!,         # Derive ShuntQ pseudo-measurements from bay current plus Vm (SE phase 2 case B).
  readMeasurementsCSV!,                   # Read a measurement CSV v1 file (atomic import; SE phase 5).
  writeMeasurementsCSV,                   # Write net.measurements as a measurement CSV v1 file (lossless roundtrip).
  runpf_from_se!,                         # Power flow starting from the last SE result (:se_state or :se_snapshot).
  takahashi_diag,                         # Diagonal of inv(A) from a UMFPACK factorization (shared Takahashi recursion).
  takahashi_selected_inverse,             # Selected inverse on the factor pattern (SE diagnostics Omega_ii).
  writeSEStateCSV,                        # Persist the last SE start state as a CSV artifact (SE -> PF chain).
  readSEStateCSV!,                        # Restore an SE start state from CSV and register it for runpf_from_se!.
  setShuntEstimation!,                    # Release a shunt susceptance as an estimator state (SE phase 2 case A).
  findPassiveBuses,
  addZeroInjectionMeasurements!,          # Add zero-injection constraints.
  SEResult,
  runse!,                                 # Run weighted least-squares state estimation.
  numeric_rank,
  runse_diagnostics,                      # Run SE diagnostics.
  validate_measurements,
  summarize_se_diagnostics,
  print_se_diagnostics,
  numerical_observable,
  structural_observable,
  numerical_row_redundant,
  structural_row_redundant,
  evaluate_observability_matrix,
  evaluate_local_observability_matrix,
  evaluate_global_observability,          # Evaluate global SE observability.
  measurement_jacobian,                   # Labeled measurement Jacobian H (rows/cols described) for matrix reports.
  evaluate_local_observability,           # Evaluate local SE observability.

  # rectangular network solver
  runpf_rectangular!,                     # Run rectangular complex-state NR PF.
  mismatch_rectangular,                   # Compute rectangular PF mismatch vector.
  project_rectangular_start,              # Build/select rectangular start candidate.

  # External solver interface
  PFModel,
  PFSolution,
  AbstractExternalSolver,
  buildPfModel,                           # Build external PF solver model.
  mismatchInf,                            # Compute infinity-norm mismatch.
  applyPfSolution!,                       # Apply an external PF solution to Net.
  solvePf,                                # Solve through an external solver interface.
  runpf_external!,                        # Run PF via external solver interface.
  apslf_solver,                           # Constructor for the AnalyticLoadFlow.jl-backed analytic power-series solver.
  ApslfSolver,                            # The solver type itself (analytic power series with optional Pade evaluation and NR polishing).
  ensure_casefile                         # Resolve or fetch a case file.


# ---------------------------------------------------------------------------
# Include order (stage 5): four blocks with a fixed direction of
# dependency, core -> adapters -> api -> webui. No file of an earlier
# block references a later block at DEFINITION time; runtime calls flow
# forward only through function names resolved at call time. Within each
# block the relative order of the previous layout is preserved, and the
# load-bearing adjacency comments travel with their lines.
# ---------------------------------------------------------------------------

# --- core: types, configuration, network model, solvers, SE, SC, N-1 ------
include("performance_profile.jl")
include("config/yamlparams.jl")
include("controller/control_framework.jl")
include("config/configuration.jl")
# configuration resolution (D5 precedence chain) lives with the
# configuration it resolves (stage 5; moved from src/api)
include("config/config_overrides.jl")
include("component.jl")
include("lines.jl")
include("transformer.jl")
include("prosumer.jl")
include("node.jl")
include("branch.jl")
include("link.jl")
include("shunt.jl")
# Native short-circuit source container (issue #299) — must precede
# network.jl because Net carries a typed sc_sources field.
include("shortcircuit/native_sc_data.jl")
# HVDC link record: must precede network.jl for the same reason (typed
# hvdcLinks field); the pair controller itself loads later.
include("hvdc_link.jl")
include("network.jl")
include("synthetic/synthetic_grids.jl")
include("controller/tap_control.jl")
include("controller/machine_control.jl")
include("controller/shunt_control.jl")
include("controller/series_reactance_control.jl")
# the UPFC composite registers a series plus a machine controller, so it
# loads after both device files
include("controller/upfc_control.jl")
# the full UPFC (#326) is a single multi-actuator controller; addUpfcControl!
# forwards to it for model = :full (resolved at call time, so include order
# after upfc_control.jl is fine)
include("controller/upfc_full_control.jl")
include("controller/hvdc_pair_control.jl")
include("controller/controller_config.jl")
include("busdata.jl")
include("equicircuit.jl")
include("limits.jl")
include("losses.jl")
include("numerics/condition_number.jl")
include("results.jl")
include("adapters/scf/scf.jl")
include("acpflow/acpflow.jl")
include("acpflow/island_diagnostics.jl")
include("remove_functions.jl")
include("acpflow/solver_core.jl")
# Rectangular power-flow helper layers. Keep dependency order:
# factorized linear-solver backend -> core equations -> voltage helpers -> Jacobian builders -> Newton step -> diagnostics/start/result helpers -> solver loop.
# The factorized backend comes first: it only depends on solve_sparse_system from
# solver_core.jl (its fallback chain) and is referenced by the Newton step.
include("powerflow_rectangular/rectangular_factorized_solver.jl")
include("powerflow_rectangular/rectangular_distributed_slack.jl")
include("powerflow_rectangular/rectangular_core_equations.jl")
include("powerflow_rectangular/rectangular_voltage_helpers.jl")
include("powerflow_rectangular/rectangular_jacobian_builders.jl")
include("powerflow_rectangular/rectangular_newton_step.jl")
include("powerflow_rectangular/rectangular_standalone_solver.jl")
include("powerflow_rectangular/rectangular_wrong_branch.jl")
include("powerflow_rectangular/rectangular_start_projection.jl")
include("powerflow_rectangular/rectangular_result_updates.jl")
include("powerflow_rectangular/rectangular_voltage_setpoints.jl")
include("powerflow_rectangular/rectangular_qlimit_trace.jl")
include("powerflow_rectangular/rectangular_qlimit_trace_logging.jl")
include("powerflow_rectangular/rectangular_qlimit_vset_adjustment.jl")
include("powerflow_rectangular/rectangular_qlimit_guard.jl")
include("powerflow_rectangular/rectangular_qlimit_iteration.jl")
include("powerflow_rectangular/rectangular_qlimit_outer_loop.jl")
include("powerflow_rectangular/current_iteration_start.jl")
include("powerflow_rectangular/rectangular_status_workspace.jl")
include("powerflow_rectangular/rectangular_finalization.jl")
include("powerflow_rectangular/rectangular_final_status.jl")
include("powerflow_rectangular/rectangular_diagnostics.jl")
include("powerflow_rectangular/ac_islands.jl")
include("powerflow_rectangular/rectangular_network_solver.jl")
# Standalone DC power flow (rundcpf!). Reuses solver_core.jl's shared
# slack-reduction helper and powerflow_rectangular/ac_islands.jl's island
# detection; does not hook into the AC Newton orchestration loop above.
include("powerflow_dc/dc_bmatrix.jl")
include("powerflow_dc/dc_solve.jl")
include("powerflow_dc/dc_status_workspace.jl")
include("powerflow_dc/dc_network_solver.jl")
include("powerflow_dc/dc_report.jl")
include("acpflow/solver_interface.jl")
# the AnalyticLoadFlow.jl adapter: a package extension while that package was
# a weak dependency, a normal part of the module since it became a required one
include("acpflow/apslf_solver.jl")
include("stateestimation/measurements.jl")
include("numerics/takahashi.jl")
include("stateestimation/tap_estimation.jl")
include("stateestimation/state_estimation.jl")
include("stateestimation/topology_validation.jl")
# IEC 60909 balanced short-circuit evaluation (issue #277). Consumes
# Net + the CGMES short-circuit harvest; reuses the equicircuit branch
# helpers and the AC-island decomposition — the PF solver stays untouched.
include("shortcircuit/short_circuit.jl")
# N-1 contingency batch API (multi-core Phase 4): needs the rectangular
# solver, remove functions, island detection, and the parallel runtime
# helpers included above.
include("contingency/contingency.jl")
include("scenario/patch.jl")
include("scenario/apply.jl")
include("scenario/engine.jl")

# --- adapters: format readers, converters, the registry --------------------
include("adapters/matpower/MatpowerIO.jl")
include("adapters/matpower/createnet_powermat.jl")
include("adapters/dtf/DTFImporter.jl")
include("adapters/cgmes/CGMESImporter.jl")
import .CGMESImporter: summarizeCGMES, createNetFromCGMES, importCGMES, compareWithSV, analyzeCGMES, shortCircuitCoverage, printShortCircuitCoverage, writeCGMESFiles, CGMESLineShortCircuit, cgmesLineShortCircuitData
include("adapters/cgmes/glue.jl")
include("adapters/matpower/exportMatPower.jl")
include("adapters/adapters.jl")
include("import/case_import.jl")
include("adapters/matpower/FetchMatpowerCase.jl")
using .FetchMatpowerCase: ensure_casefile, large_cases_dir
include("adapters/matpower/matpower_runner.jl")

# --- api: run services and their metadata/artifact plumbing ----------------
include("api/api_types.jl")
include("api/serialization.jl")
include("api/artifacts.jl")
include("api/run_metadata.jl")
include("api/run_api.jl")
include("api/powerflow_service.jl")
# Web UI/service short-circuit run: needs the ShortCircuitResult
# type above and the API result helpers included earlier.
include("api/run_short_circuit_service.jl")
include("api/run_import_analysis_service.jl")
# Web UI/service N-1 contingency run (#331 Phase 5): needs the contingency
# batch API above and the shared config-driven import + API result helpers.
include("api/run_contingency_service.jl")
include("api/run_state_estimation_service.jl")

# --- webui: the local browser UI plus the sysimage build tooling -----------
include("webui/webui.jl")
include("build/sysimage_builder.jl")
include("build/precompile.jl")
#! format: on
end # module Sparlectra
